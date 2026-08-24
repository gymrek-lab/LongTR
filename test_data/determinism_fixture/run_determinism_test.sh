#!/bin/bash
# ============================================================================
# LongTR CI determinism test — the same input must produce the same VCF.
#
# Guards against a POA subsampling regression: HaplotypeGenerator::poa() subsamples large read
# clusters, and used a std::random_device-seeded RNG. Identical runs therefore produced different
# candidate haplotypes, different ALT alleles and different genotypes — silently, at ~0.3% of loci.
# The bug survived from 2024-02 to 2026-07 because nothing checked reproducibility.
#
# Fixture (this directory): 2 public 1000G ONT samples over ONE chr21 locus known to exercise the
# subsampling path, plus a chr21 reference with real bases only where the reads land (N elsewhere,
# full length kept so coordinates need no shifting).
#
# FOUR assertions, and the third matters as much as the first two:
#   1. same seed, twice            -> byte-identical            (the regression itself)
#   2. default settings, twice     -> byte-identical            (any other nondeterminism)
#   3. two different seeds         -> MUST DIFFER               (proves the fixture still reaches the
#                                                                subsampling path; without this the
#                                                                test passes vacuously if the fixture,
#                                                                the cluster threshold or the
#                                                                admission logic ever changes)
#   4. reversed --bams order       -> byte-identical            (input ordering must not reach the
#                                                                calls; see below)
#
# poa() is only reached when a sample has >25% of its reads unadmitted, and only subsamples when a
# cluster reaches --poa-cluster-limit (default 30). This fixture is small, so the test lowers the
# limit to 5 to reach the path; assertion 3 is what verifies that is still true.
#
# On assertion 4. LongTR indexed read groups by the order reads drained out of a std::map whose key is
# prefixed with a per-BAM label, so the internal sample order was a function of --bams order, cohort
# size and BAM layout rather than of the data. Downstream (HaplotypeGenerator::gen_candidate_seqs and
# beyond) iterates samples in that order, so calls moved: on chr21 N=10, 14 of 35,703 records differed
# between per-file and one-multi-sample-BAM input, and reversing --bams reproduced the group-BAM
# result exactly. Fixed by sorting read groups by sample name in BamProcessor::read_and_filter_reads.
#
# CAVEAT, so nobody over-trusts this: unlike assertion 3, assertion 4 is NOT proven non-vacuous on
# this fixture — its single locus did not flip on reversal even with the bug present. It is here as a
# cheap guard on a real end-to-end invariant. The fixture that does reproduce the bug is
# test_data/poa_order_fixture (10 samples, chr21:14161568-14161856); a faithful single-locus subset of
# it is ~15 MB, which is why it is not wired in here. If you touch read-group ordering or the POA
# admission path, validate against that fixture, not this assertion.
#
# Usage:  run_determinism_test.sh /path/to/LongTR
# ============================================================================
set -euo pipefail
LONGTR="${1:?usage: run_determinism_test.sh /path/to/LongTR}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT

BAMS_FWD=(--bams "$DIR/HG00097.bam,$DIR/HG00102.bam" --bam-samps HG00097,HG00102 --bam-libs HG00097,HG00102)
BAMS_REV=(--bams "$DIR/HG00102.bam,$DIR/HG00097.bam" --bam-samps HG00102,HG00097 --bam-libs HG00102,HG00097)
COMMON=(--fasta "$DIR/chr21.fa.gz" --regions "$DIR/loci.bed"
        --min-reads 5 --phased-bam --output-gls --output-pls --output-phased-gls --output-filters
        --alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1)

# run <tag> <fwd|rev> [extra flags...]   -- VCF sample columns are always sorted by sample name
# (GenotyperBamProcessor sorts samples_to_genotype_), so fwd and rev bodies are directly comparable.
run(){ local tag=$1 layout=$2; shift 2
       local bams=(); case "$layout" in
         fwd) bams=("${BAMS_FWD[@]}") ;; rev) bams=("${BAMS_REV[@]}") ;;
         *) echo "internal error: bad layout '$layout'" >&2; exit 2 ;; esac
       "$LONGTR" "${bams[@]}" "${COMMON[@]}" "$@" --tr-vcf "$OUT/$tag.vcf.gz" > "$OUT/$tag.log" 2>&1 \
         || { echo "DETERMINISM FAIL: LongTR exited non-zero for '$tag'" >&2; tail -20 "$OUT/$tag.log" >&2; exit 1; }
       zcat "$OUT/$tag.vcf.gz" | grep -v '^#' > "$OUT/$tag.body"
       [ -s "$OUT/$tag.body" ] || { echo "DETERMINISM FAIL: '$tag' produced 0 VCF records" >&2; exit 1; } }

run sub_a   fwd --poa-cluster-limit 5 --seed 42
run sub_b   fwd --poa-cluster-limit 5 --seed 42
run sub_c   fwd --poa-cluster-limit 5 --seed 7
run ord_rev rev --poa-cluster-limit 5 --seed 42
run def_a   fwd
run def_b   fwd

fail=0

if cmp -s "$OUT/sub_a.body" "$OUT/sub_b.body"; then
  echo "  ok   1/4  same seed twice -> identical ($(wc -l < "$OUT/sub_a.body") records)"
else
  echo "  FAIL 1/4  same seed twice -> DIFFERS: LongTR is not reproducible" >&2
  diff "$OUT/sub_a.body" "$OUT/sub_b.body" | head -4 | cut -c1-160 >&2
  fail=1
fi

if cmp -s "$OUT/def_a.body" "$OUT/def_b.body"; then
  echo "  ok   2/4  default settings twice -> identical"
else
  echo "  FAIL 2/4  default settings twice -> DIFFERS: nondeterminism outside the POA subsampling" >&2
  diff "$OUT/def_a.body" "$OUT/def_b.body" | head -4 | cut -c1-160 >&2
  fail=1
fi

if cmp -s "$OUT/sub_a.body" "$OUT/sub_c.body"; then
  echo "  FAIL 3/4  --seed 42 and --seed 7 agree -> the fixture no longer reaches the POA" >&2
  echo "            subsampling path, so assertions 1-2 prove nothing. Rebuild the fixture" >&2
  echo "            with a locus that exercises it, or lower --poa-cluster-limit." >&2
  fail=1
else
  echo "  ok   3/4  different seeds -> differ (subsampling path is exercised)"
fi

if cmp -s "$OUT/sub_a.body" "$OUT/ord_rev.body"; then
  echo "  ok   4/4  reversed --bams order -> identical"
else
  echo "  FAIL 4/4  reversed --bams order -> DIFFERS: input ordering is reaching the calls." >&2
  echo "            Read groups must be indexed in a canonical order, not in BAM arrival order" >&2
  echo "            (BamProcessor::read_and_filter_reads). Validate against test_data/poa_order_fixture." >&2
  diff "$OUT/sub_a.body" "$OUT/ord_rev.body" | head -4 | cut -c1-160 >&2
  fail=1
fi

[ "$fail" -eq 0 ] || { echo "DETERMINISM TEST FAILED" >&2; exit 1; }
echo "DETERMINISM OK"
