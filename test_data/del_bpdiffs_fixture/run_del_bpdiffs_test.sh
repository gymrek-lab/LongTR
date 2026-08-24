#!/bin/bash
# ============================================================================
# LongTR CI test — BPDIFFS for the symbolic <DEL> allele must be tract-relative.
#
# <DEL> is emitted only when a read deletes the ENTIRE repeat tract: the flag comes from
# BamAlignment::TrimAlignment, whose accumulator is bounded by the tract span and is then tested
# against that same span, so the condition reduces to "deletion_size == span". LongTR therefore
# cannot emit a partial <DEL>, and every <DEL> means an allele of length 0.
#
# The bug this guards against: BPDIFFS for <DEL> was computed as -len(REF), the length of the
# EMITTED reference string. That string is the tract plus whatever flanking sequence the record
# happened to retain, and how much flank is retained depends on the reads -- so two samples that
# deleted the very same tract emitted different BPDIFFS. Consumers compute
#
#     allele_bp = (END - START + 1) + BPDIFFS
#
# and so read a clean full deletion as a NEGATIVE length, equal to minus the retained flank width.
# Downstream merging tools drop such alleles as unrecoverable, and a locus where several samples
# carry <DEL> gets a heterogeneous set of "lengths" for what is one and the same allele.
#
# THREE assertions, and the second is what keeps the third from passing vacuously:
#   1. a <DEL> allele is emitted             (the fixture still reaches the deletion path)
#   2. len(REF) > END-START+1                (the record retains flank -- the only condition under
#                                             which the old and new formulas differ at all; without
#                                             this the test would pass against the unfixed code)
#   3. (END-START+1) + BPDIFFS == 0          (the invariant, for every <DEL> allele emitted)
#
# Fixture (this directory): 2 public 1000G ONT samples over ONE chr21 locus at which reads fully
# delete the tract, plus a chr21 reference carrying real bases only where the reads land (N
# elsewhere, full length kept so coordinates need no shifting).
#
# On --indel-flank-len 20. At the default of 5 this locus emits a record whose REF is exactly the
# tract, which satisfies assertion 3 under the old code too -- the two formulas coincide whenever no
# flank is retained. Raising the flank length pulls a nearby indel into the haplotype block so the
# emitted REF extends 13 bp past the tract, which is the ~3%-of-records condition seen genome-wide.
# The flag is only the lever that makes a one-locus fixture exhibit it; nothing here is specific to
# non-default settings. Assertion 2 fails loudly if a future change stops producing the condition.
#
# Rebuilding the fixture: take the reads overlapping chr21:10405044-10405149 from the two samples,
# keep those <= 10 kb (the megabase-scale ONT reads dominate the file size and are not needed to
# reach the deletion path), and slice the reference over chr21:10380000-10425000.
#
# Usage:  run_del_bpdiffs_test.sh /path/to/LongTR
# ============================================================================
set -euo pipefail
LONGTR="${1:?usage: run_del_bpdiffs_test.sh /path/to/LongTR}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT

"$LONGTR" \
  --bams "$DIR/HG00097.bam,$DIR/HG00102.bam" --bam-samps HG00097,HG00102 --bam-libs HG00097,HG00102 \
  --fasta "$DIR/chr21.fa.gz" --regions "$DIR/loci.bed" \
  --min-reads 5 --indel-flank-len 20 \
  --tr-vcf "$OUT/del.vcf.gz" > "$OUT/del.log" 2>&1 \
  || { echo "BPDIFFS FAIL: LongTR exited non-zero" >&2; tail -20 "$OUT/del.log" >&2; exit 1; }

zcat "$OUT/del.vcf.gz" | grep -v '^#' > "$OUT/body"
[ -s "$OUT/body" ] || { echo "BPDIFFS FAIL: no VCF records emitted" >&2; exit 1; }

# Locate <DEL> among the ALT alleles and check the BPDIFFS entry at the same index. BPDIFFS is
# Number=A, so entry i corresponds to ALT i -- do not assume <DEL> is first.
awk -F'\t' '
BEGIN { del_records = 0; flank_records = 0; bad = 0 }
{
  n = split($8, info, ";"); start = ""; end = ""; bpd = ""
  for (i = 1; i <= n; i++) {
    if (info[i] ~ /^START=/)   { split(info[i], a, "="); start = a[2] }
    if (info[i] ~ /^END=/)     { split(info[i], a, "="); end   = a[2] }
    if (info[i] ~ /^BPDIFFS=/) { split(info[i], a, "="); bpd   = a[2] }
  }
  if (start == "" || end == "") {
    printf "  FAIL     record at POS=%s has no INFO/START or INFO/END\n", $2 > "/dev/stderr"; bad = 1; next
  }
  span = end - start + 1
  n_alt = split($5, alt, ",")
  n_bpd = split(bpd, diffs, ",")
  if (length($4) > span) flank_records++
  for (i = 1; i <= n_alt; i++) {
    if (alt[i] != "<DEL>") continue
    del_records++
    if (i > n_bpd) {
      printf "  FAIL     POS=%s has <DEL> as ALT %d but only %d BPDIFFS entries\n", $2, i, n_bpd > "/dev/stderr"
      bad = 1; continue
    }
    allele_bp = span + diffs[i]
    if (allele_bp != 0) {
      printf "  FAIL 3/3 POS=%s <DEL> has allele_bp = (%d-%d+1) + %d = %d, expected 0\n", \
             $2, end, start, diffs[i], allele_bp > "/dev/stderr"
      printf "           BPDIFFS is measured against the emitted REF (len %d, i.e. %d bp of retained\n", \
             length($4), length($4) - span > "/dev/stderr"
      printf "           flank) instead of the tract. It must be -(END-START+1) = %d.\n", -span > "/dev/stderr"
      bad = 1
    }
  }
}
END {
  if (del_records == 0) {
    print "  FAIL 1/3 no <DEL> allele emitted -- the fixture no longer reaches the full-deletion"  > "/dev/stderr"
    print "           path, so assertion 3 proves nothing. Rebuild the fixture (see header)."      > "/dev/stderr"
    exit 1
  }
  printf "  ok   1/3  <DEL> emitted (%d allele(s))\n", del_records
  if (flank_records == 0) {
    print "  FAIL 2/3 no record retains flank (len(REF) > END-START+1), so the old and new BPDIFFS" > "/dev/stderr"
    print "           formulas coincide and assertion 3 would pass even against the bug."          > "/dev/stderr"
    print "           Rebuild the fixture or raise --indel-flank-len (see header)."                > "/dev/stderr"
    exit 1
  }
  printf "  ok   2/3  %d record(s) retain flank in REF (the condition that exposes the bug)\n", flank_records
  if (bad) exit 1
  print  "  ok   3/3  (END-START+1) + BPDIFFS == 0 for every <DEL> allele"
}' "$OUT/body" || { echo "DEL BPDIFFS TEST FAILED" >&2; exit 1; }

echo "DEL BPDIFFS OK"
