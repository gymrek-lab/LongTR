#!/bin/bash
# ============================================================================
# LongTR CI smoke test — genotype a tiny, scrubbed, PUBLIC 1000G ONT fixture and require the
# emitted VCF to be VALID: every record must survive `bcftools sort`. A plain `--help` check
# cannot catch VCF-validity regressions (e.g. an undeclared FORMAT field such as DFLANKINDEL,
# which passed `--help` but broke `bcftools sort` at the first record); this does.
#
# Fixture (this directory):
#   HG00097.bam / HG00102.bam (+ .bai)  2 public 1000G ONT samples, header-scrubbed, a few-kb slice
#   chr21.fa.gz (+ .fai/.gzi)           chr21 with real bases only in the test window, N elsewhere
#                                       (kept full-length so read/locus coordinates need no shifting)
#   loci.bed                            3 short TR loci in the window
#
# Usage:  run_smoke_test.sh /path/to/LongTR
# ============================================================================
set -euo pipefail
LONGTR="${1:?usage: run_smoke_test.sh /path/to/LongTR}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT

"$LONGTR" \
  --bams "$DIR/HG00097.bam,$DIR/HG00102.bam" --bam-samps HG00097,HG00102 --bam-libs HG00097,HG00102 \
  --fasta "$DIR/chr21.fa.gz" --regions "$DIR/loci.bed" --tr-vcf "$OUT/out.vcf.gz" \
  --min-reads 5 --phased-bam --output-gls --output-pls --output-phased-gls --output-filters \
  --alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1 \
  --log-locus-signals --log-alt-alleles \
  --aln-work-median 6.3e6 --skip-aln-work-factor 1000 --max-locus-sec 900

n=$(zcat "$OUT/out.vcf.gz" | grep -vc '^#' || true)
[ "${n:-0}" -ge 1 ] || { echo "SMOKE FAIL: LongTR produced 0 VCF records" >&2; exit 1; }

# The real guard: a strict parser must accept every record (bcftools exits non-zero otherwise).
bcftools sort "$OUT/out.vcf.gz" -Oz -o "$OUT/sorted.vcf.gz"
echo "SMOKE OK: $n records genotyped; VCF passes bcftools sort"
