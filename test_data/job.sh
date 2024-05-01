./LongTR --bams HG002_sample_reads.bam,HG003_sample_reads.bam,HG004_sample_reads.bam  \
--fasta hg38.analysisSet.fa \
--regions test_regions_hg38.bed \
--tr-vcf test.vcf.gz \
--bam-samps HG002,HG003,HG004 --bam-libs HG002,HG003,HG004 \
--min-reads 5 \
--max-tr-len 10000 \
--skip-assembly \
--phased-bam
