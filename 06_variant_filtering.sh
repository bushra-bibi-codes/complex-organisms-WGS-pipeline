#!/usr/bin/env bash
#
# 06_variant_filtering.sh
# Split SNPs/indels -> GATK hard-filtering -> PASS extraction
#
# Input : variants/cohort.vcf.gz, reference/reference.fasta
# Output: filtered/cohort.snps.pass.vcf.gz, filtered/cohort.indels.pass.vcf.gz
#         (+ filtered/unfiltered intermediates and bcftools stats)
#
# Requirements: gatk4, bcftools, tabix, multiqc

set -euo pipefail

REF_DIR="reference"
REF="${REF_DIR}/reference.fasta"
VARIANTS_DIR="variants"
FILTERED_DIR="filtered"
COHORT_VCF="${VARIANTS_DIR}/cohort.vcf.gz"

JAVA_HEAP="200g"
GC_THREADS=10
BCFTOOLS_THREADS=40

mkdir -p "$FILTERED_DIR"

# ── Split into SNPs and indels ──────────────────────────────
echo "✂  Selecting SNPs..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" SelectVariants \
    -R "$REF" \
    -V "$COHORT_VCF" \
    --select-type-to-include SNP \
    -O "${FILTERED_DIR}/cohort.snps.vcf.gz"

echo "✂  Selecting indels..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" SelectVariants \
    -R "$REF" \
    -V "$COHORT_VCF" \
    --select-type-to-include INDEL \
    -O "${FILTERED_DIR}/cohort.indels.vcf.gz"

# ── Hard-filter SNPs (GATK best-practices thresholds) ───────
echo "🔍 Filtering SNPs..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" VariantFiltration \
    -R "$REF" \
    -V "${FILTERED_DIR}/cohort.snps.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${FILTERED_DIR}/cohort.snps.filtered.vcf.gz"

# ── Hard-filter indels (GATK best-practices thresholds) ─────
echo "🔍 Filtering indels..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" VariantFiltration \
    -R "$REF" \
    -V "${FILTERED_DIR}/cohort.indels.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "${FILTERED_DIR}/cohort.indels.filtered.vcf.gz"

# ── Summary stats ────────────────────────────────────────────
bcftools stats --threads "$BCFTOOLS_THREADS" "${FILTERED_DIR}/cohort.snps.filtered.vcf.gz" | grep "^SN" > "${FILTERED_DIR}/cohort.snps.filtered.stats.txt"
bcftools stats --threads "$BCFTOOLS_THREADS" "${FILTERED_DIR}/cohort.indels.filtered.vcf.gz" | grep "^SN" > "${FILTERED_DIR}/cohort.indels.filtered.stats.txt"

# ── Keep only PASS variants ──────────────────────────────────
echo "✅ Extracting PASS SNPs..."
bcftools view --threads "$BCFTOOLS_THREADS" -f PASS -O z \
    -o "${FILTERED_DIR}/cohort.snps.pass.vcf.gz" \
    "${FILTERED_DIR}/cohort.snps.filtered.vcf.gz"

tabix -p vcf "${FILTERED_DIR}/cohort.snps.pass.vcf.gz"

echo "✅ Extracting PASS indels..."
bcftools view --threads "$BCFTOOLS_THREADS" -f PASS -O z \
    -o "${FILTERED_DIR}/cohort.indels.pass.vcf.gz" \
    "${FILTERED_DIR}/cohort.indels.filtered.vcf.gz"

tabix -p vcf "${FILTERED_DIR}/cohort.indels.pass.vcf.gz"

# ── Final variant summary statistics ────────────────────────
bcftools stats "${FILTERED_DIR}/cohort.snps.pass.vcf.gz" > "${FILTERED_DIR}/cohort.snps.pass.stats.txt"
multiqc "$FILTERED_DIR" -n vcf_report --exclude rsem --exclude gatk -o "$FILTERED_DIR"

echo "Variant filtering complete. PASS-only VCFs:"
echo "  ${FILTERED_DIR}/cohort.snps.pass.vcf.gz"
echo "  ${FILTERED_DIR}/cohort.indels.pass.vcf.gz"
