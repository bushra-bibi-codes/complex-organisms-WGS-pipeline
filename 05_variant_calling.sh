#!/usr/bin/env bash
#
# 05_variant_calling.sh
# Per-sample variant calling (HaplotypeCaller, GVCF mode) -> joint genotyping
#
# Input : bqsr/<sample>.recal.bam, reference/reference.fasta,
#         SRR_Acc_List.txt (sample list, one accession per line)
# Output: variants/<sample>.g.vcf.gz, variants/cohort.g.vcf.gz, variants/cohort.vcf.gz
#
# Requirements: gatk4, parallel (only for parallel mode)

set -euo pipefail

REF_DIR="reference"
REF="${REF_DIR}/reference.fasta"
BQSR_DIR="bqsr"
VARIANTS_DIR="variants"

# JVM tuning (GATK). Sized for a 450 GB workstation:
# PARALLEL_JOBS x JAVA_HEAP ≈ 12 x 25g = 300G, leaving ~150G headroom for the
# OS, page cache and I/O buffers. Lower JAVA_HEAP or PARALLEL_JOBS if you
# run other jobs on the same machine at the same time.
JAVA_HEAP="25g"
GC_THREADS=2
PARALLEL_JOBS=12

mkdir -p "$VARIANTS_DIR"

# ── Sequential mode: calls one sample at a time ─────────────────────────
for bam in "${BQSR_DIR}"/*.recal.bam; do
    sample=$(basename "$bam" .recal.bam)
    echo "🧬 Calling variants for ${sample}..."
    gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" HaplotypeCaller \
        -R "$REF" \
        -I "$bam" \
        -O "${VARIANTS_DIR}/${sample}.g.vcf.gz" \
        -ERC GVCF \
        --native-pair-hmm-threads "$GC_THREADS"

    gatk IndexFeatureFile -I "${VARIANTS_DIR}/${sample}.g.vcf.gz"
    echo "✅ ${sample} done"
done

# ── Parallel mode: calls multiple samples at the same time ─────────────
# Requires GNU parallel and SRR_Acc_List.txt.
# PARALLEL_JOBS x JAVA_HEAP must fit your available RAM (see note above).
#
# parallel --ungroup -j ${PARALLEL_JOBS} "
#     gatk --java-options '-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}' HaplotypeCaller \
#         -R ${REF} \
#         -I ${BQSR_DIR}/{}.recal.bam \
#         -O ${VARIANTS_DIR}/{}.g.vcf.gz \
#         -ERC GVCF \
#         --native-pair-hmm-threads ${GC_THREADS}
#     gatk IndexFeatureFile -I ${VARIANTS_DIR}/{}.g.vcf.gz
# " :::: SRR_Acc_List.txt

# ── Combine per-sample GVCFs into one cohort GVCF ───────────────────────
VARIANTS=$(awk -v dir="$VARIANTS_DIR" '{printf " -V %s/%s.g.vcf.gz", dir, $0}' SRR_Acc_List.txt)

echo "🔗 Combining GVCFs..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" CombineGVCFs \
    -R "$REF" \
    $VARIANTS \
    -O "${VARIANTS_DIR}/cohort.g.vcf.gz"

# ── Joint genotyping ─────────────────────────────────────────────────────
echo "🧬 Joint genotyping..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" GenotypeGVCFs \
    -R "$REF" \
    -V "${VARIANTS_DIR}/cohort.g.vcf.gz" \
    -O "${VARIANTS_DIR}/cohort.vcf.gz"

echo "Variant calling complete. Cohort VCF: ${VARIANTS_DIR}/cohort.vcf.gz"
