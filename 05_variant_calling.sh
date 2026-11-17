#!/usr/bin/env bash
#
# 05_variant_calling.sh
# Per-sample variant calling (HaplotypeCaller, GVCF mode) -> joint genotyping
#
# Input : bqsr/<sample>.recal.bam, reference/reference.fasta,
#         SRR_Acc_List.txt (sample list, one accession per line)
# Output: variants/<sample>.g.vcf.gz, variants/cohort.g.vcf.gz, variants/cohort.vcf.gz
#         (if any reference contig is >= 2^29 bp (~537 Mb), the .tbi index cannot be created and
#          GATK cannot read .csi, so no .gz files are created for GATK steps: they use uncompressed
#          .vcf files with .idx indexes, and cohort.vcf.gz is additionally written with a .csi index)
#
# Requirements: gatk4, parallel (only for parallel mode)

set -euo pipefail

REF_DIR="reference"
REF="${REF_DIR}/reference.fasta"
BQSR_DIR="bqsr"
VARIANTS_DIR="variants"
THREADS=40

# JVM tuning (GATK). Sized for a 450 GB workstation:
# PARALLEL_JOBS x JAVA_HEAP ≈ 12 x 25g = 300G, leaving ~150G headroom for the
# OS, page cache and I/O buffers. Lower JAVA_HEAP or PARALLEL_JOBS if you
# run other jobs on the same machine at the same time.
JAVA_HEAP="25g"
GC_THREADS=2
PARALLEL_JOBS=12

mkdir -p "$VARIANTS_DIR"

# Large-contig check: .tbi cannot index contigs >= 2^29 bp (~537 Mb) and GATK cannot read .csi,
# so for such references GATK steps use uncompressed VCFs (.vcf + .idx) instead of .vcf.gz.
if awk '$2 >= 536870912 { found = 1 } END { exit !found }' "${REF}.fai"; then
    LARGE_CONTIGS=true
    EXT="vcf"
    echo "Contig >= 2^29 bp detected: using uncompressed VCFs with GATK .idx indexes"
else
    LARGE_CONTIGS=false
    EXT="vcf.gz"
fi

# ── Sequential mode: calls one sample at a time ─────────────────────────
for bam in "${BQSR_DIR}"/*.recal.bam; do
    sample=$(basename "$bam" .recal.bam)
    echo "🧬 Calling variants for ${sample}..."
    gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" HaplotypeCaller \
        -R "$REF" \
        -I "$bam" \
        -O "${VARIANTS_DIR}/${sample}.g.${EXT}" \
        -ERC GVCF \
        --native-pair-hmm-threads "$GC_THREADS"

    gatk IndexFeatureFile -I "${VARIANTS_DIR}/${sample}.g.${EXT}"
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
#         -O ${VARIANTS_DIR}/{}.g.${EXT} \
#         -ERC GVCF \
#         --native-pair-hmm-threads ${GC_THREADS}
#     gatk IndexFeatureFile -I ${VARIANTS_DIR}/{}.g.${EXT}
# " :::: SRR_Acc_List.txt

# ── Combine per-sample GVCFs into one cohort GVCF ───────────────────────
VARIANTS=$(awk -v dir="$VARIANTS_DIR" -v ext="$EXT" '{printf " -V %s/%s.g.%s", dir, $0, ext}' SRR_Acc_List.txt)

echo "🔗 Combining GVCFs..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" CombineGVCFs \
    -R "$REF" \
    $VARIANTS \
    -O "${VARIANTS_DIR}/cohort.g.${EXT}"

# ── Joint genotyping ─────────────────────────────────────────────────────
echo "🧬 Joint genotyping..."
gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" GenotypeGVCFs \
    -R "$REF" \
    -V "${VARIANTS_DIR}/cohort.g.${EXT}" \
    -O "${VARIANTS_DIR}/cohort.${EXT}"

if [[ "$LARGE_CONTIGS" == true ]]; then
    echo "📦 bgzip + CSI index..."
    bgzip -f -@ "$THREADS" -c "${VARIANTS_DIR}/cohort.vcf" > "${VARIANTS_DIR}/cohort.vcf.gz"
    tabix -f -C -p vcf "${VARIANTS_DIR}/cohort.vcf.gz"
fi

echo "Variant calling complete. Cohort VCF: ${VARIANTS_DIR}/cohort.vcf.gz"