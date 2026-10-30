#!/usr/bin/env bash
#
# 04_mark_duplicates_bqsr.sh
# Download known-sites VCF -> Mark PCR duplicates (Picard) -> Base Quality Score Recalibration (GATK)
#
# Input : align/<sample>.bam, reference/reference.fasta,
#         SRR_Acc_List.txt (sample list, one accession per line)
# Output: reference/known_variants.vcf.gz, markdup/<sample>.markdup.bam (+ metrics), bqsr/<sample>.recal.bam
#
# Requirements: wget, picard, samtools, gatk4, tabix, multiqc, (parallel — only for parallel mode)

set -euo pipefail

REF_DIR="reference"
REF="${REF_DIR}/reference.fasta"
KNOWN_SITES="${REF_DIR}/known_variants.vcf.gz"
ALIGN_DIR="align"
MARKDUP_DIR="markdup"
BQSR_DIR="bqsr"

# URL of a known-sites VCF (bgzipped, .vcf.gz) used for BQSR. Replace with your own.
VCF_URL="https://example.org/path/to/known_variants.vcf.gz"

# JVM tuning (Picard/GATK). Sized for a 450 GB workstation:
# PARALLEL_JOBS x JAVA_HEAP ≈ 12 x 25g = 300G, leaving ~150G headroom for the
# OS, page cache and BAM I/O buffers. Lower JAVA_HEAP or PARALLEL_JOBS if you
# run other jobs on the same machine at the same time.
JAVA_HEAP="25g"
GC_THREADS=2
PARALLEL_JOBS=12

mkdir -p "$MARKDUP_DIR" "$BQSR_DIR"

# Large-contig check: .bai indexes cannot cover contigs >= 2^29 bp (~537 Mb), so a .csi index is used instead.
if awk '$2 >= 536870912 { found = 1 } END { exit !found }' "${REF}.fai"; then
    BAM_INDEX_OPTS="-c"
    echo "Contig >= 2^29 bp detected: indexing BAMs with .csi"
else
    BAM_INDEX_OPTS=""
fi

# ── Download and index known-sites VCF (required by GATK BaseRecalibrator) ──
wget "$VCF_URL" -O "$KNOWN_SITES" --no-check-certificate
tabix -p vcf "$KNOWN_SITES"

# ── Sequential mode: processes one sample at a time ─────────────────────
for bam in "${ALIGN_DIR}"/*.bam; do
    sample=$(basename "$bam" .bam)
    markdup_bam="${MARKDUP_DIR}/${sample}.markdup.bam"
    recal_table="${BQSR_DIR}/${sample}.recal.table"
    recal_bam="${BQSR_DIR}/${sample}.recal.bam"

    echo "🔁 Marking duplicates for ${sample}..."
    picard -Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS} MarkDuplicates \
        I="$bam" \
        O="$markdup_bam" \
        M="${MARKDUP_DIR}/${sample}.markdup.metrics" \
        VALIDATION_STRINGENCY=SILENT

    samtools index ${BAM_INDEX_OPTS} "$markdup_bam"

    echo "📊 Computing BQSR recalibration table for ${sample}..."
    gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" BaseRecalibrator \
        -I "$markdup_bam" \
        -R "$REF" \
        --known-sites "$KNOWN_SITES" \
        -O "$recal_table"

    echo "🛠  Applying BQSR for ${sample}..."
    gatk --java-options "-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}" ApplyBQSR \
        -I "$markdup_bam" \
        -R "$REF" \
        --bqsr-recal-file "$recal_table" \
        -O "$recal_bam"

    samtools index ${BAM_INDEX_OPTS} "$recal_bam"
    echo "✅ ${sample} done"
done

# ── Parallel mode: processes multiple samples at the same time ─────────
# Requires GNU parallel and SRR_Acc_List.txt. Recommended for large sample
# sets / large BAMs, where MarkDuplicates and BQSR are the slowest stages.
# PARALLEL_JOBS x JAVA_HEAP must fit your available RAM (see note above).
#
# parallel --ungroup -j ${PARALLEL_JOBS} "
#     picard -Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS} MarkDuplicates \
#         I=${ALIGN_DIR}/{}.bam \
#         O=${MARKDUP_DIR}/{}.markdup.bam \
#         M=${MARKDUP_DIR}/{}.markdup.metrics \
#         VALIDATION_STRINGENCY=SILENT
#     samtools index ${BAM_INDEX_OPTS} ${MARKDUP_DIR}/{}.markdup.bam
#     gatk --java-options '-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}' BaseRecalibrator \
#         -I ${MARKDUP_DIR}/{}.markdup.bam \
#         -R ${REF} \
#         --known-sites ${KNOWN_SITES} \
#         -O ${BQSR_DIR}/{}.recal.table
#     gatk --java-options '-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=${GC_THREADS}' ApplyBQSR \
#         -I ${MARKDUP_DIR}/{}.markdup.bam \
#         -R ${REF} \
#         --bqsr-recal-file ${BQSR_DIR}/{}.recal.table \
#         -O ${BQSR_DIR}/{}.recal.bam
#     samtools index ${BAM_INDEX_OPTS} ${BQSR_DIR}/{}.recal.bam
# " :::: SRR_Acc_List.txt

# ── QC report ────────────────────────────────────────────────
multiqc "$MARKDUP_DIR" "$BQSR_DIR" -n markdup_bqsr_report --exclude rsem --exclude gatk -o "$MARKDUP_DIR"

echo "Mark duplicates and BQSR complete. Recalibrated BAMs are in ${BQSR_DIR}/"
