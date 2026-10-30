#!/usr/bin/env bash
#
# 02_preprocessing.sh
# Pre-trim QC -> fastp trimming -> post-trim QC
#
# Input : raw/<sample>_1.fastq.gz, raw/<sample>_2.fastq.gz
# Output: qc/pre/, trim/, qc/post/
#
# Requirements: fastqc, fastp, multiqc

set -euo pipefail

FASTQ_DIR="raw"
TRIM_DIR="trim"
PREQC_DIR="qc/pre"
POSTQC_DIR="qc/post"
THREADS=20
# fastp is run with 20 threads for stable and efficient processing; higher thread counts may cause slow or stalled output on some datasets.

mkdir -p "$TRIM_DIR" "$PREQC_DIR" "$POSTQC_DIR"

# ── Pre-trim QC ──────────────────────────────────────────────
fastqc -t "$THREADS" -o "$PREQC_DIR" "${FASTQ_DIR}"/*.fastq.gz
multiqc "$PREQC_DIR" --exclude rsem --exclude gatk -o "$PREQC_DIR" -n pre_trim_multiqc

# ── Trim with fastp ──────────────────────────────────────────
for r1 in "${FASTQ_DIR}"/*_1.fastq.gz; do
    sample=$(basename "$r1" _1.fastq.gz)
    r2="${FASTQ_DIR}/${sample}_2.fastq.gz"
    out1="${TRIM_DIR}/${sample}_1_trimmed.fastq.gz"
    out2="${TRIM_DIR}/${sample}_2_trimmed.fastq.gz"

    [[ -f "$out1" && -f "$out2" ]] && { echo "⏭  ${sample} already trimmed"; continue; }

    echo "✂  Trimming ${sample}..."
    fastp \
        --in1 "$r1" --in2 "$r2" \
        --out1 "$out1" --out2 "$out2" \
        --detect_adapter_for_pe \
        --trim_poly_x \
        --correction \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread "$THREADS" \
        --html "${TRIM_DIR}/${sample}_fastp.html" \
        --json "${TRIM_DIR}/${sample}_fastp.json"
done

# ── Post-trim QC ─────────────────────────────────────────────
fastqc -t "$THREADS" -o "$POSTQC_DIR" "${TRIM_DIR}"/*_trimmed.fq.gz
multiqc "${TRIM_DIR}"/ -n trim_report

multiqc "$POSTQC_DIR" --exclude rsem --exclude gatk -o "$POSTQC_DIR" -n post_trim_multiqc

echo "Preprocessing complete. Trimmed reads: ${TRIM_DIR}/"