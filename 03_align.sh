#!/usr/bin/env bash
#
# 03_align.sh
# Reference download -> indexing -> alignment (bwa-mem2) -> sort/index -> QC
#
# Input : trim/<sample>_1_trimmed.fq.gz, trim/<sample>_2_trimmed.fq.gz
# Output: reference/reference.fasta (+ indices), align/<sample>.bam (+ index/stats)
#
# Requirements: wget, bwa-mem2, samtools, picard, multiqc

set -euo pipefail

REF_DIR="reference"
TRIM_DIR="trim"
ALIGN_DIR="align"
THREADS=40

# URL of the reference genome FASTA (gzip-compressed, .fasta.gz). Replace with your own.
REF_URL="https://example.org/path/to/reference_genome.fasta.gz"

mkdir -p "$REF_DIR" "$ALIGN_DIR"

# ── Reference genome ────────────────────────────
wget "$REF_URL" -O "${REF_DIR}/reference.fasta.gz" --no-check-certificate
gunzip "${REF_DIR}/reference.fasta.gz"

REF="${REF_DIR}/reference.fasta"

# ── Reference indexing ──────────────────────────────────────
bwa-mem2 index "$REF"
samtools faidx "$REF"
picard CreateSequenceDictionary R="$REF" O="${REF_DIR}/reference.dict"

# Large-contig check: .bai indexes cannot cover contigs >= 2^29 bp (~537 Mb), so a .csi index is used instead.
if awk '$2 >= 536870912 { found = 1 } END { exit !found }' "${REF}.fai"; then
    BAM_INDEX_OPTS="-c"
    echo "Contig >= 2^29 bp detected: indexing BAMs with .csi"
else
    BAM_INDEX_OPTS=""
fi

# ── Alignment ────────────────────────────────────────────────
for r1 in "${TRIM_DIR}"/*_1_trimmed.fastq.gz; do
    sample=$(basename "$r1" _1_trimmed.fastq.gz)
    r2="${TRIM_DIR}/${sample}_2_trimmed.fastq.gz"
    bam="${ALIGN_DIR}/${sample}.bam"

    echo "🧬 Aligning ${sample}..."
    bwa-mem2 mem \
        -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
        "$REF" "$r1" "$r2" \
        | samtools sort -@ "$THREADS" -o "$bam"

    samtools index ${BAM_INDEX_OPTS} -@ "$THREADS" "$bam"
    samtools flagstat -@ "$THREADS" "$bam" > "${ALIGN_DIR}/${sample}_flagstat.txt"
    samtools stats -@ "$THREADS" "$bam" > "${ALIGN_DIR}/${sample}.stats"
    echo "✅ ${sample} done"
done

# ── QC report ────────────────────────────────────────────────
multiqc "$ALIGN_DIR" -n alignment_report --exclude rsem --exclude gatk -o "$ALIGN_DIR"

echo "Alignment complete. BAM files are in ${ALIGN_DIR}/"
