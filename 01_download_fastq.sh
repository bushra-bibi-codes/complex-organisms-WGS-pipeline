#!/usr/bin/env bash
#
# 01_download_fastq.sh
# Download paired-end FASTQ files from NCBI SRA using an accession list.
#
# Input : SRR_Acc_List.txt (one SRR accession per line)
# Output: raw/<SRR>_1.fastq.gz, raw/<SRR>_2.fastq.gz
#
# Requirements: sra-tools (prefetch, fasterq-dump)

set -euo pipefail

mkdir -p sra_cache raw

# ── Sequential mode: downloads one sample at a time ─────────────────────
for srr in $(cat SRR_Acc_List.txt); do
    echo "⬇  Downloading ${srr}..."
    prefetch "$srr" -O ./sra_cache
    fasterq-dump "./sra_cache/${srr}/${srr}.sra" --split-files -O ./raw/ -e 8
    gzip "./raw/${srr}_1.fastq" "./raw/${srr}_2.fastq"
    echo "✅ ${srr} done"
done

# ── Parallel mode: downloads multiple samples at the same time ─────────
# cat SRR_Acc_List.txt | xargs -P 10 -I {} bash -c '
#     echo "⬇  Downloading {}..."
#     prefetch --max-size 50G -O ./sra_cache {}
#     echo "🔄  Converting {} to FASTQ..."
#     fasterq-dump ./sra_cache/{}/{}.sra --split-files -O ./raw/ -e 4
#     gzip ./raw/{}_1.fastq ./raw/{}_2.fastq
#     echo "✅ {} done"
# '

echo "All downloads complete. FASTQ files are in raw/"
