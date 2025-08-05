#!/bin/bash

sample="$1"

r1="${sample}_1_pooled.fastq.gz"
r2="${sample}_2_pooled.fastq.gz"
outdir="trinity_results/${sample}_trinity"

# Skip if already completed
if [[ -s "$fasta" || -s "$logfile" ]]; then
    echo "[$(date)] Skipping $sample — already completed." >&2
    exit 0
fi

mkdir -p "$outdir"

if [[ -f "$r1" && -f "$r2" ]]; then
        Trinity --seqType fq \
                --left "$r1" \
                --right "$r2" \
                --CPU 16 \
                --max_memory 50G \
                --output "$outdir" > trinity_logs/${sample}.log 2>&1
else
        echo "[$(date)] Missing input files for $sample" >&2
fi
