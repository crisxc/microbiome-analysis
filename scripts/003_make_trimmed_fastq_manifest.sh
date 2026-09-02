#!/usr/bin/env bash
set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
fastp="$p/combined_samples/fastp_paired_trimmed"
out="$p/manifests/x1_trimmed_paired_fastqs.tsv"

echo -e "sample_id\ttimepoint\tr1\tr2" > "$out"

for year in 2022 2025; do
  find "$fastp/$year" -mindepth 1 -maxdepth 1 -type d | sort | while read -r dir; do
    sample=$(basename "$dir"); r1="$dir/${sample}_R1.trimmed.fastq.gz"; r2="$dir/${sample}_R2.trimmed.fastq.gz"
    [[ -f "$r1" && -f "$r2" ]] || { echo "missing pair for $sample"; exit 1; }
    echo -e "$sample\t$year\t$r1\t$r2" >> "$out"
  done
done

echo "wrote: $out"
column -t "$out"
