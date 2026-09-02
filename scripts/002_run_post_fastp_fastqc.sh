#!/usr/bin/env bash
set -euo pipefail

p="/mnt/researchdrive/mzanderson/ccarri/soil/crst_cattle_ranch_soil"
fastp="$p/combined_samples/fastp_paired_trimmed"
out="$p/combined_samples/fastqc_post_fastp"
log="$p/logs/x1_run_post_fastp_fastqc_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$out"/{2022,2025} "$p/logs"
exec > >(tee -a "$log") 2>&1

echo "started: $(date)"
echo "fastqc: $(fastqc --version 2>&1 | head -1); multiqc: $(multiqc --version 2>&1 | head -1 || true)"

find "$fastp/2022" -type f -name "*.trimmed.fastq.gz" | sort > "$out/2022_trimmed_fastqs.txt"
find "$fastp/2025" -type f -name "*.trimmed.fastq.gz" | sort > "$out/2025_trimmed_fastqs.txt"

# post-trim FastQC for both years
fastqc -t 8 -o "$out/2022" $(cat "$out/2022_trimmed_fastqs.txt")
fastqc -t 8 -o "$out/2025" $(cat "$out/2025_trimmed_fastqs.txt")

if command -v multiqc >/dev/null 2>&1; then
  multiqc "$out/2022" -o "$out/2022" -n x1_2022_post_fastp_multiqc_report.html
  multiqc "$out/2025" -o "$out/2025" -n x1_2025_post_fastp_multiqc_report.html
  multiqc "$out/2022" "$out/2025" -o "$out" -n x1_2022_2025_post_fastp_multiqc_report.html
fi
