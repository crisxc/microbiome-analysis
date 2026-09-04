#!/usr/bin/env bash
set -euo pipefail

p="/mnt/researchdrive/mzanderson/ccarri/soil/crst_cattle_ranch_soil"
in22="$p/may2022_samples/raw_named_fastqs"
in25="$p/oct2025_samples/raw_named_fastqs/combo"
out="$p/combined_samples/fastp_paired_trimmed"
log="$p/logs/x1_run_fastp_paired_trim_$(date +%Y%m%d_%H%M%S).log"
threads=8

mkdir -p "$out" "$p/logs"
exec > >(tee -a "$log") 2>&1

echo "started: $(date)"
echo "fastp: $(fastp --version 2>&1 | head -1)"

run_fastp () {
  sample="$1"; r1="$2"; r2="$3"; dir="$4"
  mkdir -p "$dir"
  [[ -e "$r1" && -e "$r2" ]] || { printf "missing input: %s\n" "$sample" >&2; exit 1; }

  fastp -i "$r1" -I "$r2" \
    -o "$dir/${sample}_R1.trimmed.fastq.gz" -O "$dir/${sample}_R2.trimmed.fastq.gz" \
    --unpaired1 "$dir/${sample}_R1.unpaired.fastq.gz" --unpaired2 "$dir/${sample}_R2.unpaired.fastq.gz" \
    --detect_adapter_for_pe -f 10 -F 10 -t 5 -T 5 -l 50 --thread "$threads" \
    --html "$dir/${sample}.fastp.html" --json "$dir/${sample}.fastp.json" \
    > "$dir/${sample}.fastp.stdout.log" 2> "$dir/${sample}.fastp.stderr.log"
}

# same trimming for 2022 main and 2025 combined runs
for code in GC GE GN GS GW HC HE HN HS HW; do
  sample="x1_2022_${code}_main"
  run_fastp "$sample" "$in22/${sample}_R1.fastq.gz" "$in22/${sample}_R2.fastq.gz" "$out/2022/$sample"
done

for code in GC GE GN GS GW HC HE HN HS HW; do
  sample="x1_2025_${code}_combo"
  run_fastp "$sample" "$in25/${sample}_R1.fastq.gz" "$in25/${sample}_R2.fastq.gz" "$out/2025/$sample"
done

find "$out" -type f \( -name "*.trimmed.fastq.gz" -o -name "*.unpaired.fastq.gz" \) -print0 | xargs -0 -n1 gzip -t
