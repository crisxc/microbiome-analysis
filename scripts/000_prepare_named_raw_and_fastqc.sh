#!/usr/bin/env bash
set -euo pipefail

p="/mnt/researchdrive/mzanderson/ccarri/soil/crst_cattle_ranch_soil"
raw22="/mnt/researchdrive/mzanderson/UWBC-Dropbox/Bioinformatics_Resource_Center/M008150/20250131_LH00283_0177_A22J3TYLT4"
raw25a="/mnt/researchdrive/mzanderson/UWBC-Dropbox/Bioinformatics_Resource_Center/M010264/20251223_LH00283_0285_B23GTHLLT3"
raw25b="/mnt/researchdrive/mzanderson/UWBC-Dropbox/Bioinformatics_Resource_Center/M010264/20251226_LH00283_0287_B23GN5WLT3"

n22="$p/may2022_samples/raw_named_fastqs"; qc22="$p/may2022_samples/fastqc_pre_fastp"
n25="$p/oct2025_samples/raw_named_fastqs"; qc25="$p/oct2025_samples/fastqc_pre_fastp"
log="$p/logs/x1_prepare_named_raw_and_fastqc_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$p/logs" "$n22" "$qc22" "$n25"/{main,extra,combo} "$qc25"/{main,extra,combo}
exec > >(tee -a "$log") 2>&1

echo "started: $(date)"
echo "fastqc: $(fastqc --version 2>&1 | head -1); multiqc: $(multiqc --version 2>&1 | head -1 || true)"

declare -A map22=(["CC24S001"]="GC" ["CE24S001"]="GE" ["CN24S001"]="GN" ["CS24S001"]="GS" ["CW24S001"]="GW"
                 ["HC24S001"]="HC" ["HE24S001"]="HE" ["HN24S001"]="HN" ["HS24S001"]="HS" ["HW24S001"]="HW")
declare -A map25=(["CC1"]="GC" ["CE1"]="GE" ["CN1"]="GN" ["CS1"]="GS" ["CW1"]="GW"
                 ["HC1"]="HC" ["HE1"]="HE" ["HN1"]="HN" ["HS1"]="HS" ["HW1"]="HW")

# clean names for the original runs
for old in "${!map22[@]}"; do
  code="${map22[$old]}"
  for read in R1 R2; do
    src=$(find "$raw22" -maxdepth 1 -type f -name "${old}_*_${read}_001.fastq.gz" | sort | head -1)
    [[ -n "$src" ]] || { printf "missing 2022 file: %s %s\n" "$old" "$read" >&2; exit 1; }
    ln -sfn "$src" "$n22/x1_2022_${code}_main_${read}.fastq.gz"
  done
done

for old in "${!map25[@]}"; do
  code="${map25[$old]}"
  for read in R1 R2; do
    a=$(find "$raw25a" -maxdepth 1 -type f -name "${old}_*_${read}_001.fastq.gz" | sort | head -1)
    b=$(find "$raw25b" -maxdepth 1 -type f -name "${old}_*_${read}_001.fastq.gz" | sort | head -1)
    [[ -n "$a" && -n "$b" ]] || { printf "missing 2025 file: %s %s\n" "$old" "$read" >&2; exit 1; }
    ln -sfn "$a" "$n25/main/x1_2025_${code}_main_${read}.fastq.gz"
    ln -sfn "$b" "$n25/extra/x1_2025_${code}_extra_${read}.fastq.gz"
  done
done

# main + extra were both parts of the 2025 sequencing
for old in "${!map25[@]}"; do
  code="${map25[$old]}"
  for read in R1 R2; do
    cat "$n25/main/x1_2025_${code}_main_${read}.fastq.gz" "$n25/extra/x1_2025_${code}_extra_${read}.fastq.gz" \
      > "$n25/combo/x1_2025_${code}_combo_${read}.fastq.gz"
  done
done
find "$n25/combo" -type f -name "*.fastq.gz" -print0 | xargs -0 -n1 gzip -t

# raw fastqc before trimming
fastqc -t 8 -o "$qc22" "$n22"/*.fastq.gz
fastqc -t 8 -o "$qc25/main" "$n25/main"/*.fastq.gz
fastqc -t 8 -o "$qc25/extra" "$n25/extra"/*.fastq.gz
fastqc -t 8 -o "$qc25/combo" "$n25/combo"/*.fastq.gz

if command -v multiqc >/dev/null 2>&1; then
  multiqc "$qc22" -o "$qc22" -n x1_2022_main_pre_fastp_multiqc_report.html
  multiqc "$qc25/main" -o "$qc25/main" -n x1_2025_main_pre_fastp_multiqc_report.html
  multiqc "$qc25/extra" -o "$qc25/extra" -n x1_2025_extra_pre_fastp_multiqc_report.html
  multiqc "$qc25/combo" -o "$qc25/combo" -n x1_2025_combo_pre_fastp_multiqc_report.html
fi
