#!/bin/bash
#SBATCH --job-name=x1_map_new
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-19
#SBATCH --output=logs/x1_map_new_%A_%a.out
#SBATCH --error=logs/x1_map_new_%A_%a.err

set -eo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
manifest="$p/manifests/x1_trimmed_paired_fastqs.tsv"
index="$p/read_mapping/bowtie2_index/x1_2022_2025_paired_global_coassembly"
out="$p/read_mapping/bam"; log="$p/read_mapping/logs"; threads="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$out" "$log" "$p/logs"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env

echo "started: $(date)"
echo "bowtie2: $(bowtie2 --version | head -1); samtools: $(samtools --version | head -1)"

line=$((SLURM_ARRAY_TASK_ID+2))
sample=$(awk -F'\t' -v n="$line" 'NR==n{print $1}' "$manifest")
r1=$(awk -F'\t' -v n="$line" 'NR==n{print $3}' "$manifest")
r2=$(awk -F'\t' -v n="$line" 'NR==n{print $4}' "$manifest")
[[ -n "$sample" && -s "$r1" && -s "$r2" ]] || exit 1

bam="$out/${sample}.sorted.bam"

# map each sample to the shared coassembly
bowtie2 -x "$index" -1 "$r1" -2 "$r2" -p "$threads" --very-sensitive 2> "$log/${sample}.bowtie2.log" |
  samtools view -@ "$threads" -bS - |
  samtools sort -@ "$threads" -o "$bam" -

samtools index "$bam"
samtools flagstat "$bam" > "$log/${sample}.flagstat.txt"
