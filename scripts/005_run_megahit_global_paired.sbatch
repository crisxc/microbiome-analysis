#!/usr/bin/env bash
#SBATCH --job-name=x1_megahit_global
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=450G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/carreracarri/crst_cattle_ranch_soil/logs/x1_megahit_global_%j.out
#SBATCH --error=/scratch/carreracarri/crst_cattle_ranch_soil/logs/x1_megahit_global_%j.err

set -eo pipefail

project="/scratch/carreracarri/crst_cattle_ranch_soil"
manifest="$project/manifests/x1_trimmed_paired_fastqs.tsv"
outdir="$project/megahit/x1_2022_2025_global_coassembly_paired_trimmed"
threads="${SLURM_CPUS_PER_TASK:-64}"; mega_mem=350000000000

mkdir -p "$project/megahit"
echo "started: $(date)"

eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/megahit-env
echo "megahit: $(megahit --version | head -1)"

# dont overwrite an old assembly
[[ -d "$outdir" ]] && exit 1

# MEGAHIT wants comma-separated R1 and R2 lists
r1=$(awk -F'\t' 'NR>1 {print $3}' "$manifest" | paste -sd "," -)
r2=$(awk -F'\t' 'NR>1 {print $4}' "$manifest" | paste -sd "," -)

megahit -1 "$r1" -2 "$r2" -t "$threads" --min-contig-len 1000 --presets meta-large \
  -m "$mega_mem" --mem-flag 0 -o "$outdir"
