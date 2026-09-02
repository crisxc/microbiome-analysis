#!/bin/bash
#SBATCH --job-name=x1_bt2_index
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=logs/x1_bt2_index_new_%j.out
#SBATCH --error=logs/x1_bt2_index_new_%j.err

set -eo pipefail

project="/scratch/carreracarri/crst_cattle_ranch_soil"
assembly="$project/megahit/x1_2022_2025_global_coassembly_paired_trimmed/final.contigs.fa"
out="$project/read_mapping/bowtie2_index"
prefix="$out/x1_2022_2025_paired_global_coassembly"
threads="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$out" "$project/logs"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env

echo "started: $(date)"
echo "bowtie2: $(bowtie2-build --version | head -1); seqkit: $(seqkit version 2>&1 | head -1)"

[[ -s "$assembly" ]] || exit 1
seqkit stats "$assembly"

# make the reusable index for all read mapping
bowtie2-build --threads "$threads" "$assembly" "$prefix"
ls -lh "${prefix}"*.bt2l
