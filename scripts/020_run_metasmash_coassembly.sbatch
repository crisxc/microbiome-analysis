#!/bin/bash
#SBATCH --job-name=x1_metasmash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/x1_metasmash_%j.out
#SBATCH --error=logs/x1_metasmash_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
fa="$p/megahit/x1_2022_2025_global_coassembly_paired_trimmed/final.contigs.fa"
db="$p/functional_annotation/metasmash/databases"
out="$p/functional_annotation/metasmash/x1_global_coassembly"
threads="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "$p/logs"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /scratch/carreracarri/mamba_envs/metasmash-env

echo "started: $(date)"
echo "metasmash: $(metasmash --version 2>&1 | head -1)"

metasmash "$fa" \
  --databases "$db" \
  --output-dir "$out" \
  --cpus "$threads" \
  --workers 8 \
  --streaming on \
  --genefinding-tool prodigal-m
