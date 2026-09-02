#!/bin/bash
#SBATCH --job-name=x1_conjscan_all
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/x1_conjscan_all_%j.out
#SBATCH --error=logs/x1_conjscan_all_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
prot="$p/functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.proteins.faa"
models="$p/functional_annotation/conjscan/models"
out="$p/functional_annotation/conjscan/all_prodigal_chromosome"
threads="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$p/logs"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/amrfinder-env

echo "started: $(date)"
echo "macsyfinder: $(macsyfinder --version 2>&1 | head -1)"

macsyfinder \
  --db-type gembase \
  --replicon-topology linear \
  --sequence-db "$prot" \
  --models-dir "$models" \
  --models CONJScan/Chromosome all \
  --worker "$threads" \
  --out-dir "$out"
