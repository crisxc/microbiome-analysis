#!/bin/bash
#SBATCH --job-name=x1_genomad_new
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/x1_genomad_new_%j.out
#SBATCH --error=logs/x1_genomad_new_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
assembly=$(find "$p/megahit" -name "final.contigs.fa" -type f | sort | tail -1)
out="$p/genomad/output/x1_2022_2025_paired_global_coassembly_genomad"
db="/mnt/cephfs/scratch/carreracarri/x1/functional_annotation/genomad/db/genomad_db"
threads="${SLURM_CPUS_PER_TASK:-32}"

eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/genomad-1.7.6-env

echo "started: $(date)"
echo "genomad: $(genomad --version 2>&1 | head -1)"

[[ -f "$assembly" && -d "$db" ]] || exit 1
mkdir -p "$(dirname "$out")"

# run virus and plasmid detection on the coassembly
genomad end-to-end --threads "$threads" "$assembly" "$out" "$db"
