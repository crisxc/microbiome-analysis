#!/bin/bash
#SBATCH --job-name=x1_prodigal_new
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/x1_prodigal_new_%j.out
#SBATCH --error=logs/x1_prodigal_new_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
assembly="$p/megahit/x1_2022_2025_global_coassembly_paired_trimmed/final.contigs.fa"
out="$p/functional_annotation/prodigal"
prefix="x1_2022_2025_paired_global_coassembly"

eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/gene-prediction-env

echo "started: $(date)"
echo "prodigal: $(prodigal -v 2>&1 | head -1); seqkit: $(seqkit version 2>&1 | head -1); pigz: $(pigz --version 2>&1 | head -1)"

[[ -s "$assembly" ]] || exit 1
mkdir -p "$out"; cd "$out"

# predict genes in metagenomic mode
prodigal -i "$assembly" -p meta -a "${prefix}.proteins.faa" -d "${prefix}.genes.fna"   -f gff -o "${prefix}.genes.gff"

pigz -p 8 -k "${prefix}.proteins.faa" "${prefix}.genes.fna" "${prefix}.genes.gff"
