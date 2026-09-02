#!/bin/bash
#SBATCH --job-name=x1_eggnog_prod
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/x1_eggnog_prodigal_%j.out
#SBATCH --error=logs/x1_eggnog_prodigal_%j.err

set -eo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
prot="$p/functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.proteins.faa"
out="$p/functional_annotation/eggnog/prodigal"
db="/mnt/cephfs/scratch/carreracarri/x1/functional_annotation/eggnog/db"
prefix="x1_2022_2025_paired_global_coassembly_prodigal"
threads="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "$out"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/eggnog-env

echo "started: $(date)"
echo "emapper: $(emapper.py --data_dir "$db" --version 2>&1 | head -1)"

[[ -s "$prot" && -d "$db" ]] || exit 1

# annotate the Prodigal protein catalog
emapper.py -i "$prot" --itype proteins --data_dir "$db" -o "$prefix"   --output_dir "$out" --cpu "$threads" --override
