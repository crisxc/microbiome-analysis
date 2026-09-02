[carreracarri@spark-login crst_cattle_ranch_soil]$ cat 014_run_amrfinder_all_proteins.sbatch
#!/bin/bash
#SBATCH --job-name=x1_amrfinder_all
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/x1_amrfinder_all_%j.out
#SBATCH --error=logs/x1_amrfinder_all_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
prot="$p/functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.proteins.faa"
out="$p/functional_annotation/amrfinder"
tsv="$out/x1_all_prodigal_amrfinder.tsv"
threads="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "$out" "$p/logs"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/amrfinder-env

echo "started: $(date)"
echo "amrfinder: $(amrfinder --version 2>&1 | head -1); db: $(amrfinder --database_version 2>&1 | head -1)"

[[ -s "$prot" ]] || exit 1

# screen the Prodigal protein catalog for AMR
amrfinder -p "$prot" --threads "$threads" > "$tsv"
