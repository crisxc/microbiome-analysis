#!/bin/bash
#SBATCH --job-name=x1_fc_fragments
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/x1_fc_fragments_%j.out
#SBATCH --error=logs/x1_fc_fragments_%j.err

set -eo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
gff="$p/functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.genes.gff"
out="$p/read_mapping/featurecounts/x1_prodigal_gene_fragment_counts.tsv"
threads="${SLURM_CPUS_PER_TASK:-16}"

cd "$p"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env

echo "started: $(date)"
echo "featureCounts: $(featureCounts -v 2>&1 | head -1)"

# count paired fragments against Prodigal CDS calls
featureCounts -T "$threads" -p --countReadPairs -B -C -t CDS -g ID -a "$gff" -o "$out"   read_mapping/bam/*.sorted.bam
