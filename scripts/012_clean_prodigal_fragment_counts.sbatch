#!/bin/bash
#SBATCH --job-name=x1_clean_fragments
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/x1_clean_fragments_%j.out
#SBATCH --error=logs/x1_clean_fragments_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
in="$p/read_mapping/featurecounts/x1_prodigal_gene_fragment_counts.tsv"
out="$p/read_mapping/clean_tables/x1_prodigal_gene_fragment_counts_clean.tsv"
tmp="${out}.tmp"

mkdir -p "$(dirname "$out")" "$p/logs"
[[ -s "$in" ]] || exit 1

echo "started: $(date)"
echo "awk: $(awk -W version 2>&1 | head -1)"

# clean featureCounts headers for the notebook
awk 'BEGIN{FS=OFS="\t"} /^#/ {next}
$1=="Geneid" {$1="gene_id"; for(i=7;i<=NF;i++){sub(/^.*\//,"",$i); sub(/\.sorted\.bam$/,"",$i)}}
{print}' "$in" > "$tmp"

mv "$tmp" "$out"
