#!/bin/bash
#SBATCH --job-name=x1_antismash_fc
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/x1_antismash_fc_%j.out
#SBATCH --error=logs/x1_antismash_fc_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
ann="$p/functional_annotation/metasmash/x1_global_coassembly/x1_antismash_regions.saf"
out="$p/read_mapping/featurecounts/x1_antismash_region_fragment_counts.tsv"
threads="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$p/read_mapping/featurecounts" "$p/logs"

fc=$(find /home/carreracarri/mamba_envs -type f -path '*/bin/featureCounts' -print -quit)
[[ -n "$fc" ]] || { echo "featureCounts not found" >&2; exit 1; }

echo "started: $(date)"
echo "featureCounts: $("$fc" -v 2>&1 | head -1)"

bams=("$p"/read_mapping/bam/*.sorted.bam)

"$fc" \
  -T "$threads" \
  -p --countReadPairs -B -C \
  -F SAF \
  -a "$ann" \
  -o "$out" \
  "${bams[@]}"
