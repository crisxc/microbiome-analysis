[carreracarri@spark-login crst_cattle_ranch_soil]$ cat 016_featurecounts_proviruses.sbatch
#!/bin/bash
#SBATCH --job-name=x1_provirus_counts
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/x1_provirus_counts_%j.out
#SBATCH --error=logs/x1_provirus_counts_%j.err

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
prov="$p/genomad/output/x1_2022_2025_paired_global_coassembly_genomad/final.contigs_find_proviruses/final.contigs_provirus.tsv"
saf="$p/read_mapping/featurecounts/x1_genomad_proviruses.saf"
out="$p/read_mapping/featurecounts/x1_genomad_provirus_fragment_counts.tsv"
threads="${SLURM_CPUS_PER_TASK:-8}"

cd "$p"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env

echo "started: $(date)"
echo "featureCounts: $(featureCounts -v 2>&1 | head -1)"

awk 'BEGIN{FS=OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
NR>1 {print $1,$2,$3,$4,"."}' "$prov" > "$saf"

featureCounts -T "$threads" -p --countReadPairs -B -C -F SAF \
  -a "$saf" -o "$out" read_mapping/bam/*.sorted.bam
