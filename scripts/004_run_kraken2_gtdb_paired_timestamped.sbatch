#!/usr/bin/env bash
#SBATCH --job-name=x1_k2_gtdb_clean
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4096
#SBATCH --time=72:00:00
#SBATCH --array=1-20%4

set -euo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
manifest="$p/manifests/x1_trimmed_paired_fastqs.tsv"
db="/scratch/carreracarri/kraken_dbs/GTDBv226"
: "${RUN_DIR:?RUN_DIR is not set}"

out="$RUN_DIR"; conf="0.10"; tag="c0100"
mkdir -p "$out"/{classification,reports,unclassified,logs}

. ~/spack/share/spack/setup-env.sh
spack env activate k2env212

IFS=$'\t' read -r sample _ r1 r2 < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$manifest")
c="$out/classification/${sample}.gtdb.${tag}.classification.txt"
r="$out/reports/${sample}.gtdb.${tag}.report.tsv"
u="$out/unclassified/${sample}.gtdb.${tag}.unclassified#.fastq"
u1="$out/unclassified/${sample}.gtdb.${tag}.unclassified_1.fastq.gz"
u2="$out/unclassified/${sample}.gtdb.${tag}.unclassified_2.fastq.gz"

echo "started: $(date)"
echo "kraken2: $(kraken2 --version 2>&1 | head -1)"

[[ -s "$c" && -s "$r" && -s "$u1" && -s "$u2" ]] && exit 0

# classify paired reads and keep the unclassified pairs
kraken2 --db "$db" --paired "$r1" "$r2" --threads "$SLURM_CPUS_PER_TASK" --memory-mapping   --report-zero-counts --confidence "$conf" --report "$r" --output "$c" --unclassified-out "$u" --gzip-compressed

for f in "$out/unclassified/${sample}.gtdb.${tag}.unclassified"_*.fastq; do
  [[ -s "$f" ]] && gzip -f "$f"
done
