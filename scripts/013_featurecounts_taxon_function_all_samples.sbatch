#!/bin/bash
#SBATCH --job-name=x1_fc_taxfun_all
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/carreracarri/crst_cattle_ranch_soil/read_mapping/taxon_function_attribution_v2/logs/x1_fc_taxfun_all_%j.out
#SBATCH --error=/scratch/carreracarri/crst_cattle_ranch_soil/read_mapping/taxon_function_attribution_v2/logs/x1_fc_taxfun_all_%j.err

set -eo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
out="$p/read_mapping/taxon_function_attribution_v2/featurecounts_assignments"
gff="$p/functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.genes.gff"
threads="${SLURM_CPUS_PER_TASK:-16}"

eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env
cd "$out"

echo "started: $(date)"
echo "featureCounts: $(featureCounts -v 2>&1 | head -1)"

# save fragment-to-gene assignments without keeping the huge CORE files
for bam in "$p"/read_mapping/bam/*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam); map="${sample}.fragment_to_gene.tsv.gz"
    detail="${sample}.sorted.bam.featureCounts"; counts="${sample}.gene_fragment_counts.tsv"; summary="${counts}.summary"
    [[ -s "$map" ]] && continue

    if [[ ! -s "$detail" || ! -s "$summary" ]]; then
        ln -sfn "$bam" "${sample}.sorted.bam"
        featureCounts -T "$threads" -p --countReadPairs -B -C -R CORE -t CDS -g ID -a "$gff" -o "$counts" "${sample}.sorted.bam"
        unlink "${sample}.sorted.bam"
    fi

    tmp="${map}.tmp"; nfile="${sample}.assigned_rows.tmp"
    awk -F'\t' -v OFS='\t' -v nfile="$nfile" '$2=="Assigned" && $4!="" {print $1,$4; n++} END {print n+0 > nfile}' "$detail" | gzip -1 > "$tmp"
    gzip -t "$tmp"

    expected=$(awk -F'\t' '$1=="Assigned"{print $2}' "$summary"); observed=$(<"$nfile")
    [[ "$expected" == "$observed" ]] || { printf "ERROR: %s expected %s assigned fragments but found %s\n" "$sample" "$expected" "$observed" >&2; exit 1; }

    mv "$tmp" "$map"; rm -f "$detail" "$nfile"
done

# small check table for all samples
printf "Sample\tAssignedFragments\n" > all_samples.fragment_to_gene_summary.tsv
for summary in *.gene_fragment_counts.tsv.summary; do
    sample=${summary%.gene_fragment_counts.tsv.summary}
    assigned=$(awk -F'\t' '$1=="Assigned"{print $2}' "$summary")
    printf "%s\t%s\n" "$sample" "$assigned"
done >> all_samples.fragment_to_gene_summary.tsv
