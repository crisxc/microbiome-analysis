#!/bin/bash
#SBATCH --job-name=x1_taxfun_all
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/carreracarri/crst_cattle_ranch_soil/read_mapping/taxon_function_attribution_v2/logs/x1_taxfun_all_%j.out
#SBATCH --error=/scratch/carreracarri/crst_cattle_ranch_soil/read_mapping/taxon_function_attribution_v2/logs/x1_taxfun_all_%j.err

set -eo pipefail

p="/scratch/carreracarri/crst_cattle_ranch_soil"
cd "$p"
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env

echo "started: $(date)"
echo "python: $(python3 --version 2>&1)"

python3 - <<'PY'
import csv,gc,gzip
from collections import defaultdict
from pathlib import Path

p=Path("/scratch/carreracarri/crst_cattle_ranch_soil")
base=p/"read_mapping/taxon_function_attribution_v2"; run=base/"cog_all_samples"
maps=base/"featurecounts_assignments"; out=run/"sample_tables"; ref=run/"reference"
out.mkdir(parents=True,exist_ok=True); ref.mkdir(parents=True,exist_ok=True)

bed=p/"read_mapping/overlap_checks/prodigal_genes.bed"
egg=p/"read_mapping/clean_tables/x1_prodigal_eggnog_annotations_with_mobile_flags.tsv"
md=p/"read_mapping/clean_tables/x1_sample_metadata.tsv"
lin_file=p/"read_mapping/taxon_function_sankey_all_samples/GTDBv226.taxid_lineage.tsv"
kraken_dir=p/"kraken2/gtdb_c0100_20260530_202648/classification"
gene_ref=ref/"prodigal_gene_to_cog.tsv.gz"

with md.open() as f: metadata={x["sample_id"]:x for x in csv.DictReader(f,delimiter="	")}
with lin_file.open() as f: lin={x["taxid"]:(x["domain"],x["phylum"]) for x in csv.DictReader(f,delimiter="	")}

# cache the eggNOG lookup so this part only happens once
gene_cog={}
if gene_ref.exists():
    with gzip.open(gene_ref,"rt") as f:
        for line in f:
            gene,cogs=line.rstrip("
").split("	"); gene_cog[gene]=cogs
else:
    query_cog={}
    with egg.open() as f:
        header=f.readline().rstrip("
").split("	"); q=header.index("query"); c=header.index("COG_category")
        for line in f:
            x=line.rstrip("
").split("	")
            if len(x)<=c: continue
            cogs="".join(dict.fromkeys(v for v in x[c] if "A"<=v<="Z"))
            if cogs: query_cog[x[q]]=cogs

    tmp=Path(str(gene_ref)+".tmp")
    with bed.open() as f,gzip.open(tmp,"wt",compresslevel=1) as g:
        for line in f:
            full=line.rstrip("
").split("	")[3]; contig,gene=full.split("|",1)
            query=f"{contig}_{gene.rsplit('_',1)[1]}"; cogs=query_cog.get(query,"")
            if cogs: gene_cog[gene]=cogs; g.write(f"{gene}	{cogs}
")
    tmp.replace(gene_ref); del query_cog; gc.collect()

suffix=".fragment_to_gene.tsv.gz"
for frag_file in sorted(maps.glob(f"*{suffix}")):
    sample=frag_file.name[:-len(suffix)]; done=out/f"{sample}.done"
    if done.exists(): continue

    meta=metadata[sample]; treatment=meta["Treatment"].title()
    kraken=kraken_dir/f"{sample}.gtdb.c0100.classification.txt"

    frag_cog={}; n_assigned=0
    with gzip.open(frag_file,"rt") as f:
        for line in f:
            frag,gene=line.rstrip("
").split("	",1)
            frag_cog[frag]=gene_cog.get(gene,""); n_assigned+=1

    cov=defaultdict(lambda:[0,0,0]); status=defaultdict(int); cog_count=defaultdict(int)
    total=classified=unclassified=gene_n=cog_n=0; missing=object()

    # match each Kraken fragment back to its Prodigal gene and COG
    with kraken.open() as f:
        for line in f:
            x=line.rstrip("
").split("	",3)
            if len(x)<3: continue
            call,frag,taxid=x[:3]; total+=1

            if call=="U" or taxid=="0":
                group,taxon="Kraken2_unclassified","Kraken2 unclassified"; unclassified+=1
            else:
                classified+=1; domain,phylum=lin.get(taxid,("",""))
                if domain=="Bacteria" and phylum: group,taxon="Bacterial_named_phylum",phylum
                elif domain=="Archaea" and phylum: group,taxon="Archaeal_named_phylum",phylum
                else: group,taxon="Classified_above_phylum","Classified above phylum"

            cogs=frag_cog.pop(frag,missing); cov[(group,taxon)][0]+=1
            if cogs is missing: function="Not_assigned_to_Prodigal_gene"
            else:
                gene_n+=1; cov[(group,taxon)][1]+=1
                if cogs:
                    function="Assigned_gene_with_COG"; cog_n+=1; cov[(group,taxon)][2]+=1
                    for cog in cogs: cog_count[(group,taxon,cog)]+=1
                else: function="Assigned_gene_without_COG"
            status[(group,taxon,function)]+=1

    if frag_cog: raise RuntimeError(f"{sample}: {len(frag_cog)} featureCounts IDs absent from Kraken2")
    info=[sample,treatment,meta["Time"],meta["PlotID"],meta["Site"]]

    with (out/f"{sample}.overall_summary.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="	")
        w.writerow(["Sample","Treatment","Time","PlotID","Site","TotalFragments","KrakenClassified","KrakenUnclassified","GeneAssigned","COGAnnotated","AssignedWithoutCOG","FeatureCountsRows"])
        w.writerow(info+[total,classified,unclassified,gene_n,cog_n,gene_n-cog_n,n_assigned])

    with (out/f"{sample}.taxon_function_coverage.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="	")
        w.writerow(["Sample","Treatment","Time","PlotID","Site","taxon_state","taxon","total_fragments","gene_assigned_fragments","cog_annotated_fragments","gene_assignment_pct","cog_coverage_pct"])
        for (group,taxon),x in sorted(cov.items()):
            w.writerow(info+[group,taxon,x[0],x[1],x[2],round(100*x[1]/x[0],6),round(100*x[2]/x[0],6)])

    with (out/f"{sample}.fragment_status.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="	")
        w.writerow(["Sample","Treatment","Time","PlotID","Site","taxon_state","taxon","function_state","fragments"])
        for key,n in sorted(status.items()): w.writerow(info+list(key)+[n])

    with (out/f"{sample}.taxon_cog_assignments.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="	")
        w.writerow(["Sample","Treatment","Time","PlotID","Site","taxon_state","taxon","COG_category","fragment_cog_assignments"])
        for key,n in sorted(cog_count.items()): w.writerow(info+list(key)+[n])

    done.write_text("complete
"); del frag_cog,cov,status,cog_count; gc.collect()

# combine the sample tables for the notebook
for suffix in ["overall_summary","taxon_function_coverage","fragment_status","taxon_cog_assignments"]:
    files=sorted(out.glob(f"x1_*.{suffix}.tsv")); dest=run/f"all_samples.{suffix}.tsv"; tmp=Path(str(dest)+".tmp")
    with tmp.open("w") as g:
        for i,file in enumerate(files):
            with file.open() as f:
                if i: next(f,None)
                for line in f: g.write(line)
    tmp.replace(dest)
PY
