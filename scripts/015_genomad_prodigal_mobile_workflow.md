# geNomAD + Prodigal mobile-overlap workflow

Final 2022/2025 workflow for geNomAD/Prodigal overlaps, eggNOG ID matching, mobile flags, and cleaned tables.

```bash
cd /scratch/carreracarri/crst_cattle_ranch_soil
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/mapping-env
mkdir -p read_mapping/{featurecounts,overlap_checks,clean_tables}
```

The corrected Prodigal → eggNOG conversion is `contig|462_1 → contig_1`.

```bash
# count Prodigal genes
featureCounts -T 16 -p -B -C -t CDS -g ID \
  -a functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.genes.gff \
  -o read_mapping/featurecounts/x1_prodigal_gene_counts.tsv read_mapping/bam/*.sorted.bam
```

```bash
# make geNomAD SAF files
python3 - <<'PY'
import pandas as pd
from pathlib import Path

base=Path("genomad/output/x1_2022_2025_paired_global_coassembly_genomad/final.contigs_summary")
out=Path("read_mapping/featurecounts")

def find_col(cols,names):
    lower={c.lower():c for c in cols}
    for x in names:
        if x.lower() in lower: return lower[x.lower()]
    raise SystemExit(f"missing column: {names}")

def make_saf(src,dest,prefix):
    df=pd.read_csv(src,sep="\t"); c=df.columns
    gene=find_col(c,["gene","gene_id","id"]); contig=find_col(c,["contig","sequence","seq_name","seqid"])
    start=find_col(c,["start","gene_start"]); end=find_col(c,["end","gene_end"]); strand=find_col(c,["strand"])
    saf=pd.DataFrame({"GeneID":prefix+"_"+df[gene].astype(str),"Chr":df[contig].astype(str),
                      "Start":df[start].astype(int),"End":df[end].astype(int),"Strand":df[strand].astype(str)})
    saf["Strand"]=saf["Strand"].replace({"1":"+","-1":"-","0":".","nan":"."})
    saf.to_csv(dest,sep="\t",index=False)

make_saf(base/"final.contigs_virus_genes.tsv",out/"x1_genomad_virus_genes.saf","virus")
make_saf(base/"final.contigs_plasmid_genes.tsv",out/"x1_genomad_plasmid_genes.saf","plasmid")
PY
```

```bash
# count geNomAD virus/plasmid genes
featureCounts -T 16 -p -B -C -F SAF -a read_mapping/featurecounts/x1_genomad_virus_genes.saf \
  -o read_mapping/featurecounts/x1_genomad_virus_gene_counts.tsv read_mapping/bam/*.sorted.bam

featureCounts -T 16 -p -B -C -F SAF -a read_mapping/featurecounts/x1_genomad_plasmid_genes.saf \
  -o read_mapping/featurecounts/x1_genomad_plasmid_gene_counts.tsv read_mapping/bam/*.sorted.bam
```

```bash
# make BED files
awk -F'\t' '$0 !~ /^#/ && $3=="CDS" {
  split($9,a,";"); id="NA";
  for(i in a) if(a[i]~/^ID=/){id=a[i]; sub(/^ID=/,"",id)}
  print $1"\t"$4-1"\t"$5"\t"id"\t.\t"$7
}' functional_annotation/prodigal/x1_2022_2025_paired_global_coassembly.genes.gff \
> read_mapping/overlap_checks/prodigal_genes.bed

awk -F'\t' 'NR>1 {
  gene=$1; contig=gene; sub(/\*[^*]+$/,"",contig); strand=$5;
  if(strand=="1") strand="+"; else if(strand=="-1") strand="-"; else strand=".";
  print contig"\t"$2-1"\t"$3"\tvirus_"gene"\t.\t"strand
}' genomad/output/x1_2022_2025_paired_global_coassembly_genomad/final.contigs_summary/final.contigs_virus_genes.tsv \
> read_mapping/overlap_checks/genomad_virus_genes.bed

awk -F'\t' 'NR>1 {
  gene=$1; contig=gene; sub(/\*[^*]+$/,"",contig); strand=$5;
  if(strand=="1") strand="+"; else if(strand=="-1") strand="-"; else strand=".";
  print contig"\t"$2-1"\t"$3"\tplasmid_"gene"\t.\t"strand
}' genomad/output/x1_2022_2025_paired_global_coassembly_genomad/final.contigs_summary/final.contigs_plasmid_genes.tsv \
> read_mapping/overlap_checks/genomad_plasmid_genes.bed
```

```bash
# overlap Prodigal genes with geNomAD calls
bedtools intersect -a read_mapping/overlap_checks/prodigal_genes.bed \
  -b read_mapping/overlap_checks/genomad_virus_genes.bed -wa -wb \
  > read_mapping/overlap_checks/prodigal_overlap_genomad_virus.tsv

bedtools intersect -a read_mapping/overlap_checks/prodigal_genes.bed \
  -b read_mapping/overlap_checks/genomad_plasmid_genes.bed -wa -wb \
  > read_mapping/overlap_checks/prodigal_overlap_genomad_plasmid.tsv
```

```bash
# make eggNOG-style mobile flags
python3 - <<'PY'
from pathlib import Path
import pandas as pd

base=Path("read_mapping/overlap_checks")
out=base/"prodigal_eggnog_mobile_overlap_flags.tsv"

def ids(path):
    x=set()
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            full=line.rstrip("\n").split("\t")[3]
            contig,gff=full.split("|",1)
            x.add(f"{contig}_{gff.rsplit('_',1)[1]}")
    return x

virus=ids(base/"prodigal_overlap_genomad_virus.tsv")
plasmid=ids(base/"prodigal_overlap_genomad_plasmid.tsv")
all_ids=sorted(virus|plasmid)

pd.DataFrame({
    "query":all_ids,
    "overlaps_genomad_virus":[x in virus for x in all_ids],
    "overlaps_genomad_plasmid":[x in plasmid for x in all_ids],
    "overlaps_genomad_mobile":[True]*len(all_ids)
}).to_csv(out,sep="\t",index=False)
PY
```

```bash
# optional overlap check against eggNOG
python3 - <<'PY'
from pathlib import Path
import pandas as pd

base=Path("read_mapping/overlap_checks")
egg=Path("functional_annotation/eggnog/prodigal/x1_2022_2025_paired_global_coassembly_prodigal.emapper.annotations")

def ids(path):
    x=set()
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            full=line.rstrip("\n").split("\t")[3]
            contig,gff=full.split("|",1)
            x.add(f"{contig}_{gff.rsplit('_',1)[1]}")
    return x

virus=ids(base/"prodigal_overlap_genomad_virus.tsv")
plasmid=ids(base/"prodigal_overlap_genomad_plasmid.tsv")
query=set(pd.read_csv(egg,sep="\t",comment="#",header=None,low_memory=False).iloc[:,0].astype(str))

print("eggnog_annotated",len(query))
print("virus_overlap",len(query&virus))
print("plasmid_overlap",len(query&plasmid))
print("mobile_overlap",len(query&(virus|plasmid)))
PY
```

Switch to the pandas environment for cleaning:

```bash
eval "$(/home/carreracarri/micromamba shell hook -s bash)"
micromamba activate /home/carreracarri/mamba_envs/kraken-biom-env
```

```bash
# clean featureCounts tables
python3 - <<'PY'
import pandas as pd
from pathlib import Path

fc=Path("read_mapping/featurecounts"); out=Path("read_mapping/clean_tables")
files={"prodigal":fc/"x1_prodigal_gene_counts.tsv",
       "genomad_virus":fc/"x1_genomad_virus_gene_counts.tsv",
       "genomad_plasmid":fc/"x1_genomad_plasmid_gene_counts.tsv"}

for name,path in files.items():
    df=pd.read_csv(path,sep="\t",comment="#",low_memory=False)
    meta=["Geneid","Chr","Start","End","Strand","Length"]
    df=df.rename(columns={c:Path(c).name.replace(".sorted.bam","") for c in df.columns if c not in meta})
    df=df.rename(columns={"Geneid":"gene_id"})
    df.to_csv(out/f"x1_{name}_gene_counts_clean.tsv",sep="\t",index=False)
PY
```

```bash
# add mobile flags to eggNOG annotations
python3 - <<'PY'
import pandas as pd
from pathlib import Path

egg=Path("functional_annotation/eggnog/prodigal/x1_2022_2025_paired_global_coassembly_prodigal.emapper.annotations")
flags=Path("read_mapping/overlap_checks/prodigal_eggnog_mobile_overlap_flags.tsv")
out=Path("read_mapping/clean_tables/x1_prodigal_eggnog_annotations_with_mobile_flags.tsv")

header=None
with open(egg) as f:
    for line in f:
        if line.startswith("#query"):
            header=line.lstrip("#").rstrip("\n").split("\t"); break

ann=pd.read_csv(egg,sep="\t",comment="#",names=header,low_memory=False)
ann=ann.merge(pd.read_csv(flags,sep="\t"),how="left",on="query")
for col in ["overlaps_genomad_virus","overlaps_genomad_plasmid","overlaps_genomad_mobile"]:
    ann[col]=ann[col].fillna(False).astype(bool)
ann.to_csv(out,sep="\t",index=False)
PY
```

```bash
# clean geNomAD annotation tables
python3 - <<'PY'
import pandas as pd
from pathlib import Path

base=Path("genomad/output/x1_2022_2025_paired_global_coassembly_genomad/final.contigs_summary")
out=Path("read_mapping/clean_tables")

for branch in ["virus","plasmid"]:
    df=pd.read_csv(base/f"final.contigs_{branch}_genes.tsv",sep="\t",low_memory=False)
    df.insert(0,"branch",branch)
    df.insert(1,"gene_id",branch+"_"+df["gene"].astype(str))
    df.insert(2,"contig",df["gene"].astype(str).str.rsplit("_",n=1).str[0])
    df.to_csv(out/f"x1_genomad_{branch}_gene_annotations_clean.tsv",sep="\t",index=False)
PY
```

```bash
# make sample metadata
python3 - <<'PY'
import pandas as pd
from pathlib import Path

m=pd.read_csv("manifests/x1_trimmed_paired_fastqs.tsv",sep="\t")
meta=m[["sample_id","timepoint"]].copy()
meta["Treatment"]=meta["sample_id"].str.extract(r"x1_\d+_([GH])")[0].map({"G":"grazed","H":"hay"})
meta["PlotID"]=meta["sample_id"].str.extract(r"x1_\d+_([GH][A-Z])")[0]
meta["Site"]=meta["PlotID"].str[1]
meta["Time"]=meta["timepoint"].map({2022:"t1",2025:"t4"}).fillna(meta["timepoint"].astype(str))
meta.to_csv(Path("read_mapping/clean_tables/x1_sample_metadata.tsv"),sep="\t",index=False)
PY
```

Final cleaned outputs:

```text
x1_prodigal_gene_counts_clean.tsv
x1_genomad_virus_gene_counts_clean.tsv
x1_genomad_plasmid_gene_counts_clean.tsv
x1_prodigal_eggnog_annotations_with_mobile_flags.tsv
x1_genomad_virus_gene_annotations_clean.tsv
x1_genomad_plasmid_gene_annotations_clean.tsv
x1_sample_metadata.tsv
```
