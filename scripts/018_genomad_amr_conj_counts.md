# 018 geNomAD AMR + conjugation marker counts

Purpose: extract geNomAD AMR and conjugation marker genes from the existing plasmid and virus gene annotations, then join them to the existing geNomAD gene count tables. No remapping or new featureCounts run is needed.

Project:
```bash
p="/scratch/carreracarri/crst_cattle_ranch_soil"
```

The geNomAD plasmid and virus gene tables contain these validated annotation columns:
```text
gene
start
end
length
strand
gc_content
genetic_code
rbs_motif
marker
evalue
bitscore
uscg
plasmid_hallmark
virus_hallmark
taxid
taxname
annotation_conjscan
annotation_amr
annotation_accessions
annotation_description
```

First confirm the annotation files and exact headers:
```bash
cd "$p"

for f in \
  functional_annotation/genomad/*/final.contigs_plasmid_genes.tsv \
  functional_annotation/genomad/*/final.contigs_virus_genes.tsv
do
  echo "$f"
  head -1 "$f" | tr '\t' '\n' | nl
done
```

Locate the existing geNomAD plasmid/virus gene count tables:
```bash
find read_mapping -maxdepth 2 -type f \( \
  -iname '*genomad*virus*count*.tsv' -o \
  -iname '*genomad*plasmid*count*.tsv' \
\) -print
```

Inspect their headers before joining:
```bash
for f in $(find read_mapping -maxdepth 2 -type f \( \
  -iname '*genomad*virus*count*.tsv' -o \
  -iname '*genomad*plasmid*count*.tsv' \
\) | sort)
do
  echo "$f"
  head -1 "$f" | tr '\t' '\n' | nl | head -30
done
```

The count tables came from the existing geNomAD gene featureCounts workflow. Keep the exact existing count key. Historical geNomAD gene IDs were created from the raw `gene` value with an element prefix:
```text
plasmid_<gene>
virus_<gene>
```

Use the validated annotation columns to make AMR and conjugation marker tables and join them to the existing counts. The exact count filenames should be filled from the discovery step above rather than guessed.

Template:
```python
import pandas as pd
from pathlib import Path

p=Path("/scratch/carreracarri/crst_cattle_ranch_soil")
ann=p/"functional_annotation/genomad"
counts=p/"read_mapping/clean_tables"
out=p/"functional_annotation/genomad/marker_counts"
out.mkdir(parents=True,exist_ok=True)

# Fill these with the exact existing files found above.
files={
    "plasmid":{
        "ann":ann/"<genomad_output_dir>/final.contigs_plasmid_genes.tsv",
        "counts":counts/"<existing_plasmid_gene_count_table>.tsv"
    },
    "virus":{
        "ann":ann/"<genomad_output_dir>/final.contigs_virus_genes.tsv",
        "counts":counts/"<existing_virus_gene_count_table>.tsv"
    }
}

amr=[]; conj=[]
for src,f in files.items():
    a=pd.read_csv(f["ann"],sep="\t",low_memory=False)
    c=pd.read_csv(f["counts"],sep="\t",low_memory=False)

    a["gene_id"]=src+"_"+a["gene"].astype(str)
    a["element_type"]=src

    keep=["gene_id","element_type","gene","start","end","length","strand",
          "marker","evalue","bitscore","taxid","taxname",
          "annotation_conjscan","annotation_amr",
          "annotation_accessions","annotation_description"]

    x=a[keep].merge(c,on="gene_id",how="left")

    amr.append(x[x["annotation_amr"].notna() & x["annotation_amr"].astype(str).ne("")])
    conj.append(x[x["annotation_conjscan"].notna() & x["annotation_conjscan"].astype(str).ne("")])

pd.concat(amr,ignore_index=True).to_csv(out/"x1_genomad_amr_marker_counts.tsv",sep="\t",index=False)
pd.concat(conj,ignore_index=True).to_csv(out/"x1_genomad_conj_marker_counts.tsv",sep="\t",index=False)
```

Some points to note:
- `annotation_amr` gives geNomAD AMR marker annotations on plasmid/virus genes.
- `annotation_conjscan` gives geNomAD conjugation-related marker annotations.
- These are gene-level marker calls. They do not by themselves establish a complete conjugation system.
- Whole-community AMR is handled separately with AMRFinderPlus.
- Whole-community conjugation systems are handled separately with MacSyFinder/CONJScan.
