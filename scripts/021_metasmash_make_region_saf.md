# 021 metaSMASH region TSV + SAF

Purpose: extract each predicted BGC region from the metaSMASH `final.contigs.json` file into a metadata TSV and a SAF annotation file for featureCounts.

Run from:
```bash
cd /scratch/carreracarri/crst_cattle_ranch_soil/functional_annotation/metasmash/x1_global_coassembly
```

metaSMASH output validation:
```bash
grep '^LOCUS' *.region*.gbk | wc -l
cut -d. -f1 < <(printf '%s\n' *.region*.gbk) | sort -u | wc -l
```

Observed:
```text
10877
10875
```

Two contigs contain two predicted regions:
```bash
printf '%s\n' *.region*.gbk | sed -E 's/\.region[0-9]+\.gbk$//' | sort | uniq -c | awk '$1>1'
```

Observed:
```text
2 k127_17486739
2 k127_17758400
```

Inspect the JSON structure:
```bash
python - <<'PY'
import json
x=json.load(open("final.contigs.json"))
print(x.keys())
print(type(x.get("records")),len(x.get("records",[])))
print(x["records"][0].keys())
PY
```

Observed:
```text
dict_keys(['version', 'input_file', 'records', 'timings', 'taxon', 'schema'])
<class 'list'> 10875
dict_keys(['id', 'seq', 'features', 'name', 'description', 'dbxrefs', 'annotations', 'letter_annotations', 'areas', 'gc_content', 'modules'])
```

Create the region metadata and SAF files:
```bash
python - <<'PY'
import json,csv

x=json.load(open("final.contigs.json"))
meta="x1_antismash_regions.tsv"
saf="x1_antismash_regions.saf"

rows=[]
for r in x["records"]:
    for i,a in enumerate(r["areas"],1):
        rid=f"{r['id']}.region{i:03d}"
        cats=sorted({p["category"] for p in a["protoclusters"].values()})
        rows.append([
            rid,r["id"],a["start"]+1,a["end"],a["end"]-a["start"],
            ";".join(a["products"]),";".join(cats),len(a["protoclusters"])
        ])

with open(meta,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["region_id","contig","Start","End","Length","products","categories","n_protoclusters"])
    w.writerows(rows)

with open(saf,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["GeneID","Chr","Start","End","Strand"])
    for r in rows:
        w.writerow([r[0],r[1],r[2],r[3],"."])

print("regions:",len(rows))
print("metadata:",meta)
print("SAF:",saf)
PY
```

Coordinate conversion:
```text
metaSMASH JSON: 0-based, end-exclusive
SAF:            1-based, inclusive
Start = JSON start + 1
End   = JSON end
Length = JSON end - JSON start
```

Validate:
```bash
wc -l x1_antismash_regions.tsv x1_antismash_regions.saf
head -3 x1_antismash_regions.tsv
```

Observed:
```text
10878 x1_antismash_regions.tsv
10878 x1_antismash_regions.saf
21756 total
```

First rows:
```text
region_id	contig	Start	End	Length	products	categories	n_protoclusters
k127_1815923.region001	k127_1815923	11248	33220	21973	thioamitides	RiPP	1
k127_14704048.region001	k127_14704048	53680	78346	24667	triceptide	RiPP	1
```

Outputs:
```text
/scratch/carreracarri/crst_cattle_ranch_soil/functional_annotation/metasmash/x1_global_coassembly/x1_antismash_regions.tsv
/scratch/carreracarri/crst_cattle_ranch_soil/functional_annotation/metasmash/x1_global_coassembly/x1_antismash_regions.saf
```

BTW, SAF file is the annotation used by the next featureCounts step to count paired fragments overlapping each predicted BGC region.
