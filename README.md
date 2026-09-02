# Environmental Microbial Metagenomics Workflow

Analysis workflow for a shotgun metagenomic study of prairie soils from South Dakota.

This project asks how land use shapes soil microbial communities, using cattle-grazed pasture and hay/reference plots sampled repeatedly across two timepoints.

## Study questions

I organized this analysis around two main questions:

1. What microbes are present in the soil?
2. What are the microbes doing?

Kraken2 with GTDB is used for community profiling, while MEGAHIT assembly and gene annotation are used to characterize functional potential.

## Sequencing and read processing

Each soil sample was shotgun sequenced (150 bp paired-end reads, and aiming for 40 million reads per sample). Raw reads are inspected with FastQC, trimmed with fastp, and checked again after trimming. 

The resulting trimmed FASTQ files are used as input for both Kraken2 classification and Megahit assembly.

**Scripts:** `000–003`

## Taxonomic profiling

This part of the workflow determines which bacterial and archaeal taxa are represented in each sample.

Kraken2 classifies the trimmed reads against GTDB V226 (really big prokaryotic reference database).I found that the this database is more comprehensive than standard NCBI RefSeq, so its boosts our classification rates in these complex communities.

I wrote a custom script to convert my Kraken2-GTDB output into a BIOM format (standard for microbiome data). The resulting taxonomic count tables can then be summarized at the domain (Bacteria and Archaea), phylum, and species levels and compared across treatments.

**Scripts:** `004`

## Functional metagenomics

Functional analysis begins with a global coassembly of the trimmed reads using MEGAHIT. The resulting contigs are then annotated with Prodigal, which predicts protein-coding genes, and eggNOG-mapper, which annotates those predicted genes with functional information from multiple sources (PFAM, KEGG, etc.).

**Scripts:** `005–007`


## Read mapping and abundance

The global coassembly provides a common reference sequence set shared across all samples.

Bowtie2 maps the trimmed paired-end reads from each sample back to this coassembly. featureCounts then counts paired fragments associated with predicted genes or other defined genomic regions.

So assembly and annotation determine what sequences and functions exist in the dataset, and this read-mapping step allows me to quantify and compare those features across plot treatments.

**Scripts:** `009–014`, `017–018`, `022`

## Mobile genetic elements

I used geNomAD to analyze the same assembled contigs to identify sequences associated with plasmids and viruses, and if those sequences carry genes associated with antimicrobial resistance or conjugation.

**Scripts:** `008`, `015–019`

## Connecting taxonomic and functional information

Since its likely that a read can map to both a predicted eggNOG function and a Kraken2 taxonomic marker, I can use the read mapping data to link functional annotations to particular taxa. This allows me to ask which microbes are carrying which functions, and whether that changes across treatments.

**Scripts:** `013–015`

## Detecting markers of antimicrobial resistance and conjugation

I was inspired by geNomAD's AMR and conjugation annotations, so I decided to implement a similar approach using the entire predicted gene catalog.

AMRFinderPlus is used to screen the full set of Prodigal-predicted proteins for antimicrobial resistance markers. CONJScan is used to search for genes associated with conjugative transfer. This allows AMR and conjugation potential to be examined both across the full microbial community and within mobile genetic elements.

**Scripts:** `016`, `018–019`

## Biosynthetic gene clusters

I also came across metaSMASH, which identifies predicted biosynthetic gene cluster regions from the global coassembly. This could help identify what types of toxins and other secondary metabolites the community has the genetic potential to produce, and if that potential differs between treatments.

**Scripts:** `020–022`


## Repository structure



`envs/` contains environment specifications used to recreate the software environments. Will be added soon.

`scripts/` contains the numbered workflow from raw-read preparation through taxonomic profiling, assembly, annotation, read mapping, and feature quantification.

`notebooks/` contains Jupyter notebooks used for data analysis and figure-making.



