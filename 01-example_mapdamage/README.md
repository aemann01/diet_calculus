# Example MapDamage Figures

## Introduction

Archaeologists may not be familiar with characteristic ancient DNA patterns.
This analysis will provide simple examples of typical ancient and modern
damage profiles.

## Set up

First we will set up the conda environment containing the tools in the root
directory of this repository.

```bash

conda env create -f environment.yml

## Once installed
conda activate mann2020_examplemapdamage
```

## Data Preparation

To generate mapDamage figures, we downloaded two samples representing an ancient
and modern human and mapped them against a human reference genome.

To download the raw sequencing data of the two samples:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/002/ERR2900752/ERR2900752_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR290/002/ERR2900752/ERR2900752_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR330/004/ERR3307054/ERR3307054_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR330/004/ERR3307054/ERR3307054_2.fastq.gz
```
Reads for both samples were quality filtered and merged using AdapterRemoval before running mapdamage

```bash
AdapterRemoval --file1 ERR2900752_1.fastq.gz --file2 ERR2900752_2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename ancienthum --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
AdapterRemoval --file1 ERR3307054_1.fastq.gz --file2 ERR3307054_2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename modhuman --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
```

To download the human reference genome and unzip:

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gzip -d GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```

## Mapping

Then we can make the Bowtie2 index for mapping:

```bash
bowtie2-build GCF_000001405.39_GRCh38.p13_genomic.fna human.db
```

Next we run Bowtie2 on the two samples mapping to the reference genome:

```bash
bowtie2 -x human.db -U ERR2900752.collapsed.gz -S ancienthum.sam --end-to-end --no-unal
bowtie2 -x human.db -U modhuman.collapsed.gz -S modhuman.sam --end-to-end --no-unal
```

## DamageProfiles

Finally, we can run mapDamage on the resulting SAM files:

```bash
mapDamage -i neanderthal.sam -r GCF_000001405.39_GRCh38.p13_genomic.fna -d ancienthum
mapDamage -i modhuman.sam -r GCF_000001405.39_GRCh38.p13_genomic.fna -d modhuman
```

To generate example mapdamage plots: scripts/mapdamage_figure.ipynb
