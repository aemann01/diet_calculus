# Example MapDamage Figures

## Introduction

Archaeologists may not be familiar with characteristic ancient DNA patterns.
This analysis will provide simple examples of typical ancient and modern
damage profiles.

## Set up

First we will set up the conda environment containing the tools in the root
directory of this repository.

```bash

conda create env -f mann2020_damagedownscaling.yml

## Once installed
conda activate mann2020_damagedownscaling.yml
```

## Data Preparation

To generate mapDamage figures, we downloaded two samples representing an ancient
and modern human and mapped them against a human reference genome.

To download the raw sequencing data of the two samples:

```bash
## Instructions missing!
```

To download the human reference genome and unzip:

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gzip -d GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```

## Mapping

Then we can make the Bowtie2 index for mapping:

```bash
bowtie2-build GCF_000001405.39_GRCh38.p13_genomic.fna .db
```

Next we run Bowtie2 on the two samples mapping to the reference genome:

```bash
bowtie2 -x .db -U ERR2900752.collapsed.gz -S err2900752.sam --end-to-end --no-unal
bowtie2 -x .db -U modHuman.collapsed.gz -S modHuman.sam --end-to-end --no-unal
```

## DamageProfiles

Finally, we can run mapDamage on the resulting SAM files:

```bash
mapDamage -i err2900752.sam -r GCF_000001405.39_GRCh38.p13_genomic.fna -d err2900752
mapDamage -i modHuman.sam -r GCF_000001405.39_GRCh38.p13_genomic.fna -d modHuman
```
