# Analysis - Minimum Reads Required for Damage Signal

By James A. Fellows Yates

## Introduction

One of the challenges in validating putative ancient dietary sequences from
dental calculus is often that there are very few molecules, making
authentication of aDNA characteristics difficult.

To illustrate this, we will map known well-preserved ancient samples of a human,
microbial pathogen, microbial commensal and a possible dietary species against
their corresponding reference genome. Then we will take the BAM files and
downscale them to small numbers of molecules that are within the ranges of
previously reported dietary studies (e.g. Warinner et al. 2014, Nat. Genet. and
Weyrich et al. 2017, Nature).

## Set up

First we will set up the conda environment containing the tools in the root
directory of this repository.

```bash

conda create env -f mann2020_damagedownscaling.yml

## Once installed
conda activate mann2020_damagedownscaling.yml
```

## Data

The data used in this analysis is as follows:

| Type | Publication | Genome | Sample | ENA Accession |
|------|-------------|--------|-------|---------------|
| Fish | Star 2017 | G_morhua | COD003 | ERR1943572
| Pathogen | Devault 2017 | G_vaginalis | Nod1-1h-nonU | SRR4885939
| Human | Gamba 2014 | H_sapiens | NE2 | SRR1187682
| Microbial Commensal | Mann et al. 2018 | T_forsythia | AO86 | SRR6877313

Sequencing Reads:

- COD003:
  [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/002/ERR1943572/ERR1943572_{1,2}.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/002/ERR1943572/ERR1943572_{1,2}.fastq.gz)
- Nod1-1h-nonU:
  [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/000/SRR4885940/SRR4885939_{1,2}.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/000/SRR4885940/SRR4885939_{1,2}.fastq.gz)
- NE2:
  [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/002/SRR1187682/SRR1187682.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/002/SRR1187682/SRR1187682.fastq.gz)
- AO86:
  [ftp.sra.ebi.ac.uk/vol1/fastq/SRR687/003/SRR6877313/SRR6877313_{1,2}.fastq.gz](ftp.sra.ebi.ac.uk/vol1/fastq/SRR687/003/SRR6877313/SRR6877313_{1,2}.fastq.gz)

Reference Genomes:

- G_morhua:
  [https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gadus_morhua/representative/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gadus_morhua/representative/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz)
- G_vaginalis:
  [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/155/GCF_000159155.2_ASM15915v2/GCF_000159155.2_ASM15915v2_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/155/GCF_000159155.2_ASM15915v2/GCF_000159155.2_ASM15915v2_genomic.fna.gz)
- H_sapiens:
  [https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
- T_forsythia:
  [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz)

## Mapping

To perform mapping of reads against a reference genome, we will use
`nf-core/eager` to do (pre)processing of raw-reads. We will then use the
resulting BAM for downstream analysis.

> Note that while we are the conda environment that  contains nextflow, we are
> using an `shh` profile that uses singularity as the container for the all the
> software. :warning: This is NOT included in this conda environment for this
> analysis! You will need to configure this accordingly for your own environment
> e.g with `conda` or `docker` as described in the documentation at
> [nf-co.re](https://nf-co.re).

```bash
mkdir -p analysis/min_reads_damage/
cd analysis/min_reads_damage/
mkdir t_forsythia g_vaginalis g_morhua h_sapiens


cd g_morhua
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/002/ERR1943572/ERR1943572_{1,2}.fastq.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gadus_morhua/representative/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'ERR1943572_{1,2}.fastq.gz' \
--fasta 'GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz'   \
--paired_end \
-name 'qi-g_morhua'

cd ../g_vaginalis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/000/SRR4885940/SRR4885939_{1,2}.fastq.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/155/GCF_000159155.2_ASM15915v2/GCF_000159155.2_ASM15915v2_genomic.fna.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'SRR4885939_{1,2}.fastq.gz' \
--fasta 'GCF_000159155.2_ASM15915v2_genomic.fna.gz'   \
--paired_end \
-name 'qi-g_vaginalis'

cd ../h_sapiens
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/002/SRR1187682/SRR1187682.fastq.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'SRR1187682.fastq.gz' \
--fasta 'hg19.fa.gz'   \
--single_end \
-name 'qi-h_sapiens'

cd ../t_forsythia
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR687/003/SRR6877313/SRR6877313_{1,2}.fastq.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'SRR6877313_{1,2}.fastq.gz' \
--fasta 'GCF_000238215.1_ASM23821v1_genomic.fna.gz' \
--paired_end \
--save_reference \
-name 'qi-t_forsythia'

cd ..

```

This will take some time.

Once complete we will place the deduplicated files in a new directory. Again
starting in the root directory of this repository.

```bash
cd analysis/min_reads_damage/
for i in g_morhua g_vaginalis h_sapiens t_forsythia; do
    mkdir -p "$i"/downsampling/raw_bam "$i"/downsampling/downsampled_bams "$i"/downsampling/damageprofiles
    cd "$i"/downsampling/raw_bam
    ln -s ../../results/deduplication/*/*bam .
    cd ../../../
done

```

## Downsampling

To use our BAMs we will use `samtools` and a GNU CoreUtils utility `shuf` (not
in the environment due to depdendncy conflicts) to randomly extract aligned
reads, and then use `DamageProfiler` to generate the damage patterns.

As the lowest number of mapped reads is in _G. vaginalis_, we will not go higher
than 10k reads. We will also do each one 100 times, so we can show how regular
the resulting pattern could occur i.e. that we weren't just unlucky to pick up a
random selection of reads that show no damage pattern (for example).

```bash
## Sampling command taken from https://www.biostars.org/p/145820/#145830
for i in g_morhua/downsampling/raw_bam/*bam; do
    for j in 25 50 100 200 500 1000 2000 5000 10000; do
        for k in {1..100}; do
        ## Extract header then the reads themselves which are subsampled
            echo "$i $j $k"
            cat <(samtools view -F 4 -SH "$i") <(samtools view -S "$i" | shuf -n "$j") | samtools view -S -b - > "${i%%.bam}"-"$j"reads-"$k"rep.bam;
        done
    done
done

for i in g_vaginalis/downsampling/raw_bam/*rmdup.bam; do
    for j in 25 50 100 200 500 1000 2000 5000 10000; do
        for k in {1..100}; do
        ## Extract header then the reads themselves which are subsampled
            echo "$i $j $k"
            cat <(samtools view -F 4 -SH "$i") <(samtools view -S "$i" | shuf -n "$j") | samtools view -S -b - > "${i%%.bam}"-"$j"reads-"$k"rep.bam;
        done
    done
done

for i in h_sapiens/downsampling/raw_bam/*rmdup.bam; do
    for j in 25 50 100 200 500 1000 2000 5000 10000; do
        for k in {1..100}; do
        ## Extract header then the reads themselves which are subsampled
            echo "$i $j $k"
            cat <(samtools view -F 4 -SH "$i") <(samtools view -S "$i" | shuf -n "$j") | samtools view -S -b - > "${i%%.bam}"-"$j"reads-"$k"rep.bam;
        done
    done
done

for i in t_forsythia/downsampling/raw_bam/*rmdup.bam; do
    for j in 25 50 100 200 500 1000 2000 5000 10000; do
        for k in {1..100}; do
        ## Extract header then the reads themselves which are subsampled
            echo "$i $j $k"
            cat <(samtools view -F 4 -SH "$i") <(samtools view -S "$i" | shuf -n "$j") | samtools view -S -b - > "${i%%.bam}"-"$j"reads-"$k"rep.bam;
        done
    done
done

```

Once generated, we can put each in a fresh folder (couldn't do that in the
previous command due to weird symlink issues))

```bash
for i in $(find -type f -name '*reads*.bam'); do
    mv "$i" "${i/raw_bam/downsampled_bams}";
done

```

Now we can generate damage profiles for each BAM file.

```bash

for i in $(find g_morhua/ -type f -name '*2000reads*rep.bam' | sort); do
    dir="$(dirname $i)"
    file="$(basename $i)"
    echo "$file"
    damageprofiler -i "$i" -o "${dir/downsampled_bams/damageprofiles}"/ -r GCF_902167405.1_gadMor3.0_genomic.fna.gz -s g_morhua
done


for i in $(find g_vaginalis/ -type f -name '*50reads*rep.bam' | sort); do
    dir="$(dirname $i)"
    file="$(basename $i)"
    echo "$file"
    damageprofiler -i "$i" -o "${dir/downsampled_bams/damageprofiles}"/ -r GCF_000159155.2_ASM15915v2_genomic.fna.gz -s g_vaginalis
done


for i in $(find h_sapiens -type f -name '*rep.bam' | sort); do
    dir="$(dirname $i)"
    file="$(basename $i)"
    echo "$file"
    damageprofiler -i "$i" -o "${dir/downsampled_bams/damageprofiles}"/ -r SRR1187682.fastq.gz -s h_sapiens
done

for i in $(find t_forsythia/ -type f -name '*28reads*rep.bam' | sort); do
    dir="$(dirname $i)"
    file="$(basename $i)"
    echo "$file"
    damageprofiler -i "$i" -o "${dir/downsampled_bams/damageprofiles}"/ -r GCF_000238215.1_ASM23821v1_genomic.fna.gz -s t_forsythia
done


```

Now the you can view the R notebook
`results/04-DamageProfilerDownsampling_Experiment_Visualisation.Rmd` to see the
commands used for the visualisation of these.
