# Analysis - Minimum Reads Required for Damage Signal

## Rationale

One of the challenges in validating putative ancient dietary sequences from
dental calculus is often that there are very few molecules, making
authentication of aDNA characteristics difficult.

To illustrate this, we will map known well-preserved ancient samples of a human,
microbial pathogen, microbial commensal and a possible dietary species against
their corresponding reference genome. Then we will take the BAM files and
downscale them to small numbers of molecules that are within the ranges of
previously reported dietary studies (e.g. Warinner et al. 2014, Nat. Genet. and
Weyrich et al. 2017, Nature).


## Data

| Type | Publication | Genome | Sample | ENA Accession |
|------|-------------|--------|-------|---------------|
| Microbial Commensal | Mann et al. 2018 | T_forsythia | AO86 | SRR6877313
| Pathogen | Devault 2017 | G_vaginalis | Nod1-1h-nonU | SRR4885939 
| Fish | Star 2017 | G_morhua | COD003 | ERR1943572
| Human | Gamba 2014 | H_sapiens | NE2 | SRR1187682 

Links AO86:
ftp.sra.ebi.ac.uk/vol1/fastq/SRR687/003/SRR6877313/SRR6877313_{1,2}.fastq.gz
Nod1-1h-nonU:
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/000/SRR4885940/SRR4885939_{1,2}.fastq.gz
COD003:
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/002/ERR1943572/ERR1943572_{1,2}.fastq.gz
NE2:
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/002/SRR1187682/SRR1187682.fastq.gz

Reference Genomes: T_forsythia:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz
G_vaginalis:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/155/GCF_000159155.2_ASM15915v2/GCF_000159155.2_ASM15915v2_genomic.fna.gz
G_morhua:
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gadus_morhua/representative/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz
H_sapiens: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

## Mapping

To perform mapping of reads against a reference genome, we will use
nf-core/eager to do (pre)processing of raw-reads. We will then use the resulting
BAM for downstream analysis.

> Note that while we are the conda environment that  contains nextflow, we are
>  using an `shh` profile that uses singularity as the container for the all the
>  software. You will need to configure this accordingly for your own
>  environment e.g with `conda` or `docker` as described in teh documentation at
>  [nf-co.re](https://nf-co.re).

```bash
mkdir -p analysis/min_reads_damage
mkdir t_forsythia g_vaginalis g_morhua h_sapiens


cd t_forsythia
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR687/003/SRR6877313/SRR6877313_{1,2}.fastq.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'SRR6877313_{1,2}.fastq.gz' \
--fasta 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz' \
--paired_end \
--save_reference \
-name 'qi-t_forsythia'

cd ../g_vaginalis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/000/SRR4885940/SRR4885939_{1,2}.fastq.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'SRR4885939_{1,2}.fastq.gz' \
--fasta 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/155/GCF_000159155.2_ASM15915v2/GCF_000159155.2_ASM15915v2_genomic.fna.gz'   \
--paired_end \
-name 'qi-g_vaginalis'

cd ../g_morhua
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/002/ERR1943572/ERR1943572_{1,2}.fastq.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'ERR1943572_{1,2}.fastq.gz' \
--fasta 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gadus_morhua/representative/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz'   \
--paired_end \
-name 'qi-g_morhua'

cd ../h_sapiens
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/002/SRR1187682/SRR1187682.fastq.gz
nextflow run nf-core/eager -r 2.1.0 \
-profile shh,sdag \
-with-tower \
--email fellows@shh.mpg.de \
--reads 'SRR1187682.fastq.gz' \
--fasta 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'   \
--single_end \
-name 'qi-h_sapiens'

cd ..

```

Once complete we will place the deduplicated files in a new directory.

```bash
cd analysis/min_reads_damage/
for i in t_forsythia g_vaginalis g_morhua h_sapiens; do
    mkdir -p "$i"/downsampling/raw_bam "$i"/downsampling/downsampled_bams "$i"/downsampling/damageprofiles
    cd "$i"/downsampling/raw_bam
    ln -s ../../results/deduplication/*/*bam .
    cd ../../../
done


f
```


## Downsampling

To use our BAMs we will use `samtools` and a GNU CoreUtils utility `shuf` (not
in the environment due to depdendncy conflicts) to randomly extract reads, and
then use `DamageProfiler` to generate the damage patterns. 

```
## Taken from
cat <(samtools view -SH input.sam) <(samtools view -S input.sam | shuf -n 1000000) > output.sam
```