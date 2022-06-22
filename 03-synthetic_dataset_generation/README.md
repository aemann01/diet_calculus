# Synthetic Dataset Generation

By Allison E. Mann

## Introduction

To demonstrate the challenges regarding metagenomic identification of dietary
reads from dental calculus, we will generate a synthetic calculus dataset.

## Setup

First we can set up the conda environment

```
conda activate mann2020_syntheticdatasetgeneration
```

Set up analysis working directories

```bash
mkdir sim mock_euks mock_oral sim/bact sim/endo sim/cont
```

## Data Downloading and Preparation

We will download the genomes of the taxa will be putting into the synthetic
dataset.

```bash
echo "Downloading data"
cd mock_euks
wget -i euk.ids -q &
cd ../examples 
wget -i examples.ids -q & 
cd ../mock_oral
wget -i bac.ids -q
echo "Download complete"
```

Next we will remove any contigs that may be contaminated.

```bash
cd ../mock_euks
ls *gz | parallel 'gzip -d {}'
grep ">" *fna | awk -F":" '{print $2}' | sed 's/>//' | awk '{print $1}' > query.txt
comm -12 <(sort query.txt) <(sort ../contaminated_genomes.Steinegger2020.txt) > contamination_hits.txt
```

We can now ask a few questions about the resulting file

How many contaminated contigs?

```bash
wc -l contamination_hits.txt
```

How many contigs do we have per organism? (Results indicated here as comments)

```bash
grep ">" *fna -c
# GCA_002197005.1_CerEla1.0_genomic.fna:11479
# GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna:22
# GCF_000001405.39_GRCh38.p13_genomic.fna:639
# GCF_000002315.6_GRCg6a_genomic.fna:464
# GCF_000003025.6_Sscrofa11.1_genomic.fna:613
# GCF_000003625.3_OryCun2.0_genomic.fna:3241
# GCF_000005005.2_B73_RefGen_v4_genomic.fna:267
# GCF_000188115.4_SL3.0_genomic.fna:3150
# GCF_000233375.1_ICSASG_v2_genomic.fna:232155
# GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna:192
# GCF_003086295.2_arahy.Tifrunner.gnm1.KYV3_genomic.fna:385
```

To remove these possibly contaminated seqs


```bash
ls *fna | sed 's/.fna//' | while read line; do awk 'BEGIN{while((getline<"contamination_hits.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $line.fna > $line.fix.fa; done
```

How much/who lost contigs?

```bash
grep ">" *fix.fa -c
# GCA_002197005.1_CerEla1.0_genomic.fix.fa:11460
# GCA_900519105.1_iwgsc_refseqv1.0_genomic.fix.fa:19
# GCF_000001405.39_GRCh38.p13_genomic.fix.fa:639
# GCF_000002315.6_GRCg6a_genomic.fix.fa:464
# GCF_000003025.6_Sscrofa11.1_genomic.fix.fa:613
# GCF_000003625.3_OryCun2.0_genomic.fix.fa:3240
# GCF_000005005.2_B73_RefGen_v4_genomic.fix.fa:267
# GCF_000188115.4_SL3.0_genomic.fix.fa:3099
# GCF_000233375.1_ICSASG_v2_genomic.fix.fa:232111
# GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fix.fa:192
# GCF_003086295.2_arahy.Tifrunner.gnm1.KYV3_genomic.fix.fa:385
```

Finally we can clean up our working directory.

```bash
#clean up
rm *fna
rename 's/.fix.fa/.fna/' *fa
```

And do the same but with the bacterial genomes

```bash
cd ../mock_oral 
ls *gz | parallel 'gzip -d {}'
grep ">" *fna | awk -F":" '{print $2}' | sed 's/>//' | awk '{print $1}' > query.txt
comm -12 <(sort query.txt) <(sort ../contaminated_genomes.Steinegger2020.txt) > contamination_hits.txt
```
This time we find no contaminated contigs.

Now we want to to label each contig of each genome is, so we know what the name
of each species is.

```bash
tagSeq()
{
    set $file
    for i in $name; do
        sed "s/^>/>${i}_/" ${1} > ${1}_fix
        shift
    done
}
name=$(awk '{print $2}' ../bac.rename)
file=$(awk '{print $1}' ../bac.rename)
tagSeq
rm *fna
cat *fix > mock_oral.fa
rm *fix
cd ../mock_euks
name=$(awk '{print $2}' ../euk.rename)
file=$(awk '{print $1}' ../euk.rename)
tagSeq

#wheat has to be processed separately (12GB file)
mv GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna_fix ../sim/endo/
rm *fna
cat *fix > mock_euks.fa
rm *fix
cd ..
```

## Simulation

Now we can use Gargammel to simulate aDNA fragmentation and metagenome creation.

Running Gargammel for generating aDNA from Eukaryotic genomes

First generate reads for wheat


```bash
cp mock_oral/mock_oral.fa sim/bact/
cp mock_oral/mock_oral.fa sim/cont/
mv sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna_fix sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna
gargammel -n 1000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/wheat_sim sim
mv sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna* mock_euks/
cp mock_euks/mock_euks.fa sim/endo/
```

Now generate reads for all other Euk. genomes

```bash
gargammel -n 5005000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_euks_sim sim
```

And finally generating aDNA from Bacterial/Archaeal genomes

```bash
rm sim/endo/mock_euks.fa
cp mock_oral/mock_oral.fa sim/endo
gargammel -n 6000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_oral_sim sim
```

Now we can move onto the next step of the analysis, which is under 04-synthetic_real_dataset_processing
