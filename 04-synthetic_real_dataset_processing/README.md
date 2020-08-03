# Synthetic and Real Dataset Processing

## Introduction

To demonstrate the challenges regarding metagenomic identification of dietary
reads from dental calculus, we will now pre-process both the synthetic dataset
and the real dataset as you would in a typical aDNA study.

## Set Up

Create working folders

```bash
mkdir adapterremoval examples samples
```

## Data Downloading

```bash
wget -i euk.ids
wget -i examples.ids
wget -i bac.ids
```

## Read Adapter Removal, Merging and Quality Filtering

Adapter removal, quality filter, merge reads of synthetic dataset.

```bash
AdapterRemoval --file1 ../03-synthetic_dataset_generation/sim/wheat_sim_s1.fq.gz --file2 ../03-synthetic_dataset_generation/sim/wheat_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_wheat --minlength 25 &
AdapterRemoval --file1 ../03-synthetic_dataset_generation/sim/mock_euks_sim_s1.fq.gz --file2 ../03-synthetic_dataset_generation/sim/mock_euks_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_euks --minlength 25 &
AdapterRemoval --file1 ../03-synthetic_dataset_generation/sim/mock_oral_sim_s1.fq.gz --file2 ../03-synthetic_dataset_generation/sim/mock_oral_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_oral --minlength 25
```

For real samples, we first have to determine the adapter sequence

```bash
AdapterRemoval --identify-adapters --file1 examples/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq.gz --file2 examples/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq.gz
AdapterRemoval --identify-adapters --file1 examples/JAE016.A0101_R1_humfilt.fastq.gz --file2 examples/JAE016.A0101_R2_humfilt.fastq.gz
AdapterRemoval --identify-adapters --file1 examples/H10b_calc_shotgun_R1.fastq.gz --file2 examples/H10b_calc_shotgun_R2.fastq.gz
```

They all have same adapters but already started run, could do this in loop to
clean up

```bash
AdapterRemoval --file1 examples/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq.gz --file2 examples/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename neanderthal --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT &
AdapterRemoval --file1 examples/JAE016.A0101_R1_humfilt.fastq.gz --file2 examples/JAE016.A0101_R2_humfilt.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename modHuman --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
ls examples/*calc*R1* | sed 's/_R1.fastq.gz//' | sed 's/examples\///' | parallel 'AdapterRemoval --file1 examples/{}_R1.fastq.gz --file2 examples/{}_R2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename {} --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
mv *gz *settings adapterremoval
rm -r sim
```

## Dereplication

Now we can clean up the output from AdapterRemoval

```bash
#keeping collapsed, singletons, and high qualtiy forward reads
ls adapterremoval/*pair1*gz | parallel 'gzip -d {}' &
ls adapterremoval/*singleton*gz | parallel 'gzip -d {}' &
ls adapterremoval/*collapsed*gz | parallel 'gzip -d {}' 
cd adapterremoval
ls | grep -v "gz" | grep -v "settings" | parallel 'fastq_to_fasta -i {} -o {}.fna'
ls *.collapsed.fna | sed 's/.collapsed.fna//' | while read line; do cat $line.collapsed.fna $line.collapsed.truncated.fna $line.pair1.fna $line.pair1.truncated.fna $line.singletons.fna $line.singletons.truncated.fna > $line.fa; done 
rm *fna

```
And then remove unnecessary sequencing duplicates

```bash
ls *fa | sed 's/.fa//' | parallel 'vsearch --derep_fulllength {}.fa --output ../{}.uniq.fa'
cd ..
rm -r adapterremoval
mv *calc* neanderthal.uniq.fa modHuman.uniq.fa samples/
```

## Trimming

Now we can ensure there are no poly-A tail sequencing artefacts that can cause
problems in repetitive regions of genomes.

```bash
ls *uniq* | sed 's/.uniq.fa//' | while read line; do cutadapt -a "A{100}" -o $line.trim.fa $line.uniq.fa; done
```

## Dataset Spiking

Now we can spike into our synthetic microbial samples, our different synthetic
eukaryotic aDNA at different levels.

```bash
## "Splitting mock Eukaryotes by species"
rm -r temp
mkdir temp
```

First fix headers so you can retrieve after kraken

```bash
sed 's/:.*//' mock_euks.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_euks.trim.fa
sed 's/:.*//' mock_oral.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_oral.trim.fa
sed 's/:.*//' mock_wheat.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_wheat.trim.fa
```

Next split out by species

```bash
cat mock_wheat.trim.fa mock_euks.trim.fa > temp.seqs
mv temp.seqs mock_euks.trim.fa
awk '{print $2}' euk.rename | parallel 'grep {} -A 1 mock_euks.trim.fa > temp/{}.fa'
echo "How many dietary reads post quality filter and dereplication?"
grep ">" temp/*fa -c
```

Now we can make the mock oral commmunities

```bash
#5k depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4995000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_5k.fa
#500 depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999500 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_500.fa
#50 depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999950 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50.fa
```

Mock oral files are now generated. Should be: 4995000, 499500, 4999950

```bash
grep ">" *oral*5*.fa -c
```

Mock Eukaryotic community preparation

```bash
cd temp
ls *fa | sed 's/.fa$//' > diet.ids
#5000 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 5000 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.5k.fa; done &
#500 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 500 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.500.fa; done &
#50 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 50 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.50.fa; done

Finally, concatenate appropriate sets, and assign sample names

```bash
cat diet.ids | while read line; do cat $line.5k.fa ../mock_oral_5k.fa > ../samples/$line.5k.spike.fa; done &
cat diet.ids | while read line; do cat $line.500.fa ../mock_oral_500.fa > ../samples/$line.500.spike.fa; done &
cat diet.ids | while read line; do cat $line.50.fa ../mock_oral_50.fa > ../samples/$line.50.spike.fa; done

```bash
echo "Generate mock samples, should all be 5000000"
grep ">" samples/*spike*fa -c
```
And cleanup

```
rm -r temp
mv neanderthal.trim.fa *calc*trim* modHuman.trim.fa samples/
```
