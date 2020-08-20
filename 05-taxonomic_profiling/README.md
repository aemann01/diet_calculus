## Taxonomic Profiling

By Allison E. Mann

## Introduction

Now the fake and real samples have been 'spiked' with dietary DNA of known
sources, we can see how well the common taxonomic profile Kraken can find them.

## Database Preparation

First we must get and prepare the Kraken database. This may take a long time.

```bash
mkdir ~/reference_databases/kraken_nt 
kraken2-build --download-taxonomy --db kraken_nt
kraken2-build --download-library nt --db kraken_nt
kraken2-build --build --db kraken_nt --kmer-len 45 --threads 10
```

## Profiling

Now we can run each sample

```bash
ln -s 04-synthetic_real_dataset_processing/samples . && cd samples
ls *fa | parallel 'gzip {}'
ls *fa.gz | sed  's/.fa.gz//' | while read line; do kraken2 --db ~/reference_databases/kraken_nt/ --threads 8 --use-names --gzip-compressed --output ../../05-taxonomic_profiling/kraken/$line.out $line.fa.gz; done
```

## Summaries

Get summary of dietary hits for each spike level output

```bash
cd ../05-taxonomic_profiling/kraken
ls *5k* | sed 's/.5k.spike.out//' | while read line; do grep $line $line.5k.spike.out | awk -F"\t" '{print $2, "\t", $3}' | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' >> 5k.spike.summary; done
ls *5k* | sed 's/.5k.spike.out//' | while read line; do grep $line $line.500.spike.out | awk -F"\t" '{print $2, "\t", $3}' | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' >> 500.spike.summary; done
ls *5k* | sed 's/.5k.spike.out//' | while read line; do grep $line $line.50.spike.out | awk -F"\t" '{print $2, "\t", $3}' | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' >> 50.spike.summary; done
```

Bacterial summary

```bash
awk -F"\t" '{print $2, "\t", $3}' Taestivum.50.spike.out | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' > bacterial.summary
```

Summary for real samples

```bash
awk -F"\t" '{print $3}' modHuman.trim.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)//g' > modHuman.summary &
awk -F"\t" '{print $3}' neanderthal.trim.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)//g' > neanderthal.summary &
ls *calc*trim* | sed 's/.out//' | while read line; do awk -F"\t" '{print $3}' $line.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)//g' > $line.summary; done
-cd ..
ls kraken/*out | parallel 'gzip {}' & 
```

Add header line

```bash
ls kraken/5* | while read line; do sed -i '1 i\count\tquery\tsubject\ttaxid' $line; done &
sed -i '1 i\count\tquery\tsubject\ttaxid' kraken/bacterial.summary
sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' kraken/modHuman.summary
sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' kraken/neanderthal.summary
ls kraken/*calc*summary | while read line; do sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' $line; done
```


Merge with taxonomy string. taxid_taxonomystr.txt should be a two column tab-separated file where the first column is labeled taxid and the second is taxonomy. First column should be the NCBI taxonomy number, second the taxonomy string. This file can be generated from NCBI Taxonomy database [fullnamelineage.dmp](ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)

|taxid|taxonomy|
|------|-------|
|1|root|
|131567|cellular_organisms|
|2157|cellular_organisms;Archaea|
|1935019|cellular_organisms;Archaea|
|2250275|cellular_organisms;Archaea;Candidatus_Hydrothermarchaeota|

```bash
cd kraken
gzip -d taxid_taxonomystr.txt.gz 
ls *summary | sed 's/.summary//' | parallel 'python scripts/merge_tax.py -i {}.summary -o {}.merged -t taxid_taxonomystr.txt'
# rm *summary
```

Taxonomic analyses and figure generation scripts can be found in: `scripts/taxonomy_figures.ipynb`.
