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
ls *fa.gz | sed  's/.fa.gz//' | while read line; do kraken2 --db ../05-taxonomic_profiling/reference_databases/kraken_nt/ --threads 8 --use-names --gzip-compressed --output ../../05-taxonomic_profiling/kraken/$line.out $line.fa.gz; done
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


Merge with taxonomy string

```bash
cd kraken
gzip -d ../taxid_taxonomystr.txt.gz 
ls *summary | sed 's/.summary//' | parallel 'python ../scripts/merge_tax.py -i {}.summary -o {}.merged -t ../taxid_taxonomystr.txt'
# rm *summary
```

## Visualisation

To visualise, we can run the following R code in `scripts/figures.r`.

Make sure to load R first, and install and load the following packages

```r
install.packages("ggplot2")
install.packages("dplyr")

library(ggplot2)
library(dplyr)
```

Load the Kraken results

```r
dat <- read.table("kraken/counts_forR.txt", header=T)
dat$Level <- factor(dat$Level, levels=c("species", "genus", "other", "unassigned"))
samp.list <- levels(dat$Sample)
```

Eukaryotic lollipop charts

```r
for(i in 1:length(samp.list)){
	x <- dat %>%
		filter(Sample == samp.list[i]) %>%
		select(Level, Value)
	xlimit <- sum(x$Value)
	pdf(paste(samp.list[i], "_lolipop.pdf", sep=""))
	p <- ggplot(x, aes(Value, Level)) +
        geom_segment(aes(x=0, y=Level, xend=Value, yend=Level), size=3) +
        geom_point(size=10) + xlim(c(0, xlimit)) + theme_minimal()
    print(p)
    dev.off()
}
```

Eukaryotic pie charts

```r
for(i in 1:length(samp.list)){
	x <- dat %>%
		filter(Sample == samp.list[i]) %>%
		select(Level, Value)
	xlimit <- sum(x$Value)
	pdf(paste(samp.list[i], "_piechart.pdf", sep=""))
	p <- pie(x$Value, labels=c("Species", "Genus", "Unclassified", "Other"), col=c("#0868ac", "#43a2ca", "#bae4bc", "#7bccc4"))
    print(p)
    dev.off()
}
```

Genome avaliability correlation

```r
dat <- read.table("kraken/correlation.txt", header=T)
pdf("genome_size_species_reads.pdf")
ggplot(dat, aes(x=Num_genomes, y=log10(Value), 
	color=Spike_level, 
	shape=Species)) + 
	geom_point(size=3) + 
	scale_shape_manual(values=c(1:11)) + 
	theme_minimal() + 
	xlab("Number of genome assemblies") + 
	ylab("log10(Read count at species)")
dev.off()
cor(dat$Value, dat$Num_genomes)
# with humans removed
dat.r <- dat[-c(10,11,12),]
ggplot(dat.r, aes(x=Num_genomes, y=log10(Value), 
	color=Spike_level, 
	shape=Species)) + 
	geom_point(size=3) + 
	scale_shape_manual(values=c(1:11)) + 
	theme_minimal() + 
	xlab("Number of genome assemblies") + 
	ylab("log10(Read count at species)")
# correlation across just the 5k reads (since that might be jacking it up)
dat <- read.table("kraken/correlation_5k.txt", header=T)
cor(dat$Value, dat$Num_genomes, method="spearman")
# without human
dat.r <- dat[-c(4),]
cor(dat.r$Value, dat.r$Num_genomes, method="spearman")
# 500 reads
dat <- read.table("kraken/correlation_500.txt", header=T)
cor(dat$Value, dat$Num_genomes, method="spearman")
# without human
dat.r <- dat[-c(4),]
cor(dat.r$Value, dat.r$Num_genomes, method="spearman")
# 50 reads
dat <- read.table("kraken/correlation_50.txt", header=T)
cor(dat$Value, dat$Num_genomes, method="spearman")
# without human
dat.r <- dat[-c(4),]
cor(dat.r$Value, dat.r$Num_genomes, method="spearman")
```

Genome size corelation pplots

```r
dat <- read.table("kraken/correlation_5k_size.txt", header=T)
cor(dat$Value, dat$Genome_size.Mb.)
# [1] -0.1442868
dat <- read.table("kraken/correlation_500_size.txt", header=T)
cor(dat$Value, dat$Genome_size.Mb.)
# [1] -0.1055102
dat <- read.table("kraken/correlation_50_size.txt", header=T)
cor(dat$Value, dat$Genome_size.Mb.)
# [1] -0.1004729
```


