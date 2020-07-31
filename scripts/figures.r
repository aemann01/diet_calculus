##################
#ENVIRONMENT SETUP
##################
library(ggplot2) # v.3.3.0
library(dplyr) # v.0.8.5
setwd("~/diet_calculus/results")
###################
#LOAD & FORMAT DATA
###################
dat <- read.table("counts_forR.txt", header=T)
dat$Level <- factor(dat$Level, levels=c("species", "genus", "other", "unassigned"))
samp.list <- levels(dat$Sample)
##############################
#EUK PROPORTION LOLIPOP CHARTS
##############################
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
##########################
#EUK PROPORTION PIE CHARTS
##########################
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
################################
#GENOME AVAILABILITY CORRELATION
################################
dat <- read.table("correlation.txt", header=T)
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
dat <- read.table("correlation_5k.txt", header=T)
cor(dat$Value, dat$Num_genomes, method="spearman")
# without human
dat.r <- dat[-c(4),]
cor(dat.r$Value, dat.r$Num_genomes, method="spearman")
# 500 reads
dat <- read.table("correlation_500.txt", header=T)
cor(dat$Value, dat$Num_genomes, method="spearman")
# without human
dat.r <- dat[-c(4),]
cor(dat.r$Value, dat.r$Num_genomes, method="spearman")
# 50 reads
dat <- read.table("correlation_50.txt", header=T)
cor(dat$Value, dat$Num_genomes, method="spearman")
# without human
dat.r <- dat[-c(4),]
cor(dat.r$Value, dat.r$Num_genomes, method="spearman")
########################
#GENOME SIZE CORRELATION
########################
dat <- read.table("correlation_5k_size.txt", header=T)
cor(dat$Value, dat$Genome_size.Mb.)
# [1] -0.1442868
dat <- read.table("correlation_500_size.txt", header=T)
cor(dat$Value, dat$Genome_size.Mb.)
# [1] -0.1055102
dat <- read.table("correlation_50_size.txt", header=T)
cor(dat$Value, dat$Genome_size.Mb.)
# [1] -0.1004729

################
#MAPDAMAGE PLOTS
################
dat <- read.table("modHuman/5pCtoT_freq.txt", header=T)
dat2 <- read.table("err2900752/5pCtoT_freq.txt", header=T)
pdf("modhuman.plot.pdf")
ggplot(dat, aes(x=pos, y=X5pC.T)) + geom_line() + ylim(c(0.0, 0.3)) + theme_minimal()
dev.off()
pdf("ancient.plot.pdf")
ggplot(dat2, aes(x=pos, y=X5pC.T)) + geom_line() + ylim(c(0.0, 0.3)) + theme_minimal()
dev.off()
