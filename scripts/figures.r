###################
#ENVIRONMENT SETUP
###################
library(ggplot2)
library(dplyr)
setwd("/Users/mann/github/diet_calculus/results/reports")

####################
#LOAD & FORMAT DATA
####################
dat <- read.table("counts_forR.txt", header=T)
dat$Level <- factor(dat$Level, levels=c("Species", "Genus", "Other", "Unclassified"))
samp.list <- levels(dat$Sample)

###############################
#EUK PROPORTION LOLIPOP CHARTS
###############################
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

###########################
#EUK PROPORTION PIE CHARTS
###########################
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






