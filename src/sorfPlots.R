rm(list=ls())
setwd("~/Documents/Slavoff_sPEP")
library(ggplot2)
library(ade4)
library(ggrepel)
source('src/multiplot.R')

f <- read.table('data/sORFs10000.txt',sep='\t',h=T)
g <- read.table('sEPs_Annotated.csv',sep=',',h=T,quote="")
levels(f$Annotation) = c(levels(f$Annotation),'Other')
levels(g$Annotation) = c(levels(g$Annotation),'lncrna')
f$Annotation[f$Annotation=='NMD'] = 'Other'

f$Annotation[f$Annotation=='NSD'] = 'Other'
f$Annotation[f$Annotation=='TEC'] = 'Other'
g$Annotation[g$Annotation=='5UTR and CDS'] = '5UTR'
g$Annotation[g$Annotation=='lncRNA'] = 'lncrna'
g2 <- g[1:90,]

fcount <- data.frame(cbind(levels(f$Start.codon), round(sapply(levels(f$Start.codon),function(x){
    sum(f$Start.codon==x)}, USE.NAMES = F)/length(f$Start.codon),3)))
gcount <- data.frame(cbind(levels(f$Start.codon),round(sapply(levels(f$Start.codon),function(x){
  sum(g2$Start.Type==x)},USE.NAMES = F)/length(g2$Start.Type),3)))
colnames(fcount) <- c('Start.codon','Proportion')
colnames(gcount) <- colnames(fcount)

quartz(title="sORFs",bg="white",height=10,width=15)
p1 <- ggplot(data=fcount,aes(x=Start.codon,y=Proportion)) +
  geom_bar(stat='identity',fill="black",alpha=0.6) + labs(title='sORFs Start codon')
p2 <- ggplot(data=gcount,aes(x=Start.codon,y=Proportion)) + 
  geom_bar(stat='identity',fill="black",alpha=0.6) + labs(title='sEPs Start codon')

p3 <- ggplot(data=f,aes(Sorf.length)) + 
  geom_histogram(aes(y =..density..),breaks=seq(10, 100, by = 5),col="black",fill="grey",alpha=0.6) + 
  labs(title="sORFs Histogram of peptides length",x="Length",y="Density")
p4 <- ggplot(data=g2,aes(Length)) + 
  geom_histogram(aes(y =..density..),breaks=seq(10, 150, by = 5),col="black",fill="grey",alpha=0.6) + 
  labs(title="sEPs Histogram of peptides length",x="Length",y="Density")

p5 <- ggplot(data=f,aes(Annotation)) + geom_bar(fill="black",alpha=0.6) + labs(title='sORFs Annotation')
p6 <- ggplot(data=g2,aes(Annotation)) + geom_bar(fill="black",alpha=0.6) + labs(title='sEPs Annotation')
multiplot(p1,p3,p5,p2,p4,p6,cols=2)
quartz.save(file="sORFs_subset-sEPs.pdf",type="pdf",device=dev.cur())

quartz(title="sORFs",bg="white",height = 9,width=11)
p1 <- ggplot(data=f) +
  geom_bar(aes(Start.codon, fill=factor(1),weight=rep(1/length(f$Start.codon),length(f$Start.codon))),alpha=0.5) +
  geom_bar(data=g2,aes(Start.Type,fill=factor(2),weight=rep(1/length(Start.Type),length(Start.Type))),
           alpha=0.5,inherit.aes = F) + scale_fill_discrete(name="Origin",h=c(1,150),labels=c('RP','MS')) +
  labs(title='Start codon proportions comparison',ylab='Proportion',xlab='Start codon')
p2 <- ggplot(data=f) +
  geom_bar(aes(Annotation, fill=factor(1),weight=rep(1/length(f$Start.codon),length(f$Start.codon))),alpha=0.5) +
  geom_bar(data=g2,aes(Annotation,fill=factor(2),weight=rep(1/length(Start.Type),length(Start.Type))),
           alpha=0.5,inherit.aes = F) + scale_fill_discrete(name="Origin",h=c(1,150),labels=c('RP','MS')) +
  labs(title='Location comparison',ylab='Proportion',xlab='Start codon')
p3 <- ggplot(data=f) + 
  geom_histogram(aes(x=Sorf.length,y=..density..,fill=factor(1)),breaks=seq(10, 100, by=5),alpha=0.6) + 
  geom_histogram(data=g2,aes(x=Length,y=..density..,fill=factor(2)),breaks=seq(10, 150, by=5),alpha=0.6, inherit.aes = F) + 
  labs(title="Histogram of peptides length",x="Length",y="Density") + 
  scale_fill_discrete(name="Origin",h=c(1,150),labels=c('RP','MS'))
multiplot(p1,p2,p3,cols = 1)
quartz.save(file="sORFs_subset-sEPs.pdf",type="pdf",device=dev.cur())



