rm(list=ls())
library(ggplot2)
library(ade4)
library(ggrepel)
library(grid)
library(gridExtra)
args <- commandArgs(trailingOnly = T)

# If not used on command line, modify the following line
args <- c('~/Documents/Slavoff_sPEP/ALL','~/Documents/Slavoff_sPEP/ALL/sEP_table.csv',
          '~/Documents/Slavoff_sPEP/data/aaproperties.tsv')

###############################################################################################
#     Open files
###############################################################################################

# Main files
file = data.frame(read.table(args[2],sep=',',h=T))
annot = vector(length=length(file$Annotation))
annot[file$Annotation=='NMD' | file$Annotation=='NSD'] = 'Other'
annot[is.na(file$Annotation)] = 'Unknown'
annot[file$Annotation=='exonic' | file$Annotation=='CDS' |
        file$Annotation=='intronic' | 
        file$Annotation=='CDS in frame' | 
        file$Annotation=='Antisense'] = 'CDS'
annot[file$Annotation=='5UTR and CDS' | file$Annotation=='5UTR'] = '5UTR'
annot[file$Annotation=='3UTR'] = '3UTR'
annot[file$Annotation=='lncRNA' | file$Annotation=='ncrna' | file$Annotation=='lncrna' | 
        file$Annotation=='ncRNA' | file$Annotation=='Non-coding'] = 'lncRNA'
annot[file$Annotation=='intergenic'] = 'intergenic'
annot[file$Annotation==' ' | file$Annotation=='Not known' ] = 'Unknown'
file$Annotation=as.factor(annot)

file$Subcellular.loc.[file$Subcellular.loc.=='None' | file$Subcellular.loc.==''] = '-'

levels(file$Start.codon) = c(levels(file$Start.codon),'Unknown')
file$Start.codon[file$Start.codon==' ' | file$Start.codon=='None' | is.na(file$Start.codon) ] = 'Unknown'

f <- file[substr(file$sPEP.id,1,3)!='SEP',]
human <- file[file$Organism=='Human',]
cf <- file[file$AA!='',]

# ProtParam file
prot <- as.data.frame(read.table(paste(args[1],"protParams.tsv",sep="/"),sep="\t",h=T))
prot = cbind(prot,ori=file[match(prot$id,file$sPEP.id),]$Origin)

# NetSurfP file
surf <- as.data.frame(read.table(paste(args[1],"NetSurfP.txt",sep='/'),skip=15))
colnames(surf) <- c("class","aa","seqname","aanumber","rsa","asa","rsascore","phelix","pshand","pcoil")

# Anchor file
anchor <- data.frame(read.table(paste(args[1],"Anchor.txt",sep='/')))
colnames(anchor) <- c('number', 'aa', 'anchorp', 'output', 'iupred', 'score', 'S', 'eint', 'egain', 'name')

# Motif file
motifs <- data.frame(read.table(paste(args[1],"Motifs.tsv",sep='/'),sep='\t',h=T))
nmotifs <- sapply(file[file$AA!='',]$sPEP.id, function(name) { return(length(motifs[motifs$sPEP.id==as.character(name),1]))})

# AA properties file
aa <- as.data.frame(read.table(args[3],sep=',',h=T))
factors <- sapply(file$sPEP.id, function(x) {
  as.numeric(substr(x,1,3)=='Sla' | substr(x,1,3)=='SEP' | substr(x,1,7)=='Sendoel')})

###############################################################################################
#     Edit tables and compute summary
###############################################################################################

# Add column of output value
addcol <- function(f) { 
  res <- apply(f[,8:10],1,function(x) {
    if (x[1]>=0.5) return("H")
    else if (x[2]>=0.5) return("S")
    else if (x[3]>=0.5) return("C")
    else return("-")
  })
  return(cbind(f,res))
}
surf <- addcol(surf)

# Calculate structure length
prednum <- function(v) {
  le <- length(v)
  phel <- sum(rle(v)$lengths[rle(v)$values=="H" & rle(v)$lengths>7])/length(v)
  nhel <- length(rle(v)$lengths[rle(v)$values=="H" & rle(v)$lengths>7])
  psh <- sum(rle(v)$lengths[rle(v)$values=="S" & rle(v)$lengths>2])/length(v)
  nsh <- length(rle(v)$lengths[rle(v)$values=="S" & rle(v)$lengths>2])
  pco <- sum(rle(v)$lengths[rle(v)$values=="C"])/length(v)
  nco <- length(rle(v)$lengths[rle(v)$values=="C"])
  return(cbind(le,phel,psh,pco,nhel,nsh,nco))
}

out=data.frame()
for (name in file$sPEP.id) {
  out <- rbind(out,prednum(as.vector(surf$res[surf$seqname==name])))
}
row.names(out)=file$sPEP.id
out <- cbind(out,origin=file$Origin)
out <- na.omit(out)

# Get Anchor disorder binary output
dis <- as.factor(apply(anchor,1, function(x) {if (x[5]>=0.5) {return(1)} else {return(0)}}))
anchor <- cbind(anchor,dis)
mdis <- sapply(unique(anchor$name), function(x) {
  return(mean(anchor[anchor$name==x,5]))
})
# rownames(out[out$origin=='Mass Spec' & out$nhel==0 & out$nsh==0,])

# Checks if motifs are in disordered regions
dismot <- apply(motifs,1,function(x) {
  return(mean(anchor[anchor$name==gsub(':','_',x[1]),][(as.numeric(x[3])+1):x[4],5]))
})
# length(dismot[dismot<0.5])/length(motifs$sPEP.id)
# length(dismot[dismot>=0.5])/length(motifs$sPEP.id)
print ('Motifs not associated to disordered regions (to discard) :')
print(motifs[dismot<0.5,])
print ('Motifs associated to disordered regions with proba > 0.5 :')
print(motifs[dismot>=0.5,])

###############################################################################################
#     Plots
###############################################################################################

# General features

quartz(title="sORFs",bg="white",height = 9,width=11)
p1<-ggplot(data=f) +
  geom_bar(aes(Start.codon, fill=Origin,y=..count../sum(..count..)),position=position_dodge(),alpha=0.7) +
  scale_fill_discrete(name="Origin") +
  labs(title='Start codon',y='Proportion',x='Start codon')
p2<-ggplot(data=file) + 
  geom_histogram(aes(x=Length,y=..density..,fill=Origin),breaks=seq(10, 150, by=5),alpha=0.7) + 
  labs(title="Peptides length",x="Length",y="Density") + 
  scale_fill_discrete(name="Origin")
p3 <- ggplot(data=f) +
  geom_bar(aes(Annotation, fill=factor(Origin),weight=rep(1/length(f$Start.codon),length(f$Start.codon))),position=position_dodge(),alpha=0.7) +
  scale_fill_discrete(name="Origin") +
  labs(title='Location',y='Proportion',x='Annotation')
p4 <- ggplot(data=file) + geom_boxplot(aes(x=Origin,y=Length,fill=Origin))
grid.arrange(p1,p4,p3,ncol=1)
quartz.save(file=paste(args[1],"Rplots/sEP_attr.pdf",sep='/'),type="pdf",device=dev.cur())

# Secondary structure predictions

quartz(title="Structure prediction",bg="white")
p1<-ggplot(out, aes(y = phel, x = le, color=origin)) +
  geom_point() + 
  geom_text_repel(data=out[out$phel>0.75,],aes(label=substr(row.names(out[out$phel>0.75,]),1,6))) +
  labs(title="Helical content vs length",x="Peptide length",y="Helix content") +
  scale_colour_discrete(name="Origin")
p2<-ggplot(out, aes(y = psh, x = le, color=origin)) +
  geom_point() + 
  geom_text_repel(data=out[out$psh>0.65,],aes(label=substr(row.names(out[out$psh>0.65,]),1,6))) +
  labs(title="Beta Sheet content vs length",x="Peptide length",y="Sheet content") +
  scale_colour_discrete(name="Origin")
grid.arrange(p1,p2,ncol=1)
quartz.save(file=paste(args[1],"Rplots/sEP_content.pdf",sep='/'),type="pdf",device=dev.cur())

quartz(title="Structure prediction",bg="white")
p1<-ggplot(data = out, aes(x=rep("hel",length(origin)), y=phel)) + geom_boxplot(aes(fill=origin)) +
  labs(title="Helical content comparison",x="Helix",y="Content") +
  scale_fill_discrete(name="Origin") 
p2<-ggplot(data = out, aes(x=rep("sheet",length(origin)), y=psh)) + geom_boxplot(aes(fill=origin)) +
  labs(title="Sheet content comparison",x="sheets",y="Content") +
  scale_fill_discrete(name="Origin") 
grid.arrange(p1,p2,ncol=2)
quartz.save(file=paste(args[1],"Rplots/Structure_pred.pdf",sep='/'),type="pdf",device=dev.cur())

quartz(title="Sheet content vs helix content",bg="white")
p1 <- ggplot(data = out, aes(x=phel, y=psh)) + geom_point(aes(color=origin),size=2,alpha=0.8) +
  labs(title="Sheet vs helical content",x="Helix",y="Sheet") +
  scale_fill_discrete(name="Origin") 
p2 <- ggplot(data = out, aes(x=nhel, y=nsh)) + geom_point(aes(color=origin),size=2,alpha=0.8) +
  labs(title="Sheet vs helical content",x="Helix",y="Sheet") +
  scale_fill_discrete(name="Origin") 
grid.arrange(p1,ncol=1)
quartz.save(file=paste(args[1],"Rplots/Sheet_vs_helix.pdf",sep='/'),type="pdf",device=dev.cur())

# Amino acid frequencies

Y1 <- apply(prot[prot$ori=='Mass Spec',6:25],2,sum)/sum(prot[prot$ori=='Mass Spec',]$len)*100
Y2 <- apply(prot[prot$ori=='Ribo Prof',6:25],2,sum)/sum(prot[prot$ori=='Ribo Prof',]$len)*100
quartz(title="AA content",bg="white")
ggplot(data=aa,aes(x=frequency)) + 
  # geom_point(aes(y=Y1,color=as.factor('MS'))) + geom_point(aes(y=Y2,color=as.factor('RP'))) +
  geom_label(aes(y=Y1,color=as.factor('MS'),label=names(Y1)),label.size = 0.15,alpha=0.8) + 
  geom_label(aes(y=Y2,color=as.factor('RP'),label=names(Y2)),label.size = 0.15,alpha=0.8) +
  geom_abline(slope = 1, intercept = 0,alpha=0.4) + xlim(0,12) + ylim(0,12) + scale_color_discrete(name='Origin') +
  labs(title='AA frequencies',x='Expected from theorical values',y='Observed in sEPs')
quartz.save(file=paste(args[1],"Rplots/AA_freqs.pdf",sep='/'),type="pdf",device=dev.cur())

# quartz(title="AA content",bg="white")
# ggplot(data=aa) + 
#   geom_label(aes(y=Y1,x=Y2), label=names(Y1)) +
#   geom_abline(slope = 1, intercept = 0,alpha=0.4) + xlim(0,12) + ylim(0,12) +
#   labs(title='AA frequencies',x='Frequency in RP peptides',y='Frequency in MS peptides')
# 
# Y1 <- apply(prot[match(file[nmotifs>0,]$sPEP.id,prot$id),6:25],2,sum)/sum(prot[match(file[nmotifs>0,]$sPEP.id,prot$id),2])*100
# Y2 <- apply(na.omit(prot[match(file[nmotifs<1,]$sPEP.id,prot$id),6:25]),2,sum)/sum(na.omit(prot[match(file[nmotifs<1,]$sPEP.id,prot$id),2]))*100
# quartz(title="AA content",bg="white")
# ggplot(data=aa,aes(x=frequency)) + 
#   geom_label(aes(y=Y1, color=as.factor(1), label=names(Y1)),label.size = 0.15,alpha=0.8) + 
#   geom_label(aes(y=Y2, color=as.factor(2), label=names(Y2)),label.size = 0.15,alpha=0.8) + 
#   geom_abline(slope = 1, intercept = 0,alpha=0.4) + xlim(0,12) + ylim(0,12)  + 
#   scale_color_discrete(name='Number of motifs', labels=c('N >= 1', 'N = 0')) +
#   labs(title='AA frequencies',x='Expected',y='Observed in sEPs')


# Factorial Correspondence Analysis on sEPs amino acid composition

afc <- dudi.coa(prot[,6:25],scannf = F,nf=2)

gravity <- function(vec, fca) {
  X <- vector(length=length(vec))
  Y <- vector(length=length(vec))
  f <- vector(length=length(vec))
  for (i in 1:length(levels(vec))) {
    X[vec==levels(vec)[i]] <- mean(fca$li$Axis1[vec==levels(vec)[i]])
    Y[vec==levels(vec)[i]] <- mean(fca$li$Axis2[vec==levels(vec)[i]])
    f[vec==levels(vec)[i]] <- levels(vec)[i]
  }
  return(data.frame(X,Y,f))
}

vec <- file[match(prot$id,file$sPEP.id),]$Origin
plotfca <- function(fca, vec, cofactor, scaleco, scaleli) {
  Z <- gravity(vec, fca)
  quartz(bg='white',width=13,height=9)
  p1 <- ggplot(data=fca$co,aes(x=Comp1, y=Comp2, color=cofactor)) + geom_point() +
    geom_vline(xintercept=0, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5) +
    geom_label_repel(aes(label=row.names(fca$co))) + xlim(-1,1) +
    scale_color_discrete(name=scaleco) +
    labs(title='Amino acids coordinates on FCA principal axis')
  p2 <- ggplot(data=fca$li) +
    geom_vline(xintercept=0, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5) +
    geom_point(aes(x=rev(Axis1),y=rev(Axis2),color=as.factor(rev(Z$f))),size=0.8) +
    geom_segment(aes(x=rev(Z$X),y=rev(Z$Y),xend=rev(Axis1),yend=rev(Axis2),color=as.factor(rev(Z$f))),alpha=0.3) + 
    labs(title='sEPs coordinates on FCA axis', x='Axis1', y='Axis2') + scale_color_discrete(name=scaleli)
  grid.arrange(p1,p2,ncol=1)
}
plotfca(afc, vec, aa$hydropathy, 'Hydropathy', 'Origin')
quartz.save(file=paste(args[1],"Rplots/FCA1.pdf",sep='/'),type="pdf",device=dev.cur())
plotfca(afc, vec, aa$charge, 'Charge', 'Origin')
quartz.save(file=paste(args[1],"Rplots/FCA2.pdf",sep='/'),type="pdf",device=dev.cur())
plotfca(afc, vec, aa$class, 'Class', 'Origin')
quartz.save(file=paste(args[1],"Rplots/FCA3.pdf",sep='/'),type="pdf",device=dev.cur())

scatterafc <- function(afc, fac, nc, nr) {
  Z <- gravity(fac, afc)
  pl <- lapply(1:length(levels(fac)), function(i) {
    ggplot(data=afc$li[Z$f==levels(fac)[i],]) +
      geom_vline(xintercept=0, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5) +
      geom_point(aes(x=Axis1,y=Axis2,color=as.factor(Z$f[Z$f==levels(fac)[i]]))) +
      geom_segment(aes(x=Z$X[Z$f==levels(fac)[i]],y=Z$Y[Z$f==levels(fac)[i]],
                       xend=Axis1,yend=Axis2,
                       color=as.factor(Z$f[Z$f==levels(fac)[i]])),alpha=0.3) +
      scale_color_discrete(name='Origin',h=c(exp(i),255))
  })
  grid.arrange(rectGrob())
  marrangeGrob(pl,ncol=nc,nrow=nr)
}
scatterafc(afc, file[match(prot$id,file$sPEP.id),]$Annotation, 3, 3)
quartz.save(file=paste(args[1],"Rplots/FCA4.pdf",sep='/'),type="pdf",device=dev.cur())


# Clustering with kmeans

cluster1 <- kmeans(as.matrix(cbind(out[,2:4],prot[,6:25])), 2)
cluster1$size

mat <- matrix(cbind(mdis,as.numeric(cf$Annotation), as.numeric(cf$Subcellular.loc.), 
                    as.numeric(cf$Signal.Peptide)),ncol=4)
cluster2 <- kmeans(mat,2)
comp <- as.numeric(cf$Origin)

quartz(bg='white')
p1 <- ggplot(data=afc$li) +
  geom_vline(xintercept=0, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5) +
  geom_point(aes(x=Axis1,y=Axis2,color=as.factor(2-cluster2$cluster+1))) +
  scale_color_discrete(name='Cluster',labels=c(1,2)) +
  # geom_text_repel(data=afc$li[(3-cluster2$cluster)!=comp,], aes(x=Axis1,y=Axis2,label=cf[(3-cluster2$cluster)!=comp,1])) +
  xlim(-1,1) + labs(title='sEPs coordinates on FCA axis')
p2 <- ggplot(data=afc$li) +
  geom_vline(xintercept=0, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5) +
  geom_point(aes(x=Axis1,y=Axis2,color=as.factor(prot$ori))) +
  xlim(-1,1) + labs(title='sEPs coordinates on FCA axis')
grid.arrange(p1,p2,ncol=1)



# Subcellular localization

# quartz(bg='white')
# ggplot(data=human[human$Origin=='Mass Spec',]) + 
#   geom_bar(aes(Subcellular.loc.,y=..count../sum(..count..))) + ylim(0,0.6) +
#   labs(title='Subcellular localization prediction', y='Proportion of peptides')

quartz(bg='white')
ggplot(data=human) + 
  geom_bar(aes(Subcellular.loc., fill=Origin,y=..count..),position=position_dodge()) +
  labs(title='Subcellular localization prediction', y='Number of peptides')
quartz.save(file=paste(args[1],"Rplots/Subcellular_loc.pdf",sep='/'),type="pdf",device=dev.cur())





