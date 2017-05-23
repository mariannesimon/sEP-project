rm(list=ls())
setwd("~/Documents/Slavoff_sPEP")
f = read.table("NetSurfP.out",skip=15,sep=" ",dec=".")
colnames(f) = c("class","aa","seqname","aanumber","rsa",
                "asa","rsascore","phelix","pstrand","pcoil")
head(f,100)
par(mfrow=c(4,4))
for (i in 83:90) {
  if (i!=51 && i != 54 && i!=89) { 
    name = paste("Slavoff_",i,sep="")
    plot(f[f$seqname==name,]$aanumber,f[f$seqname==name,]$phelix,col="red",
         type="l",ylim=c(0,1),main=name,
         xlab = "Amino acid",
         ylab="Probability"
         )
    lines(f[f$seqname==name,]$aanumber,f[f$seqname==name,]$pstrand,col="blue",type="l")
    lines(f[f$seqname==name,]$aanumber,f[f$seqname==name,]$pcoil,col="green",type="l")
  }
}

