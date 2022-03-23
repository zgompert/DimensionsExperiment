library(data.table)
## L. melissa
G<-fread("lmel_geno.txt",header=FALSE,sep=",")
Gmat<-as.matrix(G)
SNPs<-as.matrix(read.table("lmelScafPos.txt",header=FALSE))

lgScafs<-as.numeric(names(table(SNPs[,1])[which(table(SNPs[,1]) > 100)]))
Nlg<-24
LDmat<-vector("list",Nlg)
Dmat<-vector("list",Nlg)
for(i in 1:Nlg){
	cat("Working on ",i,"\n")
	Gsub<-Gmat[,SNPs[,1]==lgScafs[i]]
	Psub<-SNPs[SNPs[,1]==lgScafs[[i]],2]
	LDmat[[i]]<-cor(Gsub)^2
	Dmat[[i]]<-LDmat[[i]]
	for(k in 1:length(Psub)){
		Dmat[[i]][k,]<-abs(Psub[k]-Psub)
	}
}

## vectorize
LDvec<-vector("list",Nlg)
Dvec<-vector("list",Nlg)
for(i in 1:Nlg){
	LDvec[[i]]<-as.vector(LDmat[[i]][upper.tri(LDmat[[i]])])
	Dvec[[i]]<-as.vector(Dmat[[i]][upper.tri(Dmat[[i]])])
}
LD<-unlist(LDvec)
D<-unlist(Dvec)

bnds<-c(0,100,1000,5000,seq(10000,100000,10000))
LDest<-matrix(NA,nrow=3,ncol=14)
for(i in 1:13){
	xx<-which(D > bnds[i] & D <= bnds[i+1])
	LDest[,i]<-quantile(LD[xx],probs=c(0.5,0.95,0.99))
}
xx<-which(D > bnds[14])
LDest[,14]<-quantile(LD[xx],probs=c(0.5,0.95,0.99))

## Alfalfa
G<-fread("msat_geno.txt",header=FALSE,sep=",")
Gmat<-as.matrix(G)
SNPs<-as.matrix(read.table("msatScafPos.txt",header=FALSE))

lgScafs<-1:8
Nlg<-8
LDmat<-vector("list",Nlg)
Dmat<-vector("list",Nlg)
for(i in 1:Nlg){
	cat("Working on ",i,"\n")
	Gsub<-Gmat[,SNPs[,1]==lgScafs[i]]
	Psub<-SNPs[SNPs[,1]==lgScafs[[i]],2]
	LDmat[[i]]<-cor(Gsub)^2
	Dmat[[i]]<-LDmat[[i]]
	for(k in 1:length(Psub)){
		Dmat[[i]][k,]<-abs(Psub[k]-Psub)
	}
}

## vectorize
LDvec<-vector("list",Nlg)
Dvec<-vector("list",Nlg)
for(i in 1:Nlg){
	LDvec[[i]]<-as.vector(LDmat[[i]][upper.tri(LDmat[[i]])])
	Dvec[[i]]<-as.vector(Dmat[[i]][upper.tri(Dmat[[i]])])
}
LD<-unlist(LDvec)
D<-unlist(Dvec)

bnds<-c(0,100,1000,5000,seq(10000,100000,10000))
LDestMsat<-matrix(NA,nrow=3,ncol=14)
for(i in 1:13){
	xx<-which(D > bnds[i] & D <= bnds[i+1])
	LDestMsat[,i]<-quantile(LD[xx],probs=c(0.5,0.95,0.99))
}
xx<-which(D > bnds[14])
LDestMsat[,14]<-quantile(LD[xx],probs=c(0.5,0.95,0.99))

save(list=ls(),file="LD.rdat")

## plots
library(RColorBrewer)
cs<-rev(brewer.pal(n=4,"YlOrBr")[-1])
cl<-1.38;ca<-1.1;cm=1.38
pdf("sfig_LD.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1))
plot(bnds,LDestMsat[3,],ylim=c(0,1),xlab="Distance (bp)",ylab="LD (r2)",type='b',pch=19,col=cs[1],cex.lab=cl,cex.axis=ca)
points(bnds,LDestMsat[2,],type='b',pch=19,col=cs[2])
points(bnds,LDestMsat[1,],type='b',pch=19,col=cs[3])
legend(5e04,1,c("99th percentile","95th percentile","median"),fill=cs,bty='n')
title(main="(A) LD in M. sativa",cex.main=cm)
plot(bnds,LDest[3,],ylim=c(0,1),xlab="Distance (bp)",ylab="LD (r2)",type='b',pch=19,col=cs[1],cex.lab=cl,cex.axis=ca)
points(bnds,LDest[2,],type='b',pch=19,col=cs[2])
points(bnds,LDest[1,],type='b',pch=19,col=cs[3])
title(main="(B) LD in L. melissa",cex.main=cm)
dev.off()

round(LDest,3)
#      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#[1,] 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002
#[2,] 0.052 0.059 0.048 0.045 0.045 0.043 0.044 0.043 0.044 0.044 0.043 0.044
#[3,] 0.625 0.349 0.207 0.188 0.174 0.164 0.164 0.163 0.165 0.158 0.156 0.157
#     [,13] [,14]
#[1,] 0.002 0.002
#[2,] 0.045 0.041
#[3,] 0.159 0.141
round(LDestMsat,3)
#      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#[1,] 0.050 0.017 0.014 0.012 0.011 0.011 0.010 0.010 0.010 0.009 0.009 0.009
#[2,] 0.862 0.193 0.157 0.142 0.134 0.131 0.129 0.126 0.125 0.125 0.124 0.123
#[3,] 0.992 0.358 0.287 0.244 0.220 0.213 0.204 0.201 0.197 0.196 0.196 0.192
#     [,13] [,14]
#[1,] 0.009 0.005
#[2,] 0.122 0.109
#[3,] 0.192 0.167

