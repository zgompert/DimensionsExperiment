library(data.table)

## read in genotype matrix for M. sativa
Gmsat<-fread("msat_geno.txt",sep=",",header=FALSE)
ms_ids<-read.table("MsativaIds.txt",header=FALSE)

## read in genotype matrix for L. melissa
Glmel<-fread("lmel_geno.txt",sep=",",header=FALSE)
lm_ids<-read.table("LmelIds.txt",header=FALSE)

## fst (Nei Gst)
fst<-function(P=NA){
	## P = allele freq. matrix, cols=pops, rows=SNPs
	H<-P * (1-P) * 2
	Hw<-apply(H,1,mean)## mean across pops for each locus
	pbar<-apply(P,1,mean)
	Ht<-pbar * (1-pbar) * 2
	Fst<-1 - mean(Hw)/mean(Ht)
	return(Fst)
}

GlmelM<-as.matrix(Glmel)
lm_pops<-read.table("LmelPops.txt",header=FALSE)
Np<-6
Nl<-dim(Glmel)[2]
Plmel<-matrix(NA,nrow=Nl,ncol=Np)
for(i in 1:Nl){
	Plmel[i,]<-tapply(X=GlmelM[,i],INDEX=lm_pops[,1],mean)/2
}

fst(Plmel)
#[1] 0.02944834

GmsatM<-as.matrix(Gmsat)
ms_pops<-read.table("MsativaPops.txt",header=FALSE)
drop<-which(is.na(ms_pops))
ms_pops<-ms_pops[-drop,]
GmsatM<-GmsatM[-drop,]

Np<-11
Nl<-dim(Gmsat)[2]
Pmsat<-matrix(NA,nrow=Nl,ncol=Np)
for(i in 1:Nl){
	Pmsat[i,]<-tapply(X=GmsatM[,i],INDEX=ms_pops,mean)/2
}
fst(Pmsat)
#[1] 0.04454983

save(list=ls(),file="fst.rdat")


## PCA, center 
pcmsat<-prcomp(Gmsat,center=TRUE,scale=FALSE)
pclmel<-prcomp(Glmel,center=TRUE,scale=FALSE)

cm<-1.5;cl<-1.5;ca<-1.1

css<-read.table("colors.txt",header=FALSE,comment.char="%")
ms_css<-css[7:17,]
ms_index<-rep(NA,dim(Gmsat)[1])
for(i in 1:11){
	pop<-strsplit(as.character(ms_css[i,1]),"-")[[1]][2]
	ms_index[grep(pop,ms_ids[,1])]<-i
}
ms_index[c(22,1141)]<-2 ## was miss-spelled

lm_css<-css[1:6,]
lm_index<-rep(NA,dim(Glmel)[1])
for(i in 1:6){
	pop<-strsplit(as.character(lm_css[i,1]),"-")[[1]][2]
	lm_index[grep(pop,lm_ids[,1])]<-i
}
ms_index[c(22,1141)]<-2 ## was miss-spelled


xx<-sample(1:dim(Gmsat)[1],dim(Gmsat)[1],replace=FALSE)
yy<-sample(1:dim(Glmel)[1],dim(Glmel)[1],replace=FALSE)

pdf("fig_MapPca.pdf",width=10,height=10)
layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=FALSE),widths=c(5,5),heights=c(5,5))
par(mar=c(4,5,2.5,1.5))

plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
title(main="(a) Map of sources",cex.main=cm)

plot(pcmsat$x[xx,1],pcmsat$x[xx,2],pch=19,col=as.character(ms_css[,2])[ms_index][xx],xlab="PC1 (4.5%)",ylab="PC2 (1.8%)",cex.lab=cl,cex.axis=ca)
title(main="(b) PCA of M. sativa",cex.main=cm)

plot(pclmel$x[yy,1],pclmel$x[yy,2],pch=19,col=as.character(lm_css[,2])[lm_index][yy],xlab="PC1 (7.8%)",ylab="PC2 (3.0%)",cex.lab=cl,cex.axis=ca)
title(main="(c) PCA of L. melissa",cex.main=cm)
dev.off()

