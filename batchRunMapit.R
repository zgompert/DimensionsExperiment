## get command line arguments
argv<-as.numeric(commandArgs(trailingOnly=TRUE))
## first is phenotype number
phn<-argv[1]
## model, 1 = lmel, 2 = msat, 3 = lmel x msat
mod<-argv[2]

### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)
library(scales)
library(data.table)

### Load in functions to make QQ-plot plots ###
source("MAPIT/QQPlot.R")

#NOTE: This code assumes that the basic C++ functions are set up on the computer in use. If not, the MAPIT functions and Rcpp packages will not work properly. Mac users please refer to the homebrew applications and install the gcc commands listed in the README.md file before running the rest of the code [Warning: This step may take about an hour...].

### Load in the C++ MAPIT wrapper functions ###
source("MAPIT/Standard Version/MAPIT.R"); sourceCpp("MAPIT/Standard Version/MAPIT.cpp")

G<-fread("../Gemma/comb_geno",header=FALSE)
Gmat<-as.matrix(G[,-c(1:3)])
msat<-1:161008## snp numbers
lmel<-161009:224202


maf<-apply(Gmat,1,mean)/2
z<-maf>0.5
maf[z]<-1-maf[z]
common<-which(maf > 0.05)
msat<-msat[msat %in% common]
lmel<-lmel[lmel %in% common]


## standardize genotypes
X<-Gmat
Xmn<-apply(X,1,mean,na.rm=TRUE)
Xsd<-apply(X,1,sd,na.rm=TRUE)
L<-dim(Gmat)[1]
for(i in 1:L){
	X[i,]<-(X[i,]-Xmn[i])/Xsd[i]
}

## traits
#ph<-read.table("../Gemma/subrandom_pheno_residTraits.txt",header=FALSE)
ph<-read.table("../Gemma/subgemma_pheno_residTraits.txt",header=FALSE)

## grab relevant phenotype
Y<-ph[,phn]
mn<-mean(Y,na.rm=TRUE)
sig<-sd(Y,na.rm=TRUE)
Y<-(Y-mn)/sig
drop<-which(is.na(Y)==TRUE)

if(length(drop)>0){
	Gx<-X[,-drop]
	Y<-Y[-drop]
}
if(length(drop)==0){
	Gx<-X
	Y<-Y
}
## run model
if(mod==1){## lmel
	mapit_out<-MAPIT(Gx[lmel,],as.matrix(Y))
}
if(mod==2){## msativa
	mapit_out<-MAPIT(Gx[msat,],as.matrix(Y))
}
if(mod==3){## all
	mapit_out<-MAPIT(Gx,as.matrix(Y))
}
#o<-paste("mapitOutRan_ph",phn,"_mod",mod,".rdat",sep="")
o<-paste("mapitOutMaf_ph",phn,"_mod",mod,".rdat",sep="")
save(list=ls(),file=o)
