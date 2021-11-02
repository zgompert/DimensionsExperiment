## misc analyses looking at generalization and prediction 
library(RColorBrewer)
library(scales)
library(data.table)


### L melissa ##

## compute polygenic scores
G<-fread("../Entropy/lmel_geno.txt",header=FALSE)
G<-as.matrix(G)
N<-dim(G)[1]
L<-dim(G)[2]
## center G
for(i in 1:L){
        G[,i]<-G[,i]-mean(G[,i])
}

files_cat<-list.files(path="files_mav",pattern="mav_o_lmel_fit_gemmalmel_pheno_residTraits_ph",full.names=TRUE)

Np<-9
psCat<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
        cat(j,"\n")
        ph<-fread(files_cat[j],header=FALSE)
        ph<-unlist(as.vector(ph[,2]))# get trait
        psCat[,j]<-G %*% ph
}


ids<-read.table("../Entropy/CleanLmelIds.txt",header=FALSE)


## phenotypic data
dat<-read.table("../Pheno/dimensions_lmel_dat.csv",header=TRUE,sep=",")
matchIds<-rep(NA,N)
for(i in 1:N){
	a<-which(as.character(dat$Cat_ID) == as.character(ids[i,1]))
	if(length(a)==1){
		matchIds[i]<-a
	}
}
sortDat<-dat[matchIds,]
mean(as.character(sortDat$Cat_ID)==as.character(ids[,1]),na.rm=TRUE)
# 1 
# that worked

sg<-grep(pattern="SG",sortDat$Plant_ID)
psSG<-psCat[sg,]
sdatSG<-sortDat[sg,]
plt<-as.factor(unlist(strsplit(x=as.character(sdatSG$Plant_ID),split="-"))[seq(2,468,3)])
residSdatSG<-matrix(NA,nrow=156,ncol=9)
for(i in 1:9){
	xx<-which(is.na(sdatSG[,i+5])==FALSE)
	o<-lm(sdatSG[,i+5] ~ sdatSG$hatchD + plt)
	residSdatSG[xx,i]<-o$residuals
}

catSgCors<-matrix(NA,nrow=9,ncol=4)
for(i in 1:9){
	o<-cor.test(psSG[,i],residSdatSG[,i])
	catSgCors[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
}

save(list=ls(),file="sgPredict.rdat")

### Msativa ##

## compute polygenic scores
G<-fread("../Entropy/msat_geno.txt",header=FALSE)
G<-as.matrix(G)
N<-dim(G)[1]
L<-dim(G)[2]
## center G
for(i in 1:L){
        G[,i]<-G[,i]-mean(G[,i])
}

files_cat<-list.files(path="files_mav",pattern="mav_o_msat_fit_gemma_pheno_residTraits_ph",full.names=TRUE)

Np<-9
psCat<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
        cat(j,"\n")
        ph<-fread(files_cat[j],header=FALSE)
        ph<-unlist(as.vector(ph[,2]))# get trait
        psCat[,j]<-G %*% ph
}

ids<-read.table("../Entropy/MsativaIds.txt",header=FALSE)


## phenotypic data
dat<-read.table("../Pheno/dimensions_lmel_dat.csv",header=TRUE,sep=",")
dat$Plant_ID<-sub(pattern="-$",replacement="",x=dat$Plant_ID)
matchIds<-rep(NA,N)
for(i in 1:N){
	a<-which(as.character(dat$Plant_ID) == as.character(ids[i,1]))
	if(length(a)==1){
		matchIds[i]<-a
	}
}
sortDat<-dat[matchIds,]
mean(as.character(sortDat$Plant_ID)==as.character(ids[,1]),na.rm=TRUE)
# 1 
# that worked

sg<-grep(pattern="SG",sortDat$Plant_ID)
psSG_plant<-psCat[sg,]
sdatSG_plant<-sortDat[sg,]
plt<-as.factor(unlist(strsplit(x=as.character(sdatSG_plant$Plant_ID),split="-"))[seq(2,468,3)])
residSdatSG_plant<-matrix(NA,nrow=156,ncol=9)
for(i in 1:9){
	xx<-which(is.na(sdatSG_plant[,i+5])==FALSE)
	o<-lm(sdatSG_plant[,i+5] ~ sdatSG_plant$hatchD + plt)
	residSdatSG_plant[xx,i]<-o$residuals
}

plantSgCors<-matrix(NA,nrow=9,ncol=4)
for(i in 1:9){
	o<-cor.test(psSG_plant[,i],residSdatSG_plant[,i])
	plantSgCors[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
}

save(list=ls(),file="sgPredict.rdat")
## difference between caterpillar and plant is striking

### Combined ##

## compute polygenic scores
G<-fread("comb_sg_geno.txt",header=FALSE)
G<-as.matrix(G)
N<-dim(G)[1]
L<-dim(G)[2]
## center G
for(i in 1:L){
        G[,i]<-G[,i]-mean(G[,i])
}

files_cat<-list.files(path="files_mav",pattern="mav_o_comb_fit_subgemma_pheno_residTraits_ph",full.names=TRUE)

Np<-9
psComb<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
        cat(j,"\n")
        ph<-fread(files_cat[j],header=FALSE)
        ph<-unlist(as.vector(ph[,2]))# get trait
        psComb[,j]<-G %*% ph
}


## phenotypic data, already sorted
dat<-read.table("sg_pheno.txt",header=TRUE)

plt<-as.factor(unlist(strsplit(x=as.character(dat$Plant_ID),split="-"))[seq(2,408,3)])
residSdatSG_comb<-matrix(NA,nrow=136,ncol=9)
for(i in 1:9){
	xx<-which(is.na(dat[,i+5])==FALSE)
	o<-lm(dat[,i+5] ~ dat$hatchD + plt)
	residSdatSG_comb[xx,i]<-o$residuals
}

combSgCors<-matrix(NA,nrow=9,ncol=4)
for(i in 1:9){
	o<-cor.test(psComb[,i],residSdatSG_comb[,i])
	combSgCors[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
}

save(list=ls(),file="sgPredict.rdat")
## difference between caterpillar and plant is striking
