## comput polygenic scores and genetic correlations between plant traits and performance traits
library(data.table)
G<-fread("../../Entropy/msat_geno.txt",header=FALSE)
G<-as.matrix(G)
N<-dim(G)[1]
L<-dim(G)[2]
## center G
for(i in 1:L){
	G[,i]<-G[,i]-mean(G[,i])
}


## chemistry
Np<-1750 ## note in alpha order not strict numeric
psChem<-matrix(NA,nrow=N,ncol=Np)
files_chem<-list.files(pattern="residChem_ph")
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_chem[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psChem[,j]<-G %*% ph
}

## plant traits
files_plant<-list.files(pattern="emma_pheno_residPlant_ph")
Np<-10 ## note in alpha order not strict numeric
psPlant<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_plant[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psPlant[,j]<-G %*% ph
}

## cat traits
files_cat<-list.files(pattern="_gemma_pheno_residTraits_ph")
Np<-9
psCat<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_cat[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psCat[,j]<-G %*% ph
}

## ps correlations, fp only
ids<-read.table("../../Entropy/MsativaIds.txt")
fp<-grep(pattern="FP",ids[,1])
chem_num<-scan("chemPhNum.txt") ## ph chem order
psCors<-matrix(NA,nrow=9,ncol=10+1750)
combPs<-cbind(psChem[fp,order(chem_num)],psPlant[fp,c(1,3:10,2)])
for(i in 1:9){for(j in 1:1760){
	cat(i,"\n")
	psCors[i,j]<-cor(psCat[fp,i],combPs[,j])
}}

## read in h2
chemH2<-read.table("../files_post_chem/comb_ord_post_chem.txt",header=FALSE)
plantH2<-read.table("../output/post_gemma_resid_plant.txt",header=FALSE)
h2<-c(chemH2[,1],plantH2[,1])
catH2<-read.table("../output/post_gemma_resid_traits.txt",header=FALSE)
cTraits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")

save(list=ls(),file="polyScores.rdat")

## plot
pdf("MsatGenCors.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(4,5,2.5,1.5))
for(i in 1:9){
	plot(h2,psCors[i,],pch=19,col=c(rep(alpha("black",.3),1750),rep(alpha("forestgreen",.9),10)),xlab="Trait PVE",ylab="Genetic correlation",cex.lab=1.5)
	title(main=paste(cTraits[i],", PVE =",round(catH2[i,1],3)),cex.main=1.5)
}
dev.off()


