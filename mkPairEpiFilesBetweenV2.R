## makes pairwise epistasis genotype files
#For between species, 75 per species, all pairwise 
library(data.table)


G<-fread("comb_geno",header=FALSE)
epifiles_lm<-list.files(pattern="^comb_lmel_epi",path="../mappit",full.names=TRUE)
epifiles_ms<-list.files(pattern="^comb_msat_epi",path="../mappit",full.names=TRUE)

msat<-1:161008## snp numbers
lmel<-161009:224202


for(k in 1:3){
	ef_lm<-read.table(epifiles_lm[k],header=FALSE)
	ef_ms<-read.table(epifiles_ms[k],header=FALSE)
	cc<-unique(c(lmel[ef_lm[,1]],msat[ef_ms[,1]]))
	Gcc<-as.matrix(G[cc,-c(1:3)])
	## center and standardize
	for(j in 1:dim(Gcc)[1]){
		Gcc[j,]<-(Gcc[j,]-mean(Gcc[j,]))/sd(Gcc[j,])
	}
	Nl<-dim(Gcc)[1]
	Gepi<-matrix(NA,nrow=(Nl * (Nl-1))/2,ncol=dim(Gcc)[2])
	n<-1
	for(i in 1:(Nl-1)){
		for(j in (i+1):Nl){
			Gepi[n,]<-Gcc[i,]*Gcc[j,]
			n<-n+1
		}
	}
	hdr<-cbind(paste("epi",1:(n-1),sep=""),rep("XA",n-1),rep("Xa",n-1))
	Gdf<-data.frame(hdr,round(Gepi,6))
	NN<-dim(G)[2]
	G<-as.data.frame(G)
	colnames(G)<-1:NN
	colnames(Gdf)<-1:NN

	Gout<-rbind(G,Gdf)
	write.table(file=paste("epix_comb_geno_",k,sep=""),Gout,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

