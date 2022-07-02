## makes pairwise epistasis genotype files
#For within species,  top 150 marginal epistasis, all pairwise = ~ 11k additional terms
library(data.table)


## L. meliss
G<-fread("lmel_geno",header=FALSE)
epifiles<-list.files(pattern="^lmel_epi",path="../mappit",full.names=TRUE)

for(k in 1:3){
	ef<-read.table(epifiles[k],header=FALSE)
	cc<-ef[,1]
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
	G<-as.data.frame(G)
	colnames(G)<-1:860
	colnames(Gdf)<-1:860

	Gout<-rbind(G,Gdf)
	write.table(file=paste("epix_lmel_geno_",k,sep=""),Gout,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

## M. sativa
G<-fread("msat_geno",header=FALSE)
epifiles<-list.files(pattern="^msat_epi",path="../mappit",full.names=TRUE)

for(k in 1:3){
	ef<-read.table(epifiles[k],header=FALSE)
	cc<-ef[,1]
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
	Nj<-dim(G)[2]
	G<-as.data.frame(G)
	colnames(G)<-1:Nj
	colnames(Gdf)<-1:Nj

	Gout<-rbind(G,Gdf)
	write.table(file=paste("epix_msat_geno_",k,sep=""),Gout,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
