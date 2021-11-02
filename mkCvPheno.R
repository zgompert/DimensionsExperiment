## make 10-fold CV input files for caterpillar performance

## plant genetcs
ph<-read.table("gemma_pheno_residTraits.txt",header=FALSE)
cvg<-sample(1:10,1055,replace=TRUE)
phCv<-matrix(NA,nrow=1055,ncol=90)
for(i in 1:9){for(j in 1:10){
	x<-j + (i-1)*10
	phCv[,x]<-ph[,i]
	phCv[which(cvg==j),x]<-NA
}}
write.table(x=phCv,file="cv_gemma_pheno_residTraits.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

## caterpillar genetics
ph<-read.table("gemmalmel_pheno_residTraits.txt",header=FALSE)
cvg<-sample(1:10,857,replace=TRUE)
phCv<-matrix(NA,nrow=857,ncol=90)
for(i in 1:9){for(j in 1:10){
	x<-j + (i-1)*10
	phCv[,x]<-ph[,i]
	phCv[which(cvg==j),x]<-NA
}}
write.table(x=phCv,file="cv_gemmalmel_pheno_residTraits.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

## combined
ph<-read.table("subgemma_pheno_residTraits.txt",header=FALSE)
cvg<-sample(1:10,849,replace=TRUE)
phCv<-matrix(NA,nrow=849,ncol=90)
for(i in 1:9){for(j in 1:10){
	x<-j + (i-1)*10
	phCv[,x]<-ph[,i]
	phCv[which(cvg==j),x]<-NA
}}
write.table(x=phCv,file="cv_subgemma_pheno_residTraits.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
