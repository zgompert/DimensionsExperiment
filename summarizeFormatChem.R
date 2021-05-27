library(ape) ## for pcoa

## summarize chemistry data and format in preparation for mapping

cdat<-read.table("MESA_UT_ZG_2019_normalized_NAs.csv",header=TRUE,sep=",")

## Assuming column 7 (LPC.16.0.0.0) is the first compound and NA = 0 (not missing), but will check with Casey
## actually instead replacing with mean, seems like that is what Casey did
fp<-which(cdat$site=="FP") ## field plot
fpids<-as.character(cdat$Envelope[fp])

chem<-as.matrix(cdat[fp,-c(1:6)])
N<-dim(chem)[2]
for(i in 1:N){
	chem[is.na(chem[,i])==TRUE,i]<-0.0004205657## from casey "I used 0.0004205657 to replace NAs, which is 2500 (half of the minimum peak area I used in processing) divided by the mean ISD peak area (5944374)"
}

## log2
l2chem<-log2(chem)

## pca
pc<-prcomp(l2chem,center=TRUE,scale=FALSE)

## compare to spatial configuration in the field plot common garden
fmap<-read.table("fpMap.txt",header=TRUE)

## pcoa on truncated (at 1) spatial distance matrix (using truncation point that ensures all nearest neighbors connected
d<-dist(fmap[,1:2],upper=TRUE,diag=TRUE,method="euclidean")
dt<-as.matrix(d)
dt[dt>1]<-4 ## repalce d > t = 1 with 4
diag(dt)<-4
pcoad<-pcoa(as.dist(dt))

## sort chem
chemPCA<-matrix(NA,nrow=1080,ncol=6)
for(i in 1:1080){
	a<-which(fpids==as.character(fmap$ind_plant)[i])
	if(length(a)==1){
		chemPCA[i,]<-pc$x[a,1:6]
	}
}
chemSort<-matrix(NA,nrow=1080,ncol=1750)
for(i in 1:1080){
	a<-which(fpids==as.character(fmap$ind_plant)[i])
	if(length(a)==1){
		chemSort[i,]<-chem[a,]
	}
}

## playing with choosing MEMs for PC1 of chem
sum(pcoad$values[,1]>0)
# 758 with positive eigenvalujes
r2<-rep(NA,758)
for(i in 1:758){
	o<-lm(chemPCA[,1] ~ pcoad$vectors[,i])
	oo<-summary(o)
	r2[i]<-oo$adj.r.squared
}
newp<-which.max(r2)
o<-lm(chemPCA[,1] ~ pcoad$vectors[,newp])

## chose variable with highest r2 for model, stop when exceeds global model or P > 0.05
## need to write function to do this in a nice way, keep adding variables
fregMem<-function(Y=NA,X=NA,alpha=0.05){
	N<-dim(X)[2]
	j<-0
	## first fit global
	o<-lm(Y~X)
        oo<-summary(o)
	pv<-pf(oo$fstatistic[1],oo$fstatistic[2],oo$fstatistic[3],lower.tail=FALSE)
	globR2<-oo$adj.r.squared
	cvars<-NULL## indexes of cvars
	if(pv <= alpha){ ## global sign
		done<-0
		r2last<-0
		while(done==0){
			r2<-rep(0,N)
			for(k in 1:N){
				if((k %in% cvars)==FALSE){
				        o<-lm(Y~X[,c(cvars,k)])
					oo<-summary(o)
	       	 			r2[k]<-oo$adj.r.squared
				}
			}
			newp<-which.max(r2)
			o<-lm(Y~X[,c(cvars,newp)])
	     	        oo<-summary(o)
			nk<-dim(oo$coefficients)[1]
			pv<-oo$coefficients[nk,4]
			#print(oo$coefficients)
			if(pv < alpha & r2[newp] < globR2 & r2[newp] > r2last & j<200){ ## add and continue
				cvars<-c(cvars,newp)
				j<-j+1
				r2last<-r2[newp]
				cat(length(cvars)," ",round(r2last,3),"\n")
			}
			else{
				done<-1
			}
		}
	}
	if(is.null(cvars)==FALSE){
		o<-lm(Y~X[,cvars])
		oo<-summary(o)
		yy<-which(is.na(Y)==FALSE)
		Y[yy]<-o$residuals
	}
	txt<-paste("done r2 = ",oo$adj.r.squared,", reatined ",length(cvars)," covariates\n",sep="")
	cat(txt)
	return(list(Yresid=Y,Vars=cvars))
}

out<-fregMem(Y=chemPCA[,1],X=pcoad$vectors[,1:758])
## this seems to work
## now put it in a loop, use parallel processing to make this faster
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores-4) #not to overload your computer
registerDoParallel(cl)

finalMatrix<-foreach(i=1:1750, .combine=cbind) %dopar%{
	out<-fregMem(Y=chemSort[,i],X=pcoad$vectors[,1:758])
	oy<-out$Yresid
	oy
}
stopCluster(cl)

## another run to save specific variables used
cl <- makeCluster(cores-4) #not to overload your computer
registerDoParallel(cl)

varMatrix<-foreach(i=1:1750, .combine=list) %dopar%{
	out<-fregMem(Y=chemSort[,i],X=pcoad$vectors[,1:758])
	oy<-out$Vars
	oy
}
stopCluster(cl)

## and now PCs only
cl <- makeCluster(cores-4) #not to overload your computer
registerDoParallel(cl)

pcMatrix<-foreach(i=1:6, .combine=cbind) %dopar%{ ## just 1st 6
	out<-fregMem(Y=chemPCA[,i],X=pcoad$vectors[,1:758])
	oy<-out$Yresid
	oy
}
stopCluster(cl)

## now standardize and write
rChem<-matrix(NA,nrow=1080,ncol=1750)
for(i in 1:1750){
	rChem[,i]<-scale(finalMatrix[,i])
}
write.table(rChem,"pheno_residChem.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

sChem<-matrix(NA,nrow=1080,ncol=1750)
for(i in 1:1750){
	sChem[,i]<-scale(chemSort[,i])
}
write.table(sChem,"pheno_rawChem.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

pcChem<-matrix(NA,nrow=1080,ncol=6)
for(i in 1:6){
	pcChem[,i]<-scale(chemPCA[,i])
}
write.table(pcChem,"pheno_rawPcChem.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

rpcChem<-matrix(NA,nrow=1080,ncol=6)
for(i in 1:6){
	rpcChem[,i]<-scale(pcMatrix[,i])
}
write.table(rpcChem,"pheno_residPcChem.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
save(list=ls(),file="chemFormat.rdat")


