library(ape) ## for pcoa

## summarize caterpillar data and format in preparation for mapping

dat<-read.table("dimensions_lmel_dat.csv",header=TRUE,sep=",")

fp<-grep(pattern="FP",dat$Plant_ID) ## field plot
fpids<-as.character(dat$Plant_ID[fp])

traits<-as.matrix(dat[fp,-c(1:5)])
covar<-dat[fp,1:5]


## compare to spatial configuration in the field plot common garden
fmap<-read.table("fpMap.txt",header=TRUE)

## pcoa on truncated (at 1) spatial distance matrix (using truncation point that ensures all nearest neighbors connected
d<-dist(fmap[,1:2],upper=TRUE,diag=TRUE,method="euclidean")
dt<-as.matrix(d)
dt[dt>1]<-4 ## repalce d > t = 1 with 4
diag(dt)<-4
pcoad<-pcoa(as.dist(dt))

## trait order is
## [1] "w8day"      "w14day"     "wPupa"      "s8day"      "s14day"    
## [6] "sPupa"      "sEclose"    "stime"      "stimeTrunc"
## truncated includes eclosed at 1 day past max time

## sort traits by plant IDs 
traitSort<-matrix(NA,nrow=1080,ncol=9)
covarSort<-covar
for(i in 1:1080){
	a<-which(fpids==as.character(fmap$ind_plant)[i])
	if(length(a)==1){
		traitSort[i,]<-traits[a,]
		covarSort[i,]<-covar[a,]
	}
}

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

## first remove hatch and popualtion effects, then spatial
finalMatrix<-traitSort
traitSortCatResid<-finalMatrix
for(i in 1:9){
	xx<-which(is.na(traitSort[,i])==FALSE)
	o<-lm(traitSort[,i] ~ covarSort$hatchD+covarSort$pop)
	traitSortCatResid[xx,i]<-o$residuals
	d<-dist(fmap[xx,1:2],upper=TRUE,diag=TRUE,method="euclidean")
	dt<-as.matrix(d)
	dt[dt>1]<-4 ## repalce d > t = 1 with 4
	diag(dt)<-4
	pcoad<-pcoa(as.dist(dt))
	eigen<-1:(min(which(pcoad$values[,1] < 0))-2)
	out<-fregMem(Y=traitSortCatResid[xx,i],X=pcoad$vectors[,eigen])
	oy<-out$Yresid
	finalMatrix[xx,i]<-oy
}


## now standardize and write
rTraits<-matrix(NA,nrow=1080,ncol=9)
for(i in 1:9){
	rTraits[,i]<-scale(finalMatrix[,i])
}
write.table(rTraits,"pheno_residTraits.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

sTraits<-matrix(NA,nrow=1080,ncol=9)
for(i in 1:9){
	sTraits[,i]<-scale(traitSort[,i])
}
write.table(sTraits,"pheno_rawTraits.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
save(list=ls(),file="TraitsFormat.rdat")


