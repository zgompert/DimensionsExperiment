library(ape) ## for pcoa

## summarize basic plant trait data and format in preparation for mapping

## functions

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

## compare to spatial configuration in the field plot common garden
fmap<-read.table("fpMap.txt",header=TRUE)

## pcoa on truncated (at 1) spatial distance matrix (using truncation point that ensures all nearest neighbors connected
d<-dist(fmap[,1:2],upper=TRUE,diag=TRUE,method="euclidean")
dt<-as.matrix(d)
dt[dt>1]<-4 ## repalce d > t = 1 with 4
diag(dt)<-4
pcoad<-pcoa(as.dist(dt))

################### leaf toughness #############################
dat<-read.table("may19-toughness.csv",sep=",",header=TRUE)
fpdat<-dat[grep(pattern="FP",dat$ind_plant),]

## sort 
toughSort<-matrix(NA,nrow=1080,ncol=1)
for(i in 1:1080){
	a<-which(fpdat[,1]==as.character(fmap$ind_plant)[i])
	if(length(a)==1){
		toughSort[i,]<-fpdat[a,2]
	}
}

out<-fregMem(Y=toughSort,X=pcoad$vectors[,1:758])
#done r2 = 0.410991775945724, reatined 49 covariates
toughOut<-scale(out$Yresid)

################## plant height ####################################
dat<-read.table("may19-height.csv",sep=",",header=TRUE)
fpdat<-dat[grep(pattern="FP",dat$ind_plant),]

## sort 
heightSort<-matrix(NA,nrow=1080,ncol=1)
for(i in 1:1080){
	a<-which(fpdat[,1]==as.character(fmap$ind_plant)[i])
	if(length(a)==1){
		heightSort[i,]<-fpdat[a,2]
	}
}

out<-fregMem(Y=heightSort,X=pcoad$vectors[,1:758])
#done r2 = 0.226589018826811, reatined 46 covariates
heightOut<-scale(out$Yresid)

################# herbivory ###########################
dat<-read.table("may19-herbivory.csv",sep=",",header=TRUE)
fpdat<-dat[grep(pattern="FP",dat$ind_plant),]

## sort 
herbSort<-matrix(NA,nrow=1080,ncol=1)
for(i in 1:1080){
	a<-which(fpdat[,1]==as.character(fmap$ind_plant)[i])
	if(length(a)==1){
		herbSort[i,]<-sum(fpdat[a,2:3])
	}
}

out<-fregMem(Y=herbSort,X=pcoad$vectors[,1:758])
#done r2 = 0.353732487872706, reatined 62 covariates
herbOut<-scale(out$Yresid)

################### leaf dimensions, trichomes, weight ####
dat<-read.table("may19-leafarea-trichomes-dryweight.csv",sep=",",header=TRUE)
stands<-read.table("stands.txt",header=FALSE)
## turn trichome number to density
for(k in 1:8){
	a<-which(dat$measurer==stands[k,1])
	dat$avg_trichomes[a]<-dat$avg_trichomes[a]/stands[k,2]
}
larea<-dat$avg_leaf_length*dat$avg_leaf_width
lshape<-dat$avg_leaf_length/dat$avg_leaf_width
dat<-cbind(dat,avg_leaf_area=larea,avg_leaf_shape=lshape)

fpdat<-dat[grep(pattern="FP",dat$ind_plant),]

## sort 
leafSort<-matrix(NA,nrow=1080,ncol=7)
for(i in 1:1080){
        a<-which(fpdat[,1]==as.character(fmap$ind_plant)[i])
        if(length(a)==1){
                leafSort[i,]<-as.numeric(fpdat[a,3:9])
        }
}

leafOut<-leafSort
for(i in 1:7){
	out<-fregMem(Y=leafSort[,i],X=pcoad$vectors[,1:758])
	leafOut[,i]<-scale(out$Yresid)
}
# length, width,  weight, trichomes, SLA, area, shape 
#done r2 = 0.187210320669472, reatined 22 covariates
#done r2 = 0.340808425965393, reatined 54 covariates
#done r2 = 0.507521945715951, reatined 72 covariates
#done r2 = 0.452947358826808, reatined 57 covariates
#done r2 = 0.43690804068477, reatined 77 covariates
#done r2 = 0.365984743962929, reatined 56 covariates
#done r2 = 0.179720157726493, reatined 31 covariates

rawTrt<-cbind(leafSort,heightSort,toughSort,herbSort)
residTrt<-cbind(leafOut,heightOut,toughOut,herbOut)
for(i in 1:10){
	rawTrt[,i]<-scale(rawTrt[,i])
	residTrt[,i]<-scale(residTrt[,i])
}

write.table(residTrt,file="pheno_residPlant.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(rawTrt,file="pheno_rawPlant.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

save(list=ls(),file="plantFormat.rdat")
