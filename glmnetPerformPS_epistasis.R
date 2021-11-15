library(data.table)
library(glmnet)

## load pre-computed polygenic scores for plant traits
load("polyScores.rdat")

ids<-read.table("../../Entropy/MsativaIds.txt",header=FALSE)
fp<-grep(pattern="FP",ids[,1])
fpids<-as.character(ids[fp,1])

psCat<-psCat[fp,] ## already saved wit subset for combPs
for(i in 1:9){
        o<-scale(psCat[,i])
        psCat[,i]<-o[,1]
}

## read in caterpillar genotype and trait data
LmelG<-fread("../lmel_geno",header=FALSE)
LmelP<-read.table("../gemmalmel_pheno_residTraits.txt",header=FALSE)

LmelGM<-as.matrix(LmelG[,-c(1:3)])
pco<-prcomp(t(LmelGM),center=TRUE,scale=FALSE)

## PCA (covariance matrix) on melissa genotypes
#Importance of components:
#                            PC1      PC2     PC3     PC4     PC5     PC6
#Standard deviation     15.67133 10.36591 8.67647 8.48077 7.64412 7.46810
#Proportion of Variance  0.07259  0.03176 0.02225 0.02126 0.01727 0.01649
#Cumulative Proportion   0.07259  0.10435 0.12660 0.14786 0.16513 0.18162
# plan to use first 4,  all > 2% of variance, or first 6, seem to still be pulling pops out

## match plants and caterpillars
cnames<-read.table("../../Entropy/SortedCaterpillarSamples.txt",header=FALSE)

pcMatrix<-matrix(NA,nrow=length(fp),ncol=6)
phMatrix<-matrix(NA,nrow=length(fp),ncol=9)

for(i in 1:length(fp)){
	a<-which(as.character(cnames[,1])==as.character(fpids[i]))
	if(length(a)==1){
		pcMatrix[i,]<-pco$x[a,1:6]
		phMatrix[i,]<-as.numeric(LmelP[a,])
	}
}

keep<-which(is.na(pcMatrix[,1])==FALSE)

#Y<-phMatrix[keep,]
Y<-psCat[keep,]
XPC<-pcMatrix[keep,]
XPS<-combPs[keep,]

## scale everything
for(i in 1:9){
	Y[,i]<-scale(Y[,i])
}
for(i in 1:6){
	XPC[,i]<-scale(XPC[,i])
}
for(i in 1:1760){
	XPS[,i]<-scale(XPS[,i])
}

## make interactions
XINT_LIST<-vector("list",6)
for(i in 1:6){
	XINT_LIST[[i]]<-XPS
	for(j in 1:1760){
		XINT_LIST[[i]][,j]<-XPS[,j]*XPC[,i]
	}
}

XX<-as.matrix(cbind(XPC,XPS,XINT_LIST[[1]],XINT_LIST[[2]],XINT_LIST[[3]],XINT_LIST[[4]]))
type<-c(rep(1,6),rep(2:6,each=1760))
#XX<-as.matrix(XPS)

#Next, we’ll use the glmnet() function to fit the lasso regression model and specify alpha=1.

out<-vector("list",9)
lambda<-rep(NA,9)
r2<-rep(NA,9)
cvr2<-rep(NA,9)
N<-dim(XX)[1]
for(i in 1:9){
	## perform k-vold CV to find optimal value of lambda
	full<-which(is.na(Y[,i])==FALSE)
	cv_model<-cv.glmnet(x=XX[full,],y=Y[full,i],alpha=1,nfolds=10)
	lambda[i]<-cv_model$lambda.min
	out[[i]]<-glmnet(x=XX[full,],y=Y[full,i],alpha=1,lambda=lambda[i])
	y_pred<-predict(out[[i]],s=lambda[i],newx=XX[full,])
	sst<-sum((Y[full,i]-mean(Y[full,i]))^2)
	sse<-sum((y_pred-Y[full,i])^2)
	r2[i]<-1-sse/sst
	## 10-fold CV predictive error
	kset<-sample(1:10,N,replace=TRUE)
	for(k in 1:10){
		kk<-which(kset!=k & is.na(Y[,i])==FALSE)
		ocv<-glmnet(x=XX[kk,],y=Y[kk,i],alpha=1,lambda=lambda[i])
		y_pred[kset==k]<-predict(ocv,s=lambda[i],newx=XX[kset==k,])
	}
	sst<-sum((Y[,i]-mean(Y[,i],na.rm=TRUE))^2,na.rm=TRUE)
	sse<-sum((y_pred-Y[,i])^2,na.rm=TRUE)
	cvr2[i]<-1-sse/sst
}
## out has the coefficients, which is the main thing I want, key is the includes some interaction terms, but greater prop that arent'
barplot(tapply(X=out[[2]]$beta != 0,type,mean))
## note, r2 is generally lower than for the models with only PS

traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
library(RColorBrewer)
cs<-c("brown","black",brewer.pal(n=6,"BuPu")[-c(1:2)])
## PC alone never stands out
cm<-1.5;cl<-1.5;ca<-1.2
pdf("sfig_glmnetPSCS.pdf",width=9,height=7)
par(mfrow=c(3,3))
par(mar=c(4.5,5,2.5,1))
for(i in 1:9){
	plot(out[[i]]$beta,type='h',col=cs[type],xlab="Model covariate",ylab="Coefficient",cex.lab=cl,cex.axis=ca)
	title(main=paste("(",letters[i],") ",traits[i],sep=""),cex.main=cm)
	if(i==1){
		legend(5500,0.14,c("PC","PGS","PC1xPGS","PC2xPGS","PC3xPGS","PC4xPGS"),fill=cs,bty='n',cex=.9)
	}
}
dev.off()

## prop retained PGS vs interaction
props<-matrix(NA,nrow=9,ncol=2)
for(i in 1:9){
	props[i,1]<-mean(out[[i]]$beta[7:1767]>0)
	props[i,2]<-mean(out[[i]]$beta[-c(1:1767)]>0)
}

pdf("sfig_lassoRetained.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
plot(props[,1],props[,2],xlim=c(0,0.06),ylim=c(0,0.06),pch=19,xlab="Prop. main retained",ylab="Prop. interactions retained",cex.lab=1.4,cex.axis=1.1)
abline(a=0,b=1)
dev.off()
## checking with randomized data, just to be safe
ranR2<-matrix(NA,nrow=10,ncol=9)
for(i in 1:9){for(j in 1:10){
	aran<-sample(1:N,N,replace=FALSE)
	Yran<-Y[aran,]
	full<-which(is.na(Yran[,i])==FALSE)
	cv_model<-cv.glmnet(x=XX[full,],y=Yran[full,i],alpha=1,nfolds=10)
	ran_lambda<-cv_model$lambda.min
	ran_out<-glmnet(x=XX[full,],y=Yran[full,i],alpha=1,lambda=ran_lambda)
	y_pred<-predict(ran_out,s=ran_lambda,newx=XX[full,])
	sst<-sum((Yran[full,i]-mean(Yran[full,i]))^2)
	sse<-sum((y_pred-Yran[full,i])^2)
	ranR2[j,i]<-1-sse/sst
}}## 14 day weight is solid, rest less so but this is probably the best we got
save(list=ls(),file="glmnetLasso_epi.rdat")

