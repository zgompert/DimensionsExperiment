##trying standard (non-Bayes) lasso
library(glmnet)
library(RColorBrewer)
library(fields)
load("polyScores.rdat") ## pre-computed polygenic scores

#Next, we’ll use the glmnet() function to fit the lasso regression model and specify alpha=1.
K<-dim(combPs)

## scale covariates
for(i in 1:K[2]){
	o<-scale(combPs[,i])
	combPs[,i]<-o[,1]
}

## scale response
ids<-read.table("../../Entropy/MsativaIds.txt")
fp<-grep(pattern="FP",ids[,1])
psCat<-psCat[fp,] ## already saved wit subset for combPs
for(i in 1:9){
	o<-scale(psCat[,i])
	psCat[,i]<-o[,1]
}

## G matrix
gs<-cor(as.matrix(cbind(psCat,combPs)))
gcs<-rev(brewer.pal(n=10,"RdGy"))
pdf("sfig_gmat.pdf",width=7,height=6)
par(mar=c(1,1,1,1))
image.plot(gs,breaks=c(-1,-.75,-.5,-.25,-.1,0,.1,.25,.5,.75,1),col=gcs,axes=FALSE)
box()
dev.off()

## all of the correlations
psCors<-matrix(NA,nrow=1760,ncol=9)
for(i in 1:9){for(j in 1:1760){
	psCors[j,i]<-cor(psCat[,i],combPs[,j])
}}
## order is 1750 chem followed by 10 plant traits

## h2
chemH2<-read.table("../files_post_chem/comb_ord_post_chem.txt",header=FALSE)
plantH2<-read.table("../output/post_gemma_resid_plant.txt",header=FALSE)
h2<-c(chemH2[,1],plantH2[,1])
catH2<-read.table("../output/post_gemma_resid_traits.txt",header=FALSE)



## scale plant alone
ids<-read.table("../../Entropy/MsativaIds.txt")
fp<-grep(pattern="FP",ids[,1])
psPlant<-psPlant[fp,] ## already saved wit subset for combPs
for(i in 1:9){
	o<-scale(psPlant[,i])
	psPlant[,i]<-o[,1]
}

out<-vector("list",9)
lambda<-rep(NA,9)
r2<-rep(NA,9)
cvr2<-rep(NA,9)
N<-dim(psCat)[1]
for(i in 1:9){
	## perform k-vold CV to find optimal value of lambda
	cv_model<-cv.glmnet(x=combPs,y=psCat[,i],alpha=1,nfolds=10)
	lambda[i]<-cv_model$lambda.min
	out[[i]]<-glmnet(x=combPs,y=psCat[,i],alpha=1,lambda=lambda[i])
	y_pred<-predict(out[[i]],s=lambda[i],newx=combPs)
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	r2[i]<-1-sse/sst
	## 10-fold CV predictive error
	kset<-sample(1:10,N,replace=TRUE)
	for(k in 1:10){
		kk<-which(kset==i)
		ocv<-glmnet(x=combPs[-kk,],y=psCat[-kk,i],alpha=1,lambda=lambda[i])
		y_pred[kk]<-predict(ocv,s=lambda[i],newx=combPs[kk,])
	}
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	cvr2[i]<-1-sse/sst
}

## checking with randomized data, just to be safe
ranR2<-matrix(NA,nrow=10,ncol=9)
for(i in 1:9){for(j in 1:10){
	aran<-sample(1:1064,1064,replace=FALSE)
	cv_model<-cv.glmnet(x=combPs,y=psCat[aran,i],alpha=1,nfolds=10)
	ran_lambda<-cv_model$lambda.min
	ran_out<-glmnet(x=combPs,y=psCat[aran,i],alpha=1,lambda=ran_lambda)
	y_pred<-predict(ran_out,s=ran_lambda,newx=combPs)
	sst<-sum((psCat[aran,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[aran,i])^2)
	ranR2[j,i]<-1-sse/sst
}}## wow! this is real!!
save(list=ls(),file="glmnetLasso.rdat")


############ chosing top chemicals for annotation ####
## matrix of abs betas
betas<-matrix(NA,nrow=1750,ncol=9)
betasRaw<-matrix(NA,nrow=1750,ncol=9)
for(i in 1:9){
	betas[,i]<-abs(out[[i]]$beta[1:1750]) ## chem only
	betasRaw[,i]<-out[[i]]$beta[1:1750] 
}
mnBeta<-apply(betas,1,mean)
wBetas<-betas
for(i in 1:9){ ## weighted by h2 for chemical and pve for perform trait
	wBetas[,i]<-catH2[i,1] * (betas[,i] * h2[1:1750])
}
mnWgBeta<-apply(wBetas,1,mean)

## get top 20 
q<-20/1750
qv<-quantile(mnWgBeta,probs=1-q)
tchem<-which(mnWgBeta >= qv)


cdat<-read.table("../../Pheno/MESA_UT_ZG_2019_normalized.csv",header=TRUE,sep=",")
cnms<-colnames(cdat)[-c(1:6)]
cnms[tchem]
# [1] "PC.38.9"                                                                  
# [2] "Apigenin.4...p.coumaroyl...2..glucuronyl..1.2..glucuronide..7.glucuronide"
# [3] "Quillaic.acid.3..rhamnosyl..1.3...galactosyl..1.2...glucuronide."         
# [4] "MESA.1112"                                                                
# [5] "MESA.122"                                                                 
# [6] "MESA.124"                                                                 
# [7] "MESA.1305"                                                                
# [8] "MESA.1342"                                                                
# [9] "MESA.143"                                                                 
#[10] "MESA.339"                                                                 
#[11] "MESA.438"                                                                 
#[12] "MESA.545"                                                                 
#[13] "MESA.615"                                                                 
#[14] "MESA.730"                                                                 
#[15] "MESA.784"                                                                 
#[16] "MESA.849"                                                                 
#[17] "MESA.583"                                                                 
#[18] "MESA.584"                                                                 
#[19] "MESA.972"                                                                 
#[20] "MESA.1185"
#######################################################



## trying ps from randomized chem and plant traits to predict actual
## cat performance bvs
load("ranPs.rdat")

ids<-read.table("../../Entropy/MsativaIds.txt")
fp<-grep(pattern="FP",ids[,1])
psRanPlant<-psRanPlant[fp,] 
for(i in 1:10){
	o<-scale(psRanPlant[,i])
	psRanPlant[,i]<-o[,1]
}
psRanChem<-psRanChem[fp,] 
for(i in 1:1750){
	o<-scale(psRanChem[,i])
	psRanChem[,i]<-o[,1]
}

combRanPs<-as.matrix(cbind(psRanChem,psRanPlant))

outRan<-vector("list",9)
lambdaRan<-rep(NA,9)
r2Ran<-rep(NA,9)
cvr2Ran<-rep(NA,9)
N<-dim(psCat)[1]
for(i in 1:9){
	## perform k-vold CV to find optimal value of lambda
	cv_model<-cv.glmnet(x=combRanPs,y=psCat[,i],alpha=1,nfolds=10)
	lambdaRan[i]<-cv_model$lambda.min
	outRan[[i]]<-glmnet(x=combRanPs,y=psCat[,i],alpha=1,lambda=lambdaRan[i])
	y_pred<-predict(outRan[[i]],s=lambdaRan[i],newx=combRanPs)
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	r2Ran[i]<-1-sse/sst
	## 10-fold CV predictive error
	kset<-sample(1:10,N,replace=TRUE)
	for(k in 1:10){
		kk<-which(kset==i)
		ocv<-glmnet(x=combRanPs[-kk,],y=psCat[-kk,i],alpha=1,lambda=lambdaRan[i])
		y_pred[kk]<-predict(ocv,s=lambdaRan[i],newx=combRanPs[kk,])
	}
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	cvr2Ran[i]<-1-sse/sst
}




############## main figure ##############
traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
cm<-1.5;ca<-1.1;cl<-1.5

cs<-c(alpha("darkslategray4",.5),alpha("orange4",.9))

pdf("fig_plantPs2perform.pdf",width=7,height=11.375)

layout(matrix(c(1,2,3,3,4,4,5,5),nrow=4,ncol=2,byrow=TRUE),widths=c(4,4),heights=c(4,3,3,3))
par(mar=c(4.5,5.5,2.5,1.5))

plot(h2,psCors[,2],xlab="Plant trait PVE",ylab="Genetic correlation",pch=19,col=c(rep(cs[1],1750),rep(cs[2],10)),cex.lab=cl)
legend(.54,-.22,c("chemistry","other"),pch=19,cex=1.2,col=cs,bty='n')
abline(h=0,lty=2)
title(main="(a) 14-day weight (PVE = 0.36)",cex.main=cm)

plot(h2,psCors[,6],xlab="Plant trait PVE",ylab="Genetic correlation",pch=19,col=c(rep(cs[1],1750),rep(cs[2],10)),cex.lab=cl)
abline(h=0,lty=2)
title(main="(b) Surv. to pupa. (PVE = 0.21)",cex.main=cm)


plot(1:9,r2,ylim=c(0,1),pch=19,xlab="",ylab="Variance explained (r2)",cex.lab=cl,cex.axis=ca,axes=FALSE)
axis(1,at=1:9,traits,cex.axis=ca)
axis(2,cex.axis=ca)
box()
for(i in 1:10){points(jitter(1:9),ranR2[i,],pch=19,col="gray")}
title(main="(c) Variation explained by plant trait polygenic scores",cex.main=cm)

plot(out[[2]]$beta,col=c(rep(cs[1],1750),rep(cs[2],10)),type='h',xlab="Plant trait",ylab="Regression coef.",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
box()
title(main="(d) Effect of plant trait polygenic scores on 14-day weight",cex.main=cm)

plot(out[[6]]$beta,col=c(rep(cs[1],1750),rep(cs[2],10)),type='h',xlab="Plant trait",ylab="Regression coef.",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
box()
title(main="(e) Effect of plant trait polygenic scores on survival to pupation",cex.main=cm)
legend(-25,-.055,c("chemistry","other"),cex=1.2,fill=cs,bty='n')

dev.off()


##########################################

############## some supp. figures #######

## all the PS cor. plots
traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
cm<-1.5;ca<-1.1;cl<-1.5

cs<-c(alpha("darkslategray4",.5),alpha("orange4",.9))

pdf("sfig_MsatGenCors.pdf",width=10,height=10)

par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1:9){
    plot(h2,psCors[i,],xlab="Plant trait PVE",ylab="Genetic correlation",pch=19,col=c(rep(cs[1],1750),rep(cs[2],10)),cex.lab=cl)
    if(i==1){
    legend(.54,-.22,c("chemistry","other"),pch=19,cex=1.2,col=cs,bty='n')}
    abline(h=0,lty=2)
    title(main=paste("(",letters[i],") ",traits[i],sep=""),cex.main=cm)
}
dev.off()

## rest of the model-averaged effect plots
traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
cm<-1.5;ca<-1.1;cl<-1.5

cs<-c(alpha("darkslategray4",.5),alpha("orange4",.9))

pdf("sfig_lasso.pdf",width=10,height=10)

par(mfrow=c(4,2))
par(mar=c(4.5,5.5,2.5,1.5))

xx<-c(1,3:5,7:9)

for(i in 1:7){
    plot(out[[xx[i]]]$beta,col=c(rep(cs[1],1750),rep(cs[2],10)),type='h',xlab="Plant trait",ylab="Regression coef.",cex.lab=cl,axes=FALSE)
    axis(2,cex.axis=ca)
    box()
    title(main=paste("(",letters[i],") ",traits[xx[i]],sep=""),cex.main=cm)


}

plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
legend(.2,.8,c("chemistry","other"),cex=1.3,fill=cs,bty='n')

dev.off()


## plants only below?

## trying with only plant data
out_plant<-vector("list",9)
lambda_plant<-rep(NA,9)
r2_plant<-rep(NA,9)
cvr2_plant<-rep(NA,9)
N<-dim(psCat)[1]
o_simple<-vector("list",9)
for(i in 1:9){
	## perform k-vold CV to find optimal value of lambda
	o_simple[[i]]<-lm(psCat[,i] ~ psPlant)
	cv_model<-cv.glmnet(x=psPlant,y=psCat[,i],alpha=1,nfolds=10)
	lambda_plant[i]<-cv_model$lambda.min
	out_plant[[i]]<-glmnet(x=psPlant,y=psCat[,i],alpha=1,lambda=lambda_plant[i])
	y_pred<-predict(out_plant[[i]],s=lambda_plant[i],newx=psPlant)
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	r2_plant[i]<-1-sse/sst
	## 10-fold CV predictive error
	kset<-sample(1:10,N,replace=TRUE)
	for(k in 1:10){
		kk<-which(kset==i)
		ocv<-glmnet(x=psPlant[-kk,],y=psCat[-kk,i],alpha=1,lambda=lambda_plant[i])
		y_pred[kk]<-predict(ocv,s=lambda_plant[i],newx=psPlant[kk,])
	}
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	cvr2_plant[i]<-1-sse/sst
}## cvr2, r2 and simple lm r2 all about the same
## and all notably less than with plant and chem, still they do something
