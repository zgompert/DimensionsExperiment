##trying standard (non-Bayes) lasso
library(glmnet)
library(RColorBrewer)
library(fields)
load("polyScores.rdat") ## pre-computed polygenic scores

load("glmnetLasso.rdat")

pco<-prcomp(combPs,center=TRUE,scale=TRUE)
summary(pco)
#Importance of components:
#                           PC1     PC2     PC3     PC4     PC5    PC6    PC7
#Standard deviation     13.4792 11.0834 9.38518 8.42750 7.73141 7.1569 6.3626
#Proportion of Variance  0.1032  0.0698 0.05005 0.04035 0.03396 0.0291 0.0230
#Cumulative Proportion   0.1032  0.1730 0.22308 0.26343 0.29739 0.3265 0.3495
oo<-summary(pco)

## PCs, only 1064
pcPs<-pco$x
for(i in 1:1064){
	pcPs[,i]<-pcPs[,i]/sd(pcPs[,i])
}


out<-vector("list",9)
lambda<-rep(NA,9)
r2<-rep(NA,9)
cvr2<-rep(NA,9)
N<-dim(psCat)[1]
for(i in 1:9){
	## perform k-vold CV to find optimal value of lambda
	cv_model<-cv.glmnet(x=pcPs,y=psCat[,i],alpha=1,nfolds=10)
	lambda[i]<-cv_model$lambda.min
	out[[i]]<-glmnet(x=pcPs,y=psCat[,i],alpha=1,lambda=lambda[i])
	y_pred<-predict(out[[i]],s=lambda[i],newx=pcPs)
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	r2[i]<-1-sse/sst
	## 10-fold CV predictive error
	kset<-sample(1:10,N,replace=TRUE)
	for(k in 1:10){
		kk<-which(kset==i)
		ocv<-glmnet(x=pcPs[-kk,],y=psCat[-kk,i],alpha=1,lambda=lambda[i])
		y_pred[kk]<-predict(ocv,s=lambda[i],newx=pcPs[kk,])
	}
	sst<-sum((psCat[,i]-mean(psCat[,i]))^2)
	sse<-sum((y_pred-psCat[,i])^2)
	cvr2[i]<-1-sse/sst
}
save(list=ls(),file="glmnetLasso_PCA.rdat")

summary(r2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2657  0.4746  0.5660  0.5510  0.6373  0.7591 
 summary(cvr2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2531  0.4537  0.5391  0.5238  0.6050  0.7176 



############## figure ##############

traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
cm<-1.5;ca<-1.1;cl<-1.5


pdf("sfig_lasso.pdf_pca",width=10,height=10)

par(mfrow=c(4,2))
par(mar=c(4.5,5.5,2.5,1.5))

xx<-c(1,3:5,7:9)

for(i in 1:7){
    plot(out[[xx[i]]]$beta,col="darkgray",type='h',xlab="PC",ylab="Regression coef.",cex.lab=cl,axes=FALSE)
    axis(2,cex.axis=ca)
    box()
    title(main=paste("(",letters[i],") ",traits[xx[i]],sep=""),cex.main=cm)


}

dev.off()

