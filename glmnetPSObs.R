## comput polygenic scores and realte to observed caterpillar phenotypes 
library(data.table)
library(glmnet)
library(scales)
G<-fread("../msat_geno",header=FALSE)
G<-t(as.matrix(G[,-c(1:3)]))
N<-dim(G)[1]
L<-dim(G)[2]
## center G
for(i in 1:L){
	G[,i]<-G[,i]-mean(G[,i])
}


## chemistry
Np<-1750 ## note in alpha order not strict numeric
psChem<-matrix(NA,nrow=N,ncol=Np)
files_chem<-list.files(pattern="mav_o_msat_fit_gemma_pheno_residChem_ph")
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_chem[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psChem[,j]<-G %*% ph
}

## plant traits
files_plant<-list.files(pattern="mav_o_msat_fit_gemma_pheno_residPlant_ph")
Np<-10 ## note in alpha order not strict numeric
psPlant<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_plant[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psPlant[,j]<-G %*% ph
}

## cat traits
files_cat<-list.files(pattern="mav_o_msat_fit_gemma_pheno_residTraits_ph")
Np<-9
psCat<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_cat[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psCat[,j]<-G %*% ph
}

## ran chemistry
Np<-1750 ## note in alpha order not strict numeric
psRanChem<-matrix(NA,nrow=N,ncol=Np)
files_chem<-list.files(pattern="mav_o_msat_fit_random_pheno_residChem_ph")
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_chem[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psRanChem[,j]<-G %*% ph
}

## ran plants
files_cat<-list.files(pattern="mav_o_msat_fit_random_pheno_residPlant_ph")
Np<-10
psRanPlant<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
	cat(j,"\n")
	ph<-fread(files_cat[j],header=FALSE)
	ph<-unlist(as.vector(ph[,2]))# get trait
	psRanPlant[,j]<-G %*% ph
}


ph<-read.table("../gemma_pheno_residTraits.txt",header=FALSE)

## ps correlations, fp only
chem_num<-scan("chemPhNum.txt") ## ph chem order
psCors<-matrix(NA,nrow=9,ncol=10+1750)
phCors<-matrix(NA,nrow=9,ncol=10+1750)
combPs<-cbind(psChem[,order(chem_num)],psPlant[,c(1,3:10,2)])
combRanPs<-cbind(psRanChem[,order(chem_num)],psRanPlant[,c(1,3:10,2)])
for(i in 1:9){for(j in 1:1760){
	cat(i,"\n")
	psCors[i,j]<-cor(psCat[,i],combPs[,j])
	o<-cor.test(ph[,i],combPs[,j])
	phCors[i,j]<-o$estimate
}}

## read in h2
chemH2<-read.table("../files_post_chem/comb_ord_post_chem.txt",header=FALSE)
plantH2<-read.table("../output/post_gemma_resid_plant.txt",header=FALSE)
h2<-c(chemH2[,1],plantH2[,1])
catH2<-read.table("../output/post_gemma_resid_traits.txt",header=FALSE)
cTraits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")

save(list=ls(),file="polyScoresV2.rdat")

## plot

pdf("MsatGenPhCors.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(4,5,2.5,1.5))
for(i in 1:9){
	plot(h2,phCors[i,],pch=19,col=c(rep(alpha("black",.3),1750),rep(alpha("forestgreen",.9),10)),xlab="Trait PVE",ylab="Genetic correlation",cex.lab=1.5)
	title(main=paste(cTraits[i],", PVE =",round(catH2[i,1],3)),cex.main=1.5)
}
dev.off()

for(i in 1:dim(combPs)[2]){
	combPs[,i]<-as.numeric(scale(combPs[,i]))
}

ph<-read.table("../gemma_pheno_residTraits.txt",header=FALSE)


## glmnet scaling, on pheno
out<-vector("list",9)
lambda<-rep(NA,9)
r2<-rep(NA,9)
cvr2<-rep(NA,9)
N<-dim(ph)[1]
for(i in 1:9){
        ## perform k-vold CV to find optimal value of lambda
        nona<-which(is.na(ph[,i])==FALSE)
        nas<-which(is.na(ph[,i])==TRUE)
	cv_model<-cv.glmnet(x=combPs[nona,],y=ph[nona,i],alpha=1,nfolds=10)
        lambda[i]<-cv_model$lambda.min
        out[[i]]<-glmnet(x=combPs[nona,],y=ph[nona,i],alpha=1,lambda=lambda[i])
        y_pred<-predict(out[[i]],s=lambda[i],newx=combPs)
        sst<-sum((ph[,i]-mean(ph[,i],na.rm=TRUE))^2,na.rm=TRUE)
        sse<-sum((y_pred-ph[,i])^2,na.rm=TRUE)
        r2[i]<-1-sse/sst
        ## 10-fold CV predictive error
        kset<-sample(1:10,N,replace=TRUE)
        for(k in 1:10){
                kk<-which(kset==k)
                ocv<-glmnet(x=combPs[-unique(c(kk,nas)),],y=ph[-unique(c(kk,nas)),i],alpha=1,lambda=lambda[i])
                y_pred[kk]<-predict(ocv,s=lambda[i],newx=combPs[kk,])
        }
        sst<-sum((ph[,i]-mean(ph[,i],na.rm=TRUE))^2,na.rm=TRUE)
        sse<-sum((y_pred-ph[,i])^2,na.rm=TRUE)
        cvr2[i]<-1-sse/sst
}

## with randomization ##
ranR2<-matrix(NA,nrow=10,ncol=9)
for(i in 1:9){for(j in 1:10){
        aran<-sample(1:1055,1055,replace=FALSE)
        nona<-which(is.na(ph[aran,i])==FALSE)
	cv_model<-cv.glmnet(x=combPs[nona,],y=ph[aran,i][nona],alpha=1,nfolds=10)
        ran_lambda<-cv_model$lambda.min
        ran_out<-glmnet(x=combPs[nona,],y=ph[aran,i][nona],alpha=1,lambda=ran_lambda)
        y_pred<-predict(ran_out,s=ran_lambda,newx=combPs)
        sst<-sum((ph[aran,i]-mean(ph[,i],na.rm=TRUE))^2,na.rm=TRUE)
        sse<-sum((y_pred-ph[aran,i])^2,na.rm=TRUE)
        ranR2[j,i]<-1-sse/sst
}}## 


save(list=ls(),file="polyScoresV2.rdat")

############# plots #################
traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
cm<-1.5;ca<-1.1;cl<-1.5

cs<-c(alpha("darkslategray4",.5),alpha("orange4",.9))

pdf("sfig_lasso_pheno.pdf",width=10,height=8)

par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))

xx<-c(1:4,6:7)

for(i in 1:6){
    plot(out[[xx[i]]]$beta,col=c(rep(cs[1],1750),rep(cs[2],10)),type='h',xlab="Plant trait",ylab="Regression coef.",cex.lab=cl,axes=FALSE)
    axis(2,cex.axis=ca)
    box()
    title(main=paste("(",letters[i],") ",traits[xx[i]],sep=""),cex.main=cm)
    if(i==6){legend(.2,.8,c("chemistry","other"),cex=1.3,fill=cs,bty='n')}


}

legend(1,-.02,c("chemistry","other"),cex=1.5,fill=cs,bty='n')

dev.off()

## repeated for more precise p-values
ranR2_100<-matrix(NA,nrow=100,ncol=9)
for(i in 1:9){for(j in 1:100){
        aran<-sample(1:1055,1055,replace=FALSE)
        nona<-which(is.na(ph[aran,i])==FALSE)
	cv_model<-cv.glmnet(x=combPs[nona,],y=ph[aran,i][nona],alpha=1,nfolds=10)
        ran_lambda<-cv_model$lambda.min
        ran_out<-glmnet(x=combPs[nona,],y=ph[aran,i][nona],alpha=1,lambda=ran_lambda)
        y_pred<-predict(ran_out,s=ran_lambda,newx=combPs)
        sst<-sum((ph[aran,i]-mean(ph[,i],na.rm=TRUE))^2,na.rm=TRUE)
        sse<-sum((y_pred-ph[aran,i])^2,na.rm=TRUE)
        ranR2_100[j,i]<-1-sse/sst
}}## 


ps<-rep(NA,9)
for(i in 1:9){
	ps[i]<-mean(c(ranR2_100[,i],r2[i]) >= r2[i])
}
ps
#[1] 0.00990099 0.03960396 0.13861386 0.40594059 0.86138614 0.02970297 0.05940594
#[8] 1.00000000 1.00000000
 -2 * sum(log(ps))
#[1] 34.42128
pchisq(q=34.42128,df=18,lower.tail=FALSE)
#[1] 0.01116817

library(scales)

pdf("sfig_phVarExplain.pdf",width=7,height=5)

par(mar=c(4.5,5.5,2.5,1.5))

plot(1:9,r2,typ='n',xlab="",ylab="Variance explained (r2)",cex.lab=cl,cex.axis=ca,axes=FALSE)
x<-1:9
axis(1,at=1:9,traits,cex.axis=ca)
axis(2,cex.axis=ca)
box()
for(i in 1:100){
points(jitter(x),ranR2_100[i,],pch=19,col=alpha("gray",.6))
}
points(x,r2,pch=19)
dev.off()
