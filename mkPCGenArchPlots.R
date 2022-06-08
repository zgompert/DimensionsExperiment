## plots summarizing genetic architecture of caterpillar performance, with PCS
library(RColorBrewer)
library(scales)
library(data.table)
library(xtable)

post_plant<-read.table("output_main/post_gemma_resid_traits.txt",header=FALSE)
post_cat<-read.table("output_cat/post_gemma_resid_lmel.txttxt",header=FALSE)
post_pc_plant<-read.table("output/post_gemma_resid_msat_traits.txt",header=FALSE)
post_pc_cat<-read.table("output/post_gemma_resid_lmel_traits.txt",header=FALSE)

traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
pve<-cbind(post_plant[,1],post_cat[,1])
lb<-cbind(post_plant[,2],post_cat[,2])
ub<-cbind(post_plant[,3],post_cat[,3])

pve_pc<-cbind(post_pc_plant[,1],post_pc_cat[,1])
lb_pc<-cbind(post_pc_plant[,2],post_pc_cat[,2])
ub_pc<-cbind(post_pc_plant[,3],post_pc_cat[,3])

## comparing genetic arch estimates
cl<-1.45;ca<-1.1;cm<-1.5
pdf("sfig_pve_comps.pdf",width=4.5,height=9)
par(mfrow=c(2,1))
par(mar=c(4,5,2,1))
plot(pve[,1],pve_pc[,1],xlab="PVE without PCs",ylab="PVE with PCs",cex.lab=cl,cex.axis=ca,pch=19)
abline(a=0,b=1)
title(main="(a) M. sativa genetics",cex.main=cm)
plot(pve[,2],pve_pc[,2],xlab="PVE without PCs",ylab="PVE with PCs",cex.lab=cl,cex.axis=ca,pch=19)
abline(a=0,b=1)
title(main="(b) L. melissa genetics",cex.main=cm)
dev.off()

cor.test(pve[,1],pve_pc[,1])
#t = 54.483, df = 7, p-value = 1.84e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9941821 0.9997623
#sample estimates:
#     cor
#0.998823
cor.test(pve[,2],pve_pc[,2])
#t = 97.165, df = 7, p-value = 3.223e-12
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9981656 0.9999252
#sample estimates:
#      cor 
#0.9996295 


L_ms<-161008
L_lm<-63194

pf_ms<-list.files(path="output",full.names=TRUE,pattern="pips_o_msat")
pf_lm<-list.files(path="output",full.names=TRUE,pattern="pips_o_lmel")

pips_ms<-matrix(NA,nrow=L_ms,ncol=9)
pips_lm<-matrix(NA,nrow=L_lm,ncol=9)
for(i in 1:9){
	pdat<-read.table(pf_ms[i],header=FALSE)
	pdat<-pdat[-c(1:20),]
	pips_ms[,i]<-as.vector(pdat[,2])
	pdat<-read.table(pf_lm[i],header=FALSE)
	pdat<-pdat[-c(1:20),]
	pips_lm[,i]<-as.vector(pdat[,2])
}


sum(as.vector(pips_ms))
#[1] 433.8863 
sum(as.vector(pips_lm))
#[1] 205.0008
sum(as.vector(pips_lm) > 0.01)
#[1] 328
sum(as.vector(pips_ms) > 0.01)
#[1] 121

apply(pips_ms,2,sum)
#[1] 32.02648 80.28788 93.17584 54.58061 24.44304 31.77762 53.26780 17.60433
#[9] 46.72273
apply(pips_lm,2,sum)
#4.143801 40.564364 60.131311 10.253278  6.874558 53.082943 12.221606
#[8]  6.065982  1.662981
apply(pips_ms > 0.01,2,sum)
#9 60  2  7  5 17 18  2  1
apply(pips_lm > 0.01,2,sum)
#1] 98 28  0 60 58  7 34 35  8

### PLOTS ####
traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")

cl<-1.55;ca<-1.1;cm<-1.6
pdf("sfig_pca_genetics.pdf",width=9*1.05,height=5*1.05)
layout(matrix(c(1,1,2,2),nrow=2,ncol=2,byrow=TRUE),widths=c(5,5),heights=c(3,3))
par(mar=c(4.5,4.5,2.5,2.5))
## panel A
cs<-brewer.pal(n=3,"Dark2")


## alfalfa
px<-c(19,15,17,19,15,17,18,19,15)
csx<-alpha(c(rep("orange",3),rep("cadetblue",4),rep("darkorchid",2)),.7)
msat_snps<-read.table("msat_lg_snps",header=FALSE)
bnds<-c(1,which(msat_snps[-1,1] != msat_snps[-L_ms,1]),L_ms)
mids<-tapply(X=1:L_ms,INDEX=msat_snps[,1],mean)
cl<-1.5;ca<-1.1;cm<-1.5

plot(1:L_ms,rep(1,L_ms),type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
axis(1,at=mids,c(1:8),cex.axis=ca)
axis(2,at=c(0,0.5,1),cex.axis=ca)
box()
for(i in seq(2,8,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("grey",.25),border=NA)
}
for(i in 1:9){
a<-which(pips_ms[,i] > 0.01)
points(c(1:L_ms)[a],pips_ms[a,i],pch=px[i],col=csx[i])}
legend(115000,0.98,traits,pch=px,col=csx,bty='n',ncol=3)
abline(h=c(.1,.5),lty=2)

title(main="(a) M. sativa genotype-performance associations",cex.main=cm)

# melissa
mel_snps<-read.table("lmel_lg_snps",header=FALSE)
oo<-order(mel_snps[,1],mel_snps[,2],mel_snps[,3])
o_pips_lm<-pips_lm[oo,]

px<-c(19,15,17,19,15,17,18,19,15)
csx<-alpha(c(rep("orange",3),rep("cadetblue",4),rep("darkorchid",2)),.7)
lgs<-mel_snps[oo,1]
lgs[is.na(lgs)]<-24
bnds<-c(1,which(lgs[-1] != lgs[-L_lm]),L_lm)
mids<-tapply(X=1:L_lm,INDEX=lgs,mean)


plot(1:L_lm,rep(1,L_lm),type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
axis(1,at=mids,c(1:22,"Z","NA"),cex.axis=ca)
axis(2,at=c(0,0.5,1),cex.axis=ca)
box()
for(i in seq(2,24,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("grey",.25),border=NA)
}

for(i in 1:9){
a<-which(o_pips_lm[,i] > 0.01)
points(c(1:L_lm)[a],o_pips_lm[a,i],pch=px[i],col=csx[i])}
#legend(35000,0.98,traits,pch=px,col=csx,bty='n',ncol=3)
title(main="(b) L. melissa genotype-performance associations",cex.main=cm)
abline(h=c(.1,.5),lty=2)

dev.off()

pf_ms<-list.files(path="files_pips",full.names=TRUE,pattern="pips_o_msat_fit_gemma_pheno_residTraits_ph")
pf_lm<-list.files(path="files_pips",full.names=TRUE,pattern="pips_o_lmel_fit_gemmalmel_pheno_residTraits_ph")


pips_ms_s<-matrix(NA,nrow=L_ms,ncol=9)
pips_lm_s<-matrix(NA,nrow=L_lm,ncol=9)
for(i in 1:9){
	pdat<-read.table(pf_ms[i],header=FALSE)
	pips_ms_s[,i]<-as.vector(pdat[,2])
	pdat<-read.table(pf_lm[i],header=FALSE)
	pips_lm_s[,i]<-as.vector(pdat[,2])
}
cor.test(as.vector(pips_lm),as.vector(pips_lm_s))
#t = 3414.4, df = 568744, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9763436 0.9765854
#sample estimates:
#      cor 
#0.9764648 
cor.test(as.vector(pips_ms),as.vector(pips_ms_s))
#data:  as.vector(pips_ms) and as.vector(pips_ms_s)
#t = 1266.5, df = 1449070, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7240434 0.7255890
#sample estimates:
#      cor 
#0.7248171 


