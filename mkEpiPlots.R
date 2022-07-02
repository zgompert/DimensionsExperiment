## plots summarizing genetic architecture of caterpillar performance
library(RColorBrewer)
library(scales)
library(data.table)
library(xtable)

## all residuals, comparing plant, caterpillar and both
post_plant<-read.table("output/post_gemma_resid_traits.txt",header=FALSE)
post_cat<-read.table("output_cat/post_gemma_resid_lmel.txttxt",header=FALSE)
post_comb<-read.table("output_comb/post_resid_combined.txt",header=FALSE)
traits<-c("W8d","W14d","WPup","S8d","S14d","SPup","SAdu","Stot","Stime")
pve<-cbind(post_plant[,1],post_cat[,1],post_comb[,1])
lb<-cbind(post_plant[,2],post_cat[,2],post_comb[,2])
ub<-cbind(post_plant[,3],post_cat[,3],post_comb[,3])

pge<-cbind(post_plant[,1]*post_plant[,4],post_cat[,1]*post_cat[,4],post_comb[,1]*post_comb[,4])

L_ms<-161008
L_lm<-63194

pf_ms<-list.files(path="files_pips",full.names=TRUE,pattern="pips_o_msat_fit_gemma_pheno_residTraits_ph")
pf_lm<-list.files(path="files_pips",full.names=TRUE,pattern="pips_o_lmel_fit_gemmalmel_pheno_residTraits_ph")

pips_ms<-matrix(NA,nrow=L_ms,ncol=9)
pips_lm<-matrix(NA,nrow=L_lm,ncol=9)
for(i in 1:9){
	pdat<-read.table(pf_ms[i],header=FALSE)
	pips_ms[,i]<-as.vector(pdat[,2])
	pdat<-read.table(pf_lm[i],header=FALSE)
	pips_lm[,i]<-as.vector(pdat[,2])
}

## tables
xtable(round(cbind(post_plant[,1:3],post_cat[,1:3],post_comb[,1:3]),3))
xtable(round(cbind(post_plant[,4:6],post_cat[,4:6],post_comb[,4:6]),3))


sum(as.vector(pips_ms))
#[1] 369.5415
sum(as.vector(pips_lm))
#[1] 181.1797
sum(as.vector(pips_lm) > 0.01)
#[1] 326
sum(as.vector(pips_ms) > 0.01)
#[1] 134

apply(pips_ms,2,sum)
#[1] 51.93839 93.98602 42.92918 26.28521 25.07629 18.98927 40.35579 12.07030
#[9] 57.91110
apply(pips_lm,2,sum)
#[1] 13.871191 33.926015 40.937045  9.992100  6.718866 57.014031 11.037565
#[8]  6.000561  1.682349

apply(pips_ms > 0.01,2,sum)
#[1] 13 73  2  5  6 15 12  2  6
apply(pips_lm > 0.01,2,sum)
#[1] 86 27  1 73 58  6 34 34  7

## compute polygenic scores
G<-fread("../Entropy/msat_geno.txt",header=FALSE)
G<-as.matrix(G)
N<-dim(G)[1]
L<-dim(G)[2]
## center G
for(i in 1:L){
        G[,i]<-G[,i]-mean(G[,i])
}

files_cat<-list.files(path="files_mav",pattern="msat_fit_gemma_pheno_residTraits_ph",full.names=TRUE)

Np<-9
psCat<-matrix(NA,nrow=N,ncol=Np)
for(j in 1:Np){
        cat(j,"\n")
        ph<-fread(files_cat[j],header=FALSE)
        ph<-unlist(as.vector(ph[,2]))# get trait
        psCat[,j]<-G %*% ph
}

LmelG<-fread("lmel_geno",header=FALSE)
LmelG<-t(as.matrix(LmelG[,-c(1:3)]))
LmelN<-dim(LmelG)[1]
LmelL<-dim(LmelG)[2]
## center G
for(i in 1:LmelL){
        LmelG[,i]<-LmelG[,i]-mean(LmelG[,i])
}

files_lm<-list.files(path="files_mav",pattern="mav_o_lmel_fit_gemmalmel_pheno_residTraits_ph",full.names=TRUE)

Np<-9
psLm<-matrix(NA,nrow=LmelN,ncol=Np)
for(j in 1:Np){
        cat(j,"\n")
        ph<-fread(files_lm[j],header=FALSE)
        ph<-unlist(as.vector(ph[,2]))# get trait
        psLm[,j]<-LmelG %*% ph
}

## genetic cors
gM<-matrix(NA,nrow=9,ncol=9)
for(i in 1:9){for(j in 1:9){
if(i > j){
	gM[i,j]<-cor(psCat[,i],psCat[,j])
}
if(i < j){
	gM[i,j]<-cor(psLm[,i],psLm[,j])
}
}}
gMlm<-matrix(NA,nrow=9,ncol=9)
gMms<-matrix(NA,nrow=9,ncol=9)
for(i in 1:9){for(j in 1:9){
if(i != j){
	gMms[i,j]<-cor(psCat[,i],psCat[,j])
	gMlm[i,j]<-cor(psLm[,i],psLm[,j])
}
}}
cor.test(gMms[upper.tri(gMms)],gMlm[upper.tri(gMlm)])

#	Pearson's product-moment correlation

#data:  gMms[upper.tri(gMms)] and gMlm[upper.tri(gMlm)]
#t = 7.7001, df = 34, p-value = 5.923e-09
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.6349918 0.8920979
#sample estimates:
#      cor 
#0.7972146 


save(list=ls(),file="plots.rdat")
### PLOTS ####
cl<-1.55;ca<-1.1;cm<-1.6
pdf("fig_genetics.pdf",width=9*1.05,height=10*1.05)
layout(matrix(c(1:2,1,3,4,4,5,5),nrow=4,ncol=2,byrow=TRUE),widths=c(5,5),heights=c(5,1,3,3))
par(mar=c(4.5,4.5,2.5,2.5))
## panel A
cs<-brewer.pal(n=3,"Dark2")

rownames(pve)<-traits
dotchart(t(pve),labels=c("Ms","Lm","Ms+Lm"),cex.lab=cl,col=cs,pch=NA,xlim=c(0,.7),xlab="Variance explained",ylab="")

i<-43
for(j in 1:9){
	for(k in 1:3){
		a<-4-k
		segments(lb[j,a],i,ub[j,a],i)
		points(pve[j,a],i,pch=19,col=cs[a])
		i<-i-1
	}
	i<-i-2
}
title(main="(a) Caterpillar performance heritability",cex.main=cm)

gcs<-rev(brewer.pal(n=8,"RdGy"))
image(gM,breaks=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),col=gcs,axes=FALSE)
xa<-seq(from=0,to=9,length.out=9)/9
axis(1,at=xa,traits,las=2)
axis(2,at=xa,traits,las=2)
box()
title(main="(b) Genetic correlations",cex.main=cm)

par(mar=c(2.2,5.5,2.2,5.5))
z<-matrix(c(-0.8,-.6,-.3,-0.1,.1,.3,.6,.8),nrow=1)
image(t(z),axes=FALSE,breaks=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),col=gcs)
axis(1,at=seq(from=-0.071,to=1.071,length.out=9),c(-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0))
box()
par(mar=c(4.5,4.5,2.5,2.5))
title(main="correlation")

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

title(main="(c) M. sativa genotype-performance associations",cex.main=cm)

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
title(main="(d) L. melissa genotype-performance associations",cex.main=cm)
abline(h=c(.1,.5),lty=2)

dev.off()


