## plant and caterpillar population effects
## needs R 4.0 or >
## using 4.0.2
library("lme4")
library("RLRsim")
library("RColorBrewer")

cat<-read.csv("caterpillars_cleaned_work.csv")

## plot
tab<-rbind(tapply(X=cat$Binary_8day,INDEX=cat$Cat_Pop,mean),
tapply(X=cat$Binary_14day,INDEX=cat$Cat_Pop,mean),
tapply(X=cat$Binary_Pupation,INDEX=cat$Cat_Pop,mean),
tapply(X=cat$Binary_Eclosion,INDEX=cat$Cat_Pop,mean))
tab<-tab[,c(1:2,4:5,3,6)]

pdf("sfig_surv2018.pdf",width=7,height=5)
par(mar=c(5,5,1,1))
barplot(tab,beside=TRUE,names.arg=c("Lm-BST","Lm-BWP","Lm-HWR","Lm-JJT","Ce","Vc"),ylab="Survival proportion",xlab="Species or population",cex.lab=1.5,col=brewer.pal(n=4,"BuGn"),ylim=c(0,1))
legend(26,.99,c("8d","14d","Pup.","Eclos."),fill=brewer.pal(n=4,"BuGn"),bty='n')
box()
dev.off()

mat8d<-tapply(cat$X8_day_weight_.mg.,INDEX=list(cat$Plant,cat$Cat_Pop),mean,na.rm=TRUE)
mat14d<-tapply(cat$X14_day_weight_.mg.,INDEX=list(cat$Plant,cat$Cat_Pop),mean,na.rm=TRUE)
mat8d<-mat8d[,c(1:2,4:5,3,6)]
mat14d<-mat14d[,c(1:2,4:5,3,6)]

gcs<-rev(brewer.pal(n=8,"RdGy"))
xa<-seq(0,6,length.out=6)/6
nms<-c("Lm-BST","Lm-BWP","Lm-HWR","Lm-JJT","Ce","Vc")
cm<-1.5

pdf("sfig_cors2018.pdf",width=7,height=12)
par(mfrow=c(2,1))
image.plot(cor(mat8d),breaks=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),col=gcs,axes=FALSE)
axis(1,xa,nms,las=2)
axis(2,xa,nms,las=2)
box()
title(main="(a) 8 day weight",cex.main=cm)

image(cor(mat14d),breaks=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),col=gcs,axes=FALSE)
axis(1,xa,nms,las=2)
axis(2,xa,nms,las=2)
box()
title(main="(b) 14 day weight",cex.main=cm)
dev.off()

### variance components
melissa<-c("BST","BWP","HWR","JJT")
Lm<-which(cat$Cat_Pop %in% melissa) ## 672 melissa, 133 colias and 196 painted ladies
mdat<-cat[Lm,]

## hatch date, turn to numeric and center
hd<-as.Date(mdat$Hatch_date,format="%m/%d/%y")
hdn<-as.numeric(hd-mean(hd))

w8d<-as.numeric(mdat$X8_day_weight_.mg.)
w14d<-as.numeric(mdat$X14_day_weight_.mg.)
wPup<-as.numeric(mdat$Pupal_weight_.mg.) ## only 16 non-NA? too few pupated to really use this

mean(mdat$Binary_8day)## most made it to 8 and 14 days
#[1] 0.9136905
mean(mdat$Binary_14day)
#[1] 0.8258929
mean(mdat$Binary_Pupation) ## but very few pupated
#[1] 0.02678571
mean(mdat$Binary_Eclosion)
#[1] 0.005952381

df<-data.frame(hdn=hdn,w8d=w8d,w14d=w14d,ppop=mdat$Plant,cpop=mdat$Cat_Pop)

## 8dw
est<-lmer(w8d ~ 1 + hdn + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.04207419
out<-exactRLRT(est)
#RLRT = 13.556, p-value < 2.2e-16
estP<-est

est<-lmer(w8d ~ 1 + hdn + (1|cpop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
# 0.0797906
out<-exactRLRT(est)
#RLRT = 15.595, p-value < 2.2e-16
estC<-est

est<-lmer(w8d ~ 1 + hdn + (1|cpop) + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)## cat
#[1] 0.03987918
vc$vcov[2]/sum(vc$vcov)## plant
#[1] 0.07956096
out<-exactRLRT(mA=est,m=estC,m0=estP)## test cat
#RLRT = 16.128, p-value < 2.2e-16
out<-exactRLRT(mA=est,m=estP,m0=estC)## test plant
#RLRT = 14.089, p-value < 1e-04

## 14dw
est<-lmer(w14d ~ 1 + hdn + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.1408897
out<-exactRLRT(est)
#RLRT = 56.88, p-value < 2.2e-16
estP<-est

est<-lmer(w14d ~ 1 + hdn + (1|cpop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#boundary (singular) fit: see ?isSingular
out<-exactRLRT(est)
#NA

est<-lmer(w14d ~ 1 + hdn + (1|cpop) + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
#boundary (singular) fit: see ?isSingular

## ended up with singular fits for caterpillar family as well

## colias
Co<-which(cat$Cat_Pop == "COL") ## 133 colias
cdat<-cat[Co,]

## hatch date, turn to numeric and center
hd<-as.Date(cdat$Hatch_date,format="%m/%d/%y")
hdn<-as.numeric(hd-mean(hd))

w8d<-as.numeric(cdat$X8_day_weight_.mg.)
w14d<-as.numeric(cdat$X14_day_weight_.mg.)
wPup<-as.numeric(cdat$Pupal_weight_.mg.) 

mean(cdat$Binary_8day)## most at least  made it pupation
#[1] 0.924812
mean(cdat$Binary_14day)
#[1] 0.8796992
mean(cdat$Binary_Pupation) 
#[1] 0.8120301
mean(cdat$Binary_Eclosion)
#[1] 0.6992481
df<-data.frame(hdn=hdn,w8d=w8d,w14d=w14d,wPup=wPup,ppop=cdat$Plant)

## 8dw
est<-lmer(w8d ~ 1 + hdn + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.02969577
out<-exactRLRT(est)
#RLRT = 0.70178, p-value = 0.1637

## 14dw
est<-lmer(w14d ~ 1 + hdn + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.09398847
out<-exactRLRT(est)
#RLRT = 4.2894, p-value = 0.0141

## wPup
est<-lmer(wPup ~ 1 + hdn + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.05287054
out<-exactRLRT(est)
#RLRT = 1.251, p-value = 0.0969

## painted lady 
Vc<-which(cat$Cat_Pop == "PLY") ## 196 V. cardui
vdat<-cat[Vc,]

## hatch date, turn to numeric and center
hd<-as.Date(vdat$Hatch_date,format="%m/%d/%y")
hdn<-as.numeric(hd-mean(hd)) ## no real variation

w8d<-as.numeric(vdat$X8_day_weight_.mg.)
w14d<-as.numeric(vdat$X14_day_weight_.mg.)
wPup<-as.numeric(vdat$Pupal_weight_.mg.) 

mean(vdat$Binary_8day)## about half made it to 8 or 14, none to pupation
#0.4489796
mean(vdat$Binary_14day)
#0.4030612
mean(vdat$Binary_Pupation) 
#0
mean(vdat$Binary_Eclosion)
#0
df<-data.frame(w8d=w8d,w14d=w14d,ppop=vdat$Plant)

## 8dw
est<-lmer(w8d ~ 1 + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.1008093
out<-exactRLRT(est)
#RLRT = 3.4369, p-value = 0.0227

## 14dw
est<-lmer(w14d ~ 1 + (1|ppop),data=df,verbose=1,REML=TRUE)
vc<-as.data.frame(VarCorr(est))
vc$vcov[1]/sum(vc$vcov)
#0.1297467
out<-exactRLRT(est)
#RLRT = 4.4349, p-value = 0.0119
