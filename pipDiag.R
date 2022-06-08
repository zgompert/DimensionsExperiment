## how consistent are pips from subsets of runs
library(RColorBrewer)

ff<-list.files(pattern="sub")

N<-length(ff)
cors<-rep(NA,N)
top01<-rep(NA,N)
top10<-rep(NA,N)
for(i in 1:N){
	dat<-read.table(ff[i],header=FALSE)
	cors[i]<-cor(dat[,2],dat[,3])
	top01[i]<-sum(dat[,2] > 0.01 & dat[,3] > 0.01)/sum(dat[,2] > 0.01 | dat[,3] > 0.01)
	top10[i]<-sum(dat[,2] > 0.1 & dat[,3] > 0.1)/sum(dat[,2] > 0.1 | dat[,3] > 0.1)
}

summary(cors)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03296 0.30722 0.90020 0.66862 0.98824 0.99980 
summary(top10)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.6429  0.7143  0.8333  0.8547  1.0000  1.0000      10 

cl<-1.5;ca<-1.1;cm<-1.5
pdf("sfig_pipConsistency.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(4,5,2,1))
cs<-rep(brewer.pal(n=3,"Dark2"),each=9)
plot(cors,pch=19,col=cs,ylim=c(0,1),xlab="Performance traits",ylab="PIP correlation",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
box()
title(main="(a) Correlation in SNP PIPs",cex.main=cm)

plot(top10,pch=19,col=cs,ylim=c(0,1),xlab="Performance traits",ylab="Proportion PIP > 0.1",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
box()
legend(20,.35,c("M. sativa","L. melissa","Combined"),pch=19,col=unique(cs),bty='n',cex=ca*1.25)

title(main="(b) Consistency in high PIP SNPs",cex.main=cm)
dev.off()
