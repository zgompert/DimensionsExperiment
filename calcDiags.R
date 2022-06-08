## gelman-rubin diagonstic for main gemma runs
library(coda)
library(RColorBrewer)
library(scales)

## read list of trait files, without chain numbers
ff<-read.table("diagFiles.txt",header=FALSE)

N<-dim(ff)[1]

pve<-matrix(NA,nrow=N,ncol=2)
pge<-matrix(NA,nrow=N,ncol=2)
ess<-matrix(NA,nrow=N,ncol=2)

for(k in 1:N){
	cat(ff[k,],"\n")

	clist<-vector("list",10)
	md<-100000
	for(i in 1:10){
		dat<-read.table(paste(ff[k,],i-1,".hyp.txt",sep=""),header=TRUE)
	        if(dim(dat)[1] < md){md<-dim(dat)[1]}
		clist[[i]]<-mcmc(dat[,c(2,4)])
	}
    	if(md<100000){
		for(i in 1:10){
			clist[[i]]<-as.mcmc(clist[[i]][1:md,])
		}
	}
	mlist<-as.mcmc.list(clist)
	o<-gelman.diag(mlist)
	pve[k,]<-o$psrf[1,]
	pge[k,]<-o$psrf[2,]
	o<-effectiveSize(mlist)
	ess[k,]<-o
}

cl<-1.5;ca<-1.1;cm<-1.5
png("sfig_mixing.png",width=9,height=8,units="in",res=480)
par(mfrow=c(2,1))
par(mar=c(4,5,2,1))
cs<-alpha(brewer.pal("BrBG",n=10),.5)
k<-1
plot(1:100000,rep(0,100000),type='n',ylim=c(0,1),ylab="PVE",xlab="MCMC sample",cex.lab=cl,cex.axis=ca)
for(i in 1:10){
	dat<-read.table(paste(ff[k,],i-1,".hyp.txt",sep=""),header=TRUE)
	lines(dat[,2],col=cs[i])
}
title(main=substitute(paste("(a) Trace ",italic("M. sativa")," W8d",sep="")),cex.main=cm)
k<-13
plot(1:100000,rep(0,100000),type='n',ylim=c(0,1),ylab="PVE",xlab="MCMC sample",cex.lab=cl,cex.axis=ca)
for(i in 1:10){
	dat<-read.table(paste(ff[k,],i-1,".hyp.txt",sep=""),header=TRUE)
	lines(dat[,2],col=cs[i])
}
title(main=substitute(paste("(b) Trace ",italic("L. melissa")," S8d",sep="")),cex.main=cm)
dev.off()


cl<-1.5;ca<-1.1;cm<-1.5
pdf("sfig_mixingSummary.pdf",width=9,height=8)
par(mfrow=c(2,2))
par(mar=c(4,5,2,1))
cs<-rep(brewer.pal(n=3,"Dark2"),each=9)

plot(1:27,pve[1:27,1],pch=19,col=cs,ylim=c(1,1.3),xlab="Performance traits",ylab="Scale reduction factor",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
segments(1:27,pve[1:27,1],1:27,pve[1:27,2],col=cs)
box()
title(main="(a) Diagnostic PVE",cex.main=cm)
legend(1,1.3,c("M. sativa","L. melissa","Combined"),fill=unique(cs),bty='n',cex=ca*1.3)

plot(1:27,pge[1:27,1],pch=19,col=cs,ylim=c(1,1.3),xlab="Performance traits",ylab="Scale reduction factor",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
segments(1:27,pge[1:27,1],1:27,pge[1:27,2],col=cs)
box()
title(main="(b) Diagnostic PGE",cex.main=cm)

plot(1:27,log10(ess[1:27,1]),pch=19,col=cs,ylim=c(0,6),xlab="Performance traits",ylab="Effective sample size (log10)",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
box()
title(main="(c) Effective sample size PVE",cex.main=cm)

plot(1:27,log10(ess[1:27,2]),pch=19,col=cs,ylim=c(0,6),xlab="Performance traits",ylab="Effective sample size (log10)",cex.lab=cl,axes=FALSE)
axis(2,cex.axis=ca)
box()
title(main="(d) Effective sample size PGE",cex.main=cm)

dev.off()



