## more data wrangling, compute actual response variables from raw data and fix anything that is still screwy

dat<-read.csv("../Scratch/captive-butterflies-clean-2019.csv")
## all plant and cat IDs are unique
## butterfly DNA ids for 3 appear twice:
## Lyc-048-F11 Lyc-048-H12  Lyc-056-B4
## will either drop or figure out what is what

wgt8day<-(as.numeric(as.character(dat$X8_day_weight_.mg._1)) + as.numeric(as.character(dat$X8_day_weight_.mg._2)))/2

wgt14day<-(as.numeric(as.character(dat$X14_day_weight_1.mg)) + as.numeric(as.character(dat$X14_day_weight_2.mg)))/2

wgtPup<-(as.numeric(as.character(dat$Pupal_weight_1.mg)) + as.numeric(as.character(dat$Pupal_weight_2.mg)))/2


s8day<-as.numeric(is.na(dat$X8_day_weight_.mg._1)==FALSE)
## 88% made it to day 8
s14day<-as.numeric(is.na(dat$X14_day_weight_1.mg.)==FALSE)
## 82% made it to day 14
sPupa<-as.numeric(is.na(dat$Pupation_date)==FALSE)
## 25% made it to pupation
sEclose<-as.numeric(is.na(dat$Eclosion_date)==FALSE)
## 12% made it to adult

hatchDay<-as.Date(dat$Hatch_Date,format="%m/%d")
deathDay<-as.Date(dat$Death_date,format="%m/%d")
stime<-as.numeric(deathDay-hatchDay)
stimeTrunc<-stime
stimeTrunc[sEclose==1]<-max(stime,na.rm=TRUE)+1

day<-as.numeric(as.Date(dat$Hatch_Date,format="%m/%d")-min(as.Date(dat$Hatch_Date,format="%m/%d")))
pops<-gsub(dat$ind_cat,pattern="[0-9-]+",perl=TRUE,replacement="")

## almost certainly need to control for hatch date, it seems to really matter
## and population matters some (maybe confounded with hatch date... FRP seems to really
## rock, how does controlling work here, need to think, really depends on whether
## we are capturing maternal effects or polygenic variation, probably worth doing both ways
## and it looks like really might be hatch date more than population per se

df<-data.frame(Plant_ID=dat$Plant_ID,Cat_ID=dat$ind_cat,DNA=dat$Extraction.Plate,pop=pops,hatchD=day,w8day=wgt8day,w14day=wgt14day,wPupa=wgtPup,s8day=s8day,s14day=s14day,sPupa=sPupa,sEclose=sEclose,stime=stime,stimeTrunc=stimeTrunc)

write.table(df,file="dimensions_lmel_dat.csv",quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")
