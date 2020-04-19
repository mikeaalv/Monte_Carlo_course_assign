rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(ggplot2)
require(xml2)
require(cowplot)
require(pracma)
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/Monte_Carlo_course_assign/assign6/result/"
setwd(dir)
fondsize=10
L=16
##wanglandau sampling for estimating E and Cv
locfiles=list.files(pattern="output\\.wanglandau\\.rand.*\\.gener\\.tab")
tablist=vector(mode="list");
for(file in locfiles){
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> rand
  datatab=read.table(file)
  colnames(datatab)=c("E","glog","H_last")
  tablist[[paste(rand)]]=datatab
}
Tseq=c(2.3,5.0)
taball=as.data.frame(matrix(NA,nrow=3*2,ncol=4))
colnames(taball)=c("T","rand","E","Cv")
rowi=1
for(Tele in Tseq){
  for(rand in names(tablist)){
    taball[rowi,"T"]=Tele
    taball[rowi,"rand"]=rand
    loctab=tablist[[rand]]
    glog_cent=loctab[,"glog"]
    glog_cent=glog_cent-mean(glog_cent)
    Zvec=exp(glog_cent)*exp(-loctab[,"E"]/Tele)
    Prob=Zvec/sum(Zvec)
    taball[rowi,"E"]=sum(loctab[,"E"]*Prob)/(L^2)
    taball[rowi,"Cv"]=(sum(loctab[,"E"]^2*Prob)-sum(loctab[,"E"]*Prob)^2)/(Tele^2*L^2)
    rowi=rowi+1
  }
}
##plotting E
temptab=taball[,c("T","rand","E")]
meantab=aggregate(E~T,data=temptab,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(E~T,data=temptab,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp1=cbind(sdtab,rep("metropolis",times=dim(sdtab)[1]))
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   labs(title=expression(paste("E vs T",)),x="T",y="E")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"E_vs_T_wanglandau.pdf"))

##plotting Cv
temptab=taball[,c("T","rand","Cv")]
meantab=aggregate(Cv~T,data=temptab,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv~T,data=temptab,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp1=cbind(sdtab,rep("metropolis",times=dim(sdtab)[1]))
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_wanglandau.pdf"))

write.table(taball,file="wanglandau_est.tab",sep="\t")

##plot one density of state
rand="1"
temptab=tablist[[rand]]
p<-ggplot(data=temptab,aes(x=E,y=glog))+
   geom_line()+
   labs(title=expression(paste("log density vs E",)),x="E",y="log density")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"logden_vs_E_wanglandau.pdf"))

##plot histogram
rand="1"
temptab=tablist[[rand]]
p<-ggplot(data=temptab,aes(x=E,y=H_last))+
   geom_line()+
   labs(title=expression(paste("histogram vs E",)),x="E",y="histogram")+
   ylim(c(0,100000))+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"histogram_vs_E_wanglandau.pdf"))

##metropolis estimation
stepsize=5
samplesize=10000
heatup=100
locfiles=list.files(pattern="output\\.metropois\\.T\\d+\\.*\\d*\\.rand[123]*\\.gener\\.tab")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=4))
colnames(tabsumm)=c("T","rand","E","Cv")
sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.T\\d+\\.*\\d*+") %>%
           str_replace_all(.,pattern="\\.T",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"T"]
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"rand"]
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  tabsumm[rowi,"E"]=mean(datatab[sampleseq,"H"])/L^2
  tabsumm[rowi,"Cv"]=var(datatab[sampleseq,"H"])/(tabsumm[rowi,"T"]^2*L^2)
  rowi=rowi+1
}
##plotting E
temptab=tabsumm[,c("T","rand","E")]
meantab=aggregate(E~T,data=temptab,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(E~T,data=temptab,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp1=cbind(sdtab,rep("metropolis",times=dim(sdtab)[1]))
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   labs(title=expression(paste("E vs T",)),x="T",y="E")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"E_vs_T_metropolis.pdf"))

##plotting Cv
temptab=tabsumm[,c("T","rand","Cv")]
meantab=aggregate(Cv~T,data=temptab,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv~T,data=temptab,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp1=cbind(sdtab,rep("metropolis",times=dim(sdtab)[1]))
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_metropolis.pdf"))

write.table(tabsumm,file="metropolis_est.tab",sep="\t")
