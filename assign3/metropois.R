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
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/Monte_Carlo_course_assign/assign3/result/"
setwd(dir)
## plot test run data
# ### check converged(plateau)
# endregion=5
# tab=read.table("output.L20.T2.4.rand1.tab")
# len=dim(tab)[1]
# par(mfrow=c(2,2))
# plot(tab[,1],type='l')
# plot(tab[,2],type='l')
# for(i=1:2){
#   Mend=mean(tab[(len-endregion):len,i])
#   phi=(tab[1:(len-endregion),i]-Mend)/(tab[1,i]-Mend)
#   plot(phi,type='l')
# }
# ### simple check of relaxation
# acf(tab[50:1000,1])
# acf(tab[50:1000,2])
fondsize=10
## analysis of acommulated data
### general run of different L and T
stepsize=5
samplesize=1000
heatup=100
locfiles=list.files(pattern="gener\\.tab$")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=9))
colnames(tabsumm)=c("L","T","rand","H","H_S","M","M_S","Cv","sus")
sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.L\\d+") %>%
           str_replace_all(.,pattern="\\.L",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"L"]
  file %>% str_extract_all(.,pattern="\\.T\\d+\\.*\\d*+") %>%
           str_replace_all(.,pattern="\\.T",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"T"]
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"rand"]
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  tabsumm[rowi,"H"]=mean(datatab[sampleseq,"H"])
  tabsumm[rowi,"M"]=mean(datatab[sampleseq,"M"])
  tabsumm[rowi,"H_S"]=tabsumm[rowi,"H"]/tabsumm[rowi,"L"]^2
  tabsumm[rowi,"M_S"]=tabsumm[rowi,"M"]/tabsumm[rowi,"L"]^2
  tabsumm[rowi,"Cv"]=var(datatab[sampleseq,"H"])/tabsumm[rowi,"T"]^2
  tabsumm[rowi,"sus"]=var(datatab[sampleseq,"M"])/tabsumm[rowi,"T"]
  rowi=rowi+1
}
tabsumm[,"L"]=as.factor(tabsumm[,"L"])
##H per spin plot
temptab=tabsumm[,c("L","T","rand","H_S")]
meantab=aggregate(H_S~L+T,data=temptab,mean)
colnames(meantab)=c("L","T","mean")
sdtab=aggregate(H_S~L+T,data=temptab,sd)
colnames(sdtab)=c("L","T","sd")
summarytab=merge(meantab,sdtab,by=c("L","T"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=T,y=mean,color=L,group=L))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("U vs T",)),x="T",y="U(per spin)")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"U_vs_T.pdf"))

##M per spin plot
temptab=tabsumm[,c("L","T","rand","M_S")]
meantab=aggregate(M_S~L+T,data=temptab,mean)
colnames(meantab)=c("L","T","mean")
sdtab=aggregate(M_S~L+T,data=temptab,sd)
colnames(sdtab)=c("L","T","sd")
summarytab=merge(meantab,sdtab,by=c("L","T"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=T,y=mean,color=L,group=L))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("M vs T",)),x="T",y="M(per spin)")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"M_vs_T.pdf"))

##Cv plot
temptab=tabsumm[,c("L","T","rand","Cv")]
meantab=aggregate(Cv~L+T,data=temptab,mean)
colnames(meantab)=c("L","T","mean")
sdtab=aggregate(Cv~L+T,data=temptab,sd)
colnames(sdtab)=c("L","T","sd")
summarytab=merge(meantab,sdtab,by=c("L","T"),all=TRUE)
Cvesti_var=summarytab
p<-ggplot(data=summarytab,aes(x=T,y=mean,color=L,group=L))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T.pdf"))

##susceptibility plot
temptab=tabsumm[,c("L","T","rand","sus")]
meantab=aggregate(sus~L+T,data=temptab,mean)
colnames(meantab)=c("L","T","mean")
sdtab=aggregate(sus~L+T,data=temptab,sd)
colnames(sdtab)=c("L","T","sd")
summarytab=merge(meantab,sdtab,by=c("L","T"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=T,y=mean,color=L,group=L))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("susceptibility vs T",)),x="T",y="susceptibility")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"sus_vs_T.pdf"))

##numeric differential for Cv
# calc_point=c(2.5,19.5)
temptab=tabsumm[,c("L","T","rand","H")]
Cvest_der=as.data.frame(matrix(NA,nrow=length(unique(temptab[,"L"]))*2,ncol=4))
colnames(Cvest_der)=c("L","T","mean","sd")
rowi=1
for(Lele in unique(temptab[,"L"])){
  cvvec=list("2.5"=c(),"19.5"=c())
  for(randele in unique(temptab[,"rand"])){
    loctab=temptab[temptab[,"L"]==Lele&temptab[,"rand"]==randele,]
    rownames(loctab)=loctab[,"T"]
    cvvec[["2.5"]]=c(cvvec[["2.5"]],(loctab["3","H"]-loctab["2","H"])/1)
    cvvec[["19.5"]]=c(cvvec[["19.5"]],(loctab["20","H"]-loctab["19","H"])/1)
  }
  Cvest_der[rowi,"L"]=Lele
  Cvest_der[rowi,"T"]=2.5
  Cvest_der[rowi,"mean"]=mean(cvvec[["2.5"]])
  Cvest_der[rowi,"sd"]=sd(cvvec[["2.5"]])
  Cvest_der[rowi+1,"L"]=Lele
  Cvest_der[rowi+1,"T"]=19.5
  Cvest_der[rowi+1,"mean"]=mean(cvvec[["19.5"]])
  Cvest_der[rowi+1,"sd"]=sd(cvvec[["19.5"]])
  rowi=rowi+2
}
p<-ggplot(data=Cvest_der,aes(x=T,y=mean,color=L,group=L))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_deriv.pdf"))

Cvesti_var

## run with different length
stepsize=5
samplesizevec=c(100,10000)
heatup=100
locfiles=list.files(pattern="mccomp\\.tab$")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=10))
colnames(tabsumm)=c("L","T","rand","Nsample","H","H_S","M","M_S","Cv","sus")
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.L\\d+") %>%
           str_replace_all(.,pattern="\\.L",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"L"]
  file %>% str_extract_all(.,pattern="\\.T\\d+\\.*\\d*+") %>%
           str_replace_all(.,pattern="\\.T",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"T"]
  file %>% str_extract_all(.,pattern="\\.mc\\d+") %>%
           str_replace_all(.,pattern="\\.mc",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"Nsample"]
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"rand"]
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  if(tabsumm[rowi,"Nsample"]==50099){
    samplesize=10000
  }else if(tabsumm[rowi,"Nsample"]==599){
    samplesize=100
  }
  sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
  tabsumm[rowi,"H"]=mean(datatab[sampleseq,"H"])
  tabsumm[rowi,"M"]=mean(datatab[sampleseq,"M"])
  tabsumm[rowi,"H_S"]=tabsumm[rowi,"H"]/tabsumm[rowi,"L"]^2
  tabsumm[rowi,"M_S"]=tabsumm[rowi,"M"]/tabsumm[rowi,"L"]^2
  tabsumm[rowi,"Cv"]=var(datatab[sampleseq,"H"])/tabsumm[rowi,"T"]^2
  tabsumm[rowi,"sus"]=var(datatab[sampleseq,"M"])/tabsumm[rowi,"T"]
  rowi=rowi+1
}
tabsumm[,"Nsample"]=as.factor(tabsumm[,"Nsample"])
##H per spin plot
temptab=tabsumm[,c("Nsample","T","rand","H_S")]
meantab=aggregate(H_S~Nsample+T,data=temptab,mean)
colnames(meantab)=c("Nsample","T","mean")
sdtab=aggregate(H_S~Nsample+T,data=temptab,sd)
colnames(sdtab)=c("Nsample","T","sd")
summarytab=merge(meantab,sdtab,by=c("Nsample","T"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=T,y=mean,color=Nsample,group=Nsample))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~Nsample,nrow=3,scales="free_y")+
   labs(title=expression(paste("U vs T",)),x="T",y="U(per spin)")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"U_vs_T_mc.pdf"))
temptab=tabsumm[,c("Nsample","T","rand","Cv")]
meantab=aggregate(Cv~Nsample+T,data=temptab,mean)
colnames(meantab)=c("Nsample","T","mean")
sdtab=aggregate(Cv~Nsample+T,data=temptab,sd)
colnames(sdtab)=c("Nsample","T","sd")
summarytab=merge(meantab,sdtab,by=c("Nsample","T"),all=TRUE)
Cvesti_var=summarytab
p<-ggplot(data=summarytab,aes(x=T,y=mean,color=Nsample,group=Nsample))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   facet_wrap(~Nsample,nrow=3,scales="free_y")+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_mc.pdf"))

## use dU/dT to estimate Cv
stepsize=5
samplesize=1000
heatup=100
locfiles=list.files(pattern="Tc\\.tab$")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=9))
colnames(tabsumm)=c("L","T","rand","H","H_S","M","M_S","Cv","sus")
sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.L\\d+") %>%
           str_replace_all(.,pattern="\\.L",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"L"]
  file %>% str_extract_all(.,pattern="\\.T\\d+\\.*\\d*+") %>%
           str_replace_all(.,pattern="\\.T",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"T"]
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> tabsumm[rowi,"rand"]
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  tabsumm[rowi,"H"]=mean(datatab[sampleseq,"H"])
  tabsumm[rowi,"M"]=mean(datatab[sampleseq,"M"])
  tabsumm[rowi,"H_S"]=tabsumm[rowi,"H"]/tabsumm[rowi,"L"]^2
  tabsumm[rowi,"M_S"]=tabsumm[rowi,"M"]/tabsumm[rowi,"L"]^2
  tabsumm[rowi,"Cv"]=var(datatab[sampleseq,"H"])/tabsumm[rowi,"T"]^2
  tabsumm[rowi,"sus"]=var(datatab[sampleseq,"M"])/tabsumm[rowi,"T"]
  rowi=rowi+1
}
Llist=unique(tabsumm[,"L"])
randlist=unique(tabsumm[,"rand"])
tabtc=as.data.frame(matrix(NA,nrow=length(Llist)*length(randlist),ncol=3))
colnames(tabtc)=c("L","rand","Tc")
rowi=1
for(Lele in Llist){
  for(randele in randlist){
    loctab=tabsumm[tabsumm[,"L"]==Lele&tabsumm[,"rand"]==randele,]
    uvec=loctab[order(loctab[,"T"]),"H"]
    len=length(uvec)
    num_der=(uvec[3:len]-uvec[1:(len-2)])/0.4
    tabtc[rowi,"L"]=Lele
    tabtc[rowi,"rand"]=randele
    tabtc[rowi,"Tc"]=loctab[which(num_der==max(num_der))+1,"T"]
    rowi=rowi+1
  }
}
meantab=aggregate(Tc~L,data=tabtc,mean)
colnames(meantab)=c("L","mean")
sdtab=aggregate(Tc~L,data=tabtc,sd)
colnames(sdtab)=c("L","sd")
summarytab=merge(meantab,sdtab,by=c("L"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=L,y=mean))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   # facet_wrap(~Nsample,nrow=3,scales="free_y")+
   labs(title=expression(paste("Tc vs L",)),x="L",y="Tc")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Tc_vs_L.pdf"))
