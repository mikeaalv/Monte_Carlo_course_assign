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
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/Monte_Carlo_course_assign/assign5/result/"
setwd(dir)
fondsize=10
L=8
##histogram based on infinite temperature
heatup=100
stepsize=1
samplesize=5000000
locfiles=list.files(pattern="output\\.simplesampling\\.hist\\.rand.*\\.gener\\.tab")
taball=c()
sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> rand
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  hvec=datatab[sampleseq,"H"]
  taball=rbind(taball,cbind(rep(rand,times=length(hvec)),hvec))
  rowi=rowi+1
}
colnames(taball)=c("rand","H")
Trange=seq(from=2.0,to=4,by=0.1)
randseq=unique(taball[,"rand"])
Cvec_col_list=vector(mode="list")
plist=vector(mode="list")
for(randi in randseq){
  loc_histogram=table(taball[taball[,"rand"]==randi,"H"])
  loc_hami=as.numeric(names(loc_histogram))
  ## histogram to calculate Cv
  Cvvec=c()
  for(Tele in Trange){
    delK= -1/Tele;
    Pele=loc_histogram*exp(delK*loc_hami)
    P=Pele/sum(Pele)
    Cv=(sum(P*loc_hami^2)-sum(P*loc_hami)^2)/(L^2*Tele^2)
    Cvvec=c(Cvvec,Cv)
    if(Tele==2||Tele==4){
      plist[[paste0(Tele)]]=P
    }
  }
  Cvec_col_list[[randi]]=cbind(cbind(Trange,rep(randi,times=length(Trange))),Cvvec)
}
Cvectab_col=as.data.frame(Reduce('rbind',Cvec_col_list))
colnames(Cvectab_col)=c("T","randi","Cv")
meantab=aggregate(Cv~T,data=Cvectab_col,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv~T,data=Cvectab_col,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_ribbon(aes(ymin=mean-2*sd,ymax=mean+2*sd),fill="grey70",colour="grey70")+
   geom_line(size=1)+
   geom_point()+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_histogram_Tinfi.pdf"))
summarytab[summarytab[,"mean"]==max(summarytab[,"mean"]),"T"]
##plot P
pdf(paste0(dir,"P_histogram_Tinfi_2.pdf"))
plot(plist[["2"]],xlab="H",ylab="P")
dev.off()
pdf(paste0(dir,"P_histogram_Tinfi_4.pdf"))
plot(plist[["4"]],xlab="H",ylab="P")
dev.off()

##Cv from metropois sampling
stepsize=5
samplesize=10000
heatup=100
locfiles=list.files(pattern="output\\.metropois\\.T\\d+\\.*\\d*\\.rand[123]*\\.gener\\.tab")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=5))
colnames(tabsumm)=c("T","rand","H","H_S","Cv")
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
  tabsumm[rowi,"H"]=mean(datatab[sampleseq,"H"])
  tabsumm[rowi,"H_S"]=tabsumm[rowi,"H"]/L^2
  tabsumm[rowi,"Cv"]=var(datatab[sampleseq,"H"])/(tabsumm[rowi,"T"]^2*L^2)
  rowi=rowi+1
}
##plotting
temptab=tabsumm[,c("T","rand","Cv")]
meantab=aggregate(Cv~T,data=temptab,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv~T,data=temptab,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp1=cbind(sdtab,rep("metropolis",times=dim(sdtab)[1]))
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_metropolis.pdf"))

##histogram based on Tc
heatup=100
stepsize=1
samplesize=300000
locfiles=list.files(pattern="output\\.metropois\\.Tc\\.hist\\.rand[123]+.gener\\.tab")
taballtc=c()
sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
Tc=3.0
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> rand
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  hvec=datatab[sampleseq,"H"]
  taballtc=rbind(taballtc,cbind(rep(rand,times=length(hvec)),hvec))
  rowi=rowi+1
}
colnames(taballtc)=c("rand","H")
Trange=seq(from=2.0,to=4,by=0.1)
randseq=unique(taballtc[,"rand"])
Cvec_col_tc_list=vector(mode="list")
for(randi in randseq){
  loc_histogram=table(taballtc[taballtc[,"rand"]==randi,"H"])
  loc_hami=as.numeric(names(loc_histogram))
  ## histogram to calculate Cv
  Cvvec=c()
  for(Tele in Trange){
    delK= -1/Tele+1/Tc;
    Pele=loc_histogram*exp(delK*loc_hami)
    P=Pele/sum(Pele)
    Cv=(sum(P*loc_hami^2)-sum(P*loc_hami)^2)/(L^2*Tele^2)
    Cvvec=c(Cvvec,Cv)
  }
  Cvec_col_tc_list[[randi]]=cbind(cbind(Trange,rep(randi,times=length(Trange))),Cvvec)
}
Cvectab_tc_col=as.data.frame(Reduce('rbind',Cvec_col_tc_list))
colnames(Cvectab_tc_col)=c("T","randi","Cv")
meantab=aggregate(Cv~T,data=Cvectab_tc_col,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv~T,data=Cvectab_tc_col,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp2=cbind(sdtab,rep("hist",times=dim(sdtab)[1]))
p<-ggplot(data=summarytab,aes(x=T,y=mean))+
   geom_ribbon(aes(ymin=mean-2*sd,ymax=mean+2*sd),fill="grey70",colour="grey70")+
   geom_line(size=1)+
   geom_point()+
   labs(title=expression(paste("Cv vs T",)),x="T",y="Cv")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cv_vs_T_histogram_Tc.pdf"))
summarytab[summarytab[,"mean"]==max(summarytab[,"mean"]),"T"]


##plot sd for metropolis and histogram method
##metropolis
stepsize=5
samplesize=10000
heatup=100
locfiles=list.files(pattern="output\\.metropois\\.T\\d+\\.*\\d*\\.rand\\d*\\.gener\\.tab")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=5))
colnames(tabsumm)=c("T","rand","H","H_S","Cv")
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
  tabsumm[rowi,"H"]=mean(datatab[sampleseq,"H"])
  tabsumm[rowi,"H_S"]=tabsumm[rowi,"H"]/L^2
  tabsumm[rowi,"Cv"]=var(datatab[sampleseq,"H"])/(tabsumm[rowi,"T"]^2*L^2)
  rowi=rowi+1
}
temptab=tabsumm[,c("T","rand","Cv")]
tabsd_summ=as.data.frame(matrix(NA,nrow=dim(temptab)[1]/3,ncol=3))
colnames(tabsd_summ)=c("T","group","Cv_sd")
irow=1
for(Tele in unique(temptab[,"T"])){
  for(randgroup in c(0,1,2)){
    subtab=temptab[temptab[,"T"]==Tele&temptab[,"rand"]%%3==randgroup,]
    tabsd_summ[irow,"T"]=Tele
    tabsd_summ[irow,"Cv_sd"]=sd(subtab[,"Cv"])
    tabsd_summ[irow,"group"]=randgroup
    irow=irow+1
  }
}
meantab=aggregate(Cv_sd~T,data=tabsd_summ,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv_sd~T,data=tabsd_summ,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp1=cbind(summarytab,rep("metropolis",times=dim(summarytab)[1]))
##histogram
heatup=100
stepsize=1
samplesize=300000
locfiles=list.files(pattern="output\\.metropois\\.Tc\\.hist\\.rand\\d+.gener\\.tab")
taballtc=c()
sampleseq=seq(from=heatup,by=stepsize,length.out=samplesize)
Tc=3.0
rowi=1
for(file in locfiles){
  # print(file)
  file %>% str_extract_all(.,pattern="\\.rand\\d+") %>%
           str_replace_all(.,pattern="\\.rand",replacement="") %>%
           as.numeric(.) -> rand
  datatab=read.table(file)
  colnames(datatab)=c("H","M")
  hvec=datatab[sampleseq,"H"]
  taballtc=rbind(taballtc,cbind(rep(rand,times=length(hvec)),hvec))
  rowi=rowi+1
}
colnames(taballtc)=c("rand","H")
Trange=seq(from=2.0,to=4,by=0.1)
randseq=unique(taballtc[,"rand"])
Cvec_col_tc_list=vector(mode="list")
for(randi in randseq){
  loc_histogram=table(taballtc[taballtc[,"rand"]==randi,"H"])
  loc_hami=as.numeric(names(loc_histogram))
  ## histogram to calculate Cv
  Cvvec=c()
  for(Tele in Trange){
    delK= -1/Tele+1/Tc;
    Pele=loc_histogram*exp(delK*loc_hami)
    P=Pele/sum(Pele)
    Cv=(sum(P*loc_hami^2)-sum(P*loc_hami)^2)/(L^2*Tele^2)
    Cvvec=c(Cvvec,Cv)
  }
  Cvec_col_tc_list[[randi]]=cbind(cbind(Trange,rep(randi,times=length(Trange))),Cvvec)
}
Cvectab_tc_col=as.data.frame(Reduce('rbind',Cvec_col_tc_list))
colnames(Cvectab_tc_col)=c("T","randi","Cv")
tabsd_summ=as.data.frame(matrix(NA,nrow=dim(Cvectab_tc_col)[1]/3,ncol=3))
colnames(tabsd_summ)=c("T","group","Cv_sd")
irow=1
for(Tele in unique(Cvectab_tc_col[,"T"])){
  for(randgroup in c(0,1,2)){
    subtab=Cvectab_tc_col[Cvectab_tc_col[,"T"]==Tele&Cvectab_tc_col[,"randi"]%%3==randgroup,]
    tabsd_summ[irow,"T"]=Tele
    tabsd_summ[irow,"Cv_sd"]=sd(subtab[,"Cv"])
    tabsd_summ[irow,"group"]=randgroup
    irow=irow+1
  }
}
meantab=aggregate(Cv_sd~T,data=tabsd_summ,mean)
colnames(meantab)=c("T","mean")
sdtab=aggregate(Cv_sd~T,data=tabsd_summ,sd)
colnames(sdtab)=c("T","sd")
summarytab=merge(meantab,sdtab,by=c("T"),all=TRUE)
sdtab_comp2=cbind(summarytab,rep("hist",times=dim(summarytab)[1]))

colnames(sdtab_comp1)=c("T","mean","sd","meth")
colnames(sdtab_comp2)=c("T","mean","sd","meth")
sdtab_comp=rbind(sdtab_comp1,sdtab_comp2)
# dplyr::bind_rows(sdtab_comp,sdtab_comp2)
sdtab_comp[,"meth"]=as.factor(sdtab_comp[,"meth"])
p<-ggplot(data=sdtab_comp,aes(x=T,y=mean,color=meth))+
   geom_line(size=1)+
   geom_point()+
   geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
   labs(title=expression(paste("Cv_sd vs T",)),x="T",y="Cv_sd")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"Cvsd_vs_T_histogram_Tc_vs_metropolis2.pdf"))
