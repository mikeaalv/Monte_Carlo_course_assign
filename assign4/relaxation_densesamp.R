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
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/Monte_Carlo_course_assign/assign4/result/"
setwd(dir)
fondsize=10
## analysis of acommulated data
### general run of different L and T
##for relaxation sampling function: moving product calculation
mp_cal<-function(x,tseq,meanval,varval){
  dim(x)=c(tseq,length(x)%/%tseq)
  x=t(x)
  refx=x[,1]
  samp_mat=x*as.vector(refx)
  samp_mat=(samp_mat-meanval^2)/varval
  returnlist=list(mean=colMeans(samp_mat),sd=apply(samp_mat,2,sd))
  return(returnlist)
}
##general sampling
heatup=100
endind=600099
stepsize=1
##relaxation sampling
tsteppre=10
tstrep_alt1=200
tstrep_alt2=500
tstrep_alt3=3000
locfiles=list.files(pattern="gener\\.tab$")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=9))
colnames(tabsumm)=c("L","T","rand","H","H_var","M","M_var","H_relaxseq_ind","M_relaxseq_ind")
sampleseq_balk=seq(from=heatup,by=stepsize,to=endind)
rowi=1
corrlist=vector(mode="list")
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
  corr_vec_list=vector(mode="list")
  for(varele in c("H","M")){
    t0step=tsteppre#the distance between each t0
    tseq=tsteppre#the length of each sampling
    if(varele=="M"&&tabsumm[rowi,"T"]==2.3){
      t0step=tstrep_alt1
      tseq=tstrep_alt1
    }
    if(varele=="H"&&tabsumm[rowi,"T"]==2.3&&tabsumm[rowi,"L"]==32){
      t0step=tsteppre
      tseq=tstrep_alt2
    }
    if(varele=="M"&&tabsumm[rowi,"T"]==2.3&&tabsumm[rowi,"L"]==32){
      t0step=tsteppre
      tseq=tstrep_alt3
    }
    tottimeseq=0:(tseq-1)
    ngroup_t0=floor((endind-heatup-tseq)/t0step)
    sampleseq_relax=sampleseq_balk[rep(seq(from=1,by=1,to=tseq),times=ngroup_t0)+rep(t0step*seq(from=0,by=1,to=ngroup_t0-1),each=tseq)]
    tabsumm[rowi,varele]=mean(datatab[sampleseq_balk,varele])
    tabsumm[rowi,paste0(varele,"_var")]=var(datatab[sampleseq_balk,varele])
    # H_corr_infor=mp_cal(datatab[sampleseq_relax,"H"],tseq)
    # M_corr_infor=mp_cal(datatab[sampleseq_relax,"M"],tseq)
    # corr_vec_list[["H"]]=list(ind=tottimeseq,
    #                           mean=(H_corr_infor[["mean"]]-tabsumm[rowi,"H"]^2)/tabsumm[rowi,"H_var"],
    #                           sd=H_corr_infor[["sd"]]/tabsumm[rowi,"H_var"])
    # corr_vec_list[["M"]]=list(ind=tottimeseq,
    #                           mean=(M_corr_infor[["mean"]]-tabsumm[rowi,"M"]^2)/tabsumm[rowi,"M_var"],
    #                           sd=M_corr_infor[["sd"]]/tabsumm[rowi,"M_var"])
    corr_infor=mp_cal(datatab[sampleseq_relax,varele],tseq,tabsumm[rowi,varele],tabsumm[rowi,paste0(varele,"_var")])
    corr_vec_list[[varele]]=list(ind=tottimeseq,
                              mean=corr_infor[["mean"]],
                              sd=corr_infor[["sd"]])
    tabsumm[rowi,paste0(varele,"_relaxseq_ind")]=rowi
  }
  corrlist[[rowi]]=corr_vec_list
  rowi=rowi+1
}
##test plot for mean curve on relaxation
# L_plot=32#4 32
# T_plot=2.3#2.3 5.0
# rand_plot=1
for(T_plot in unique(tabsumm[,"T"])){
  for(L_plot in unique(tabsumm[,"L"])){
    inds=tabsumm[tabsumm[,"L"]==L_plot&tabsumm[,"T"]==T_plot,"H_relaxseq_ind"]#&tabsumm[,"rand"]==rand_plot
    for(varele in c("H","M")){
      plotlist=vector(mode="list")
      for(ind in inds){
        storelist=corrlist[[ind]][[varele]]
        storedataframe=as.data.frame(storelist)
        storedataframe[,"rep"]=rep(tabsumm[ind,"rand"],times=dim(storedataframe)[1])
        plotlist[[ind]]=storedataframe
      }
      plottab=cbind(Reduce("rbind",plotlist))
      plottab[,"rep"]=as.factor(plottab[,"rep"])
      temptab=plottab[,c("ind","mean")]
      meantab=aggregate(mean~ind,data=temptab,mean)
      colnames(meantab)=c("ind","mean")
      sdtab=aggregate(mean~ind,data=temptab,sd)
      colnames(sdtab)=c("ind","sd")
      summarytab=merge(meantab,sdtab,by=c("ind"),all=TRUE)

      p<-ggplot(data=summarytab,aes(x=ind,y=mean))+
         geom_ribbon(aes(ymin=mean-2*sd,ymax=mean+2*sd),fill="grey70",colour="grey70")+
         geom_line(size=0.3)+
         # geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
         # facet_wrap(~rep,nrow=3,scales="free_y")+
         labs(title=paste0("correlation vs t"),x="t",y="correlation")+
         theme(axis.text=element_text(size=fondsize))
      ggsave(plot=p,filename=paste0(dir,"correlation_vs_delt_",varele,"_T_",T_plot,"_L_",L_plot,".pdf"))
    }
  }
}
##replicate 1 curve for relaxation
for(T_plot in unique(tabsumm[,"T"])){
  for(L_plot in unique(tabsumm[,"L"])){
    inds=tabsumm[tabsumm[,"L"]==L_plot&tabsumm[,"T"]==T_plot,"H_relaxseq_ind"]#&tabsumm[,"rand"]==rand_plot
    for(varele in c("H","M")){
      plotlist=vector(mode="list")
      for(ind in inds){
        storelist=corrlist[[ind]][[varele]]
        storedataframe=as.data.frame(storelist)
        storedataframe[,"rep"]=rep(tabsumm[ind,"rand"],times=dim(storedataframe)[1])
        plotlist[[ind]]=storedataframe
      }
      plottab=cbind(Reduce("rbind",plotlist))
      plottab[,"rep"]=as.factor(plottab[,"rep"])
      plottab2=plottab[plottab[,"rep"]==1,]
      p<-ggplot(data=plottab2,aes(x=ind,y=mean))+
         # geom_ribbon(aes(ymin=mean-2*sd,ymax=mean+2*sd),fill="grey70",colour="grey70")+
         geom_line(size=0.3)+
         # geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
         # facet_wrap(~rep,nrow=3,scales="free_y")+
         labs(title=paste0("correlation vs t"),x="t",y="correlation")+
         theme(axis.text=element_text(size=fondsize))
      ggsave(plot=p,filename=paste0(dir,"correlation_vs_delt_",varele,"_T_",T_plot,"_L_",L_plot,"_rep1.pdf"))
    }
  }
}

for(T_plot in unique(tabsumm[,"T"])){
  for(L_plot in unique(tabsumm[,"L"])){
    inds=tabsumm[tabsumm[,"L"]==L_plot&tabsumm[,"T"]==T_plot,"H_relaxseq_ind"]#&tabsumm[,"rand"]==rand_plot
    for(varele in c("H","M")){
      plotlist=vector(mode="list")
      for(ind in inds){
        storelist=corrlist[[ind]][[varele]]
        storedataframe=as.data.frame(storelist)
        storedataframe[,"rep"]=rep(tabsumm[ind,"rand"],times=dim(storedataframe)[1])
        plotlist[[ind]]=storedataframe
      }
      plottab=cbind(Reduce("rbind",plotlist))
      plottab[,"rep"]=as.factor(plottab[,"rep"])
      temptab=plottab[,c("ind","mean")]
      meantab=aggregate(mean~ind,data=temptab,mean)
      colnames(meantab)=c("ind","mean")
      sdtab=aggregate(mean~ind,data=temptab,sd)
      colnames(sdtab)=c("ind","sd")
      summarytab=merge(meantab,sdtab,by=c("ind"),all=TRUE)

      p<-ggplot(data=summarytab,aes(x=ind,y=mean))+
         geom_ribbon(aes(ymin=mean-2*sd,ymax=mean+2*sd),fill="grey70",colour="grey70")+
         geom_line(size=0.3)+
         scale_y_continuous(trans='log2')+
         # geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
         # facet_wrap(~rep,nrow=3,scales="free_y")+
         labs(title=paste0("correlation vs t"),x="t",y="correlation")+
         theme(axis.text=element_text(size=fondsize))
      ggsave(plot=p,filename=paste0(dir,"correlation_vs_delt_",varele,"_T_",T_plot,"_L_",L_plot,"log2.pdf"))
    }
  }
}
##replicate 1 curve for relaxation
for(T_plot in unique(tabsumm[,"T"])){
  for(L_plot in unique(tabsumm[,"L"])){
    inds=tabsumm[tabsumm[,"L"]==L_plot&tabsumm[,"T"]==T_plot,"H_relaxseq_ind"]#&tabsumm[,"rand"]==rand_plot
    for(varele in c("H","M")){
      plotlist=vector(mode="list")
      for(ind in inds){
        storelist=corrlist[[ind]][[varele]]
        storedataframe=as.data.frame(storelist)
        storedataframe[,"rep"]=rep(tabsumm[ind,"rand"],times=dim(storedataframe)[1])
        plotlist[[ind]]=storedataframe
      }
      plottab=cbind(Reduce("rbind",plotlist))
      plottab[,"rep"]=as.factor(plottab[,"rep"])
      plottab2=plottab[plottab[,"rep"]==1,]
      p<-ggplot(data=plottab2,aes(x=ind,y=mean))+
         # geom_ribbon(aes(ymin=mean-2*sd,ymax=mean+2*sd),fill="grey70",colour="grey70")+
         geom_line(size=0.3)+
         scale_y_continuous(trans='log2')+
         # geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
         # facet_wrap(~rep,nrow=3,scales="free_y")+
         labs(title=paste0("correlation vs t"),x="t",y="correlation")+
         theme(axis.text=element_text(size=fondsize))
      ggsave(plot=p,filename=paste0(dir,"correlation_vs_delt_",varele,"_T_",T_plot,"_L_",L_plot,"_rep1_log2.pdf"))
    }
  }
}
##exponent estimation (relaxation time)
##relaxation exponent approximation
lnseq1=2:10
lnseq2=2:150
lnseq3=10:150
lnseq4=10:3000
lnseq5=2:4
tabtau=as.data.frame(matrix(NA,nrow=floor(length(locfiles)/3),ncol=10))
colnames(tabtau)=c("L","T","tau_HH","tau_MM","tau_HH_sd","tau_MM_sd","tau_HH_integ","tau_MM_integ","tau_HH_integ_sd","tau_MM_integ_sd")
exporange=c()
rowi=1
for(Tele in unique(tabsumm[,"T"])){
  for(Lele in unique(tabsumm[,"L"])){
    tabtau[rowi,"L"]=Lele
    tabtau[rowi,"T"]=Tele
    inds=tabsumm[tabsumm[,"L"]==Lele&tabsumm[,"T"]==Tele,"H_relaxseq_ind"]
    for(varele in c("H","M")){
      timvec=c()
      corrvec=c()
      lnseq=lnseq1
      if(varele=="M"&&tabtau[rowi,"T"]==2.3&&tabtau[rowi,"L"]==4){
        lnseq=lnseq2
      }
      if(varele=="H"&&tabtau[rowi,"T"]==2.3&&tabtau[rowi,"L"]==32){
        lnseq=lnseq3
      }
      if(varele=="M"&&tabtau[rowi,"T"]==2.3&&tabtau[rowi,"L"]==32){
        lnseq=lnseq4
      }
      if(varele=="H"&&tabtau[rowi,"T"]==5&&tabtau[rowi,"L"]==32){
        lnseq=lnseq5
      }
      for(ind in inds){
        storelist=corrlist[[ind]][[varele]]
        timvec=c(timvec,lnseq)
        corrvec=c(corrvec,storelist[["mean"]][lnseq])
      }
      corrvec[corrvec<0]=min(abs(corrvec))/10##the negative value are placed by 0.1* mini_absolute
      xvec=1/timvec
      yvec=1/log(corrvec)
      modelsummary=summary(lm(yvec~xvec))
      varstr=paste0("tau_",varele,varele)
      tabtau[rowi,varstr]= -modelsummary[["coefficients"]]["xvec","Estimate"]
      tabtau[rowi,paste0(varstr,"_sd")]=modelsummary[["coefficients"]]["xvec","Std. Error"]
    }
    rowi=rowi+1
  }
}

tabtau[,"L"]=as.factor(tabtau[,"L"])
tabtau[,"T"]=as.factor(tabtau[,"T"])
##exponent for H and M plot
p<-ggplot(data=tabtau,aes(x=T,y=tau_HH,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_HH-tau_HH_sd*2,ymax=tau_HH+tau_HH_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=paste0("tau_HH vs T and L"),x="T",y="tau_HH")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_HH_vs_TandL.pdf"))
p<-ggplot(data=tabtau,aes(x=T,y=tau_MM,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_MM-tau_MM_sd*2,ymax=tau_MM+tau_MM_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=paste0("tau_MM vs T and L"),x="T",y="tau_MM")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_MM_vs_TandL.pdf"))

##integrated tau
for(Tele in unique(tabsumm[,"T"])){
  for(Lele in unique(tabsumm[,"L"])){
    for(varele in c("H","M")){
      tauvec=c()
      for(randele in unique(tabsumm[,"rand"])){
        ind=tabsumm[tabsumm[,"L"]==Lele&tabsumm[,"T"]==Tele&tabsumm[,"rand"]==randele,"H_relaxseq_ind"]
        templist=corrlist[[ind]][[varele]]
        x=templist[["ind"]]
        y=templist[["mean"]]
        tauvec=c(tauvec,trapz(x,y))
      }
      tauname=paste0("tau_",varele,varele,"_integ")
      rowind=which(tabtau[,"L"]==Lele&tabtau[,"T"]==Tele)
      tabtau[rowind,tauname]=mean(tauvec)
      tabtau[rowind,paste0(tauname,"_sd")]=sd(tauvec)
    }
  }
}
##exponent for H and M plot
p<-ggplot(data=tabtau,aes(x=T,y=tau_HH_integ,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_HH_integ-tau_HH_integ_sd*2,ymax=tau_HH_integ+tau_HH_integ_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=paste0("tau_HH_integ vs T and L"),x="T",y="tau_HH_integ")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_HH_integ_vs_TandL.pdf"))
p<-ggplot(data=tabtau,aes(x=T,y=tau_MM_integ,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_MM_integ-tau_MM_integ_sd*2,ymax=tau_MM_integ+tau_MM_integ_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=paste0("tau_MM_integ vs T and L"),x="T",y="tau_MM_integ")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_MM_integ_vs_TandL.pdf"))

write.table(tabtau,file="relaxation_time.tab",sep="\t",row.names=FALSE)
