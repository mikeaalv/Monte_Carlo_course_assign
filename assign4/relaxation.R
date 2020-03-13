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
mp_cal<-function(x,tseq,mean,var){
  dim(x)=c(tseq,length(x)%/%tseq)
  x=t(x)
  refx=x[,1]
  samp_mat=x*as.vector(refx)
  samp_mat=(samp_mat-mean^2)/var
  returnlist=list(mean=colMeans(samp_mat),sd=apply(samp_mat,2,sd))
  return(returnlist)
}
##general sampling
heatup=100
endind=50099
stepsize=1
##relaxation sampling
t0step=10
tseq=10
ngroup_t0=floor((endind-heatup)/t0step)
##relaxation exponent approximation
ln_range=c(10,20)
lnseq=ln_range[1]:ln_range[2]
tottimeseq=0:(tseq-1)
locfiles=list.files(pattern="gener\\.tab$")
tabsumm=as.data.frame(matrix(NA,nrow=length(locfiles),ncol=9))
colnames(tabsumm)=c("L","T","rand","H","H_var","M","M_var","H_relaxseq_ind","M_relaxseq_ind")
sampleseq_balk=seq(from=heatup,by=stepsize,to=endind)
sampleseq_relax=sampleseq_balk[rep(seq(from=1,by=1,to=tseq),times=ngroup_t0)+rep(t0step*seq(from=0,by=1,to=ngroup_t0-1),each=tseq)]
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
  tabsumm[rowi,"H"]=mean(datatab[sampleseq_balk,"H"])
  tabsumm[rowi,"M"]=mean(datatab[sampleseq_balk,"M"])
  tabsumm[rowi,"H_var"]=var(datatab[sampleseq_balk,"H"])
  tabsumm[rowi,"M_var"]=var(datatab[sampleseq_balk,"M"])
  corr_vec_list=vector(mode="list")
  # H_corr_infor=mp_cal(datatab[sampleseq_relax,"H"],tseq)
  # M_corr_infor=mp_cal(datatab[sampleseq_relax,"M"],tseq)
  # corr_vec_list[["H"]]=list(ind=tottimeseq,
  #                           mean=(H_corr_infor[["mean"]]-tabsumm[rowi,"H"]^2)/tabsumm[rowi,"H_var"],
  #                           sd=H_corr_infor[["sd"]]/tabsumm[rowi,"H_var"])
  # corr_vec_list[["M"]]=list(ind=tottimeseq,
  #                           mean=(M_corr_infor[["mean"]]-tabsumm[rowi,"M"]^2)/tabsumm[rowi,"M_var"],
  #                           sd=M_corr_infor[["sd"]]/tabsumm[rowi,"M_var"])
  H_corr_infor=mp_cal(datatab[sampleseq_relax,"H"],tseq,tabsumm[rowi,"H"],tabsumm[rowi,"H_var"])
  M_corr_infor=mp_cal(datatab[sampleseq_relax,"M"],tseq,tabsumm[rowi,"M"],tabsumm[rowi,"M_var"])
  corr_vec_list[["H"]]=list(ind=tottimeseq,
                            mean=H_corr_infor[["mean"]],
                            sd=H_corr_infor[["sd"]])
  corr_vec_list[["M"]]=list(ind=tottimeseq,
                            mean=M_corr_infor[["mean"]],
                            sd=M_corr_infor[["sd"]])
  corrlist[[rowi]]=corr_vec_list
  tabsumm[rowi,"H_relaxseq_ind"]=rowi
  tabsumm[rowi,"M_relaxseq_ind"]=rowi
  rowi=rowi+1
}
##test plot
# L_plot=32#4 32
# T_plot=2.3#2.3 5.0
# rand_plot=1
for(T_plot in unique(tabsumm[,"T"])){
  for(L_plot in unique(tabsumm[,"L"])){
    inds=tabsumm[tabsumm[,"L"]==L_plot&tabsumm[,"T"]==T_plot,"H_relaxseq_ind"]#&tabsumm[,"rand"]==rand_plot
    for(var in c("H","M")){
      plotlist=vector(mode="list")
      for(ind in inds){
        storelist=corrlist[[ind]][[var]]
        storedataframe=as.data.frame(storelist)
        storedataframe[,"rep"]=rep(tabsumm[ind,"rand"],times=dim(storedataframe)[1])
        plotlist[[ind]]=storedataframe
      }
      plottab=cbind(Reduce("rbind",plotlist))
      plottab[,"rep"]=as.factor(plottab[,"rep"])
      p<-ggplot(data=plottab,aes(x=ind,y=mean,color=rep,group=rep))+
         geom_line()+
         geom_point()+
         # geom_errorbar(aes(ymin=mean-sd*2,ymax=mean+sd*2))+
         # facet_wrap(~rep,nrow=3,scales="free_y")+
         labs(title=expression(paste("correlation vs t",)),x="t",y="correlation")+
         theme(axis.text=element_text(size=fondsize))
      ggsave(plot=p,filename=paste0(dir,"correlation_vs_delt_",var,"_T_",T_plot,"_L_",L_plot,".pdf"))
    }
  }
}

##exponent estimation (relaxation time)
tabtau=as.data.frame(matrix(NA,nrow=floor(length(locfiles)/3),ncol=10))
colnames(tabtau)=c("L","T","tau_HH","tau_MM","tau_HH_sd","tau_MM_sd","tau_HH_integ","tau_MM_integ","tau_HH_integ_sd","tau_MM_integ_sd")
rowi=1
for(Tele in unique(tabsumm[,"T"])){
  for(Lele in unique(tabsumm[,"L"])){
    tabtau[rowi,"L"]=Lele
    tabtau[rowi,"T"]=Tele
    inds=tabsumm[tabsumm[,"L"]==Lele&tabsumm[,"T"]==Tele,"H_relaxseq_ind"]
    for(var in c("H","M")){
      timvec=c()
      corrvec=c()
      for(ind in inds){
        storelist=corrlist[[ind]][[var]]
        timvec=c(timvec,lnseq)
        corrvec=c(corrvec,storelist[["mean"]])
      }
      xvec=1/timvec
      yvec=1/log(corrvec)
      modelsummary= -summary(lm(yvec~xvec))
      varstr=paste("tau_",var,var)
      tabtau[rowi,varstr]=modelsummary[["coefficients"]]["xvec","Estimate"]
      tabtau[rowi,paste(varstr,"_sd")]=modelsummary[["coefficients"]]["xvec","Std. Error"]
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
   labs(title=expression(paste("tau_HH vs T and L",)),x="T",y="tau_HH")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_HH_vs_TandL.pdf"))
p<-ggplot(data=tabtau,aes(x=T,y=tau_MM,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_MM-tau_MM_sd*2,ymax=tau_MM+tau_MM_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("tau_MM vs T and L",)),x="T",y="tau_MM")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_MM_vs_TandL.pdf"))
##integrated tau
for(Tele in unique(tabsumm[,"T"])){
  for(Lele in unique(tabsumm[,"L"])){
    for(var in c("H","M")){
      tauvec=c()
      for(randele in tabsumm[,"rand"]){
        inds=tabsumm[tabsumm[,"L"]==Lele&tabsumm[,"T"]==Tele&tabsumm[,"rand"]==randele,"H_relaxseq_ind"]
        templist=corrlist[[ind]][[var]]
        x=tottimeseq
        y=templist[["mean"]]
        tauvec=c(tauvec,trapz(x,y))
      }
      tauname=paste("tau_",var,var,"_integ")
      rowind=tabtau[,"L"]==Lele&tabtau[,"T"]==Tele
      tabtau[rowind,tauname]=mean(tauvec)
      tabtau[rowind,paste(tauname,"_sd")]=sd(tauvec)
    }
  }
}
##exponent for H and M plot
p<-ggplot(data=tabtau,aes(x=T,y=tau_HH_integ,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_HH_integ-tau_HH_integ_sd*2,ymax=tau_HH_integ+tau_HH_integ_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("tau_HH_integ vs T and L",)),x="T",y="tau_HH_integ")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_HH_integ_vs_TandL.pdf"))
p<-ggplot(data=tabtau,aes(x=T,y=tau_MM_integ,color=L,group=L))+
   geom_point()+
   geom_errorbar(aes(ymin=tau_MM_integ-tau_MM_integ_sd*2,ymax=tau_MM_integ+tau_MM_integ_sd*2))+
   facet_wrap(~L,nrow=3,scales="free_y")+
   labs(title=expression(paste("tau_MM_integ vs T and L",)),x="T",y="tau_MM_integ")+
   theme(axis.text=element_text(size=fondsize))
ggsave(plot=p,filename=paste0(dir,"tau_MM_integ_vs_TandL.pdf"))
