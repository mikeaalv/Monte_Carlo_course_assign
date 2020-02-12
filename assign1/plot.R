### this is the script for plotting
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
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/assignment/assign1/result/"
setwd(dir)
n_indep_chain=10
datafilelist=list.files(".",pattern="*.tab")
# column: sample size; estimated pi
list_pi_array=vector(mode="list")
for(datafile in datafilelist){
  pi_array_temp=matrix(NA,nrow=n_indep_chain,ncol=2)
  pivec=scan(paste0(dir,datafile))
  datafile %>% str_extract(string=.,pattern="\\_\\d+\\.") %>%
               str_replace_all(string=.,pattern="\\_|\\.",replacement="") %>%
               as.numeric(.) -> sample_size
  pi_array_temp[,1]=rep(sample_size,times=n_indep_chain)
  pi_array_temp[,2]=pivec
  list_pi_array[[datafile]]=pi_array_temp
}
pi_array=as.data.frame(Reduce('rbind',list_pi_array))
colnames(pi_array)=c("sample_size","esti_pi")
# plot pi ~ sample size
p<-ggplot(data=pi_array,aes(x=sample_size,y=esti_pi))+
   geom_point()+
   stat_summary(aes(y=esti_pi,group=1),fun.y=mean,colour="red",geom="line",group=1)+
   scale_x_log10()+
   labs(title=expression(paste("estimated ",pi," vs N",)),x="N",y=expression(paste("estimated ",pi)))
ggsave(plot=p,filename=paste0(dir,"est_pi_vs_samplesize.pdf"))
p<-ggplot(data=pi_array[pi_array[,1]>100000,],aes(x=sample_size,y=esti_pi))+
   geom_point()+
   stat_summary(aes(y=esti_pi,group=1),fun.y=mean,colour="red",geom="line",group=1)+
   stat_summary(aes(label=round(..y..,6)),fun.y=mean,geom="text",size=3,colour="blue",vjust=-2.5)+
   scale_x_log10()+
   labs(title=expression(paste("estimated ",pi," vs N",)),x="N",y=expression(paste("estimated ",pi)))
ggsave(plot=p,filename=paste0(dir,"est_pi_vs_samplesize2.pdf"))
# plot sd(pi) ~ sample size
pi_array_temp=pi_array
pi_array_temp[,"sample_size"]=as.factor(pi_array_temp[,"sample_size"])
pi_array_sd=aggregate(.~sample_size,data=pi_array_temp,sd)
pi_array_sd[,"sample_size"]=as.numeric(as.character(pi_array_sd[,"sample_size"]))
p<-ggplot(data=pi_array_sd,aes(x=sample_size,y=esti_pi))+
   geom_line()+
   scale_x_log10()+
   scale_y_log10()+
   labs(title="estimation precision",x="N",y="estimation standard error")
ggsave(plot=p,filename=paste0(dir,"esti_precision.pdf"))
# plot accuracy(pi) ~ sample size
# by mse: mean squared error
pi_array_accuracy=aggregate(.~sample_size,data=pi_array_temp,function(x){
  sqrt(sum((x-pi)^2)/length(x))
})
pi_array_accuracy[,"sample_size"]=as.numeric(as.character(pi_array_accuracy[,"sample_size"]))
p<-ggplot(data=pi_array_accuracy,aes(x=sample_size,y=esti_pi))+
   geom_line()+
   scale_x_log10()+
   scale_y_log10()+
   labs(title="estimation accuracy",x="N",y="RMSE")
ggsave(plot=p,filename=paste0(dir,"esti_accuracy.pdf"))
