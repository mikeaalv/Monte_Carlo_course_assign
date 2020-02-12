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
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/assignment/assign2/result/"
setwd(dir)
##plot infinite cluster probability with error bar
tab=read.table(paste0(dir,"rangeplot.tab"),header=FALSE,sep=" ")
colnames(tab)=c("L","p","exit","prop")
tabuse=tab[tab[,"L"]>50,]
tabuse_count=aggregate(exit~L+p,data=tabuse,function(x){
  sum(x)/length(x)
})
errtab=read.table(paste0(dir,"rangeplot_uncernty.tab"),header=FALSE,sep=" ")
colnames(errtab)=c("sampi",colnames(tab))
errtab_count=aggregate(exit~L+p+sampi,data=errtab,function(x){
  sum(x)/length(x)
})
sdtab=aggregate(exit~L+p,data=errtab_count,sd)
tabuse_count=merge(tabuse_count,sdtab,by=c("L","p"),all.x=TRUE)
colnames(tabuse_count)=c("L","p","exit","sd")
tabuse_count[,"L"]=as.factor(tabuse_count[,"L"])
tabuse_count[is.na(tabuse_count[,"sd"]),"sd"]=0
tabuse_count=tabuse_count[tabuse_count[,"p"]>0.25&tabuse_count[,"p"]<0.8,]
p<-ggplot(data=tabuse_count,aes(x=p,y=exit,group=L,color=L))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=exit-sd*2,ymax=exit+sd*2))+
   facet_wrap(~L,nrow=3)+
   labs(title=expression(paste("P(inifite cluster) vs p",)),x="p",y="P(inf)")
ggsave(plot=p,filename=paste0(dir,"Prob_vs_p.pdf"))

##example matrix
mat=as.matrix(read.table(paste0(dir,"examp_plot1000.5.tab"),header=FALSE,sep=" ",skip=1))
seqind=seq(dim(mat)[1])
pdf(paste0(dir,"example_100_0.5.pdf"))
p=image(mat,x=seqind,y=seqind,col=c("white","red"))
dev.off()
mat=as.matrix(read.table(paste0(dir,"examp_plot1000.8.tab"),header=FALSE,sep=" ",skip=1))
seqind=seq(dim(mat)[1])
pdf(paste0(dir,"example_100_0.8.pdf"))
image(mat,x=seqind,y=seqind,col=c("white","red"))
dev.off()

##estimate pc
## M=B(p-pc)^beta
## use x axis cross to extrapolate
tabMcal=read.table(paste0(dir,"save_pc_calc.tab"),header=FALSE,sep=" ")
colnames(tabMcal)=c("randi","L","p","exit","prop")
# tabMcal=tabMcal[tabMcal[,"L"]>50,]
tabMcaluse_mean=aggregate(prop~L+p+randi,data=tabMcal,mean)
Lseq=unique(tabMcaluse_mean[,"L"])
randseq=unique(tabMcaluse_mean[,"randi"])
m_estmat=as.data.frame(matrix(NA,nrow=length(Lseq),ncol=3))
colnames(m_estmat)=c("L","pc","sd")
for(Lele_i in seq(length(Lseq))){
  Lele=Lseq[Lele_i]
  m_estmat[Lele_i,"L"]=Lele
  pcvec=c()
  for(randi in randseq){
    loctab=tabMcaluse_mean[tabMcaluse_mean[,"L"]==Lele&tabMcaluse_mean[,"randi"]==randi,]
    pcvec=c(pcvec,max(loctab[loctab[,"prop"]==0,"p"]))
  }
  m_estmat[Lele_i,"pc"]=mean(pcvec)
  m_estmat[Lele_i,"sd"]=sd(pcvec)
}
p<-ggplot(data=m_estmat,aes(x=L,y=pc))+
   geom_line()+
   geom_point()+
   geom_errorbar(aes(ymin=pc-sd*2,ymax=pc+sd*2))+
   labs(title=expression(paste("pc_estimate vs L",)),x="L",y="pc")
ggsave(plot=p,filename=paste0(dir,"pc_vs_L.pdf"))


## estimate tao
## cn=an^(-tao)
## use log transformed and slope to extrapolate
vectao=readLines(paste0(dir,"tao_esti100.tab"))
vectao=vectao[seq(length(vectao))%%2==0]
vectao %>% paste(.,collapse=" ") %>%
           str_split(string=.,pattern=" ") %>%
           extract2(.,1) %>%
           as.numeric(.) -> nvec
cn=table(nvec)/length(vectao)
x=log(as.numeric(names(cn)))
y=log(cn)
pdf(paste0(dir,"logcn_vs_logn100.pdf"))
plot(x,y,xlab="log n",ylab="log cn",main="log(cn) vs log(n) L=100")
dev.off()
summary(lm(y~x))
vectao=readLines(paste0(dir,"tao_esti500.tab"))
vectao=vectao[seq(length(vectao))%%2==0]
vectao %>% paste(.,collapse=" ") %>%
           str_split(string=.,pattern=" ") %>%
           extract2(.,1) %>%
           as.numeric(.) -> nvec
cn=table(nvec)/length(vectao)
x=log(as.numeric(names(cn)))
y=log(cn)
pdf(paste0(dir,"logcn_vs_logn500.pdf"))
plot(x,y,xlab="log n",ylab="log cn",main="log(cn) vs log(n) L=500")
dev.off()
summary(lm(y~x))
