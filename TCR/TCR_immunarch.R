setwd('G:/Chimera/VDJ/')

library(immunarch)
library(ggpubr)
library(ggforce)

combined_single<-repLoad('G:/Chimera/VDJ/TCR_filter',.mode = 'single')
combined_single$meta$sample<-rep(c('CHI','129'),c(4,4))

TRAV<-geneUsage(combined_single$data[c(1,3,5,7)], "musmus.trav",.norm = T)
TRAV[is.na(TRAV)] <- 0
TRAJ<-geneUsage(combined_single$data[c(1,3,5,7)], "musmus.traj",.norm = T)

TRBV<-geneUsage(combined_single$data[c(2,4,6,8)], "musmus.trbv",.norm = T)
TRBJ<-geneUsage(combined_single$data[c(2,4,6,8)], "musmus.trbj",.norm = T)

#VDJ gene usage merge in one plot
TRA_merge<-rbind(TRAV,TRAJ)
TRA_merge[,2:5]<-scale(TRA_merge[,2:5])
rownames(TRA_merge)<-TRA_merge$Names
TRA_merge<-data.frame(t(TRA_merge),check.rows = F,check.names = F)
TRA_merge<-TRA_merge[-1,]
TRA_merge$Sample<-rep(c('CHI','CON'),c(2,2))
TRA_merge<-melt(TRA_merge,measure.vars = colnames(TRA_merge)[1:113],variable.name = "gene",value.name = "value")
TRA_merge$value<-as.numeric(TRA_merge$value)
mean<-aggregate(TRA_merge$value,by=list(TRA_merge$Sample,TRA_merge$gene),FUN=mean) #计算均值
sd<-aggregate(TRA_merge$value,by=list(TRA_merge$Sample,TRA_merge$gene),FUN=sd) #计算标准差
N<-aggregate(TRA_merge$value,by=list(TRA_merge$Sample,TRA_merge$gene),FUN=length) #计算个数
TRA_data<-data.frame(mean,sd=sd$x,N=N$x) #合并数据框
colnames(TRA_data)=c("Sample","gene","value","sd","N")
TRA_data$se <- TRA_data$sd / sqrt(TRA_data$N) #计算标准误
pdf("./TRA_usage_merge.pdf",width = 18,height = 6)
ggplot(TRA_data,aes(x=gene,y=value,color=Sample,group=Sample))+
  geom_point(shape=16,size=3)+
  theme_bw() + theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=0,size=1)+
  theme(axis.text.x = element_text(size = 8, color = "black",angle = 90,vjust = 0.5,hjust = 1))+
  theme(axis.text.y = element_text(size = 8, color = "black"),axis.title.y=element_text(size=10))+
  theme(plot.title = element_text(size=10,hjust = 0.5),legend.title=element_text(size=10),
        legend.key = element_rect(colour = "transparent", fill = "white"),legend.text=element_text(size=10))+
  scale_color_manual(name='Sample',values=c("CON"='#00BFC4',"CHI"='#F8766D'))+
  labs(x = "",y = "Relative usage (z-score)", title = "TRAV genes                            TRAJ genes")#+
dev.off()

TRB_merge<-rbind(TRBV,TRBJ)
TRB_merge[,2:5]<-scale(TRB_merge[,2:5])
rownames(TRB_merge)<-TRB_merge$Names
TRB_merge<-data.frame(t(TRB_merge),check.rows = F,check.names = F)
TRB_merge<-TRB_merge[-1,]
TRB_merge$Sample<-rep(c('CHI','CON'),c(2,2))
TRB_merge<-melt(TRB_merge,measure.vars = colnames(TRB_merge)[1:34],variable.name = "gene",value.name = "value")
TRB_merge$value<-as.numeric(TRB_merge$value)
mean<-aggregate(TRB_merge$value,by=list(TRB_merge$Sample,TRB_merge$gene),FUN=mean) #计算均值
sd<-aggregate(TRB_merge$value,by=list(TRB_merge$Sample,TRB_merge$gene),FUN=sd) #计算标准差
N<-aggregate(TRB_merge$value,by=list(TRB_merge$Sample,TRB_merge$gene),FUN=length) #计算个数
TRB_data<-data.frame(mean,sd=sd$x,N=N$x) #合并数据框
colnames(TRB_data)=c("Sample","gene","value","sd","N")
TRB_data$se <- TRB_data$sd / sqrt(TRB_data$N) #计算标准误
pdf("./TRB_usage_merge.pdf",width = 6,height = 6)
ggplot(TRB_data,aes(x=gene,y=value,color=Sample,group=Sample))+
  geom_point(shape=16,size=3)+
  theme_bw() + theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=0,size=1)+
  theme(axis.text.x = element_text(size = 8, color = "black",angle = 90,vjust = 0.5,hjust = 1))+
  theme(axis.text.y = element_text(size = 8, color = "black"),axis.title.y=element_text(size=10))+
  theme(plot.title = element_text(size=10,hjust = 0.5),legend.title=element_text(size=10),
        legend.key = element_rect(colour = "transparent", fill = "white"),legend.text=element_text(size=10))+
  scale_color_manual(name='Sample',values=c("CON"='#00BFC4',"CHI"='#F8766D'))+
  labs(x = "",y = "Relative usage (z-score)", title = "TRBV genes                 TRBJ genes")
dev.off()
