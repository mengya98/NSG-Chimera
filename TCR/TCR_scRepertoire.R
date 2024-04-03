setwd('G:/Chimera/VDJ/')

library(Seurat)
library(scRepertoire)
library(ggsci)
library(ggpubr)
library(patchwork)
library(reshape2)
library(viridis)
library(RColorBrewer)

tcell<-readRDS('G:/Chimera/RNA/All/T/tcell_seurat.rds')
meta<-data.frame(tcell@meta.data)

meta$barcode<-rownames(meta)
meta<-meta[meta$celltype!='Gamma delta T',]
meta<-meta[meta$tissue!='BM',]

#S1<-read.csv('G:/Chimera/VDJ/filtered/CON_BM/vdj_t/filtered_contig_annotations.csv',stringsAsFactors = F)
S2<-read.csv('G:/Chimera/VDJ/filtered/CON_LN/vdj_t/filtered_contig_annotations.csv',stringsAsFactors = F)
S3<-read.csv('G:/Chimera/VDJ/filtered/CON_TM/vdj_t/filtered_contig_annotations.csv',stringsAsFactors = F)
#S4<-read.csv('G:/Chimera/VDJ/filtered/CHI_BM/vdj_t/filtered_contig_annotations.csv',stringsAsFactors = F)
S5<-read.csv('G:/Chimera/VDJ/filtered/CHI_LN/vdj_t/filtered_contig_annotations.csv',stringsAsFactors = F)
S6<-read.csv('G:/Chimera/VDJ/filtered/CHI_TM/vdj_t/filtered_contig_annotations.csv',stringsAsFactors = F)
S2$sample_barcode<-paste0('CON_LN','_',S2$barcode)
S2<-S2[S2$sample_barcode%in%c(meta$barcode),]
S3$sample_barcode<-paste0('CON_TM','_',S3$barcode)
S3<-S3[S3$sample_barcode%in%c(meta$barcode),]
S5$barcode<-paste0(sapply(strsplit(as.character(S5$barcode),'-'),'[[',1),'.1')
S5$sample_barcode<-paste0('CHI_LN','_',S5$barcode)
S5<-S5[S5$sample_barcode%in%c(meta$barcode),]
S6$barcode<-paste0(sapply(strsplit(as.character(S6$barcode),'-'),'[[',1),'.1')
S6$sample_barcode<-paste0('CHI_TM','_',S6$barcode)
S6<-S6[S6$sample_barcode%in%c(meta$barcode),]

write.table(S2,file = "./TCR_filter/CON_LN_filtered_contig_annotations.csv", row.names = F, col.names = T, quote = F, sep = ",")
write.table(S3,file = "./TCR_filter/CON_TM_filtered_contig_annotations.csv", row.names = F, col.names = T, quote = F, sep = ",")
write.table(S5,file = "./TCR_filter/CHI_LN_filtered_contig_annotations.csv", row.names = F, col.names = T, quote = F, sep = ",")
write.table(S6,file = "./TCR_filter/CHI_TM_filtered_contig_annotations.csv", row.names = F, col.names = T, quote = F, sep = ",")

meta_celltype<-meta[,c(10,13)]

#CON CHI-BM/LN/TM
contig_list <- list(S2, S5, S3, S6)
combined <- combineTCR(contig_list,
                       samples = c("CON","CHI","CON","CHI"),
                       ID = c("LN","LN","TM","TM"),
                       cells ="T-AB")
combined_S2<-data.frame(combined$CON_LN)
combined_S2<-merge(combined_S2,meta_celltype,by = 'barcode')
combined$CON_LN<-combined_S2
combined_S3<-data.frame(combined$CON_TM)
combined_S3<-merge(combined_S3,meta_celltype,by = 'barcode')
combined$CON_TM<-combined_S3
combined_S5<-data.frame(combined$CHI_LN)
combined_S5<-merge(combined_S5,meta_celltype,by = 'barcode')
combined$CHI_LN<-combined_S5
combined_S6<-data.frame(combined$CHI_TM)
combined_S6<-merge(combined_S6,meta_celltype,by = 'barcode')
combined$CHI_TM<-combined_S6

#quantContig
quantContig_output<-quantContig(combined, cloneCall="gene+nt", scale = T,group.by = 'sample',exportTable = T) 
quantContig_output$sample<-factor(quantContig_output$sample,levels = c('CON','CHI'))

pdf('quantContig(gene+nt)_origin.pdf',width = 4.5,height =4)
ggplot(quantContig_output,aes(x=sample,y=scaled,fill=sample))+
  stat_summary(fun='mean',geom = 'bar', width = 0.7,position = position_dodge(), colour = "black")+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge(.6))+
  scale_fill_manual(name='Sample',values=c("CON"='#00BFC4',"CHI"='#F8766D')) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 12, color = "black"),axis.title.y=element_text(size=14, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(plot.title = element_text(size=16,hjust = 0.5,color = "black"),legend.title=element_text(size=14,color = "black"),
        legend.text=element_text(size=14,color = "black"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,110),breaks = c(0,50,100))+
  geom_signif(comparisons = list(c("CON","CHI")),
              map_signif_level=T, 
              tip_length=c(0.04,0.04),
              y_position = c(104), 
              size=0.3, 
              textsize = 3, 
              test = "t.test")+
  labs(x = "",y = "Percent of Unique Clonotypes", title = "Combined chain visualization")
dev.off()
write.table(quantContig_output,file = "quantContig_output.xls", row.names = F, col.names = T, quote = F, sep = "\t")

##Length of Clonotypes##
lengthContig_celltype<-lengthContig(combined, cloneCall="aa", chain = "both",group.by = 'celltype',exportTable = T)
lengthContig_celltype$sample<-sapply(strsplit(lengthContig_celltype$values,"_"),"[",1)
lengthContig_celltype$tissue<-sapply(strsplit(lengthContig_celltype$values,"_"),"[",2)
lengthContig_celltype$celltype<-factor(lengthContig_celltype$celltype,levels = c('DP','Naive CD4','Naive CD8','Treg'))
lengthContig_celltype$TRA<-sapply(strsplit(as.character(lengthContig_celltype$CT),'_'),'[[',1)
lengthContig_celltype$TRB<-sapply(strsplit(as.character(lengthContig_celltype$CT),'_'),'[[',2)
lengthContig_part_celltype<-lengthContig_celltype
lengthContig_part_celltype<-lengthContig_part_celltype %>%tidyr::separate(col = "TRA", sep = ";", into = paste0("TRA_", 1:2))
lengthContig_part_celltype<-lengthContig_part_celltype %>%tidyr::separate(col = "TRB", sep = ";", into = paste0("TRB_", 1:2))
lengthContig_part_celltype<-lengthContig_part_celltype[,-c(1:2)]
lengthContig_part_celltype<-melt(lengthContig_part_celltype,id.var=c('celltype','values','sample','tissue'),value.name = 'CT')
lengthContig_part_celltype<-lengthContig_part_celltype[lengthContig_part_celltype$CT!='NA',]
lengthContig_part_celltype<-na.omit(lengthContig_part_celltype)
lengthContig_part_celltype$length<-nchar(lengthContig_part_celltype$CT)
lengthContig_part_celltype$count<-1
lengthContig_part_origin_length<-aggregate(lengthContig_part_celltype$count,by = list(lengthContig_part_celltype$length,lengthContig_part_celltype$sample),FUN = sum)
colnames(lengthContig_part_origin_length)<-c('length','sample','count')
lengthContig_part_origin_length<-lengthContig_part_origin_length[-19,]
lengthContig_part_origin_length_supp<-data.frame(length=c(10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5),
                                                 sample=c('CHI','CHI','CHI','CHI','CHI','CHI','CHI','CHI','CHI','CHI','CON','CON','CON','CON','CON','CON','CON','CON','CON','CON'),
                                                 count=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
lengthContig_part_origin_length<-rbind(lengthContig_part_origin_length,lengthContig_part_origin_length_supp)

pdf('CDR3number_lengthContig_origin(aa).pdf',width = 5.1,height = 4)
ggplot(lengthContig_part_origin_length)+
  geom_line(aes(x=length,y=count,color=sample),size=0.5)+
  scale_color_manual(name='Sample',values=c("CON"='#00BFC4',"CHI"='#F8766D'))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12, color = "black"),axis.title.y=element_text(size=14, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(plot.title = element_text(size=16,hjust = 0.5,color = "black"),legend.title=element_text(size=14,color = "black"),
        legend.text=element_text(size=14,color = "black"))+
  scale_y_continuous(limits = c(0,10000))+
  labs(x = "CDR3 length (aa)",y = "Number of CDR3", title = "Combined chain visualization")
dev.off()

lengthContig_part_celltype_length<-aggregate(lengthContig_part_celltype$count,by = list(lengthContig_part_celltype$length,lengthContig_part_celltype$celltype,lengthContig_part_celltype$sample,lengthContig_part_celltype$tissue),FUN = sum)
colnames(lengthContig_part_celltype_length)<-c('length','celltype','sample','tissue','count')
lengthContig_part_celltype_length<-lengthContig_part_celltype_length[-90,]
lengthContig_part_celltype_length$sample<-factor(lengthContig_part_celltype_length$sample,levels = c('CON','CHI'))
pdf('lengthContig_celltype(aa).pdf',width = 6,height = 4)
ggplot(lengthContig_part_celltype_length,aes(x=sample,y=length,fill=sample))+
  stat_summary(fun='mean',geom = 'bar', width = 0.7,position = position_dodge(), colour = "black")+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge(.6))+
  scale_fill_manual(name='Sample',values=c("CON"='#00BFC4',"CHI"='#F8766D')) +
  theme_classic()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_text(size=10, color = "black"))+
  theme(axis.text.y = element_text(size = 10, color = "black"))+
  theme(plot.title = element_text(size=14,hjust = 0.5,color = "black"),legend.title=element_text(size=10,color = "black"),
        legend.text=element_text(size=10,color = "black"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,30))+
  geom_signif(comparisons = list(c("CON","CHI")),
              map_signif_level=T, 
              tip_length=c(0.04,0.04), 
              y_position = c(25),
              size=0.3, 
              textsize = 3, 
              test = "t.test")+
  facet_grid(tissue~celltype)+
  theme(strip.placement = "outside",
        panel.spacing.y = unit(0.5, "cm"))+
  labs(x = "Celltype",y = "CDR3 length (aa)", title = "")
dev.off()

#Diversity Analysis
pdf('F1E, clonal_diversity(gene).pdf',width = 8,height = 4)
clonalDiversity(combined,cloneCall = "gene",group.by = "sample",x.axis = "ID",n.boots = 100)+ 
  scale_color_manual(values=c('#F8766D','#00BFC4'))+
  theme(axis.text = element_text(color = "black"))
dev.off()

##Interacting with Single-Cell Objects##
seurat <- combineExpression(combined, tcell, 
                            cloneCall="gene", 
                            group.by = "sample", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=2, Large=500))
Idents(seurat)<-seurat$celltype
slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType,
                                              levels = c("Large (2 < X <= 500)",
                                                         "Small (1 < X <= 2)",
                                                         "Single (0 < X <= 1)", NA))
colorblind_vector <- c('#FC800F','#289A6A','#1A72B4')
pdf('seurat_cloneType_3levels.pdf',width = 8.1,height = 6)
DimPlot(seurat, group.by = "cloneType",pt.size = 0.5) +
  scale_color_manual(values = colorblind_vector, na.value="grey") + 
  theme(plot.title = element_blank())
dev.off()

#T subset clonotypes
seurat_CON<-subset(seurat,origin=='CON')
proportion_CON_3levels<-occupiedscRepertoire(seurat_CON, x.axis = "ident",exportTable =T)
seurat_CHI<-subset(seurat,origin=='CHI')
proportion_CHI_3levels<-occupiedscRepertoire(seurat_CHI, x.axis = "ident",exportTable =T)

#put 3levels CON/CHI in a plot
proportion_CON_3levels$origin<-'CON'
proportion_CHI_3levels$origin<-'CHI'
proportion_merge_3levels<-rbind(proportion_CON_3levels,proportion_CHI_3levels)
proportion_merge_3levels$origin_celltype<-paste0(proportion_merge_3levels$origin,'-',proportion_merge_3levels$ident)
proportion_merge_3levels$origin_celltype<-factor(proportion_merge_3levels$origin_celltype,
                                                 levels = rev(c('CON-DP','CHI-DP','CON-Naive CD4','CHI-Naive CD4',
                                                                'CON-Naive CD8','CHI-Naive CD8','CON-Treg','CHI-Treg')))
pdf("Tsubset_proportion_merge_3levels.pdf",width = 7,height = 4)
ggplot(proportion_merge_3levels,aes(x=origin_celltype,fill=cloneType,y=value))+
  geom_bar(stat = 'identity',
           position = 'fill',width = 0.7)+ 
  labs(x = "",y = "Proportion", title = "",fill='T subsets')+
  theme_classic()+
  scale_fill_manual(values = c('#FC800F','#289A6A','#1A72B4'))+
  #guides(fill=guide_legend(title=),size)+
  theme(axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        legend.text = element_text(color="black", size = 12))+
  theme(axis.text = element_text(colour = 'black'))+
  coord_flip()
dev.off()
