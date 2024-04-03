setwd("G:/Chimera/RNA/All/")

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)
library(harmony)

con_bm<-Read10X(data.dir = "G:/Chimera/RNA/CON/BM/raw_feature_bc_matrix/")
con_bm<-CreateSeuratObject(counts = con_bm,Project="con_bm",min.cells = 3,min.features = 200)
con_ln<-Read10X(data.dir = "G:/Chimera/RNA/CON/LN/raw_feature_bc_matrix/")
con_ln<-CreateSeuratObject(counts = con_ln,Project="con_ln",min.cells = 3,min.features = 200)
con_tm<-Read10X(data.dir = "G:/Chimera/RNA/CON/TM/raw_feature_bc_matrix/")
con_tm<-CreateSeuratObject(counts = con_tm,Project="con_tm",min.cells = 3,min.features = 200)
chi_bm<-Read10X(data.dir = "G:/Chimera/RNA/CHI/BM/raw_feature_bc_matrix/")
chi_bm<-CreateSeuratObject(counts = chi_bm,Project="chi_bm",min.cells = 3,min.features = 200)
chi_ln<-Read10X(data.dir = "G:/Chimera/RNA/CHI/LN/raw_feature_bc_matrix/")
chi_ln<-CreateSeuratObject(counts = chi_ln,Project="chi_ln",min.cells = 3,min.features = 200)
chi_tm<-Read10X(data.dir = "G:/Chimera/RNA/CHI/TM/raw_feature_bc_matrix/")
chi_tm<-CreateSeuratObject(counts = chi_tm,Project="chi_tm",min.cells = 3,min.features = 200)

##merge 6 samples into one##
all <- merge(con_bm, y = c(con_ln,con_tm,chi_bm,chi_ln,chi_tm), add.cell.ids = c("CON_BM", "CON_LN", "CON_TM",'CHI_BM','CHI_LN','CHI_TM'), project = "all")
all[["percent.mito"]] <- PercentageFeatureSet(object = all, pattern = "^mt-")
all[["origin"]] <- sapply(strsplit(colnames(all),'_'),'[',1)
all[["tissue"]] <- sapply(strsplit(colnames(all),'_'),'[',2)
all$sample<-paste0(all$origin,'_',all$tissue)
all <- subset(x = all, subset = nFeature_RNA > 500 & nCount_RNA < 30000 & percent.mito < 10 )
saveRDS(all,'all_merge_seurat.rds')
remove(con_bm,con_ln,con_tm,chi_bm,chi_ln,chi_tm)

pdf("nFeature_RNA.pdf",width = 4.5,height = 4)
VlnPlot(object = all, features = c("nFeature_RNA"),pt.size = 0,group.by = 'sample')
dev.off()
pdf("nCount_RNA.pdf",width = 4.6,height = 4)
VlnPlot(object = all, features = c("nCount_RNA"), pt.size = 0,group.by = 'sample')
dev.off()
pdf("nFeature_RNA_nCount_RNA.pdf",width = 7,height = 6)
FeatureScatter(object = all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample') 
dev.off()

all <- NormalizeData(object = all, normalization.method = "LogNormalize", scale.factor = 1e5)
all <- FindVariableFeatures(object = all,selection.method = 'vst', nfeatures = 3000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = VariableFeatures(object = all))

#harmony
all <- all %>%
  RunHarmony("origin", plot_convergence = TRUE)

#Visualization
all <- FindNeighbors(object = all, reduction="harmony",dims = 1:30)
all <- FindClusters(object = all, resolution = 0.8)
all <- RunTSNE(object = all,reduction="harmony",dims = 1:30)

all$celltype<-''
all$celltype[all$seurat_clusters%in%c(0,1,2,4,6,8,9,12,14,16,23)]<-'T'
all$celltype[all$seurat_clusters%in%c(3,5,13,21,26,30)]<-'B'
all$celltype[all$seurat_clusters%in%c(7,10,11,22,28)]<-'Neutrophil'
all$celltype[all$seurat_clusters%in%c(15,17,27)]<-'Monocyte'
all$celltype[all$seurat_clusters%in%c(25)]<-'Macrophage'
all$celltype[all$seurat_clusters%in%c(18)]<-'Erythrocyte'
all$celltype[all$seurat_clusters%in%c(19)]<-'pDC'
all$celltype[all$seurat_clusters%in%c(20)]<-'cDC'
all$celltype[all$seurat_clusters%in%c(24)]<-'NK'
all$celltype[all$seurat_clusters%in%c(29)]<-'Basophil'

cols1<-c('darkgoldenrod2','dodgerblue','cornflowerblue','lightblue','lightsalmon4',
         'lemonchiffon4','darkorange','darkorange3','deepskyblue','lightsalmon')
pdf("all.Celltype.pdf",width = 7.5,height = 6)
DimPlot(object = all,group.by = 'celltype',reduction = 'tsne',pt.size = 0.5,cols = cols1)
dev.off()

pdf("all.sample.pdf",width = 7.5,height = 6)
DimPlot(object = all,group.by = 'sample',reduction = 'tsne',pt.size = 0.5)
dev.off()

all.markers <- FindAllMarkers(all, only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
all.markers_order<-arrange(all.markers,all.markers$cluster,-all.markers$avg_log2FC)
all.markers_order<-arrange(all.markers_order,all.markers_order$cluster,-all.markers_order$avg_log2FC)

pdf('Dotplot_celltype.pdf',width=11,height=4)
all$celltype<-factor(all$celltype,levels = rev(c('B','Basophil','cDC','Erythrocyte','Macrophage','Monocyte','Neutrophil','NK','pDC','T')))
DotPlot(all,features = c("Cd79a","Cd79b","Cd19",'Ms4a1',
                         "Mcpt8","Fcer1a","Cd200r3","Prss34",
                         "Fscn1","Cst3",'H2-M2','Apol7c',
                         "Tfrc","Gypa","Rhag","Car1",
                         "C1qb","C1qa",'Fcer1g',"Lyz2",
                         "Ms4a6c","F13a1","Ctsc",
                         "S100a8","S100a9","Lcn2","Ly6g",
                         'Nkg7',"Ncr1",'Klrb1a','Gzma',
                         "Irf8","Bst2","Siglech","Cox6a2",
                         "Cd3e","Cd3g","Cd3d","Thy1"),
        dot.scale = 6,group.by = 'celltype')+
  theme_bw()+
  theme(strip.background = element_blank())+
  scale_colour_gradientn(colours = c("#9E0142","#D53E4F","#F46D43",'orange',"#FDAE61","#FEE08B",'#E6F598',"#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))+
  theme(axis.text = element_text(color = "black"))+
  theme(axis.title = element_text(color = "black"))+
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent"))+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
dev.off()

meta<-data.frame(all@meta.data)
meta$sample<-factor(meta$sample,levels = rev(c('CON_BM','CHI_BM','CON_LN','CHI_LN','CON_TM','CHI_TM')))
meta$celltype<-factor(meta$celltype,levels = rev(c('B','Basophil','cDC','Erythrocyte','Macrophage','Monocyte','Neutrophil','NK','pDC','T')))
pdf("Celltype.sample.proportion.pdf",width = 8, height = 4)
ggplot(meta,aes(x=sample,fill=celltype))+
  geom_bar(stat = 'count',position = 'fill')+ 
  labs(x = "",y = "Proportion", title = "",fill='Celltype')+
  theme_classic()+
  scale_fill_manual(values = rev(cols1))+
  theme(axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        legend.text = element_text(color="black", size = 12))+
  theme(axis.text = element_text(colour = 'black'))+
  coord_flip()
dev.off()

#######################t-SNE based on origin of CON/CHI_tomato+/CHI_tomato-#######################
CHI_BM_tomato<-read.table('G:/Chimera/RNA/SNP/CHI_BM/bm.ratio.xls',header = T,sep = '\t')
CHI_BM_tomato<-CHI_BM_tomato[,c(1,5,9)]
colnames(CHI_BM_tomato)[1]<-'identity'
CHI_BM_tomato$identity<-paste0('CHI_BM_',CHI_BM_tomato$identity)
CHI_LN_tomato<-read.table('G:/Chimera/RNA/SNP/CHI_LN/ln.ratio.xls',header = T,sep = '\t')
CHI_LN_tomato<-CHI_LN_tomato[,c(1,7,9)]
colnames(CHI_LN_tomato)[1]<-'identity'
CHI_LN_tomato$identity<-paste0('CHI_LN_',CHI_LN_tomato$identity)
CHI_TM_tomato<-read.table('G:/Chimera/RNA/SNP/CHI_TM/tm.ratio.xls',header = T,sep = '\t')
CHI_TM_tomato<-CHI_TM_tomato[,c(1,7,9)]
colnames(CHI_TM_tomato)[1]<-'identity'
CHI_TM_tomato$identity<-paste0('CHI_TM_',CHI_TM_tomato$identity)
CHI_tomato<-rbind(CHI_BM_tomato,CHI_LN_tomato,CHI_TM_tomato)

meta<-merge(meta,CHI_tomato,by='identity',all.x = T)
meta[is.na(meta)]<-'CON'
rownames(meta)<-meta$identity
meta<-meta[,c(2:13)]
meta_original<-data.frame(all@meta.data)
meta<-meta[rownames(meta_original),]
write.table(data.frame(identity=rownames(meta),meta),file = "./meta_SNP_RFP.xls", row.names = F, col.names = T, quote = F, sep = "\t")
all@meta.data<-meta
cols1<-c('darkgoldenrod2','dodgerblue','cornflowerblue','lightblue','lightsalmon4',
         'lemonchiffon4','darkorange','darkorange3','deepskyblue','lightsalmon')

#group by SNP-RFP
SNP_CON<-subset(all,snp=='CON')
pdf("SNP_CON.tsne.pdf",width = 7.5,height = 6)
DimPlot(object = SNP_CON,reduction = 'tsne',pt.size = 0.5,group.by = 'celltype',cols = cols1)+
  labs(title = 'CON')
dev.off()

SNP_CHI_tomato_plus<-subset(all,snp=='129')
pdf("SNP_CHI_tomato+.tsne.pdf",width = 7.5,height = 6)
DimPlot(object = SNP_CHI_tomato_plus,reduction = 'tsne',pt.size = 0.5,group.by = 'celltype',cols = cols1)+
  labs(title = 'CHI_tomato+')
dev.off()

SNP_CHI_tomato_minus<-subset(all,snp=='nod')
pdf("SNP_CHI_tomato-.tsne.pdf",width = 7.5,height = 6)
DimPlot(object = SNP_CHI_tomato_minus,reduction = 'tsne',pt.size = 0.5,group.by = 'celltype',cols = cols1)+
  labs(title = 'CHI_tomato-')
dev.off()
remove(SNP_CON,SNP_CHI_tomato_plus,SNP_CHI_tomato_minus)

###########################PCA based on markers+CON/CHI_tomato+/CHI_tomato-##########################
meta_SNP<-read.table('meta_SNP_RFP.xls',row.names = 1,header = T)
rna<-data.frame(all@assays$RNA@data[c("Cd19","Cd79a","Cd79b",'Ms4a1','Ly6d',"Ighd",
                                       "Mcpt8","Fcer1a","Cd200r3","Prss34",
                                       "S100a4","Fscn1","Batf3","Cst3",'H2-M2','Apol7c',
                                       "Tfrc","Gypa","Rhag","Car1","Hba-a1","Hbb-bs",
                                       "C1qb","C1qa","Lyz2",'Fcer1g','Csf1r','Fcgr3',
                                       "S100a4","Klf4","F13a1","Csf1r","Ms4a6c","Ctsc",
                                       "S100a8","S100a9","Ly6g","Lyz2","Lcn2","Cd177",
                                       'Klrb1a',"Klrb1c","Ncr1",'Gzma','Nkg7','Ccl5',
                                       "Irf8","Bst2","Siglech","Cox6a2","Cd209d","Pltp",
                                       "Cd3e","Cd3g","Cd3d","Thy1","Cd4","Cd8a"),],check.rows = F,check.names = F)
rna<-data.frame(t(rna),check.rows = F,check.names = F)
meta_SNP$snp_celltype<-paste0(meta_SNP$snp,'_',meta_SNP$celltype)
rna$snp_celltype<-meta_SNP$snp_celltype
attach(rna)
rna<-aggregate(rna[,1:58],by=list(snp_celltype),FUN = mean)
colnames(rna)[1]<-'snp_celltype'
rna<-rna[c(11,1,21,12,2,22,13,3,23,14,4,24,15,5,25,16,6,26,17,7,27,18,8,28,19,9,29,20,10),]
rownames(rna)<-rna$snp_celltype
rna<-rna[,-1]

pca<-prcomp((rna))
plot<-as.data.frame(pca$x[,1:2])
plot$origin<-rep(c("CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-",
                   "CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+","CHI_tomato-","CON","CHI_tomato+"),c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
plot$sample<-rownames(plot)
plot$celltype<-rep(c('B','Basophil','cDC','Erythrocyte','Macrophage','Monocyte','Neutrophil','NK','pDC','T'),
                    c(3,3,3,3,3,3,3,3,3,2))

plot$origin<-factor(plot$origin,levels = c("CON","CHI_tomato+","CHI_tomato-"))
#计算方差,round(x,0)取整,round(x,1)保留一位小数
percentVar<-pca$sdev^2/sum(pca$sdev^2)
pdf("PCA.snp_celltype_markers.pdf",width = 5.2,height = 4)
ggplot(plot,aes(x=PC1,y=PC2,color=celltype,shape=origin))+geom_point(size=4)+
  theme_classic()+
  theme(axis.text= element_text(size=10,color = "black"),axis.title=element_text(size=10))+
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10))+
  theme(plot.title = element_text(size = 12))+
  labs(x = paste0("PC1(", round(percentVar[1]*100,1),"%)"), y = paste0("PC2(", round(percentVar[2]*100,1),"%)"), title = "")+
  scale_color_manual(name='Sample',values=cols1)
dev.off()

##################################Clustering of T subsets############################################
tcell<-subset(all,celltype=='T')
meta_T_old<-data.frame(tcell@meta.data)
meta_T_old<-meta_T_old[,1:7]
remove(all)

tcell.data<-data.frame(tcell@assays$RNA@counts,check.names = F)
tcell<-CreateSeuratObject(counts = tcell.data,Project="tcell",min.cells = 3,min.features = 200)
tcell@meta.data<-meta_T_old
remove(tcell.data)

tcell <- NormalizeData(object = tcell, normalization.method = "LogNormalize", scale.factor = 1e5)
tcell <- FindVariableFeatures(object = tcell,selection.method = 'vst', nfeatures = 3000)
tcell <- ScaleData(object = tcell)
tcell <- RunPCA(object = tcell, features = VariableFeatures(object = tcell))

#harmony
tcell <- tcell %>%
  RunHarmony("origin", plot_convergence = TRUE)

tcell <- FindNeighbors(object = tcell, reduction="harmony",dims = 1:30)
tcell <- FindClusters(object = tcell, resolution = 1.5)
tcell <- RunTSNE(object = tcell,reduction="harmony",dims = 1:30)

pdf("tcell.sample.pdf",width =7,height = 6)
DimPlot(tcell, reduction = "tsne", group.by = "sample", pt.size = 0.5)
dev.off()

tcell.markers <- FindAllMarkers(tcell, only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
tcell.markers_order<-arrange(tcell.markers,tcell.markers$cluster,-tcell.markers$avg_log2FC)
write.table(tcell.markers_order,'./02.Cluster/tcell.marker.xls',row.names = F, col.names = T, quote = F, sep = "\t")

##正式marker##
DP_markers<-c("Cd3g","Cd8a","Cd4")
naive_markers<-c('Sell','Il7r','Ccr7')
Treg_markers<-c("Il2ra","Foxp3","Ikzf2",'Ctla4')
Gamma_delta_T_markers<-c('Tcrg-C1','Tcrg-C2','Tcrg-C4','Trgv2')

tcell$celltype<-''
tcell$celltype[tcell$seurat_clusters%in%c(0,1,3,4,6,7,9,11,12,13,15,17,18,20,27)]='DP'
tcell$celltype[tcell$seurat_clusters%in%c(2,5,19,22,24,28,30)]='Naive CD4'
tcell$celltype[tcell$seurat_clusters%in%c(8,14)]='Naive CD8'
tcell$celltype[tcell$seurat_clusters%in%c(10,16,21,26)]='Treg'
tcell$celltype[tcell$seurat_clusters%in%c(23,25,29)]='Gamma delta T'

pdf("tcell.tsne.pdf",width = 7.9,height = 6)
DimPlot(object = tcell,reduction = 'tsne',group.by = 'celltype',pt.size = 0.5,
        cols = c('#9FDAA4','#468141','#D3ECC2','#A98467','#9ABF6F'))
dev.off()

pdf('tcell_subset_dotplot.pdf',width=14,height=8)
DotPlot(tcell,features = c(DP_markers,naive_markers,Treg_markers,Gamma_delta_T_markers),dot.scale = 6)+
  theme_bw()+
  theme(strip.background = element_blank())+
  scale_colour_gradientn(colours = c("#9E0142","#D53E4F","#F46D43",'orange',"#FDAE61","#FEE08B",'#E6F598',"#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),
                         values = c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))+
  theme(axis.text = element_text(color = "black"))+
  theme(axis.title = element_text(color = "black"))+
  RotatedAxis()
dev.off()

rna1<-data.frame(tcell@assays$RNA@data,check.rows = F,check.names = F)
rna1<-data.frame(t(rna1),check.rows = F,check.names = F)
rna1$origin_celltype<-meta$origin_celltype
rna1_origin_celltype<-rna1
m<-colSums(rna1_origin_celltype)
rna1_origin_celltype.use<-rna1_origin_celltype[,which(m!=0)]
rna1_origin_celltype.use<-rna1_origin_celltype.use[,c(DP_markers,naive_markers,Treg_markers,Gamma_delta_T_markers)]
pca3<-prcomp((rna1_origin_celltype.use))
plot3<-as.data.frame(pca3$x[,1:2])
plot3$origin<-rep(c("CON","CHI"),c(5,5))
plot3$sample<-rownames(plot3)
plot3$origin<-factor(plot3$origin,levels = c('CON','CHI'))
plot3<-plot3[c(1,6,2,7,3,8,4,9,5,10),]
plot3$celltype<-rep(c('DP','Gamma delta T','Naive CD4','Naive CD8','Treg'),
                    c(2,2,2,2,2))
percentVar3<-pca3$sdev^2/sum(pca3$sdev^2)
pdf("PCA.10origin_celltype.pdf",width = 5.6,height = 4)
ggplot(plot3,aes(x=PC1,y=PC2,color=celltype,shape=origin))+geom_point(size=5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12, color = "black"),axis.title.x=element_text(size=14))+
  theme(axis.text.y = element_text(size = 12, color = "black"),axis.title.y=element_text(size=14))+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=12))+
  theme(plot.title = element_text(size = 16))+
  labs(x = paste0("PC1(", round(percentVar3[1]*100,1),"%)"), y = paste0("PC2(", round(percentVar3[2]*100,1),"%)"), title = "")+
  scale_color_manual(name='Sample',values=c('#9FDAA4','#468141','#D3ECC2','#A98467','#9ABF6F'))
dev.off()

#correlation of T subsets
matrix<-data.frame(cor(data.frame(t(rna1_origin_celltype[,c(DP_markers,naive_markers,Treg_markers,Gamma_delta_T_markers)]),check.rows = F,check.names = F)))
pdf('T.origin_celltype.correlation.pdf',width = 5,height = 4.7)
pheatmap(matrix,cluster_cols = T,cluster_rows=T,show_colnames = T,show_rownames =T,
         angle_col = '90',
         border_color = NA,
         main = 'Correlation')
dev.off()

#proportion of T subsets in 6 samples#
meta$sample<-factor(meta$sample,levels = (c('CON_BM','CHI_BM','CON_LN','CHI_LN','CON_TM','CHI_TM')))
meta$celltype<-factor(meta$celltype,levels = rev(c('DP','Gamma delta T','Naive CD4','Naive CD8','Treg')))
pdf("T.subsets.sample.proportion.pdf",width = 8, height = 4)
ggplot(meta,aes(x=sample,fill=celltype))+
  geom_bar(stat = 'count',position = 'fill')+  
  labs(x = "",y = "Proportion", title = "",fill='T subsets')+
  theme_classic()+
  scale_fill_manual(values = rev(c('#9FDAA4','#468141','#D3ECC2','#A98467','#9ABF6F')))+
  #guides(fill=guide_legend(title=),size)+
  theme(axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        legend.text = element_text(color="black", size = 12))+
  theme(axis.text = element_text(colour = 'black'))+
  coord_flip()
dev.off()

