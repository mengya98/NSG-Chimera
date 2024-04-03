#Clustering cells of chimera mice into 129mTmG or NSG origin

library(dplyr)
library(Seurat)
library(ggplot2) 
library(reshape2)

setwd('G:/Chimera/RNA/SNP/CHI_BM')

bm.alt.data<-Read10X(data.dir = "./raw_allele_bc_matrices_mex/bm_alt/")
bm_alt<-data.frame(bm.alt.data)
bm.ref.data<-Read10X(data.dir = "./raw_allele_bc_matrices_mex/bm_ref/")
bm_ref<-data.frame(bm.ref.data)
snps_freq<-bm_alt/(bm_alt+bm_ref+0.1)

freq <- lapply(as.list(snps_freq), function(x) rownames(snps_freq)[which(x>0.5)])
snp129<-read.table('../129snp.txt')
rownames(snp129)<-paste0(snp129$V1,',',snp129$V2,',',snp129$V5,',',snp129$V7)
snp_loc_129<-paste0(snp129$V1,',',snp129$V2,',',snp129$V5,',',snp129$V7)
frq_129<-data.frame(lapply(freq, function(x) length(intersect(snp_loc_129,x))))

snpnod<-read.table('../nodsnp.txt')
rownames(snpnod)<-paste0(snpnod$V1,',',snpnod$V3,',',snpnod$V6,',',snpnod$V8)
snp_loc_nod<-paste0(snpnod$V1,',',snpnod$V3,',',snpnod$V6,',',snpnod$V8)
frq_nod<-data.frame(lapply(freq, function(x) length(intersect(snp_loc_nod,x))))

frq_all<-rbind(frq_129,frq_nod)
frq_all<-data.frame(t(frq_all))
for (i in 1:nrow(frq_all)){
  if(frq_all[i,1]>=frq_all[i,2]){
    frq_all[i,3]<-"129"
    if(frq_all[i,2]==0){
      frq_all[i,4]<-frq_all[i,1]/(frq_all[i,2]+1)
    }
    if(frq_all[i,2]>0){
      frq_all[i,4]<-frq_all[i,1]/frq_all[i,2]
    }
  }
  if(frq_all[i,1]<frq_all[i,2]){
    frq_all[i,3]<-"nod"
    frq_all[i,4]<-frq_all[i,1]/frq_all[i,2]
  }
}
tsne<-read.table('G:/Chimera/RNA/CHI_BM/tsne.info.txt',sep='\t',header = T,row.names = 1)
tsne$snp<-frq_all$V3
tsne$ratio<-frq_all$V4
write.table(data.frame(rownames(tsne),tsne),file='CHI_BM.ratio.xls',row.names = F, col.names = T, quote = F, sep = "\t")


setwd('G:/Chimera/RNA/SNP/CHI_LN')
ln.alt.data<-Read10X(data.dir = "./raw_allele_bc_matrices_mex/ln_alt/")
ln_alt<-data.frame(ln.alt.data)
ln.ref.data<-Read10X(data.dir = "./raw_allele_bc_matrices_mex/ln_ref/")
ln_ref<-data.frame(ln.ref.data)
snps_freq<-ln_alt/(ln_alt+ln_ref+0.1)

freq <- lapply(as.list(snps_freq), function(x) rownames(snps_freq)[which(x>0.5)])
frq_129<-data.frame(lapply(freq, function(x) length(intersect(snp_loc_129,x))))
frq_nod<-data.frame(lapply(freq, function(x) length(intersect(snp_loc_nod,x))))

frq_all<-rbind(frq_129,frq_nod)
frq_all<-data.frame(t(frq_all))
for (i in 1:nrow(frq_all)){
  if(frq_all[i,1]>=frq_all[i,2]){
    frq_all[i,3]<-"129"
    if(frq_all[i,2]==0){
      frq_all[i,4]<-frq_all[i,1]/(frq_all[i,2]+1)
    }
    if(frq_all[i,2]>0){
      frq_all[i,4]<-frq_all[i,1]/frq_all[i,2]
    }
  }
  if(frq_all[i,1]<frq_all[i,2]){
    frq_all[i,3]<-"nod"
    frq_all[i,4]<-frq_all[i,1]/frq_all[i,2]
  }
}
tsne<-read.table('G:/Chimeria/RNA/CHI_LN/tsne.info.txt',sep='\t',header = T,row.names = 1)
tsne$snp<-frq_all$V3
tsne$ratio<-frq_all$V4
write.table(data.frame(rownames(tsne),tsne),file='ln.ratio.xls',row.names = F, col.names = T, quote = F, sep = "\t")


setwd('G:/Chimera/RNA/SNP/CHI_TM')
tm.alt.data<-Read10X(data.dir = "./raw_allele_bc_matrices_mex/tm_alt/")
tm_alt<-data.frame(tm.alt.data)
tm.ref.data<-Read10X(data.dir = "./raw_allele_bc_matrices_mex/tm_ref/")
tm_ref<-data.frame(tm.ref.data)
snps_freq<-tm_alt/(tm_alt+tm_ref+0.1)

freq <- lapply(as.list(snps_freq), function(x) rownames(snps_freq)[which(x>0.5)])
frq_129<-data.frame(lapply(freq, function(x) length(intersect(snp_loc_129,x))))
frq_nod<-data.frame(lapply(freq, function(x) length(intersect(snp_loc_nod,x))))

frq_all<-rbind(frq_129,frq_nod)
frq_all<-data.frame(t(frq_all))
for (i in 1:nrow(frq_all)){
  if(frq_all[i,1]>=frq_all[i,2]){
    frq_all[i,3]<-"129"
    if(frq_all[i,2]==0){
      frq_all[i,4]<-frq_all[i,1]/(frq_all[i,2]+1)
    }
    if(frq_all[i,2]>0){
      frq_all[i,4]<-frq_all[i,1]/frq_all[i,2]
    }
  }
  if(frq_all[i,1]<frq_all[i,2]){
    frq_all[i,3]<-"nod"
    frq_all[i,4]<-frq_all[i,1]/frq_all[i,2]
  }
}
tsne<-read.table('G:/Chimerisma/RNA/CHI_TM/tsne.info.txt',sep='\t',header = T,row.names = 1)
tsne$snp<-frq_all$V3
tsne$ratio<-frq_all$V4
write.table(data.frame(rownames(tsne),tsne),file='tm.ratio.xls',row.names = F, col.names = T, quote = F, sep = "\t")
