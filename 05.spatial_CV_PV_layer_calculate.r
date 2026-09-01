library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table) 
library(Seurat)
library(dplyr)

sp_path<-Sys.glob("../*_rctd_sp.rds")
shot_path<-basename(sp_path)
sample_slc<-sub("_rctd_sp.*","",shot_path)

cv_genes<-c("Akr1c6","Aldh1a1","Car3","Ces1c","Cyp1a2","Cyp2c37","Cyp2c50","Cyp2e1","Gsta3","Mgst1","Mup11","Mup17","Mup18","Oat","Pon1","Rgn")
pv_genes<-c("Alb","Aldob","Arg1","Asl","Cyp2f2","Fbp1","Hal","Hpx","Hsd17b13","Mup20","Pck1","Slc25a47","Trf","Ass1")

AddCPscore<-function(rna){
    rna<-NormalizeData(rna)
   rna<-AddModuleScore(rna,features = list(cv_genes),name = "cv_score")
rna<-AddModuleScore(rna,features = list(pv_genes),name = "pv_score")
rna$z_score<-rna$pv_score1-rna$cv_score1
rna@meta.data$rank_cp=cut_number(rna@meta.data$pv_score1 - rna@meta.data$cv_score1,n =5,label=1:5) 
    return(rna)
}

for (i in c(1:length(sample_slc))){
sp_rna<-readRDS(sp_path[i])
sp_rna<-AddCPscore(sp_rna)
saveRDS(sp_rna,paste0("../",sample_slc[i],"_rctd_cp_sp.rds"))
}


