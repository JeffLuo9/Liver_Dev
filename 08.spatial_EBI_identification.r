library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table) 
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(circlize)
library(org.Mm.eg.db)
library(clusterProfiler)

rna_1<-readRDS("../Stereo_seq/new_stereo_filter_EK_niche/FP200002512_G2_rctd_sp.rds")  #P0
rna_2<-readRDS("../Stereo_seq/new_stereo_filter_EK_niche/FP200002512_G5_rctd_sp.rds")   #p3
rna_3<-readRDS("../Stereo_seq/new_stereo_filter_EK_niche/FP200000525TR_B2_rctd_sp.rds")   #p7
rna_4<-readRDS("../Stereo_seq/new_stereo_filter_EK_niche/SS200000158BR_E5_rctd_sp.rds")   #p14


for (i in c(1:4)){
  rna<-get(paste0("rna_",i))
  rna$Ery_cluster<-kmeans(rna$Erythroblast,centers = 2)$cluster
  rna$Kup_cluster<-kmeans(rna$Kupffer,centers = 2)$cluster
  assign(paste0("rna_",i),rna)
}
k_Ery<-function(rna){
    num_ery1<-mean(rna$Erythroblast[rna$Ery_cluster  == 1])
    num_ery2<-mean(rna$Erythroblast[rna$Ery_cluster  == 2])
    rna$Ery_pos<-rna$Ery_cluster
    ifelse(num_ery1 > num_ery2,rna$Ery_pos[rna$Ery_pos == 1]<-"Ery",rna$Ery_pos[rna$Ery_pos == 2]<-"Ery")
    return(rna)
}
k_Kup<-function(rna){
    num_kup1<-mean(rna$Kupffer[rna$Kup_cluster  == 1])
    num_kup2<-mean(rna$Kupffer[rna$Kup_cluster  == 2])
    rna$Kup_pos<-rna$Kup_cluster
    ifelse(num_kup1 > num_kup2,rna$Kup_pos[rna$Kup_pos == 1]<-"Kup",rna$Kup_pos[rna$Kup_pos == 2]<-"Kup")
    return(rna)
}
for (i in c(1:4)){
  rna<-get(paste0("rna_",i))
  rna<-k_Ery(rna)
  rna<-k_Kup(rna)
  rna$EBI_niche<-"other"
  rna@meta.data$EBI_niche[rna$Kup_cluster == "Kup" & rna$Ery_cluster == "Ery"]<-"niche"
  assign(paste0("rna_",i),rna)
}
saveRDS(rna_1,"../Stereo_seq/new_stereo_filter_EK_niche/FP200002512_G2_rctd_sp.rds")
saveRDS(rna_2,"../Stereo_seq/new_stereo_filter_EK_niche/FP200002512_G5_rctd_sp.rds")
saveRDS(rna_3,"../Stereo_seq/new_stereo_filter_EK_niche/FP200000525TR_B2_rctd_sp.rds")
saveRDS(rna_4,"../Stereo_seq/new_stereo_filter_EK_niche/SS200000158BR_E5_rctd_sp.rds")

df_EBI_genes<-read.csv("../EBI_DEGs.csv",row.names = 1)
EBI_genes<-rownames(df_EBI_genes)[df_EBI_genes$significant == "yes" & df_EBI_genes$log2.fold_change. > 0]

rna_m<-merge(x = rna_1,y = c(rna_2,rna_3,rna_4))
rna_m<-NormalizeData(rna_m)
rna_m<-AddModuleScore(rna_m,features = list(EBI_genes),name = "EBI_s")

p <- ggplot(meta_long, aes(x = time, y = exp, fill = E_c)) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.85,
    color = "black",
    linewidth = 0.6,
    outlier.shape = 21,
    outlier.size = 1,
    outlier.alpha = 0.4,
    outlier.fill = "black",
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = supplement_colors) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Expression Level",
    title = genes[1]
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )
pdf("../EBI_score_bt_sp_time.pdf",width = 5,height = 4)
p
dev.off()







