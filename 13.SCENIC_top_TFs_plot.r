library(limma)
library(SCENIC)
library(AUCell)
library(data.table)
library(Seurat)
library(SCopeLoomR)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

rna<-readRDS(""../scRNA_HSCs_subset.rds"")
loom <- open_loom("../data/work/pySCENIC/HSC_out_SCENICn.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 

df_rna<-rna@meta.data
auc_mat <- as.matrix(regulonAUC@assays@data$AUC)
common_cells <- intersect(colnames(auc_mat), rownames(df_rna))
auc_mat <- auc_mat[, common_cells]
df_rna <- df_rna[common_cells, ]
cell_groups <- df_rna$time
regulon_avg <- t(apply(auc_mat, 1, function(x) {
  tapply(x, cell_groups, mean, na.rm = TRUE)
}))
head(regulon_avg)

library(tibble)
top_tfs <- lapply(colnames(regulon_avg), function(g) {
  regulon_avg[, g] %>%
    sort(decreasing = TRUE) %>%
    head(15) %>%
    names()
}) %>% unlist() %>% unique()
regulon_top <- regulon_avg[top_tfs, ]
#The top 15 transcription‑factor genes were extracted for each time point.

colors<-c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30")
df_tg<-read.csv("../data/work/pySCENIC/LSEC_adj.samplen.tsv",sep = "\t")
df_tgf<-df1[df_tg$target %in% c("Wnt2","Wnt9b","Rspo3","Dll4","Jag1"),]
slc_tf<-df_tgf$TF[df_tgf$importance > 0.5]
slc_tf<-paste0(slc_tf,"(+)")
df_reg<-regulon_top[rownames(regulon_top) %in% slc_tf,]
df_reg_s<-t(scale(t(df_reg)))

cnn<-colorRampPalette(colors)(100)
p<-pheatmap::pheatmap(reg_s,cluster_cols = F,border_color = NA,cnn)
pdf("../HSCs_WNT_NOTCH_TFs.pdf",width = 5,height = 8)
p
dev.off()


rna<-readRDS("../scRNA_LSECs_subset.rds")
loom <- open_loom("../data/work/pySCENIC/LSEC_out_SCENICn.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 

df_rna<-rna@meta.data
auc_mat <- as.matrix(regulonAUC@assays@data$AUC)
common_cells <- intersect(colnames(auc_mat), rownames(df_rna))
auc_mat <- auc_mat[, common_cells]
df_rna <- df_rna[common_cells, ]
cell_groups <- df_rna$time
regulon_avg <- t(apply(auc_mat, 1, function(x) {
  tapply(x, cell_groups, mean, na.rm = TRUE)
}))
head(regulon_avg)

library(tibble)
top_tfs <- lapply(colnames(regulon_avg), function(g) {
  regulon_avg[, g] %>%
    sort(decreasing = TRUE) %>%
    head(15) %>%
    names()
}) %>% unlist() %>% unique()
regulon_top <- regulon_avg[top_tfs, ]

df_tg<-read.csv("../data/work/pySCENIC/LSEC_adj.samplen.tsv",sep = "\t")
df_tgf<-df1[df_tg$target %in% c("Wnt2","Wnt9b","Rspo3","Dll4","Jag1"),]
slc_tf<-df_tgf$TF[df_tgf$importance > 0.5]
slc_tf<-paste0(slc_tf,"(+)")
df_reg<-regulon_top[rownames(regulon_top) %in% slc_tf,]
df_reg_s<-t(scale(t(df_reg)))

cnn<-colorRampPalette(colors)(100)
p<-pheatmap::pheatmap(reg_s,cluster_cols = F,border_color = NA,cnn)
pdf("../LSECs_WNT_NOTCH_TFs.pdf",width = 5,height = 8)
p
dev.off()
