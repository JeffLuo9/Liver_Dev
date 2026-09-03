library(limma)
library(SCENIC)
library(AUCell)
library(data.table)
library(Seurat)
library(SCopeLoomR)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

rna<-readRDS("/data/work/pySCENIC_NPC0520/HSC1_0427_1.rds")
loom <- open_loom("/data/work/pySCENIC_NPC0520/HSC_out_SCENICn.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 

rna<-readRDS("/data/work/pySCENIC_NPC0520/EC1_0427_1.rds")
loom <- open_loom("/data/work/pySCENIC_NPC0520/LSEC_out_SCENICn.loom")
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


