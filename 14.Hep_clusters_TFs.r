library(limma)
library(SCENIC)
library(AUCell)
library(data.table)
library(Seurat)
library(SCopeLoomR)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

loom <- open_loom("../data/work/pySCENIC/hep_out_SCENICn.loom")
rna<-readRDS("../SC_hep_harmony_knn5.rds")

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 

auc_mat <- as.matrix(regulonAUC@assays@data$AUC)
common_cells <- intersect(colnames(auc_mat), rownames(df_rna))
auc_mat <- auc_mat[, common_cells]
df_rna <- df_rna[common_cells, ]
cell_groups <- df_rna$anno_knn5
regulon_avg <- t(apply(auc_mat, 1, function(x) {
  tapply(x, cell_groups, mean, na.rm = TRUE)
}))
head(regulon_avg)

library(tibble)
top_tfs <- lapply(colnames(regulon_avg), function(g) {
  regulon_avg[, g] %>%
    sort(decreasing = TRUE) %>%
    head(20) %>%
    names()
}) %>% unlist() %>% unique()
regulon_top <- regulon_avg[top_tfs, ]

mtx_s<-t(scale(t(regulon_top)))
colors<-brewer.pal(name = "RdBu",n = 8)
cnn<-colorRampPalette(colors)(100)
p<-pheatmap::pheatmap(mtx_s,cluster_cols = F,border_color = NA,rev(cnn))
pdf("../hep_cluster_SigTFs.heatmap.pdf",width = 5,height = 9)
p
dev.off()

regulon_filter <- regulon_avg[apply(regulon_avg, 1, function(row){
  row["Hepatocyte 2"] == max(row)
}), ]
#Transcription‑factor regulators with high regulatory activity in Tra‑Hep cells were screened.
df_tg<-read.csv("../data/work/pySCENIC/hep_adj.samplen.tsv",sep = "\t")
df_slc<-df_tg[df_tg$TF %in% rownames(regulon_filter),]
df_tf_counts <- df_slc %>% 
  group_by(TF) %>% 
  summarise(gene_counts = n(), .groups = "drop") %>% 
  rename(TF_name = TF)
df_tf_counts<-df_tf_counts[order(df_tf_counts$gene_counts,decreasing = T),]
#Gene‑regulatory‑network plots were visualized using Cytoscape software (https://cytoscape.org/)





