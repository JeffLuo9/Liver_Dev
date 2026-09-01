library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(clusterProfiler)

df_lr<-read.csv("../mouse_lr_pair_from_celltalk.csv",row.names = 1)
gene_l<-unique(df_lr$ligand_gene_symbol)
gene_r<-unique(df_lr$receptor_gene_symbol)

mtx<-AverageExpression(sp_rna,features = c(gene_r,gene_l),group.by = "rank_cp")$RNA
exp_mat<-mtx
lr_product <- t(apply(df_lr, 1, function(row) {
  lig <- row["ligand_gene_symbol"]
  rec <- row["receptor_gene_symbol"]
  if (lig %in% rownames(exp_mat) & rec %in% rownames(exp_mat)) {
    exp_mat[lig, ] * exp_mat[rec, ]
  } else {
    rep(NA, 5)
  }
}))
lr_product <- as.data.frame(lr_product)
colnames(lr_product) <- paste0("Layer_", 1:5)
rownames(lr_product) <- paste(df_lr$ligand_gene_symbol,
                              df_lr$receptor_gene_symbol,
                              sep = "_")
lr_product <- lr_product[!apply(lr_product, 1, function(x) all(is.na(x))), ]

lr_s<-t(scale(t(lr_product)))

lr_df_s<-as.data.frame(lr_s)

result_pv <- lr_df_s %>%
  filter(Layer_2 > Layer_1 &
         Layer_3 > Layer_2 &
         Layer_4 > Layer_3 &
         Layer_5 > Layer_4)
result_cv <- lr_df_s %>%
  filter(Layer_2 < Layer_1 &
         Layer_3 < Layer_2 &
         Layer_4 < Layer_3 &
         Layer_5 < Layer_4)

                                


                                


