library("data.table")
library(Seurat)

rna_hep<-readRDS("../SC_hep_harmony_knn5.rds")
rna_hep<-subset(x = rna_hep, downsample = 50000)
write.csv(t(as.matrix(rna_hep@assays$RNA@counts)),file = "../data/work/Hep_for.scenic.data.csv")

rna_EC<-readRDS("../scRNA_LSECs_subset.rds")
write.csv(t(as.matrix(rna_EC@assays$RNA@counts)),file = "../data/work/EC_for.scenic.data.csv")

rna_HSC<-readRDS("../scRNA_HSCs_subset.rds")
write.csv(t(as.matrix(rna_HSC@assays$RNA@counts)),file = "../data/work/HSC_for.scenic.data.csv")

