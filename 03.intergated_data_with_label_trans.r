library("data.table")
library(Seurat)
library("harmony")
library(RANN)

rna_sc<-readRDS("../Stereo_cell_data_with_annotation.rds")  
#Clustering analysis and cell‑type annotation were performed on the Stereo‑cell dataset, identifying four distinct hepatocyte populations.

rna_all<-readRDS("../scRNA_and_Stereo_cell_harmony_UMAP.rds")
#This dataset was obtained by merging the Stereo‑cell dataset and the single‑cell C4 hepatocyte dataset.

rna_all@assays$RNA@var.features<-rna_sc@assays$RNA@var.features
Harmony_run<-function(rna_all,dims=1:15,features=3000){
  rna_all<-NormalizeData(rna_all)
  rna_all <- ScaleData(rna_all)
  rna_all<-RunPCA(rna_all)
  rna_all<-RunHarmony(rna_all,"sample")
  rna_all <- rna_all %>%
  FindNeighbors(reduction = "harmony",dims = dims) %>%
  FindClusters(resolution = 0.9) %>%
  RunUMAP(reduction = "harmony",dims = dims) 
  return(rna_all)
}
rna_all<-Harmony_run(rna_all)  #The Stereo‑cell and C4 scRNA‑seq datasets were integrated using Harmony.

winner_kNN <- function(initial_df,query_df,sm_vector,knn=5,threhold=1){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    sm_vector <- matrix(sm_vector[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE)
    sm_vector <- sapply(1:nrow(sm_vector),function(x){
        tmp_table <- table(sm_vector[x,])
        name <- names(sort(tmp_table,decreasing = T))[1]
        if(tmp_table[name]>threhold){return(name)}else{return(sm_vector[x,1])}
    })
    return(sm_vector)
}

rna_all$label_unify <- rna_all$anno_hep    #The annotation column in the Stereo‑Cell dataset is labelled "anno_hep"
rna_all@meta.data[is.na(rna_all$label_unify),'label_unify'] <- winner_kNN(
                                       initial_df = rna_all@reductions$harmony@cell.embeddings[!is.na(rna_all$label_unify),1:15],
                                       query_df = rna_all@reductions$harmony@cell.embeddings[is.na(rna_all$label_unify),1:15],
                                        sm_vector = rna_all@meta.data[!is.na(rna_all$label_unify),'label_unify'],
                                         knn = 5)


