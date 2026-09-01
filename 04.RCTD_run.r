library("spacexr")
library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

gem_path<-Sys.glob("F:/stereo_gem/*.gem.gz")
shot_path<-basename(gem_path)
sample_name<-sub("\\.gem\\.gz","",shot_path)

gem2rds<-function(gem_path){
gem=fread(gem_path)
gem$cell<-paste0(round(gem$x /50,digits = 0) ,"_",round(gem$y / 50 ,digits = 0))
gene=1:length(unique(gem$geneID))
names(gene)=unique(gem$geneID)
cell=1:length(unique(gem$cell))
names(cell)=unique(gem$cell)
mat1=sparseMatrix(i = gene[gem$geneID],j=cell[ as.character(gem$cell) ], x=gem$MIDCount)
rownames(mat1)=names(gene)
colnames(mat1)=names(cell)

obj1=CreateSeuratObject(counts = mat1)
obj1@meta.data$x<-sub("_.*","",rownames(obj1@meta.data))
obj1@meta.data$y<-sub(".*_","",rownames(obj1@meta.data))
return(obj1)    
}

for  (i in c(1:length(gem_path))){
rna <-gem2rds(gem_path[i])
    assign(sample_name[i],rna)
    saveRDS(get(sample_name[i]),paste0("F:/stereo_gem/stereo_rds/",sample_name[i],"_bin50.orig.rds"))
}         
#The raw Stereo‑seq spatial transcriptomic GEM file was converted into a Seurat object.

rna_all<-readRDS("../SC_hep_harmony_knn5.rds")
rna_npc<-readRDS("../rna_npc.rds")   

rna_all$annotation<-rna_all$label_unify
rna_merge<-merge(x = rna_all,y = rna_npc)

sc <- rna_merge    #The HEPs and NPCs libraries were merged to build the reference index for RCTD.
rm(rna_merge)

sp_path<-Sys.glob("F:/stereo_gem/stereo_rds/*_bin50.orig.rds")
shot_path<-basename(sp_path)
sample_slc<-sub("_bin50.*","",shot_path)

for (i in c(1:length(sp_path))){
spatial<-readRDS(sp_path[i])
    message("readRDS done")
genes<-rownames(spatial@assays$RNA@counts)[rownames(spatial@assays$RNA@counts) %in% rownames(sc@assays$RNA@counts)]

DefaultAssay(sc)="RNA"
sc$celltype=as.character(sc$annotation)
Idents(sc)=as.character(sc$celltype)

celltype_sc=sc$celltype
celltype_sc=as.factor(celltype_sc)
nUMI_sc=sc@meta.data[,"nCount_RNA"]
names(nUMI_sc)=rownames(sc@meta.data)

exp_sc <- GetAssayData(sc, slot = "counts")
exp_sc<-exp_sc[rownames(exp_sc) %in% genes,]
reference <- Reference(exp_sc, celltype_sc, nUMI_sc)
message("ref conduct")

coord_spatial=spatial@meta.data[,c("x","y")]
coord_spatial$x<-as.numeric(coord_spatial$x)
coord_spatial$y<-as.numeric(coord_spatial$y)
nUMI_spatial=spatial@meta.data[,"nCount_RNA"]
names(nUMI_spatial)=rownames(spatial@meta.data)
exp_spatial <- GetAssayData(spatial, slot = "counts")
exp_spatial<-exp_spatial[rownames(exp_spatial) %in% genes,]
puck <- SpatialRNA(coord_spatial, exp_spatial, nUMI_spatial)

rm(spatial)
myRCTD <- create.RCTD(puck, reference, max_cores = 4)
    message("rctd start")
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
saveRDS(myRCTD,paste0("F:/stereo_gem/RCTD260105/"sample_slc[i],".anno_hepRCTD_20260105.rds"))
    message(paste0(sample_name[i],"done"))
}
#The RCTD pipeline was run and the resulting files were saved.

for (i in c(1:length(sample_slc))){
rctd<-readRDS(paste0("F:/stereo_gem/RCTD260105/"sample_slc[i],".anno_hepRCTD_20260105.rds"))
sp_rna<-readRDS("F:/stereo_gem/stereo_rds/",sample_slc[i],"_bin50.orig.rds")
weight=lapply(rctd@results,function(x){
    as.vector(x$all_weights)})
    weight<-Reduce(rbind,weight)
weight<-as.data.frame(weight)
rownames(weight)<-rownames(rctd@spatialRNA@coords)
colnames(weight)<-rctd@cell_type_info$info[[2]]
colnames(weight)<- sub(" ","_",colnames(weight))
sp_rna<-AddMetaData(sp_rna,weight)
saveRDS(sp_rna,paste0("../",sample_slc[i],"_rctd_sp.rds"))
}
#The weight matrix derived from RCTD computations was added to the meta.data of the spatial transcriptomic dataset.
