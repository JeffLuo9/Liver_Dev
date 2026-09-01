


sc <- readRDS("../SC_scRNA_hep_npc.rds")
for (i in c(1:length(rna_path_s))){
spatial<-readRDS(rna_path_s[i])
    message("readRDS done")
genes<-rownames(spatial@assays$RNA@counts)[rownames(spatial@assays$RNA@counts) %in% rownames(sc@assays$RNA@counts)]

DefaultAssay(sc)="RNA"
sc$celltype=as.character(sc$anno_hep)
Idents(sc)=as.character(sc$celltype)

celltype_sc=sc$celltype
celltype_sc=as.factor(celltype_sc)
nUMI_sc=sc@meta.data[,"nCount_RNA"]
names(nUMI_sc)=rownames(sc@meta.data)


exp_sc <- GetAssayData(sc, slot = "counts")
exp_sc<-exp_sc[rownames(exp_sc) %in% genes,]
reference <- Reference(exp_sc, celltype_sc, nUMI_sc)
message("ref conduct")
#空间转录组操作
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
