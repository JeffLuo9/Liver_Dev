library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(RColorBrewer)
library(RANN)
library(ggpubr)

sp_rna<-readRDS("../SS200000158BR_E5_rctd_sp.rds")

smooth_kNN <- function(initial_df,query_df,sm_vector,round=100,knn=30){
    result <- nn2(initial_df,query_df,k=knn)
    for(i in 1:round){
        if(i==1){initail_col <- sm_vector}else{initail_col <- sm_vector_sm}
        sm_vector_sm <- rowMeans(matrix(initail_col[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE))
    }
    return(sm_vector_sm)
}

df <- sp_rna@meta.data
sp_rna$CP_score_sm <- smooth_kNN(df[,c('x','y')],df[,c('x','y')], df[,'CP_score'], knn=20,round=3)
sp_rna$layer=cut_number(sp_rna$CP_score_sm,n = 5)
sp_rna$layer=as.numeric(factor(sp_rna$layer))
#Contour lines were smoothed.

df_cvlr<-read.csv("../P14_cv_zonation_LR.csv",row.names = 1)

for (n in c(1:nrow(df_cvlr))){  #Ligand‑receptor pairs with gradient‑changing CV characteristics were iteratively plotted.
gene_n<-strsplit(rownames(df_cvlr)[n],split = "_",fixed = T)[[1]]
sp_rna@meta.data[,gene_n[1]]<-sp_rna@assays$RNA@data[gene_n[1],]
sp_rna@meta.data[,gene_n[2]]<-sp_rna@assays$RNA@data[gene_n[2],] 
df <- sp_rna@meta.data
sp_rna@meta.data[,paste0(gene_n[1],"_sm")] <- smooth_kNN(df[,c('x','y')],df[,c('x','y')], df[,gene_n[1]], knn=20,round=3)
sp_rna@meta.data[,paste0(gene_n[2],"_sm")] <- smooth_kNN(df[,c('x','y')],df[,c('x','y')], df[,gene_n[2]], knn=20,round=3)
#Gene‑expression values were smoothed using the smooth_KNN function.
sp_rna@meta.data[,rownames(df_cvlr)[n]]<-sp_rna@meta.data[,paste0(gene_n[1],"_sm")] * sp_rna@meta.data[,paste0(gene_n[2],"_sm")]

options(repr.plot.width=8,repr.plot.height=8)
obj1<-sp_rna
cp <- paste0(rownames(df_cvlr)[n])
tmp <- obj1@meta.data[,c("x","y","layer",ng)]
names(tmp) <- c("x","y","z","gene")
p1 <- ggplot()+geom_tile(data = tmp,aes(x=x,y=y),fill = "#808080",size=0.1)+
    geom_tile(data = tmp %>% arrange(gene),aes(x=x,y=y,fill=gene),alpha=1,size=0.1)+
    scale_fill_viridis_c(option="B",direction = -1)+labs(fill=cp)+
    theme_void()+coord_fixed()+geom_contour(data = tmp,mapping = aes(x,y,z= z),breaks=3.05,color="white",alpha=1)
print(p1)

p1n<-p1 + theme(legend.position = "none")
png(paste0("../P14_",rownames(df_cvlr)[n],"_sp_Plot.no_leg.png"), width=2200, height=2200, res=500,bg = "transparent")
print(p1n)
dev.off()

leg <- get_legend(p1 + theme(legend.position="right"))
pdf("../P14_",rownames(df_cvlr)[n],"legend.pdf", width=2, height=3)
print(as_ggplot(leg))
dev.off()
}











