library("ArchR")
library(patchwork)
addArchRThreads(threads = 8)
addArchRGenome("mm10")
library(data.table)

proj1<-loadArchRProject("../ATAC_hep_0226")

proj1 <- addIterativeLSI(
    ArchRProj = proj1,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    force =T
)
proj1 <- addHarmony(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force =T
)
proj1 <- addClusters(
    input = proj1,
    reducedDims = "Harmony", #IterativeLSI
    method = "Seurat",
    name = "HarmonyClusters",
    resolution = 0.8,
    force =T
)
proj1 <- addUMAP(
    ArchRProj = proj1,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force =T
)
proj1 <- addImputeWeights(proj1)

p1 <- plotEmbedding(ArchRProj = proj1 , colorBy = "GeneScoreMatrix" , name = c("Afp","H19","Mki67","Nfia","Egfr","Zbtb20","Immp2l","Cyp2e1","Cyp2f2f") , embedding = "UMAPHarmony")
p2 <- lapply(p1, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
pdf("../HarmonyClusters_featurePlot.pdf",width = 21,height = 21)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()

proj1$anno_hep=as.vector(mapvalues(proj1$HarmonyClusters,from = paste0("C",c(1:13)),to=c("Hep1","Hep1","Hep4","Hep1","Hep1",
                                                                                       "Hep3","Hep3","Hep3","Hep4","Hep2",
                                                                                      "Hep2","Hep1","Hep4")))
markersGS <- getMarkerFeatures(
    ArchRProj = proj1, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "anno_hep",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "anno_hep")
pathToMacs2 <- findMacs2()

proj1 <- addReproduciblePeakSet(
    ArchRProj = proj1, 
    groupBy = "anno_hep", 
    pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj1 <- addPeakMatrix(proj1)
proj1<-saveArchRProject(proj1,outputDirectory = "../ATAC_hep_0226",overwrite = T,load = T)

proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")

markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "anno_hep",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Hep2"
)
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj1,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

proj1 <- addBgdPeaks(proj1)
proj1 <- addDeviationsMatrix(
  ArchRProj = proj1, 
  peakAnnotation = "Motif",
  force = TRUE
)
proj1 <- addCoAccessibility(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI"
)



