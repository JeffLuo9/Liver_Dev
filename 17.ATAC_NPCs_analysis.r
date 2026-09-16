library("ArchR")
library(patchwork)
addArchRThreads(threads = 8)
addArchRGenome("mm10")
library(data.table)
library(ggseqlogo)

#Due to the poor cross‑platform portability of ArchR project files, we regenerated the ArchR project after changing the analysis directory path.

df_name<-read.csv("../All_cellColData.csv",row.names =1 )
LSEC_names<-rownames(df_name)[df_name$annotation == "Endo"]


ArrowFile1<-Sys.glob("../data/work/P0_P28_HEP/hep_P0P280129/ArrowFiles/*.arrow")
ArrowFile2<-Sys.glob("../data/work/P0_P28_HEP/new_arrowFiles/ArrowFiles/*.arrow")

proj <- ArchRProject(
  ArrowFiles = c(ArrowFile1,ArrowFile2),
  outputDirectory = "../data/work/LSEC_HSC_arrowfile/new_all_atac/",
  copyArrows = FALSE 
)

proj1<-proj[rownames(proj@cellColData) %in% LSEC_names]

proj1$Time<-gsub("ATAC-|ATAC-XJS-|_.*|NPC|npc|HEP|Hep|NCP","",proj1$Sample)
proj1$Time<-gsub("-.*","",proj1$Time)
proj1$Time<-gsub("P351","P35",proj1$Time)
proj1$Time<-gsub("P352","P35",proj1$Time)
proj1$Time<-gsub("p","P",proj1$Time)

proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Time",force = TRUE)
pathToMacs2 <- findMacs2()

proj1 <- addReproduciblePeakSet(plot = FALSE,
    ArchRProj = proj1, 
    groupBy = "Time", 
    pathToMacs2 = pathToMacs2,
  force = TRUE,
)

proj1 <- addPeakMatrix(proj1,force = TRUE)
proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")

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
    dimsToUse = 1:30
)

proj1 <- addCoAccessibility(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI"
)

proj1 <- addBgdPeaks(proj1,force = TRUE)

proj1 <- addDeviationsMatrix(
  ArchRProj = proj1, 
  peakAnnotation = "Motif",
  force = TRUE
)

proj1<-saveArchRProject(proj1,outputDirectory = "../data/work/LSEC_HSC_arrowfile/Erythrobalst/",overwrite = T,load = T)
rna_npc<-readRDS("../rna_npc.rds") 
rna_npc<-subset(rna_npc,subset = annotation == "LSEC")

groupList <- SimpleList(
    P0 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P0"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P0"]
    ),
    P3 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P3"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P3"]
    ),
     P7 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P7"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P7"]
    ),
     P10 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P10"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P10"]
    ),
     P14 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P14"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P14"]
    ),
     P21 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P21"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P21"]
    ),
     P28 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P28"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P28"]
    ),
     P35 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P35"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P35"]
    ),
     P42 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P42"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P42"]
    ),
     P56 = SimpleList(
        ATAC =rownames(proj1@cellColData)[proj1$Time == "P56"],
        RNA = rownames(rna_npc@meta.data)[rna_npc$time == "P56"]
    )
    )    
proj1 <- addGeneIntegrationMatrix(
    ArchRProj = proj1, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = rna_npc,
    addToArrow = TRUE, 
    groupList = groupList,
    groupRNA = "time",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co",
    force = TRUE
)
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "Time",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "P7",
  bgdGroups = "P10"
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

ggUp

seGroupMotif <- getGroupSE(ArchRProj = proj1, useMatrix = "MotifMatrix", groupBy = "Time")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
    ArchRProj = proj1,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.6& corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )+ggrepel::geom_label_repel(
        data = as.data.frame(corGSM_MM[corGSM_MM$TFRegulator == "YES",]), aes(x = cor, y = maxDelta, label = MotifMatrix_name), 
        size = 1.5,
        nudge_x = 0.2,
        max.overlaps = 25, 
        color = "black"
  )

p

pdf("../data/work/HEP/LSEC_GSM_MM_order_plot.pdf",width = 7.5,height = 6)
p
dev.off()

proj1 <- addPeak2GeneLinks(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj1,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
p2g

markerGenes  <- c("Klf2","Klf4"
  )

p <- plotBrowserTrack(
    ArchRProj = proj1, 
    groupBy = "Time", 
    geneSymbol = markerGenes, 
    upstream = 30000,
    downstream = 15000,
    loops = getPeak2GeneLinks(proj1)
)
#VlnPlot
p2 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Time", 
    colorBy = "GeneScoreMatrix", 
  name ="Klf2",
  plotAs = "violin",
  alpha = 0.7,
    baseSize = 10,
  addBoxPlot = F,
)

p2 <- p2 +

  geom_boxplot(width = 0.3, fill = "white", colour = "black", linewidth = 0.5, outlier.shape = NA) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

p2
#Erythroblast, HSC and LSEC were processed using the same analytical workflow.






