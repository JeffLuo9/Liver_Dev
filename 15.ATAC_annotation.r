library("ArchR")
library(patchwork)
addArchRThreads(threads = 8)
addArchRGenome("mm10")
library(data.table)

input_path<-Sys.glob("../data/*fragments.tsv.gz")
sample_name<-basename(input_path)
sample_name<-gsub("fragments.*","",sample_name)
names(input_path)<-sample_name

ArrowFiles <- createArrowFiles(
  inputFiles = input_path,
  sampleNames = names(input_path),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

projI_L <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/mnt/gandan/luojiahui/ArchR_project/",
  copyArrows = FALSE 
)

projI_L$Time<-gsub("-.*|_.*|NPC|npc|HEP|Hep|NCP","",projI_L$Sample)
projI_L$Time<-gsub("P351","P35",projI_L$Time)
projI_L$Time<-gsub("P352","P35",projI_L$Time)

projI_L=filterDoublets(projI_L,filterRatio = 8)

projI_L <- addIterativeLSI(
    ArchRProj = projI_L,
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
projI_L <- addHarmony(
    ArchRProj = projI_L,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force =T
)
projI_L <- addClusters(
    input = projI_L,
    reducedDims = "Harmony", #IterativeLSI
    method = "Seurat",
    name = "HarmonyClusters",
    resolution = 0.8,
    force =T
)
projI_L <- addUMAP(
    ArchRProj = projI_L,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force =T
)
projI_L <- addImputeWeights(projI_L)

markerGenes  <- c(
"Cyp2c37","Cyp2c50","Gulo","Cyp1a2","Cyp2a5","Cyp2e1",#Central hep
"Cdh1","Aldh1b1","Cyp2f2","Hal",#Portal hep
"Vsig4","Clec4f","Cd5l","Cd163",#Kuffer
"Cxcr2","S100a8","S100a9",#Inf.macro
"Runx2","Siglech","Flt3",#DC
"Ccr2","Cx3cr1","Fcgr3",#Mono
"Lrat","Dcn","Des","Cygb",#Stellate
"Krt19","Krt7","Epcam",#Cholangio
"Pecam1","Cdh5","Clec4g","Dnase1l3","Aqp1",#Endo
"Cd3d","Nkg7",#NKs_TCs
"Cd79b","Cd19","Pax5",#B cell
"Hba-a1","Hba-a2","Mki67",
"Epor","Gata1"    #Ery
  )
p <- plotEmbedding(
    ArchRProj = projI_L, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projI_L)
)
p2 <- lapply(p, function(x){
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
pdf("../celltype_marker_featureplot.pdf",width = 21,height = 90)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()

p <- plotGroups(
    ArchRProj = projI_L, 
    groupBy = "HarmonyClusters", 
    colorBy = "cellColData", 
  name = c("TSSEnrichment","nFrags","DoubletScore"),
  plotAs = "violin",
  alpha = 0.4,
    baseSize = 10,
  addBoxPlot = TRUE,
)

p2 <- lapply(p, function(x){
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
pdf("../HarmonyClusters_Vlnplot_QC.pdf",width = 18,height = 80)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()

f_sample<-c("P10NPC_2_1","P10NPC_2_2","P10NPC_3_1","P10NPC_3_2","P21NPC_2_1","P21NPC_2_2","P28HEP_2_1","P28HEP_2_2","P21HEP_1_2L","P21HEP_1_3L",
                                                    "P21HEP_1_1L","P35Hep1-2","P21NPC_2_3")
f_clusters<-c("C3","C43","C26","C35","C24","C46")

projI_L<-projI_L[! projI_L@cellColData$Sample %in% f_sample]
projI_L<-projI_L[! projI_L@cellColData$HarmonyClusters %in% f_clusters]

projI_L@cellColData$anno1=as.vector(mapvalues(proJI_L@cellColData$HarmonyClusters,from = paste0("C",c(1:43)),to=c("Hep","Hep" , "Hep","Hep","Hep",
                                                                                                                              "Hep","Hep","Hep","Hep",
                                                                                         "Ery","Ery","Hep","Hep","Hep","Hep","Ery","Hep",
                                                                                         "Hep","Hep","Endo","Endo","DC","Monocle","Kupffer",
                                                                                            "DC","Kupffer","Kupffer","Stellate","Cholangiocyte",
                                                                                            "Cholangiocyte","Monocle","Monocle","DC","DC",
                                                                                            "Ery","Nk_T","B_cell","B_cell","Nk_T","Nk_T","Nk_T",
                                                                                            "B_cell","Nk_T")))

proj_hep<-subsetArchRProject(projI_L,cells =rownames(projI_L@cellColData)[which(projI_L$anno1 %in% c("Hep"))],outputDirectory = "ATAC_hep_0226",force = T)
proj_npc<-subsetArchRProject(projI_L,cells =rownames(projI_L@cellColData)[which(! projI_L$anno1 %in% c("Hep"))],outputDirectory = "ATAC_npc_0226",force = T)


