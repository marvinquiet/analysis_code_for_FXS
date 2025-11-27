## conda activate Seurat

library(ArchR)
set.seed(2022)
addArchRThreads(threads=1)
addArchRGenome("mm10")

project_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"
## get inputFiles
inputFiles = c("cortex_WT1"=file.path(project_dir, "Cortex_multiome_WT1_counts/outs/atac_fragments.tsv.gz"),
               "cortex_WT2"=file.path(project_dir, "Cortex_multiome_WT2_counts/outs/atac_fragments.tsv.gz"),
               "cortex_FXTAS1"=file.path(project_dir, "Cortex_multiome_FXTAS1_counts/outs/atac_fragments.tsv.gz"),
               "cortex_FXTAS2"=file.path(project_dir, "Cortex_multiome_FXTAS2_counts/outs/atac_fragments.tsv.gz"))

## create ArrowFiles indicated by ArchR
ArrowFiles = createArrowFiles(
        inputFiles = inputFiles,
        sampleNames = names(inputFiles),
        filterTSS = 4,
        filterFrags = 1000, 
        addTileMat = TRUE,
        addGeneScoreMat = TRUE
)

doubScores = addDoubletScores(
        input = ArrowFiles,
        k = 10,
        knnMethod = "UMAP", #
        LSIMethod = 1
)

proj = ArchRProject(
        ArrowFiles = ArrowFiles,
        outputDirectory = "Cortex_multiome_ArchR",
        copyArrows = TRUE
)

## plot TSS enrichment scores
p1 = plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
)
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "TSS_enrichment.png"))

p1 = plotFragmentSizes(ArchRProj = proj)
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "Fragsize_distribution.png"))
p2 = plotTSSEnrichment(ArchRProj = proj)
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "TSS_enrichment_distribution.png"))

saveArchRProject(ArchRProj = proj, outputDirectory = "Cortex_multiome_ArchR", load = FALSE)


## load ArchR Project
proj = loadArchRProject(path=file.path(project_dir, "Cortex_multiome_ArchR"))

## filter doublet
proj_filtered = filterDoublets(proj)
proj_filtered = addIterativeLSI(
    ArchRProj = proj_filtered,
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

## This does not seem to have batch effect.. we might not need the harmony
#proj_filtered = addHarmony(
#    ArchRProj = proj_filtered,
#    reducedDims = "IterativeLSI",
#    name = "Harmony",
#    groupBy = "Sample"
#)


proj_filtered <- addClusters(
    input = proj_filtered,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

proj_filtered = addUMAP(
    ArchRProj = proj_filtered, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
cluster_plot = plotEmbedding(ArchRProj = proj_filtered, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "Cluster_umap.png"))

sample_plot = plotEmbedding(ArchRProj = proj_filtered, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "Sample_umap.png"))

## load cell types
pred_celltypes = read.csv(file.path(project_dir, "Cortex_multiome_ArchR", "dscATAC_predict_celltypes.csv"), header=T, row.names=1)
proj_filtered$pred_celltype = pred_celltypes[rownames(proj_filtered), "pred_celltype"]
pred_celltype_plot = plotEmbedding(ArchRProj = proj_filtered, colorBy = "cellColData", name = "pred_celltype", embedding = "UMAP")
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "dscATAC_celltyping_umap.png"))

pred_celltypes = read.csv(file.path(project_dir, "Cortex_multiome_ArchR", "atlas_predict_celltypes.csv"), header=T, row.names=1)
proj_filtered$pred_celltype = pred_celltypes[rownames(proj_filtered), "pred_celltype"]
pred_celltype_plot = plotEmbedding(ArchRProj = proj_filtered, colorBy = "cellColData", name = "pred_celltype", embedding = "UMAP")
ggsave(file.path(project_dir, "Cortex_multiome_ArchR", "atlas_celltyping_umap.png"))



## find sample-specific clusters
cM = confusionMatrix(paste0(proj_filtered$Clusters), paste0(proj_filtered$Sample))
library(pheatmap)
cM = cM / Matrix::rowSums(cM)
png(file.path(project_dir, "Cortex_multiome_ArchR", "Cluster_heatmap.png"))
p = pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
dev.off()

## get gene score matrices from ArchR project to do Cellcano annotation
#GS_mat = getMatrixFromProject(
#    ArchRProj = proj_filtered,
#    useMatrix = "GeneScoreMatrix",
#    useSeqnames = NULL,
#    verbose = TRUE,
#    binarize = FALSE,
#    threads = getArchRThreads(),
#    logFile = createLogFile("getMatrixFromProject")
#)
#GS_df = as.matrix(assays(GS_mat)[[1]])
#rownames(GS_df) = rowData(GS_mat)$name
#write.csv(GS_df, file.path(project_dir, "Cortex_multiome_ArchR", "FMR_Multiome_GS.csv"), quote=F)


## get marker genes
markersGS = getMarkerFeatures(
    ArchRProj = proj_filtered, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
## save GS for manual annotation
mean_markerGS = assays(markersGS)[['Mean']]
rownames(mean_markerGS) = rowData(markersGS)$name
write.csv(mean_markerGS, file.path(project_dir, "Cortex_multiome_ArchR", "FMR_scATAC_mean_genescores.csv"), quote=F)

FDR_markerGS = assays(markersGS)[['FDR']]
rownames(FDR_markerGS) = rowData(markersGS)$name
write.csv(FDR_markerGS, file.path(project_dir, "Cortex_multiome_ArchR", "FMR_scATAC_FDR_genescores.csv"), quote=F)

markerList = getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markergeneList = lapply(markerList, function(x) {x$name})
markergeneDf = data.frame(lapply(markergeneList, function(x) {
    x <- unlist(x)
    length(x) <- max(lengths(markergeneList))
    return(x)
}))
write.csv(markergeneDf, file.path(project_dir, "Cortex_multiome_ArchR", "FMR_scATAC_markergenes_FDR_0.01_log2FC_1.25.csv"), quote=F)

markerGenes = c('Snap25', 'Syt1', ## Neurons
                'Slc17a7', 'Satb2', 'rorb', 'Neurod2', ## excitatory neurons
                'Gad1', 'Gad2', 'Nxph1', ## inhibitory neurons
                'Gfap', 'Aqp4', 'Slc1a2', ## Astrocytes
                'Csf1r', 'Cd74', 'P2ry12', 'Ptprc', 'Tmem119', ## microglia
                'Mobp', 'Mbp', 'St18', 'Klk6', 'Slc5a11', ## Oligodendrocytes
                'Pdgfra', 'Cspg4', ## Oligo precursor
                'Flt1', 'Cldn5', 'Abcb1', 'Ebf1' ## endothelial
)


heatmapGS = markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

pdf(file.path(project_dir, "Cortex_multiome_ArchR", "Markers_heatmap.pdf"))
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
