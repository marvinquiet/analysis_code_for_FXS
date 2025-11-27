## conda activate Seurat
library(Seurat)
library(Matrix)
library(ggplot2)
set.seed(2022)

project_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"
## get inputFiles
inputFiles = c("cortex_WT1"="Cortex_multiome_WT1_counts",
               "cortex_WT2"="Cortex_multiome_WT2_counts",
               "cortex_FXTAS1"="Cortex_multiome_FXTAS1_counts",
               "cortex_FXTAS2"="Cortex_multiome_FXTAS2_counts")

sample_obj_list = list()
for (sample in names(inputFiles)) {
    input_dir = file.path(project_dir, inputFiles[sample], 'outs', 'filtered_feature_bc_matrix')
    input_mat = readMM(file.path(input_dir, 'matrix.mtx.gz'))
    input_barcodes = scan(file.path(input_dir, 'barcodes.tsv.gz'), what=character())
    input_features = read.csv(file.path(input_dir, 'features.tsv.gz'), header=F, sep='\t')

    ## filter gene count matrix and create Seurat matrix
    gene_expr_id = which(input_features$V3 == 'Gene Expression')
    gene_expr_mat = input_mat[gene_expr_id, ]
    rownames(gene_expr_mat) = input_features[gene_expr_id, 'V2']
    colnames(gene_expr_mat) = input_barcodes

    ## create Seurat object
    sample_obj = CreateSeuratObject(counts=gene_expr_mat, project=sample, 
                                    min.cells=0, min.features=0)
    sample_obj_list[[sample]] = sample_obj
}
sample_combined_obj = merge(sample_obj_list[[1]], 
                     y=c(sample_obj_list[[2]], sample_obj_list[[3]], sample_obj_list[[4]]), 
                     add.cell.ids=names(sample_obj_list))
sample_combined_obj[["percent.mt"]] = PercentageFeatureSet(sample_combined_obj, pattern = "^mt-")  ## already filtered

result_dir = file.path(project_dir, 'Cortex_multiome_Seurat_pFC_200_10000')
dir.create(result_dir, showWarnings = FALSE)

VlnPlot(sample_combined_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(result_dir, 'QuanlityViolinplot.png'))
VlnPlot(sample_combined_obj, features='nFeature_RNA') & geom_hline(yintercept=200) & geom_hline(yintercept=10000)
ggsave(file.path(result_dir, 'QuanlityViolinplot_nFeatureRNA.png'))
VlnPlot(sample_combined_obj, features='percent.mt') & geom_hline(yintercept=5)
ggsave(file.path(result_dir, 'QuanlityViolinplot_mt-percent.png'))

nFeature_cutoff = 200  ## set 200 to keep more available cells
result_dir = file.path(project_dir, 'Cortex_multiome_Seurat_pFC_200_10000')
dir.create(result_dir, showWarnings = FALSE)

combined_obj = subset(sample_combined_obj, subset = nFeature_RNA > nFeature_cutoff & nFeature_RNA < 10000)
## subset with ArchR scATAC-seq data
GS_df = read.csv(file.path(project_dir, 'Cortex_multiome_ArchR', 'FMR_Multiome_GS.csv'), header=T, row.names=1, sep=',', check.names=F)
scATAC_id = gsub('#', '_', colnames(GS_df))
common_cells = intersect(colnames(combined_obj), scATAC_id)
combined_obj = combined_obj[, common_cells]
print(combined_obj)

combined_obj = NormalizeData(combined_obj)
combined_obj = FindVariableFeatures(combined_obj, selection.method = "vst", nfeatures = 2000)
combined_obj = ScaleData(combined_obj, features=rownames(combined_obj))
combined_obj = RunPCA(combined_obj, features = VariableFeatures(object=combined_obj))
g = ElbowPlot(combined_obj, ndims=50)
ggsave(file.path(result_dir, 'PC_elbowplot.png')) ## will select 1:40 PCs

combined_obj = FindNeighbors(combined_obj, dims = 1:40)
resolutions = c(0.3, 0.4, 0.5, 0.6)
for (resolution in resolutions) {
    combined_obj = FindClusters(combined_obj, resolution = resolution)
    combined_obj = RunUMAP(combined_obj, dims=1:40)
    DimPlot(combined_obj, reduction = "umap")
    ggsave(file.path(result_dir, paste0('UMAP_louvain', resolution, '.png')))
    DimPlot(combined_obj, reduction = "umap", group.by='orig.ident')
    ggsave(file.path(result_dir, paste0('UMAP_louvain', resolution, '_sample.png')))
}

## table of summary statistics
table(combined_obj@meta.data$orig.ident)
table(combined_obj@meta.data$orig.ident, Idents(combined_obj))

## boxplot of number of genes for each cluster and predicted cell types
pred_celltypes = read.csv(file.path(project_dir, 'Cortex_multiome_Seurat_lessstringent', 'mousepFC_lessstringent_predcelltypes.csv'), header=T, row.names=1)
rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
combined_obj@meta.data$mousepFC_pred_celltype = pred_celltypes[colnames(combined_obj), 'pred_celltype']

ggplot(combined_obj@meta.data, aes(x=nFeature_RNA, y=seurat_clusters, fill=seurat_clusters)) + geom_boxplot() + 
    geom_vline(xintercept=800, linetype="dashed", color='blue', size=1.5) + 
    geom_vline(xintercept=1500, linetype="dashed", color='red', size=1.5) 
ggsave(file.path(result_dir, 'nFeature_based_on_clusters.png'))

ggplot(combined_obj@meta.data, aes(x=nFeature_RNA, y=mousepFC_pred_celltype, fill=mousepFC_pred_celltype)) + geom_boxplot() + 
    geom_vline(xintercept=800, linetype="dashed", color='blue', size=1.5) + 
    geom_vline(xintercept=1500, linetype="dashed", color='red', size=1.5) 
ggsave(file.path(result_dir, 'nFeature_based_on_mousepFC_predcelltypes.png'))

DimPlot(combined_obj, reduction = "umap", group.by='mousepFC_pred_celltype')
ggsave(file.path(result_dir, 'UMAP_louvain0.5_Cellcano_pFC_majorcelltype.png'))

## === different cutoff for neuronal clusters
clusters = c(0,2,4,5,6,7,10,13,14,15,17,18,19,21,22,23) ## decided by nFeature_RNA per cluster
neuronal_cells = colnames(combined_obj)[Idents(combined_obj) %in% clusters]
neuronal_obj = combined_obj[, neuronal_cells]
filtered_neuronal_cells = colnames(neuronal_obj)[neuronal_obj@meta.data$nFeature_RNA >= 1000] ## change to 1000

## === filtered results
nonneuronal_cells = colnames(combined_obj)[!(Idents(combined_obj) %in% clusters)]
all_cells = c(nonneuronal_cells, filtered_neuronal_cells)
filtered_obj = combined_obj[, all_cells]
#writeMM(filtered_obj[['RNA']]@counts, file.path(result_dir, 'filtered_pFC_GE.mtx'))
#write(rownames(filtered_obj[['RNA']]@counts), file.path(result_dir, 'filtered_pFC_GE_genes.tsv'))
#write(colnames(filtered_obj[['RNA']]@counts), file.path(result_dir, 'filtered_pFC_GE_barcodes.tsv'))
## after writing out, use Cellcano to predict, the following are the command lines
# Cellcano predict -i filtered_pFC_GE --trained_model mousepFC_Cellcano_output/majorcelltypes_predMLP_model/ -o mousepFC_Cellcano_output --prefix filtered_pFC_GE
# Cellcano predict -i filtered_pFC_GE --trained_model Saunders_Cellcano_output/majorcelltypes_predMLP_model/ -o Saunders_Cellcano_output --prefix filtered_pFC_GE_Saunders

## === after filtering
filtered_obj = NormalizeData(filtered_obj)
filtered_obj = FindVariableFeatures(filtered_obj, selection.method = "vst", nfeatures = 2000)
filtered_obj = ScaleData(filtered_obj, features=rownames(filtered_obj))
filtered_obj = RunPCA(filtered_obj, features = VariableFeatures(object=filtered_obj))
g = ElbowPlot(filtered_obj, ndims=50)
ggsave(file.path(result_dir, 'filtered_PC_elbowplot.png')) ## will select 1:40 PCs

filtered_obj = FindNeighbors(filtered_obj, dims = 1:40)
filtered_obj = FindClusters(filtered_obj, resolution = 0.5)
filtered_obj = RunUMAP(filtered_obj, dims=1:40)
DimPlot(filtered_obj, reduction = "umap", label=T)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5.png'))
DimPlot(filtered_obj, reduction = "umap", group.by='orig.ident')
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_sample.png'))
filtered_obj@meta.data$condition = ifelse(grepl('WT', filtered_obj@meta.data$orig.ident), 'control', 'FXTAS2')
DimPlot(filtered_obj, reduction = "umap", split.by='condition', label=T)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_condition.png'), width=10, height=5)

## plot with mouse pFC as ref
pred_celltypes = read.csv(file.path(result_dir, 'filtered_pFC_GEcelltypes.csv'), header=T, row.names=1)
rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
filtered_obj@meta.data$mousepFC_pred_celltype = pred_celltypes[colnames(filtered_obj), 'pred_celltype']
DimPlot(filtered_obj, reduction = "umap", group.by='mousepFC_pred_celltype')
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_Cellcano_pFC_majorcelltype.png'))

## plot with Saunders as ref
pred_celltypes = read.csv(file.path(result_dir, 'filtered_pFC_GE_Saunderscelltypes.csv'), header=T, row.names=1)
rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
filtered_obj@meta.data$Saunders_pred_celltype = pred_celltypes[colnames(filtered_obj), 'pred_celltype']
DimPlot(filtered_obj, reduction = "umap", group.by='Saunders_pred_celltype')
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_Cellcano_pFC_Saunders_majorcelltype.png'))


## new marker genes available
### === marker genes provided by Dr. Yulin Jin
#markerGenes = c('Snap25',  ## Neurons
#                'Slc17a7', ## excitatory neurons
#                'Gad2',  ## inhibitory neurons
#                'Gja1',  ## Astrocytes
#                'C1qa',  ## microglia
#                'Aspa',  ## Oligodendrocytes
#                'Pdgfra',  ## Oligo precursor
#                'Flt1',  ## endothelial
#                'Bmp4' ## NF oligo
#)
#FeaturePlot(filtered_obj, features = markerGenes)
#ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_markergenes.png'))

## === curate cell types
cluster_celltype_table = table(Idents(filtered_obj), filtered_obj@meta.data$mousepFC_pred_celltype)
cluster_celltype_table = cluster_celltype_table/rowSums(cluster_celltype_table)
cluster_celltype_df = as.data.frame(cluster_celltype_table)
colnames(cluster_celltype_df) = c('cluster', 'celltype', 'prop')

ggplot(cluster_celltype_df, aes(x=celltype, y=cluster, fill=prop)) + geom_tile() + 
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'pFC_cluster_vs_celltype_heatmap.png'))

## cluster proportion in each sample
cluster_table = table(Idents(filtered_obj), filtered_obj@meta.data$orig.ident)
cluster_table = cluster_table / rowSums(cluster_table)
cluster_df = as.data.frame(cluster_table)
colnames(cluster_df) = c('cluster', 'sample', 'prop')
ggplot(cluster_df, aes(x=cluster, y=prop, fill=sample)) + geom_bar(position="fill", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'pFC_sample_vs_cluster_barplot.png'))

## pFC cell type proportion in each sample
mousepFC_celltype_table = table(filtered_obj@meta.data$mousepFC_pred_celltype, filtered_obj@meta.data$orig.ident)
mousepFC_celltype_table = mousepFC_celltype_table / rowSums(mousepFC_celltype_table)
mousepFC_celltype_df = as.data.frame(mousepFC_celltype_table)
colnames(mousepFC_celltype_df) = c('celltype', 'sample', 'prop')
ggplot(mousepFC_celltype_df, aes(x=celltype, y=prop, fill=sample)) + geom_bar(position="dodge", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'pFC_sample_vs_celltype_barplot_dodge.png'))

## Saunders cell type proportion in each sample
saunders_celltype_table = table(filtered_obj@meta.data$Saunders_pred_celltype, filtered_obj@meta.data$orig.ident)
saunders_celltype_table = saunders_celltype_table / rowSums(saunders_celltype_table)
saunders_celltype_df = as.data.frame(saunders_celltype_table)
colnames(saunders_celltype_df) = c('celltype', 'sample', 'prop')
ggplot(saunders_celltype_df, aes(x=celltype, y=prop, fill=sample)) + geom_bar(position="dodge", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'pFC_sample_vs_Saunders_celltype_barplot_dodge.png'))

## use marker genes to annotate
FeaturePlot(filtered_obj, features = c('Snap25'))
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_neuronal_markergenes.png'))
FeaturePlot(filtered_obj, features = c('Gja1', 'Mlc1', 'Gib6', 'Aqp4', 'Acsbg1'))
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_astrocytes_markergenes.png'), width=10, height=10)
FeaturePlot(filtered_obj, features=c('Slc17a7', 'Slc17a6', 'Neurod6', 'Nrn1',
        'Tshz2', 'Foxp2', 'Pou3f1', 'Adam33', 'Nnat', 'Penk'), ncol=3)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_excitatory_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Gad1', 'Gad2', 'Pnoc', 'Slc32a', 'Pbx3'), ncol=2)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_inhibitory_markergenes.png'))
FeaturePlot(filtered_obj, features=c('C1qa', 'C1qb', 'C1qc', 'Ctss'), ncol=2)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_microglia_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Bmp4', 'Park4', 'Enpp6', 'Brca1'), ncol=2)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_NFOligo_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Mog', 'Aspa', 'Ermn', 'Opalin', 'Mobp'), ncol=3)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_Oligo_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Pdgfra', 'Cacng4', 'Matn4'), ncol=2)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_OPCs_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Anpep', 'Mcam', 'Cspg4', 'Pdgfrb', 'Acta2', 'Des', 'Aldha1'), ncol=3)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_mural_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Cux2', 'Calb1', 'Adm33'), ncol=1)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_L23_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Etv1', 'Pou3f1', 'Tshz2'), ncol=2)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_L5_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Syt6', 'Sla', 'Foxp2', 'Oprk1'), ncol=2)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_L6_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Sst', 'Pavlb', 'Stim2', 'Cul4a', 'B3gat2', 'Cartpt', 'Kcnmb4','Crhbp'))
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_inhibMGE_markergenes.png'))
FeaturePlot(filtered_obj, features=c('Vip', 'Cplx3', 'Synpr'))
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.5_inhibCGE_markergenes.png'))


## DE genes across conditions using mouse pFC as ref
for (celltype in unique(filtered_obj@meta.data$mousepFC_pred_celltype)) {
    cell_selection = subset(filtered_obj, cells=colnames(filtered_obj)[filtered_obj@meta.data[, 'mousepFC_pred_celltype'] == celltype])
    condition = ifelse(grepl('WT', cell_selection@meta.data$orig.ident), 'control', 'FXTAS')
    cell_selection@meta.data$condition = condition
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, logfc.threshold = 0, test.use = "wilcox",
                                        min.pct = 0, only.pos = F)
    write.csv(DEG_cell_selection, file.path(result_dir, paste0(gsub('/', '', celltype), '_DEGs.csv')), quote=F)
    DEG_df = DEG_cell_selection[DEG_cell_selection$p_val_adj < 0.1, ]
    cat(celltype, dim(DEG_df))
    VlnPlot(cell_selection, features = c('Fmr1'), group.by = "condition",
                assay = "RNA", pt.size = 0.1)
    ggsave(file.path(result_dir, paste0(celltype, '_FMR1_GE_diff.png')))

    VlnPlot(cell_selection, features = c('Prkcg'), group.by = "condition",
                assay = "RNA", pt.size = 0.1)
    ggsave(file.path(result_dir, paste0(celltype, '_Prkcg_GE_diff.png')))

    if ('Prkcg' %in% rownames(DEG_cell_selection)) { print('Prkcg in the DEGs')}
}

## DE genes across conditions using Saunders as ref
for (celltype in unique(filtered_obj@meta.data$Saunders_pred_celltype)) {
    cell_selection = subset(filtered_obj, cells=colnames(filtered_obj)[filtered_obj@meta.data[, 'Saunders_pred_celltype'] == celltype])
    condition = ifelse(grepl('WT', cell_selection@meta.data$orig.ident), 'control', 'FXTAS')
    cell_selection@meta.data$condition = condition
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, logfc.threshold = 0, test.use = "wilcox",
                                        min.pct = 0, only.pos = F)
    write.csv(DEG_cell_selection, file.path(result_dir, paste0(gsub('/', '', celltype), '_Sauders_DEGs.csv')), quote=F)
    DEG_df = DEG_cell_selection[DEG_cell_selection$p_val_adj < 0.1, ]
    cat(celltype, dim(DEG_df))
    VlnPlot(cell_selection, features = c('Fmr1'), group.by = "condition",
                assay = "RNA", pt.size = 0.1)
    ggsave(file.path(result_dir, paste0(gsub('/', '', celltype), '_FMR1_Saunders_GE_diff.png')))

    VlnPlot(cell_selection, features = c('Prkcg'), group.by = "condition",
                assay = "RNA", pt.size = 0.1)
    ggsave(file.path(result_dir, paste0(gsub('/', '', celltype), '_Prkcg_Saunders_GE_diff.png')))
}





