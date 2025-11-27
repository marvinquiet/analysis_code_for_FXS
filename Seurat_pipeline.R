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

combined_obj = merge(sample_obj_list[[1]], 
                     y=c(sample_obj_list[[2]], sample_obj_list[[3]], sample_obj_list[[4]]), 
                     add.cell.ids=names(sample_obj_list))
combined_obj[["percent.mt"]] = PercentageFeatureSet(combined_obj, pattern = "^MT-")  ## already filtered
VlnPlot(combined_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'QuanlityViolinplot.png'))

combined_obj = subset(combined_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
## subset with ArchR scATAC-seq data
GS_df = read.csv(file.path(project_dir, 'Cortex_multiome_ArchR', 'FMR_Multiome_GS.csv'), header=T, row.names=1, sep=',', check.names=F)
scATAC_id = gsub('#', '_', colnames(GS_df))
common_cells = intersect(colnames(combined_obj), scATAC_id)
combined_obj = combined_obj[, common_cells]

combined_obj = NormalizeData(combined_obj)
combined_obj = FindVariableFeatures(combined_obj, selection.method = "vst", nfeatures = 2000)
combined_obj = ScaleData(combined_obj, features=rownames(combined_obj))
combined_obj = RunPCA(combined_obj, features = VariableFeatures(object=combined_obj))

g = ElbowPlot(combined_obj, ndims=50)
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'PC_elbowplot.png')) ## will select 1:30 PCs

combined_obj = FindNeighbors(combined_obj, dims = 1:50)
combined_obj = FindClusters(combined_obj, resolution = 0.5)  ## go with 20 clusters for now
combined_obj = RunUMAP(combined_obj, dims=1:50)
DimPlot(combined_obj, reduction = "umap")
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5.png')) ## will select 1:30 PCs

## plot identity
DimPlot(combined_obj, reduction = "umap", group.by='orig.ident')
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_sample.png'))

## table of summary statistics
table(combined_obj@meta.data$orig.ident)
table(combined_obj@meta.data$orig.ident, Idents(combined_obj))

## any doublets?
require(DoubletFinder) ## https://github.com/chris-mcginnis-ucsf/DoubletFinder
homotypic.prop = modelHomotypic(combined_obj@meta.data$seurat_clusters)
nExp_poi = round(0.05*nrow(combined_obj@meta.data))  ## Assuming 5% doublet formation
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
combined_obj = doubletFinder_v3(combined_obj, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN=F, sct = FALSE)
DimPlot(combined_obj, reduction='umap', group.by='DF.classifications_0.25_0.09_804')
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_doublet.png'))

saveRDS(combined_obj, file=file.path(project_dir, 'Cortex_multiome_Seurat', 'Seurat_obj.RDS'))
write.csv(combined_obj[['RNA']]@counts, file.path(project_dir, 'Cortex_multiome_Seurat', 'GE.csv'), quote=F)
#writeMM(combined_obj[['RNA']]@counts, file.path(project_dir, 'Cortex_multiome_Seurat', 'GE.mtx'))
#write(rownames(combined_obj[['RNA']]@counts), file.path(project_dir, 'Cortex_multiome_Seurat', 'genes.tsv'))
#write(colnames(combined_obj[['RNA']]@counts), file.path(project_dir, 'Cortex_multiome_Seurat', 'barcodes.tsv'))

## ==== cell type prediction

pred_celltypes = read.csv(file.path(project_dir, 'Cortex_multiome_Seurat', 'Saunders_majorcelltypes_predcelltypes.csv'), header=T, row.names=1)
rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
combined_obj@meta.data$pred_celltype = pred_celltypes[colnames(combined_obj), 'pred_celltype']
DimPlot(combined_obj, reduction = "umap", group.by='pred_celltype')
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_Cellcano_Saunders_majorcelltype.png'))

combined_obj@meta.data$first_pred_celltype = pred_celltypes[colnames(combined_obj), 'firstround_pred_celltype']
DimPlot(combined_obj, reduction = "umap", group.by='first_pred_celltype')
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_MLP_Saunders_majorcelltype.png'))

pred_celltypes = read.csv(file.path(project_dir, 'Cortex_multiome_Seurat', 'pFC_majorcelltypes_predcelltypes.csv'), header=T, row.names=1)
rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
combined_obj@meta.data$pred_celltype = pred_celltypes[colnames(combined_obj), 'pred_celltype']
DimPlot(combined_obj, reduction = "umap", group.by='pred_celltype')
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_Cellcano_pFC_majorcelltype.png'))

combined_obj@meta.data$first_pred_celltype = pred_celltypes[colnames(combined_obj), 'firstround_pred_celltype']
DimPlot(combined_obj, reduction = "umap", group.by='first_pred_celltype')
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_MLP_pFC_majorcelltype.png'))


## ==== marker genes information from Dr. Peng Jin's paper
markerGenes = c('Snap25', 'Syt1', ## Neurons
                'Slc17a7', 'Satb2', 'Rorb', 'Neurod2', ## excitatory neurons
                'Gad1', 'Gad2', 'Nxph1', ## inhibitory neurons
                'Gfap', 'Aqp4', 'Slc1a2', ## Astrocytes
                'Csf1r', 'Cd74', 'P2ry12', 'Ptprc', 'Tmem119', ## microglia
                'Mobp', 'Mbp', 'St18', 'Klk6', 'Slc5a11', ## Oligodendrocytes
                'Pdgfra', 'Cspg4', ## Oligo precursor
                'Flt1', 'Cldn5', 'Abcb1', 'Ebf1' ## endothelial
)
FeaturePlot(combined_obj, features = markerGenes[1:9])
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_markergenes_neurons.png'))
FeaturePlot(combined_obj, features = markerGenes[10:12])
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_markergenes_astro.png'))
FeaturePlot(combined_obj, features = markerGenes[13:17])
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_markergenes_microglia.png'))
FeaturePlot(combined_obj, features = markerGenes[18:24])
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_markergenes_oligos.png'))
FeaturePlot(combined_obj, features = markerGenes[25:28])
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_markergenes_endothelial.png'))

## ==== marker genes information from Cell paper 
mousebrain_markergenes = c('Cplx3', 'Synpr', 'Sst', 'Pvalb', 'Syt6', 
                           'Fezf2', 'Nr4a2', 'Nptxr', 'Parm1', 'Gja1', 
                           'Tfr', 'Tnr', 'C1qb', 'Flt1', 'Rgs5', 'Acta2',
                           'Dcn')
FeaturePlot(combined_obj, features = mousebrain_markergenes, ncol=6)
ggsave(file.path(project_dir, 'Cortex_multiome_Seurat', 'UMAP_louvain0.5_mousebrain_markergenes.png'), width=20, height=10)


