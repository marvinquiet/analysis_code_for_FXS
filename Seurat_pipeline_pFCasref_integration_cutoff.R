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
combined_obj = subset(sample_combined_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000) ## set to 200 to 3000

result_dir = file.path(project_dir, 'Cortex_multiome_Seurat_pFC_separatecutoff_integration')
dir.create(result_dir, showWarnings = FALSE)

## subset with ArchR scATAC-seq data
GS_df = read.csv(file.path(project_dir, 'Cortex_multiome_ArchR', 'FMR_Multiome_GS.csv'), header=T, row.names=1, sep=',', check.names=F)
scATAC_id = gsub('#', '_', colnames(GS_df))
common_cells = intersect(colnames(combined_obj), scATAC_id)
combined_obj = combined_obj[, common_cells]
print(combined_obj)

library(harmony)
combined_obj = Seurat::NormalizeData(combined_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = combined_obj@var.genes, npcs = 30, verbose = FALSE)
combined_obj@meta.data$condition = ifelse(Idents(combined_obj) %in% names(inputFiles)[1:2], 'control', 'disease')
combined_obj = combined_obj %>% 
        RunHarmony("condition", plot_convergence = TRUE)
DimPlot(object = combined_obj, reduction = "harmony", pt.size = .1, group.by = "condition")
ggsave(file.path(result_dir, 'Harmony_integrated_UMAP_PCA_condition.png'),
    width=5, height=5)

combined_obj = combined_obj %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20)
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) {
    combined_obj = combined_obj %>% 
        FindClusters(resolution = resolution) %>% 
        identity()
    DimPlot(combined_obj, reduction = "umap", pt.size = .1, label=T)
    ggsave(file.path(result_dir, paste0('Harmony_integrated_UMAP_louvain', resolution, '.png')),
        width=5, height=5)
    DimPlot(combined_obj, reduction = "umap", pt.size = .1, split.by = 'condition', label=T)
    ggsave(file.path(result_dir, paste0('Harmony_integrated_UMAP_louvain', resolution, '_condition_split.png')),
        width=10, height=5)
}

## set cell types based on resolution 0.3 and endothelials by 0.5
combined_obj@meta.data$annotated_celltype = NA
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(1, 5, 8, 14)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Excitatory neuron'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(3, 4, 12) ## 17 unsure, in my understanding, cluster 13 also unsure
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Inhibitory neuron'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(0, 11)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Neuronal cells'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(6)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Astrocytes'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(7)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Microglia'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(2)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Oligo'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(10)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'OPC'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.3 %in% c(9)
#combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'NF Oligo'
#cell_idx = combined_obj@meta.data$RNA_snn_res.0.5 %in% c(20)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Endothelial'
cell_idx = is.na(combined_obj@meta.data$annotated_celltype)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Unknown'
unknown_cells = combined_obj@meta.data$annotated_celltype == 'Unknown'
annotated_combined_obj = combined_obj[, !unknown_cells]
## save first round object
saveRDS(annotated_combined_obj, file.path(result_dir, 'firstround_annotated_object.RDS'))

neuronal_cells = colnames(annotated_combined_obj)[annotated_combined_obj$annotated_celltype %in% c('Excitatory neuron', 'Inhibitory neuron', 'Neuronal cells')]
neuronal_obj = annotated_combined_obj[, neuronal_cells]
filtered_neuronal_cells = colnames(neuronal_obj)[neuronal_obj@meta.data$nFeature_RNA >= 800]

nonneuronal_cells = colnames(annotated_combined_obj)[!annotated_combined_obj$annotated_celltype %in% c('Excitatory neuron', 'Inhibitory neuron', 'Neuronal cells')]
all_cells = c(nonneuronal_cells, filtered_neuronal_cells)
filtered_obj = annotated_combined_obj[, all_cells]

## === after filtering
filtered_obj = NormalizeData(filtered_obj)
filtered_obj = FindVariableFeatures(filtered_obj, selection.method = "vst", nfeatures = 2000)
filtered_obj = ScaleData(filtered_obj, features=rownames(filtered_obj))
filtered_obj = RunPCA(filtered_obj, features = VariableFeatures(object=filtered_obj))
g = ElbowPlot(filtered_obj, ndims=50)
ggsave(file.path(result_dir, 'filtered_PC_elbowplot.png')) ## will select 1:40 PCs

filtered_obj = FindNeighbors(filtered_obj, dims = 1:40)
filtered_obj = FindClusters(filtered_obj, resolution = 0.3)
filtered_obj = RunUMAP(filtered_obj, dims=1:40)
DimPlot(filtered_obj, reduction = "umap", label=T)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.3.png'))
DimPlot(filtered_obj, reduction = "umap", group.by='orig.ident')
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.3_sample.png'))
filtered_obj@meta.data$condition = ifelse(grepl('WT', filtered_obj@meta.data$orig.ident), 'control', 'FXTAS2')
DimPlot(filtered_obj, reduction = "umap", split.by='condition', label=T)
ggsave(file.path(result_dir, 'filtered_UMAP_louvain0.3_condition.png'), width=10, height=5)

## use marker genes to annotate
marker_list = list(
                   neuronal=c('Stmn2'),
                   excitatory=c('Slc17a7', 'Satb2', 'Tbr1', 'Neurod6', 'Nrn1', 
                                'Sox5'),
                   inhibitory=c('Gad1', 'Gad2', 'Dlx1', 'Dlx2', 'Sox2', 'Galb2'),
                   microglia=c('Ctss', 'Csf1r', 'Aif1'),
                   astrocytes=c('Aqp4', 'Gia1', 'Acsbg1'),
                   oligodendrocytes=c('Mog', 'Mobp', 'Ermn', 'Aspa'),
                   OPC=c('Pdgfra', 'Olig1', 'Matn4'),
                   endothelial=c('Flt1', 'Cldn5', 'Itm2a'),
                   ependymal=c('Foxj1', 'Tuba1a'),
                   choroid_plexus=c('Ttr', 'Folr1', 'Prlr'),
                   mural=c('Pdgfrb', 'Cspg4')
                   )
marker_df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(marker_df) = c('celltype', 'marker')
for (celltype in names(marker_list)) {
    markers = marker_list[[celltype]]
    for (marker in markers) {
        marker_df[nrow(marker_df)+1,] = c(celltype, marker)
    }
}

## Featureplot per marker
for (i in 1:nrow(marker_df)) {
    celltype = marker_df[i, 'celltype']
    marker = marker_df[i, 'marker']
    celltype_dir = file.path(result_dir, 'markers', celltype)
    dir.create(celltype_dir, showWarnings = FALSE)
    if (!marker %in% rownames(filtered_obj)) {
        cat(marker, 'not found!\n')
        next()
    }

    FeaturePlot(filtered_obj, features=marker, label=T) + 
        ggtitle(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '.png')),
           width=5, height=5)
    FeaturePlot(filtered_obj, features=marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '_split.png')),
           width=10, height=5)
}


## dot plot
features = marker_df$marker
features = features[!duplicated(features)]
for (resolution in c(0.3)) {
    Idents(filtered_obj) = filtered_obj@meta.data[, paste0('RNA_snn_res.', resolution)]
    DotPlot(filtered_obj, features=features) & coord_flip() & theme_bw()
    ggsave(file.path(result_dir, paste0('Harmony_integrated_markers_louvain', resolution, '_dotplot.png')), height=20, width=10) 
}


ggplot(filtered_obj@meta.data, aes(x=nFeature_RNA, y=RNA_snn_res.0.3, fill=RNA_snn_res.0.3)) + geom_boxplot() + 
    geom_vline(xintercept=800, linetype="dashed", color='red', linewidth=1.5) 
ggsave(file.path(result_dir, 'nFeature_based_on_clusters.png'))

## === save filtered object
saveRDS(filtered_obj, file.path(result_dir, 'filtered_object.RDS'))

## === annotation
## set cell types based on resolution 0.3
filtered_obj@meta.data$annotated_celltype = NA
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(2, 3, 7, 9, 20)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Excitatory neuron'
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(0, 12, 14, 8, 4, 13) ## remove 16
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Inhibitory neuron'
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(5)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Astrocytes'
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(6)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Microglia'
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(1)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Oligo'
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(10)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'OPC'
cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(17)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Endothelial'
cell_idx = is.na(filtered_obj@meta.data$annotated_celltype)
filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Unknown'



## filter out unknown cells
unknown_cells = filtered_obj@meta.data$annotated_celltype == 'Unknown'
annotated_filtered_obj = filtered_obj[, !unknown_cells]
DimPlot(annotated_filtered_obj, reduction = "umap", group.by='condition', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'annotated_filtered_Harmony_integrated_UMAP_condition.png'),
    width=6, height=5)
DimPlot(annotated_filtered_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'annotated_filtered_Harmony_integrated_UMAP.png'),
    width=6, height=5)
DimPlot(annotated_filtered_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'annotated_filtered_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)


## check coverage for annotated combined object
ggplot(annotated_filtered_obj@meta.data, aes(x=nFeature_RNA, y=annotated_celltype, fill=annotated_celltype)) + geom_boxplot() + 
    geom_vline(xintercept=400, linetype="dashed", color='blue', size=1.5) + 
    geom_vline(xintercept=1000, linetype="dashed", color='red', size=1.5) 
ggsave(file.path(result_dir, 'nFeature_based_on_annotated_filtered_celltypes.png'))
## cell type composition
annotated_celltype_table = table(annotated_filtered_obj@meta.data$annotated_celltype, annotated_filtered_obj@meta.data$orig.ident)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
annotated_celltype_df = as.data.frame(annotated_celltype_table)
colnames(annotated_celltype_df) = c('celltype', 'sample', 'prop')
ggplot(annotated_celltype_df, aes(x=sample, y=prop, fill=celltype)) + geom_bar(position="fill", stat="identity", color = "black") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'annotated_filtered_Harmony_integrated_celltype_proportion.png'))

annotated_celltype_table = table(annotated_filtered_obj@meta.data$annotated_celltype, annotated_filtered_obj@meta.data$condition)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
annotated_celltype_df = as.data.frame(annotated_celltype_table)
colnames(annotated_celltype_df) = c('celltype', 'condition', 'prop')
ggplot(annotated_celltype_df, aes(x=condition, y=prop, fill=celltype)) + geom_bar(position="fill", stat="identity", color = "black") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'annotated_filtered_Harmony_integrated_celltype_proportion_per_condition.png'))

annotated_cluster_table = table(as.character(annotated_filtered_obj@meta.data$RNA_snn_res.0.3), annotated_filtered_obj@meta.data$condition)
annotated_cluster_table = t(t(annotated_cluster_table) / colSums(annotated_cluster_table))
annotated_cluster_df = as.data.frame(annotated_cluster_table)
colnames(annotated_cluster_df) = c('cluster', 'condition', 'prop')
ggplot(annotated_cluster_df, aes(x=condition, y=prop, fill=cluster)) + geom_bar(position="fill", stat="identity", color = "black") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'annotated_filtered_Harmony_integrated_cluster_proportion_per_condition.png'))


## ==== Cellcano pipeline
#writeMM(annotated_filtered_obj[['RNA']]@counts, file.path(result_dir, 'annotated_filtered_pFC_GE.mtx'))
#write(rownames(annotated_filtered_obj[['RNA']]@counts), file.path(result_dir, 'annotated_filtered_pFC_GE_genes.tsv'))
#write(colnames(annotated_filtered_obj[['RNA']]@counts), file.path(result_dir, 'annotated_filtered_pFC_GE_barcodes.tsv'))
## Cellcano predict -i annotated_filtered_pFC_GE --trained_model mousepFC_Cellcano_output/majorcelltypes_predMLP_model/ -o mousepFC_Cellcano_output --prefix annotated_filtered_pFC_GE
## Cellcano predict -i annotated_filtered_pFC_GE --trained_model Saunders_Cellcano_output/majorcelltypes_predMLP_model/ -o Saunders_Cellcano_output --prefix annotated_filtered_pFC_GE_Saunders
#
### load predicted cell type as validation
### plot with mouse pFC as ref
#pred_celltypes = read.csv(file.path(result_dir, 'annotated_filtered_pFC_GEcelltypes.csv'), header=T, row.names=1)
#rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
#annotated_filtered_obj@meta.data$mousepFC_pred_celltype = pred_celltypes[colnames(annotated_filtered_obj), 'pred_celltype']
#DimPlot(annotated_filtered_obj, reduction = "umap", group.by='mousepFC_pred_celltype')
#ggsave(file.path(result_dir, 'annotated_filtered_UMAP_louvain0.3_Cellcano_pFC_majorcelltype.png'))
#
### plot with Saunders as ref
#pred_celltypes = read.csv(file.path(result_dir, 'annotated_filtered_pFC_GE_Saunderscelltypes.csv'), header=T, row.names=1)
#rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
#annotated_filtered_obj@meta.data$Saunders_pred_celltype = pred_celltypes[colnames(annotated_filtered_obj), 'pred_celltype']
#DimPlot(annotated_filtered_obj, reduction = "umap", group.by='Saunders_pred_celltype')
#ggsave(file.path(result_dir, 'annotated_filtered_UMAP_louvain0.3_Cellcano_pFC_Saunders_majorcelltype.png'))
## ==== 

## === Number of cells for cluster and major cell types
annotated_celltype_table = table(annotated_filtered_obj@meta.data$annotated_celltype, annotated_filtered_obj@meta.data$condition)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
print(annotated_celltype_table)
annotated_cluster_table = table(annotated_filtered_obj@meta.data$RNA_snn_res.0.3, annotated_filtered_obj@meta.data$condition)
annotated_cluster_table = t(t(annotated_cluster_table) / colSums(annotated_cluster_table))
print(annotated_cluster_table)

## === DEGs per cluster and major cell types
marker_dir = file.path(result_dir, 'DEGs')
dir.create(marker_dir, showWarnings = FALSE)
for (celltype in unique(annotated_filtered_obj@meta.data$annotated_celltype)) {
    cell_selection = subset(annotated_filtered_obj, cells=colnames(annotated_filtered_obj)[annotated_filtered_obj@meta.data[, 'annotated_celltype'] == celltype])
    condition = ifelse(grepl('WT', cell_selection@meta.data$orig.ident), 'control', 'FXTAS')
    cell_selection@meta.data$condition = condition
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, logfc.threshold = 0, test.use = "wilcox",
                                        min.pct = 0, only.pos = F)
    write.csv(DEG_cell_selection, file.path(marker_dir, paste0(gsub('/', '', celltype), '_celltype_DEGs.csv')), quote=F)
    DEG_df = DEG_cell_selection[DEG_cell_selection$p_val_adj < 0.1, ]
    cat(celltype, dim(DEG_df))
}


for (cluster in unique(annotated_filtered_obj@meta.data$RNA_snn_res.0.3)) {
    cell_selection = subset(annotated_filtered_obj, cells=colnames(annotated_filtered_obj)[annotated_filtered_obj@meta.data[, 'RNA_snn_res.0.3'] == cluster])
    condition = ifelse(grepl('WT', cell_selection@meta.data$orig.ident), 'control', 'FXTAS')
    cell_selection@meta.data$condition = condition
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, logfc.threshold = 0, test.use = "wilcox",
                                        min.pct = 0, only.pos = F)
    write.csv(DEG_cell_selection, file.path(marker_dir, paste0(gsub('/', '', cluster), '_cluster_DEGs.csv')), quote=F)
    DEG_df = DEG_cell_selection[DEG_cell_selection$p_val_adj < 0.1, ]
    cat(cluster, dim(DEG_df))
}


