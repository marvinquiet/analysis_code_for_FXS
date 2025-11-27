## conda activate Seurat
library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
set.seed(2022)

project_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"
result_dir = file.path(project_dir, 'Cortex_multiome_Seurat_pFC_separatecutoff_integration_neuronal')
dir.create(result_dir, showWarnings = FALSE)
## === save filtered object
filtered_obj = readRDS(file.path(result_dir, 'Seurat_pFC_200_10000_integration_annotated_object.RDS'))
#filtered_obj = readRDS(file.path(result_dir, 'filtered_object.RDS')) # from separate cutoff
#filtered_obj@meta.data$annotated_celltype = NA
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(2, 3, 7, 9, 20)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Excitatory neuron'
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(0, 12, 14, 8, 4, 13) ## remove 16
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Inhibitory neuron'
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(5)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Astrocytes'
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(6)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Microglia'
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(1)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Oligo'
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(10)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'OPC'
#cell_idx = filtered_obj@meta.data$RNA_snn_res.0.3 %in% c(17)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Endothelial'
#cell_idx = is.na(filtered_obj@meta.data$annotated_celltype)
#filtered_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Unknown'
#unknown_cells = filtered_obj@meta.data$annotated_celltype == 'Unknown'
#filtered_obj = filtered_obj[, !unknown_cells] ## 13971 cells

#' select neuronal cells out and perform clustering on those marker genes
neuronal_obj = subset(filtered_obj, subset=annotated_celltype %in% c("Inhibitory neuron", "Excitatory neuron", "Neuronal cells")) ## 11,038 cells
neuronal_obj@meta.data$original_clusters = neuronal_obj@meta.data$RNA_snn_res.0.3
neuronal_obj = NormalizeData(neuronal_obj)
neuronal_obj = ScaleData(neuronal_obj)
neuronal_obj = RunPCA(neuronal_obj, npcs=5)
neuronal_obj = FindNeighbors(neuronal_obj, dims=1:5)
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) {
    neuronal_obj = FindClusters(neuronal_obj, resolution=resolution)
    neuronal_obj = RunUMAP(neuronal_obj, dims=1:5)
    DimPlot(neuronal_obj, reduction = "umap", pt.size = .1, label=T)
    ggsave(file.path(result_dir, 'clusters', paste0('sub_neuronal_UMAP_louvain', resolution, '.png')),
        width=5, height=5)
    DimPlot(neuronal_obj, reduction = "umap", pt.size = .1, split.by = 'condition', label=T)
    ggsave(file.path(result_dir, 'clusters', paste0('sub_neuronal_UMAP_louvain', resolution, '_condition_split.png')),
        width=10, height=5)
}

# check original cluster and current cluster correspondence
table(neuronal_obj@meta.data$original_clusters, neuronal_obj@meta.data$RNA_snn_res.0.8)

## check marker genes
marker_list = list(
                   neuronal=c('Stmn2'),
                   excitatory=c('Slc17a7', 'Satb2', 'Tbr1', 'Neurod6', 'Nrn1', 
                                'Sox5'),
                   inhibitory=c('Gad1', 'Gad2', 'Dlx1', 'Dlx2', 'Sox2', 'Galb2')
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
    if (!marker %in% rownames(neuronal_obj)) {
        cat(marker, 'not found!\n')
        next()
    }
    FeaturePlot(neuronal_obj, features=marker, label=T) + 
        ggtitle(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '.png')),
           width=5, height=5)
    FeaturePlot(neuronal_obj, features=marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '_split.png')),
           width=10, height=5)
}
## dot plot
features = marker_df$marker
features = features[!duplicated(features)]
for (resolution in c(0.8)) {
    Idents(neuronal_obj) = neuronal_obj@meta.data[, paste0('RNA_snn_res.', resolution)]
    DotPlot(neuronal_obj, features=features) & coord_flip() & theme_bw()
    ggsave(file.path(result_dir, paste0('Harmony_integrated_markers_louvain', resolution, '_dotplot.png')), height=20, width=10) 
}


#' Derive marker genes from pFC dataset




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


