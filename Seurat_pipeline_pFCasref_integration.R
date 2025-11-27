## conda activate Seurat
library(Seurat)
library(Matrix)
library(ggplot2)
library(harmony)
set.seed(2022)

## filter criteria
lower_bound = 200
upper_bound = 4000 ## 10000 too loose, has a lot of noises; 3000 is too stringent

#project_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"
project_dir = "/beegfs/home/wma36/collaborations/Yulin/FMRpolyG_Cortex_10Xmultiomics"
result_dir = file.path(project_dir, paste0('Cortex_multiome_Seurat_pFC_', lower_bound, '_', upper_bound, '_integration'))
dir.create(result_dir, showWarnings = FALSE)

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

combined_obj = subset(sample_combined_obj, subset = nFeature_RNA > lower_bound & nFeature_RNA < upper_bound)

## subset with ArchR scATAC-seq data
GS_df = read.csv(file.path(project_dir, 'Cortex_multiome_ArchR', 'FMR_Multiome_GS.csv'), header=T, row.names=1, sep=',', check.names=F)
scATAC_id = gsub('#', '_', colnames(GS_df))
common_cells = intersect(colnames(combined_obj), scATAC_id)
combined_obj = combined_obj[, common_cells]
print(combined_obj)

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

## boxplot of number of genes for each cluster and predicted cell types
#pred_celltypes = read.csv(file.path(project_dir, 'Cortex_multiome_Seurat_lessstringent', 'mousepFC_lessstringent_predcelltypes.csv'), header=T, row.names=1)
#rownames(pred_celltypes) = gsub('\\.', '-', rownames(pred_celltypes))
#combined_obj@meta.data$mousepFC_pred_celltype = pred_celltypes[colnames(combined_obj), 'pred_celltype']
#DimPlot(combined_obj, reduction = "umap", group.by='mousepFC_pred_celltype')
#ggsave(file.path(result_dir, 'Harmony_integrated_UMAP_Cellcano_pFC_majorcelltype.png'))

## use marker genes to annotate
marker_list = list(
                   neuronal=c('Snap25', 'Stmn2'),
                   astrocytes=c('Gja1', 'Mlc1', 'Gib6', 'Aqp4', 'Acsbg1'),
                   endothelial=c('Flt1', 'Cldn5', 'Tm4sf1', 'Pglyrp1', 'Rgs5',
                                 'Slc38a5', 'Mgp', 'Cp'),
                   excitatory=c('Slc17a7', 'Slc17a6', 'Neurod6', 'Nrn1',
                                'Tshz2', 'Foxp2', 'Pou3f1', 'Adam33', 'Nnat', 'Penk',
                                'Slc17a8', 'Grin2b', 'Gls', 'Glul', 'Dlg4',
                                'Tbr1', 'Satb2', 'Ctp2', 'Brn2', 'Cux1', 'Cux2', 'Bcl11a', 'Bcl11b', 'Sox5'),
                   inhibitory=c('Gad1', 'Gad2', 'Pnoc', 'Slc32a', 'Pbx3', 'Adrb2', 
                                'Htr3a', 'Lhx6', 'Pvalb', 'Sst', 'Vip', 'Lamp5',
                                'Gabbr1', 'Gabbr2', 'Slc32a1', 'Gphn',
                                'Dlx1', 'Dlx2', 'Nkx2-1', 'Sox2', 'Calb2'),
                   microglia=c('Csf1r', 'C1qa', 'C1qb', 'C1qc', 'Ctss'),
                   NFoligo=c('Bmp4', 'Park4', 'Enpp6', 'Brca1'),
                   oligo=c('Mog', 'Aspa', 'Ermn', 'Opalin', 'Mobp'),
                   OPC=c('Pdgfra', 'Cacng4', 'Matn4'),
                   mural=c('Anpep', 'Mcam', 'Cspg4', 'Pdgfrb', 'Acta2', 'Des', 'Aldha1'),
                   L2_3=c('Cux2', 'Calb1', 'Adm33'),
                   L5=c('Etv1', 'Pou3f1', 'Tshz2'),
                   L6=c('Syt6', 'Sla', 'Foxp2', 'Oprk1'),
                   inhibitory_MGE=c('Sst', 'Pavlb', 'Stim2', 'Cul4a', 'B3gat2', 'Cartpt', 'Kcnmb4','Crhbp'),
                   inhibitory_CGE=c('Vip', 'Cplx3', 'Synpr'))
marker_df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(marker_df) = c('celltype', 'marker')
for (celltype in names(marker_list)) {
    markers = marker_list[[celltype]]
    for (marker in markers) {
        marker_df[nrow(marker_df)+1,] = c(celltype, marker)
    }
}

## Featureplot per marker
Idents(combined_obj) = combined_obj@meta.data$RNA_snn_res.0.5
for (i in 1:nrow(marker_df)) {
    celltype = marker_df[i, 'celltype']
    marker = marker_df[i, 'marker']
    celltype_dir = file.path(result_dir, 'markers', celltype)
    dir.create(celltype_dir, showWarnings = FALSE)
    if (!marker %in% rownames(combined_obj)) {
        cat(marker, 'not found!\n')
        next()
    }

    FeaturePlot(combined_obj, features=marker, label=T) + 
        ggtitle(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '.png')),
           width=5, height=5)
    FeaturePlot(combined_obj, features=marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '_split.png')),
           width=10, height=5)
}

## dot plot
features = marker_df$marker
features = features[!duplicated(features)]
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) {
    Idents(combined_obj) = combined_obj@meta.data[, paste0('RNA_snn_res.', resolution)]
    DotPlot(combined_obj, features=features) & coord_flip() & theme_bw()
    ggsave(file.path(result_dir, paste0('Harmony_integrated_markers_louvain', resolution, '_dotplot.png')), height=20, width=10) 
}


region_markers = c('Omp', 'Syt1', 'Tubb3', 'Snap25', 
                   'Calb2', 'Apold1', 'Sox11', 'Th', 'Dcx', 
                   'Calb1', 'Rprm', 'Lgr4', 'Stxbp6', 'Crhr1', 'Tpbg', 
                   'S100b', 'Kcnj4', 'Cdhr1' , 'Coro6', 'Eomes', 
                   'Spp1', 'Slc32a1', 'Tbx21')
region_marker_dir = file.path(result_dir, 'region_markers')
dir.create(region_marker_dir, showWarnings = FALSE)
for (region_marker in region_markers) {
    if (!region_marker %in% rownames(combined_obj)) {
        cat(region_marker, 'not found!\n')
        next()
    }
    FeaturePlot(combined_obj, features=region_marker, label=T) + 
        ggtitle(region_marker)
    ggsave(file.path(region_marker_dir, 
                     paste0('Harmony_integrated_', region_marker, '.png')),
           width=5, height=5)
    FeaturePlot(combined_obj, features=region_marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(region_marker)
    ggsave(file.path(region_marker_dir, 
                     paste0('Harmony_integrated_', region_marker, '_split.png')),
           width=10, height=5)
}


## dot plot
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) {
    Idents(combined_obj) = combined_obj@meta.data[, paste0('RNA_snn_res.', resolution)]
    DotPlot(combined_obj, features=region_markers) & coord_flip() & theme_bw()
    ggsave(file.path(region_marker_dir, paste0('Harmony_integrated_markers_louvain', resolution, '_dotplot.png')), height=20, width=10) 
}

#' another marker list
mouse_marker_info = "Neurons: Tubb3
Excitatory: Nrgn, Lpl
Immature interneurons: Htr3a
Ganglionic eminence inhibitory precursors: Ascl1, Dlx2, Top2a
Cajal Retzius cells: Reln, Calb2, Lhx5
Subventricular zone migrating neurons: Eomes, Sema3c, Neurod1
Astrocytes_1: Aqp4, Aldoc, Apoe 
Astrocytes_2: Ptx3, Igfbp5, C4b
T Cells: Cd3d, Ccr7, Cd28
Border macrophages: Lyz2, Msr1, Fcgr4
Mature oligodendrocytes: Mbp, Sirt2, Plp1
OPC: Pdgfra, Olig2, Sox10
Neural stem cells: Btg2, Dll1, Dbx2
Endothelial: Cldn5
Pericytes: Vtn"

mouse_marker_list = list()
items = unlist(strsplit(mouse_marker_info, split='\n'))
for (item in items) {
    res = unlist(strsplit(item, split=':'))
    genes = unlist(lapply(strsplit(res[2], split=','), trimws))
    mouse_marker_list[[res[1]]] = genes
}
marker_df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(marker_df) = c('celltype', 'marker')
for (celltype in names(mouse_marker_list)) {
    markers = mouse_marker_list[[celltype]]
    for (marker in markers) {
        marker_df[nrow(marker_df)+1,] = c(celltype, marker)
    }
}

## Featureplot per marker
Idents(combined_obj) = combined_obj@meta.data$RNA_snn_res.0.6
for (i in 1:nrow(marker_df)) {
    celltype = marker_df[i, 'celltype']
    marker = marker_df[i, 'marker']
    celltype_dir = file.path(result_dir, 'mouse_markers', celltype)
    dir.create(celltype_dir, showWarnings = FALSE)
    if (!marker %in% rownames(combined_obj)) {
        cat(marker, 'not found!\n')
        next()
    }

    FeaturePlot(combined_obj, features=marker, label=T) + 
        ggtitle(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '.png')),
           width=5, height=5)
    FeaturePlot(combined_obj, features=marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '_split.png')),
           width=10, height=5)
}




#' cell counts in resolution 0.6
write.csv(as.data.frame(table(combined_obj$condition, combined_obj$RNA_snn_res.0.6)),
          file.path(result_dir, 'cell_counts_resolution0.6.csv'), quote=F)

#' Generate cell markers for resolution 0.6
Idents(combined_obj) = combined_obj@meta.data$RNA_snn_res.0.6
Seurat_marker_dir = file.path(result_dir, 'Seurat_derived_markers')
dir.create(Seurat_marker_dir, showWarnings = FALSE)
for (cluster in levels(Idents(combined_obj))) {
    Seurat_markers = FindMarkers(combined_obj, ident.1=cluster, min.pct=0.25)
    Seurat_markers = Seurat_markers[Seurat_markers$p_val_adj < 0.05,]
    write.csv(Seurat_markers, file.path(Seurat_marker_dir, paste0(cluster, '_cluster_markers_adjp_0.05.csv')), quote=F)
}


#' Annotate cell clusters
combined_obj@meta.data$annotated_celltype = NA
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(1, 5, 7, 9, 14)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Excitatory'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(3, 4) 
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Inhibitory'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(6, 15)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Astrocytes'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(8)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Microglia'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(2)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Oligodendrocytes'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(11)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'OPC'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.6 %in% c(10)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Endothelial'
cell_idx = is.na(combined_obj@meta.data$annotated_celltype)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Unknown'

## filter out unknown cells
unknown_cells = combined_obj@meta.data$annotated_celltype == 'Unknown'
annotated_combined_obj = combined_obj[, !unknown_cells]
DimPlot(annotated_combined_obj, reduction = "umap", group.by='condition', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP_condition.png'),
    width=6, height=5)
DimPlot(annotated_combined_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP.png'),
    width=6, height=5)
DimPlot(annotated_combined_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)
#saveRDS(annotated_combined_obj, file=file.path(result_dir, 'Seurat_integration_annotated_object.RDS'))

annotated_combined_obj = readRDS(file.path(result_dir, 'Seurat_integration_annotated_object.RDS'))
# re-run clustering
annotated_combined_obj = Seurat::NormalizeData(annotated_combined_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE, features=rownames(annotated_combined_obj)) %>% 
    RunPCA(pc.genes = annotated_combined_obj@var.genes, npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)
DimPlot(annotated_combined_obj, reduction = "umap", group.by='condition', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'reUMAP_annotated_Harmony_integrated_UMAP_condition.png'),
    width=6, height=5)

DimPlot(annotated_combined_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'reUMAP_annotated_Harmony_integrated_UMAP.png'),
    width=6, height=5)
DimPlot(annotated_combined_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'reUMAP_annotated_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)

annotated_combined_obj@meta.data$annotated_final = paste0(annotated_combined_obj$RNA_snn_res.0.6, ':', annotated_combined_obj$annotated_celltype)
annotated_combined_obj@meta.data$annotated_final = factor(annotated_combined_obj@meta.data$annotated_final,
                                                          levels=c('1:Excitatory', '5:Excitatory', '7:Excitatory', '9:Excitatory', '14:Excitatory',
                                                                   '3:Inhibitory', '4:Inhibitory',
                                                                   '2:Oligodendrocytes', '6:Astrocytes', '15:Astrocytes', 
                                                                   '8:Microglia', '11:OPC', '10:Endothelial'))
DimPlot(annotated_combined_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, label=F, repel=T, label.size=4)
ggsave(file.path(result_dir, 'reUMAP_annotated_final_Harmony_integrated_UMAP_unlabeled.png'),
    width=6, height=5)
DimPlot(annotated_combined_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, split.by = 'condition', label=F, repel=T, label.size=4)
ggsave(file.path(result_dir, 'reUMAP_annotated_final_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)

## Featureplot per marker
new_marker_info = "Excitatory: Slc17a7, Satb2, Neurod6, Nrn1
Inhibitory: Gad1, Gad2, Dlx1
Oligo: Mobp, Mog
Astrocytes: Acsbg1, Aqp4
Microglia: Csf1r, Ctss
Endothelia: Flt1
OPC: Pdgfra
"

new_marker_list = list()
items = unlist(strsplit(new_marker_info, split='\n'))
for (item in items) {
    res = unlist(strsplit(item, split=':'))
    genes = unlist(lapply(strsplit(res[2], split=','), trimws))
    new_marker_list[[res[1]]] = genes
}
marker_df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(marker_df) = c('celltype', 'marker')
for (celltype in names(new_marker_list)) {
    markers = new_marker_list[[celltype]]
    for (marker in markers) {
        marker_df[nrow(marker_df)+1,] = c(celltype, marker)
    }
}


Idents(annotated_combined_obj) = annotated_combined_obj$annotated_final
for (i in 1:nrow(marker_df)) {
    celltype = marker_df[i, 'celltype']
    marker = marker_df[i, 'marker']
    celltype_dir = file.path(result_dir, 'mouse_markers_new', celltype)
    dir.create(celltype_dir, showWarnings = FALSE)
    if (!marker %in% rownames(annotated_combined_obj)) {
        cat(marker, 'not found!\n')
        next()
    }
    FeaturePlot(annotated_combined_obj, features=marker, label=T) + 
        ggtitle(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '.png')),
           width=5, height=5)
    FeaturePlot(annotated_combined_obj, features=marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(paste(marker, '(', celltype, ')'))
    ggsave(file.path(celltype_dir, 
                     paste0('Harmony_integrated_', marker, '_split.png')),
           width=10, height=5)
}

features = marker_df$marker
features = features[!duplicated(features)]
annotated_combined_obj$annotated_celltype = factor(annotated_combined_obj$annotated_celltype,
                                                  levels=c('Excitatory', 'Inhibitory', 'Oligodendrocytes',
                                                           'Astrocytes', 'Microglia', 'Endothelial', 'OPC'))
Idents(annotated_combined_obj) = annotated_combined_obj$annotated_celltype
g = DotPlot(annotated_combined_obj, features=features) & coord_flip() & theme_bw()
g = g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file.path(result_dir, 'mouse_markers_new', 'celltype_dotplot.pdf'), height=6, width=8)

# --- generate heatmap
write.csv(annotated_combined_obj@meta.data, file.path(result_dir, 'Figure2D_metadata.csv'))
write.csv(annotated_combined_obj@assays$RNA@scale.data[features, ], file.path(result_dir, 'Figure2D_heatmap.csv'))

g = DoHeatmap(annotated_combined_obj, features=features, size=2, draw.lines=T, raster=F) + scale_fill_gradient2(low='blue', high='red')
ggsave(file.path(result_dir, 'mouse_markers_new', 'celltype_heatmap.pdf'), height=6, width=10)

# umap of neurons
neuron_obj = subset(annotated_combined_obj, subset=annotated_celltype %in% c('Excitatory', 'Inhibitory'))
neuron_obj = Seurat::NormalizeData(neuron_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = neuron_obj@var.genes, npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)
DimPlot(neuron_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, label=T, repel=T, label.size=2)
ggsave(file.path(result_dir, 'reUMAP_annotated_neurons_Harmony_integrated_UMAP.png'),
    width=6, height=5)
DimPlot(neuron_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'reUMAP_annotated_neurons_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)

# umap of inhibitory neurons
inh_neuron_obj = subset(annotated_combined_obj, subset=annotated_celltype %in% c('Inhibitory'))
inh_neuron_obj = Seurat::NormalizeData(inh_neuron_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = neuron_obj@var.genes, npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)
DimPlot(inh_neuron_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, label=F, repel=T, label.size=2)
ggsave(file.path(result_dir, 'reUMAP_annotated_inhibitory_neurons_Harmony_integrated_UMAP.png'),
    width=6, height=5)
DimPlot(inh_neuron_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'reUMAP_annotated_inhibitory_neurons_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)






# umap of non-neurons
nonneuron_obj = subset(annotated_combined_obj, subset=annotated_celltype %in% c('Oligodendrocytes', 'Astrocytes', 'Microglia', 'Endothelial', 'OPC'))
nonneuron_obj = Seurat::NormalizeData(nonneuron_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = neuron_obj@var.genes, npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)
DimPlot(nonneuron_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, label=T, repel=T, label.size=2)
ggsave(file.path(result_dir, 'reUMAP_annotated_nonneurons_Harmony_integrated_UMAP.png'),
    width=6, height=5)
DimPlot(nonneuron_obj, reduction = "umap", group.by='annotated_final', 
        pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'reUMAP_annotated_nonneurons_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)


## check coverage for annotated combined object
ggplot(annotated_combined_obj@meta.data, aes(x=nFeature_RNA, y=annotated_celltype, fill=annotated_celltype)) + geom_boxplot() + 
    geom_vline(xintercept=400, linetype="dashed", color='blue', size=1.5) + 
    geom_vline(xintercept=1000, linetype="dashed", color='red', size=1.5) 
ggsave(file.path(result_dir, 'nFeature_based_on_annotated_celltypes.png'))
ggplot(annotated_combined_obj@meta.data, aes(x=nFeature_RNA, y=RNA_snn_res.0.6, fill=RNA_snn_res.0.6)) + geom_boxplot() + 
    geom_vline(xintercept=400, linetype="dashed", color='blue', size=1.5) + 
    geom_vline(xintercept=1000, linetype="dashed", color='red', size=1.5) 
ggsave(file.path(result_dir, 'nFeature_based_on_clusters.png'))


## cell type composition
annotated_celltype_table = table(annotated_combined_obj@meta.data$annotated_celltype, annotated_combined_obj@meta.data$orig.ident)
write.csv(as.data.frame(annotated_celltype_table), file.path(result_dir, 'annotated_Harmony_integrated_celltype_number.csv'), quote=F)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
annotated_celltype_df = as.data.frame(annotated_celltype_table)
colnames(annotated_celltype_df) = c('celltype', 'sample', 'prop')
ggplot(annotated_celltype_df, aes(x=sample, y=prop, fill=celltype)) + geom_bar(position="fill", stat="identity", color = "black") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_celltype_proportion.png'))

# cluster composition
annotated_celltype_table = table(annotated_combined_obj@meta.data$annotated_final, annotated_combined_obj@meta.data$orig.ident)
write.csv(as.data.frame(annotated_celltype_table), file.path(result_dir, 'annotated_Harmony_integrated_cluster_number.csv'), quote=F)

## per condition
annotated_celltype_table = table(annotated_combined_obj@meta.data$annotated_celltype, annotated_combined_obj@meta.data$condition)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
annotated_celltype_df = as.data.frame(annotated_celltype_table)
colnames(annotated_celltype_df) = c('celltype', 'condition', 'prop')
ggplot(annotated_celltype_df, aes(x=condition, y=prop, fill=celltype)) + geom_bar(position="fill", stat="identity", color = "black") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_celltype_proportion_per_condition.png'))

## === DEG analysis
Idents(annotated_combined_obj) = annotated_combined_obj@meta.data$condition
#condition_DEGs = FindMarkers(annotated_combined_obj, logfc.threshold = 0, ident.1="disease", ident.2="control")
#condition_DEGs = condition_DEGs[condition_DEGs$p_val_adj < 0.05, ]
#write.csv(condition_DEGs, file.path(result_dir, 'condition_DEGs.csv'), quote=F)
#up_condition_DEGs = condition_DEGs[condition_DEGs$avg_log2FC > 0, ]
#write.csv(up_condition_DEGs, file.path(result_dir, 'condition_DEGs_UP.csv'), quote=F)
#down_condition_DEGs = condition_DEGs[condition_DEGs$avg_log2FC < 0, ]
#write.csv(down_condition_DEGs, file.path(result_dir, 'condition_DEGs_DOWN.csv'), quote=F)
# direcitons: https://github.com/satijalab/seurat/issues/5127
DEG_dir = file.path(result_dir, 'DEGs_nothreshold')
dir.create(DEG_dir, showWarnings=F)
for (cluster_num in unique(annotated_combined_obj$RNA_snn_res.0.6)) {
    celltype_combined_obj = subset(annotated_combined_obj, subset = RNA_snn_res.0.6 == cluster_num)
    celltype_DEGs = FindMarkers(celltype_combined_obj, logfc.threshold = 0, ident.1="disease", ident.2="control")
    write.csv(celltype_DEGs, file.path(DEG_dir, paste0(cluster_num, '_condition_DEGs.csv')), quote=F)
    up_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC > 0, ]
    write.csv(up_celltype_DEGs, file.path(DEG_dir, paste0(cluster_num, '_condition_DEGs_UP.csv')), quote=F)
    down_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC < 0, ]
    write.csv(down_celltype_DEGs, file.path(DEG_dir, paste0(cluster_num, '_condition_DEGs_DOWN.csv')), quote=F)
}


Idents(annotated_combined_obj) = annotated_combined_obj@meta.data$condition
#condition_DEGs = FindMarkers(annotated_combined_obj, logfc.threshold = 0, ident.1="disease", ident.2="control")
#condition_DEGs = condition_DEGs[condition_DEGs$p_val_adj < 0.05, ]
#write.csv(condition_DEGs, file.path(result_dir, 'condition_DEGs.csv'), quote=F)
#up_condition_DEGs = condition_DEGs[condition_DEGs$avg_log2FC > 0, ]
#write.csv(up_condition_DEGs, file.path(result_dir, 'condition_DEGs_UP.csv'), quote=F)
#down_condition_DEGs = condition_DEGs[condition_DEGs$avg_log2FC < 0, ]
#write.csv(down_condition_DEGs, file.path(result_dir, 'condition_DEGs_DOWN.csv'), quote=F)
# direcitons: https://github.com/satijalab/seurat/issues/5127
DEG_dir = file.path(result_dir, 'DEGs_nothreshold_per_celltype')
dir.create(DEG_dir, showWarnings=F)
for (celltype in unique(annotated_combined_obj$annotated_celltype)) {
    celltype_combined_obj = subset(annotated_combined_obj, subset = annotated_celltype == celltype)
    celltype_DEGs = FindMarkers(celltype_combined_obj, logfc.threshold = 0, ident.1="disease", ident.2="control")
    write.csv(celltype_DEGs, file.path(DEG_dir, paste0(celltype, '_condition_DEGs.csv')), quote=F)
    up_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC > 0, ]
    write.csv(up_celltype_DEGs, file.path(DEG_dir, paste0(celltype, '_condition_DEGs_UP.csv')), quote=F)
    down_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC < 0, ]
    write.csv(down_celltype_DEGs, file.path(DEG_dir, paste0(celltype, '_condition_DEGs_DOWN.csv')), quote=F)
}
#for (celltype in unique(annotated_combined_obj@meta.data$annotated_celltype)) {
#    celltype_combined_obj = subset(annotated_combined_obj, subset = annotated_celltype == celltype)
#    celltype_DEGs = FindMarkers(celltype_combined_obj, logfc.threshold = 0, ident.1="disease", ident.2="control")
#    celltype_DEGs = celltype_DEGs[celltype_DEGs$p_val_adj < 0.05, ]
#    write.csv(celltype_DEGs, file.path(result_dir, paste0(celltype, '_condition_DEGs.csv')), quote=F)
#    up_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC > 0, ]
#    write.csv(up_celltype_DEGs, file.path(result_dir, paste0(celltype, '_condition_DEGs_UP.csv')), quote=F)
#    down_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC < 0, ]
#    write.csv(down_celltype_DEGs, file.path(result_dir, paste0(celltype, '_condition_DEGs_DOWN.csv')), quote=F)
#}


#for (celltype in unique(filtered_obj@meta.data$Saunders_pred_celltype)) {
#    cell_selection = subset(filtered_obj, cells=colnames(filtered_obj)[filtered_obj@meta.data[, 'Saunders_pred_celltype'] == celltype])
#    condition = ifelse(grepl('WT', cell_selection@meta.data$orig.ident), 'control', 'FXTAS')
#    cell_selection@meta.data$condition = condition
#    cell_selection = SetIdent(cell_selection, value = "condition")
#    DEG_cell_selection = FindAllMarkers(cell_selection, logfc.threshold = 0, test.use = "wilcox",
#                                        min.pct = 0, only.pos = F)
#    write.csv(DEG_cell_selection, file.path(result_dir, paste0(gsub('/', '', celltype), '_Sauders_DEGs.csv')), quote=F)
#    DEG_df = DEG_cell_selection[DEG_cell_selection$p_val_adj < 0.1, ]
#    cat(celltype, dim(DEG_df))
#    VlnPlot(cell_selection, features = c('Fmr1'), group.by = "condition",
#                assay = "RNA", pt.size = 0.1)
#    ggsave(file.path(result_dir, paste0(gsub('/', '', celltype), '_FMR1_Saunders_GE_diff.png')))
#
#    VlnPlot(cell_selection, features = c('Prkcg'), group.by = "condition",
#                assay = "RNA", pt.size = 0.1)
#    ggsave(file.path(result_dir, paste0(gsub('/', '', celltype), '_Prkcg_Saunders_GE_diff.png')))
#}



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





