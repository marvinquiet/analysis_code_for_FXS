# conda activate Seurat
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
set.seed(2022)

project_dir = "/beegfs/home/wma36/collaborations/Yulin/FMRpolyG_Cortex_10Xmultiomics"
result_dir = file.path(project_dir, 'Cortex_GFPonly')
dir.create(result_dir, showWarnings = FALSE)
## get inputFiles
inputFiles = c("cortex_WT1"="Cortex_GFPonly_WT1_counts",
               "cortex_WT2"="Cortex_GFPonly_WT2_counts",
               "cortex_FXTAS1"="Cortex_GFPonly_FXTAS1_counts",
               "cortex_FXTAS2"="Cortex_GFPonly_FXTAS2_counts")

#sample_obj_list = list()
#for (sample in names(inputFiles)) {
#    input_dir = file.path(project_dir, inputFiles[sample], 'outs', 'filtered_feature_bc_matrix')
#    input_mat = readMM(file.path(input_dir, 'matrix.mtx.gz'))
#    input_barcodes = scan(file.path(input_dir, 'barcodes.tsv.gz'), what=character())
#    input_features = read.csv(file.path(input_dir, 'features.tsv.gz'), header=F, sep='\t')
#
#    ## filter gene count matrix and create Seurat matrix
#    gene_expr_id = which(input_features$V3 == 'Gene Expression')
#    gene_expr_mat = input_mat[gene_expr_id, ]
#    rownames(gene_expr_mat) = input_features[gene_expr_id, 'V2']
#    colnames(gene_expr_mat) = input_barcodes
#
#    ## create Seurat object
#    sample_obj = CreateSeuratObject(counts=gene_expr_mat, project=sample, 
#                                    min.cells=0, min.features=0)
#    sample_obj_list[[sample]] = sample_obj
#}
#sample_combined_obj = merge(sample_obj_list[[1]], 
#                     y=c(sample_obj_list[[2]], sample_obj_list[[3]], sample_obj_list[[4]]), 
#                     add.cell.ids=names(sample_obj_list))
#
#saveRDS(sample_combined_obj, file.path(result_dir, 'combined_obj.RDS'))
sample_combined_obj = readRDS(file.path(result_dir, 'combined_obj.RDS'))

# --- load the annotated object from original results
annotated_obj = readRDS(file.path(project_dir, 'Cortex_multiome_Seurat_pFC_separatecutoff_integration', 'annotated_filtered_obj.RDS')) # 13971 cells
common_cells = intersect(colnames(sample_combined_obj), colnames(annotated_obj)) # 13968 cells
sample_combined_obj = sample_combined_obj[, common_cells]
annotated_obj = annotated_obj[, common_cells]
sample_combined_obj$annotated_celltype = annotated_obj$annotated_celltype
sample_combined_obj = NormalizeData(sample_combined_obj)
all.genes = rownames(sample_combined_obj)
sample_combined_obj = ScaleData(sample_combined_obj, features = all.genes)
VlnPlot(sample_combined_obj, features = c("GFP"))
ggsave(file.path(result_dir, 'GFP_per_sample.png'))

# --- set identities
VlnPlot(sample_combined_obj, features = c("GFP"), group.by='annotated_celltype', split.by='orig.ident')
ggsave(file.path(result_dir, 'GFP_per_sample_celltype.png'))

VlnPlot(sample_combined_obj, features = c("GFP"), slot="counts", log=TRUE, group.by='annotated_celltype', split.by='orig.ident')
ggsave(file.path(result_dir, 'GFP_per_sample_celltype_log.png'))

# --- compute log2FC
sample_combined_obj@meta.data$condition = ifelse(Idents(sample_combined_obj) %in% names(inputFiles)[1:2], 'control', 'disease')
sample_combined_obj$celltype_condition = paste0(sample_combined_obj$condition, '-', sample_combined_obj$annotated_celltype)
Idents(sample_combined_obj) = "celltype_condition"
celltype_list = list()
for (celltype in unique(sample_combined_obj$annotated_celltype)) {
    celltype_marker = FindMarkers(sample_combined_obj, ident.1=paste0('control-', celltype), ident.2=paste0('disease-', celltype),
				 logfc.threshold=0, min.pct=0, verbose = FALSE)
    celltype_list[[celltype]] = celltype_marker['GFP', ]
}
saveRDS(celltype_list, file.path(result_dir, "GFP_celltype_log2FC.RDS"))

# --- compute sparsity and average expression for each sample per cell type
sample_combined_obj$GFP_exprs = FetchData(object = sample_combined_obj, vars = "GFP", layer = "counts")
sample_combined_obj$GFP_expressing = sample_combined_obj$GFP_exprs > 0
write.csv(sample_combined_obj@meta.data, file.path(result_dir, 'GFP_metadata.csv'))

frac_df = sample_combined_obj@meta.data %>% group_by(annotated_celltype,orig.ident) %>% summarize(frac = mean(GFP_expressing), n=n()) %>% as.data.frame()
write.csv(frac_df, file.path(result_dir, 'GFP_faction.csv'))
summary_df = sample_combined_obj@meta.data %>% group_by(annotated_celltype, orig.ident) %>% summarize(sum_expr=sum(GFP_exprs), mean_expr=mean(GFP_exprs)) %>% as.data.frame()
write.csv(summary_df, file.path(result_dir, 'GFP_avgexprs.csv'))

nonzero_summary_df = sample_combined_obj@meta.data %>%
	    group_by(annotated_celltype, orig.ident) %>%
	    summarize(mean_nonzero = mean(GFP_exprs[GFP_exprs > 0], na.rm=TRUE)) %>%
	    as.data.frame()
write.csv(nonzero_summary_df, file.path(result_dir, 'GFP_avgexprs_nonzero.csv'))

