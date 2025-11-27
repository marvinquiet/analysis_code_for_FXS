## conda activate Seurat
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(Seurat)
library(Matrix)
library(ggplot2)
set.seed(2022)

#source('Manuscript_figures_utils.R')

#project_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"
project_dir = "/beegfs/home/wma36/collaborations/Yulin/FMRpolyG_Cortex_10Xmultiomics"
result_dir = file.path(project_dir, 'Cortex_multiome_Seurat_pFC_separatecutoff_integration')
dir.create(result_dir, showWarnings = FALSE)
figure_dir = file.path(project_dir, 'manuscript_figure')
dir.create(figure_dir, showWarnings = FALSE)

annotated_filtered_obj = readRDS(file.path(result_dir, 'annotated_filtered_obj.RDS'))

## --- Figure 3.1
split_plot = DimPlot(annotated_filtered_obj, reduction = "umap", pt.size = .1, group.by='RNA_snn_res.0.3', label=T)
g = DimPlot(annotated_filtered_obj, reduction = "umap", group.by='condition', pt.size = .1, label=F)
g = g + split_plot$layers[[2]]
ggsave(file.path(figure_dir, 'Figure3.1.b_condition_UMAP.pdf'),
    width=5, height=5)
g = DimPlot(annotated_filtered_obj, reduction = "umap", group.by='annotated_celltype', pt.size = .1, label=F)
g = g + split_plot$layers[[2]]
ggsave(file.path(figure_dir, 'Figure3.1.a_celltype_UMAP.pdf'),
    width=6, height=5)

## --- Figure 3.2
genes = c('Stmn2', 'Gad1', 'Gad2', 'Dlx1',
          'Slc17a7', 'Nrn1', 'Satb2', 'Tbr1',
          'Csf1r', 'Ctss', 'Aif1',
          'Acsbg1', 'Gja1', 'Aqp4',
          'Aspa', 'Mobp', 'Mog',
          'Pdgfra', 'Matn4', 'Flt1')
annotated_filtered_obj$annotated_celltype = factor(annotated_filtered_obj$annotated_celltype, 
                                                   levels=c('Inhibitory neuron',
                                                            'Excitatory neuron',
                                                            'Microglia', 'Astrocytes',
                                                            'Oligo', 'OPC', 'Endothelial'))
DoHeatmap(annotated_filtered_obj, 
          features=genes,
          group.by='annotated_celltype',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5,
          raster=F) + 
          theme(legend.position='bottom') +
          guides(color=guide_legend(ncol=7))
ggsave(file.path(figure_dir, 'Figure3.2_heatmap.pdf'), height=4, width=8)

# --- revision Figure 3.2
annotated_result_dir = file.path(project_dir, "Cortex_multiome_Seurat_pFC_200_4000_integration")
annotated_combined_obj = readRDS(file.path(annotated_result_dir, 'Seurat_integration_annotated_object.RDS'))
genes = c('Slc17a7', 'Satb2', 'Neurod6', 'Nrn1',
	  'Gad1', 'Gad2', 'Dlx1', 'Pvalb', 'Sst', 'Vip', 'Adarb2',
	  'Mobp', 'Mog',
	  'Acsbg1', 'Aqp4',
	  'Csf1r', 'Ctss',
	  'Flt1', 'Pdgfra'
	  )
annotated_combined_obj@meta.data$annotated_final = paste0(annotated_combined_obj$RNA_snn_res.0.6, ':', annotated_combined_obj$annotated_celltype)
annotated_combined_obj@meta.data$annotated_final = factor(annotated_combined_obj@meta.data$annotated_final,
                                                          levels=c('1:Excitatory', '5:Excitatory', '7:Excitatory', '9:Excitatory', '14:Excitatory',
                                                                   '3:Inhibitory', '4:Inhibitory',
                                                                   '2:Oligodendrocytes', '6:Astrocytes', '15:Astrocytes', 
                                                                   '8:Microglia', '11:OPC', '10:Endothelial'))
Idents(annotated_combined_obj) = annotated_combined_obj$annotated_final
DoHeatmap(annotated_combined_obj, 
          features=genes,
	  size=2, draw.lines=T, raster=F) + 
scale_fill_gradient2(low='blue', high='red')
ggsave(file.path(figure_dir, 'Figure3.2_heatmap_subcelltype.pdf'), height=6, width=10)



# based on cluster
annotated_filtered_obj$RNA_snn_res.0.3 = factor(annotated_filtered_obj$RNA_snn_res.0.3, 
                                                levels=c(0, 12, 14, 8, 4, 13,
                                                         2, 3, 7, 9, 20,
                                                         6, 5, 1, 10, 17))
DoHeatmap(annotated_filtered_obj, features=genes,
          group.by='RNA_snn_res.0.3', 
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5,
          raster=F) + 
          theme(legend.position='bottom') 
ggsave(file.path(figure_dir, 'Figure3.2_heatmap_cluster.pdf'), height=4, width=8)



Idents(annotated_filtered_obj) = factor(annotated_filtered_obj$RNA_snn_res.0.3,
                                       levels=c(2, 3, 7, 9, 20,
                                                0, 12, 14, 8, 4, 13,
                                                5, 6, 1, 10, 17))
VlnPlot(annotated_filtered_obj, features=genes, idents=levels(annotated_filtered_obj), fill.by='ident',
        stack=T) + 
        theme(legend.position='none')
ggsave(file.path(figure_dir, 'Figure3.2_violin.pdf'), height=6, width=10)

# violin plot
#Idents(annotated_filtered_obj) = factor(annotated_filtered_obj$RNA_snn_res.0.3,
#                                        levels=c(2, 3, 7, 9, 20,
#                                                 0, 12, 14, 8, 4, 13,
#                                                 5, 6, 1, 10, 17))
#sub_counts = GetAssayData(annotated_filtered_obj)[genes,]
#sub_counts = as.matrix(sub_counts)
#sub_counts_df = melt(sub_counts)
#colnames(sub_counts_df) = c('genes', 'barcodes', 'exprs')
#sub_counts_df$cluster = annotated_filtered_obj@meta.data[sub_counts_df$barcodes, 'RNA_snn_res.0.3']
#sub_counts_df$cluster = factor(sub_counts_df$cluster, 
#                               levels=c(2, 3, 7, 9, 20,
#                                        0, 12, 14, 8, 4, 13,
#                                        5, 6, 1, 10, 17))
#g = ggviolin(sub_counts_df, x="genes", y="exprs", fill='cluster') + 
#    facet_wrap('cluster', strip.position="right", ncol=1) + 
#    theme(legend.position='none')
#ggsave(file.path(figure_dir, 'Figure3.2_violin.pdf'), height=6, width=10)


## --- Figure 3.3
marker_dir = file.path(result_dir, 'DEGs')
marker_genes = c()
for (celltype in unique(annotated_filtered_obj@meta.data$annotated_celltype)) {
    DEG_df = read.csv(file.path(marker_dir,  paste0(gsub('/', '', celltype), '_celltype_DEGs.csv')),
                      header=T, row.names=1) 
    #DEgenes = DEG_df[DEG_df$p_val_adj < 0.05, ]$gene
    DEgenes = DEG_df[1:5, 'gene']
    marker_genes = c(marker_genes, unique(DEgenes))
}
DoHeatmap(annotated_filtered_obj, 
          features=marker_genes,
          group.by='annotated_celltype',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5) + 
          theme(legend.position='bottom')
ggsave(file.path(figure_dir, 'Figure3.3_DEGs_celltype_top5_heatmap.pdf'), height=4, width=8)

marker_genes = c()
for (celltype in unique(annotated_filtered_obj@meta.data$annotated_celltype)) {
    DEG_df = read.csv(file.path(marker_dir,  paste0(gsub('/', '', celltype), '_celltype_DEGs.csv')),
                      header=T, row.names=1) 
    DEgenes = DEG_df[DEG_df$p_val_adj < 0.05, ]$gene
    #DEgenes = DEG_df[1:5, 'gene']
    marker_genes = c(marker_genes, unique(DEgenes))
}
DoHeatmap(annotated_filtered_obj, 
          features=marker_genes,
          group.by='annotated_celltype',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5) + 
          theme(legend.position='bottom')
ggsave(file.path(figure_dir, 'Figure3.3_DEGs_celltype_FDR_heatmap.pdf'), height=8, width=8)

marker_genes = c()
for (cluster in unique(annotated_filtered_obj@meta.data$RNA_snn_res.0.3)) {
    DEG_df = read.csv(file.path(marker_dir,  paste0(gsub('/', '', cluster), '_cluster_DEGs.csv')),
                      header=T, row.names=1) 
    #DEgenes = DEG_df[DEG_df$p_val_adj < 0.05, ]$gene
    DEgenes = DEG_df[1:5, 'gene']
    marker_genes = c(marker_genes, unique(DEgenes))
}
DoHeatmap(annotated_filtered_obj, 
          features=marker_genes,
          group.by='RNA_snn_res.0.3',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5) + 
          theme(legend.position='bottom')
ggsave(file.path(figure_dir, 'Figure3.3_DEGs_cluster_top5_heatmap.pdf'), height=6, width=8)

marker_genes = c()
for (cluster in unique(annotated_filtered_obj@meta.data$RNA_snn_res.0.3)) {
    DEG_df = read.csv(file.path(marker_dir,  paste0(gsub('/', '', cluster), '_cluster_DEGs.csv')),
                      header=T, row.names=1) 
    DEgenes = DEG_df[DEG_df$p_val_adj < 0.05, ]$gene
    marker_genes = c(marker_genes, unique(DEgenes))
}
DoHeatmap(annotated_filtered_obj, 
          features=marker_genes,
          group.by='RNA_snn_res.0.3',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5) + 
          theme(legend.position='bottom')
ggsave(file.path(figure_dir, 'Figure3.3_DEGs_cluster_FDR_heatmap.pdf'), height=8, width=8)


# inhibitory neurons - subclusters
inhibitory_obj = subset(annotated_filtered_obj, subset=annotated_celltype == 'Inhibitory neuron')
inhibitory_obj$RNA_snn_res.0.3 = factor(inhibitory_obj$RNA_snn_res.0.3, levels = unique(inhibitory_obj$RNA_snn_res.0.3))
inhibitory_obj$condition = factor(inhibitory_obj$condition, levels=unique(inhibitory_obj$condition))
marker_genes = c()
for (cluster in unique(inhibitory_obj@meta.data$RNA_snn_res.0.3)) {
    DEG_df = read.csv(file.path(marker_dir,  paste0(gsub('/', '', cluster), '_cluster_DEGs.csv')),
                      header=T, row.names=1) 
    DEG_df = DEG_df[DEG_df$cluster == 'control' & DEG_df$pct.1 > 0.1 & DEG_df$pct.2 > 0.1, ]
    DEG_df = DEG_df[order(-abs(DEG_df$avg_log2FC)), ]
    DEgenes = DEG_df[DEG_df$p_val_adj < 0.05, ]$gene
    cat('Cluster:', cluster, 'DEGs:', length(DEgenes), '\n')
    num_of_DEGs = ifelse(cluster == 0, 20, 10)
    if (length(DEgenes) < num_of_DEGs) {
        marker_genes = c(marker_genes, DEgenes[!duplicated(DEgenes)])
    } else {
        top10_genes = DEgenes[1:num_of_DEGs]
        marker_genes = c(marker_genes, top10_genes[!duplicated(top10_genes)])
    }
}
dedup_marker_genes = marker_genes[!duplicated(marker_genes)]
mat = GetAssayData(inhibitory_obj)
mat = as.matrix(mat[dedup_marker_genes, ])
melted_mat = melt(mat)
colnames(melted_mat) = c('genes', 'barcodes', 'exprs')
melted_mat$condition = inhibitory_obj@meta.data[melted_mat$barcodes, 'condition']
melted_mat$cluster = inhibitory_obj@meta.data[melted_mat$barcodes, 'RNA_snn_res.0.3']
avg_exprs_df = melted_mat %>% 
    group_by(genes, cluster, condition) %>% 
    summarise(mean_exprs=mean(exprs)) %>%
    as.data.frame()
avg_exprs_mat = acast(avg_exprs_df, genes ~ cluster + condition)
scaled_mat = t(apply(avg_exprs_mat, 1, scale))
colnames(scaled_mat) = colnames(avg_exprs_mat)
cols = brewer.pal(10, 'Set3')
cluster_cols = cols[1:length(unique(inhibitory_obj$RNA_snn_res.0.3))]
names(cluster_cols) = unique(inhibitory_obj$RNA_snn_res.0.3)
condition_cols = c('red', 'blue')
names(condition_cols) = unique(inhibitory_obj$condition)
col_df = data.frame(cluster=sapply(strsplit(colnames(avg_exprs_mat), split='_'), '[', 1),
                    condition=sapply(strsplit(colnames(avg_exprs_mat), split='_'), '[', 2))
rownames(col_df) = colnames(avg_exprs_mat)
col_df$cluster = as.numeric(col_df$cluster)

pdf(file.path(figure_dir, 'Figure3.3_inhibitory_DEGs_cluster_top10_aggregated_heatmap.pdf'), width=6, height=8)
col_ha = HeatmapAnnotation(df=col_df,
                           col=list(cluster=cluster_cols,
                                    condition=condition_cols))
h = Heatmap(scaled_mat, name='Avg Exprs',
            show_row_names=T, show_column_names=F,
            cluster_rows=T, cluster_columns=F, column_split=col_df$cluster,
            top_annotation=col_ha)
draw(h)
dev.off()

#cols = brewer.pal(10, 'Set3')
#cluster_cols = cols[1:length(unique(inhibitory_obj$RNA_snn_res.0.3))]
#names(cluster_cols) = unique(inhibitory_obj$RNA_snn_res.0.3)
#condition_cols = c('red', 'blue')
#names(condition_cols) = unique(inhibitory_obj$condition)
#pdf(file.path(figure_dir, 'Figure3.3_inhibitory_DEGs_cluster_top10_original_scale_heatmap.pdf'), width=8, height=8)
#col_df = inhibitory_obj@meta.data[colnames(mat), c('RNA_snn_res.0.3', 'condition')]
#colnames(col_df) = c('cluster', 'condition')
#col_df = col_df[order(col_df[, 'cluster'], col_df[, 'condition']), ]
#col_ha = HeatmapAnnotation(df=col_df,
#                           col=list(cluster=cluster_cols,
#                                    condition=condition_cols))
#h = Heatmap(mat[, rownames(col_df)], show_row_names=T, show_column_names=F,
#            cluster_rows=F, cluster_columns=F, column_split=col_df$cluster,
#            top_annotation=col_ha)
#draw(h)
#dev.off()



DoMultiBarHeatmap(inhibitory_obj, 
          features=marker_genes[!duplicated(marker_genes)],
          group.by='RNA_snn_res.0.3',
          group.bar=T,
          additional.group.by = 'condition',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5, raster=F) + 
          theme(legend.position='bottom')
ggsave(file.path(figure_dir, 'Figure3.3_inhibitory_DEGs_cluster_top10_heatmap.pdf'), height=4, width=8)


marker_genes = c()
for (cluster in unique(inhibitory_obj@meta.data$RNA_snn_res.0.3)) {
    DEG_df = read.csv(file.path(marker_dir,  paste0(gsub('/', '', cluster), '_cluster_DEGs.csv')),
                      header=T, row.names=1) 
    DEgenes = DEG_df[DEG_df$p_val_adj < 0.05, ]$gene
    marker_genes = c(marker_genes, DEgenes[!duplicated(DEgenes)])
}
DoHeatmap(inhibitory_obj, 
          features=marker_genes[!duplicated(marker_genes)],
          group.by='RNA_snn_res.0.3',
          slot = "scale.data",
          label=F, disp.min=-2.5, disp.max=2.5) + 
          theme(legend.position='bottom')
ggsave(file.path(figure_dir, 'Figure3.3_inhibitory_DEGs_cluster_FDR_heatmap.pdf'), height=12, width=8)

## --- Figure 3.4
inhibitory_obj = subset(annotated_filtered_obj, subset=annotated_celltype == 'Inhibitory neuron')
split_plot = DimPlot(inhibitory_obj, reduction = "umap", pt.size = .1, group.by='RNA_snn_res.0.3', split.by='condition', label=T) + 
    xlim(-4, 4) + ylim(-3, 10)
g = DimPlot(inhibitory_obj, reduction = "umap", group.by='RNA_snn_res.0.3', split.by='condition', pt.size = .1, label=F) + 
    xlim(-4, 4) + ylim(-3, 10)
g = g + split_plot$layers[[2]]
ggsave(file.path(figure_dir, 'Figure3.4_inhibitory_UMAP.pdf'),
    width=10, height=5)





