# conda activate Seurat
library(Seurat)
library(LRcell)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)

#result_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics/Cortex_multiome_Seurat_pFC_200_4000_integration"
result_dir = "/beegfs/home/wma36/collaborations/Yulin/FMRpolyG_Cortex_10Xmultiomics/Cortex_multiome_Seurat_pFC_200_4000_integration"
annotated_combined_obj = readRDS(file.path(result_dir, 'Seurat_integration_annotated_object.RDS'))
annotated_combined_obj@meta.data$annotated_final = paste0(annotated_combined_obj$RNA_snn_res.0.6, ':', annotated_combined_obj$annotated_celltype)
annotated_combined_obj@meta.data$annotated_final = factor(annotated_combined_obj@meta.data$annotated_final,
                                                          levels=c('1:Excitatory', '5:Excitatory', '7:Excitatory', '9:Excitatory', '14:Excitatory',
                                                                   '3:Inhibitory', '4:Inhibitory',
                                                                   '2:Oligodendrocytes', '6:Astrocytes', '15:Astrocytes', 
                                                                   '8:Microglia', '11:OPC', '10:Endothelial'))

# control markers
control_obj = subset(annotated_combined_obj, subset=condition=='control')
control_exprs = GetAssayData(object = control_obj, assay = "RNA", slot = "counts")
control_annotations = control_obj$annotated_final
control_enriched_res = LRcell_gene_enriched_scores(expr = control_exprs,
                annot = control_annotations, parallel = FALSE)
saveRDS(control_enriched_res, file.path(result_dir, 'LRcell', 'control_enriched_genes.RDS'))

control_major_annotations = control_obj$annotated_celltype
control_major_enriched_res = LRcell_gene_enriched_scores(expr = control_exprs,
                annot = control_major_annotations, parallel = FALSE)
saveRDS(control_major_enriched_res, file.path(result_dir, 'LRcell', 'control_majorcelltype_enriched_genes.RDS'))

control_markers = get_markergenes(control_enriched_res, method="LR", topn=100)
DEG_df = read.csv(file.path(result_dir, 'bulkDEGs', '3m_CTX_Male_diffexpr_results.csv'))
DEG_pvalues = DEG_df$pvalue
names(DEG_pvalues) = DEG_df$Gene

control_res = LRcell(gene.p = DEG_pvalues, marker.g = control_markers,
                     method='LR')
plot_df = control_res$user
plot_df$cell_type = factor(plot_df$cell_type, 
                           levels=c('1:Excitatory', '5:Excitatory', '7:Excitatory', '9:Excitatory', '14:Excitatory',
                                    '3:Inhibitory', '4:Inhibitory',
                                    '2:Oligodendrocytes', '6:Astrocytes', '15:Astrocytes', 
                                    '8:Microglia', '11:OPC', '10:Endothelial'))
plot_manhattan_enrich(plot_df, sig.cutoff = .05, label.topn = 5) + ggplot2::theme(axis.text.x = element_text(angle = 90)) + 
    ggplot2::labs(x = "Annotated celltypes")
ggsave(file.path(result_dir, 'LRcell', 'LRcell_enriched.png'))

LiR_control_markers = get_markergenes(control_enriched_res, method="LiR", topn=100)
LiR_control_res = LRcell(gene.p = DEG_pvalues, marker.g = LiR_control_markers,  method='LiR')
plot_df = LiR_control_res$user
plot_df$cell_type = factor(plot_df$cell_type, 
                           levels=c('1:Excitatory', '5:Excitatory', '7:Excitatory', '9:Excitatory', '14:Excitatory',
                                    '3:Inhibitory', '4:Inhibitory',
                                    '2:Oligodendrocytes', '6:Astrocytes', '15:Astrocytes', 
                                    '8:Microglia', '11:OPC', '10:Endothelial'))
plot_manhattan_enrich(plot_df, sig.cutoff = .05, label.topn = 5) + ggplot2::theme(axis.text.x = element_text(angle = 90)) + 
    ggplot2::labs(x = "Annotated celltypes")
ggsave(file.path(result_dir, 'LRcell', 'LRcell_enriched_linear_regression.png'))

# control major markers
control_major_markers = get_markergenes(control_major_enriched_res, method="LR", topn=100)
control_major_res = LRcell(gene.p = DEG_pvalues, marker.g = control_major_markers,
                     method='LR')
plot_df = control_major_res$user
plot_df$cell_type = factor(plot_df$cell_type, levels=c('Excitatory', 'Inhibitory',
                                                       'Oligodendrocytes', 'Astrocytes',
                                                       'Microglia', 'OPC', 'Endothelial'))
plot_manhattan_enrich(plot_df, sig.cutoff = .05, label.topn = 5) + ggplot2::theme(axis.text.x = element_text(angle = 90)) + 
    ggplot2::labs(x = "Annotated celltypes")
ggsave(file.path(result_dir, 'LRcell', 'LRcell_major_celltype_enriched.png'))

LiR_control_major_markers = get_markergenes(control_major_enriched_res, method="LiR", topn=100)
LiR_control_major_res = LRcell(gene.p = DEG_pvalues, marker.g = LiR_control_major_markers,  method='LiR')
plot_df = LiR_control_major_res$user
plot_df$cell_type = factor(plot_df$cell_type, levels=c('Excitatory', 'Inhibitory',
                                                       'Oligodendrocytes', 'Astrocytes',
                                                       'Microglia', 'OPC', 'Endothelial'))
plot_manhattan_enrich(plot_df, sig.cutoff = .05, label.topn = 5) + ggplot2::theme(axis.text.x = element_text(angle = 90)) + 
    ggplot2::labs(x = "Annotated celltypes")
ggsave(file.path(result_dir, 'LRcell', 'LRcell_major_celltype_enriched_linear_regression.png'))

# DEGs log2FC correlation
bulk_DEG_df = read.csv(file.path(result_dir, 'bulkDEGs', '3m_CTX_Male_diffexpr_results.csv'))
bulk_DEG_pvalues = bulk_DEG_df$pvalue
names(bulk_DEG_pvalues) = bulk_DEG_df$Gene

DEG_dir = file.path(result_dir, 'DEGs_nothreshold')
DEG_files = list.files(DEG_dir, pattern='*_DEGs.csv')
DEG_markers = list()
for (DEG_file in DEG_files) {
    DEG_cluster = gsub('_condition_DEGs.csv', '', DEG_file)
    celltype = as.character(unique(annotated_combined_obj$annotated_final[annotated_combined_obj$RNA_snn_res.0.6 == DEG_cluster]))
    DEG_df = read.csv(file.path(DEG_dir, DEG_file))
    DEG_markers[[celltype]] = -log10(DEG_df$p_val[1:100])
    names(DEG_markers[[celltype]]) = DEG_df$X[1:100]
}
DEG_LRcell_res = LRcell(gene.p = bulk_DEG_pvalues, marker.g = DEG_markers, method='LiR')
plot_df = DEG_LRcell_res$user
plot_df$cell_type = factor(plot_df$cell_type, 
                           levels=c('1:Excitatory', '5:Excitatory', '7:Excitatory', '9:Excitatory', '14:Excitatory',
                                    '3:Inhibitory', '4:Inhibitory',
                                    '2:Oligodendrocytes', '6:Astrocytes', '15:Astrocytes', 
                                    '8:Microglia', '11:OPC', '10:Endothelial'))
plot_manhattan_enrich(plot_df, sig.cutoff = .05, label.topn = 5) + ggplot2::theme(axis.text.x = element_text(angle = 90)) + 
    ggplot2::labs(x = "Annotated celltypes")
ggsave(file.path(result_dir, 'LRcell', 'LRcell_DEG_enriched.png'))


bulk_DEG_df = read.csv(file.path(result_dir, 'bulkDEGs', '3m_CTX_Male_diffexpr_results.csv'))
DEG_dir = file.path(result_dir, 'DEGs_nothreshold')
DEG_files = list.files(DEG_dir, pattern='*_DEGs.csv')
subtype_DEG_correlation_df = data.frame(matrix(ncol=4, nrow=0))
colnames(subtype_DEG_correlation_df) = c('bulk_lfc', 'sclfc', 'scFDR', 'celltype')
for (DEG_file in DEG_files) {
    DEG_cluster = gsub('_condition_DEGs.csv', '', DEG_file)
    celltype = as.character(unique(annotated_combined_obj$annotated_final[annotated_combined_obj$RNA_snn_res.0.6 == DEG_cluster]))
    DEG_df = read.csv(file.path(DEG_dir, DEG_file))
    match_idx = match(bulk_DEG_df$Gene, DEG_df$X)
    avg_log2FC = DEG_df$avg_log2FC[match_idx]
    FDR = DEG_df$p_val_adj[match_idx]
    df = data.frame(bulk_lfc=bulk_DEG_df$log2FoldChange, sclfc=avg_log2FC,
                    scFDR=FDR, celltype=celltype)
    subtype_DEG_correlation_df = rbind(subtype_DEG_correlation_df, df)
}
subtype_DEG_correlation_df$significant = ifelse(subtype_DEG_correlation_df$scFDR < 0.05, 'sig', 'non-sig')
subtype_DEG_correlation_df_cleaned = na.omit(subtype_DEG_correlation_df)
subtype_DEG_correlation_df_cleaned$logscFDR = -log10(subtype_DEG_correlation_df_cleaned$scFDR)
colors = c('non-sig'='gray', 'sig'='red')
#cor_results = subtype_DEG_correlation_df_cleaned %>%
#    group_by(celltype) %>%
#    summarize(
#              cor = cor(bulk_lfc, sclfc, method = "spearman"),
#              pvalue = format(cor.test(bulk_lfc, sclfc, method = "spearman")$p.value, scientific=T, digits=3),
#              n=n()
#    )
g = ggscatter(subtype_DEG_correlation_df_cleaned, x='bulk_lfc', y='sclfc', color='significant', size=1, shape=20, alpha=0.5) + 
    stat_cor(method='spearman', label.y.npc="top", label.x.npc = "left") + #, label.sep = sprintf(", n = %s, ", nrow(.)))  + 
    geom_hline(yintercept=0, linetype="dotted") + 
    geom_vline(xintercept=0,linetype="dotted") + 
    scale_color_manual(values=colors) + 
    facet_wrap(~celltype, scales = "free", nrow=3) + 
    scale_x_continuous(limits = c(-0.5, 0.5)) + 
    geom_text(data = subtype_DEG_correlation_df_cleaned %>% count(celltype), 
            aes(label = paste("N =",n), y = -Inf, x  = Inf, vjust=-1, hjust=1))
ggsave(file.path(result_dir, 'LRcell', 'bulk_sc_lfc_subtype_correlation.png'), height=9, width=15)


celltype_DEG_dir = file.path(result_dir, 'DEGs_nothreshold_per_celltype')
celltype_DEG_files = list.files(celltype_DEG_dir, pattern='*_condition_DEGs.csv')
celltype_DEG_correlation_df = data.frame(matrix(ncol=4, nrow=0))
colnames(celltype_DEG_correlation_df) = c('bulk_lfc', 'sclfc', 'scFDR', 'celltype')
for (celltype_DEG_file in celltype_DEG_files) {
    celltype = gsub('_condition_DEGs.csv', '', celltype_DEG_file)
    celltype_DEG_df = read.csv(file.path(celltype_DEG_dir, celltype_DEG_file))
    match_idx = match(bulk_DEG_df$Gene, celltype_DEG_df$X)
    avg_log2FC = celltype_DEG_df$avg_log2FC[match_idx]
    #avg_log2FC[is.na(avg_log2FC)] = 0
    #print(celltype_DEG_file)
    #print(cor.test(avg_log2FC, bulk_DEG_df$log2FoldChange, method='spearman'))
    FDR = celltype_DEG_df$p_val_adj[match_idx]
    #FDR[is.na(FDR)] = 1
    df = data.frame(bulk_lfc=bulk_DEG_df$log2FoldChange, sclfc=avg_log2FC,
                    scFDR=FDR, celltype=celltype)
    #print(celltype_DEG_file)
    #print(cor.test(pvalue, bulk_DEG_df$pvalue, method='spearman'))
    #celltype_DEG_markers[[celltype]] = celltype_DEG_df$X
    #celltype_DEG_markers[[celltype]] = -log10(celltype_DEG_df$p_val_adj)
    #names(celltype_DEG_markers[[celltype]]) = celltype_DEG_df$X
    celltype_DEG_correlation_df = rbind(celltype_DEG_correlation_df, df)
}   
celltype_DEG_correlation_df$significant = ifelse(celltype_DEG_correlation_df$scFDR < 0.05, 'sig', 'non-sig')
celltype_DEG_correlation_df_cleaned = na.omit(celltype_DEG_correlation_df)
write.csv(celltype_DEG_correlation_df_cleaned, file.path(result_dir, 'LRcell', 'bulk_sc_lfc_correlation.csv'))
colors = c('non-sig'='gray', 'sig'='red')
for (celltype in unique(celltype_DEG_correlation_df_cleaned$celltype)) {
    celltype_correlation_df = celltype_DEG_correlation_df_cleaned[celltype_DEG_correlation_df_cleaned$celltype == celltype,]
    ggscatter(celltype_correlation_df, x='bulk_lfc', y='sclfc', color='significant', size=1, shape=20, alpha=0.5) + 
        stat_cor(method='spearman', label.y.npc="top", label.x.npc = "left") + #, label.sep = sprintf(", n = %s, ", nrow(.)))  + 
        geom_hline(yintercept=0, linetype="dotted") + 
        geom_vline(xintercept=0,linetype="dotted") + 
        scale_color_manual(values=colors) + 
        #facet_wrap(~celltype, scales = "free", nrow=2) + 
        scale_x_continuous(limits = c(-0.5, 0.5))  + 
        #geom_text(data = celltype_DEG_correlation_df_cleaned %>% count(celltype), 
        #        aes(label = paste("N =",n), y = -Inf, x  = Inf, vjust=-1, hjust=1))
        geom_text(aes(label = paste("N =",nrow(celltype_correlation_df)), y = -Inf, x  = Inf, vjust=-1, hjust=1))
    ggsave(file.path(result_dir, 'LRcell', paste0(celltype, '_bulk_sc_lfc_correlation.png')), height=3, width=3)
}
#ggsave(file.path(result_dir, 'LRcell', 'bulk_sc_lfc_correlation.png'), height=3, width=3)


# heatmap
bulk_DEG_df = read.csv(file.path(result_dir, 'bulkDEGs', '3m_CTX_Male_diffexpr_results.csv'))
rownames(bulk_DEG_df) = bulk_DEG_df$Gene
DEG_dir = file.path(result_dir, 'DEGs_nothreshold')
DEG_files = list.files(DEG_dir, pattern='*_DEGs.csv')
for (DEG_file in DEG_files) {
    DEG_cluster = gsub('_condition_DEGs.csv', '', DEG_file)
    if (!DEG_cluster %in% c('1', '5', '7', '14', '3', '4')) {next}
    celltype = as.character(unique(annotated_combined_obj$annotated_final[annotated_combined_obj$RNA_snn_res.0.6 == DEG_cluster]))
    DEG_df = read.csv(file.path(DEG_dir, DEG_file), row.names=1)
    DEG_df = DEG_df[DEG_df$p_val_adj < 0.05, ]
    DEG_df = DEG_df[order(-DEG_df$avg_log2FC), ]
    gene_list = rownames(DEG_df)
    sub_obj = subset(annotated_combined_obj, subset=annotated_final==celltype)
    sub_obj = sub_obj[gene_list]
    sub_counts = sub_obj@assays$RNA@counts
    # calculate average gene expression
    avg_mat = matrix(nrow=length(gene_list), ncol=length(unique(sub_obj$orig.ident)))
    rownames(avg_mat) = gene_list
    colnames(avg_mat) = unique(sub_obj$orig.ident)
    for (ident in unique(sub_obj$orig.ident)) {
        ident_idx = which(sub_obj$orig.ident == ident)
        avg_counts = rowMeans(sub_counts[, ident_idx])
        avg_mat[, ident] = avg_counts
    }
    avg_mat = t(scale(t(avg_mat)))
    bulk_log2FC = bulk_DEG_df[gene_list, 'log2FoldChange']
    bulk_log2FC[is.na(bulk_log2FC)] = 0
    row_ha = rowAnnotation(bulkFC=bulk_log2FC)
    column_ha = HeatmapAnnotation(condition=as.factor(ifelse(grepl('WT', colnames(avg_mat)), 'WT', 'FXTAS')))
    h = Heatmap(avg_mat, name='sn_avgCounts', top_annotation = column_ha, right_annotation = row_ha,
                cluster_rows = F, cluster_columns = F)
    pdf(file.path(result_dir, paste0(celltype, '_FCheatmap.pdf')))
    draw(h, merge_legend=T)
    dev.off()
}

# heatmap for inhibitory and excitatory
#result_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics/Cortex_multiome_Seurat_pFC_200_4000_integration"
result_dir = "/beegfs/home/wma36/collaborations/Yulin/FMRpolyG_Cortex_10Xmultiomics/Cortex_multiome_Seurat_pFC_200_4000_integration"
annotated_combined_obj = readRDS(file.path(result_dir, 'Seurat_integration_annotated_object.RDS'))
DEGs_dir = file.path(result_dir, 'DEGs_nothreshold_per_celltype')
for (celltype in c('Inhibitory', 'Excitatory')) {
    sub_obj = subset(annotated_combined_obj, subset=annotated_celltype==celltype)
    DEGs_df = read.csv(file.path(DEGs_dir, paste0(celltype, '_condition_DEGs.csv')), header=T, row.names=1)
    selected_DEGs_df = as.data.frame(DEGs_df[1:200, ])
    selected_DEGs_df = selected_DEGs_df[order(-selected_DEGs_df$avg_log2FC), ]
    selected_DEGs = rownames(selected_DEGs_df)
    sub_obj = sub_obj[selected_DEGs]
    sub_counts = sub_obj@assays$RNA@counts
    # calculate average gene expression
    avg_mat = matrix(nrow=length(selected_DEGs), ncol=length(unique(sub_obj$orig.ident)))
    rownames(avg_mat) = selected_DEGs
    colnames(avg_mat) = unique(sub_obj$orig.ident)
    for (ident in unique(sub_obj$orig.ident)) {
        ident_idx = which(sub_obj$orig.ident == ident)
        avg_counts = rowMeans(sub_counts[, ident_idx])
        avg_mat[, ident] = avg_counts
    }
    avg_mat = t(scale(t(avg_mat)))
    write.csv(avg_mat, file.path(result_dir, paste0(celltype, '_FCheatmap_data.csv')))
    #column_ha = HeatmapAnnotation(condition=as.factor(ifelse(grepl('WT', colnames(avg_mat)), 'WT', 'FXTAS')),
    #                              col=list(condition=c('WT'='blue', 'FXTAS'='red')))
    #h = Heatmap(avg_mat, name=paste(celltype, 'sn_avgCounts'), top_annotation = column_ha,
    #            cluster_rows = F, cluster_columns = F, row_names_gp = grid::gpar(fontsize = 6))
    #pdf(file.path(result_dir, paste0(celltype, '_FCheatmap.pdf')), width=6, height=12)
    #draw(h, merge_legend=T)
    #dev.off()
}
