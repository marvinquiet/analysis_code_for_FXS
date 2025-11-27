library(Seurat)
library(Matrix)
library(ggplot2)
set.seed(2022)

project_dir = "/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"

dsc_snATAC_annotation = read.csv(file.path(project_dir, 'Cortex_multiome_ArchR', 'dscATAC_predict_celltypes.csv'),
                                 header=T, row.names=1, sep=',')
rownames(dsc_snATAC_annotation) = gsub('#', '_', rownames(dsc_snATAC_annotation))
pFC_snRNA_annotation = read.csv(file.path(project_dir, 'Cortex_multiome_Seurat', 'pFC_majorcelltypes_predcelltypes.csv'),
                                header=T, row.names=1, sep=',', check.names=F)
rownames(pFC_snRNA_annotation) = gsub('\\.', '-', rownames(pFC_snRNA_annotation))
common_cells = intersect(rownames(dsc_snATAC_annotation), rownames(pFC_snRNA_annotation))

concord_df = data.frame(table(pFC_snRNA_annotation[common_cells, 'pred_celltype'], 
      dsc_snATAC_annotation[common_cells, 'pred_celltype']))
colnames(concord_df) = c('snRNA', 'snATAC', 'Freq')
ggplot(concord_df, aes(x=snRNA, y=snATAC, fill=Freq)) + geom_tile() + 
    geom_text(aes(label=Freq), color = "white", size = 4) +
    coord_fixed() 
ggsave(file.path(project_dir, 'Cortex_multiome_joint', 'pFC_snRNA_vs_dsc_snATAC_Cellcano_heatmap.png'))

concord_df = data.frame(table(pFC_snRNA_annotation[common_cells, 'firstround_pred_celltype'], 
      dsc_snATAC_annotation[common_cells, 'pred_celltype']))
colnames(concord_df) = c('snRNA', 'snATAC', 'Freq')
ggplot(concord_df, aes(x=snRNA, y=snATAC, fill=Freq)) + geom_tile() + 
    geom_text(aes(label=Freq), color = "white", size = 4) +
    coord_fixed() 
ggsave(file.path(project_dir, 'Cortex_multiome_joint', 'pFC_snRNA_MLP_vs_dsc_snATAC_heatmap.png'))



