#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1
#SBATCH --partition=largemem
#SBATCH -o SCENIC-aucell-combined.log

source ~/.bashrc
AUX_DATA_DIR="/panfs/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics/SCENIC_output/auxiliary_data"
PROJECT_DIR="/panfs/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics/SCENIC_output/"
conda activate Seurat

#pyscenic grn $PROJECT_DIR/filtered_scenic.loom $AUX_DATA_DIR/allTFs_mm.txt -o $PROJECT_DIR/adj.csv --num_workers 24

#pyscenic ctx $PROJECT_DIR/adj.csv \
#    $AUX_DATA_DIR/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather $AUX_DATA_DIR/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather $AUX_DATA_DIR/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather $AUX_DATA_DIR/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
#    --annotations_fname $AUX_DATA_DIR/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
#    --expression_mtx_fname $PROJECT_DIR/filtered_scenic.loom \
#    --output $PROJECT_DIR/reg.csv \
#    --mask_dropouts \
#    --num_workers 24

pyscenic aucell \
    $PROJECT_DIR/filtered_scenic.loom \
    $PROJECT_DIR/reg.csv \
    --output $PROJECT_DIR/pyscenic_output.loom \
    --num_workers 24
