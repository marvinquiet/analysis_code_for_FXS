#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --partition=largemem
#SBATCH -o Cortex_FXTAS2.log

source ~/.bashrc
PROJECT_DIR="/projects/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"


## obtain snATAC-seq data
#cellranger-arc mkfastq --id=snATAC \
#    --run=$PROJECT_DIR/BCL/220418_A00327_0240_BHYCYTDRXY \
#    --samplesheet=$PROJECT_DIR/BCL/220418_A00327_0240_BHYCYTDRXY/SampleSheet_corrected.csv   ## there are some problems with the original sample sheet

## obtain snRNA-seq data
#cellranger-arc mkfastq --id=snRNA \
#    --run=$PROJECT_DIR/BCL/220419_A00327_0241_BHYFN7DRXY \
#    --samplesheet=$PROJECT_DIR/BCL/220419_A00327_0241_BHYFN7DRXY/SampleSheet.csv 


## === create corresponding foldes in snATAC-seq according to snRNA-seq output because the generated fastqs are all returned under the fastq_path without organized

## run cellranger arc to count gene expression and generate fragment files
## note that when creating the library files, the Gene expression samples are named as E3 (WT1), G3 (WT2), F3 (cortex_FXTAS1), H3 (cortex_FXTAS2) 
#cellranger-arc count --id=Cortex_multiome_WT1_counts \
#    --reference=$PROJECT_DIR/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
#    --libraries=$PROJECT_DIR/Cortex_WT1_library.csv \
#    --localcores=8 


#cellranger-arc count --id=Cortex_multiome_WT2_counts \
#    --reference=$PROJECT_DIR/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
#    --libraries=$PROJECT_DIR/Cortex_WT2_library.csv \
#    --localcores=8 

#cellranger-arc count --id=Cortex_multiome_FXTAS1_counts \
#    --reference=$PROJECT_DIR/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
#    --libraries=$PROJECT_DIR/Cortex_FXTAS1_library.csv \
#    --localcores=8 

cellranger-arc count --id=Cortex_multiome_FXTAS2_counts \
    --reference=$PROJECT_DIR/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
    --libraries=$PROJECT_DIR/Cortex_FXTAS2_library.csv \
    --localcores=8 


