#!/bin/bash
#SBATCH --job-name=cell-ranger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=100G
#SBATCH --output=cell-rang.out

export PATH=/mnt/scratch/wisnoskilab/shared/bioinformatics/sc-RNAseq/cellranger-10.1.0:$PATH

cellranger count \
    --id=run_count_1kpbmcs \
    --fastqs=/mnt/scratch/wisnoskilab/dc2484/GitHub/Single-Cell_RNA_Seq/run_cellranger_count/pbmc_1k_v3_fastqs \
    --sample=pbmc_1k_v3 \
    --transcriptome=/mnt/scratch/wisnoskilab/dc2484/GitHub/Single-Cell_RNA_Seq/run_cellranger_count/Ref_genome \
    --create-bam=true \
    --localcores=80 \
    --localmem=90