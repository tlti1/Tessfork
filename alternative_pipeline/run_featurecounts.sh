#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/alice/t/tlti1/hisat2_pipeline/logs/featurecounts_%j.out

source ~/.bashrc
conda activate rnaseq

featureCounts \
    -T 8 \
    -p \
    -a /scratch/alice/t/tlti1/hisat2_pipeline/genome/gencode.v22.protein_coding_only.gtf \
    -o /scratch/alice/t/tlti1/hisat2_pipeline/counts/all_cells_counts.txt \
    -g gene_name \
    --extraAttributes gene_id \
    /scratch/alice/t/tlti1/hisat2_pipeline/aligned/*_sorted.bam

echo "featureCounts complete"
