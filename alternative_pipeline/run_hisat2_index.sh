#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --time=06:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/alice/t/tlti1/hisat2_pipeline/logs/hisat2_index_%j.out

module load hisat2/2.2.1-eovlo7b

mkdir -p /scratch/alice/t/tlti1/hisat2_pipeline/index

hisat2-build -p 8 /scratch/alice/t/tlti1/hisat2_pipeline/genome/GRCh38.primary_assembly.genome.fa /scratch/alice/t/tlti1/hisat2_pipeline/index/grch38_hisat2_index

echo "Index building complete"
