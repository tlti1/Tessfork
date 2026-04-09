#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/alice/t/tlti1/hisat2_pipeline/logs/fastqc_%j.out

module load fastqc/0.12.1-hkgpcde

fastqc /scratch/alice/t/tlti1/hisat2_pipeline/fastq/*.fastq \
    -o /scratch/alice/t/tlti1/hisat2_pipeline/fastqc_out/ \
    -t 8
