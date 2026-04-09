#!/bin/bash
#SBATCH --job-name=trimming
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/alice/t/tlti1/hisat2_pipeline/logs/trimming_%j.out

module load trimmomatic/0.39-tfj436w

mkdir -p /scratch/alice/t/tlti1/hisat2_pipeline/trimmed

for R1 in /scratch/alice/t/tlti1/hisat2_pipeline/fastq/*_1.fastq; do

    SAMPLE=$(basename "$R1" _1.fastq)
    R2=/scratch/alice/t/tlti1/hisat2_pipeline/fastq/${SAMPLE}_2.fastq

    trimmomatic PE \
        -threads 8 \
        "$R1" "$R2" \
        /scratch/alice/t/tlti1/hisat2_pipeline/trimmed/${SAMPLE}_1_paired.fastq \
        /scratch/alice/t/tlti1/hisat2_pipeline/trimmed/${SAMPLE}_1_unpaired.fastq \
        /scratch/alice/t/tlti1/hisat2_pipeline/trimmed/${SAMPLE}_2_paired.fastq \
        /scratch/alice/t/tlti1/hisat2_pipeline/trimmed/${SAMPLE}_2_unpaired.fastq \
        ILLUMINACLIP:/cm/shared/spack/opt/spack/linux-rocky9-x86_64_v3/gcc-12.3.0/trimmomatic-0.39-tfj436wtywk6jxr35a3pzf3qtnmmf2d6/share/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done

