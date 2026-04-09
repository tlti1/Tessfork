#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/alice/t/tlti1/hisat2_pipeline/logs/hisat2_align_%j.out

module load hisat2/2.2.1-eovlo7b
module load samtools/1.18

for R1 in /scratch/alice/t/tlti1/hisat2_pipeline/trimmed/*_1_paired.fastq; do

    SAMPLE=$(basename "$R1" _1_paired.fastq)
    R2=/scratch/alice/t/tlti1/hisat2_pipeline/trimmed/${SAMPLE}_2_paired.fastq

    hisat2 -p 8 -x /scratch/alice/t/tlti1/hisat2_pipeline/index/grch38_hisat2_index --dta -1 "$R1" -2 "$R2" 2>/scratch/alice/t/tlti1/hisat2_pipeline/logs/${SAMPLE}_hisat2.log | samtools view -bS | samtools sort -o /scratch/alice/t/tlti1/hisat2_pipeline/aligned/${SAMPLE}_sorted.bam -@ 8

    samtools index /scratch/alice/t/tlti1/hisat2_pipeline/aligned/${SAMPLE}_sorted.bam

done

echo "Alignment complete"
