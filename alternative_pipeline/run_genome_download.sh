#!/bin/bash
#SBATCH --job-name=genome_download
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/alice/t/tlti1/hisat2_pipeline/logs/genome_download_%j.out

cd /scratch/alice/t/tlti1/hisat2_pipeline/genome/

# Download genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/GRCh38.primary_assembly.genome.fa.gz

# Download GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.primary_assembly.annotation.gtf.gz

# Decompress both
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v22.primary_assembly.annotation.gtf.gz

# Filter GTF to protein coding genes only
grep 'gene_type "protein_coding"' gencode.v22.primary_assembly.annotation.gtf > gencode.v22.protein_coding_only.gtf

echo "Download and filtering complete"
