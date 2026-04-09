## HISAT2-based scRNA-seq Pipeline
## Replication of Camp et al. (2015) Figure 3D
### University of Leicester — CO7093 Bioinformatics Coursework

---

## Overview

This repository contains an alternative RNA-seq pipeline to replicate Figure 3D from:

> Camp, J.G. et al. (2015). Human cerebral organoids recapitulate gene expression programs of fetal neocortex development. *PNAS*, 112(51), 15672-15677. https://doi.org/10.1073/pnas.1520760112

The original paper used TopHat + Bowtie2 for alignment and Cufflinks for quantification. This pipeline replaces those tools with HISAT2 (alignment) and featureCounts (quantification), followed by Scanpy for clustering analysis.

## Data

Raw sequencing data: GEO accession [GSE75140](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75140)

- 734 single cells (508 organoid, 226 fetal neocortex)
- Paired-end reads, 96-100bp
- Illumina HiSeq 2500
- Fluidigm C1 platform

Reference genome: GENCODE GRCh38 release 22
- Genome FASTA: `GRCh38.primary_assembly.genome.fa`
- GTF annotation: `gencode.v22.primary_assembly.annotation.gtf`
- Available at: https://www.gencodegenes.org/human/release_22.html

---

## Pipeline Overview

    Raw FASTQ (GEO/SRA)
    Quality control (FastQC + MultiQC) ↓
    Adapter trimming (Trimmomatic) ↓
    Genome indexing (HISAT2) ↓
    Alignment (HISAT2) ↓
    BAM sorting and indexing (SAMtools) ↓
    Quantification (featureCounts) ↓
    Clustering and visualisation (Scanpy) ↓
    ROC test validation (sklearn) ↓
    Pipeline comparison (vs paper's data)


---

## Scripts

Run in this order on an HPC cluster (SLURM):

| Script | Description |
| `run_fastqc.sh` | Quality control on raw FASTQ files |
| `run_trimming.sh` | Adapter trimming with Trimmomatic |
| `run_genome_download.sh` | Download GRCh38 genome and GENCODE v22 GTF |
| `run_hisat2_index.sh` | Build HISAT2 genome index |
| `run_alignment.sh` | Align trimmed reads to genome |
| `run_featurecounts.sh` | Quantify gene expression |
| `scanpy_analysis.py` | Clustering, t-SNE, and figure generation |
| `run_roc_test.py` | ROC test for cluster validation |
| `compare_pipelines.py` | Compare our pipeline to paper's expression data |

---

## Software Versions

| Tool | Version |
| FastQC | 0.12.1 |
| Trimmomatic | 0.39 |
| HISAT2 | 2.2.1 |
| SAMtools | 1.18 |
| featureCounts (Subread) | 2.1.1 |
| Scanpy | latest |
| Python | 3.14 |
| sklearn | latest |

---
