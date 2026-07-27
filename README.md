# RNA-Seq Differential Gene Expression Analysis

## Project

This project analyzes RNA-seq data from human airway smooth muscle (HASM) cells treated with dexamethasone to identify differentially expressed genes and enriched biological pathways. The analysis was performed using a standard RNA-seq workflow in Linux and R.

## Objectives

- Perform quality assessment of RNA-seq data
- Generate gene count matrix using featureCounts
- Identify differentially expressed genes using DESeq2
- Visualize expression patterns using PCA, heatmap, and volcano plot
- Perform Gene Ontology (GO) enrichment analysis

## Tools

- Linux
- FastQC
- HISAT2
- Samtools
- featureCounts
- R
- DESeq2
- clusterProfiler

## Workflow

FASTQ files
↓
Quality Control (FastQC)
↓
Alignment (HISAT2)
↓
Gene Quantification (featureCounts)
↓
DE analysis (DESeq2)
↓
PCA
↓
Heatmap
↓
Volcano Plot
↓
GO Enrichment

## Results

- Significant DEGs identified
- PCA separates treated and untreated samples
- Volcano plot visualizes differential expression
- GO enrichment highlights biological processes affected by dexamethasone

## Reference Files

The reference genome (GRCh38) and GTF annotation files were used during alignment and gene quantification but are not included in this repository because of their large size.

## Author

Sheetal




















