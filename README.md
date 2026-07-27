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

-  Total of 62,754 genes were quantified across the samples.
-  Differential gene expression analysis was performed using DESeq2.
-  Identified **1,041 significantly differentially expressed genes** using an adjusted p-value < 0.05 and |log2FoldChange| > 1.
-  **491 genes were upregulated**.
-  **550 genes were downregulated**.
-  Principal Component Analysis (PCA) showed clear separation between the control and dexamethasone-treated samples.
-  The volcano plot identified significantly differentially expressed genes based on both fold change and statistical significance.
-  A heatmap of the most variable genes revealed clear clustering of samples according to their experimental groups.
-  Gene Ontology (GO) enrichment analysis was performed on the significant differentially expressed genes using the Biological Process ontology.

## Reference Files

The reference genome (GRCh38) and GTF annotation files were used during alignment and gene quantification but are not included in this repository because of their large size.

## Author

Sheetal
