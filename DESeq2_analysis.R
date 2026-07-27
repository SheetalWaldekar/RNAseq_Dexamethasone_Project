library(DESeq2)
counts <- read.table(
  "//wsl.localhost/Ubuntu-24.04/home/shital/my_project/RNASeq2/data/counts/gene_counts.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)
dim(counts)
colnames(counts)
colSums(counts[,7:14])
dim(counts)
countData <- counts[,7:ncol(counts)]
rownames(countData) <- counts$Geneid 
countData <- as.matrix(countData)
dim(countData)

sampleInfo <- data.frame(
  row.names = colnames(countData),
  condition=c(
    "untreated","dex-treated","untreated","dex-treated",
    "untreated","dex-treated","untreated","dex-treated"
  )
)
sampleInfo
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = sampleInfo,
  design = ~ condition
)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds)
summary(res)
dim(res)
resultsNames(dds)

sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nSignificant DE genes:", nrow(sig), "\n")
cat("upregulated:", sum(sig$log2FoldChange > 0), "\n")
cat("downregulated:", sum(sig$log2FoldChange < 0), "\n")
head(sig)

sum(res$padj<0.05, na.rm=TRUE)
sum(!is.na(res$padj))

levels(dds$condition)
write.csv(as.data.frame(sig), "significant_gene.csv")

library(org.Hs.eg.db)
library(AnnotationDbi)

ensembl_ids <- rownames(sig)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
sig$GeneSymbol <- gene_symbols
sig
nrow(sig)

write.csv(as.data.frame(sig),"sig_genesymbol.csv")

sig_symbol<-sig[!is.na(sig$GeneSymbol), ]
sig_symbol<-sig_symbol[order(sig_symbol$padj), ]
top20 <- sig_symbol[order(sig_symbol$padj), ][1:20, ]
top20_df <- as.data.frame(top20[, c("GeneSymbol","log2FoldChange","padj")])
print(top20_df)
head(rownames(sig), 20)

upregulated <- subset(
  sig,
  padj < 0.05 & log2FoldChange < -1
)
upregulated
dex_up <- upregulated[order(upregulated$log2FoldChange), ]
head(dex_up, 40)
top100 <- dex_up[order(dex_up$log2FoldChange), ][100:200, ]
top100_df <- as.data.frame(top100[, c("GeneSymbol","log2FoldChange","padj")])
print(top100_df)

downregulated <- subset(
  sig,
  padj < 0.05 & log2FoldChange > 0
)
downregulated

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

genes<-rownames(sig)
gene.df <- bitr(
  genes,
  fromType = "ENSEMBL",
  toType = c("ENTREZID","SYMBOL"),
  OrgDb = org.Hs.eg.db
)
length(genes)
go <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
#check the result
go
summary(go)
slotNames(go)
length(gene.df$ENTREZID)
length(unique(gene.df$ENTREZID))
head(go)
length(go@gene)

barplot(go, showCategory = 20)
dotplot(go, showCategory = 20)
cnetplot(go)

write.csv(as.data.frame(go),"go_result.csv")

library(ggplot2)
#volcano plot
df <- as.data.frame(res)
df$sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1
df

plot(res$log2FoldChange,
     -log10(res$pvalue),
     pch = 20,
     xlab = "log2FoldChange",
     ylab = "-Log10(P-value)"
)

# MA plot
plotMA(res, ylim=c(-5,+5))

# Histogram of p-values 
hist(res$pvalue,
     breaks=30,
     col="skyblue",
     main= "P-value distribution"
)

#PCA plot
vsd <- vst(dds,blind = FALSE)
vsd
plotPCA(vsd, intgroup="condition")

#heatmap
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)
pheatmap(sampleDistMatrix)



















