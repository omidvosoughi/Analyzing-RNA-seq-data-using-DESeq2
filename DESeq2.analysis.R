# install and load libraries
install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

library(ggplot2)

BiocManager::install("pheatmap")
BiocManager::install("ggsci")
BiocManager::install("apeglm")
BiocManager::install("vsn")

library(pheatmap)
library(ggsci)
library(apeglm)
library(vsn)

# reading data
# RNA-seq count data per gene for 8 samples
counts <- read.delim(file="liver_forebrain.tsv", header = TRUE, row.names = 1)
# information about samples
annotation <- read.csv(file = "annotation.table.csv", row.names = 1, header = TRUE)

dim(counts)

head(counts)

head(annotation)

# transform the data to define the DESeq2 object called dds.
# transform counts into matrix
counts.matrix <- as.matrix(counts)

# the conditions and age should be factors
annotation$condition <- as.factor(annotation$condition)
annotation$age <- as.factor(annotation$age)

# define DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts.matrix, colData = annotation, design = ~condition)

# keep genes with at least 10 reads
keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

# differential expressed genes (and other calculations)
dds <- DESeq(dds)
# the results of DESeq
res <- results(dds)
res

#summarize some basic tallies using the summary function
summary(res)

# Exporting results

# How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm = TRUE)
# change the cuttof for adjusted p-value for the summary
res05 <- results(dds, alpha = 0.05)
summary(res05)

# order the results table by the smallest adjusted p-value
resOrdered <- res[order(res$padj),]
resOrdered

# export the results as a plain-text file
write.csv(as.data.frame(resOrdered), file = "results/liver_vs_forebrain_results.csv")

# export only the results which pass an adjusted p-value threshold
resSig <- subset(resOrdered, padj<0.05)
write.csv(as.data.frame(resSig), file = "results/liver_vs_forebrain_results_significant.csv")

# Log fold change shrinkage for visualization and ranking
# first get the names of the coefficients in the dds object
resultsNames(dds)
# LFC shrinkage
resLFCa <- lfcShrink(dds, coef = "condition_liver_vs_forebrain", type = "apeglm")

summary(resLFCa)

# MA plot

levels(annotation$condition)

plotMA(res05, main="MA plot")

# creates a display with2 windows
par(mfrow=c(1,2))
# set the limits for plotting
limits <- c(-5,5)
plotMA(res05, ylim=limits, main="No shrinkage")
plotMA(resLFCa, ylim=limits, main="Apeglm shrinkage")

# Multiple factor design

head(annotation)

# we can easily include the age variable into the design
dds.mf <- DESeqDataSetFromMatrix(countData = counts.matrix, colData = annotation, design = ~age+condition)
keep <- rowSums(counts(dds.mf)) >= 10
dds.mf <- dds.mf[keep,]

# Data transformation and visualization
# Data normalization and transformation

# this gives log2(n + 1)
ntd <- normTransform(dds.mf)
# variance stabilizing transformation
vsd <- vst(dds.mf, blind = TRUE)
# export transformed counts
nt.counts <- assay(ntd)
vsd.counts <- assay(vsd)
dim(vsd.counts)
head(vsd.counts)

# plot mean vs variance
msd <- meanSdPlot(assay(ntd))


msd.df <- data.frame(x=msd$px, y=msd$py)
mds.line <- data.frame(rank=msd$rank, sd=msd$sd)

ggplot(msd.df, aes(x=x, y=y)) + geom_bin_2d(bins=50) + 
  geom_line(data = mds.line, aes(rank, sd), color="red4") + 
  xlab("Mean rank") + ylab("St.deviation")+
  labs(title="Log 2 transform")+theme_classic()

msd <- meanSdPlot(assay(vsd))

# Data visualization with heatmap
# select genes with larges Log2 Fold Changes
select <- order(abs(res$log2FoldChange), decreasing = TRUE)[1:100]
# plot a heatmap without clustering for the log2-normalized data
pheatmap(nt.counts[select,], 
         cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation, 
         show_rownames = FALSE, show_colnames = TRUE, fontsize_col = 8)

# plot a heatmap without clustering for the vsd data
pheatmap(vsd.counts[select,], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         show_colnames = TRUE, annotation_col = annotation, main = "VSD normalized counts", fontsize_col = 8)

# cluster the rows and columns in the heatmap and cut them into clusters
pheatmap(vsd.counts[select,], cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         annotation_col = annotation, main = "VSD normalized counts", fontsize_col = 6, cutree_cols = 4, cutree_rows = 2)

ann_colors <- list(
  condition=c(liver="palegreen3", forebrain="gold"),
  age=c(day11.5="seashell2", day12.5="seashell4")
)


pheatmap(vsd.counts[select,], cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = TRUE,
         annotation_col = annotation, main = "VSD normalized counts", fontsize_col = 6,
         cutree_cols = 4, cutree_rows = 2,
         color = colorRampPalette(c("royalblue3", "white", "violetred"))(20),
         annotation_colors = ann_colors, filename = "plots/heatmap.vsd.100.pdf")

# PCA Plot
# internal functino in DESeq2
plotPCA(vsd, intgroup=c("condition"), ntop=20000)
ggsave("plots/qq.pdf")

pcaData <- plotPCA(ntd, intgroup=c("condition","age"), ntop=15000, returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, shape=age)) +
  geom_point(size=4, alpha=0.8) +
  xlab(paste0("PC1:", percentVar[1], "% variance")) +
  ylab(paste0("PC2:", percentVar[2], "% variance")) + 
  scale_color_npg() +
  theme_minimal()
ggsave("plots/ss.pdf")
