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

counts <- read.delim(file="liver_forebrain.tsv", header = TRUE, row.names = 1)

annotation <- read.csv(file = "annotation.table.csv", row.names = 1, header = TRUE)

dim(counts)

head(counts)

head(annotation)

counts.matrix <- as.matrix(counts)

annotation$condition <- as.factor(annotation$condition)
annotation$age <- as.factor(annotation$age)

dds <- DESeqDataSetFromMatrix(countData = counts.matrix, colData = annotation, design = ~condition)

keep <- rowSums(counts(dds))>=10

dds <- dds[keep,]

dds <- DESeq(dds)

res <- results(dds)
res

summary(res)

sum(res$padj < 0.05, na.rm = TRUE)

res05 <- results(dds, alpha = 0.05)
summary(res05)

resOrdered <- res[order(res$padj),]
resOrdered

write.csv(as.data.frame(resOrdered), file = "results/liver_vs_forebrain_results.csv")

resSig <- subset(resOrdered, padj<0.05)
write.csv(as.data.frame(resSig), file = "results/liver_vs_forebrain_results_significant.csv")

resultsNames(dds)

resLFCa <- lfcShrink(dds, coef = "condition_liver_vs_forebrain", type = "apeglm")

summary(resLFCa)

levels(annotation$condition)

plotMA(res05, main="MA plot")

par(mfrow=c(1,2))
limits <- c(-5,5)
plotMA(res05, ylim=limits, main="No shrinkage")
plotMA(resLFCa, ylim=limits, main="Apeglm shrinkage")

head(annotation)

dds.mf <- DESeqDataSetFromMatrix(countData = counts.matrix, colData = annotation, design = ~age+condition)
keep <- rowSums(counts(dds.mf)) >= 10
dds.mf <- dds.mf[keep,]


ntd <- normTransform(dds.mf)

vsd <- vst(dds.mf, blind = TRUE)

nt.counts <- assay(ntd)
vsd.counts <- assay(vsd)
dim(vsd.counts)
head(vsd.counts)

msd <- meanSdPlot(assay(ntd))


msd.df <- data.frame(x=msd$px, y=msd$py)
mds.line <- data.frame(rank=msd$rank, sd=msd$sd)

ggplot(msd.df, aes(x=x, y=y)) + geom_bin_2d(bins=50) + 
  geom_line(data = mds.line, aes(rank, sd), color="red4") + 
  xlab("Mean rank") + ylab("St.deviation")+
  labs(title="Log 2 transform")+theme_classic()

msd <- meanSdPlot(assay(vsd))



select <- order(abs(res$log2FoldChange), decreasing = TRUE)[1:100]
pheatmap(nt.counts[select,], 
         cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation, 
         show_rownames = FALSE, show_colnames = TRUE, fontsize_col = 8)

pheatmap(vsd.counts[select,], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
         show_colnames = TRUE, annotation_col = annotation, main = "VSD normalized counts", fontsize_col = 8)

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
