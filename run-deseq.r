## Setting up

setwd("Z:/DataAnalysis/temp requests/deseq_test")
library('DESeq2')
sampleTable<-read.table("merged.counts",header=TRUE, row.names = 1)
SampleNames=colnames(sampleTable)
connames=as.factor(c(rep("2089",2),rep("GFP",2)))
samples <- data.frame(row.names=SampleNames, condition=connames)

## DEseq analysis
ddsHTSeq <- DESeqDataSetFromMatrix(countData = sampleTable, colData=samples, design=~condition)
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref="GFP")

dds <- DESeq(ddsHTSeq)
res <- results(dds)

## Order the results
res<-res[order(res$padj),]
res

## PCA plot using ggplot
rld <- rlogTransformation(dds, blind=TRUE)
data <- plotPCA(rld, intgroup=c("condition"), returnData=T)
percentVar <- round(100 * attr(data, "percentVar"))

library(ggplot2)
pdf('deseq2_pca.pdf', height=8, width=8)
ggplot(data, aes(PC1, PC2, color=condition)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()

## Intersample distance
library("pheatmap")
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

## MA plot
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()


## Plot individual genes:
pdf('ITCH_ENSG00000078747_Exp.pdf', height=4, width=4)
plotCounts(dds, gene="ENSG00000078747", intgroup="condition")
dev.off()


## Save Results
write.table(file='deseq-results.txt', res, quote=F, col.names=T, sep='\t')
norm.counts <-as.data.frame(counts(dds, normalized=T))
write.table(norm.counts, file='norm-counts-all-data.txt', quote=F, col.names=T, sep='\t')

save(res,file="deseq_res.rd")

