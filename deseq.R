args <- commandArgs(TRUE)

# counts
raw.counts <- read.table("genematrix_rawcounts.txt", header=T, row.names=1, sep='\t')

#conditions
col.data <- read.table("sample.txt", header=T, row.names=1, sep='\t')

#erase NA
raw.counts[is.na(raw.counts)] <- 0


####### DeSeq2 #########


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = raw.counts, colData = col.data, design =~Condition)


dds <- DESeq(dds)

res <- results(dds, alpha=0.05)

res


# MA plot (padj < 0.05) saving it as PDF
pdf('MAplot.pdf')
plotMA(res, main="MA plot - FDR at .5%")
dev.off()

# Volcano plot (padj < 0.05)
pdf('volcanoPlot.pdf')
## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  require(calibrate)
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=rownames(res), cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

volcanoplot(res, main="Volcano plot - FDR at .5%", legendpos = 'topright')
dev.off()
