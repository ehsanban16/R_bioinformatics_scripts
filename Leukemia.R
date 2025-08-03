setwd("~/Desktop/gene/Leukemia/")

library(limma)
library(GEOquery)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)

#setRepositories() # add 2 BioC software (1 2)
#install.packages(c("GEOquerry", limma, gplots,reshape2,plyr)
res <- read.delim("GSE9476.top.table.tsv")
aml.up <-subset(res, logFC > 1 & adj.P.Val < 0.05)
aml.up.gene <- unique(aml.up$Gene.symbol)

series <- "GSE9476" 
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE)

if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gr <- c("CD34", rep("BM", 10), rep("CD34", 7), rep("AML", 26), rep("PB", 10), rep("CD34", 10))

ex <- exprs(gset) #ex<-log2(ex+1) #ex<-normalizeQuantiles(ex)

#QC box plot comparison of groups
pdf("Results/boxplot.pdf", width=64)
boxplot(ex)
dev.off()

#cor heatmap
pdf("Results/CorHeatMap.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row = gr, labels_col = gr, color=bluered(256), border_color = NA)
dev.off()

#PCA
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

#PCA on change instead of absolute val
ex.scale <- t(scale(t(ex), scale=F))
pc2 <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc2)
plot(pc2$x[,1:2])
dev.off()

pcr <- data.frame(pc2$r[,1:3], group=gr)
pdf("Results/PC_scaled_rots.pdf")
ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=3) + theme_bw()
dev.off()

# Differential expression analysis ------------------------------------- 
# assign samples to groups and set up design matrix
gr <- factor(gr)
gset$description <- factor(gr)
design <- model.matrix(~description + 0, gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design)  # fit linear model
cont.matrix <- makeContrasts(AML-CD34, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val","logFC"))
write.table(tT, "Results/ALM_CD34.txt", row.names=F, sep="\t", quote=F)
## -------

## pathway ----
aml.up <- subset(tT, logFC >1 & adj.P.Val > 0.05)
#aml.up.genes <- unique(aml.up$Gene.symbol)
#aml.up.genes <- sub('///.*' , "", aml.up.gene) #replace w/ method not removing genes
aml.up.genes <- unique(as.character( strsplit2( aml.up$Gene.symbol, "///" ) ))
write.table(aml.up.genes, file="Results/AML_CD34_Up.txt", quote=F, row.names = F, col.names = F)

aml.down <- subset(tT, logFC > -1 & adj.P.Val > 0.05)
aml.down.genes <- unique(as.character( strsplit2( aml.down$Gene.symbol, "///" ) ))
write.table(aml.down.genes, file="Results/AML_CD34_Down.txt", quote=F, row.names = F, col.names = F)

