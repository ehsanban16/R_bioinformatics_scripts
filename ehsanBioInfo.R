library(pheatmap)
library(ggplot2)
library(reshape)
library(grid)

x <- read.csv("Endoderm.csv")
rownames(x) <- x[,1]
x <- x[,-1]
x <- log2(x+1)
x <- na.omit(x)

any(is.na(x))
pdf("heatmap2.pdf")
pheatmap(x, border_color = NA)
dev.off()
t.test(x[1:3, 5], x[4:6, 5])

library(ggplot2)
x <- t(x)
x <- data.frame(x)
x$Gene <- rownames(x) #same as cbind(x, Gene=rownames(x))
ggplot(x, aes( x= DE.1, y= DE.2 )) + geom_point() + geom_text(label=x$Gene)
p <- ggplot(x) + geom_col( aes(y = DE.1, x=Gene, fill= DE.2) )
p <-p + ylab("Expression of gene in Definitive Endoderm (log2)")
pdf("ExpressionBarPlot.pdf")
p
dev.off()
x<-x[,!names(x) %in% "Gene"]

#correlation
y <- data.frame( BP = c( 100, 110, 120) , Weight = c( 160, 170, 100))
cor(y$BP, y$Weight)
ggplot(y, aes(x=BP, y=Weight)) + geom_point() + geom_smooth(method=lm)


#violin
library(reshape)
x.m <- melt(x)
ggplot(x.m , aes(x=variable, y=value)) + geom_boxplot(outlier.size = -1)

#pearson corr heatmap
x <- na.omit(x)
x <- x[-2:-3,]
pheatmap(x, border_color = NA, clustering_distance_cols = "correlation", clustering_distance_rows = "correlation")
ggplot(x , aes( ISL1, HLXb9) ) + geom_point() + geom_smooth(method = "lm")
cor(x[, c("ISL1", "HLXb9")])

#apply
apply(x, 1, var)
var(as.numeric(x[2,]))
apply(x, 2, min, na.rm = T) # 2 is row
Mytest <- function(i){
  t.test(x[4:10, i], x[11:32, i])$p.value
}
sapply(1:9, Mytest)

#standard error
se <- function(x){
  sd(x)/sqrt(length(x))
}

#PCA
pc <- prcomp(x)
pcx <- data.frame(pc$x)
pcr <- data.frame(pc$rotation)
pcr$Gene <- rownames(pcr)
ggplot(pcx, aes(x=PC1, y = PC2)) + geom_point()
pcx$Sample <- rownames(pcx)
pcx$Sample = substr(pcx$Sample, 1, nchar(pcx$Sample) - 2)
ggplot(pcx, aes(x=PC1, y = PC2, color = Sample)) + geom_point()
ggplot(pcr, aes(PC1, PC2, label=Gene)) + geom_text()
x[,3,drop=F]

#ANOVA
x$Sample <- substr(rownames(x) , 1, nchar(rownames(x)) - 2)
y <- x[, c("HLXb9", "Sample")]
ggplot(y, aes(x=Sample, y= HLXb9)) + geom_point()
t.test(y[ y$Sample == "SC", "HLXb9"], y[ y$Sample == "DE", "HLXb9"])
subset(y, y$HLXb9>4)
a <- aov(HLXb9 ~Sample, y)
anova(a)
AOV <- function(Gene){
  y <- x[, c(Gene, "Sample")]
  colnames(y)[1] = "Gene"
  anova(aov(Gene ~Sample, y))[1,5]
}
genes <- genes[genes!="Sample"]
sapply(genes, AOV)
ggplot(x, aes( x= Sample, y= HLXb9 )) + geom_boxplot()

x.t <- t(x)
x.t <- x.t[, -2:-3]
xc <- cor(x.t)
pheatmap(xc)

pc <- prcomp(x)  
pcx <- data.frame(pc$x)
pcx$Sample <- rownames(pcx)
pcx$Sample <- substr( pcx$Sample, 1, nchar(pcx$Sample) -2)
ggplot(pcx, aes(PC1, PC2, color=Sample)) + geom_point()

x.m <- melt(as.matrix(x[3:10,]))
colnames(x.m) <- c("Sample", "Gene", "Exp")
anova( aov(Exp~Gene+Sample, x.m ) )


#nonparam tests
xx <- data.frame(Person=1:7, Pressure1=runif(7,11,12) , Pressure2 = runif(7,10,11))
wilcox.test(xx$Pressure1, xx$Pressure2, paired=T) #Mann-Whitney
y <- data.frame(Mean=colMeans(xx[,-1]), SD=apply(xx[,-1],2,sd))
y$Group <- rownames(y)
se <- function(x) sqrt(sd(x))/ (length(x)-1)
y$SE <- apply(xx[,-1], 2, se)
ggplot(y, aes(Group, Mean, ymin= Mean-SE, ymax=Mean+SE, fill=Group)) + geom_bar(stat="identity") + geom_errorbar()

