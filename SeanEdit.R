#library(BiocInstaller)
#biocLite("curatedOvarianData")
library(curatedOvarianData)
library(Biobase)
data(package="curatedOvarianData") #pulls all the datasets

data(GSE9891_eset)
head(GSE9891_eset)
phenoData(GSE9891_eset)$sample_type
table(GSE9891_eset$sample_type, GSE9891_eset$histological_type)

head(exprs(GSE9891_eset))[,1:5]

colnames(phenoData(GSE9891_eset))

library(matrixStats)
rv <- rowVars(exprs(GSE9891_eset))
o <- order(rv,decreasing=TRUE)
dists <- dist(t(exprs(GSE9891_eset)[head(o,500),]))

hc <- hclust(dists)
dend <- as.dendrogram(hc)

library(dendextend)
library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
o.dend <- order.dendrogram(dend)
#labels(dend) <- GSE9891_eset$sample[o.dend]
#labels_colors(dend) <- as.integer(GSE9891_eset$sample[o.dend])
plot(dend)