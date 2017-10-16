#library(BiocInstaller)
#biocLite("curatedOvarianData")
library(curatedOvarianData)
library(Biobase)
data(package="curatedOvarianData")
data(TCGA_eset)
head(TCGA_eset)

data("GSE2109_eset")
data("GSE26712_eset")

table(GSE2109_eset$sample_type)

table(GSE26712_eset$sample_type)
