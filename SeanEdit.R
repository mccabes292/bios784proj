#library(BiocInstaller)
#biocLite("curatedOvarianData")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(curatedOvarianData)
library(Biobase)
library(matrixStats)
data(package="curatedOvarianData")

data("GSE9891_eset")

rM=colMeans(log2(t(exprs(GSE9891_eset)+1)))
rVar=rowVars(log(exprs(GSE9891_eset)+1))


#Grade, Stage (Summary vs. tumor), substage, age, tax, recurrence status, site of first recurrence,primary therapy outcome success
#Debulking, % normal stromal tumor cells , batch


outMeans=kmeans(t(exprs(GSE9891_eset)),centers=6)
head(outMeans$size)


outPC=prcomp(t(exprs(GSE9891_eset)))


d=exprs(GSE9891_eset)

results = ConsensusClusterPlus(d,maxK=6,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
table(results[[6]][["consensusClass"]])
icl = calcICL(results)

