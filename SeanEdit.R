#library(BiocInstaller)
#biocLite("curatedOvarianData")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(curatedOvarianData)
library(Biobase)
library(matrixStats)
library(dplyr)
data(package="curatedOvarianData")

data("GSE9891_eset")

rM=rowMeans((exprs(GSE9891_eset)))
rVar=rowVars((exprs(GSE9891_eset)))


#Grade, Stage (Summary vs. tumor), substage, age, tax, recurrence status, site of first recurrence,primary therapy outcome success
#Debulking, % normal stromal tumor cells , batch





d=(exprs(GSE9891_eset)+1)
rv=rowVars(d)
d2=d[order(rv,decreasing = T)[1:5000],]

results = ConsensusClusterPlus(d2,maxK=6,reps=1000,pItem=0.8,pFeature=1,clusterAlg="kmdist",distance="pearson",seed=1262118388.71279)
table(results[[6]][["consensusClass"]])
icl = calcICL(results)
icl6=(icl$itemConsensus)[(icl$itemConsensus)$k==6,]

icl6%>%group_by(item)%>%filter( itemConsensus==max(itemConsensus))->outICL6
sum(outICL6$itemConsensus<0.8)

#Subset to consensus set
consenSet=outICL6$item[outICL6$itemConsensus>0.8]
#Obtain clusters for consensus set
yClust=outICL6$cluster[outICL6$itemConsensus>0.8]

library(class)
trainSet=d[,consenSet]
crossVal=knn.cv(t(trainSet),yClust,k=3)






outPC=prcomp(t(exprs(GSE9891_eset)))

t=pData(GSE9891_eset)
t2=t[consenSet,]
table(t2$sample_type,outICL6$cluster[outICL6$itemConsensus>0.8])

