#library(BiocInstaller)
#biocLite("curatedOvarianData")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(curatedOvarianData)
library(Biobase)

data(package="curatedOvarianData") #pulls all the datasets

data(GSE9891_eset) #pulls data we are using
head(GSE9891_eset)
phenoData(GSE9891_eset)$sample_type
table(GSE9891_eset$sample_type, GSE9891_eset$histological_type)

head(exprs(GSE9891_eset))[,1:5] #view a small bit of the count data matrix

colnames(phenoData(GSE9891_eset)) #all the clinical variables


library(matrixStats)
library(dplyr)


#gene filtering: the data was already log transformed
#Genes with log expression values of <7 and a variance of
  #<0.5 were filtered out before clustering, leaving 8,732 probe sets
rM=colMeans(t(exprs(GSE9891_eset)))
rVar=rowVars(exprs(GSE9891_eset))

length(rM)
length(rVar)

#Clinical Variables
#Grade, Stage (Summary vs. tumor), substage, age, tax, recurrence status, site of first recurrence,primary therapy outcome success
#Debulking, % normal stromal tumor cells , batch





#k means clustering
d=exprs(GSE9891_eset)
dim(d)
d2 = d[which(rM>=7 | rVar >= 0.5),] #filtering out the genes

dim(d2)

#Run Consensus Cluster Plus to obtain consensus set
results = ConsensusClusterPlus(d2,maxK=6,reps=1000,pItem=0.8,pFeature=0.8,clusterAlg="kmdist",distance="pearson",seed=1262118388.71279)
table(results[[6]][["consensusClass"]])
icl = calcICL(results)
icl6=(icl$itemConsensus)[(icl$itemConsensus)$k==6,]

icl6%>%group_by(item)%>%filter( itemConsensus==max(itemConsensus))->outICL6
sum(outICL6$itemConsensus<0.8)

#Subset to consensus set
consenSet=outICL6$item[outICL6$itemConsensus>0.8]
#Obtain clusters for consensus set
yClust=outICL6$cluster[outICL6$itemConsensus>0.8]

trainSet=d2[,consenSet]

#Find NonConsensus Set
nonConSet=outICL6$item[outICL6$itemConsensus<=0.8]
testSet=d2[,nonConSet]

#Run K Nearest Neighbor Algorithm
library(class)
crossVal=knn.cv(t(trainSet),yClust,k=3)

#Calculates error rate for specified k value on our training set. 
#low and high are the lower and upper bounds of k to be used
findBestK=function(low,high,by=1){
  error=NULL
  
for(k in seq(low,high,1)){
    crossVal=knn.cv(t(trainSet),yClust,k=k)
    error[k]<-mean(crossVal!=yClust)
  }
  return(error)
}

kLow=1
kHigh=5
kError=findBestK(kLow,kHigh)
chosenK=((kLow:kHigh)[(kError==min(kError))])[1]


##Find KNN Classification of test data set
knnOut=knn(train=t(trainSet),test=t(testSet), cl=yClust,k=chosenK)



####LDA####


library(HiDimDA)
outLda=DACrossVal(t(trainSet),as.factor(yClust),TrainAlg = Dlda,kfold=2)

out1=Dlda(t(trainSet),as.factor(yClust))
#outPC=prcomp(t(exprs(GSE9891_eset)))

#Table 1 in paper: See how many of each sample type (LMP and Malignant) are in the cluster
t=pData(GSE9891_eset)
t2=t[consenSet,]
table(t2$sample_type,outICL6$cluster[outICL6$itemConsensus>0.8])