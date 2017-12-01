#library(BiocInstaller)
#biocLite("curatedOvarianData")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(curatedOvarianData)
library(Biobase)

data(package="curatedOvarianData") #pulls all the datasets

data(GSE9891_eset) #pulls data we are using


library(matrixStats)
library(dplyr)


#gene filtering: the data was already log transformed
#Genes with log expression values of <7 and a variance of
  #<0.5 were filtered out before clustering, leaving 8,732 probe sets
rM=rowMeans((exprs(GSE9891_eset)))
rVar=rowVars(exprs(GSE9891_eset))

length(rM)
length(rVar)

#Clinical Variables
#Grade, Stage (Summary vs. tumor), substage, age, tax, recurrence status, site of first recurrence,primary therapy outcome success
#Debulking, % normal stromal tumor cells , batch


table(rM<7,rVar<0.5)


#k Means Clustering Begins Here
d=exprs(GSE9891_eset)
dim(d)
d2 = d[which(rM>=7 | rVar >= 0.5),] #filtering out the genes
dim(d2)

#Run Consensus Cluster Plus to obtain consensus set
results = ConsensusClusterPlus(d2,maxK=6,reps=1000,pItem=1,pFeature=1,clusterAlg="kmdist",distance="pearson",seed=12)
table(results[[6]][["consensusClass"]])
icl = calcICL(results)
icl6=(icl$itemConsensus)[(icl$itemConsensus)$k==6,]

icl6%>%group_by(item)%>%filter( itemConsensus==max(itemConsensus))->outICL6
sum(outICL6$itemConsensus<0.7)

#Subset to consensus set
#PROPOSAL: Subset to 70% instead
consenSet=outICL6$item[outICL6$itemConsensus>0.7]
#Obtain clusters for consensus set
yClust=outICL6$cluster[outICL6$itemConsensus>0.7]
names(yClust)=consenSet

trainSet=d2[,consenSet]

#Find NonConsensus Set
nonConSet=outICL6$item[outICL6$itemConsensus<=0.7]
testSet=d2[,nonConSet]

#Run K Nearest Neighbor Algorithm
library(class)

#Calculates error rate for specified k value on our training set. 
#low and high are the lower and upper bounds of k to be used
#finding the appropriate k to be used 
findBestK=function(low,high,by=1){
  error=NULL
  
for(k in seq(low,high,1)){
    crossVal=knn.cv(t(trainSet),yClust,k=k)
    error[k]<-mean(crossVal!=yClust) #proportion of misclassification, since crossVal! != yClust is 0 or 1
  }
  return(error)
}

kLow=1
kHigh=5
kError=findBestK(kLow,kHigh)
chosenK=((kLow:kHigh)[(kError==min(kError))])[1] #choose k to minimize the error


##Find KNN Classification of test data set
  #prediction of a value in our noncensenus set
knnOut=knn(train=t(trainSet),test=t(testSet), cl=yClust,k=chosenK)



####LDA####


library(HiDimDA)
outLda=Dlda(t(trainSet),as.factor(yClust),ldafun="classification")
p1=predict(outLda,t(testSet))




classVec=(knnOut==(as.numeric(p1$class)))
names(classVec)=colnames(testSet)

classVec2=rep(1,length(consenSet))
names(classVec2)=consenSet

classVec3=ifelse(classVec==T,2,3)



clusterMem=ifelse(classVec==TRUE, knnOut, 0)
clusterMem2=c(clusterMem,yClust)
clusterMem3=clusterMem2[colnames(GSE9891_eset)]


classification=c(classVec3,classVec2)
classification2=classification[colnames(GSE9891_eset)]


##Classification Variable is 1 if >0.7 in Classification Cluster +, 2 if classified via KNN and Dlda, and 3 if NC
finalSE=GSE9891_eset[which(rM>=7 | rVar >= 0.5),]
finalSE$classification=classification2
finalSE$clusterMem=clusterMem3

save(finalSE,file="finalSet")







#Table 1 in paper: See how many of each sample type (LMP and Malignant) are in the cluster
t=pData(GSE9891_eset)
t2=t[consenSet,]
table(t2$sample_type,outICL6$cluster[outICL6$itemConsensus>0.8])