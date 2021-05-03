#install libraries
install.packages("tidyverse")
install.packages("factoextra")
install.packages("cluster")
install.packages("clusterSim")
install.packages("ade4")

#load libraries
library(tidyverse)
library(cluster)
library(clusterSim)
library(factoextra)
library(ade4)

#Load the file
df <- read.table("R_inputfile_BlastNCBI_fecalsamples_for Enterotype_unclFirmiunclSpiro_54animals.txt",
                 header = TRUE,
                 sep = "\t",
                 row.names = 1)


#transpose the data
df <- t(df)

#function for calculating distances

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}


#claculate the distance matrix
data.dist <- dist.JSD(df)


#function for calculating medoid
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}


#run a test on the data, TAKE CARE WITH THE NUMBER OF CLUSTER
data.cluster <- pam.clustering(data.dist, k = 3)


#calculate the number of clusters
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(df),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

#calculate cluster with an alternative method
#elbow-method
fviz_nbclust(df, kmeans, method = "wss")


#Between-class analysis
obs.pca <- dudi.pca(data.frame(t(df)), scannf=F, nf=10)
obs.bet <- bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)

#put labels (animals) on top of points.
labs <- rownames(obs.bet$ls)
s.class(obs.bet$ls,fac=factor(data.cluster),grid=F, xlim=c(-4,4))
text(obs.bet$ls,labels=labs,adj=c(-.1,-.8),cex=0.8)

#lists your clusters with your rowheaders, in this case animals to clusters.
o=order(data.cluster)
data.frame(labs[o],data.cluster[o])

#PCOA
obs.pcoa <- dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)

s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

#filter data

taxa_weight <- obs.bet$tab

write.csv(taxa_weight, file = "taxa_weight.csv")

#ggplot PCoA
pcoa_data <- obs.pcoa$li

ggplot(data = pcoa_data,
       aes(x = A1,
           y = A2,
           color = as.factor(data.cluster))) +
  geom_point(size =3) +
  xlab("PCoA1") +
  xlab("PCoA2") +
  theme_classic()

apply(taxa_weight, 1, function(x) which(x==max(x)))

#gives the names of the genera with highest weight per cluster in cluster 1 the genus with higest weight is in row 9...
colnames(taxa_weight[76])
colnames(taxa_weight[102])
colnames(taxa_weight[9])

