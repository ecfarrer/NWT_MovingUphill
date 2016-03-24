library(vegan)
library(reshape)
library(bioDist)

#I don't believe this method entirely - when you have a matrix of columns species, treatment then all the samples, by definition, if you have two treatments half of the rows will be zeros (because the sample is only associated with one treatment). this creates nonsense distances in the distance matrix which becmoe NAs but then which Ryan fills in with 1 (equal to correlation of zero). During the permutation test, however, these nonsense 1's will be permuted, so it will do weird things with what your null expectation.

#use rarefied comm.data from original cooccurrence script, make sure comm.data$lomehi is a factor
comm.data

#take out doubletons and singletons
comm.data.sub<-data.frame(comm.data[,1:2],comm.data[,(which(colSums(comm.data[,28:2209])>2)+27)])

#melt data set and recast so that samples are columns and microbes are metadata are rows
comm.data.melt<-melt(comm.data.sub, id=c("lomehi","X.SampleID"))
comm.data.cast<-cast(comm.data.melt, lomehi+variable~X.SampleID, value="value",add.missing=TRUE,fun.aggregate=sum)

#make the distance matrix, this will take a while.  Make sure columns of metadata are not included
microbes.dist<-spearman.dist(data.matrix(comm_data_sub_cast[,-c(1:3)]))

#NA's will be removed and replaced with 1's. They are equal to correlations of 0 (1-0=0 for spearman correlation)
for(i in 1:length(microbes.dist)){
  
  if(is.na(microbes.dist[i])==TRUE){microbes.dist[i]<-1}
}

#run a permanova with trt as the independent variable
adonis(microbes.dist~comm_data_sub_cast$trt,data=NULL)