######Edge creation######

#input is "results" from Pairwise_Cooccurrence_Calcs

library(igraph)
library(fdrtool)
library(ggplot2)


co_occur_pairs<-function(dataset){
  final.results<-data.frame()
  rhos<-c(-.75,-.5,.5,.75)
  trts<-as.vector(unique(dataset$trt))
  
  for(t in 1:length(trts)){
    #t<-1
    dataset_trt<-subset(dataset, trt==trts[t])
    dataset_trt_no0<-subset(dataset_trt, ab1 > 2 & ab2 > 2)#Ryan had this at 0
    dataset_trt_no0$pairs<-paste(dataset_trt_no0$taxa1,dataset_trt_no0$taxa2)
    
    for(r in 1:4){
      #r<-5
      if(rhos[r] < 0){temp<-subset(dataset_trt_no0, rho <= rhos[r])}
      if(rhos[r] > 0){temp<-subset(dataset_trt_no0, rho >= rhos[r])}
      if(dim(temp)[1]>1){
        #temp.graph<-simplify(graph.edgelist(as.matrix(temp[,c(2,3)]),directed=FALSE))#Emily: this reorders the edges and the list of species, so it is not the same as temp. thus below the merge misses some species because the name is first on one list and second on the other. the only thing this does is remove loops and multiple edges, but there should be any of these things in the data set to begin with, so I'm skipping this.
        #edge_list<-data.frame(get.edgelist(temp.graph,names=TRUE))
        #edge_list$pairs<-paste(edge_list$X1,edge_list$X2)
        #edge_list_pvals<-merge(edge_list,dataset_trt_no0,by="pairs",sort=FALSE  )
        temp$rho_cut<-rhos[r]
        #edge_list_pvals$trt<-trts[t]
        #temp$qval1<-fdrtool(temp$p.value,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
        temp$qval2<-p.adjust(temp$p.value, method="fdr")#using a different method for p-value adjustment. the fdrtool wasn't working when I had a lot of p values that were significant and none around 1
        temp<-cbind(pairs=temp$pairs,temp[,1:7],temp[,9:10])
        final.results<-rbind(final.results,temp)	}
    }
    print(t/length(trts))
  }
  return(final.results)
}

edge_lists<-co_occur_pairs(results)
dim(subset(edge_lists,rho>.5))
head(edge_lists)


######Edge creation II - with a p cutoff rather than a rho cutoff######
#
co_occur_pairs_all<-function(dataset){
  final.results<-data.frame()
  trts<-as.vector(unique(dataset$trt))
  for(t in 1:length(trts)){
    dataset_trt<-subset(dataset, trt==trts[t])
    dataset_trt_no0<-subset(dataset_trt, ab1 > 0 & ab2 > 0)#Ryan had this at 0
    dataset_trt_no0$pairs<-paste(dataset_trt_no0$taxa1,dataset_trt_no0$taxa2)
    temp<-dataset_trt_no0
    if(dim(temp)[1]>1){
      temp$qval1<-fdrtool(temp$spearmanp.value,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
      temp$qval2<-p.adjust(temp$spearmanp.value, method="fdr")#using a different method for p-value adjustment. the fdrtool wasn't working when I had a lot of p values that were significant and none around 1
      temp<-cbind(pairs=temp$pairs,temp[,1:11],temp[,13:14])
      #temp<-cbind(pairs=temp$pairs,temp[,1:7],temp[,9:10])
      final.results<-rbind(final.results,temp)	}
    print(t/length(trts))
  }
  return(final.results)
}

edge_listsno4orlower<-co_occur_pairs_all(resultsno4orlower)#spearman, p values are bad
edge_listsno4orlowerK<-co_occur_pairs_all(resultsno4orlowerK)#kelmann, pvalues good
edge_listsno4abunorlowerK<-co_occur_pairs_all(resultsno4abunorlowerK)
edge_listsno4orlowerK$qval3<-p.adjust(edge_listsno4orlowerK$p.value,method="fdr")

edge_listsKSno4<-cbind(pairs=paste(resultsKSno4$taxa1,resultsKSno4$taxa2),resultsKSno4)
edge_listsKSno4$qval<-p.adjust(edge_listsKSno4$spearmanp.value,method="fdr")
edge_listsKSno4b<-subset(edge_listsKSno4,spearmanrho>0)

edge_listsKSno5<-cbind(pairs=paste(resultsKSno5$taxa1,resultsKSno5$taxa2),resultsKSno5)
edge_listsKSno5$qval<-p.adjust(edge_listsKSno5$spearmanp.value,method="fdr")
edge_listsKSno5b<-subset(edge_listsKSno5,spearmanrho>0)

edge_listsKSno3<-cbind(pairs=paste(resultsKSno3$taxa1,resultsKSno3$taxa2),resultsKSno3)
edge_listsKSno3$qval<-p.adjust(edge_listsKSno3$spearmanp.value,method="fdr")
edge_listsKSno3b<-subset(edge_listsKSno3,spearmanrho>0)

edge_listsKS32no2<-cbind(pairs=paste(resultsKS32no2$taxa1,resultsKS32no2$taxa2),resultsKS32no2)
edge_listsKS32no2$qval<-p.adjust(edge_listsKS32no2$spearmanp.value,method="fdr")
edge_listsKS32no2b<-subset(edge_listsKS32no2,spearmanrho>0)

edge_listsKS32p<-cbind(pairs=paste(resultsKS32p$taxa1,resultsKS32p$taxa2),resultsKS32p)#the frequency cutoff was >2 *I think*, this was done within the loop so not here
edge_listsKS32p$qval<-p.adjust(edge_listsKS32p$spearmanp.value,method="fdr")
edge_listsKS32pb<-subset(edge_listsKS32p,spearmanrho>0)
head(edge_listsKS32pb)


edge_listsBEP<-cbind(pairs=paste(resultsBEP$taxa1,resultsBEP$taxa2),resultsBEP)
edge_listsBEP$qval<-p.adjust(edge_listsBEP$spearmanp.value,method="fdr")

edge_listsBEPb<-subset(edge_listsBEP,ab1freq>10&ab2freq>10)
#edge_listsBEPb<-subset(edge_listsBEP,ab1freq>5&ab2freq>5)
edge_listsBEPb$qval<-p.adjust(edge_listsBEPb$spearmanp.value,method="fdr")
edge_listsBEPc<-subset(edge_listsBEPb,spearmanrho>0)
dim(edge_listsBEPc)
min(edge_listsBEPc$ab2freq)
dim(subset(edge_listsBEPc,trt=="hi"))
#dim(comm16S.spelo2)[2]+dim(commEuk.spelo2)[2]


