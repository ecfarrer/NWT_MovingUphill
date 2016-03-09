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
        
        #temp$qval<-fdrtool(temp$p.value,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
        temp$qval<-p.adjust(temp$p.value, method="fdr")#using a different method for p-value adjustment. the fdrtool wasn't working when I had a lot of p values that were significant and none around 1
        
        temp<-cbind(pairs=temp$pairs,temp[,1:7],temp[,9:10])
        final.results<-rbind(final.results,temp)	}
    }
    print(t/length(trts))
  }
  return(final.results)
}

#results<-read.csv(file.choose())
edge_lists<-co_occur_pairs(results)




######Edge creation II - with a p cutoff rather than a rho cutoff######
co_occur_pairs_all<-function(dataset){
  final.results<-data.frame()
  trts<-as.vector(unique(dataset$trt))
  for(t in 1:length(trts)){
    dataset_trt<-subset(dataset, trt==trts[t])
    dataset_trt_no0<-subset(dataset_trt, ab1 > 2 & ab2 > 2)#Ryan had this at 0
    dataset_trt_no0$pairs<-paste(dataset_trt_no0$taxa1,dataset_trt_no0$taxa2)
    temp<-dataset_trt_no0
    if(dim(temp)[1]>1){
      temp$qval1<-fdrtool(temp$p.value,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
      temp$qval2<-p.adjust(temp$p.value, method="fdr")#using a different method for p-value adjustment. the fdrtool wasn't working when I had a lot of p values that were significant and none around 1
      temp<-cbind(pairs=temp$pairs,temp[,1:7],temp[,9:10])
      final.results<-rbind(final.results,temp)	}
    print(t/length(trts))
  }
  return(final.results)
}

#results<-read.csv(file.choose())
edge_lists2<-co_occur_pairs_all(results)
edge_lists2s<-subset(edge_lists2,qval2<0.05&rho>0)