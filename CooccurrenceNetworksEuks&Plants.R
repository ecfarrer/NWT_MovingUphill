

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")


library(vegan)
library(reshape)
library(plotrix)
library(foreach)
library(doParallel)

getDoParWorkers()


######Co-occurrence pairwise correlations###### 
#Input dataset is dats6order

comm.data<-dats6order#read.csv("/Users/farrer/Desktop/Networks/total_order_info.csv")
#comm.data<-comm.data[,-1]
#comm.data.read<-subset(comm.data, reads >= 1407)

#Standardize to 100 or rarified to 1250. if I standardize, I will need to change code in the loop below because it has cutoffs of 1
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1250))
#comm.data[,27:243]<-comm.data[,27:243]/rowSums(comm.data[,27:243])*100

greater66plants<-factor(ifelse(comm.data$Plant_Dens>66,"hi","lo")) #this is stem density, including mosses
comm.data<-cbind(greater66plants,comm.data)
trts<-as.vector(unique(greater66plants))


#Read in plant data
plantcomp<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots
dim(plantcomp)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantlabels<-as.data.frame(cbind(colnames(plantcomp2)[2:56],"Plant"))
colnames(plantlabels)<-c("orders","kingdomlabels")

#merge plants with microbes
comm.data$Sample_name<-as.numeric(as.character(comm.data$Sample_name))
comm.datamp<-merge(comm.data,plantcomp2,"Sample_name",sort=F)



#number of iterations
#setup parallel backend to use 8 processors
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
#loop
results<-matrix(nrow=0,ncol=7)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.datamp, greater66plants==trt.temp)
  
  #in this case the community data started at column 28, so the loop for co-occurrence has to start at that point
  for(b in 28:(dim(temp)[2]-1)){

	results1<-foreach(c=(b+1):(dim(temp)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      #if the column is all 0's no co-occurrence will be performed
      if(species1.ab >1 & species2.ab >1){
        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        rho<-test$estimate
        p.value<-test$p.value
      }    
      if(species1.ab <=1 | species2.ab <= 1){
        rho<-0
        p.value<-1
      }	
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],rho,p.value,species1.ab,species2.ab)
      new.row		
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsoldmp<-results
print(Sys.time()-strt)#1 minute
stopCluster(cl)

rownames(resultsoldmp)<-1:dim(resultsoldmp)[1]
results2<-data.frame(resultsoldmp[,1:3],as.numeric(resultsoldmp[,4]),as.numeric(resultsoldmp[,5]),as.numeric(resultsoldmp[,6]),as.numeric(resultsoldmp[,7]))
results<-results2
names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")

#start time
strt<-Sys.time()
#loop
results<-matrix(nrow=0,ncol=7)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.datamp, greater66plants==trt.temp)
  
  #in this case the community data started at column 28, so the loop for co-occurrence has to start at that point
  for(b in 28:(dim(temp)[2]-1)){
    #every species will be compared to every other species, so there has to be another loop that iterates down the rest of the columns
    for(c in (b+1):(dim(temp)[2])){
      
      #summing the abundances of species of the columns that will be compared
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      #if the column is all 0's no co-occurrence will be performed
      if(species1.ab >1 & species2.ab >1){
        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        rho<-test$estimate
        p.value<-test$p.value
      }
      
      if(species1.ab <=1 | species2.ab <= 1){
        rho<-0
        p.value<-1
      }	
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],rho,p.value,species1.ab,species2.ab)
      results<-rbind(results,new.row)			
    }
  }
  print(a/length(trts))
}
head(results)
resultsoldmp2<-results
print(Sys.time()-strt) #11min
#results<-data.frame(data.matrix(results))
#names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")

#Emily: it made the numeric values in results factors so I had to change it
rownames(resultsoldmp2)<-1:dim(resultsoldmp2)[1]
results2<-data.frame(resultsoldmp2[,1:3],as.numeric(resultsoldmp2[,4]),as.numeric(resultsoldmp2[,5]),as.numeric(resultsoldmp2[,6]),as.numeric(resultsoldmp2[,7]))
results2<-results2
names(results2)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")

#######probably should do one more test to make sure the results are the same
#results (from parallel), results2 (not parallel)







######Edge creation######

#input is "results" from above

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
edge_listsmp<-co_occur_pairs(results)



#Plotting

plantcols<-data.frame(kingdomlabels=c("Amoebozoa","Archaeplastida","Discicristoidea","Fungi","Holozoa","Nonphotosynthetic_Alveolata","Nonphotosynthetic_Discoba","Nonphotosynthetic_Eukaryota","Photosynthetic_Alveolata","Photosynthetic_Discoba","Rhizaria","Plant"),color=c("red","blue","red","yellow","orange","red","red","red","blue","blue","red","green"))
labelsall<-rbind(labelfile,plantlabels)
labelsall2<-merge(labelsall,plantcols)
labelsall2$color<-as.character(labelsall2$color)

inputhi<-subset(edge_listsmp, rho_cut>0.4&trt=="hi")[,3:4]
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
eb<-edge.betweenness.community(graph1)
membership(eb)
plot(graph1,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"orders"
colorgraph1<-unique(merge(verticesgraph1,labelsall2,"orders",all.y=F,all.x=F,sort=F))
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,edge.curved=T,vertex.label=NA)
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,edge.curved=T)

inputlo<-subset(edge_listsmp, rho_cut>0.4&trt=="lo")[,3:4]
graph2<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
eb<-edge.betweenness.community(graph2)
membership(eb)
plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"orders"
colorgraph2<-unique(merge(verticesgraph2,labelsall2,"orders",all.y=F,all.x=F,sort=F))
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,edge.curved=T,vertex.label=NA)
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,edge.curved=T)#yellow are fungi, light blue are green algae



temp<-subset(comm.data,greater66plants=="hi")
plot(jitter(temp$Classiculales),temp[,"Unclassified Coccolithales"])
summary(lm(temp$Classiculales~temp[,"Unclassified Coccolithales"]))
cor.test(temp$Classiculales,temp[,"Unclassified Coccolithales"],method="spearman",na.action=na.rm)






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

#Plotting
inputhi<-subset(edge_lists2s,trt=="hi")[,3:4]
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
eb<-edge.betweenness.community(graph1)
membership(eb)
plot(graph1,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"orders"
colorgraph1<-unique(merge(verticesgraph1,labelfile,"orders",all.y=F,all.x=F,sort=F))
plot(graph1,vertex.size=4,vertex.color=colorgraph1$kingdomlabels,edge.curved=T,vertex.label=NA)

inputlo<-subset(edge_lists2s,trt=="lo")[,3:4]
graph2<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
eb<-edge.betweenness.community(graph2)
membership(eb)
plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"orders"
colorgraph2<-unique(merge(verticesgraph2,labelfile,"orders",all.y=F,all.x=F,sort=F))
plot(graph2,vertex.size=4,vertex.color=colorgraph2$kingdomlabels,edge.curved=T,vertex.label=NA)
#yellow are fungi, light blue are green algae





