setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")


#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata") #results has the results from OTU cooccurrence patterns

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")


library(vegan)
library(reshape)
library(plotrix)
library(foreach)
library(doParallel)



######Co-occurrence pairwise correlations###### 
#Input dataset is dats6order

comm.data<-dats6order
comm.data<-dats9otu
  
#Standardize to 100 or rarified to 1250 (order) or 1191 (otu). if I standardize, I will need to change code in the loop below because it has cutoffs of 1
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1250))
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1191))
#comm.data[,27:243]<-comm.data[,27:243]/rowSums(comm.data[,27:243])*100

#temp<-cbind(as.character(comm.data$Sample_name),comm.data$Plant_Dens)
#temp[order(as.numeric(as.character(comm.data$Sample_name))),]

greater66plants<-factor(ifelse(comm.data$Plant_Dens>66,"hi","lo")) #this is stem density, including mosses
comm.data<-cbind(greater66plants,comm.data)
trts<-as.vector(unique(greater66plants))

results<-matrix(nrow=0,ncol=7)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.data, greater66plants==trt.temp)
  
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
resultsold<-results
resultsoldotu<-results
#results<-data.frame(data.matrix(results))
#names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")

#Emily: it made the numeric values in results factors so I had to change it
rownames(resultsold)<-1:dim(resultsold)[1]
results2<-data.frame(resultsold[,1:3],as.numeric(resultsold[,4]),as.numeric(resultsold[,5]),as.numeric(resultsold[,6]),as.numeric(resultsold[,7]))
results<-results2
names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")

rownames(resultsoldotu)<-1:dim(resultsoldotu)[1]
results2<-data.frame(resultsoldotu[,1:3],as.numeric(resultsoldotu[,4]),as.numeric(resultsoldotu[,5]),as.numeric(resultsoldotu[,6]),as.numeric(resultsoldotu[,7]))
results<-results2
names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")









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
edge_lists<-co_occur_pairs(results)



#Plotting

plantcols<-data.frame(kingdomlabels=c("Amoebozoa","Archaeplastida","Discicristoidea","Fungi","Holozoa","Nonphotosynthetic_Alveolata","Nonphotosynthetic_Discoba","Nonphotosynthetic_Eukaryota","Photosynthetic_Alveolata","Photosynthetic_Discoba","Rhizaria"),color=c("red","blue","red","yellow","orange","red","red","red","blue","blue","red"))
labelsall3<-unique(merge(labelfile,plantcols))
labelsall3$color<-as.character(labelsall3$color)

#high density
inputhi<-subset(edge_lists, rho_cut>0.4&trt=="hi")[,3:4]
inputhiv<-subset(edge_lists, rho_cut>0.4&trt=="hi")
vertexsizes<-unique(data.frame(orders=c(as.character(inputhiv$taxa1),as.character(inputhiv$taxa2)),abun=c(inputhiv$ab1,inputhiv$ab2)))

graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph1)
#membership(eb)
plot(graph1,vertex.size=2,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"orders"
colorgraph1<-unique(merge(verticesgraph1,labelsall3,"orders",all.y=F,all.x=F,sort=F))
sizesgraph1<-merge(verticesgraph1,vertexsizes,"orders",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityorder.pdf")
plot(graph1,vertex.size=log(sizesgraph1$abun)*2,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,vertex.label=NA)#
#dev.off()

taxagraph1<-as.data.frame(get.edgelist(graph1))
colnames(taxagraph1)[1]<-"orders"
taxagraph1a<-unique(merge(taxagraph1,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph1a)<-c("V1","orders","V1kingdoms")
taxagraph1b<-unique(merge(taxagraph1a,labelfile,"orders",sort=F))
colnames(taxagraph1b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph1c<-data.frame(V1=taxagraph1b$V1,V2=taxagraph1b$V2,taxagraph1b[,3:4])


#low density
inputlo<-subset(edge_lists, rho_cut>0.4&trt=="lo")[,3:4]
inputlov<-subset(edge_lists, rho_cut>0.4&trt=="lo")
vertexsizes<-unique(data.frame(orders=c(as.character(inputlov$taxa1),as.character(inputlov$taxa2)),abun=c(inputlov$ab1,inputlov$ab2)))

graph2<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"orders"
colorgraph2<-unique(merge(verticesgraph2, labelsall3,"orders",all.y=F,all.x=F,sort=F))
sizesgraph2<-merge(verticesgraph2,vertexsizes,"orders",sort=F)
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityorder.pdf")
plot(graph2,vertex.size=log(sizesgraph2$abun)*2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,vertex.label=NA)#
dev.off()

taxagraph2<-as.data.frame(get.edgelist(graph2))
colnames(taxagraph2)[1]<-"orders"
taxagraph2a<-unique(merge(taxagraph2,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph2a)<-c("V1","orders","V1kingdoms")
taxagraph2b<-unique(merge(taxagraph2a,labelfile,"orders",sort=F))
colnames(taxagraph2b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph2c<-data.frame(V1=taxagraph2b$V1,V2=taxagraph2b$V2,taxagraph2b[,3:4])





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





