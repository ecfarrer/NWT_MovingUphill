#NetworkStats

library(NetIndices)

inputhiv<-subset(edge_listsno4orlower,trt=="hi")#rho_cut>0.4&rho>0
inputlov<-subset(edge_listsno4orlower,trt=="lo")#rho_cut>0.4&rho>0
inputnov<-subset(edge_listsno4orlower,trt=="no")#rho_cut>0.4&rho>0

inputhiv<-subset(edge_listsno4orlowerK,qval3<.1&trt=="hi")
inputlov<-subset(edge_listsno4orlowerK,qval3<.1&trt=="lo")
inputnov<-subset(edge_listsno4orlowerK,qval3<.1&trt=="no")

inputhiv<-subset(edge_listsno4orlowerK,p.value<.05&rho>=.45&trt=="hi")
inputlov<-subset(edge_listsno4orlowerK,p.value<.05&rho>=.45&trt=="lo")
inputnov<-subset(edge_listsno4orlowerK,p.value<.05&rho>=.45&trt=="no")

dim(inputhiv)
dim(inputlov)
dim(inputnov)

inputhiv<-subset(edge_listsno4orlower,rho>=.5&trt=="hi")
inputlov<-subset(edge_listsno4orlower,rho>=.5&trt=="lo")
inputnov<-subset(edge_listsno4orlower,rho>=.5&trt=="no")

graphhi<-simplify(graph.edgelist(as.matrix(inputhiv[,3:4]),directed=FALSE))
graphlo<-simplify(graph.edgelist(as.matrix(inputlov[,3:4]),directed=FALSE))
graphno<-simplify(graph.edgelist(as.matrix(inputnov[,3:4]),directed=FALSE))

length(V(graphhi))
length(V(graphlo))
length(V(graphno))
length(E(graphhi))
length(E(graphlo))
length(E(graphno))
graph.density(graphhi)
graph.density(graphlo)
graph.density(graphno)
transitivity(graphhi, type="global")#clustering, triangle
transitivity(graphlo, type="global")
transitivity(graphno, type="global")
modularity(graphhi, membership(walktrap.community(graphhi)))
modularity(graphlo, membership(walktrap.community(graphlo)))
modularity(graphno, membership(walktrap.community(graphno)))
#edge.connectivity(graphhi)
#edge.connectivity(graphlo)
#edge.connectivity(graphno)
#modularity(graphhi, membership(edge.betweenness.community(graphhi)))
#modularity(graphlo, membership(edge.betweenness.community(graphlo)))
#modularity(graphno, membership(edge.betweenness.community(graphno)))
diameter(graphhi)#length of longest geodesic
diameter(graphlo)
diameter(graphno)
clusters(graphhi)$no
clusters(graphlo)$no
clusters(graphno)$no


graphhi.adj<-get.adjacency(graphhi,sparse=F)
graphhi.properties<-GenInd(graphhi.adj)
graphlo.adj<-get.adjacency(graphlo,sparse=F)
graphlo.properties<-GenInd(graphlo.adj)
graphno.adj<-get.adjacency(graphno,sparse=F)
graphno.properties<-GenInd(graphno.adj)

graphhi.properties$C#connectance=proportion of possible links between species that are realized, this is the same as graph.density() above
graphlo.properties$C
graphno.properties$C


hist(inputhiv$p.value,freq=F,xlim=c(0,.1),breaks=1000)
hist(inputlov$p.value,freq=F,xlim=c(0,.1),breaks=1000)
hist(inputnov$p.value,freq=F,xlim=c(0,.1),breaks=1000)




#for low medium and hi
#Kendall is too conservative, if I do fdr correction, only 2, 12,and 4 vertices remain in my networks. Spearman is less conservative and actually seemed closer to the permutation values I was getting (however this needs to be investigated more if I want to argue it), many vertices remain after fdr correction 
#I think the most correct way of doing this is: do the qvalue correction on all spearman p values (positive and negative rhos), then subset and use only positive rho
#edge_listsKSno4b<-subset(edge_listsKSno4,spearmanrho>0)
#edge_listsKSno4b$qval<-p.adjust(edge_listsKSno4b$spearmanp.value,method="fdr")
#edge_listsKSno4b<-subset(edge_listsKSno4qbygroup,spearmanrho>0)

edge_listsKSno4b,edge_listsKS32no2b edge_listsKS32pb# fdr correction on all spearman pvalues, and subset to only positive spearman rhos

inputhiv<-subset(edge_listsKS32pb,qval<.05&spearmanrho>=.65&trt=="hi")
inputmev<-subset(edge_listsKS32pb,qval<.05&spearmanrho>=.65&trt=="me")
inputlov<-subset(edge_listsKS32pb,qval<.05&spearmanrho>=.65&trt=="lo")
#inputhiv<-subset(edge_listsKS32no2b,qval<.05&trt=="hi")
#inputmev<-subset(edge_listsKS32no2b,qval<.05&trt=="me")
#inputlov<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")

graphhi<-simplify(graph.edgelist(as.matrix(inputhiv[,3:4]),directed=FALSE))
graphme<-simplify(graph.edgelist(as.matrix(inputmev[,3:4]),directed=FALSE))
graphlo<-simplify(graph.edgelist(as.matrix(inputlov[,3:4]),directed=FALSE))

length(V(graphhi))#689
length(V(graphme))#667
length(V(graphlo))#500
length(E(graphhi))#1280
length(E(graphme))#1549
length(E(graphlo))#1022
graph.density(graphhi)
graph.density(graphme)
graph.density(graphlo)
transitivity(graphhi, type="global")#clustering, triangle
transitivity(graphme, type="global")
transitivity(graphlo, type="global")
modularity(graphhi, membership(walktrap.community(graphhi)))
modularity(graphme, membership(walktrap.community(graphme)))
modularity(graphlo, membership(walktrap.community(graphlo)))

#get random transitivity values for hi me lo, this is just one random draw, could do many
transitivity(erdos.renyi.game(length(V(graphhi)),length(E(graphhi)),type="gnm"))
transitivity(erdos.renyi.game(length(V(graphme)),length(E(graphme)),type="gnm"))
transitivity(erdos.renyi.game(length(V(graphlo)),length(E(graphlo)),type="gnm"))

statshi<-data.frame(otu=row.names((as.matrix(degree(graphhi,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphhi,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphhi))))
statsme<-data.frame(otu=row.names((as.matrix(degree(graphme,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphme,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphme))))
statslo<-data.frame(otu=row.names((as.matrix(degree(graphlo,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphlo,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphlo))))

statshi[order(statshi$betweenness,decreasing=T),][1:20,]
statsme[order(statsme$betweenness,decreasing=T),][1:20,]
statslo[order(statslo$betweenness,decreasing=T),][1:20,]



#intersection of edges of two graphs
graph.intersection(graphlo,graphhi)#union of the vertex names and an intersection of the edges
graph.intersection(graphlo,graphme)
graph.intersection(graphme,graphhi)
graph.union(graphlo,graphhi)

intersect(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo))))
intersect(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphme))))
intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo))))

union(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo))))

intersection(V(graphhi),V(graphlo))
union(V(graphhi),V(graphlo))

tax_table(dats2)[rownames(tax_table(dats2))=="denovo215842",]
tax_table(dats2)[rownames(tax_table(dats2))=="denovo12769",]


