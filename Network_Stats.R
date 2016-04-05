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




#Ror low medium and hi
#Kendall is too conservative, if I do fdr correction, only 2, 12,and 4 vertices remain in my networks. Spearman is less conservative and actually seemed closer to the permutation values I was getting (however this needs to be investigated more if I want to argue it), many vertices remain after fdr correction 
#I think the most correct way of doing this is: do the qvalue correction on all spearman p values (positive and negative rhos), then subset and use only positive rho
#edge_listsKSno4b<-subset(edge_listsKSno4,spearmanrho>0)
#edge_listsKSno4b$qval<-p.adjust(edge_listsKSno4b$spearmanp.value,method="fdr")
#edge_listsKSno4b<-subset(edge_listsKSno4qbygroup,spearmanrho>0)

edge_listsKSno4b,edge_listsKS32no2b edge_listsKS32pb edge_listsBEPc# fdr correction on all spearman pvalues, and subset to only positive spearman rhos 

inputlov<-subset(edge_listsBEPc,qval<.005&spearmanrho>.7&trt=="lo")
inputmev<-subset(edge_listsBEPc,qval<.005&spearmanrho>.7&trt=="me")
inputhiv<-subset(edge_listsBEPc,qval<.005&spearmanrho>.7&trt=="hi")
#inputhiv<-subset(edge_listsKS32no2b,qval<.05&trt=="hi")
#inputmev<-subset(edge_listsKS32no2b,qval<.05&trt=="me")
#inputlov<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")

graphlo<-simplify(graph.edgelist(as.matrix(inputlov[,3:4]),directed=FALSE))
graphme<-simplify(graph.edgelist(as.matrix(inputmev[,3:4]),directed=FALSE))
graphhi<-simplify(graph.edgelist(as.matrix(inputhiv[,3:4]),directed=FALSE))

length(V(graphlo))#500
length(V(graphme))#667
length(V(graphhi))#689
length(E(graphlo))#1022
length(E(graphme))#1549
length(E(graphhi))#1280
graph.density(graphlo)
graph.density(graphme)
graph.density(graphhi)
transitivity(graphlo, type="global")
transitivity(graphme, type="global")
transitivity(graphhi, type="global")#clustering, triangle
modularity(graphlo, membership(walktrap.community(graphlo)))
modularity(graphme, membership(walktrap.community(graphme)))
modularity(graphhi, membership(walktrap.community(graphhi)))

#get random transitivity values for hi me lo, this is just one random draw, could do many
transitivity(erdos.renyi.game(length(V(graphlo)),length(E(graphlo)),type="gnm"))
transitivity(erdos.renyi.game(length(V(graphme)),length(E(graphme)),type="gnm"))
transitivity(erdos.renyi.game(length(V(graphhi)),length(E(graphhi)),type="gnm"))

statslo<-data.frame(otu=row.names((as.matrix(degree(graphlo,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphlo,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphlo))))
statsme<-data.frame(otu=row.names((as.matrix(degree(graphme,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphme,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphme))))
statshi<-data.frame(otu=row.names((as.matrix(degree(graphhi,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphhi,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphhi))))

statslo1<-statslo[order(statslo$betweenness,decreasing=T),][1:5,]
statsme1<-statsme[order(statsme$betweenness,decreasing=T),][1:5,]
statshi1<-statshi[order(statshi$betweenness,decreasing=T),][1:5,]

labelsall[labelsall$otu%in%rownames(statslo1),]
labelsall[labelsall$otu%in%rownames(statsme1),]
labelsall[labelsall$otu%in%rownames(statshi1),]

labelsall[labelsall$otu=="denovo143776",]
tax_table(dats2)[rownames(tax_table(dats2))=="denovo143776",]#euk
tax_table(dat16Ss2)[rownames(tax_table(dat16Ss2))=="denovo462276",]



# Number of plants bacteria and euks in networks
verticesgraphlo<-as.data.frame(rownames(as.matrix(V(graphlo))))
verticesgraphlo2<-labelsall[which(labelsall$otuxy%in%verticesgraphlo[,1]),]
length(which(verticesgraphlo2$domain=="plant"))
length(which(verticesgraphlo2$domain=="bac"))
length(which(verticesgraphlo2$domain=="euk"))

verticesgraphme<-as.data.frame(rownames(as.matrix(V(graphme))))
verticesgraphme2<-labelsall[which(labelsall$otuxy%in%verticesgraphme[,1]),]
length(which(verticesgraphme2$domain=="plant"))
length(which(verticesgraphme2$domain=="bac"))
length(which(verticesgraphme2$domain=="euk"))

verticesgraphhi<-as.data.frame(rownames(as.matrix(V(graphhi))))
verticesgraphhi2<-labelsall[which(labelsall$otuxy%in%verticesgraphhi[,1]),]
length(which(verticesgraphhi2$domain=="plant"))
length(which(verticesgraphhi2$domain=="bac"))
length(which(verticesgraphhi2$domain=="euk"))


#intersection of edges of two graphs
graph.intersection(graphlo,graphme)#union of the vertex names and an intersection of the edges
graph.union(graphlo,graphme)

graph.intersection(graphlo,graphhi)
graph.union(graphlo,graphhi)

graph.intersection(graphme,graphhi)
graph.union(graphme,graphhi)


#intersection of vertices
length(intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo)))))
length(union(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo)))))

length(intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphhi)))))
length(union(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphhi)))))

length(intersect(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo)))))
length(union(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo)))))


tax_table(dats2)[rownames(tax_table(dats2))=="denovo215842",]
tax_table(dats2)[rownames(tax_table(dats2))=="denovo12769",]


