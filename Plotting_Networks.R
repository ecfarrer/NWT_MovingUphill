#Plotting
#using edge_listsKSno4b or edge_listsKS32no2b - has fdr correction on all spearman pvalues, and subset to only positive spearman rhos

plantcols<-data.frame(kingdomlabels=c("Amoebozoa","Archaeplastida","Discicristoidea","Fungi","Holozoa","Nonphotosynthetic_Alveolata","Nonphotosynthetic_Discoba","Nonphotosynthetic_Eukaryota","Photosynthetic_Alveolata","Photosynthetic_Discoba","Rhizaria","Plant"),color=c("red","blue","red","yellow","orange","red","red","red","blue","blue","red","green"))
labelsall1<-rbind(labelfile,plantlabels)
labelsall<-unique(merge(labelsall1,plantcols,"kingdomlabels"))
labelsall$color<-as.character(labelsall$color)


#high density
inputhi<-subset(edge_listsKS32pb,qval<.05&spearmanrho>.65&trt=="hi")[,3:4]#
inputhi<-subset(edge_listsKSno5b,spearmanrho>=.5&trt=="hi")[,3:4]#
inputhi<-subset(edge_listsKS32no2b,qval<.05&spearmanrho>.6&trt=="hi")[,3:4]#could use .65 if I want to simplify it
inputhiv<-subset(edge_listsKS32no2b,qval<.05&trt=="hi")#
vertexsizes1<-unique(data.frame(otu=c(as.character(inputhiv$taxa1),as.character(inputhiv$taxa2)),abun=c(inputhiv$ab1,inputhiv$ab2)))

graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph1)
#plot(graph1,vertex.size=2,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu" #change to "order" if doing things by order
colorgraph1<-unique(merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F))
sizesgraph1<-merge(verticesgraph1,vertexsizes1,"otu",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantf3q.05r.65.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2
#dev.off()



taxagraph1<-as.data.frame(get.edgelist(graph1))
colnames(taxagraph1)[1]<-"orders"
taxagraph1a<-unique(merge(taxagraph1,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph1a)<-c("V1","orders","V1kingdoms")
taxagraph1b<-unique(merge(taxagraph1a,labelfile,"orders",sort=F))
colnames(taxagraph1b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph1c<-data.frame(V1=taxagraph1b$V1,V2=taxagraph1b$V2,taxagraph1b[,3:4])


#medium density
inputme<-subset(edge_listsKS32pb,qval<.05&spearmanrho>.65&trt=="me")[,3:4]
inputme<-subset(edge_listsKSno5b,spearmanrho>=.5&trt=="me")[,3:4]
inputme<-subset(edge_listsKS32no2b,qval<.05&trt=="me")[,3:4]
inputmev<-subset(edge_listsKS32no2b,qval<.05&trt=="me")
vertexsizes2<-unique(data.frame(otu=c(as.character(inputmev$taxa1),as.character(inputmev$taxa2)),abun=c(inputmev$ab1,inputmev$ab2)))

graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-unique(merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F))
sizesgraph2<-merge(verticesgraph2,vertexsizes2,"otu",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantf3q.05r.65.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2
#dev.off()

taxagraph2<-as.data.frame(get.edgelist(graph2))
colnames(taxagraph2)[1]<-"orders"
taxagraph2a<-unique(merge(taxagraph2,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph2a)<-c("V1","orders","V1kingdoms")
taxagraph2b<-unique(merge(taxagraph2a,labelfile,"orders",sort=F))
colnames(taxagraph2b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph2c<-data.frame(V1=taxagraph2b$V1,V2=taxagraph2b$V2,taxagraph2b[,3:4])


#low density
inputlo<-subset(edge_listsKS32pb,qval<.05&spearmanrho>.65&trt=="lo")[,3:4]
inputlo<-subset(edge_listsKSno5b,spearmanrho>=.5&trt=="lo")[,3:4]
inputlo<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")[,3:4]
inputlov<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")
vertexsizes3<-unique(data.frame(otu=c(as.character(inputlov$taxa1),as.character(inputlov$taxa2)),abun=c(inputlov$ab1,inputlov$ab2)))

graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-unique(merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F))
sizesgraph3<-merge(verticesgraph3,vertexsizes3,"otu",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantf3q.05r.65.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()






temp<-subset(comm.data,greater66plants=="hi")
plot(jitter(temp$Classiculales),temp[,"Unclassified Coccolithales"])
summary(lm(temp$Classiculales~temp[,"Unclassified Coccolithales"]))
cor.test(temp$Classiculales,temp[,"Unclassified Coccolithales"],method="spearman",na.action=na.rm)





old old


















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






#Plotting
inputhi<-subset(edge_lists,trt=="hi")[,3:4]
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
eb<-edge.betweenness.community(graph1)
membership(eb)
plot(graph1,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu"
colorgraph1<-unique(merge(verticesgraph1,labelsall2,"otu",all.y=F,all.x=F,sort=F))
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotu.pdf")
plot(graph1,vertex.size=2,vertex.color=colorgraph1$color,edge.curved=T,vertex.label=NA)
dev.off()

inputlo<-subset(edge_lists,trt=="lo")[,3:4]
graph2<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
eb<-edge.betweenness.community(graph2)
membership(eb)
plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-unique(merge(verticesgraph2,labelsall2,"otu",all.y=F,all.x=F,sort=F))
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotu.pdf")
plot(graph2,vertex.size=2,vertex.color=colorgraph2$color,edge.curved=T,vertex.label=NA)
dev.off()



