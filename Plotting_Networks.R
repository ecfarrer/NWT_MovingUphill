#Plotting

plantcols<-data.frame(kingdomlabels=c("Amoebozoa","Archaeplastida","Discicristoidea","Fungi","Holozoa","Nonphotosynthetic_Alveolata","Nonphotosynthetic_Discoba","Nonphotosynthetic_Eukaryota","Photosynthetic_Alveolata","Photosynthetic_Discoba","Rhizaria","Plant"),color=c("red","blue","red","yellow","orange","red","red","red","blue","blue","red","green"))
#labelsall<-rbind(labelfile,plantlabels)
labelsall2<-merge(labelfile,plantcols)
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



