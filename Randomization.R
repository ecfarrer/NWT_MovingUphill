#Random networks


#####equalizing euks and bacteria#####

#draw 1000 pairs from euk and bac and make a network

head(edge_listsBEPb)
temp<-unique(c(as.character(edge_listsBEPb$taxa1),as.character(edge_listsBEPb$taxa2)))

tempe<-intersect(temp,labelsall$otuxy[labelsall$domain=="euk"])
length(tempe)
tempe1<-sample(tempe,200)

tempb<-intersect(temp,labelsall$otuxy[labelsall$domain=="bac"])
length(tempb)
tempb1<-sample(tempb,200)

ind<-which(edge_listsBEPb$taxa1%in%c(tempe1,tempb1)&edge_listsBEPb$taxa2%in%c(tempe1,tempb1))
temp2<-edge_listsBEPb[ind,]
dim(temp2)
head(temp2)
#edge_listsBEPb$qval<-p.adjust(edge_listsBEPb$spearmanp.value,method="fdr")
temp3<-subset(temp2,spearmanrho>0)
dim(temp3)

#Low density
inputlo<-subset(temp3,qval<.05&spearmanrho>.5&trt=="lo")[,3:4]#
dim(inputlo)
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otuxy"
colorgraph3<-merge(verticesgraph3,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbaceukf10q.005r.7.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()


#Medium density
inputme<-subset(temp3,qval<.05&spearmanrho>.5&trt=="me")[,3:4]
dim(inputme)
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otuxy"
colorgraph2<-merge(verticesgraph2,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbaceukf10q.005r.7.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2
#dev.off()


#High density
inputhi<-subset(temp3,qval<.05&spearmanrho>.5&trt=="hi")[,3:4]
dim(inputhi)
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otuxy" #change to "order" if doing things by order
colorgraph1<-merge(verticesgraph1,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbaceukf10q.005r.7.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2
#dev.off()

length(which(colorgraph3$domain=="euk"))
length(which(colorgraph3$domain=="bac"))
length(which(colorgraph2$domain=="euk"))
length(which(colorgraph2$domain=="bac"))
length(which(colorgraph1$domain=="euk"))
length(which(colorgraph1$domain=="bac"))


#conclusions: the amount of bacteria and euks in the network is more equal but still slightly skewed toward bacteria, even when starting with equal number of taxa input. there are roughtly 10% to 40% more bacteria in the network.







