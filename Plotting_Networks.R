#Plotting
#using edge_listsKSno4b or edge_listsKS32no2b - has fdr correction on all spearman pvalues, and subset to only positive spearman rhos



#Fix label files due to duplicates in euk and bac datasets for OTU names
#the duplicates for euks and bac "denovo60343"  "denovo147101" "denovo169780" "denovo188318" "denovo207652"
#not there "denovo169780.x" "denovo188318.x"
#euks is x, bacteria y
which(resultsBEP$taxa1=="denovo188318.y")

#changing the label files
head(bacterialabels)
head(labelfile)#euks

bacterialabels$otuxy<-as.character(bacterialabels$otu)
bacterialabels[which(bacterialabels$otuxy%in%c("denovo60343","denovo147101","denovo169780","denovo188318","denovo207652")),"otuxy"]<-c("denovo60343.y","denovo147101.y","denovo169780.y","denovo188318.y","denovo207652.y")
bacterialabels[8000:8030,]

labelfile$otuxy<-as.character(labelfile$otu)
labelfile[which(labelfile$otuxy%in%c("denovo60343","denovo147101","denovo169780","denovo188318","denovo207652")),"otuxy"]<-c("denovo60343.x","denovo147101.x","denovo169780.x","denovo188318.x","denovo207652.x")
labelfile[9700:9770,]

plantlabels$otuxy<-plantlabels$otu

#get sequence names from the comm.data files, match with euk or 16S label file, then rbind labelfiles
names16S<-data.frame(colnames(comm.data16S)[27:dim(comm.data16S)[2]])
names(names16S)<-"otu"
names16Sb<-merge(names16S,bacterialabels,"otu",all.y=F)
head(names16Sb)

namesEuk<-data.frame(colnames(comm.dataEuk)[27:dim(comm.dataEuk)[2]])
names(namesEuk)<-"otu"
namesEukb<-merge(namesEuk,labelfile,"otu",all.y=F)


plantcols<-data.frame(kingdomlabels=c("Amoebozoa","Archaeplastida","Discicristoidea","Fungi","Holozoa","Nonphotosynthetic_Alveolata","Nonphotosynthetic_Discoba","Nonphotosynthetic_Eukaryota","Photosynthetic_Alveolata","Photosynthetic_Discoba","Rhizaria","Plant","PhotosyntheticBacteria","NonphotosyntheticBacteria"),color=c("orangered","green","orangered","yellow","red","orangered","orangered","orangered","green","green","orangered","blue","darkgreen","#802b00"))
labelsall1<-rbind(cbind(namesEukb,domain=c("euk")),cbind(names16Sb,domain=c("bac")),cbind(plantlabels,domain=c("plant")))
labelsall<-merge(labelsall1,plantcols,"kingdomlabels")
labelsall$color<-as.character(labelsall$color)





###### Plotting ######

#Low density
inputlo<-subset(edge_listsKS32pb,qval<.05&spearmanrho>.65&trt=="lo")[,3:4]
inputlo<-subset(edge_listsKSno5b,spearmanrho>=.5&trt=="lo")[,3:4]
inputlo<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")[,3:4]
inputlo<-subset(edge_listsBEPc,qval<.005&spearmanrho>.7&trt=="lo")[,3:4]
dim(inputlo)
#inputlov<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")
#vertexsizes3<-unique(data.frame(otu=c(as.character(inputlov$taxa1),as.character(inputlov$taxa2)),abun=c(inputlov$ab1,inputlov$ab2)))

graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otuxy"
colorgraph3<-merge(verticesgraph3,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#sizesgraph3<-merge(verticesgraph3,vertexsizes3,"otu",sort=F)
sizesgraph3<-ifelse(verticesgraph3$otuxy=="denovo559741",5,2)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbaceukf10q.005r.7.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()
plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2



#Medium density
inputme<-subset(edge_listsKS32pb,qval<.05&spearmanrho>.65&trt=="me")[,3:4]
inputme<-subset(edge_listsKSno5b,spearmanrho>=.5&trt=="me")[,3:4]
inputme<-subset(edge_listsKS32no2b,qval<.05&trt=="me")[,3:4]
inputme<-subset(edge_listsBEPc,qval<.005&spearmanrho>.7&trt=="me")[,3:4]
dim(inputme)
#inputmev<-subset(edge_listsKS32no2b,qval<.05&trt=="me")
#vertexsizes2<-unique(data.frame(otu=c(as.character(inputmev$taxa1),as.character(inputmev$taxa2)),abun=c(inputmev$ab1,inputmev$ab2)))

graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otuxy"
colorgraph2<-merge(verticesgraph2,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#sizesgraph2<-merge(verticesgraph2,vertexsizes2,"otuxy",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbaceukf10q.005r.7.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2
#dev.off()



#High density
inputhi<-subset(edge_listsKS32pb,qval<.05&spearmanrho>.65&trt=="hi")[,3:4]#
inputhi<-subset(edge_listsKSno5b,spearmanrho>=.5&trt=="hi")[,3:4]#
inputhi<-subset(edge_listsKS32no2b,qval<.05&spearmanrho>.6&trt=="hi")[,3:4]#could use .65 if I want to simplify it
inputhi<-subset(edge_listsBEPc,qval<.005&spearmanrho>.7&trt=="hi")[,3:4]
dim(inputhi)
#inputhiv<-subset(edge_listsKS32no2b,qval<.05&trt=="hi")#
#vertexsizes1<-unique(data.frame(otu=c(as.character(inputhiv$taxa1),as.character(inputhiv$taxa2)),abun=c(inputhiv$ab1,inputhiv$ab2)))

graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph1)
#plot(graph1,vertex.size=2,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otuxy" #change to "order" if doing things by order
colorgraph1<-merge(verticesgraph1,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#sizesgraph1<-merge(verticesgraph1,vertexsizes1,"otuxy",sort=F)
sizesgraph1<-ifelse(verticesgraph1$otuxy=="denovo559741",6,2)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbaceukf10q.005r.7.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2
#dev.off()
plot(graph1,vertex.size=sizesgraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2





##if i want to only plot vertices with >1 edge
graph3b<-delete_vertices(graph3, which(degree(graph3)<2))
plot(graph3b,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2



temp<-subset(comm.data,greater66plants=="hi")
plot(jitter(temp$Classiculales),temp[,"Unclassified Coccolithales"])
summary(lm(temp$Classiculales~temp[,"Unclassified Coccolithales"]))
cor.test(temp$Classiculales,temp[,"Unclassified Coccolithales"],method="spearman",na.action=na.rm)


taxagraph1<-as.data.frame(get.edgelist(graph1))
colnames(taxagraph1)[1]<-"orders"
taxagraph1a<-unique(merge(taxagraph1,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph1a)<-c("V1","orders","V1kingdoms")
taxagraph1b<-unique(merge(taxagraph1a,labelfile,"orders",sort=F))
colnames(taxagraph1b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph1c<-data.frame(V1=taxagraph1b$V1,V2=taxagraph1b$V2,taxagraph1b[,3:4])

taxagraph2<-as.data.frame(get.edgelist(graph2))
colnames(taxagraph2)[1]<-"orders"
taxagraph2a<-unique(merge(taxagraph2,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph2a)<-c("V1","orders","V1kingdoms")
taxagraph2b<-unique(merge(taxagraph2a,labelfile,"orders",sort=F))
colnames(taxagraph2b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph2c<-data.frame(V1=taxagraph2b$V1,V2=taxagraph2b$V2,taxagraph2b[,3:4])





#extract only vertices that are connected to a plant in all plant densities to determine if plant correlations are still there
#figure out why there is more richness in high plant density but the same number of network vertices - is the added diversity all low abundance?


inputlo<-subset(edge_listsBEPc,qval<.05&spearmanrho>.6&trt=="lo")[,3:4]
dim(inputlo)
inputlo2<-inputlo[which(inputlo$taxa1%in%plantlabels$otu|inputlo$taxa2%in%plantlabels$otu),]
graph3<-simplify(graph.edgelist(as.matrix(inputlo2),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otuxy"
colorgraph3<-merge(verticesgraph3,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbaceukf10q.005r.7.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()


#Medium density
inputme<-subset(edge_listsBEPc,qval<.05&spearmanrho>.6&trt=="me")[,3:4]
dim(inputme)
inputme2<-inputme[which(inputme$taxa1%in%plantlabels$otu|inputme$taxa2%in%plantlabels$otu),]
graph2<-simplify(graph.edgelist(as.matrix(inputme2),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otuxy"
colorgraph2<-merge(verticesgraph2,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbaceukf10q.005r.7.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2
#dev.off()

#High density
inputhi<-subset(edge_listsBEPc,qval<.05&spearmanrho>.6&trt=="hi")[,3:4]
dim(inputhi)
inputhi2<-inputhi[which(inputhi$taxa1%in%plantlabels$otu|inputhi$taxa2%in%plantlabels$otu),]
graph1<-simplify(graph.edgelist(as.matrix(inputhi2),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otuxy" #change to "order" if doing things by order
colorgraph1<-merge(verticesgraph1,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbaceukf10q.005r.7.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2
#dev.off()





plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40")










#old old


















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



