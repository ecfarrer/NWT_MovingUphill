load("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_AnalysisOrder.Rdata")

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_AnalysisOrder.Rdata")

set.seed(10)
comm.data<-dats10otu
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1121))#by otu

nohilo<-ifelse(comm.data$Plant_Dens==0,"no","else");nohilo[which(comm.data$Plant_Dens<43&comm.data$Plant_Dens>0)]<-"lo";nohilo[which(comm.data$Plant_Dens>=43)]<-"hi";nohilo<-factor(nohilo)#14 no,41 lo, 42 hi 

comm.data<-cbind(nohilo,comm.data)
trts<-as.vector(unique(nohilo))


comm.data[,28:1035]<-ifelse(comm.data[,28:1035]>0,1,0)

head(comm.data)[,2:30]
head(dats10otu)[,1:30]

temp1<-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
temp2<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
chisq.test(temp1,temp2,rescale.p=T)
cor.test(temp1,temp2)

M <- as.table(rbind(c(10,5), c(7,8)))
chisq.test(M,rescale.p=T,correct=F,simulate.p.value = F)  #corrected is better for small sample size, as sample size increases it approaches the uncorrected value

######Co-occurrence pairwise correlations###### 
#Input dataset is dats6order or dats10otu

comm.data<-dats6order
comm.data<-dats10otu
  
#Standardize to 100 or rarified to 1250 (order) or 1121 (otu). if I standardize, I will need to change code in the loop below because it has cutoffs of 1
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1250))#by order
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1121))#by otu
#comm.data[,27:243]<-comm.data[,27:243]/rowSums(comm.data[,27:243])*100

#temp<-cbind(as.character(comm.data$Sample_name),comm.data$Plant_Dens)
#temp[order(as.numeric(as.character(comm.data$Sample_name))),]

greater66plants<-factor(ifelse(comm.data$Plant_Dens>66,"hi","lo")) #this is stem density, including mosses

nohilo<-ifelse(comm.data$Plant_Dens==0,"no","else");nohilo[which(comm.data$Plant_Dens<81&comm.data$Plant_Dens>0)]<-"lo";nohilo[which(comm.data$Plant_Dens>=81)]<-"hi";nohilo<-factor(nohilo)#14 no,41 lo, 42 hi #THIS WAS MESSED UP NEED TO REDO

comm.data<-cbind(nohilo,comm.data)
trts<-as.vector(unique(nohilo))






#by order
comm.data<-dats6order

set.seed(10)
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1250))#by order

nohilo<-ifelse(comm.data$Plant_Dens==0,"no","else");nohilo[which(comm.data$Plant_Dens<81&comm.data$Plant_Dens>0)]<-"lo";nohilo[which(comm.data$Plant_Dens>=81)]<-"hi";nohilo<-factor(nohilo)#14 no,41 lo, 42 hi 

lomehi<-ifelse(comm.data$Plant_Dens<36,"lo","else");lomehi[which(comm.data$Plant_Dens<91&comm.data$Plant_Dens>=36)]<-"me";lomehi[which(comm.data$Plant_Dens>=91)]<-"hi";lomehi<-factor(lomehi)#32 no, 33 lo, 32 hi 

comm.data<-cbind(lomehi,comm.data)
trts<-as.vector(unique(lomehi))


#Setup parallel backend to use 4 processors, for all metrics
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
#loop
results<-matrix(nrow=0,ncol=21)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.data, lomehi==trt.temp)
  
  #in this case the community data started at column 28, so the loop for co-occurrence has to start at that point
  for(b in 28:(dim(temp)[2]-1)){
    
    results1<-foreach(c=(b+1):(dim(temp)[2]),.combine=rbind,.packages="Kendall") %dopar% {
      #for(c in 29:(dim(temp)[2])){      #for debugging
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(temp[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)
      #I changed this so that it will calculate the correlation for all species pairs, even singletons, and I can subset them later if I want to remove infrequent otus for qvalue calculation
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
        test<-Kendall(temp[,b],temp[,c])
        kendallrho<-test$tau
        kendallp.value<-test$sl
        
        if(length(unique(ifelse(temp[,b]>0,1,0)))>1&length(unique(ifelse(temp[,c]>0,1,0)))>1){
        test<-cor.test(ifelse(temp[,b]>0,1,0),ifelse(temp[,c]>0,1,0),method="pearson")
        phi<-test$estimate
        phip.value<-chisq.test(ifelse(temp[,b]>0,1,0),ifelse(temp[,c]>0,1,0),correct=F,simulate.p.value = F)$p.value
        phip.valuecor<-chisq.test(ifelse(temp[,b]>0,1,0),ifelse(temp[,c]>0,1,0),correct=T,simulate.p.value = F)$p.value
        phip.valuesim<-chisq.test(ifelse(temp[,b]>0,1,0),ifelse(temp[,c]>0,1,0),correct=F,simulate.p.value = T)$p.value}else{
          phi<-0
          phip.value<-1
          phip.valuecor<-1
          phip.valuesim<-1}
        
        test<-cor.test(temp[,b],temp[,c],method="pearson")
        pearsonr<-test$estimate
        pearsonp.value<-test$p.value
        
        ind<-which(temp[,b]>0&temp[,c]>0)#still might give an NA if there is no variation
        if(length(ind)>2){
          test<-cor.test(temp[ind,b],temp[ind,c],method="pearson")
          pearsonrno0<-test$estimate
          pearsonp.valueno0<-test$p.value}
        if(length(ind)<=2){
          pearsonrno0<-0
          pearsonp.valueno0<-1}
 
        ind<-which(temp[,b]+temp[,c]>0)
        if(length(ind)>2){
          test<-cor.test(temp[ind,b],temp[ind,c],method="pearson")
          pearsonrno00<-test$estimate
          pearsonp.valueno00<-test$p.value}
        if(length(ind)<=2){
          pearsonrno00<-0
          pearsonp.valueno00<-1}
        
        #do randomization test for p value
        #ran1<-t(replicate(1000,sample(temp[,b],length(temp[,b]))))
        #ran2<-t(replicate(1000,sample(temp[,c],length(temp[,c]))))
        #ken<-apply(cbind(ran1,ran2),1,function(x){Kendall(x[1:length(temp[,b])],x[(length(temp[,b])+1):(length(temp[,b])*2)])$tau})
        #p.value.perm<-length(which(ken>rho))/1000
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
        kendallrho<-0
        kendallp.value<-1  
        phi<-0
        phip.value<-1
        phip.valuecor<-1
        phip.valuesim<-1
        pearsonr<-0
        pearsonp.value<-1
        pearsonrno0<-0
        pearsonp.valueno0<-1
        pearsonrhono00<-0
        pearsonp.valueno00<-1
      }	
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],spearmanrho,spearmanp.value,kendallrho,kendallp.value,phi,phip.value,phip.valuecor,phip.valuesim,pearsonr,pearsonp.value,pearsonrno0,pearsonp.valueno0,pearsonrno00,pearsonp.valueno00,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
      new.row		
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt)
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:21],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","kendallrho","kendallp.value","phi","phip.value","phip.valuecor","phip.valuesim","pearsonr","pearsonp.value","pearsonrno0","pearsonp.valueno0","pearsonrno00","pearsonp.valueno00","ab1","ab2","ab1freq","ab2freq")



#Setup parallel backend to use 4 processors, for spearman permutation and kendall metrics
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
#loop
results<-matrix(nrow=0,ncol=12)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.data, lomehi==trt.temp)
  
  #in this case the community data started at column 28, so the loop for co-occurrence has to start at that point
  for(b in 28:(dim(temp)[2]-1)){
    
    results1<-foreach(c=(b+1):(dim(temp)[2]),.combine=rbind,.packages="Kendall") %dopar% {
      #for(c in 29:(dim(temp)[2])){      #for debugging
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(temp[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)
      #I changed this so that it will calculate the correlation for all species pairs, even singletons, and I can subset them later if I want to remove infrequent otus for qvalue calculation
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0

        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
        #do randomization test for p value
        ran1<-t(replicate(1000,sample(temp[,b],length(temp[,b]))))
        ran2<-t(replicate(1000,sample(temp[,c],length(temp[,c]))))
        spe<-apply(cbind(ran1,ran2),1,function(x){cor.test(x[1:length(temp[,b])],x[(length(temp[,b])+1):(length(temp[,b])*2)])$estimate})
        spearmanp.valueperm<-length(which(spe>spearmanrho))/1000
        
        test<-Kendall(temp[,b],temp[,c])
        kendallrho<-test$tau
        kendallp.value<-test$sl
        
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
        spearmanp.valueperm<-1
        kendallrho<-0
        kendallp.value<-1  
      }	
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],spearmanrho,spearmanp.value,spearmanp.valueperm,kendallrho,kendallp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
      new.row		
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt)
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:12],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","spearmanp.valueperm","kendallrho","kendallp.value","ab1","ab2","ab1freq","ab2freq")

head(results)
head(resultsold)

resultstrialorder<-subset(results,ab1freq>4&ab2freq>4&rho>0)
resultsorder<-subset(results,ab1freq>4&ab2freq>4&rho>0)
resultsorderm<-subset(results,ab1freq>4&ab2freq>4)
resultsorderm2<-subset(results,ab1freq>4&ab2freq>4)
resultsordersp<-subset(results,ab1freq>4&ab2freq>4)

co_occur_pairs_all<-function(dataset){
  final.results<-data.frame()
  trts<-as.vector(unique(dataset$trt))
  for(t in 1:length(trts)){
    dataset_trt<-subset(dataset, trt==trts[t])
    dataset_trt_no0<-subset(dataset_trt, ab1 > 0 & ab2 > 0)#Ryan had this at 0
    dataset_trt_no0$pairs<-paste(dataset_trt_no0$taxa1,dataset_trt_no0$taxa2)
    temp<-dataset_trt_no0
    if(dim(temp)[1]>1){
      temp$qval1<-fdrtool(temp$p.value,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
      temp$qval2<-p.adjust(temp$p.value.perm, method="fdr")#using a different method for p-value adjustment. the fdrtool wasn't working when I had a lot of p values that were significant and none around 1
      temp<-cbind(pairs=temp$pairs,temp[,1:10],temp[,12:13])
      #temp<-cbind(pairs=temp$pairs,temp[,1:7],temp[,9:10])
      final.results<-rbind(final.results,temp)	}
    print(t/length(trts))
  }
  return(final.results)
}

edge_listsorder<-co_occur_pairs_all(resultstrialorder)
edge_listsorder$qval3<-p.adjust(edge_listsorder$p.value.perm,method="fdr")
edge_listsorderm2<-cbind(pairs=paste(resultsorderm2$taxa1,resultsorderm2$taxa2),resultsorderm2)
edge_listsordersp<-cbind(pairs=paste(resultsordersp$taxa1,resultsordersp$taxa2),resultsordersp)

head(edge_listsordersp)

#summary: even at the order level, where there should be fewer zeros, we are having the same problem of hi frequency, low rho, low p

edge_listsorderm1<-subset(edge_listsorderm,spearmanrho>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$spearmanp.value,method="fdr")

edge_listsorderm1<-subset(edge_listsorderm,kendallrho>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$kendallp.value,method="fdr")

edge_listsorderm1<-subset(edge_listsorderm,phi>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$phip.valuesim,method="fdr")

edge_listsorderm1<-subset(edge_listsorderm,pearsonr>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$pearsonp.value,method="fdr")

edge_listsorderm1<-subset(edge_listsorderm,pearsonrno0>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$pearsonp.valueno0,method="fdr")

edge_listsorderm1<-subset(edge_listsorderm,pearsonrno00>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$pearsonp.valueno00,method="fdr")

#the more zeros you delete the fewer significant interaction you get

inputhiv<-subset(edge_listsorderm1,qval<0.05&trt=="hi")
inputlov<-subset(edge_listsorderm1,qval<0.05&trt=="lo")
inputnov<-subset(edge_listsorderm1,qval<0.05&trt=="no")

dim(inputhiv)
dim(inputlov)
dim(inputnov)
mean(inputhiv$ab1freq)
mean(inputlov$ab1freq)
mean(inputnov$ab1freq)
mean(inputhiv$spearmanrho)
mean(inputlov$spearmanrho)
mean(inputnov$spearmanrho)



edge_listsorderm1<-subset(edge_listsorderm2,spearmanrho>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$spearmanp.value,method="fdr")

edge_listsorderm1<-subset(edge_listsordersp,spearmanrho>0.5)
edge_listsorderm1<-subset(edge_listsordersp,spearmanp.value<0.01)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$spearmanp.value,method="fdr")
edge_listsorderm1$qval<-fdrtool(edge_listsorderm1$spearmanp.value,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval

edge_listsorderm1<-subset(edge_listsordersp,kendallrho>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$kendallp.value,method="fdr")

edge_listsorderm1<-subset(edge_listsorderm2,pearsonr>0)
edge_listsorderm1$qval<-p.adjust(edge_listsorderm1$pearsonp.value,method="fdr")

inputhiv<-subset(edge_listsorderm1,qval<0.1&trt=="hi")
inputmev<-subset(edge_listsorderm1,qval<0.1&trt=="me")
inputlov<-subset(edge_listsorderm1,qval<0.1&trt=="lo")
inputhiv<-subset(edge_listsorderm1,trt=="hi")
inputmev<-subset(edge_listsorderm1,trt=="me")
inputlov<-subset(edge_listsorderm1,trt=="lo")

dim(inputhiv)
dim(inputmev)
dim(inputlov)
mean(inputhiv$ab1freq)
mean(inputmev$ab1freq)
mean(inputlov$ab1freq)
mean(inputhiv$spearmanrho)
mean(inputmev$spearmanrho)
mean(inputlov$spearmanrho)


#high density
inputhi<-subset(edge_listsordersp,spearmanp.value<0.01&trt=="hi")[,3:4]#
inputhiv<-subset(edge_listsordersp,spearmanp.value<0.01&trt=="hi")#
vertexsizes<-unique(data.frame(orders=c(as.character(inputhiv$taxa1),as.character(inputhiv$taxa2)),abun=c(inputhiv$ab1,inputhiv$ab2)))

graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph1)
#plot(graph1,vertex.size=2,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"orders" #change to "orders" if doing things by order
colorgraph1<-unique(merge(verticesgraph1,labelsall[,c(1,3,4)],"orders",all.y=F,all.x=F,sort=F))
sizesgraph1<-merge(verticesgraph1,vertexsizes,"orders",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityorderspearp.01.pdf") #f=frequency cutoff 5, r=rho cutoff .5
plot(graph1,vertex.size=log(sizesgraph1$abun)*2,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)
#dev.off()


#med density
inputme<-subset(edge_listsordersp,spearmanp.value<0.01&trt=="me")[,3:4]
inputmev<-subset(edge_listsordersp,spearmanp.value<0.01&trt=="me")
inputme<-subset(edge_listsordersp,spearmanrho>0.5&trt=="me")[,3:4]
inputmev<-subset(edge_listsordersp,spearmanrho>0.5&trt=="me")
vertexsizes<-unique(data.frame(orders=c(as.character(inputmev$taxa1),as.character(inputmev$taxa2)),abun=c(inputmev$ab1,inputmev$ab2)))

graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"orders"
colorgraph2<-unique(merge(verticesgraph2,labelsall[,c(1,3,4)],"orders",all.y=F,all.x=F,sort=F))
sizesgraph2<-merge(verticesgraph2,vertexsizes,"orders",sort=F)
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuspearmanp.01.pdf")
plot(graph2,vertex.size=log(sizesgraph2$abun)*2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#
dev.off()



#lo density
inputlo<-subset(edge_listsordersp,spearmanp.value<0.01&trt=="lo")[,3:4]
inputlov<-subset(edge_listsordersp,spearmanp.value<0.01&trt=="lo")
inputlo<-subset(edge_listsordersp,spearmanrho>0.5&trt=="lo")[,3:4]
inputlov<-subset(edge_listsordersp,spearmanrho>0.5&trt=="lo")
vertexsizes<-unique(data.frame(orders=c(as.character(inputlov$taxa1),as.character(inputlov$taxa2)),abun=c(inputlov$ab1,inputlov$ab2)))

graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"orders"
colorgraph3<-unique(merge(verticesgraph3,labelsall[,c(1,3,4)],"orders",all.y=F,all.x=F,sort=F))
sizesgraph3<-merge(verticesgraph3,vertexsizes,"orders",sort=F)
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuspearmanp.01.pdf")
plot(graph3,vertex.size=log(sizesgraph3$abun)*2,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#
dev.off()



inputhiv<-subset(edge_listsordersp,spearmanrho>.5&trt=="hi")
inputmev<-subset(edge_listsordersp,spearmanrho>.5&trt=="me")
inputlov<-subset(edge_listsordersp,spearmanrho>.5&trt=="lo")

inputhiv<-subset(edge_listsordersp,spearmanp.value<.01&trt=="hi")
inputmev<-subset(edge_listsordersp,spearmanp.value<.01&trt=="me")
inputlov<-subset(edge_listsordersp,spearmanp.value<.01&trt=="lo")

graphhi<-simplify(graph.edgelist(as.matrix(inputhiv[,3:4]),directed=FALSE))
graphme<-simplify(graph.edgelist(as.matrix(inputmev[,3:4]),directed=FALSE))
graphlo<-simplify(graph.edgelist(as.matrix(inputlov[,3:4]),directed=FALSE))

length(V(graphhi))
length(V(graphme))
length(V(graphlo))
length(E(graphhi))
length(E(graphme))
length(E(graphlo))
graph.density(graphhi)
graph.density(graphme)
graph.density(graphlo)
transitivity(graphhi, type="global")#clustering, triangle
transitivity(graphme, type="global")
transitivity(graphlo, type="global")
modularity(graphhi, membership(walktrap.community(graphhi)))
modularity(graphme, membership(walktrap.community(graphme)))
modularity(graphlo, membership(walktrap.community(graphlo)))
