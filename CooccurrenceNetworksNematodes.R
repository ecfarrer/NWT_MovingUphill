#correlations including the nematode data (floating from 20g subsample and then sequenced)




#Read in nematode data
nematodecomp<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Nematodes/Nematode_species.csv")
head(nematodecomp)

#Remove nematodes only present in one or two plots
dim(nematodecomp2)
nematodecomp2<-nematodecomp[,colSums(nematodecomp>0)>2]


#make label file
nematodelabels1<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Nematodes/Nematode_species_info.csv")
nematodelabels<-as.data.frame(cbind(as.character(nematodelabels1$otu),as.character(nematodelabels1$NCBI_Taxonomy),as.character(nematodelabels1$Trophic_Group),as.character(nematodelabels1$otu)))
colnames(nematodelabels)<-c("otu","orders","kingdomlabels","otuxy")
nematodelabels

#use comm.dataALL from the CooccurrenceNetworksEuksbac.R script and merge with nematode community data. There are 92 samples in common between the nematode and the bac/euk sequencing dataset
comm.dataALLn<-merge(comm.dataALL,nematodecomp2,"Sample_name",sort=F,all.y=F,all.x=F)
comm.dataALLn$Sample_name
head(comm.dataALLn)[,1:30]

#same as before
trts



#Setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
results<-matrix(nrow=0,ncol=11)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp1<-subset(comm.dataALLn, lomehi==trt.temp)####change to comm.data if not doing plants
  
  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  temp2<-cbind(temp1[,1:27],temp1[,((which(colSums(temp1[,28:dim(temp1)[2]]>0)>2))+27)])
  ind<-which(substr(colnames(temp2),start=1,stop=1)=="n")
  tempnem<-temp2[,ind]
  temp<-temp2[,-ind]
  
  #take nematode data, for each nematode otu do correlations with every bacteria/plant otu. remember bact community data started at column 28, so the loop for co-occurrence has to start at that point
  for(b in 1:(dim(tempnem)[2])){
    results1<-foreach(c=28:(dim(temp)[2]),.combine=rbind,.packages="Kendall") %dopar% {
      species1.ab<-sum(tempnem[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(tempnem[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)
      #I changed this so that it will calculate the correlation for all species pairs (except doubletons and singletons), and I can subset them later if I want to remove infrequent otus for qvalue calculation
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(tempnem[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
        test<-Kendall(tempnem[,b],temp[,c])
        kendallrho<-test$tau
        kendallp.value<-test$sl
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
        kendallrho<-0
        kendallp.value<-1
      }	
      new.row<-c(trts[a],names(tempnem)[b],names(temp)[c],spearmanrho,spearmanp.value,kendallrho,kendallp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
      new.row		
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt) #nematodes vs Euk, bact, plants, took 7min with 4 cores
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","kendallrho","kendallp.value","ab1","ab2","ab1freq","ab2freq")


head(results)

resultsN<-results

