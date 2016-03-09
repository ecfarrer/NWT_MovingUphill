

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

#Standardize to 100 or rarify to 1250. if I standardize, I will need to change code in the loop below because it has cutoffs of 1
comm.data<-cbind(comm.data[,c(1:27)],rrarefy(comm.data[,-c(1:27)],1250))
#comm.data[,27:243]<-comm.data[,27:243]/rowSums(comm.data[,27:243])*100

trts<-as.vector(unique(greater66plants))

#Merge plants with microbes, plantcomp is everything, plantcomp2 removes doubletons/singletons
comm.data<-merge(comm.data,plantcomp2,"Sample_name",sort=F)





#Running it in parallel, setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)

#Record start time
strt<-Sys.time()
results<-matrix(nrow=0,ncol=7)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.data,greater66plants==trt.temp)
  
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
rownames(results)<-1:dim(results)[1]
results<-data.frame(results[,1:3],as.numeric(results[,4]),as.numeric(results[,5]),as.numeric(results[,6]),as.numeric(results[,7]))
names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")
print(Sys.time()-strt)#1 minute
stopCluster(cl)
head(results)
resultsm<-results






#Not running it in parallel
#Note: everythign below needs to be tweaked after I reanmed and added columns to the microbial data
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



















