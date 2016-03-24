#Alpha diversity


#Dorota says to rarefy data first before calculating alpha diversity
#I can also do multiple rarefactions and then take the mean for each "treatment"

#the input file is dats2 or dats2otu (the latter has the sample/mapping data attached to the beginning)

comm.dataalpha1<-dats2otu
min(rowSums(comm.dataalpha1[,-c(1:26)]))
comm.dataalpha<-cbind(lomehi,comm.dataalpha1)

#Setup parallel backend to use 4 processors
cl<-makeCluster(6)
registerDoParallel(cl)
strt<-Sys.time()

resultsalpha<-foreach(a=1:100,.combine=cbind,.packages="vegan") %dopar% {
  comm.temp<-cbind(comm.dataalpha[,c(1:27)],rrarefy(comm.dataalpha[,-c(1:27)],1280))
  cbind(vegan::diversity(comm.temp[,-c(1:27)]),specnumber(comm.temp[,-c(1:27)]))
  }
  
head(resultsalpha)
print(Sys.time()-strt)
stopCluster(cl)
  
resultsdiv<-rowSums(resultsalpha[,seq(from=1,to=199,by=2)])/100
restulsrich<-rowSums(resultsalpha[,seq(from=2,to=200,by=2)])/100

  