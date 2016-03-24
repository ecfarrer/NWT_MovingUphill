#Alpha diversity


#Dorota says to rarefy data first before calculating alpha diversity
#I can also do multiple rarefactions and then take the mean for each "treatment"

#the input file is dats2 or dats2otu (the latter has the sample/mapping data attached to the beginning)

comm.dataalpha1<-dats2otu
min(rowSums(comm.dataalpha1[,-c(1:26)]))
comm.dataalpha<-cbind(lomehif,comm.dataalpha1)




######Rarefy data 100 times and get mean diversity and richness######
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
restultsrich<-rowSums(resultsalpha[,seq(from=2,to=200,by=2)])/100

plot(comm.dataalpha$lomehi,restulsrich)
plot(comm.dataalpha$Plant_Div,restulsrich)
plot(comm.dataalpha$Plant_Dens,restulsrich)

ggplot(cbind(comm.dataalpha,rich=restultsrich,div=resultsdiv),aes(x=Plant_Div,y=rich))+#as.numeric(fert),color=species
  #scale_y_log10() +##ylim(0,5) +#
  labs(x="Plant richness",y="Euk richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=15))+
  geom_point(size=4)+
  geom_smooth(method=lm,se=F,size=1.5) 
#geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.25,size=1.5)
#geom_line(stat = "identity", position = "identity",size=1.5)+





######Make rarefaction curves for samples######
rarecurve(comm.dataalpha[,-c(1:27)],step=100,xlab="Sample size",ylab="OTU",label=F,col=lomehi)#,sample=rowSums(comm.dataalpha[,-c(1:27)])









