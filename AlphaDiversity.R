#Alpha diversity

library(BiodiversityR) #this requires X11 and takes a while to load, you need to close the window that it opens in rcommander


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

plot(comm.dataalpha$lomehif,restultsrich)
plot(comm.dataalpha$Plant_Div,restultsrich)
plot(comm.dataalpha$Plant_Dens,restultsrich)

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
rarecurve(comm.dataalpha[,-c(1:27)],step=100,xlab="Sample size",ylab="OTU",label=F,col=lomehi)#,sample=rowSums(comm.dataalpha[,-c(1:27)]) #takes a few minutes to make




######compare species abundance distributions in lo me hi######
#comm.data16Sc and comm.dataEukc come from the BetaDiversity.R script

#RankAbun.1 <- rankabundance(dune)
#rankabunplot(RankAbun.1,scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabuncomp(comm.dataEukc[,-c(1:27)], y=comm.dataEukc[,c(1:27)], factor='lomehi', legend=F,xlim=c(1,70),ylim=c(1,2000))#, scale='proportion'

comm.data16Sc
comm16S.spe
specnumber(comm.data16Sc[,-c(1:27)])

aggregate.data.frame(rowSums(comm16S.spe>0),by=list(comm.data16Sc$lomehi),mean)

