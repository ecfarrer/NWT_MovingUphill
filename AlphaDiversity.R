#Alpha diversity

#library(BiodiversityR) #this requires X11 and takes a while to load, you need to close the window that it opens in rcommander
library(picante)

#Dorota says to rarefy data first before calculating alpha diversity. then she said maybe not rarefy on anything because we do so much standardization in the lab protocol that number differences likely mean something
#I can also do multiple rarefactions and then take the mean for each "treatment"


#Euks- the input file is dats2 or dats2otu (the latter has the sample/mapping data attached to the beginning)
#Bacteria - the input file is dat16Ss2a or dat16Ss2a_s, dat16Ss2a_o (these have the 3 samples that didn't amplify removed). the sample and otu files for bacteria could not be cbound because too large

#First take out samples that did not amplify in one or other dataset: euks sample 81 did not amplify. for bacteria, samples 126, 5, 34 did not amplify.
#then sort them both so they are in the same order

comm.dataalphaEuk1<-dats2otu
comm.dataalphaEuk1$Sample_name<-as.numeric(as.character(comm.dataalphaEuk1$Sample_name))
comm.dataalphaEuk<-comm.dataalphaEuk1%>%
  filter(Sample_name!=5&Sample_name!=34&Sample_name!=81&Sample_name!=126)%>%
  arrange(Sample_name)
head(comm.dataalphaEuk)[,1:30]

dat16Ss2a_s$Sample_name<-as.numeric(as.character(dat16Ss2a_s$Sample_name))
ind<-which(dat16Ss2a_s$Sample_name!=5&dat16Ss2a_s$Sample_name!=34&dat16Ss2a_s$Sample_name!=81&dat16Ss2a_s$Sample_name!=126)
comm.dataalpha16S_o1<-dat16Ss2a_o[ind,]
comm.dataalpha16S_s1<-dat16Ss2a_s[ind,]
comm.dataalpha16S_o<-comm.dataalpha16S_o1[order(comm.dataalpha16S_s1$Sample_name),]
comm.dataalpha16S_s<-comm.dataalpha16S_s1[order(comm.dataalpha16S_s1$Sample_name),]
head(comm.dataalpha16S_s)[,1:20]
head(comm.dataalpha16S_o)[,1:20]

lomehialpha<-ifelse(comm.dataalphaEuk$Plant_Dens<36,"lo","else");lomehialpha[which(comm.dataalphaEuk$Plant_Dens<89&comm.dataalphaEuk$Plant_Dens>=36)]<-"me";lomehialpha[which(comm.dataalphaEuk$Plant_Dens>=89)]<-"hi";lomehialpha<-factor(lomehialpha,levels=c("lo","me","hi"))

comm.dataalphaEuk<-cbind(lomehi=lomehialpha,comm.dataalphaEuk)
comm.dataalpha16S_s<-cbind(lomehi=lomehialpha,comm.dataalpha16S_s)

#how many to rarefy to
min(rowSums(comm.dataalphaEuk[,-c(1:27)]))#1280
min(rowSums(comm.dataalpha16S_o))#5697

#plantcomp2b<-plantcomp2%>%
#  filter(Sample_name%in%comm.data16Sb$Sample_name)%>%
#  arrange(Sample_name)




######Rarefy data 100 times and get mean diversity and richness######

#Euks
#Setup parallel backend to use 4 processors, takes 6 minutes
cl<-makeCluster(6)
registerDoParallel(cl)
strt<-Sys.time()

resultsalphaEuk<-foreach(a=1:100,.combine=cbind,.packages="vegan") %dopar% {
  comm.temp<-cbind(comm.dataalphaEuk[,c(1:27)],rrarefy(comm.dataalphaEuk[,-c(1:27)],1280))
  cbind(vegan::diversity(comm.temp[,-c(1:27)]),specnumber(comm.temp[,-c(1:27)]))
}

head(resultsalphaEuk)
print(Sys.time()-strt)
stopCluster(cl)

resultsdivEuk<-rowSums(resultsalphaEuk[,seq(from=1,to=199,by=2)])/100
restultsrichEuk<-rowSums(resultsalphaEuk[,seq(from=2,to=200,by=2)])/100

plot(comm.dataalphaEuk$lomehi,restultsrichEuk)
plot(comm.dataalphaEuk$Plant_Div,restultsrichEuk)
plot(comm.dataalphaEuk$Plant_Dens,restultsrichEuk)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/eukrichnessvsplantdensity.pdf")
ggplot(cbind(comm.dataalphaEuk,rich=restultsrichEuk,div=resultsdivEuk),aes(x=log10(Plant_Dens+1),y=rich))+#as.numeric(fert),color=species
  #scale_x_log10(breaks=c(1,10,100,1000)) + 
  ylim(0,500) +
  labs(x="Plant density",y="Euk richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=20))+
  geom_point(size=3)+
  geom_smooth(method=lm,se=F,size=1.5,color="black") 
dev.off()

summary(lm(restultsrichEuk~log((comm.dataalphaEuk$Plant_Dens+1),base=10)))


#Bacteria
#Setup parallel backend to use 4 processors, took 15 minutes
cl<-makeCluster(6)
registerDoParallel(cl)
strt<-Sys.time()

resultsalpha16S<-foreach(a=1:100,.combine=cbind,.packages="vegan") %dopar% {
  comm.temp<-rrarefy(comm.dataalpha16S_o,5697)
  cbind(vegan::diversity(comm.temp),specnumber(comm.temp))
}

head(resultsalpha16S)
print(Sys.time()-strt)
stopCluster(cl)

resultsdiv16S<-rowSums(resultsalpha16S[,seq(from=1,to=199,by=2)])/100
restultsrich16S<-rowSums(resultsalpha16S[,seq(from=2,to=200,by=2)])/100

plot(comm.dataalpha16S_s$lomehi,restultsrich16S)
plot(comm.dataalpha16S_s$Plant_Div,restultsrich16S)
plot(comm.dataalpha16S_s$Plant_Dens,restultsrich16S)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/bacrichnessvsplantdensity.pdf")
ggplot(cbind(comm.dataalpha16S_s,rich=restultsrich16S,div=resultsdiv16S),aes(x=log10(Plant_Dens+1),y=rich))+#as.numeric(fert),color=species
  #scale_y_log10() +##ylim(0,5) +#
  ylim(0,2000) +
  labs(x="Plant density",y="Bacteria richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=20))+
  geom_point(size=3)+
  geom_smooth(method=lm,se=F,size=1.5,color="black") 
dev.off()

summary(lm(restultsrich16S~log((comm.dataalpha16S_s$Plant_Dens+1),base=10)))


#merge data so I can graph them on panels
richdata<-as.data.frame(rbind(cbind(Plant_Dens=comm.dataalpha16S_s$Plant_Dens,rich=restultsrich16S),cbind(Plant_Dens=comm.dataalphaEuk$Plant_Dens,rich=restultsrichEuk)));rownames(richdata)<-1:dim(richdata)[1]
richdata$type<-factor(rep(c("Bac","Euk"),each=94))
richdata$lomehi<-comm.dataalpha16S_s$lomehi
richdata$Sample_name<-comm.dataalphaEuk$Sample_name
richdata$X.SampleID<-comm.dataalphaEuk$X.SampleID

#I don't think you can set different limits for each panel, for example I'd like them both to have ymin of 0
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/microberichnessvsplantdensity.pdf",width=7, height=3.5)
ggplot(richdata,aes(x=log10(Plant_Dens+1),y=rich))+#as.numeric(fert),color=species
  labs(x="Plant density",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  facet_wrap(~type,scales="free")
dev.off()


#aggregate by lomehi so I can make a barplot
richdata2<- richdata %>% group_by(type,lomehi) %>%
  summarise(mean_rich = mean(rich),se_rich=std.error(rich))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/microberichnessvsplantdensitygroups.pdf",width=7, height=3.5)
ggplot(richdata2,aes(x=lomehi,y=mean_rich,group=type))+
  labs(x = "Plant density",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_line(stat = "identity", position = "identity",size=.8)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean_rich+se_rich, ymin=mean_rich-se_rich),width=.15,size=.8) +
  facet_wrap(~type,scales="free")
dev.off()





######Phylogenetic diversity######
#file is pdEuk from 18S_DataCleaning script and pd16S
#need to add bac when server is done
pdEuk$X.SampleID<-rownames(pdEuk)
richdataEuk<-subset(richdata,type=="Euk")
richdataEuk2<-merge(richdataEuk,pdEuk,by="X.SampleID",all.y=F,sort=F)

pd16S$X.SampleID<-rownames(pd16S)
richdataBac<-subset(richdata,type=="Bac")
richdataBac2<-merge(richdataBac,pd16S,by="X.SampleID",all.y=F,sort=F)

pdN<-read.delim("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Nematodes/Euk_Nema_truncate_99_alphadiv.txt",row.names = 1)
temp<-t(as.data.frame(strsplit(rownames(pdN),"[.]")))
temp2<-paste("S",temp[,2],temp[,3],sep=".")
pdN$X.SampleID<-temp2
pdN2<-merge(richdataEuk2[,1:6],pdN,all.x = F,all.y = F,sort=F) #
pdN3<-cbind(pdN2[,1:6],PD=pdN2$PD_whole_tree,SR=NA)
pdN3$type<-"Nematode"
pdN3$rich<-NA

richdata3<-rbind(richdataEuk2,richdataBac2,pdN3)

#fig with euks and bac
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/microbeNPDvsplantdensity.pdf",width=7, height=3) #for two panels width=7, height=3.5)
ggplot(richdata3,aes(x=log10(Plant_Dens+1),y=PD))+#as.numeric(fert),color=species
  labs(x="Plant density",y="Phylogenetic diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  facet_wrap(~type,scales="free")
dev.off()

summary(lm(richdataBac2$PD~log((richdataBac2$Plant_Dens+1),base=10)))
summary(lm(richdataEuk2$PD~log((richdataEuk2$Plant_Dens+1),base=10)))
summary(lm(pdN3$PD~log((pdN3$Plant_Dens+1),base=10)))

#plot with dashed vertical lines for low med hi, just euks, trying it out, I don't really like how it looks
ggplot(richdataEuk2,aes(x=log10(Plant_Dens+1),y=PD))+#as.numeric(fert),color=species
  labs(x="Plant density",y="Phylogenetic diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  geom_vline(aes(xintercept=log10(36+1)), colour="gray30", linetype="dashed")+
  geom_vline(aes(xintercept=log10(89+1)), colour="gray30", linetype="dashed")#+facet_wrap(~type,scales="free")









#Euks
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



#Bacteria
#the input file is dat16Ss2a or dat16Ss2a_s, dat16Ss2a_o (these have the 3 samples that didn't amplify removed)

comm.data16Salpha1_o<-dat16Ss2a_o

min(rowSums(comm.data16Salpha1_o))#5697
comm.data16Salpha_o<-cbind(lomehif,comm.data16Salpha1)


######Rarefy data 100 times and get mean diversity and richness######
#Setup parallel backend to use 4 processors
cl<-makeCluster(6)
registerDoParallel(cl)
strt<-Sys.time()

results16Salpha<-foreach(a=1:100,.combine=cbind,.packages="vegan") %dopar% {
  comm.temp<-cbind(comm.data16Salpha[,c(1:27)],rrarefy(comm.data16Salpha[,-c(1:27)],1280))
  cbind(vegan::diversity(comm.temp[,-c(1:27)]),specnumber(comm.temp[,-c(1:27)]))
}

head(results16Salpha)
print(Sys.time()-strt)
stopCluster(cl)

results16Sdiv<-rowSums(results16Salpha[,seq(from=1,to=199,by=2)])/100
restults16Srich<-rowSums(results16Salpha[,seq(from=2,to=200,by=2)])/100

plot(comm.data16Salpha$lomehif,restults16Srich)
plot(comm.data16Salpha$Plant_Div,restults16Srich)
plot(comm.data16Salpha$Plant_Dens,restults16Srich)

ggplot(cbind(comm.data16Salpha,rich=restults16Srich,div=results16Sdiv),aes(x=Plant_Div,y=rich))+#as.numeric(fert),color=species
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
rankabuncomp(comm.dataEukc[,-c(1:27)], y=comm.dataEukc[,c(1:27)], factor='lomehi', legend=F,xlim=c(1,770),ylim=c(1,20))#, scale='proportion'
rankabuncomp(comm.data16Sc[,-c(1:27)], y=comm.data16Sc[,c(1:27)], factor='lomehi', legend=F,xlim=c(1,1770),ylim=c(1,100))#, scale='proportion' 

comm.data16Sc
comm16S.spe
specnumber(comm.data16Sc[,-c(1:27)])

aggregate.data.frame(rowSums(comm16S.spe>0),by=list(comm.data16Sc$lomehi),mean)

mean(specnumber(comm16S.spelo2))
mean(specnumber(comm16S.speme2))
mean(specnumber(comm16S.spehi2))
mean(specnumber(commEuk.spelo2))
mean(specnumber(commEuk.speme2))
mean(specnumber(commEuk.spehi2))



