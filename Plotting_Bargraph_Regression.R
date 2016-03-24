
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata")

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata")


library(tidyr)


######Diversity caluculations on non-rarefied data######
#I think we do need to rarefy, there is a huge impact of sequencing depth on richness (not a strong effect on diversty though)

dats2
lomehif<-factor(lomehi,levels=c("lo","me","hi"))
rich<-estimate_richness(dats2,split=T, measures=c("Observed","Shannon"))# Observed is the number of OTUs
temp2<-sample_data(dats2)
plot(log(temp2$Plant_Dens),rich$Observed)
plot(temp2$Plant_Div,rich$Observed)

m1<-aggregate.data.frame(rich$Observed, by=list(lomehif),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(lomehif),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)

unique(kingdomlabels)

fungi<-subset_taxa(dats2,kingdomgroup=="Fungi")
fungi
rich<-estimate_richness(fungi,split=T, measures=c("Observed","Shannon"))
plot(temp2$Plant_Dens,rich$Shannon)
m1<-aggregate.data.frame(rich$Observed, by=list(lomehif),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(lomehif),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)

archaeplastida<-subset_taxa(dats2,kingdomgroup=="Archaeplastida")
rich<-estimate_richness(archaeplastida,split=T, measures=c("Observed","Shannon"))
plot(temp2$Plant_Dens,rich$Shannon)
m1<-aggregate.data.frame(rich$Observed, by=list(lomehif),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(lomehif),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)

archaeplastida<-subset_taxa(dats2,kingdomgroup=="Archaeplastida")
rich<-estimate_richness(archaeplastida,split=T, measures=c("Observed","Shannon"))
plot(temp2$Plant_Dens,rich$Shannon)
m1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)









######Changes in relative abundance of kingdoms######
dats9kingdom #a file of relative abundance by kingdom output in the data cleaning r file
#dats2otu<-cbind(sample_data(dats2),t(otu_table(dats2))) ##takes a long time
#relabunfunc<-function(x){x/rowSums(x)*100}
#dats2otusprel<-dats2otu %>% select(denovo3:denovo358659) %>% relabunfunc()
dats9kingdom<-cbind(lomehi=factor(lomehi,levels=c("lo","me","hi")),dats9kingdom)

plotdata<-dats9kingdom %>% 
  select(lomehi,X.SampleID, Plant_Div, Plant_Dens,Amoebozoa:Rhizaria) %>%
  #mutate(Fungi=rowSums(dats2otur[,which(labelsall$kingdomlabels=="Fungi")])) %>%
  gather(species,abun,Amoebozoa:Rhizaria) %>%
  group_by(species,lomehi) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun))
data.frame(plotdata)

plotdatalarge<-plotdata %>%
  filter(species%in%c("Fungi","Nonphotosynthetic_Alveolata","Rhizaria","Archaeplastida"))

ggplot(plotdatalarge,aes(x=as.numeric(lomehi),y=mean_abun,color=species))+#as.numeric(fert)
  #scale_y_log10() +##ylim(0,5) +#
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=15))+
  geom_line(stat = "identity", position = "identity",size=1.5)+
  geom_point(size=4)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.25,size=1.5)

plotdatasmall<-plotdata %>%
  filter(species%in%c("Amoebozoa","Discicristoidea","Holozoa","Nonphotosynthetic_Discoba","Nonphotosynthetic_Eukaryota","Photosynthetic_Alveolata","Photosynthetic_Discoba"))

ggplot(plotdatasmall,aes(x=as.numeric(lomehi),y=mean_abun,color=species))+#as.numeric(fert)
  #scale_y_log10() +##ylim(0,5) +#
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=15))+
  geom_line(stat = "identity", position = "identity",size=1.5)+
  geom_point(size=4)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.25,size=1.5)











#Input dataset is microbplant which has no rarefaction for microbes and includes all plants including doubletons/singletons

species<-dats6order[,27:243]
diversity<-vegan::diversity(species)
richness<-rowSums(species>0)
m1<-aggregate.data.frame(diversity, by=list(greater66plants),mean)
se1<-aggregate.data.frame(diversity, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)
m1<-aggregate.data.frame(richness, by=list(greater66plants),mean)
se1<-aggregate.data.frame(richness, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)
plot(dats6order$Plant_Div,diversity)
