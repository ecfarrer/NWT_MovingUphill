
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




######Bacteria and Euks######
#input files are: dats9kingdom and dat16Ss9kingdom, raw data only single reads removed, plants and proks removed, then relative abundance was calculated
comm.datarelEuk1<-dats9kingdom[,-c(1)]
comm.datarelEuk1$Sample_name<-as.numeric(as.character(comm.datarelEuk1$Sample_name))
comm.datarelEuk<-comm.datarelEuk1%>%
  filter(Sample_name!=5&Sample_name!=34&Sample_name!=81&Sample_name!=126)%>%
  arrange(Sample_name)
head(comm.datarelEuk)[,1:30]

comm.datarel16S1<-dat16Ss9kingdom
comm.datarel16S1$Sample_name<-as.numeric(as.character(comm.datarel16S1$Sample_name))
comm.datarel16S<-comm.datarel16S1%>%
  filter(Sample_name!=5&Sample_name!=34&Sample_name!=81&Sample_name!=126)%>%
  arrange(Sample_name)

lomehirel<-ifelse(comm.datarelEuk$Plant_Dens<36,"lo","else");lomehirel[which(comm.datarelEuk$Plant_Dens<89&comm.datarelEuk$Plant_Dens>=36)]<-"me";lomehirel[which(comm.datarelEuk$Plant_Dens>=89)]<-"hi";lomehirel<-factor(lomehirel,levels=c("lo","me","hi"))

comm.datarelEuk<-cbind(lomehi=lomehirel,comm.datarelEuk)
comm.datarel16S<-cbind(lomehi=lomehirel,comm.datarel16S)
#comm.datarelALL<-cbind(comm.datarel16S,comm.datarelEuk[,28:38])

comm.datarelEukl<-comm.datarelEuk %>% 
  dplyr::select(lomehi,X.SampleID, Plant_Div, Plant_Dens,Amoebozoa:Rhizaria) %>%
  gather(Taxa,abun,Amoebozoa:Rhizaria) %>%
  mutate(type="euk")
comm.datarel16Sl<-comm.datarel16S %>% 
  dplyr::select(lomehi,X.SampleID, Plant_Div, Plant_Dens,Acidobacteria:Verrucomicrobia) %>%
  gather(Taxa,abun,Acidobacteria:Verrucomicrobia) %>%
  mutate(type="bac")

comm.datarellALL<-rbind(comm.datarel16Sl,comm.datarelEukl)
head(comm.datarellALL)

plotdata<-comm.datarellALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  group_by(Taxa,lomehi,type,typeTaxa) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) %>%
  filter(mean_abun>4)
as.data.frame(plotdata)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroups.pdf",width=7, height=3.5)
ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_line(stat = "identity", position = "identity",size=.8)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.8)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,scales="free")
dev.off()

#http://tools.medialab.sciences-po.fr/iwanthue/
mycols<-c("#4BC366",
          "#D9A125",
          "#659125",
          "#6768A3",
          "#5C426C",
          "#D185E0",
          "#6F94DE",
          "#B4405E",
          "#D063A5",
          "#C25833",
          "#555516",
          "#8AD93B")

ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_line(stat = "identity", position = "identity",size=.8)+
  geom_point(size=3)+
  scale_color_brewer(palette="Paired") 





##euk only
######Changes in relative abundance of kingdoms######
dats9kingdom #a file of relative abundance by kingdom output in the data cleaning r file
#dats2otu<-cbind(sample_data(dats2),t(otu_table(dats2))) ##takes a long time
#relabunfunc<-function(x){x/rowSums(x)*100}
#dats2otusprel<-dats2otu %>% select(denovo3:denovo358659) %>% relabunfunc()
dats9kingdom<-cbind(lomehi=factor(lomehi,levels=c("lo","me","hi")),dats9kingdom)

plotdata<-dats9kingdom %>% 
  dplyr::select(lomehi,X.SampleID, Plant_Div, Plant_Dens,Amoebozoa:Rhizaria) %>%
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

#connecting lines for a factor
ggplot(plotdatalarge,aes(x=lomehi,y=mean_abun,color=species,group=species))+#as.numeric(fert)
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



dat16Ss9kingdom
dat16Ss9kingdom<-cbind(lomehi=factor(lomehi,levels=c("lo","me","hi")),dat16Ss9kingdom)







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
