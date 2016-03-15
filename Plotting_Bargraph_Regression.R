
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata")

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")




######Diversity caluculations on non-rarefied data######

dats2
rich<-estimate_richness(dats2,split=T, measures=c("Observed","Shannon"))# Observed is the number of OTUs
temp2<-sample_data(dats2)
plot(temp2$Plant_Dens,rich$Observed)

m1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)

unique(kingdomlabels)

fungi<-subset_taxa(dats2,kingdomgroup=="Fungi")
fungi
rich<-estimate_richness(fungi,split=T, measures=c("Observed","Shannon"))
plot(temp2$Plant_Dens,rich$Shannon)
m1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)

archaeplastida<-subset_taxa(dats2,kingdomgroup=="Archaeplastida")
rich<-estimate_richness(archaeplastida,split=T, measures=c("Observed","Shannon"))
plot(temp2$Plant_Dens,rich$Shannon)
m1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)

archaeplastida<-subset_taxa(dats2,kingdomgroup=="Archaeplastida")
rich<-estimate_richness(archaeplastida,split=T, measures=c("Observed","Shannon"))
plot(temp2$Plant_Dens,rich$Shannon)
m1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),mean)
se1<-aggregate.data.frame(rich$Observed, by=list(greater66plants),std.error)
plotCI(barplot(m1$x,names.arg=m1$Group.1),m1$x,uiw=se1$x,add=T,pch=NA)






######Rarefaction######
dats2r = rarefy_even_depth(dats2)
dats2rr = transform_sample_counts(dats2r, function(x) x/sum(x))
?make_network
#taxa_sums(dats2rr)#just sums OTUs in all samples








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
