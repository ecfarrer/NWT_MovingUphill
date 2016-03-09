
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata")

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")




######Diversity caluculations######

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
