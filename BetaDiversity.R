
library(dplyr)

######Bray-curtis BetaDiversity######

#These are rarefied but nothing else
comm.dataEuk[,1:30]
comm.data16S[,1:30]

#First take out samples that did not amplify in one or other dataset: euks sample 81 did not amplify. for bacteria, samples 126, 5, 34 did not amplify.
#then sort them both so they are in the same order
comm.dataEukb<-comm.dataEuk%>%
  filter(Sample_name!=5&Sample_name!=34&Sample_name!=81&Sample_name!=126)%>%
  arrange(Sample_name)
comm.data16Sb<-comm.data16S%>%
  filter(Sample_name!=5&Sample_name!=34&Sample_name!=81&Sample_name!=126)%>%
  arrange(Sample_name)
plantcomp2b<-plantcomp2%>%
  filter(Sample_name%in%comm.data16Sb$Sample_name)%>%
  arrange(Sample_name)
  

#Calculate lo me hi and put in front of each data.frame
lomehi2<-ifelse(comm.dataEukb$Plant_Dens<36,"lo","else");lomehi2[which(comm.dataEukb$Plant_Dens<89&comm.dataEukb$Plant_Dens>=36)]<-"me";lomehi2[which(comm.dataEukb$Plant_Dens>=89)]<-"hi";lomehi2<-factor(lomehi2, levels=c("lo","me","hi"))
lomehi2

comm.dataEukc<-cbind(lomehi=lomehi2,comm.dataEukb)
comm.data16Sc<-cbind(lomehi=lomehi2,comm.data16Sb)
plantcomp2c<-cbind(lomehi=lomehi2,plantcomp2b)

#separate species
#euks & 16S
commEuk.spe<-comm.dataEukc[,28:dim(comm.dataEukc)[2]]
commEuk.spelo<-subset(commEuk.spe,lomehi2=="lo")
commEuk.speme<-subset(commEuk.spe,lomehi2=="me")
commEuk.spehi<-subset(commEuk.spe,lomehi2=="hi")
comm16S.spe<-comm.data16Sc[,28:dim(comm.data16Sc)[2]]
comm16S.spelo<-subset(comm16S.spe,lomehi2=="lo")
comm16S.speme<-subset(comm16S.spe,lomehi2=="me")
comm16S.spehi<-subset(comm16S.spe,lomehi2=="hi")
commPlant.spe<-plantcomp2c[,3:dim(plantcomp2c)[2]]
commPlant.spelo<-subset(commPlant.spe,lomehi2=="lo")
commPlant.speme<-subset(commPlant.spe,lomehi2=="me")
commPlant.spehi<-subset(commPlant.spe,lomehi2=="hi")

otuEuklo<-colnames(commEuk.spelo)[(which(colSums(commEuk.spelo)>0))]
otuEukme<-colnames(commEuk.speme)[(which(colSums(commEuk.speme)>0))]
otuEukhi<-colnames(commEuk.spehi)[(which(colSums(commEuk.spehi)>0))]
length(intersect(otuEukme,otuEuklo))/length(union(otuEukme,otuEuklo))
length(intersect(otuEukhi,otuEukme))/length(union(otuEukhi,otuEukme))
length(intersect(otuEukhi,otuEuklo))/length(union(otuEukhi,otuEuklo))

otu16Slo<-colnames(comm16S.spelo)[(which(colSums(comm16S.spelo)>0))]
otu16Sme<-colnames(comm16S.speme)[(which(colSums(comm16S.speme)>0))]
otu16Shi<-colnames(comm16S.spehi)[(which(colSums(comm16S.spehi)>0))]
length(intersect(otu16Sme,otu16Slo))/length(union(otu16Sme,otu16Slo))
length(intersect(otu16Shi,otu16Sme))/length(union(otu16Shi,otu16Sme))
length(intersect(otu16Shi,otu16Slo))/length(union(otu16Shi,otu16Slo))



aggregate.data.frame(specnumber(comm.dataEukc[,-c(1:27)]),by=list(lomehi2),mean)
aggregate.data.frame(specnumber(comm.data16Sc[,-c(1:27)]),by=list(lomehi2),mean)


#over lap in abundance/frequent species. in lo me hi datasets I removed any species that had frequency of 10 or less. similar to network analysis input.
#euks & 16S
comm16S.spelo2<-comm16S.spelo[,which(colSums(comm16S.spelo>0)>10)]
comm16S.speme2<-comm16S.speme[,which(colSums(comm16S.speme>0)>10)]
comm16S.spehi2<-comm16S.spehi[,which(colSums(comm16S.spehi>0)>10)]

commEuk.spelo2<-commEuk.spelo[,which(colSums(commEuk.spelo>0)>10)]
commEuk.speme2<-commEuk.speme[,which(colSums(commEuk.speme>0)>10)]
commEuk.spehi2<-commEuk.spehi[,which(colSums(commEuk.spehi>0)>10)]

commPlant.spelo2<-commEuk.spelo[,which(colSums(commEuk.spelo>0)>10)]
commPlant.speme2<-commEuk.speme[,which(colSums(commEuk.speme>0)>10)]
commPlant.spehi2<-commEuk.spehi[,which(colSums(commEuk.spehi>0)>10)]

otuEuklo<-colnames(commEuk.spelo2)[(which(colSums(commEuk.spelo2)>0))]
otuEukme<-colnames(commEuk.speme2)[(which(colSums(commEuk.speme2)>0))]
otuEukhi<-colnames(commEuk.spehi2)[(which(colSums(commEuk.spehi2)>0))]
length(intersect(otuEukme,otuEuklo))/length(union(otuEukme,otuEuklo))
length(intersect(otuEukhi,otuEukme))/length(union(otuEukhi,otuEukme))
length(intersect(otuEukhi,otuEuklo))/length(union(otuEukhi,otuEuklo))

otu16Slo<-colnames(comm16S.spelo2)[(which(colSums(comm16S.spelo2)>0))]
otu16Sme<-colnames(comm16S.speme2)[(which(colSums(comm16S.speme2)>0))]
otu16Shi<-colnames(comm16S.spehi2)[(which(colSums(comm16S.spehi2)>0))]
length(intersect(otu16Sme,otu16Slo))/length(union(otu16Sme,otu16Slo))
length(intersect(otu16Shi,otu16Sme))/length(union(otu16Shi,otu16Sme))
length(intersect(otu16Shi,otu16Slo))/length(union(otu16Shi,otu16Slo))



aggregate.data.frame(specnumber(comm.dataEukc[,-c(1:27)]),by=list(lomehi2),mean)
aggregate.data.frame(specnumber(comm.data16Sc[,-c(1:27)]),by=list(lomehi2),mean)






#Only euk data
head(comm.data)[,1:30] #make sure comm.data has nohilo as the first column
dim(comm.data)

comm.spe<-comm.data[,28:dim(comm.data)[2]]
comm.spehi<-subset(comm.spe,lomehi=="hi")
comm.speme<-subset(comm.spe,lomehi=="me")
comm.spelo<-subset(comm.spe,lomehi=="lo")

comm.disthi<-vegdist(comm.spehi,method="jaccard",binary=T)
comm.distme<-vegdist(comm.speme,method="jaccard",binary=T)
comm.distlo<-vegdist(comm.spelo,method="jaccard",binary=T)

mean(comm.disthi)
mean(comm.distme)
mean(comm.distlo)
std.error(comm.disthi)
std.error(comm.distme)
std.error(comm.distlo)

alpha <- with(comm.data, tapply(specnumber(comm.spe), lomehi, mean))
gamma <- with(comm.data, specnumber(comm.spe, lomehi))
gamma/alpha - 1




######Taxonomic overlap######

dim(comm.spehi)

#Jaccard SJ = a/(a + b + c)
otuhi<-colnames(comm.spehi)[(which(colSums(comm.spehi)>0))]
otume<-colnames(comm.speme)[(which(colSums(comm.speme)>0))]
otulo<-colnames(comm.spelo)[(which(colSums(comm.spelo)>0))]
length(intersect(otuhi,otulo))/length(union(otuhi,otulo))
length(intersect(otuhi,otume))/length(union(otuhi,otume))
length(intersect(otume,otulo))/length(union(otume,otulo))

#Sorenson SS = 2a/(2a + b + c), this might be the more common
(2*length(intersect(otuhi,otulo)))/(length(intersect(otuhi,otulo))+length(union(otuhi,otulo)))
(2*length(intersect(otuhi,otume)))/(length(intersect(otuhi,otume))+length(union(otuhi,otume)))
(2*length(intersect(otume,otulo)))/(length(intersect(otume,otulo))+length(union(otume,otulo)))




betameanse <- comm.data %>% group_by(lomehi) %>% select(year,class_3,anpp) %>% summarise_each(funs(mean,sd,std.error))

#from prod code: how to use dplyr

betameanse <- prod %>% filter(year!=1997) %>% group_by(year,class_3) %>% select(year,class_3,anpp) %>% summarise_each(funs(mean,sd,std.error))  

ggplot(prodmean, aes(x=year, y=anpp, color=class_3)) + geom_line()







######Ordination######
#dats10 is the phyloseq object, dats10otu is the data.frame, dats10r is the rarefied phyloseq object, doubletons and singletons removed and sequences <.2% rel abun removed
#ordination to compare with html file
ordu = ordinate(talus,method="PCoA", distance="unifrac", weighted = TRUE,normalized=FALSE)
p=plot_ordination(talus, ordu, color = "Depth", shape = "Description")
p = p + geom_point(size = 7, alpha = 0.75)
p
distmatp<-UniFrac(talus,weighted=TRUE)
distmatp2<-UniFrac(talus,weighted=TRUE,normalized=FALSE)#this produces the output from the qiime code, the axis percent varience explained are slightly different from qiime, probably b/c how they deal with negative eigenvalues. here it is eig/total


#Ordination with bray curtis
sample_data(dats10r)$Description<-factor(lomehi,levels=c("lo","me","hi"))
sample_data(dats10r)$Description<-factor(ifelse(sample_data(dats10r)$pH>5,"hi","lo"))
sample_data(dats10r)$Description<-factor(ifelse(sample_data(dats10r)$moisture>10,"hi","lo"))

ordu = ordinate(dats10r,method="PCoA", distance="bray") #weighed and normalized are passed to unifrac, so they don't matter if you use bray
plot_ordination(dats10r, ordu, color = "Description") +#, shape = "Description"
  geom_point(size = 5, alpha = 1)#alpha is the opaqueness of the point, can be set to continuous numbers

#see http://joey711.github.io/phyloseq/plot_ordination-examples for more examples on how to make graph prettier and integrate better with ggplot

range01 <- function(x){(x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))}











