


######Bray-curtis BetaDiversity######

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











