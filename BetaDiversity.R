#Bray-curtis BetaDiversity

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

