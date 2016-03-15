
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Cleaning.Rdata")
load("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Cleaning.Rdata")

#######Read in OTU data#######
library(phyloseq)
packageVersion("phyloseq")
#library(ape)
#library(ggplot2)

#I tried to use this, but it doesn't work with Euks because the number of levels are all different
#parse_taxonomy_silva<-function (char.vec) {  parse_taxonomy_greengenes(strsplit(char.vec, "_", TRUE)[[1]]) }
#otufile <-import_biom("Euk_ALL_97_OTU_filtsing.biom",parseFunction = parse_taxonomy_silva)

#Read in biom file, with singletons removed
otufile <-import_biom("Euk_ALL_97_OTU_filtsing.biom")
tax_table(otufile)
head(tax_table(otufile))

#Read in a file that lists taxonomic levels in the silva database, https://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/tax_slv_ssu_nr_119.txt
silvamap<-read.delim("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/elitalus/tax_slv_ssu_nr_119.txt",header=F)
head(silvamap)

#Assign taxonomic levels to the tax file. 
#First get tax_table in the write format with ; instead of _ between tax levels
tax_table2<-as.matrix(tax_table(otufile))
#tax_table2<-tax_table2[1:5] #trial on first five rows
dim(tax_table2)
tax_table3<-as.factor(tax_table2[,1])
levels(tax_table3)
tax_table4<-as.factor(gsub("_",";",tax_table3))
tax_table5<-as.factor(paste(tax_table4,rep(";",length(tax_table4)),sep=""))
#ind<-match(tax_table5,silvamap$V1)
#silvamap$V3[ind]

#Then make a table with the names to search against the silva map file
head(silvamap)
head(tax_table5)
tax_table6<-matrix(0,nrow=length(tax_table5),ncol=15)
for (i in 1:length(tax_table5)){
  for (j in 1:length(as.vector(gregexpr(";",tax_table5[i])[[1]]))){
    second<-as.vector(gregexpr(";",tax_table5[i])[[1]])[j]
    tax_table6[i,j]<-substr(tax_table5[i],1,second)
  }
}

#Then search against the silva map file# this takes a while, maybe 8 min
silvamapmatrix<-matrix(0,nrow=length(tax_table5),ncol=15)
for (i in 1:nrow(tax_table6)){
  for (j in 1:ncol(tax_table6)){
    ind<-match(tax_table6[i,j],silvamap$V1)
    silvamapmatrix[i,j]<-as.character(silvamap$V3[ind])
  }
}

#Make file with just last name from the table6 file
tax_table7<-matrix(0,nrow=length(tax_table5),ncol=15)
for (i in 1:length(tax_table5)){
  for (j in 1:length(as.vector(gregexpr(";",tax_table5[i])[[1]]))){
    if (j==1){
      first<-0
      second<-as.vector(gregexpr(";",tax_table5[i])[[1]])[j]
    } else {   
      first<-as.vector(gregexpr(";",tax_table5[i])[[1]])[j-1]
      second<-as.vector(gregexpr(";",tax_table5[i])[[1]])[j]
    }
    tax_table7[i,j]<-substr(tax_table5[i],first+1,second-1)
  }
}

#Then order them correctly so that each level is put in the correct column
mytaxlevels<-c("domain","major_clade","superkingdom","kingdom","subkingdom","superphylum","phylum","subphylum","class","subclass","superorder","order","family","subfamily","genus")

mytaxmatrix<-matrix(0,nrow=length(tax_table5),ncol=15)
for(i in 1:nrow(tax_table7)){
  ind<-match(silvamapmatrix[i,],mytaxlevels)
  len<-length(as.vector(gregexpr(";",tax_table5[i])[[1]]))
  ind2<-which(is.na(ind)==F)
  ind3<-ind[ind2]
  mytaxmatrix[i,ind3]<-tax_table7[i,][ind2]
}
colnames(mytaxmatrix)<-mytaxlevels
rownames(mytaxmatrix)<-rownames(tax_table2)
head(mytaxmatrix)
mytaxmatrix[which(mytaxmatrix=="0")]<-NA
#mytaxtable<-as.data.frame(mytaxmatrix)
#class(mytaxmatrix)<-class(tax_table2)
#head(mytaxmatrix)
#mytaxmatrix[10:50,]


#I should do some aggregating for orders and kingdom levels and integrate into mytaxmatrix here. however I only did it below for the 2015 data so I would need to look more closely at it if I wanted to insert it here.



#Replace the new tax matrix in the otu file
tax_table(otufile)<-mytaxmatrix

#Import mapping and tree file
map<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/EukBr_Niwot_20072015_All_MapFile.txt")

dat<-merge_phyloseq(otufile,map)


XXXXXXXXX
treefile <- read_tree("talus_euks_sinaaln_rep_set.tre")
plot(treefile)

dat<-merge_phyloseq(otufile,map,treefile)
  
talus <-merge_phyloseq(otufile,map,treefile)
#taxtable<-as.data.frame(run1@tax_table)
#taxtable2<-taxtable %>% separate(Rank1, c("r1", "r2","r3","r4","r5","r6","r7","r8","r9","r10","r11","r13"), "_", extra="merge")


myTaxa = names(sort(taxa_sums(talus), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, talus)
plot(phy_tree(ex1), show.node.label = TRUE)
plot_tree(ex1, color = "Description", label.tips = "phylum", ladderize = "left", size = "Abundance",plot.margin=.8)#,justify = "left"


#Remove OTUs that do not appear more than 5 times in more than half the samples
wh0 = genefilter_sample(talus, filterfun_sample(function(x) x > 5), A = 0.5 * nsamples(talus))
talus1 = prune_taxa(wh0, talus)
plot_tree(ex1, color = "Description", label.tips = "kingdom", ladderize = "left", justify = "left" , size = "Abundance")

#ordination to compare with html file
ordu = ordinate(talus,method="PCoA", distance="unifrac", weighted = TRUE,normalized=FALSE)
p=plot_ordination(talus, ordu, color = "Depth", shape = "Description")
p = p + geom_point(size = 7, alpha = 0.75)
p
distmatp<-UniFrac(talus,weighted=TRUE)
distmatp2<-UniFrac(talus,weighted=TRUE,normalized=FALSE)#this produces the output from the qiime code, the axis percent varience explained are slightly different from qiime, probably b/c how they deal with negative eigenvalues. here it is eig/total










######Soil data from 2015######
dats<-subset_samples(dat, SampleType=="soil"&year=="2015")

#take plants out
dats2<-subset_taxa(dats,domain=="Eukaryota"&class!="Embryophyta")

#Renaming and organizing orders and kingdoms here, then replace the tax table in dats2 to include these new groupings
dats2tax<-as.data.frame(tax_table(dats2))
dim(dats2tax)

#aggregate things by order - first for orders that are NA/uncultured/incertaesedis use class/phylum information in the order colum
sort(unique(as.character(dats2tax$order)))
orders<-as.character(dats2tax$order)
classes<-as.character(dats2tax$class)
phyla<-as.character(dats2tax$phylum)
kingdom<-as.character(dats2tax$kingdom)
indo<-which(orders%in%c("Incertae Sedis","uncultured",NA))
indc<-which(classes%in%c("Incertae Sedis","uncultured",NA))
indp<-which(phyla%in%c("Incertae Sedis","uncultured",NA))
indk<-which(kingdom%in%c("Incertae Sedis","uncultured",NA))
intersect(indp,indk)#every sample either has a phylum or kingdom
orders[indo]<-NA
classes[indc]<-NA
phyla[indp]<-NA
kingdom[indk]<-NA
ind1<-which(is.na(orders)==T&is.na(classes)==F)
ind2<-which(is.na(orders)==T&is.na(classes)==T&is.na(phyla)==F)
#ind3<-which(is.na(orders)==T&is.na(classes)==T&is.na(phyla)==T) #not needed, the things without phylum must have had a class/order
orders[ind1]<-paste("Unclassified",as.character(dats2tax$class[ind1]))
orders[ind2]<-paste("Unclassified",as.character(dats2tax$phylum[ind2]))
#orders[ind3]<-paste("Unclassified",as.character(dats4tax$kingdom[ind3]))

#Making groups for labeling graphs
#Fungi (kingdom), 
#Archaeplastida (major_clade), combine the kingdoms Chloroplastida, Stramenopiles, Rhodophyceae, Haptophyta, 
#Rhizaria (kingdom: unicellular amoeboid euks), 
#Amoebozoa(kingdom),
#Holozoa(kingdom: all animals, not fungi) - note within the animals the silva tax map file was blank for a numer of the important groups like nematodes/tardigrades/arthropods. I can either leave as is and call all "animals" or go to the silva file and annotate by hand
#Discicristoidea(kingdom: amoebas),
#photosynthetic Alveolata (phylum Dinoflagellata: mostly photosynthetic; nonphotosynthetic Protoperidinium, both SL163A10 and Pfiesteria (can if it eats an alga), unknown D244), 
#nonphotosynthetic Alveolata (phyla Ciliophora(predators), protalveolata, apicomplexa (parasites)), 
#photosynthetic Discoba (phylum Euglenozoa: mostly photosynthetic), 
#nonphotosnthetic Discoba (phylum Heterolobosea: parasites, free living, symbiotic, amoeba-like), 
#NA - are all heterotrophic protists/parasites

head(orders)
kingdomlabels<-kingdom
ind<-which(kingdomlabels%in%c("Chloroplastida","Stramenopiles","Rhodophyceae","Haptophyta"))
kingdomlabels[ind]<-"Archaeplastida"
ind<-which(kingdomlabels=="Alveolata"&phyla=="Dinoflagellata")
kingdomlabels[ind]<-"Photosynthetic_Alveolata"
ind<-which(kingdomlabels=="Alveolata")
kingdomlabels[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(kingdomlabels=="Discoba"&phyla=="Euglenozoa")
kingdomlabels[ind]<-"Photosynthetic_Discoba"
ind<-which(kingdomlabels=="Discoba"&phyla=="Heterolobosea")
kingdomlabels[ind]<-"Nonphotosynthetic_Discoba"
ind<-which(is.na(kingdomlabels))# a few are NA, they are both heterotrophs or parasites/symbiont phyla[which(is.na(kingdomlabels))]
kingdomlabels[ind]<-"Nonphotosynthetic_Eukaryota"

labelfile<-unique(as.data.frame(cbind(otu=rownames(dats2tax),orders,kingdomlabels)))
dim(labelfile)
head(labelfile)

head(dats2tax)
dats2tax$ordergroup<-orders
dats2tax$kingdomgroup<-kingdomlabels
dats2taxm<-as.matrix(dats2tax)
head(dats2taxm)
head(tax_table(dats2))
tax_table(dats2)<-dats2taxm



#remove doubletons and singletons, taxa present in only 1 or 2 samples
wh0 = genefilter_sample(dats2, filterfun_sample(function(x) x > 0), A = 3)
dats3 = prune_taxa(wh0, dats2)

#transform to relative abudance
#dats3 = transform_sample_counts(dats2, function(x) 100*x/sum(x) )

#aggregate by order, this either removes NAs or groups NAs into one group, it doesn't look at the next highest level of taxonomy
#dats5<-tax_glom(dats4,"order",NArm=FALSE)

unique(tax_table(dats3)[,"major_clade"])
unique(tax_table(dats2)[,"kingdom"])
unique(tax_table(subset_taxa(dats3,kingdom=="Discoba")))

#dat5<-psmelt(dats4)#this puts it in long format (abundance is a column)

dats3tax<-as.data.frame(tax_table(dats3))
dats3sample<-as.data.frame(sample_data(dats3))
dats3otu<-as.data.frame(otu_table(dats3))

#aggregate things by order - first for orders that are NA/uncultured/incertaesedis use class/phylum information in the order colum
sort(unique(as.character(dats3tax$order)))
orders<-as.character(dats3tax$order)
classes<-as.character(dats3tax$class)
phyla<-as.character(dats3tax$phylum)
kingdom<-as.character(dats3tax$kingdom)
indo<-which(orders%in%c("Incertae Sedis","uncultured",NA))
indc<-which(classes%in%c("Incertae Sedis","uncultured",NA))
indp<-which(phyla%in%c("Incertae Sedis","uncultured",NA))
indk<-which(kingdom%in%c("Incertae Sedis","uncultured",NA))
intersect(indp,indk)#every sample either has a phylum or kingdom
orders[indo]<-NA
classes[indc]<-NA
phyla[indp]<-NA
kingdom[indk]<-NA
ind1<-which(is.na(orders)==T&is.na(classes)==F)
ind2<-which(is.na(orders)==T&is.na(classes)==T&is.na(phyla)==F)
#ind3<-which(is.na(orders)==T&is.na(classes)==T&is.na(phyla)==T) #not needed, the things without phylum must have had a class/order
orders[ind1]<-paste("Unclassified",as.character(dats3tax$class[ind1]))
orders[ind2]<-paste("Unclassified",as.character(dats3tax$phylum[ind2]))
#orders[ind3]<-paste("Unclassified",as.character(dats4tax$kingdom[ind3]))

head(orders)
kingdomlabels<-kingdom
ind<-which(kingdomlabels%in%c("Chloroplastida","Stramenopiles","Rhodophyceae","Haptophyta"))
kingdomlabels[ind]<-"Archaeplastida"
ind<-which(kingdomlabels=="Alveolata"&phyla=="Dinoflagellata")
kingdomlabels[ind]<-"Photosynthetic_Alveolata"
ind<-which(kingdomlabels=="Alveolata")
kingdomlabels[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(kingdomlabels=="Discoba"&phyla=="Euglenozoa")
kingdomlabels[ind]<-"Photosynthetic_Discoba"
ind<-which(kingdomlabels=="Discoba"&phyla=="Heterolobosea")
kingdomlabels[ind]<-"Nonphotosynthetic_Discoba"
ind<-which(is.na(kingdomlabels))# two are NA, they are both heterotrophs phyla[which(is.na(kingdomlabels))]
kingdomlabels[ind]<-"Nonphotosynthetic_Eukaryota"

labelfile<-as.data.frame(cbind(orders,kingdomlabels))

head(dats3tax)
dats3tax$ordergroup<-orders
XXXXXXX

head(dats3otu)
head(dats3sample)

dats4<-cbind(dats3tax,dats3otu)
dats5<-aggregate.data.frame(dats4[,17:113],by=list(ordergroup=dats4$ordergroup),sum)
rownames(dats5)<-dats5$ordergroup
dats5$ordergroup<-NULL
dats6<-t(dats5)
head(dats6)

#the minimum number of reads is 1250, which is pretty high, so right now I'm not rarefying
min(rowSums(dats6))

dats6order<-cbind(dats3sample,dats6)
head(dats6order)
greater66plants<-factor(ifelse(dats6order$Plant_Dens>66,"hi","lo")) #this is stem density, including mosses
dats6order<-cbind(greater66plants,dats6order)
dats6order$Sample_name<-as.numeric(as.character(dats6order$Sample_name))


#Read in plant data to merge with order file
plantcomp<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots
dim(plantcomp)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantlabels<-as.data.frame(cbind(colnames(plantcomp2)[2:56],"Plant"))
colnames(plantlabels)<-c("orders","kingdomlabels")

#Merge plants with microbes, plantcomp is everything, plantcomp2 removes doubletons/singletons
microbplant<-merge(dats6order,plantcomp,"Sample_name",sort=F)
microbplant2<-merge(dats6order,plantcomp2,"Sample_name",sort=F)
head(microbplant)






######Grouping by kingdom#####
dat7<-data.frame(dats4[,1:16],kingdomlabels,dats4[,17:113])

dat8<-aggregate.data.frame(dat7[,18:114],by=list(kingdomlabels=dat7$kingdomlabels),sum)
rownames(dat8)<-dat8$kingdomlabels
dat8$kingdomlabels<-NULL
dat9<-t(dat8)
head(dat9)
dat9kingdom<-cbind(greater66plants,dats3sample,dat9)
species<-dat9kingdom[,28:38]
speciesrel<-species/rowSums(species)*100
m1<-aggregate.data.frame(speciesrel, by=list(greater66plants),mean)







#####Group things by genus / otu#####
head(dats3otu)
dats7otu<-t(dats3otu)
dats8otu<-cbind(dats3sample,dats7otu)
dim(dats8otu)
head(dats8otu)[1:4,1:28]#otus start at column 27
#Remove otus only present in five or fewer plots
dats8otuspe<-dats8otu[,27:5042]
dats9otu<-cbind(dats8otu[,1:26],dats8otuspe[,colSums(dats8otuspe>0)>5])#2455 taxa
head(dats9otu)[,1:35]

#second try, pruning a bit more, need to be in 10 samples
wh0 = genefilter_sample(dats2, filterfun_sample(function(x) x > 0), A = 15)#1008 taxa
dats10 = prune_taxa(wh0, dats2)
dats10
min(sample_sums(dats10))
dats10otu<-cbind(sample_data(dats2),t(otu_table(dats10)))
head(dats10otu)[,1:40]