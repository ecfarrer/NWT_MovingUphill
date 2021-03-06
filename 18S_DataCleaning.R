
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/")

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

#checking for duplicate levels of taxonmy
#none for phylum, kingdom, domain
#for class, it is fine, all dups are bacteria
#for order there are 1065, most are plants which I will ignore, 252 are from euks. most are from Collophora which has 'Incertae Sedis' duplicated for "order", which is fine. 1:7 are "Incertae Sedis" "Zoopagales" which is in the correct order. 9:10 are 'Incertae Sedis' duplicated. 10:19 are "Incertae Sedis" "Zoopagales" which is fine. so all are fine
ind<-which(rowSums(silvamapmatrix=="domain",na.rm=T)>1);ind
ind<-which((rowSums(silvamapmatrix=="order",na.rm=T)>1)&(rowSums(tax_table7=="Embryophyta",na.rm=T)==0)&(rowSums(tax_table7=="Collophora",na.rm=T)==0));ind
silvamapmatrix[ind,]
tax_table7[ind,]

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
} #this is slighly incorrect. I just realized that sometimes there is a row with two items that are called "order". this loop will put the second name in the final mytaxmatrix. I checked above and it looks like this doesnt matter, all the times order is duplicated, it is two "Incertae Sedis", or Incertae Sedis then the correct order.

colnames(mytaxmatrix)<-mytaxlevels
rownames(mytaxmatrix)<-rownames(tax_table2)
head(mytaxmatrix)
mytaxmatrix[which(mytaxmatrix=="0")]<-NA
#mytaxtable<-as.data.frame(mytaxmatrix)
#class(mytaxmatrix)<-class(tax_table2)
#head(mytaxmatrix)
#mytaxmatrix[10:50,]


#I should do some aggregating for orders and kingdom levels and integrate into mytaxmatrix here. however I only did it below for the 2015 data so I would need to look more closely at it if I wanted to insert it here. it actually might not be useful to do the aggregating here since there are still bacteria in the dataset



#Replace the new tax matrix in the otu file
tax_table(otufile)<-mytaxmatrix

#Import mapping and tree file
map<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/EukBr_Niwot_20072015_All_MapFile.txt")

dat<-merge_phyloseq(otufile,map)




#Add tree file
#This has not yet been integrated into dat
treefile <- read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_truncate_97_sinaaln_rep_set.tre")
plot(treefile) #takes a long time

dattree<-merge_phyloseq(otufile,map,treefile) #takes a long time, maybe 30 min to 1 hr. I lose 23,855 taxa when I merge the treefile, I assume this is because some samples could not be aligned


myTaxa = names(sort(taxa_sums(dattree), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, dattree)
plot(phy_tree(ex1), show.node.label = TRUE)
plot_tree(ex1, color = "Description", label.tips = "phylum", ladderize = "left", size = "Abundance",plot.margin=.8)#,justify = "left"








######Soil data from 2015######
dats<-subset_samples(dat, SampleType=="soil"&year=="2015")
datst<-subset_samples(dattree, SampleType=="soil"&year=="2015")

#take plants out
dats2<-subset_taxa(dats,domain=="Eukaryota"&class!="Embryophyta")
dats2t<-subset_taxa(datst,domain=="Eukaryota"&class!="Embryophyta")
dats2tr<-rarefy_even_depth(dats2t,sample.size=min(sample_sums(dats2t)),rngseed=10,replace=F)#rarefy to 761, this is quite low
pdEuk<-pd(as.matrix(as.data.frame(t(otu_table(dats2tr)))), phy_tree(dats2tr), include.root=TRUE) #takes a few minutes

length(which(colSums(t(otu_table(dats2)))>0))#remove zeros, count total number of OTUs
sum(t(otu_table(dats2)))

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

dats2otu<-cbind(sample_data(dats2),t(otu_table(dats2)))
head(dats2otu)[,1:10]

#dats2 is the final data file for diversity analyses. the only thing removed are plants and otus with only one read (which is likely error)





#Then remove more rare taxa which will create a datset for network analysis
#Remove doubletons and singletons, taxa present in only 1 or 2 samples
wh0 = genefilter_sample(dats2, filterfun_sample(function(x) x > 0), A = 3)
dats3 = prune_taxa(wh0, dats2)

#Transform to relative abudance
#dats3 = transform_sample_counts(dats2, function(x) 100*x/sum(x) )

#One way to aggregate by order, this either removes NAs or groups NAs into one group, it doesn't look at the next highest level of taxonomy
#dats5<-tax_glom(dats4,"order",NArm=FALSE)

unique(tax_table(dats3)[,"major_clade"])
unique(tax_table(dats2)[,"kingdom"])
unique(tax_table(subset_taxa(dats3,kingdom=="Discoba")))

#dat5<-psmelt(dats4)#this puts it in long format (abundance is a column)

#To create a datset at the order level
dats3tax<-as.data.frame(tax_table(dats3))
dats3sample<-as.data.frame(sample_data(dats3))
dats3otu<-as.data.frame(otu_table(dats3))

head(dats3otu)
head(dats3sample)

dats4<-cbind(dats3tax,dats3otu)
dats5<-aggregate.data.frame(dats4[,18:114],by=list(ordergroup=dats4$ordergroup),sum)
rownames(dats5)<-dats5$ordergroup
dats5$ordergroup<-NULL
dats6<-t(dats5)
head(dats6)

#The minimum number of reads is 1250, which is pretty high, so right now I'm not rarefying
min(rowSums(dats6))

dats6order<-cbind(dats3sample,dats6)
head(dats6order)
dats6order$Sample_name<-as.numeric(as.character(dats6order$Sample_name))








#####Group things by genus / otu#####
#First try - Remove otus only present in five or fewer plots
wh0 = genefilter_sample(dats2, filterfun_sample(function(x) x > 0), A = 6)#2455 taxa
dats10 = prune_taxa(wh0, dats2)
dats10
min(sample_sums(dats10)) #1191 reads
dats10otu<-cbind(sample_data(dats2),t(otu_table(dats10)))
head(dats10otu)[,1:40]

#Second try, pruning a bit more, need to be in 15 samples
wh0 = genefilter_sample(dats2, filterfun_sample(function(x) x > 0), A = 15)#1008 taxa
dats10 = prune_taxa(wh0, dats2)
dats10
min(sample_sums(dats10)) #need to rarefy to 1121 reads
dats10otu<-cbind(sample_data(dats2),t(otu_table(dats10)))
head(dats10otu)[,1:40]

#Third try, filter doubletons/singletons. Follwing Widder et al 2014 PNAS
wh0 = genefilter_sample(dats2, filterfun_sample(function(x) x > 0), A = 3)#5016 taxa
dats10a = prune_taxa(wh0, dats2)
dats10a
#transform to relative abundance
dats10b = transform_sample_counts(dats10a, function(x) x/sum(x) )
wh0 = filter_taxa(dats10b, function(x) sum(x) > (0.002), prune=F)
dats10 = prune_taxa(wh0, dats10a)  #2182
dats10
min(sample_sums(dats10)) #need to rarefy to 1236 reads
dats10otu<-cbind(sample_data(dats2),t(otu_table(dats10)))
head(dats10otu)[,1:40]
dats10r<-rarefy_even_depth(dats10,sample.size=min(sample_sums(dats10)),rngseed=10,replace=F)#rarefy dats10 for betadiversity metrics, 12 OTUS were removed because thye were no longer present in any sample after rarefaction
#nice, this is the same exact rarefaction as the comm.data
head(comm.data)[,28:50]
head(t(otu_table(dats10r)))[,1:30]






#I can't extract nematode samples OTUs because the holozoa do not have labels in the mapping file. tax_table2 is the original tax table with full line of taxononmy Eukaryota_Opisthokonta_Nucletmycea_Fungi_Dikarya_Ascomycota_Pezizomycotina_Sordariomycetes_Hypocreales
#all taxa in the file going into the network analysis (not rarefied but doesn't matter)
temp3<-rownames(tax_table(dats10)) 

#full taxonomies
head(tax_table2)
temptax<-as.character(tax_table2)
ind<-grep("Nematoda",temptax)
temp2<-tax_table2[ind,]
nematodeids<-rownames(temp2)

intersect(temp,nematodeids) #there are no nematodes present in the otufile going into network analysis
length(union(temp,nematodeids))

#checking that there are still nematodes in the original 2015 subset file, they are present in dat and dats but not dats2
temp<-rownames(otu_table(dats))
intersect(temp,nematodeids)

temp <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euk_ALL_97_OTU_filtsing.biom")
dim(tax_table(temp))
dim(otu_table(temp))

ind<-which(rownames(otu_table(temp))=="denovo244871")
tax_table(temp)[ind,] #it imported correctly and has holozoa in the label

#after my taxonomy rearrangement what did it do? that is correct
tax_table(otufile)[ind,]

#











######Grouping by kingdom#####
dats7<-transform_sample_counts(dats2, function(x) 100*x/sum(x) )
dats8<-otu_table(dats7)
dats9<-aggregate.data.frame(dats8,by=list(kingdomlabels=kingdomlabels),sum)
rownames(dats9)<-dats9$kingdomlabels
dats9$kingdomlabels<-NULL
dats9kingdom<-cbind(sample_data(dats2),t(dats9))







##### Nematode data #####
#import_biom gave an error when I tried to import the biom file I created on the server from Dorota's otu text table. I can't figure out why it gave the error beacuse the text file looks exactly like the other text files that have all euk data.
otufileN<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Nematodes/Nema_OTUtable_species.csv")
head(otufileN)
otufileN1<-otufileN[,2:(dim(otufileN)[2]-1)]
rownames(otufileN1)<-otufileN[,1]
head(otufileN1)
otufileN2 = otu_table(otufileN1, taxa_are_rows = TRUE)
taxfileN<-matrix(otufileN[,dim(otufileN)[2]])
rownames(taxfileN)<-otufileN[,1]
colnames(taxfileN)<-"Rank1"
taxfileN1 = tax_table(taxfileN)
datN = phyloseq(otufileN2,taxfileN1,map,treefile)
datN
pdN<-pd(as.matrix(as.data.frame(t(otu_table(datN)))),phy_tree(datN),include.root=TRUE) 

rownames(otu_table(otufile))[which(rownames(otu_table(otufile))%in%rownames(otufileN1))]
treefile$tip.label[which(treefile$tip.label%in%rownames(otufileN1))]
denovo120212 denovo68570
tax_table(otufile)[which(rownames(tax_table(otufile))=="denovo256418")]
#denovo256418

















#Haven't done this yet. I did this in the cooccurrencenetworkseuks file
######Read in plant data to merge with order file######
#plantcomp<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_comp2015.csv")
#head(plantcomp)
#names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots
#dim(plantcomp)
#plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
#plantlabels<-as.data.frame(cbind(colnames(plantcomp2)[2:56],"Plant"))
#colnames(plantlabels)<-c("orders","kingdomlabels")

#Merge plants with microbes, plantcomp is everything, plantcomp2 removes doubletons/singletons
#microbplant<-merge(dats6order,plantcomp,"Sample_name",sort=F)
#microbplant2<-merge(dats6order,plantcomp2,"Sample_name",sort=F)
#head(microbplant)









