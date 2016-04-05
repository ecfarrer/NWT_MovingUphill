#######Read in OTU data#######
library(phyloseq)
library(foreach)
library(doParallel)
packageVersion("phyloseq")

#Read in biom file, with singletons removed
otufile16S <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/16S_ALL_97_OTU_filtsing.biom")
tax_table(otufile16S)
head(tax_table(otufile16S))

#Read in a file that lists taxonomic levels in the silva database, https://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/tax_slv_ssu_nr_119.txt
silvamap<-read.delim("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/elitalus/tax_slv_ssu_nr_119.txt",header=F)
head(silvamap)

#silvamap has Bacteria;Actinobacteria; listed as a class, which is incorrect, it should be a phylum. changing this:
silvamap$V3[which(silvamap$V1=="Bacteria;Actinobacteria;")]<-"phylum"
#it also has "uncultured" listed as a class incorrectly many times, I should put is as a genus (this is only a problem for bacteria not eukarya). but I have to fix this later on

#Assign taxonomic levels to the tax file. 
#First get tax_table in the write format with ; instead of _ between tax levels
tax_table16S2<-as.matrix(tax_table(otufile16S))
#tax_table2<-tax_table2[1:5] #trial on first five rows
dim(tax_table16S2)
tax_table16S3<-as.factor(tax_table16S2[,1])
levels(tax_table16S3)
tax_table16S4<-as.factor(gsub("_",";",tax_table16S3))
tax_table16S5<-as.factor(paste(tax_table16S4,rep(";",length(tax_table16S4)),sep=""))

#Then make a table with the names to search against the silva map file
head(silvamap)
head(tax_table16S5)
tax_table16S6<-matrix(0,nrow=length(tax_table16S5),ncol=11)#was 15 for euks
for (i in 1:length(tax_table16S5)){
  for (j in 1:length(as.vector(gregexpr(";",tax_table16S5[i])[[1]]))){
    second<-as.vector(gregexpr(";",tax_table16S5[i])[[1]])[j]
    tax_table16S6[i,j]<-substr(tax_table16S5[i],1,second)
  }
}
unique(tax_table16S6[,1]) #it only has a max of 11 columns
#the taxonomy is messed up, within one taxonmy line there are multiple items that match to "class" for example, so my original code will not work. Thus I will try to just pull out kingdom, order, and otu

#Then search against the silva map file# this takes a while, maybe 2 min
silvamapmatrix16S<-matrix(0,nrow=length(tax_table16S5),ncol=11)#was 15 for euks
for (i in 1:nrow(tax_table16S6)){
  ind<-match(tax_table16S6[i,],silvamap$V1)
  silvamapmatrix16S[i,]<-as.character(silvamap$V3[ind])
}

#only 5 instances of order being duplicated, and my code would actually treat them correctly
#below I fixed the class duplication (whcih was thousands)
#there are 5families that are wrong, and my code would do it incorrectly
ind<-which(rowSums(silvamapmatrix16S=="class",na.rm=T)>1)
silvamapmatrix16S[ind[2],]
tax_table16S6[ind[2],]

#Make file with just last name from the table6 file
tax_table16S7<-matrix(0,nrow=length(tax_table16S5),ncol=15)
for (i in 1:length(tax_table16S5)){
  for (j in 1:length(as.vector(gregexpr(";",tax_table16S5[i])[[1]]))){
    if (j==1){
      first<-0
      second<-as.vector(gregexpr(";",tax_table16S5[i])[[1]])[j]
    } else {   
      first<-as.vector(gregexpr(";",tax_table16S5[i])[[1]])[j-1]
      second<-as.vector(gregexpr(";",tax_table16S5[i])[[1]])[j]
    }
    tax_table16S7[i,j]<-substr(tax_table16S5[i],first+1,second-1)
  }
}

head(tax_table16S7)

#fix the places where uncultured is a class, change it to genus
ind<-which(tax_table16S7[,1]=="Bacteria"&rowSums(tax_table16S7=="uncultured")==1)
ind
for(i in 1:length(ind)){
  ind2<-which(tax_table16S7[ind[i],]=="uncultured")
  silvamapmatrix16S[ind[i],ind2]<-"genus"
}
head(silvamapmatrix16S)


#Then order them correctly so that each level is put in the correct column
mytaxlevels<-c("domain","major_clade","superkingdom","kingdom","subkingdom","superphylum","phylum","subphylum","class","subclass","superorder","order","family","subfamily","genus")

mytaxmatrix16S<-matrix(0,nrow=length(tax_table16S5),ncol=15)
for(i in 1:nrow(tax_table16S7)){
  ind<-match(silvamapmatrix16S[i,],mytaxlevels)
  len<-length(as.vector(gregexpr(";",tax_table16S5[i])[[1]]))
  ind2<-which(is.na(ind)==F)
  ind3<-ind[ind2]
  mytaxmatrix16S[i,ind3]<-tax_table16S7[i,][ind2]
}
colnames(mytaxmatrix16S)<-mytaxlevels
rownames(mytaxmatrix16S)<-rownames(tax_table16S2)
head(mytaxmatrix16S)
mytaxmatrix16S[which(mytaxmatrix16S=="0")]<-NA
#mytaxtable<-as.data.frame(mytaxmatrix)
#class(mytaxmatrix)<-class(tax_table2)
#head(mytaxmatrix)
#mytaxmatrix[10:50,]



#Replace the new tax matrix in the otu file
tax_table(otufile16S)<-mytaxmatrix16S

#Import mapping and tree file
map16S<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/515BC_Niwot_20072015_All_MapFile.txt")

dat16S<-merge_phyloseq(otufile16S,map16S)


#XXXXXXXXX add tree file when fastree stops
#This has not yet been integrated into dat
treefile <- read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euk_ALL_truncate_97_sinaaln_rep_set.tre")
plot(treefile) #takes a long time

dattree<-merge_phyloseq(otufile,map,treefile) #takes a long time, maybe 30 min to 1 hr







######Soil data from 2015######
dat16Ss<-subset_samples(dat16S, SampleType=="soil"&year=="2015")

#take euks out
dat16Ss2<-subset_taxa(dat16Ss,domain%in%c("Archaea","Bacteria")&family!="mitochondria"&class!="Chloroplast")

#Renaming and organizing orders and kingdoms here, then replace the tax table in dats2 to include these new groupings
dat16Ss2tax<-as.data.frame(tax_table(dat16Ss2))
dim(dat16Ss2tax)


#all samples have a domain and a phylum and a class (although one level is "unknown class"). no kingdoms and there are NAs for order 
#I'm going to skip the labels for order right now
sort(unique(dat16Ss2tax$phylum))
#photosynthetic phyla are Chloroflexi and Cyanobacteria
photbac<-ifelse(tax_table(dat16Ss2)[,"phylum"]%in%c("Chloroflexi","Cyanobacteria"),"PhotosyntheticBacteria","NonphotosyntheticBacteria")
bacterialabels<-data.frame(cbind(otu=rownames(tax_table(dat16Ss2)),orders=tax_table(dat16Ss2)[,"phylum"],kingdomlables=photbac));colnames(bacterialabels)<-c("otu","orders","kingdomlabels")
head(bacterialabels)


#Then remove more rare taxa which will create a datset for network analysis

#Beginning of code for order groups
#Remove doubletons and singletons, taxa present in only 1 or 2 samples
wh0 = genefilter_sample(dat16Ss2, filterfun_sample(function(x) x > 0), A = 3)
dat16Ss3 = prune_taxa(wh0, dat16Ss2)



#####Group things by genus / otu#####
#Third try, trying to be consistent with sample size minimum in the network analysis (my cutoff there is 3 or 4 or 5). Follwing Widder et al 2014 PNAS
wh0 = genefilter_sample(dat16Ss2, filterfun_sample(function(x) x > 0), A = 3)#30635 taxa
dat16Ss10a = prune_taxa(wh0, dat16Ss2)
dat16Ss10a
#transform to relative abundance
dat16Ss10b = transform_sample_counts(dat16Ss10a, function(x) x/sum(x) )
wh0 = filter_taxa(dat16Ss10b, function(x) sum(x) > (0.002), prune=F)
dat16Ss10c = prune_taxa(wh0, dat16Ss10a)  #3019
dat16Ss10c
sort(sample_sums(dat16Ss10c)) #need to take out samples 34, 5, 126 because they had less than 30 total reads and need to rarefy the rest to 4896 reads
dat16Ss10 = prune_samples(names(which(sample_sums(dat16Ss10c) >= 200)), dat16Ss10c)

dat16Ss10otu<-cbind(sample_data(dat16Ss10),t(otu_table(dat16Ss10)))
head(dat16Ss10otu)[,1:40]
dim(dat16Ss10otu)

dat16Ss10r<-rarefy_even_depth(dat16Ss10,sample.size=min(sample_sums(dat16Ss10)),rngseed=10,replace=F)#rarefy dats10 for betadiversity metrics, 12 OTUS were removed because thye were no longer present in any sample after rarefaction
#nice, this is the same exact rarefaction as the comm.data
#head(comm.data)[,28:50]
#head(t(otu_table(dats10r)))[,1:30]









