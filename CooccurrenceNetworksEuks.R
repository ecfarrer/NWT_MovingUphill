setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")
#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis.Rdata")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata") #results has the results from OTU cooccurrence patterns

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU.Rdata")
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_AnalysisOrder.Rdata")


library(vegan)
library(reshape)
library(plotrix)
library(foreach)
library(doParallel)
library(Kendall)


######Co-occurrence pairwise correlations###### 
#Input dataset is dats6order or dats10otu, dats10otu has doubletons and singletons removed

comm.data<-dats6order
comm.data<-dats10otu
  
#Standardize to 100 or rarified to 1250 (order) or 1121 (otu) or 1236 (otu with cutoff of 3 rather than 15). if I standardize, I will need to change code in the loop below because it has cutoffs of 1
set.seed(10)
#comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1250))#by order
#comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1121))#by otu
comm.data<-cbind(comm.data[,c(1:26)],rrarefy(comm.data[,-c(1:26)],1236))#by otu
#comm.data[,27:243]<-comm.data[,27:243]/rowSums(comm.data[,27:243])*100

greater66plants<-factor(ifelse(comm.data$Plant_Dens>66,"hi","lo")) #this is stem density, including mosses

#nohilo<-ifelse(comm.data$Plant_Dens==0,"no","else");nohilo[which(comm.data$Plant_Dens<43&comm.data$Plant_Dens>0)]<-"lo";nohilo[which(comm.data$Plant_Dens>=43)]<-"hi";nohilo<-factor(nohilo)#14 no,41 lo, 42 hi ##THESE SAMPLE SIZES ARE INCORRECT
lomehi<-ifelse(comm.data$Plant_Dens<36,"lo","else");lomehi[which(comm.data$Plant_Dens<91&comm.data$Plant_Dens>=36)]<-"me";lomehi[which(comm.data$Plant_Dens>=91)]<-"hi";lomehi<-factor(lomehi)#32 no, 33 lo, 32 hi 

comm.data<-cbind(lomehi,comm.data)
trts<-as.vector(unique(lomehi))

###XXXXXXXXXXXcould also make the loop quicker if I take out doubletons and singletons here, since the rarefaction probably introduced some. there aren't very many, I'll just remove them within the loop
#which(colSums(comm.data[,28:2209])<1)


#Read in plant data if you want
plantcomp<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots; there are some plots that have plant data but not microbe data. 69 70 71 77 81 108 117 118 147 148 149 151. This is because when we started doing the surveys we were going to all plots for plants and only sample some for microbes, then we realized that that was insane!
dim(plantcomp)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantlabels<-as.data.frame(cbind(colnames(plantcomp2)[2:56],colnames(plantcomp2)[2:56],"Plant"))
colnames(plantlabels)<-c("otu","orders","kingdomlabels")

#merge plants with microbes
comm.data$Sample_name<-as.numeric(as.character(comm.data$Sample_name))
comm.datamp<-merge(comm.data,plantcomp2,"Sample_name",sort=F,all.y=F)
comm.datamp$Sample_name


#Setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
#loop
results<-matrix(nrow=0,ncol=11)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp1<-subset(comm.datamp, lomehi==trt.temp)####change to comm.data if not doing plants

  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  temp<-cbind(temp1[,1:27],temp1[,((which(colSums(temp1[,28:dim(temp1)[2]]>0)>2))+27)])#for me should be 907 from denovo946 to trispi
  
  #in this case the community data started at column 28, so the loop for co-occurrence has to start at that point
  for(b in 28:(dim(temp)[2]-1)){
    
    results1<-foreach(c=(b+1):(dim(temp)[2]),.combine=rbind,.packages="Kendall") %dopar% {
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(temp[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)
      #I changed this so that it will calculate the correlation for all species pairs, even singletons, and I can subset them later if I want to remove infrequent otus for qvalue calculation
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
        test<-Kendall(temp[,b],temp[,c])
        kendallrho<-test$tau
        kendallp.value<-test$sl
        
        #do randomization test for p value
        #ran1<-t(replicate(1000,sample(temp[,b],length(temp[,b]))))
        #ran2<-t(replicate(1000,sample(temp[,c],length(temp[,c]))))
        #ken<-apply(cbind(ran1,ran2),1,function(x){Kendall(x[1:length(temp[,b])],x[(length(temp[,b])+1):(length(temp[,b])*2)])$tau})
        #p.value.perm<-length(which(ken>rho))/1000#the p value is not correct for negative correlations
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
        kendallrho<-0
        kendallp.value<-1
      }	
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],spearmanrho,spearmanp.value,kendallrho,kendallp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
      new.row		
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt)
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","kendallrho","kendallp.value","ab1","ab2","ab1freq","ab2freq")

#from run with spearman
resultsS<-results
resultsno4orlower<-subset(results,ab1freq>4&ab2freq>4&rho>0)
resultsno5orlower<-subset(results,ab1freq>5&ab2freq>5&rho>0)

resultsK<-results
resultsno4orlowerK<-subset(resultsK,ab1freq>4&ab2freq>4&rho>0)
resultsno4abunorlowerK<-subset(resultsK,ab1>4&ab2>4&rho>0)

#results from partial run with randomization p value
resultsP<-results

#results from kendall and spearman
resultsKS<-results
resultsKSno4<-subset(resultsKS,ab1freq>4&ab2freq>4)
resultsKSno5<-subset(resultsKS,ab1freq>5&ab2freq>5)#like the fierer paper, more than 5 sequences (although he did not rarefy??? weird)
resultsKSno3<-subset(resultsKS,ab1freq>3&ab2freq>3)

#results from kendall and spearman, with a 3 sample frequency cutoff and removal of otus with less than .2% relative abundance. took 4 hrs to run
resultsKS32<-results
resultsKS32no2<-subset(resultsKS32,ab1freq>2&ab2freq>2)
head(resultsKS32no2)

#results from kendall and spearman, with a 3 sample frequency cutoff and removal of otus with less than .2% relative abundance with plant data. took 40 min to run (i did the removal of zeros 1s and 2s in the temp file within the loop)
resultsKS32p<-results
#resultsKS32pno2<-subset(resultsKS32p,ab1freq>2&ab2freq>2) #not necessary
head(resultsKS32p)












#Looking into the pvalue corrections, ignore this
b=34
c=37
sum(temp[,34])
sum(temp[,37])

cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
plot(temp[,b],temp[,c])
plot(rank(temp[,b]),rank(temp[,c]))

temp$denovo1219 
temp$denovo218063
cor.test(temp$denovo1219 ,temp$denovo218063,method="spearman",na.action=na.rm)
plot(temp[,b],temp[,c])
plot(jitter(rank(temp$denovo1219)),rank(temp$denovo218063))

p.adjust(c(.01,.02,.02,.03,.03,.7,.7,.8,.8,.9))
p.adjust(c(.01,.02,.03,.7,.8,.9))

head(edge_listsno4orlower)
temphi<-subset(tempno,trt=="hi")
sort(temphi$p.value)[400:500]
sort(p.adjust(temphi$p.value, method="fdr"))[400:500]
which(p.adjust(temphi$p.value, method="fdr")<.05)

templo<-subset(edge_listsno4orlower,trt=="lo")


tempno<-subset(edge_lists2,trt=="no")
sort(tempno$p.value)[1:100]
sort(p.adjust(tempno$p.value, method="fdr"))[1:100]
sort(p.adjust(c(tempno$p.value,tempno$p.value), method="fdr"))[1:100]

length(which(temphi$p.value<.0001))
length(which(temphi$p.value>=.0001))

length(which(templo$p.value<.0001))
length(which(templo$p.value>=.0001))

length(which(tempno$p.value<.0001))
length(which(tempno$p.value>=.0001))

#tempno and templo have proportionally fewer really low p-values compared to the temphi. it also appears that because of this (?) the pvalue cutoff for being a qval<0.05 is higher for the hi treatment
#if the p.adjust method of fdr correction, the more significant pvalues you in a certain range (say .001 to .005) the less the p value is adjusted. which is probably why the hi plant density has more significant p values after adjustment - there were more significant p values to begin with. To test this, I calculated adjusted values on all pvalues together, but it still didn't really solve the problem and I'm left with very few in no plant plots
#I'm convincing myself that the no plant plots have more zeros (the mean frequency is lower) thus the correlation coefficients are stronger because there are many points anchored at 0,0 and because of the ties the p value is adjusted so that it is not as significant


temphi[ind,][10:20,]
head(tempno)
plot(tempno$p.value,(tempno$ab1+tempno$ab2)/2)
abline(lm((tempno$ab1+tempno$ab2)/2~tempno$p.value))

tempno<-subset(comm.data,nohilo=="no")
templo<-subset(comm.data,nohilo=="lo")
temphi<-subset(comm.data,nohilo=="hi")

edge_listsno4orlowerhi<-subset(edge_listsno4orlowerK,trt=="hi")
edge_listsno4orlowerlo<-subset(edge_listsno4orlowerK,trt=="lo")
edge_listsno4orlowerno<-subset(edge_listsno4orlowerK,trt=="no")

ind<-which(edge_listsno4orlowerhi$qval1<.05)
edge_listsno4orlowerhi[ind,][1:100,]

ind<-which(edge_listsno4orlowerlo$qval1<.05)
edge_listsno4orlowerlo[ind,][1:100,]

ind<-which(edge_listsno4orlowerno$p.value<.05&edge_listsno4orlowerno$rho>.5)
edge_listsno4orlowerno[ind,][1:10,]

avgh<-(edge_listsno4orlowerhi$ab1freq+edge_listsno4orlowerhi$ab2freq)/2
plot(avgh,edge_listsno4orlowerhi$p.value)
abline(lm(edge_listsno4orlowerhi$p.value~avgh),col=2)
summary(lm(edge_listsno4orlowerhi$p.value~avgh))


otu1<-temphi$A31
otu2<-temphi$Archaeosporales
cbind(otu1,otu2)
cor.test(otu1,otu2,method="spearman",na.action=na.rm)
Kendall(otu1,otu2)
plot(jitter(rank(otu1)),rank(otu2))
plot(otu1,otu2)

#to get a permutation p value
ran1<-data.frame(matrix(nrow=0,ncol=14))
for(x in 1:1000){ran1[x,]<-sample(otu1,14)}
ran2<-data.frame(matrix(nrow=0,ncol=14))
for(x in 1:1000){ran2[x,]<-sample(otu2,14)}
ken<-data.frame(matrix(nrow=0,ncol=1))
for(x in 1:20000){ken[x,1]<-Kendall(ran1[x,],ran2[x,])$tau}
length(which(ken[,1]>Kendall(otu1,otu2)$tau))/20000

#faster way to permute
ran1<-t(replicate(20000,sample(otu1,14)))
ran2<-t(replicate(20000,sample(otu2,14)))
ken2<-apply(cbind(ran1,ran2),1,function(x){Kendall(x[1:14],x[15:28])$tau})
length(which(ken2>Kendall(otu1,otu2)$tau))/20000

strt<-Sys.time()
print(Sys.time()-strt)

cor.test(resultsno4orlowerlo$denovo3495 ,resultsno4orlowerlo$denovo24997,method="spearman",na.action=na.rm)
plot(jitter(resultsno4orlowerlo$denovo3495),resultsno4orlowerlo$denovo24997)
plot(jitter(rank(resultsno4orlowerlo$denovo3495)),rank(resultsno4orlowerlo$denovo24997))

cor.test(resultsno4orlowerhi$denovo16645 ,resultsno4orlowerhi$denovo318490,method="spearman",na.action=na.rm)
plot(jitter(resultsno4orlowerhi$denovo16645),resultsno4orlowerhi$denovo318490)
plot(jitter(rank(resultsno4orlowerhi$denovo16645)),rank(resultsno4orlowerhi$denovo318490))

ind<-which(templo$Erysiphales+templo$Boletales!=0)
plot(jitter(rank(templo$Boletales)),jitter(rank(templo$Erysiphales)))
abline(lm(rank(templo$Erysiphales)~rank(templo$Boletales)),col=2)
cor.test(templo$Erysiphales,templo$Boletales,method="pearson")
Kendall(templo$Erysiphales,templo$Boletales)
#taking zeros out will reduce rho and increase p
cor.test(ifelse(templo$Erysiphales>0,1,0),ifelse(templo$Boletales>0,1,0),method="pearson")
chisq.test(ifelse(templo$Erysiphales>0,1,0),ifelse(templo$Boletales>0,1,0),correct=T,simulate.p.value = T)
pcc(as.matrix(rbind(as.integer(ifelse(templo$Erysiphales>0,1,2)),as.integer(ifelse(templo$Boletales>0,1,2)))),corrected=F)
pcc(as.matrix(rbind(as.integer(ifelse(templo$Monoblepharidales>0,1,2)),as.integer(ifelse(templo$Boletales>0,1,2)))),corrected=F)
pcc(as.matrix(rbind(as.integer(c(2,2,ifelse(templo$Monoblepharidales>0,1,2))),as.integer(c(2,2,ifelse(templo$Boletales>0,1,2))))),corrected=F)
phi(matrix(c(7,1,23,11),nrow=2,byrow=T),digits=6)


#pcc is different from phi but overall pretty similar, the maximum value for phi is one less than the smaller of the two dimensions of the contingency table. due to this it is finefor comparison among different contingency tables of the same size (but not of different sizes). both adjust the statistic by sample size phi=sqrt(chisq/n). pcc may be a totally fine alternative in the future if the code is faster
cbind(ifelse(templo$Erysiphales>0,1,0),ifelse(templo$Boletales>0,1,0))
Phi(ifelse(templo$Erysiphales>0,1,0),ifelse(templo$Boletales>0,1,0))
Phi(rep(ifelse(templo$Erysiphales>0,1,0),2),rep(ifelse(templo$Boletales>0,1,0),2))
ContCoef(ifelse(templo$Erysiphales>0,1,0),ifelse(templo$Boletales>0,1,0))






