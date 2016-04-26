#Loading/saving/packages needed

setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill")

#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU4.Rdata") #results has the results from OTU cooccurrence patterns, 2 is with plants/bac/euk, 3 is with nematodes, 4 is with .tre files for bacteria

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_Analysis_byOTU4.Rdata")

#load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill_Workspace_AnalysisOrder.Rdata") #has old order level analyses for euk data


library(vegan)
library(reshape)
library(plotrix)
library(foreach)
library(doParallel)
library(Kendall)
library(tidyr)
library(grid)
library(phyloseq)
library(foreach)
library(doParallel)
library(data.table)