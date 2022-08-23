#Analysis explanation
#Two-way ANOVA for each gene, having as factors the group and the time.
#The data file has expression values for each gene in different time points, for two groups.
#Afer applying ANOVA, results are then corrected for multiple comparisons (multiple genes).

##Loading packages
rm(list=ls())
library(dplyr)
library(tidyr)

##Reading the raw data
totalData<-read.csv("gene1.csv")
##General process
#For each gene, data will be read from the original table, stored in a temporary input
#object, run through the ANOVA and the results will be stored in a final output ANOVA table.
#The process will be repeated for all genes.
##PREPARING THE TEMPORARY INPUT TABLE

#######33
groups<-c("0h","2h","8h","16h")
nGroups<-4
nSamples<-2  #number samples collected for each group
nReplicates<-2 #number of replicate samples for each time point in each group
#Object to store the data for a single gene in every iteration of the ANOVA
dataTr<-data.frame(group=rep(groups,each=nSamples),
                   time=rep(rep(c("a","b"),each=nReplicates),nGroups),
                   expression=rep(NA,nSamples*nGroups))




##PREPARING THE OUTPUT TABLE
#Object to store two-way ANOVA result
nGenes<-length(totalData$GeneNames)
anova2way_data<-data.frame(
  GeneNames=totalData$GeneNames,
  ANOVAF=rep(NA,nGenes),
  ANOVAP=rep(NA,nGenes)
)

##Performing two-way ANOVA for each gene, using a for loop
for(i in seq(nGenes)){   #For each gene
  dataCols<-(2:9)   #columns in the original table that have expression data
  #Gets the data in the original table and stores in the input table dataTr.
  dataTr$expression <- unlist(totalData[i,dataCols])
  
  #Performs two-way ANOVA for the given gene
  d_obj <- lm(expression ~ group*time , data=dataTr)
  d_anova <- anova(d_obj)
  #Stores the ANOVA result in the output object
 
  anova2way_data$ANOVAF[i] <- d_anova$`F value`[1]
  anova2way_data$ANOVAP[i] <-  d_anova$`Pr(>F)`[1]
}



###Qvalue

library(qvalue)
pvalues <- anova2way_data$ANOVAP

#Checking the p-value histogram
#hist(hedenfalk$p, nclass = 20)
hist(anova2way_data$ANOVAP,nclass = 20)

qobj <- qvalue(p = pvalues)
names(qobj)
summary(qobj)
hist(qobj)
plot(qobj)

qvalues <- qobj$qvalues

