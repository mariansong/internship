
#ANOVA


#######import the table 
library(readxl)
mainAll_time <- read_excel("mainAll_time.xls")
View(mainAll_time)
a<-data.frame(mainAll_time)
b<-t(a)
head(b)
write.csv(b,file = "b.csv")

#Read in the modified file "expr_max_change"
###
library(readr)
b_sample1 <- read_csv("b_time.csv")
sample1<-data.frame(b_sample1)

#Then, the time variable is set to the categorical variable.
sample1$time<-factor(sample1$time,levels = c("0","2","8","16"))
sample1$time
#Then. We can do a parametric or non-parametric ANOVA of.

#Parametric ANOVA (slow; simple for loop)

baseformula <- " ~ time"
#n Counting, initialization
n<-0
#a<-array(data = NA,dim=length(sample))
#b. Store the result of the for loop
b<-array(data = NA,dim=length(sample))
for (i in 2:ncol(sample1)) {
  formula <- paste(colnames(sample1)[i], baseformula, sep="")
  
  p <- summary(aov(as.formula(formula), data=sample1))[[1]][["Pr(>F)"]][1]
  
  #set  P<0.1
  if(p<0.1){
    n=n+1
    #a[n]<-formula
    print(paste(formula, ": p=", p, sep=""))
    b[n]<-paste(formula, ": p=", p, sep="")
  }
  
  #write.table(a,file="aa.txt")
  write.table(b,file="bb.txt")

}

print(n)
#[1] 7536
###Under the condition of P<1, a total of 7536 differential genes were obtained, with duplicates



###limma
#Bioconductor limma
BiocManager::install('limma')

#load limma
library(limma)

#Data pre-processing to remove duplicate gene names
library(readxl)
X1DE <- read_excel("mainAll_time_1.xls")
X1<-data.frame(X1DE)
head(X1)

#For the same gene, we should pick the whole row with the larger row mean and not break it up.
#Calculate row averages, sorted in descending order
expr<-X1
index=order(rowMeans(expr[,-1]),decreasing = T)
#Adjusting the gene order of expression profiles
expr_ordered=expr[index,]
#For genes with duplicates, keep the one with the first occurrence, i.e. the one with the larger row mean
keep=!duplicated(expr_ordered$Gene)
#Obtain the expression spectrum matrix after the final processing

expr_max=expr_ordered[keep,]
expr_max

#write results
write.table(expr_max, file = "expr_max.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#Gene expression matrix

#exprSet <- read.delim('expr_max _change.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#exprSet1 <- read_excel("expr_max _change.xls")

#write.table(exprSet1,file = "exprSet1.txt",sep="\t",row.names = T,col.names = F,quote = F)
#exprSet2<-read.table("exprSet1.txt")
#exprSet<-data.frame(exprSet1)

exprSet <- read.delim('expr_max_change.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

boxplot(exprSet)
plotDensities(exprSet)

#Experimental design matrix
#Be careful to ensure that the order of samples in the expression matrix and the grouping order here are one-to-one correspondence
group <- factor(rep(c('control', 'treat'), each = 4), levels = c('control', 'treat'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(group)
design

#Filter low count data, e.g. directly on selected count values
#exprSet <- exprSet[rowSums(exprSet) >= 50, ]

#voom  normalize
norm <- voom(exprSet, design, plot = TRUE)

boxplot(norm$E)
plotDensities(norm$E)

#Output the normalized factor expression matrix
#write.table(norm$E, 'gene_normalized.txt', col.names = NA, sep = '\t', quote = FALSE)

#Linear fit
fit <- lmFit(norm, design, method = 'ls')

#Identify the two groups for comparison
#The up-/down-regulation status of the gene expression values of the group with marker 1 relative to the group with marker -1 will be calculated later.
contrast <- makeContrasts('treat-control', levels = design)
contrast

#Fitting standard errors using an empirical Bayesian model

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2) 

qqt(fit2$t, df = fit2$df.prior+fit2$df.residual, pch = 16, cex = 0.2)
abline(0,1)

#p-value correction, extraction of analysis of variance results
diff_gene <- topTable(fit2, number = Inf, adjust.method = 'fdr')
head(diff_gene, 10)

write.table(diff_gene, 'gene_diff.txt', col.names = NA, sep = '\t', quote = FALSE)
#Afterwards, the exact gene name needs to be matched back again

##picture
#Example of a volcano map
#For example, here "difference" is defined according to |log2FC| >= 1 & FDR p-value < 0.01
diff_gene[which(diff_gene$logFC >= 1 & diff_gene$adj.P.Val < 0.01),'sig'] <- 'red'
diff_gene[which(diff_gene$logFC <= -1 & diff_gene$adj.P.Val < 0.01),'sig'] <- 'blue'
diff_gene[which(abs(diff_gene$logFC) < 1 | diff_gene$adj.P.Val >= 0.01),'sig'] <- 'gray'

log2FoldChange <- diff_gene$logFC
FDR <- diff_gene$adj.P.Val
plot(log2FoldChange, -log10(FDR), pch = 20, col = diff_gene$sig)
abline(v = 1, lty = 2)
abline(v = -1, lty = 2)
abline(h = -log(0.01, 10), lty = 2)




