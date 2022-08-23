###limma
BiocManager::install("plotDensities")
install.packages(plotDensities)
library(plotDensities)
#write results
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
BiocManager::install("voom")
install.packages("voom")
library(voom)
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
#For example, here "difference" is defined according to |log2FC| >= 1 & FDR p-value < 0.02
diff_gene[which(diff_gene$logFC >= 1 & diff_gene$adj.P.Val < 0.02),'sig'] <- 'red'
diff_gene[which(diff_gene$logFC <= -1 & diff_gene$adj.P.Val < 0.02),'sig'] <- 'blue'
diff_gene[which(abs(diff_gene$logFC) < 1 | diff_gene$adj.P.Val >= 0.02),'sig'] <- 'gray'

log2FoldChange <- diff_gene$logFC
FDR <- diff_gene$adj.P.Val
plot(log2FoldChange, -log10(FDR), pch = 20, col = diff_gene$sig)
abline(v = 1, lty = 2)
abline(v = -1, lty = 2)
abline(h = -log(0.01, 10), lty = 2)
