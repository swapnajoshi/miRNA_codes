
library(NanoStringNorm)
library(ggplot2)
library(limma)
library(Biobase)
library(edgeR)
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")

setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/tissueRCC/")

# count data


# edgeR -------------------------------------------------------------------


## Make design matrix
condition <- relevel(as.factor(phenoDat1$status), ref="HC")
# libType <- factor(meta$libType)
edesign <- model.matrix(~condition)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
e <- DGEList(counts=norm.data2, group = condition)
e$samples
keep <- rowSums(cpm(e)>1) >= 2
e <- e[keep, , keep.lib.sizes=FALSE]; dim(e)

# e <- calcNormFactors(e)
e <- estimateGLMCommonDisp(e, edesign)
# e <- estimateGLMTrendedDisp(e, edesign) 
# e <- estimateGLMTagwiseDisp(e, edesign)

## MDS Plot
png("MDS_edgeR.png")
plotMDS(e, main="edgeR MDS Plot")
dev.off()
plotMDS(e, main="edgeR MDS Plot")

## Biological coefficient of variation plot
png("Biological coefficient of variation plot.png")
plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")
dev.off()
plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")

##############
et <- exactTest(e)
topTags(et)
etable <- topTags(et, n=nrow(e))$table
etable <- etable[order(etable$FDR), ]

# ## Fit the model, testing the coefficient for the treated vs untreated comparison
# efit <- glmFit(e, edesign)
# efit <- glmLRT(efit, coef="conditionIBS")
# 
# ## Make a table of results
# etable <- topTags(efit, n=nrow(e))$table
# etable <- etable[order(etable$FDR), ]
# head(etable)

## ~MA Plot

with(etable, plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
with(subset(etable, PValue<0.05), points(logCPM, logFC, pch=20, col="red"))
abline(h=c(-1,1), col="blue")

png("FOLD CHANGE VS ABUNDANCE_edgeR.png")
with(etable, plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
with(subset(etable, PValue<0.05), points(logCPM, logFC, pch=20, col="red"))
abline(h=c(-1,1), col="blue")
dev.off()

write.table(etable, file = "IBSvsHC_miRNA_DE.csv", sep = ",", col.names = NA)

