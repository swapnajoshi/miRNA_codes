library(NanoStringNorm)
library(ggplot2)
library(limma)
library(Biobase)
library(DESeq)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")

setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/tissueRCC/")
rccfiles<-read.markup.RCC(rcc.pattern="*.RCC|*.rcc",  exclude=NULL, include = NULL, nprobes=-1)
#Returns a list with two components. The first is the header information which contains sample IDs
#and diagnostic information on the quality of the samples. The second is the count data and can be
#directly used in the input to NanoStringNorm.

tissueMirDf<-rccfiles$x

phenoDat<-colnames(tissueMirDf)[4:48]

phenoDat1 <- as.data.frame(matrix(unlist(strsplit(phenoDat,"_")), nrow=45, ncol=6, byrow=TRUE))
row.names(phenoDat1) <-phenoDat
colnames(phenoDat1)<-c("NDP_number","Condition","Sex","Lane","Tissue_type","Run_date")
phenoDat1$Cartridge<-substr(phenoDat1$Lane,1,1)########### however, A and B are repeated as cartridge names, hence they cannot be same as batch
phenoDat1<-phenoDat1[order(phenoDat1$Cartridge, phenoDat1$Lane),]

TS = factor(phenoDat1$Cartridge,levels=c("A","B","C","F"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design1<-as.data.frame(design)
design1$A<-gsub("1","2",design1$A)
design1$B<-gsub("1","2",design1$B)
design1$D<-gsub("1","2",design1$C)
design1$E<-gsub("1","2",design1$F)

design1$A<-gsub("0","1",design1$A)
design1$B<-gsub("0","1",design1$B)
design1$D<-gsub("0","1",design1$C)
design1$E<-gsub("0","1",design1$F)

row.names(design1)<-phenoDat
design1[,1]<-as.numeric(as.character(design1[,1]))
design1[,2]<-as.numeric(as.character(design1[,2]))
design1[,3]<-as.numeric(as.character(design1[,3]))
design1[,4]<-as.numeric(as.character(design1[,4]))


#trait2
phenoDat1$status<-phenoDat1$Condition
phenoDat1$status<-gsub("IBSC","IBS",phenoDat1$status)
phenoDat1$status<-gsub("IBSD","IBS",phenoDat1$status)
TS1 = factor(phenoDat1$status,levels=c("IBS","HC"))
cond <- model.matrix(~0+TS1)
colnames(cond) <- levels(TS1)
cond1<-as.data.frame(cond)
cond1$IBS<-gsub("1","2",cond1$IBS)
cond1$HC<-gsub("1","2",cond1$HC)
cond1$IBS<-gsub("0","1",cond1$IBS)
cond1$HC<-gsub("0","1",cond1$HC)

row.names(cond1)<-phenoDat
cond1[,1]<-as.numeric(as.character(cond1[,1]))
cond1[,2]<-as.numeric(as.character(cond1[,2]))

#positive control and negative control sanity check


posC<-tissueMirDf[1:6,]
row.names(posC)<-posC[,2]


posC<-posC[,4:48]
posC1<-t(posC[sort(row.names(posC)),])

png("boxplot.posControls.png")
boxplot(posC1) 
dev.off()
boxplot(posC1) 

# for (i in 1:45)
# {
#   
#   png(paste(i,"posCon.png"))	
#   barplot(posC1[i,]) # to look at individual samples
#   dev.off()
# }


negC<-tissueMirDf[7:12,]
row.names(negC)<-negC[,2]
negC<-t(negC[,4:48])

png("boxplot.negControls.png")
boxplot(negC) 
dev.off()

tissueMirDf1<-tissueMirDf
tissueMirDf1$"A6145_HC_F_F09_Tis_Apr15" <- NULL
tissueMirDf1$"A5813_IBSC_M_B07_Tis_Feb15" <- NULL

nano.norm<-NanoStringNorm(tissueMirDf1,anno = NA,header = NA,
                          Probe.Correction.Factor ='none',CodeCount ='geo.mean',
                          Background = 'none', SampleContent =	'top.geo.mean',	OtherNorm =	'none',
                          round.values = TRUE,
                          is.log = FALSE,
                          take.log = FALSE,
                          return.matrix.of.endogenous.probes = FALSE,
                          traits = NA,
                          predict.conc = TRUE,
                          verbose = TRUE,
                          genes.to.fit,
                          genes.to.predict,
                          guess.cartridge = TRUE
)


norm.data<-nano.norm$normalized.data
norm.data1<-subset(norm.data, norm.data$Code.Class!="Negative" & norm.data$Code.Class!="Housekeeping" & norm.data$Code.Class!="Positive" & norm.data$Code.Class!="SpikeIn")
norm.data2<-norm.data1[,4:46]

png('NanoStringNorm_Example_Plots_tissue_ibs_hc_%03d_topgeo-positivecontoldrop_2sdbckg.png', units ='in',  height = 6,   width = 6, res = 250, pointsize = 10);
Plot.NanoStringNorm( x = nano.norm,
                     label.best.guess = TRUE,
                     plot.type = c('cv','mean.sd','RNA.estimates','volcano','missing','norm.factors','positive.controls','batch.effects'))
dev.off()

Plot.NanoStringNorm( x = nano.norm,
                     label.best.guess = TRUE,
                     plot.type = c('cv','mean.sd','RNA.estimates','volcano','missing','norm.factors','positive.controls','batch.effects'))


raw.data<-nano.norm$raw.data

phenoDat1 <- phenoDat1[!row.names(phenoDat1)%in%c("A6145_HC_F_F09_Tis_Apr15","A5813_IBSC_M_B07_Tis_Feb15"),]
rawdata1<-raw.data[,4:46]
raw.data1<-rawdata1[,row.names(phenoDat1)]
mean1.sam<-apply(raw.data1,2,mean)

condition <- relevel(as.factor(phenoDat1$status), ref="HC")


# DESeq -------------------------------------------------------------------

## Make a new countDataSet
d <- newCountDataSet(norm.data2, group = condition)

## Estimate library size and dispersion
d <- estimateSizeFactors(d)
d <- estimateDispersions(d)
plotDispEsts(d, main="DESeq: Per-gene dispersion estimates")

## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
print(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition", "libType")))

## Fit full and reduced models, get p-values
dfit1 <- fitNbinomGLMs(d, count~libType+condition)
dfit0 <- fitNbinomGLMs(d, count~libType)
dpval <- nbinomGLMTest(dfit1, dfit0)
dpadj <- p.adjust(dpval, method="BH")

## Make results table with pvalues and adjusted p-values
dtable <- transform(dfit1, pval=dpval, padj=dpadj)
dtable <- dtable[order(dtable$padj), ]
head(dtable)
















# edgeR -------------------------------------------------------------------
row.names(tissueMirDf) <- tissueMirDf[,2]
tissueMirDfx <- tissueMirDf[,-c(1:3)]
## Make design matrix
condition <- relevel(as.factor(phenoDat1$status), ref="HC")
# libType <- factor(meta$libType)
edesign <- model.matrix(~condition)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
e <- DGEList(counts=tissueMirDfx, group = condition)
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

## Fit the model, testing the coefficient for the treated vs untreated comparison
efit <- glmFit(e, edesign)
efit <- glmLRT(efit, coef="conditionIBS")

## Make a table of results
etable <- topTags(efit, n=nrow(e))$table
etable <- etable[order(etable$FDR), ]
head(etable)

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

