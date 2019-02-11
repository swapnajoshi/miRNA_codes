# TODO: NanoStringNorm
#import, normalization serum, normalization tissue, separate analysis, correlation between tissue and serum, target selection 
# 
# Author: SwapnaJoshi 4102016

###############################################################################

#Funcitons and libraries:

Wilcox.test <- function(x, s1, s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  wx.out <- wilcox.test(x1, x2, alternative = "two.sided",
                        paired = FALSE)
  p.Wilcox <- as.numeric(wx.out$p.value)
  return(p.Wilcox)
}

qqplot.data <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  d <- data.frame(resids = vec)
  
  ggplot(d, aes(sample = resids)) + stat_qq() + geom_abline(slope = slope, intercept = int)
  
}

# source("http://bioconductor.org/biocLite.R")
# biocLite("vsn")
library(NanoStringNorm)
library(ggplot2)
library(limma)


####################################################################################

setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/Serum_RCC/")
rccfiles<-read.markup.RCC(rcc.pattern="*.RCC|*.rcc",  exclude=NULL, include = NULL, nprobes=-1)
#Returns a list with two components. The first is the header information which contains sample IDs
#and diagnostic information on the quality of the samples. The second is the count data and can be
#directly used in the input to NanoStringNorm.

serumMirDf<-rccfiles$x

phenoDat<-colnames(serumMirDf)[4:48]

phenoDat1<-as.data.frame(matrix(unlist(strsplit(phenoDat,"_")), nrow=45, ncol=6, byrow=TRUE))
row.names(phenoDat1)<-phenoDat
colnames(phenoDat1)<-c("NDP_number","Condition","Sex","Lane","Tissue_type","Run_date")
phenoDat1$Cartridge<-substr(phenoDat1$Lane,1,1)########### however, A and B are repeated as cartridge names, hence they cannot be same as batch
phenoDat1<-phenoDat1[order(phenoDat1$Cartridge, phenoDat1$Lane),]

# TS = factor(phenoDat1$Cartridge,levels=c("A","B","D","E"))
# design <- model.matrix(~0+TS)
# colnames(design) <- levels(TS)
# design1<-as.data.frame(design)
# design1$A<-gsub("1","2",design1$A)
# design1$B<-gsub("1","2",design1$B)
# design1$D<-gsub("1","2",design1$D)
# design1$E<-gsub("1","2",design1$E)
# 
# design1$A<-gsub("0","1",design1$A)
# design1$B<-gsub("0","1",design1$B)
# design1$D<-gsub("0","1",design1$D)
# design1$E<-gsub("0","1",design1$E)
# 
# row.names(design1)<-phenoDat
# design1[,1]<-as.numeric(as.character(design1[,1]))
# design1[,2]<-as.numeric(as.character(design1[,2]))
# design1[,3]<-as.numeric(as.character(design1[,3]))
# design1[,4]<-as.numeric(as.character(design1[,4]))

################################################################# trait2
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/Results/Plasma_miR/")

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


posC<-serumMirDf[1:6,]
row.names(posC)<-posC[,2]


posC<-posC[,4:48]
posC1<-t(posC[sort(row.names(posC)),])

png("boxplot.posControls.png")
boxplot(posC1) 
dev.off()

for (i in 1:45)
{
  
  png(paste(i,"posCon.png"))	
  barplot(posC1[i,]) # to look at individual samples
  dev.off()
}

# 34 [1] "A5765_IBSC_F_D02_Ser_Apr15" doesnt look good others are good.  Finally included them in the analyses


negC<-serumMirDf[7:12,]
row.names(negC)<-negC[,2]
negC<-t(negC[,4:47])

png("boxplot.negControls.png")
boxplot(negC) 
dev.off()
boxplot(negC) 

serumMirDf1<-serumMirDf[-7,]


nano.norm<-NanoStringNorm(serumMirDf1,anno = NA,header = NA,
                          Probe.Correction.Factor ='none',CodeCount ='geo.mean',
                          Background = 'mean',SampleContent =	'top.geo.mean',	OtherNorm =	'none',
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


# A5765_IBSC_F_D02_Ser_Apr15 needs to be removed high positive normalization factor,12  (5 others have high, but with in 5) )

norm.data<-nano.norm$normalized.data
norm.data1<-subset(norm.data, norm.data$Code.Class!="Negative" & norm.data$Code.Class!="Housekeeping" & norm.data$Code.Class!="Positive" & norm.data$Code.Class!="SpikeIn")

norm.data2<-norm.data1[,4:48]
norm.data2$A5765_IBSC_F_D02_Ser_Apr15 <- NULL

png('NanoStringNorm_Example_Plots_serum_ibs_hc_%03d_topgeo-positivecontoldrop_2sdbckg.png', units ='in',  height = 6,   width = 6, res = 250, pointsize = 10);
Plot.NanoStringNorm( x = nano.norm,
                     label.best.guess = TRUE,
                     plot.type = c('cv','mean.sd','RNA.estimates','volcano','missing','norm.factors','positive.controls','batch.effects'))
dev.off()







############################################################# q q plots
norm.data4.t<-t(norm.data4)
norm.data4.t<-as.data.frame(norm.data4.t)
norm.data4.t$Condition<-phenoDat1$status

#parametric or non parametric? generate a q-q plot
meanDat<-phenoDat1
meanDat$mean1<-mean2.sam
library(ggplot2)

qqplot1<-qqplot.data(meanDat$mean1)
ggsave(qqplot1, file="qqplotSerum.mean.IBS1.png") ## 
qqplot1

ibs<-rownames(subset(phenoDat1, phenoDat1$status=="IBS"))
hc<-rownames(subset(phenoDat1, phenoDat1$status=="HC"))
ibsc<-rownames(subset(phenoDat1, phenoDat1$Condition=="IBSC"))
ibsd<-rownames(subset(phenoDat1, phenoDat1$Condition=="IBSD"))



##########################################################################
exprs1<-norm.data4
colnames(exprs1)<-substr(colnames(exprs1),1,5)

## limma

#import phenotype data

row.names(phenoDat1)<-phenoDat1[,1]

match(row.names(phenoDat1),colnames(exprs1))
pData <- as.matrix(phenoDat1[,2]); row.names(pData)=row.names(phenoDat1); colnames(pData)="id"
pData<-(as.data.frame(pData))

phenoData <- new("AnnotatedDataFrame", data=pData)

exprDat<-new("ExpressionSet", exprs=as.matrix(exprs1), phenoData=phenoData)

#expressionset
TS = factor(pData(phenoData)$id,levels=c("IBSC","IBSD","HC"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
v <- voom(exprs(exprDat),design,plot=TRUE)
fit <- lmFit(v, design)

#contrast matrix
cont.matrix <- makeContrasts(IBS_HC = (IBSC+IBSD)/2-HC, IBSC-HC, IBSD-HC, IBSD-IBSC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
ibs_fit2 <- eBayes(fit2)

png("IbsMirPlasmaVolca1.png")
volcanoplot(ibs_fit2)
dev.off()

volcanoplot(ibs_fit2)
#Differentially expressed genes	
topTable(ibs_fit2, number=10, coef=1)
# write.fit(ibs_fit2, results = NULL, file = "AllResults_plasma.csv", method = "separate", sep = ",")
resultsIBS=topTable(ibs_fit2, number=327)
indiv.P.values <- ibs_fit2$p.value
indiv.q.values <- apply(ibs_fit2$p.value, 2, p.adjust, method="fdr")

pvalues<-cbind(indiv.P.values, indiv.q.values)
pMirRes <- merge(resultsIBS, pvalues, by = "row.names")

row.names(pMirRes) <- pMirRes[,1]
pMirRes$Row.names <- NULL
foldChange <- function(x) {ifelse (x >= 0, round(2^x,2), round(2^(-1*x),2))}

pMirRes[,c("FC_IBS_HC","FC_IBS_C_HC","FC_IBS_D_HC","FC_IBS_D_IBS_C")] <- foldChange(pMirRes[,c(1,2,3,4)])

tableSig <- function(x)
{ a1 <- subset(x, x[,1] < 0.05)
a2 <- subset(x, x[,2] < 0.05)
a3 <- subset(x, x[,3] < 0.05)
a4 <- subset(x, x[,4] < 0.05) 
y <- rbind( a1, a2, a3, a4)
z <- colnames(x[,1:4])
Comparison <- cbind(c(rep(z[1], dim(a1)[1])),rep(z[2], dim(a2)[1]), rep(z[3], dim(a3)[1]), rep(z[4], dim(a4)[1]))
y <- y$Comparison
return(y)  
}

SigRes <- tableSig(pMirRes[,c(9:12)])

save(resultsIBS, file="resultsIBSPlasma.rda")


resultsIBS$ibs_cont_sig0.05<-resultsIBS$P.Value<0.05
resultsIBS$fdr_0.1<-resultsIBS$adj.P.Val<0.1

SignP0.05<-subset(resultsIBS, resultsIBS$ibs_cont_sig0.05==TRUE); dim(SignP0.05)
Sign_q<-subset(resultsIBS, resultsIBS$fdr_0.1==TRUE); dim(Sign_q)

save(SignP0.05,Sign_q, file="sigPlasma.Rda")
write.table(SignP0.05, file="IbsMirPlasmap05.csv", sep=",", col.names=NA)
write.table(Sign_q, file="IbsMirPlasma_q0.1.csv", sep=",", col.names=NA)

table1<-SignP0.05[,c(1,4,5)]


table4

anySig<-rbind(table1,table2, table3,table4)
anySig$comparison<-c(rep("IBSvsHC",3),rep("IBSDvsHC",6),rep("IBSCvsHC",4),rep("IBSDvsIBSC",10))
write.table(anySig, "anySigPlasmaMir.csv", sep=",", col.names=NA)

#####################################################################################





