
# TODO: NanoStringNorm
# Biopsy Tissue 
# Author: SwapnaJoshi 4182016

###############################################################################

```{r, message=FALSE} 

#Funcitons and libraries:
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

library(NanoStringNorm)
library(ggplot2)
library(limma)
library(Biobase)

####################################################################################
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/tissueRCC/")
rccfiles<-read.markup.RCC(rcc.pattern="*.RCC|*.rcc",  exclude=NULL, include = NULL, nprobes=-1)
#Returns a list with two components. The first is the header information which contains sample IDs
#and diagnostic information on the quality of the samples. The second is the count data and can be
#directly used in the input to NanoStringNorm.

tissueMirDf<-rccfiles$x

phenoDat<-colnames(tissueMirDf)[4:48]

phenoDat1<-as.data.frame(matrix(unlist(strsplit(phenoDat,"_")), nrow=45, ncol=6, byrow=TRUE))
row.names(phenoDat1)<-phenoDat
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

for (i in 1:45)
{
  
  png(paste(i,"posCon.png"))	
  barplot(posC1[i,]) # to look at individual samples
  dev.off()
}


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



# positive normalization factor very high for A6145_HC_F_F09_Tis_Apr15 , remove and rerun

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


png("boxplot_mirna_raw-dropposcontrol.png", height=1000, width=1000, res=200)
boxplot(raw.data1,horizontal=FALSE,axes=FALSE,cex=0.2,outline=FALSE,xlab="x-axis label", ylab="y-axis label",las=1)
axis(1,cex.axis=0.5)
axis(2,cex.axis=0.5) 
dev.off()

boxplot(raw.data1,horizontal=FALSE,axes=FALSE,cex=0.2,outline=FALSE,xlab="x-axis label", ylab="y-axis label",las=1)
axis(1,cex.axis=0.5)
axis(2,cex.axis=0.5) 

norm.data3<-norm.data2[,row.names(phenoDat1)]
mean2.sam<-apply(norm.data3,2,mean)


png("boxplot_mirna_norm-droppositivecontrol.png", height=1000, width=1000, res=200)
boxplot(norm.data3,horizontal=FALSE,axes=TRUE,cex=0.2,outline=FALSE,xlab="x-axis label", ylab="y-axis label",las=1)
axis(1,cex.axis=0.5)
axis(2,cex.axis=0.5) 
dev.off()

boxplot(norm.data3,horizontal=FALSE,axes=TRUE,cex=0.2,outline=FALSE,xlab="x-axis label", ylab="y-axis label",las=1)
axis(1,cex.axis=0.5)
axis(2,cex.axis=0.5) 

save(norm.data3, file="norm.data3.mean.tissue.rda")
save(raw.data1, file="raw.data1.tissue.rda")
save(phenoDat1, file="phenoDat1.tissue.rda")

dim(norm.data3)

# [1] 800  43


count.dat<-matrix(NA, nrow=dim(norm.data3)[1],ncol=dim(norm.data3)[2])
for (i in 1:dim(norm.data3)[1])
  for (j in 1: dim(norm.data3)[2])
  {
    if (norm.data3[i,j]>0)#> log(25,2) [1] 4.32
      
    { count.dat[i,j]<-1}
    
    else if (norm.data3[i,j]<=0)
    { count.dat[i,j]<-0}
  }

row.names(count.dat)<-row.names(norm.data3)
colnames(count.dat)<-colnames(norm.data3)

##### 50% of the controls or 50 % of patients samples must have atlest one count, 50% pat (29) =12, 50% cont (14) = 7

count.dat<-as.data.frame(count.dat)
count.dat$Expressed<-0
for (i in 1:dim(count.dat)[1]) {
  if(apply(count.dat[i,grep("IBS",colnames(count.dat))],1,sum) >=12) {
  count.dat[i,44] <- 1
  }
    else if(apply(count.dat[i,grep("HC",colnames(count.dat)),],1,sum) >=6) {
    count.dat[i,44] <- 1
    }
    else(count.dat[i,44] <- 0)
  
   }
  
probesInt<-subset(count.dat, count.dat$Expressed==1)
norm.data4<-norm.data3[row.names(norm.data3)%in%row.names(probesInt),]

dim(norm.data4)
# [1] 333  44
save(norm.data4, file="norm.data4.rda")
write.table(norm.data4, file="norm.data4.csv", sep=",")

## q q plots
norm.data4.t<-t(norm.data4)
norm.data4.t<-as.data.frame(norm.data4.t)
norm.data4.t$Condition<-phenoDat1$status

#parametric or non parametric? generate a q-q plot
meanDat<-phenoDat1
meanDat$mean1<-mean2.sam

library(ggplot2)
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
qqplot1<-qqplot.data(meanDat$mean1)
ggsave(qqplot1, file="qqplottissue.mean.IBS1.png") ## 
qqplot1

# ibs<-rownames(subset(phenoDat1, phenoDat1$status=="IBS"))
# hc<-rownames(subset(phenoDat1, phenoDat1$status=="HC"))
# ibsc<-rownames(subset(phenoDat1, phenoDat1$Condition=="IBSC"))
# ibsd<-rownames(subset(phenoDat1, phenoDat1$Condition=="IBSD"))

###################################################################################

# Limma
setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/TissueMir/")

# load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/tissueRCC/phenoDat1.tissue.rda")

# load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/tissueRCC/norm.data4.rda")

exprs1<-norm.data4

colnames(exprs1)<-substr(colnames(exprs1),1,5)

## limma
#  import phenotype data

row.names(phenoDat1)<-phenoDat1[,1]
match(row.names(phenoDat1),colnames(exprs1))


#expressionset

TS = factor(pData(phenoData)$id,levels=c("IBS","HC"))
design <- model.matrix(~0+TS)
v <- voom(exprs(exprDat),design,plot=TRUE)
colnames(design) <- levels(TS)
fit <- lmFit(v, design)

#contrast matrix
cont.matrix <- makeContrasts(IBS_HC = IBS-HC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
ibs_fit2 <- eBayes(fit2)

png("IbsMirtissueVolca1_voom_limma.png")
volcanoplot(ibs_fit2)
dev.off()
volcanoplot(ibs_fit2)

#Differentially expressed genes	
topTable(ibs_fit2, number=10, coef=1)
resultsIBS=topTable(ibs_fit2, number=333)
save(resultsIBS, file="resultsIBS.rda")


resultsIBS$ibs_cont_sig0.05<-resultsIBS$P.Value<0.05
SignP0.05<-subset(resultsIBS, resultsIBS$ibs_cont_sig0.05==TRUE)

# SignP0.05

save(SignP0.05, file="IbsMirtissuep05.rda")
write.table(SignP0.05, file="IbsMirtissuep05.csv", sep=",", col.names=NA)

table1<-SignP0.05[,c(1,4,5)]
table1; dim(table1) #25 3

##################################################################################

#####IBSD

row.names(phenoDat1)<-phenoDat1[,1]

match(row.names(phenoDat1),colnames(exprs1))
pData <- as.matrix(phenoDat1[phenoDat1[,2]!="IBSC",2]); row.names(pData)=row.names(phenoDat1[phenoDat1[,2]!="IBSC",]); colnames(pData)="id"
pData<-(as.data.frame(pData))

phenoData <- new("AnnotatedDataFrame", data=pData)
exprs1<-norm.data4
colnames(exprs1)<-substr(colnames(exprs1),1,5)
exprs2<-exprs1[,colnames(exprs1)%in%row.names(pData)]
match(row.names(pData),colnames(exprs2))

exprDat<-new("ExpressionSet", exprs=as.matrix(exprs2), phenoData=phenoData)

#expressionset
TS = factor(pData(phenoData)$id,levels=c("IBSD","HC"))
design <- model.matrix(~0+TS)
v <- voom(exprs(exprDat),design,plot=TRUE)
colnames(design) <- levels(TS)
fit <- lmFit(v, design)

#contrast matrix
cont.matrix <- makeContrasts(IBSD_HC = IBSD-HC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
ibs_fit2 <- eBayes(fit2)

png("IbsDMirtissueVolca1.png")
volcanoplot(ibs_fit2)
dev.off()
volcanoplot(ibs_fit2)

#Differentially expressed genes	
topTable(ibs_fit2, number=10, coef=1)

resultsIBSD=topTable(ibs_fit2, number=333)
save(resultsIBSD, file="resultsIBSD.rda")
resultsIBSD$ibs_cont_sig0.05<-resultsIBSD$P.Value<0.05

SignP0.05<-subset(resultsIBSD, resultsIBSD$ibs_cont_sig0.05==TRUE)

save(SignP0.05, file="IbsDMirtissuep0.05.rda")
write.table(SignP0.05, file="IbsDMirtissuep0.05.csv", sep=",", col.names=NA)
table2<-SignP0.05[,c(1,4,5)]
table2; dim(table2) # [1] 16 3
#################################################################################

#####IBSC
row.names(phenoDat1)<-phenoDat1[,1]

match (row.names(phenoDat1),colnames(exprs1))
pData <- as.matrix(phenoDat1[phenoDat1[,2]!="IBSD",2]); row.names(pData)=row.names(phenoDat1[phenoDat1[,2]!="IBSD",]); colnames(pData)="id"
pData<-(as.data.frame(pData))

phenoData <- new("AnnotatedDataFrame", data=pData)
exprs1<-norm.data4
colnames(exprs1)<-substr(colnames(exprs1),1,5)
exprs2<-exprs1[,colnames(exprs1)%in%row.names(pData)]
match(row.names(pData),colnames(exprs2))

exprDat<-new("ExpressionSet", exprs=as.matrix(exprs2), phenoData=phenoData)

#expressionset
TS = factor(pData(phenoData)$id,levels=c("IBSC","HC"))
design <- model.matrix(~0+TS)
v <- voom(exprs(exprDat),design,plot=TRUE)
colnames(design) <- levels(TS)
fit <- lmFit(v, design)

#contrast matrix
cont.matrix <- makeContrasts(IBSC_HC = IBSC-HC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

ibs_fit2 <- eBayes(fit2)

png("IbsCMirtissueVolca1.png")
volcanoplot(ibs_fit2)
dev.off()
volcanoplot(ibs_fit2)

#Differentially expressed genes	
topTable(ibs_fit2, number=10, coef=1)
resultsIBSC=topTable(ibs_fit2, number=333)
save(resultsIBSC, file="resultsIBSC.rda")

resultsIBSC$ibs_cont_sig0.05<-resultsIBSC$P.Value<0.05
SignP0.05 <- subset(resultsIBSC, resultsIBSC$ibs_cont_sig0.05 == TRUE)
save(SignP0.05, file = "IbsCMirtissuep0.05.rda")
write.table(SignP0.05, file = "IbsCMirtissuep0.05.csv", sep = ",", col.names =    NA)
table3 <- SignP0.05[,c(1,4,5)]
table3; dim(table3) # 19  3



####################################################################################
# IBSD vs IBSC

row.names(phenoDat1)<-phenoDat1[,1]

match (row.names(phenoDat1),colnames(exprs1))
pData <- as.matrix(phenoDat1[phenoDat1[,2]!="HC",2]); row.names(pData)=row.names(phenoDat1[phenoDat1[,2]!="HC",]); colnames(pData)="id"
pData<-(as.data.frame(pData))

phenoData <- new("AnnotatedDataFrame", data=pData)
exprs1<-norm.data4
colnames(exprs1)<-substr(colnames(exprs1),1,5)
exprs2<-exprs1[,colnames(exprs1)%in%row.names(pData)]
match(row.names(pData),colnames(exprs2))

exprDat<-new("ExpressionSet", exprs=as.matrix(exprs2), phenoData=phenoData)

#expressionset
TS = factor(pData(phenoData)$id,levels=c("IBSD","IBSC"))
design <- model.matrix(~0+TS)
v <- voom(exprs(exprDat),design,plot=TRUE)
colnames(design) <- levels(TS)
fit <- lmFit(v, design)

#contrast matrix
cont.matrix <- makeContrasts(IBSD_IBSC = IBSD-IBSC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

ibs_fit2 <- eBayes(fit2)

png("IbsDCMirtissueVolca1.png")
volcanoplot(ibs_fit2)
dev.off()
volcanoplot(ibs_fit2)

#Differentially expressed genes	
topTable(ibs_fit2, number=10, coef=1)

resultsIBSCD=topTable(ibs_fit2, number=333)
save(resultsIBSCD, file="resultsIBSCD.rda")
resultsIBSCD$ibs_cont_sig0.05<-resultsIBSCD$P.Value<0.05

SignP0.05 <- subset(resultsIBSCD, resultsIBSCD$ibs_cont_sig0.05 == TRUE)

save(SignP0.05, file = "IbsCDMirtissuep0.05.rda")
write.table(SignP0.05, file = "IbsCDMirtissuep0.05.csv", sep = ",", col.names =
              NA)
table4 <- SignP0.05[,c(1,4,5)]
table4; dim(table4) # 6 3

anySig<-as.data.frame(rbind(table1, table2, table3, table4)); dim(anySig)

anySig$comparison<-c(rep("IBSvsHC",25),rep("IBSDvsHC",16),rep("IBSCvsHC",19),rep("IBSDvsIBSC",6))
write.table(anySig, "anySigtissueMir.csv", sep=",", col.names=NA)

############################ code to select the most efficient method

mat <- matrix(NA, nrow = 24, ncol = 4)
mat[,1] <- rep(c("sum","geo.mean"),12)
mat[,2] <- rep(c("mean","mean.2sd","max"),8)
mat[,3] <- rep(c("total.sum","low.cv.geo.mean","top.mean","top.geo.mean"),6)
colnames(mat) <- c("CodeCount", "Background", "SampleContent", "nExprGenes")
row.names(mat) <- c(1:24)

for (i in 1:24){
        mat[i,4] <- dim(na.omit(ifelse(NanoStringNorm (tissueMir,anno = NA,
                                                       header = NA,
                          Probe.Correction.Factor ='none',CodeCount = mat[i,1],
                          Background = mat[i,2], SampleContent =	mat[i,3],	
                          OtherNorm =	'none',
                          round.values = TRUE,
                          is.log = FALSE,
                          take.log = TRUE,
                          return.matrix.of.endogenous.probes = TRUE,
                          traits = NA,
                          predict.conc = FALSE,
                          verbose = TRUE,
                          genes.to.fit,
                          genes.to.predict,
                          guess.cartridge = FALSE) > 0, 
                          NanoStringNorm (tissueMir,anno = NA,header =NA,
                          Probe.Correction.Factor ='none',CodeCount = mat[i,1],
                          Background = mat[i,2], SampleContent =	mat[i,3],	
                          OtherNorm =	'none',
                          round.values = TRUE,
                          is.log = FALSE,
                          take.log = TRUE,
                          return.matrix.of.endogenous.probes = TRUE,
                          traits = NA,
                          predict.conc = FALSE,
                          verbose = TRUE,
                          genes.to.fit,
                          genes.to.predict,
                          guess.cartridge = FALSE), NA)))[1]
        
      }
################################################################
# trials with different DE methods

# Differentially expressed genes

```{r}
phenoDat2 <- phenoDat1[-grep("A5813",row.names(phenoDat1)),]
phenoDat2 <- phenoDat2[-grep("A6145",row.names(phenoDat2)),]

phenoDat2 <- as.data.frame(phenoDat1[colnames(biopsyNormDat),])
pData <- as.matrix(phenoDat2[,c(2,3,8)]);row.names(pData)=row.names(phenoDat2); 
colnames(pData) <- c("BowelHabit","Sex","Dx")
pData <- as.data.frame(pData)
phenoData <- new("AnnotatedDataFrame", data=pData)
exprDat <- new("ExpressionSet", exprs=as.matrix(biopsyNormDat), 
               phenoData=phenoData)

TS = factor(pData(phenoData)$Dx,levels=c("IBS","HC"))
design <- model.matrix(~0 + TS)
v <- voom(exprs(exprDat),design,plot=TRUE)
colnames(design) <- levels(TS)
fit <- lmFit(v, design)

#contrast matrix
cont.matrix <- makeContrasts(IBS_HC = IBS-HC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
ibs_fit2 <- eBayes(fit2)

png("IbsMirtissueVolca1_voom_limma.png")
volcanoplot(ibs_fit2)
dev.off()
volcanoplot(ibs_fit2)

#Differentially expressed genes	
topTable(ibs_fit2, number=30, coef=1)
resultsIBS=topTable(ibs_fit2, number=333)
save(resultsIBS, file="resultsIBS.rda")


resultsIBS$ibs_cont_sig0.05<-resultsIBS$P.Value<0.05
SignP0.05<-subset(resultsIBS, resultsIBS$ibs_cont_sig0.05==TRUE)


```

```{r}
library(edgeR)
tissueMir2 <- tissueMir1
row.names(tissueMir2) <- tissueMir1[,2]
tissueMir2 <- tissueMir2 [tissueMir2$CodeClass == "Endogenous1", -c(1,2,3)]
d <- DGEList(counts=as.matrix(tissueMir2))
phenoDat2 <- phenoDat1[-grep("A6145",row.names(phenoDat1)),]

phenoDat2 <- as.data.frame(phenoDat2[colnames(tissueMir2),])
pData <- as.matrix(phenoDat2[,c(2,3,8)]);row.names(pData)=row.names(phenoDat2); 
colnames(pData) <- c("BowelHabit","Sex","Dx")
pData <- as.data.frame(pData)
phenoData <- new("AnnotatedDataFrame", data=pData)

TS = factor(pData(phenoData)$Dx,levels=c("IBS","HC"))
design <- model.matrix(~0 + TS)
v <- voom(d, design)
colnames(design) <- levels(TS)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(IBS_HC = IBS-HC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
ibs_fit2 <- eBayes(fit2)

v <- voom(d, design, normalize="quantile")
v <- voom(d, design, normalize="cyclicloess")
```

```{r}

tissueMir2 <- tissueMir1
row.names(tissueMir2) <- tissueMir1[,2]
tissueMir2 <- tissueMir2 [tissueMir2$CodeClass == "Endogenous1", -c(1,2,3)]
count.dat<-matrix(NA, nrow=dim(tissueMir2)[1],ncol=dim(tissueMir2)[2])
for (i in 1:dim(tissueMir2)[1])
  for (j in 1: dim(tissueMir2)[2])
  {
    if (tissueMir2[i,j]>25)#> log(25,2) [1] 4.32

    { count.dat[i,j]<-1}

    else if (tissueMir2[i,j]<=25)
    { count.dat[i,j]<-0}
  }
row.names(count.dat)<-row.names(tissueMir2)
colnames(count.dat)<-colnames(tissueMir2)
##### 50% of the controls or 50 % of patients samples must have atlest one count, 50% pat (29) =15, 50% cont (14) = 7

count.dat<-as.data.frame(count.dat)
count.dat$Expressed<-0
for (i in 1:dim(count.dat)[1]) {
  if(apply(count.dat[i,grep("IBSC",colnames(count.dat))],1,sum) >=7) {
  count.dat[i,44] <- 1
  }
    else if(apply(count.dat[i,grep("IBSD",colnames(count.dat)),],1,sum) >=7) {
    count.dat[i,44] <- 1
    }
   else if(apply(count.dat[i,grep("HC",colnames(count.dat)),],1,sum) >=7) {
    count.dat[i,44] <- 1
   }
    else(count.dat[i,44] <- 0)

   }

probesInt<-subset(count.dat, count.dat$Expressed==1)
biopsyDat<-tissueMir2[row.names(tissueMir2)%in%row.names(probesInt),]



row.names(pData1)=row.names(phenoDat2)
colnames(pData1)[1] <- "BowelHabit"
phenoData <- new("AnnotatedDataFrame", data=pData1)




# Dx
TS = factor(pData(phenoData)$Dx,levels=c("IBS","HC"))
design <- model.matrix(~0 + TS)
colnames(design) <- levels(TS)

d <- DGEList(counts=as.matrix(biopsyDat))
# v <- voom(d, design)
# v <- voom(d, design, normalize="quantile")
v <- voom(d, design, normalize="cyclicloess")


fit <- lmFit(v, design)
cont.matrix <- makeContrasts(IBSC_HC = IBSC-HC,IBSD_HC = IBSD-HC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
ibs_fit2 <- eBayes(fit2)
topTable(ibs_fit2, number=50, coef=2)

resultsIBS=topTable(ibs_fit2, number=dim(biopsyDat)[1])
SignP0.05<-subset(resultsIBS, resultsIBS$P.Value<0.05)

# vMethod can be "none" "quantile" "cyclicloess" 
lmRes <- function(x, pData,vMethod) { 
  TS <- factor(pData(phenoData)$Dx,levels=c("IBS","HC")) 
  design <- model.matrix(~0 + TS)
  colnames(design) <- levels(TS)
  d <- DGEList(counts=as.matrix(x))
  v <- voom(d, design, normalize=vMethod, lib.size = NULL)

}


```

```{r}
nano.norm <- NanoStringNorm(tissueMir,anno = NA,header = NA,
                          Probe.Correction.Factor ='none',CodeCount ='none',
                          Background = 'mean',SampleContent =	'none',	
                          OtherNorm =	'none',
                          round.values = TRUE,
                          is.log = FALSE,
                          take.log = FALSE,
                          return.matrix.of.endogenous.probes = TRUE,
                          traits = NA,
                          predict.conc = FALSE,
                          verbose = TRUE,
                          genes.to.fit,
                          genes.to.predict,
                          guess.cartridge = FALSE
)
nano.norm <- as.data.frame(nano.norm)
nano.norm$mean1 <- apply(nano.norm,1,mean)
nano.norm1 <- subset(nano.norm, nano.norm$mean1 > 10)
nano.norm1$mean1 <- NULL

phenoDat2 <- phenoDat1[row.names(phenoDat1)%in%colnames(nano.norm1),]
pData1 <- as.data.frame(phenoDat2[colnames(nano.norm1),c(2,3,8)])
library(DESeq)
Dx <- as.factor(pData1$Dx)
BH <- as.factor(pData1$Condition)
cds <- newCountDataSet(nano.norm1, BH)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds ) 
plotDispEsts <- function( cds ){
plot(rowMeans(counts( cds, normalized=TRUE) ), fitInfo(cds)$perGeneDispEsts, 
     pch = '.', log ="xy", ylab ="dispersion",xlab="mean of normalized counts")
  xg = 10^seq( -.5, 5, length.out=300 )
lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
plotDispEsts(cds)
results <- nbinomTest(cds, 'IBSC', 'HC', pvals_only = FALSE, eps=NULL )
subset(results, results$pval <=0.05)

```

###################################################END
