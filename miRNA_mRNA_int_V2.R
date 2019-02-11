```{r, message=FALSE}
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rcpp")

library(Hmisc)
library(heatmap.plus)
library(matlab)
# library(survival)
.libPaths()
#IBS vs HC DE's

setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/")
mRNA <- read.delim("mRNA_IBS_HC_May2016_Beth.csv", sep = ",", row.names = 1); dim(mRNA)

sig_mrna <- subset(mRNA, mRNA$P.Value <=0.05); dim(sig_mrna)
dat_mrna <- read.delim("Dat_mrna.csv", sep = ",", row.names = 1); dim(dat_mrna)
colnames(dat_mrna)[1:30] <- substr(colnames(dat_mrna)[1:30],3,7)
annot_mrna <- dat_mrna[,c(61:76)]
dat_mrna <- dat_mrna[,c(1:30)]; dim(dat_mrna); dat_mrna[1:5,1:5]

dat_mrna_sig <- dat_mrna[row.names(dat_mrna)%in%row.names(sig_mrna),]; dim(dat_mrna_sig)

sig_mir <- read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/Results/biopsy_miR/IbsMirbiopsy_SigRes.csv", sep = ",", row.names = 1); dim(sig_mir)
load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/Results/biopsy_miR/norm.data4.rda"); dim(norm.data4)
dat_mir <- norm.data4[row.names(norm.data4)%in%row.names(sig_mir),]; dim(dat_mir); dat_mir[1:5,1:5]
colnames(dat_mir) <- substr(colnames(dat_mir), 1, 5)
dat_mrna1 <- dat_mrna_sig[,colnames(dat_mrna_sig)%in%colnames(dat_mir)]
dat_mir1 <- dat_mir[ ,colnames(dat_mir)%in%colnames(dat_mrna1)]
dat_mrna1 <- dat_mrna1[,colnames(dat_mir1)]
match(colnames(dat_mrna1),colnames(dat_mir1))
# correlation between miRNA and mRNA
library(mixOmics)
color<-color.jet
par(mfrow = c(1,2));
png("mixOmics.png", height=2000, width=3000, res=300)
liver.splsda <- imgCor(t(dat_mir1), t(dat_mrna1), X.var.names = TRUE,
                       Y.var.names = TRUE,
                       sideColors = TRUE,
                       interactive.dev = TRUE,
                       main = TRUE)
# color, row.cex, col.cex,symkey, keysize,
# xlab, ylab, margins, lhei, lwid)
dev.off()

data1<-cbind(t(dat_mir1),t(dat_mrna1)); dim(data1)
library(ellipse)
ctab <- cor(data1)
corMat<-round(ctab, 2); dim(corMat)

###### pvalue based filtering 
library(Hmisc)
ptab <- rcorr(as.matrix(data1))
corp<-ptab$P
adj.p<-matrix(p.adjust(corp, method = "BH"), nrow = 2094, ncol = 2094, byrow = FALSE)
row.names(adj.p) <- row.names(corp)
colnames(adj.p) <- colnames(corp) 
p.table <- as.data.frame(as.table(adj.p)) # matric to cor pairs
sig.p <- subset(p.table, p.table[,3]<=0.05); dim(sig.p)
sig.p <-as.data.frame(sig.p )
sig.p$unique <- ifelse(sig.p[,1]!=sig.p[,2],1,0)
sig.p1 <- subset(sig.p, sig.p$unique == 1); dim(sig.p1) #682204  
sig.p$unique <- NULL

sig.p$Pair <- 0
for ( i in 1:682204) {
  cat(paste("row",i))
  if (substr(sig.p[i,1],1,3) == substr(sig.p[i,2],1,3)) {
    sig.p[i,4] <- "within_mRNA_mRNA" }
  else if (substr(sig.p[i,1],1,3) != substr(sig.p[i,2],1,3)) {
    sig.p[i,4] <- "mRNA-miRNA correlation" 
  }
}

miR_mrna <- subset(sig.p, sig.p$Pair == "mRNA-miRNA correlation"); dim(miR_mrna) # 2884
corc<-ptab$r
c.table <- as.data.frame(as.table(corc))
miR_mrna$cor1<-c.table[row.names(c.table)%in%row.names(miR_mrna),3]
miR_mrna_neg <- subset(miR_mrna, miR_mrna$cor1 < 0); dim(miR_mrna_neg) # 1684 5
miR_mrna_pos <- subset(miR_mrna, miR_mrna$cor1 > 0); dim(miR_mrna_pos) # 1200 5
write.csv(miR_mrna_neg, file = "miR_mrna_neg.csv", col.names = NA)
# annotate mRNAs
# miR_mrna$GeneSymbol <- annot_mrna[row.names(annot_mrna)%in%miR_mrna[,1],4]
miR_mrna<- as.data.frame(miR_mrna)
colnames(miR_mrna)[1] <-"ID"
annot_mrna$ID <- row.names(annot_mrna)
miR_mrna_annot <- merge(miR_mrna, annot_mrna, by = "ID")

save (miR_mrna_annot, file= "miR_mrna_annot.rda")
write.table(miR_mrna_annot, file = "miR_mrna_annot.csv", sep = ",", col.names = NA)
save (ptab, file = "ptab.rda")
save (sig.p, file = "sig.p.rda")
load("miR_mrna_annot.rda")
load("sig.p.rda")

# validated targets of miR_mrna_pos in the list of correlated genes
miR_mrna_annot$Var2 <- gsub("\\|0.*","",miR_mrna_annot$Var2)
# miR_mrna_annot$miRNA_id <- matrix(unlist(strsplit(as.character(miR_mrna$Var2),"\\|")),nrow=1442, ncol = 2, byrow = TRUE)[,1]
val1<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/miRNA_new_analyses_IBS9232015/miRTarBase_MTI_validated targets_human.csv", sep = ",")
val1$pair.val <- paste(val1$miRNA,val1$Target.Gene, sep = " ")
miR_mrna_annot$datPair <- paste(miR_mrna_annot$Var2,miR_mrna_annot$GeneSymbol, sep = " ")

write.table(miR_mrna_annot, file = "miR_mrna_annot.csv", sep = ",", col.names = NA)


val.mirna.mrna <- val1[which(val1$pair.val%in%miR_mrna_annot$datPair),]; dim(val.mirna.mrna)
# 15   12
write.table(val.mirna.mrna, file = "val.mirna.mrna.csv", sep = ",")


# Heatmap of top 50 negatively correlated genes
row.names(dat_mir1) <- gsub("\\|0.*","",row.names(dat_mir1))

# miR_mrna_annot_val <- miR_mrna_annot[which(miR_mrna_annot$GeneSymbol %in% val.mirna.mrna1$Target.Gene),]; dim(miR_mrna_annot_val)
# miR_mrna_annot_val1 <- miR_mrna_annot_val[,c(1,2,3,4,5,9)]
miRNA1 <- merge(miR_mrna_annot, dat_mir1, by.x = "Var2", by.y = "row.names"); dim(miRNA1)
miRNA1 <- miRNA1[order(miRNA1$cor1),]
miRNA2 <- miRNA1[c(1:25,1418:1442),]; dim(miRNA2)
miRNA2 <- miRNA2[order(miRNA2$Var2, miRNA2$ID),]
miRNA3 <- miRNA2[,c(1,23:51)]; dim(miRNA3)

mRNA1 <- merge(miR_mrna_annot, dat_mrna1, by.x = "ID", by.y = "row.names"); dim(mRNA1)
mRNA1 <- mRNA1[order(mRNA1$cor1),]
mRNA2 <- mRNA1[c(1:25,1418:1442),]; dim(mRNA2)
mRNA2 <- mRNA2[order(mRNA2$Var2, mRNA2$ID),]
mRNA3 <- mRNA2[,c(9,23:51)]; dim(mRNA3)

save(mRNA3, file="mRNA3.rda")
save(miRNA3, file = "miRNA3.rda")

hm.mirna <- t(apply(miRNA3[,2:30],1,scale))
row.names(hm.mirna)<-row.names(miRNA3[,2:30])
colnames(hm.mirna)<-colnames(miRNA3[,2:30])
phenoDat2 <- phenoDat2[colnames(hm.mirna),]
hm.mirna <- as.data.frame(hm.mirna)
col.means1 <- colMeans(hm.mirna)
col.means2 <- as.data.frame(cbind(col.means1, as.character(phenoDat2$Condition)))
col.means2$col.means1 <- as.numeric(as.character(col.means2$col.means1))
col.means3 <- col.means2[rev(order(col.means2$V2, col.means2$col.means1)),]
phenoDat2 <- phenoDat2[row.names(col.means3),]
hm.mirna <- hm.mirna[,row.names(phenoDat2)]
match(row.names(phenoDat2), colnames(hm.mirna))

phenoDat2 <- phenoDat2[,c(2,1,3:8)]
#cc.sex<-apply(as.matrix(colnames(ibs)[1:25]),1,vlookup,var.ibs,8)
#cc.col5<-unlist(lapply(cc.sex, sex))

BowelHabit<-function(mol.biol) {
  if (mol.biol=="IBSC") "purple" #red
  else if (mol.biol=="IBSD") "green" #yellow
  else if (mol.biol=="HC") "orange" #Black
  #	else if (mol.biol=="6") "yellow" #Black
  else if (mol.biol=="NA") "grey" #Black
}

col.j<-jet.colors(75)
vlookup <- function(val, df, col){df[df[1] == val,col]}
cc.bh<-apply(as.matrix(colnames(hm.mirna)[1:29]),1,vlookup,phenoDat2,3)
cc.col6<-unlist(lapply(cc.bh, BowelHabit))
mat<-cbind(cc.col6,cc.col6)

png("hm_mirna_frompairs_top50_col.final.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(hm.mirna),na.rm=TRUE,scale="none",  
                  # RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=T,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("pancreatic_cancer", dim(hm.mirna)[1],"Probes; ",dim(hm.mirna)[2],"Samples"),
                  labRow=miRNA3$Var2
)
dev.off()


# hm.mrna <- merge(dat_mrna1, a1, by.x = "row.names", by.y = "ID", all.y = TRUE)
# hm.mrna <- hm.mrna[order(hm.mrna$Freq),]
# # row.names(hm.mrna) <- hm.mrna[,1]
# hm.mrna1 <- hm.mrna[,c(2:30,38)]

hm.mrna<-t(apply(mRNA3[,2:30],1,scale))
row.names(hm.mrna)<-row.names(mRNA3[,2:30])
colnames(hm.mrna)<-colnames(mRNA3[,2:30])
hm.mrna <- hm.mrna[,rownames(phenoDat2)]
match(row.names(phenoDat2), colnames(hm.mrna))
# phenoDat2 <- phenoDat2[,c(2,1,3:8)]
#cc.sex<-apply(as.matrix(colnames(ibs)[1:25]),1,vlookup,var.ibs,8)
#cc.col5<-unlist(lapply(cc.sex, sex))

BowelHabit<-function(mol.biol) {
  if (mol.biol=="IBSC") "purple" #red
  else if (mol.biol=="IBSD") "green" #yellow
  else if (mol.biol=="HC") "orange" #Black
  #	else if (mol.biol=="6") "yellow" #Black
  else if (mol.biol=="NA") "grey" #Black
}
col.j<-jet.colors(75)
vlookup <- function(val, df, col){df[df[1] == val,col]}

row.order <- hv2$rowInd
#col.order <- colnames(tum.dic.pr1[ ,hv1$colInd])
hm.mrna<-hm.mrna[rev(row.order),]
# col.order <- hv2$colInd
# hm.mrna2<-hm.mrna2[,col.order]
cc.bh<-apply(as.matrix(colnames(hm.mrna)[1:29]),1,vlookup,phenoDat2,3)
cc.col6<-unlist(lapply(cc.bh, BowelHabit))
mat<-cbind(cc.col6,cc.col6)


png("hm_mrna_frompairs_top50_col.final.png", bg="white", res=400, width=3000, height=3000)
hv3<-heatmap.plus(as.matrix(hm.mrna),na.rm=TRUE,scale="none",  
                  #		RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=NA,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("ibs", dim(hm.mrna)[1],"Probes; ",dim(hm.mrna)[2],"Samples"),
                  labRow=mRNA3$GeneSymbol[row.order]
)
dev.off()

#############################################################################################



install.packages('mRNAnormSet', type = "source", repos = NULL)

IBS<-as.factor(pheno$IBS)
mod1 <- model.matrix(~ 0 + IBS) 
colnames(mod)=c("HC","IBS")

##############################################################################

val1.mRNA <- val1[val1$Target.Gene %in% sig_mrna$GeneSymbol,]
val1.mRNA.sel <- val1.mRNA[val1.mRNA$miRNA == c("hsa-miR-548ah-5p", "hsa-miR-363-3p","hsa-miR-219a-5p", "hsa-miR-188-5p", "hsa-miR-106b-5p"),] 
write.table(val1.mRNA.sel, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/val1.mRNA.sel.csv", sep = ",", col.names = NA)

######################
# preparing files for bicluster analysis
setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/")
mRNA <- read.delim("mRNA_IBS_HC_May2016_Beth.csv", sep = ",", row.names = 1); dim(mRNA)

sig_mrna <- subset(mRNA, mRNA$P.Value <=0.05); dim(sig_mrna)
dat_mrna <- read.delim("Dat_mrna.csv", sep = ",", row.names = 1); dim(dat_mrna)
colnames(dat_mrna)[1:30] <- substr(colnames(dat_mrna)[1:30],3,7)
annot_mrna <- dat_mrna[,c(61:76)]
dat_mrna <- dat_mrna[,c(1:30,64,73)]; dim(dat_mrna); dat_mrna[1:5,1:5]
dat_mrna_sig <- dat_mrna[row.names(dat_mrna)%in%row.names(sig_mrna),]; dim(dat_mrna_sig)
write.table(dat_mrna_sig, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/Bicluster_analysis/mRNAData.txt", sep = '\t')

sig_mir <- read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/Results/biopsy_miR/IbsMirbiopsy_SigRes.csv", sep = ",", row.names = 1); dim(sig_mir)
load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/Results/biopsy_miR/norm.data4.rda"); dim(norm.data4)
dat_mir <- norm.data4[row.names(norm.data4)%in%row.names(sig_mir),]; dim(dat_mir); dat_mir[1:5,1:5]
colnames(dat_mir) <- substr(colnames(dat_mir), 1, 5)
dat_mrna1 <- dat_mrna_sig[,colnames(dat_mrna_sig)%in%colnames(dat_mir)]
dat_mir1 <- dat_mir[ ,colnames(dat_mir)%in%colnames(dat_mrna1)]
dat_mrna1 <- dat_mrna1[,colnames(dat_mir1)]
match(colnames(dat_mrna1),colnames(dat_mir1))
row.names(dat_mir1) <- matrix(unlist(strsplit(row.names(dat_mir1), "\\|")), nrow = 34, ncol = 2, byrow = TRUE)[,1]
write.table(dat_mir1, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/Bicluster_analysis/dataMir.txt", sep = '\t')

#####################################

c("LIPG", "BCL2", "STC1", "BCHE")

c("TMEM8B", "ARMC10", "SHC1")
