source("https://bioconductor.org/biocLite.R")
biocLite("survival")

library(Hmisc)
library(heatmap.plus)
library(matlab)

#IBS vs HC DE's

setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/")
mRNA <- read.delim("mRNA_IBS_HC_May2016_Beth.csv", sep = ",", row.names = 1); dim(mRNA)

sig_mrna <- subset(mRNA, mRNA$P.Value <=0.05); dim(sig_mrna)
dat_mrna <- read.delim("Dat_mrna.csv", sep = ",", row.names = 1); dim(dat_mrna)
colnames(dat_mrna)[1:30] <- substr(colnames(dat_mrna)[1:30],3,7)
annot_mrna <- dat_mrna[,c(61:76)]
dat_mrna <- dat_mrna[,c(1:30)]; dim(dat_mrna); dat_mrna[1:5,1:5]

dat_mrna_sig <- dat_mrna[row.names(dat_mrna)%in%row.names(sig_mrna),]; dim(dat_mrna_sig)

sig_mir <- read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/TissueMir/IbsMirtissuep05.csv", sep = ",", row.names = 1); dim(sig_mir)
load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/TissueRCC/norm.data4.rda"); dim(norm.data4)
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

# miRNA target pair plots
phenoDat1 <- read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/TissueMir/phenoDat1.csv", sep = ",")
row.names(phenoDat1) <- substr(phenoDat1[,1],1,5)
phenoDat2 <- phenoDat1[row.names(phenoDat1)%in%row.names(data1),]
phenoDat2 <- phenoDat2[row.names(data1),]
match(row.names(data1), row.names(phenoDat2))
data2 <- cbind(as.data.frame(data1), as.character(phenoDat2[,3]))
colnames(data2) <- gsub("-","_",colnames(data2))
colnames(data2) <- gsub("\\|","_",colnames(data2))
colnames(data2)[2095] <- "condition"

miR_mrna_annot[miR_mrna_annot$GeneSymbol %in% "MAPK9",]
p1 <- ggplot(data2, aes(x = ASHGA5P052543, y = hsa_miR_363_3p_0.005))
p2 <- p1 + geom_point(aes(color=factor(condition), size = 3)) + scale_color_manual(values = c("orange", "purple", "green"))
ggsave(p2, file = "MAPK9_hsa-miR-363-3p.png",width = 15, height = 15, units = "cm")

miR_mrna_annot[miR_mrna_annot$GeneSymbol %in% "PIK3CB",]
p1 <- ggplot(data2, aes(x = ASHGA5P004286, y = hsa_miR_363_3p_0.005))
p2 <- p1 + geom_point(aes(color=factor(condition), size = 3)) + scale_color_manual(values = c("orange", "purple", "green"))
ggsave(p2, file = "PIK3CB_hsa-miR-363-3p.png",width = 15, height = 15, units = "cm")

miR_mrna_annot[miR_mrna_annot$GeneSymbol %in% "PRKCG",]
p1 <- ggplot(data2, aes(x = ASHGA5P034761, y = hsa_miR_4531_0))
p2 <- p1 + geom_point(aes(color=factor(condition), size = 3)) + scale_color_manual(values = c("orange", "purple", "green"))
ggsave(p2, file = "PRKCG_hsa-hsa_miR_4531_0.png",width = 15, height = 15, units = "cm")

miR_mrna_annot[miR_mrna_annot$GeneSymbol %in% "CAMK2B",]
p1 <- ggplot(data2, aes(x = ASHGA5P013679, y = hsa_miR_450a_5p_0))
p2 <- p1 + geom_point(aes(color=factor(condition), size = 3)) + scale_color_manual(values = c("orange", "purple", "green"))
ggsave(p2, file = "CAMK2B_hsa_miR_450a_5p_0.png", width = 15, height = 15, units = "cm")


## TBD check predicted or validated genes-mirna pairs within these pairs.
#  mirwalk pairs
mw_pairs <- read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/miRNA_mRNA_integration/3utr-region-mirwalk-interactions-pos1.csv", sep = ",")
mw_pairs <- as.data.frame(mw_pairs)
mw_pairs$Pre_pair <-NA
for (i in 1: dim(mw_pairs)[1]){
  mw_pairs[i,11] <- paste (as.character(mw_pairs[i,1]), as.character(mw_pairs[i,3], sep = "_"))
}

miR_mrna_annot <- as.data.frame(miR_mrna_annot)
# miR_mrna_annot$miRNA <- matrix(unlist(strsplit(as.character(miR_mrna_annot[,2]), "\\|")), nrow = 1442, ncol = 2, byrow = FALSE)[1,]
miR_mrna_annot$miRNA_name <- gsub("-","_",miR_mrna_annot$Var2)
miR_mrna_annot$miRNA_name <- gsub("\\|","_",miR_mrna_annot$miRNA_name)

miR_mrna_annot$datPair <- NA
for (i in 1: dim(miR_mrna_annot)[1]){
  miR_mrna_annot[i,23] <- paste (as.character(miR_mrna_annot[i,22]), as.character(miR_mrna_annot[i,9], sep = "_"))
}

predAnddatPairs <- merge (miR_mrna_annot, mw_pairs, by.x = "datPair", by.y = "Pre_pair")

write.table(predAnddatPairs, file = "predAnddatPairs.csv", sep = ",", col.names =NA)


## perform pathway analyses with these # show some plots

# heatmap of top 50 mrna-mrna pais

ibs<-subset(phenoDat1, phenoDat1$Condition != "HC")
#hc<-rownames(subset(phenoDat1, phenoDat1$status=="HC"))
#ibsc<-rownames(subset(phenoDat1, phenoDat1$Condition=="IBSC"))
#ibsd<-rownames(subset(phenoDat1, phenoDat1$Condition=="IBSD"))

ibs<-as.data.frame(ibs)
ibs$id<-row.names(ibs)
ibs<-ibs[,c(9,1:8)]
data2.t <- t(data2)
mrna1<-merge(val.mirna.mrna1, data2.t, by.x = "Target.Gene", by.y = "GeneSymbol")
micrna1<-merge(miR_mrna_annot, data2.t, by.x = "Var2", by.y = "row.names")

# hm.mRNA <- t(data2[,colnames(data2) %in% miR_mrna_annot$ID])
# hm.mrna <-hm.mrna[miR_mrna_annot$ID,]
# hm.microRNA <- t(data2[,c(1:25)])

miR_mrna_annot1 <- miR_mrna_annot[order(miR_mrna_annot$Freq),]


a1 <-  miR_mrna_annot1[1:50,]
row.names(dat_mir1) <- gsub("-","_",row.names(dat_mir1))
row.names(dat_mir1) <- gsub("\\|","_",row.names(dat_mir1))
hm.mirna <- merge(dat_mir1, a1, by.x = "row.names", by.y = "miRNA_name", all.y = TRUE)
hm.mirna <- hm.mirna[order(hm.mirna$Freq),]
hm.mirna1 <- hm.mirna[,c(2:30,1)]
hm.mirna2 <- t(apply(hm.mirna1[,1:29],1,scale))
row.names(hm.mirna2)<-row.names(hm.mirna1[,1:29])
colnames(hm.mirna2)<-colnames(hm.mirna1[,1:29])
# phenoDat2 <- phenoDat2[colnames(hm.mirna2),]
hm.mirna2 <- as.data.frame(hm.mirna2)
col.means1 <- colMeans(hm.mirna2)
col.means2 <- as.data.frame(cbind(col.means1, as.character(phenoDat2$Condition)))
col.means2$col.means1 <- as.numeric(as.character(col.means2$col.means1))
col.means3 <- col.means2[rev(order(col.means2$V2, col.means2$col.means1)),]
phenoDat2 <- phenoDat2[row.names(col.means3),]
hm.mirna2 <- hm.mirna2[,row.names(phenoDat2)]
match(row.names(phenoDat2), colnames(hm.mirna2))



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
cc.bh<-apply(as.matrix(colnames(hm.mirna2)[1:29]),1,vlookup,phenoDat2,3)
cc.col6<-unlist(lapply(cc.bh, BowelHabit))
mat<-cbind(cc.col6,cc.col6)

png("hm_mirna_frompairs_top50_col.sorted.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(hm.mirna2),na.rm=TRUE,scale="none",  
                  # RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=T,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("pancreatic_cancer", dim(hm.mirna2)[1],"Probes; ",dim(hm.mirna2)[2],"Samples"),
                  labRow=substr(hm.mirna1$Row.names,4,20)
)
dev.off()


hm.mrna <- merge(dat_mrna1, a1, by.x = "row.names", by.y = "ID", all.y = TRUE)
hm.mrna <- hm.mrna[order(hm.mrna$Freq),]
# row.names(hm.mrna) <- hm.mrna[,1]
hm.mrna1 <- hm.mrna[,c(2:30,38)]


hm.mrna2<-t(apply(hm.mrna1[,1:29],1,scale))
row.names(hm.mrna2)<-row.names(hm.mrna1[,1:29])
colnames(hm.mrna2)<-colnames(hm.mrna1[,1:29])
hm.mrna2 <- hm.mrna2[,rownames(phenoDat2)]
match(row.names(phenoDat2), colnames(hm.mrna2))
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
hm.mrna2<-hm.mrna2[row.order,]
# col.order <- hv2$colInd
# hm.mrna2<-hm.mrna2[,col.order]
cc.bh<-apply(as.matrix(colnames(hm.mrna2)[1:29]),1,vlookup,phenoDat2,3)
cc.col6<-unlist(lapply(cc.bh, BowelHabit))
mat<-cbind(cc.col6,cc.col6)


png("hm_mrna_frompairs_top50_col.sorted.png", bg="white", res=400, width=3000, height=3000)
hv3<-heatmap.plus(as.matrix(hm.mrna2),na.rm=TRUE,scale="none",  
                  #		RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=NA,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("ibs", dim(hm.mrna2)[1],"Probes; ",dim(hm.mrna2)[2],"Samples"),
                  labRow=hm.mrna1$GeneSymbol
)
dev.off()



# validated targets of miR_mrna_pos in the list of correlated genes
miR_mrna_annot$miRNA_id <- matrix(unlist(strsplit(as.character(miR_mrna$Var2),"\\|")),nrow=1442, ncol = 2, byrow = TRUE)[,1]
val1<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/miRNA_new_analyses_IBS9232015/miRTarBase_MTI_validated targets_human.csv", sep = ",")
val.mirna.mrna <- val1[which(val1[,2]%in%miR_mrna_annot$miRNA_id),]; dim(val.mirna)
# 4661   11
val.mirna.mrna1 <- val.mirna.mrna[which(val.mirna.mrna[,4]%in%miR_mrna_annot$GeneSymbol),]; dim(val.mirna.mrna1)
write.table(val.mirna.mrna1, file = "val.mirna.mrna1.csv", sep = ",")

require("mirIntegrator")
data(kegg_pathways)
data(mirTarBase)
kegg_pathways <- kegg_pathways[18:20] #delete this for augmenting all pathways.
augmented_pathways <- integrate_mir(kegg_pathways, mirTarBase)
head(augmented_pathways)


# require(graph)
# >   require(ROntoTools)
# >   data(GSE43592_mRNA)
# >   data(GSE43592_miRNA)
# >   data(augmented_pathways)
# >   data(names_pathways)
# >   lfoldChangeMRNA <- GSE43592_mRNA$logFC
# >   names(lfoldChangeMRNA) <- GSE43592_mRNA$entrez
# >   lfoldChangeMiRNA <- GSE43592_miRNA$logFC
# >   names(lfoldChangeMiRNA) <- GSE43592_miRNA$entrez
# >   keggGenes <- unique(unlist( lapply(augmented_pathways,nodes) ) )
# >   interGMi <- intersect(keggGenes, GSE43592_miRNA$entrez)
# >   interGM <- intersect(keggGenes, GSE43592_mRNA$entrez)
# >   ## For real-world analysis, nboot should be >= 2000
#   >   peRes <- pe(x= c(lfoldChangeMRNA, lfoldChangeMiRNA ),
#                   +               graphs=augmented_pathways, nboot = 200, verbose = FALSE)
# >   message(paste("There are ", length(unique(GSE43592_miRNA$entrez)),
#                   +                 "miRNAs meassured and",length(interGMi),
#                   +                 "of them were included in the analysis."))
# >   message(paste("There are ", length(unique(GSE43592_mRNA$entrez)),
#                   +                 "mRNAs meassured and", length(interGM),
#                   +                 "of them were included in the analysis."))
# >   summ <- Summary(peRes)
# >   rankList <- data.frame(summ,path.id = row.names(summ))
# >   tableKnames <- data.frame(path.id = names(names_pathways),names_pathways)
# >   rankList <- merge(tableKnames, rankList, by.x = "path.id", by.y = "path.id")
# >   rankList <- rankList[with(rankList, order(pAcc.fdr)), ]
# >   head(rankList)