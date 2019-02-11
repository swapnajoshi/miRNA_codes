```{r, message=FALSE} 
library(matlab)
library(heatmap.plus)
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Beth_lnRNA/")
# readin array star and ns raw data
ma <- read.delim("MAraw_50by28.csv", sep = ",", row.names = 1); dim(ma)
ns <- read.delim("NSraw_50by28.csv", sep = ",", row.names = 1); dim(ns)
mRNA <- read.delim("mRNA_IBS_HC_May2016_Beth.csv", sep = ",", row.names = 1); dim(mRNA)

gs1 <- mRNA[row.names(mRNA)%in%row.names(ma),]
gs1 <- gs1[row.names(ma),]
match(row.names(gs1), row.names(ma))
row.names(ma) <- gs1$GeneSymbol

# match(row.names(ma), row.names(ns))
match(colnames(ma), colnames(ns))
ma <- ma[rownames(ns),]
match(colnames(ma), colnames(ns))

match(rownames(ma), rownames(ns))
ma <- ma[,colnames(ns)]
match(colnames(ma), colnames(ns))

# Read in phenotype data and match the rownames
phenoDat1 <- read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Beth_lnRNA/nanostring_manifest.csv", sep = ",", row.names = 2)
phenoDat2 <- phenoDat1[row.names(phenoDat1)%in%colnames(ma),]
colnames(phenoDat2)[8] <- "BH"
# phenoDat2 <- phenoDat2[,colnames(ma)]
match(colnames(ma), row.names(phenoDat2))

# Scale the data
hm.ma <- t(apply(ma[,1:28],1,scale))
row.names(hm.ma)<-row.names(ma)
colnames(hm.ma)<-colnames(ma)
hm.ma <- as.data.frame(hm.ma)
col.means1 <- colMeans(hm.ma)
col.means2 <- as.data.frame(cbind(col.means1, as.character(phenoDat2$BH)))
col.means2$col.means1 <- as.numeric(as.character(col.means2$col.means1))
col.means3 <- col.means2[rev(order(col.means2$V2, col.means2$col.means1)),]
phenoDat2 <- phenoDat2[row.names(col.means3),]
hm.ma <- hm.ma[,row.names(phenoDat2)]
match(row.names(phenoDat2), colnames(hm.ma))
phenoDat2 <- as.data.frame(phenoDat2)
phenoDat2$ID <- row.names(phenoDat2)
phenoDat2 <- phenoDat2[,c(9,1:8)]
phenoDat2$BH <- gsub("N", "HC", phenoDat2$BH)
#cc.sex<-apply(as.matrix(colnames(ibs)[1:25]),1,vlookup,var.ibs,8)
#cc.col5<-unlist(lapply(cc.sex, sex))

BowelHabit<-function(mol.biol) {
  if (mol.biol=="C") "purple" #red
  else if (mol.biol=="D") "green" #yellow
  else if (mol.biol=="HC") "orange" #Black
  #	else if (mol.biol=="6") "yellow" #Black
  else if (mol.biol=="NA") "grey" #Black
}

col.j<-jet.colors(75)
vlookup <- function(val, df, col){df[df[1] == val,col]}

# Bar for hm
cc.bh<-apply(as.matrix(colnames(hm.ma)[1:28]),1,vlookup,phenoDat2,9)
cc.col6<-unlist(lapply(cc.bh, BowelHabit))
mat<-cbind(cc.col6,cc.col6)

png("hm_ma_arraystar_raw_1.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(hm.ma),na.rm=TRUE,scale="none",  
                  # RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=T,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("arraystar_raw", dim(hm.ma)[1],"Probes; ",dim(hm.ma)[2],"Samples"),
                  labRow=row.names(hm.ma)
)
dev.off()




# hm for ns

hm.ns<-t(apply(ns,1,scale))
row.names(hm.ns)<-row.names(ns)
colnames(hm.ns)<-colnames(ns)
hm.ns <- hm.ns[,rownames(phenoDat2)]
match(row.names(phenoDat2), colnames(hm.ns))

# keep the row order 
row.order <- hv2$rowInd
hm.ns<-hm.ns[rev(row.order),]  # heatmap.plus reverses the order
# col.order <- hv2$colInd
# hm.mrna2<-hm.mrna2[,col.order]
cc.bh<-apply(as.matrix(colnames(hm.ns)[1:28]),1,vlookup,phenoDat2,9)
cc.col6<-unlist(lapply(cc.bh, BowelHabit))
mat<-cbind(cc.col6,cc.col6)


png("hm_ns_raw_1.png", bg="white", res=400, width=3000, height=3000)
hv3<-heatmap.plus(as.matrix(hm.ns),na.rm=TRUE,scale="none",  
                  #		RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=NA,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("ibs", dim(hm.ns)[1],"Probes; ",dim(hm.ns)[2],"Samples"),
                  labRow = rev(row.names(hm.ns))
)
dev.off()


# For markdown
par(mfrow = c(1,2))
hv2<-heatmap.plus(as.matrix(hm.ma),na.rm=TRUE,scale="none",  
                  # RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=T,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("arraystar_raw", dim(hm.ma)[1],"Probes; ",dim(hm.ma)[2],"Samples"),
                  labRow=row.names(hm.ma)
)


hv3<-heatmap.plus(as.matrix(hm.ns),na.rm=TRUE,scale="none",  
                  #		RowSideColor=mat,
                  ColSideColors=mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=NA,Rowv=NA,
                  cexRow=1,cexCol=1,
                  # keep.dendro = NA,
                  main = paste("ibs", dim(hm.ns)[1],"Probes; ",dim(hm.ns)[2],"Samples"),
                  labRow = rev(row.names(hm.ns))
)
