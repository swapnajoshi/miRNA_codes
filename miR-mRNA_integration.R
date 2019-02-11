```
  # output: pdf_document
```
  
#   title: "Pancreatic cancer UCLA"
# author: "Swapna Mahurkar-Joshi"
# date: "May 20, 2016"
# output: word_document


```{r, message=FALSE} 
setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/Panc_can_ucla/")
# source("https://bioconductor.org/biocLite.R")
# biocLite("FDb.InfiniumMethylation.hg19")
library("FDb.InfiniumMethylation.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("ggplot2")
library("matlab")
library("heatmap.plus")

# library("rtracklayer")
# library("ggbio")
# library("GenomicFeatures")

func.list<-list()

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

func.list<-list()
func.list$vlookup<-function(val, df, col){
  df[df[1] == val, col][1]
}
sort.data.frame <- function(x, key, ...) {
  if (missing(key)) {
    rn <- rownames(x)
    if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
    x[order(rn, ...), , drop=FALSE]
  } else {
    x[do.call("order", c(x[key], ...)), , drop=FALSE]
  }
}

Wilcox.test.paired <- function(x, s1, s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  wx.out <- wilcox.test(x1, x2, alternative = "two.sided",paired = TRUE)
  p.Wilcox <- as.numeric(wx.out$p.value)
  return(p.Wilcox)
}
setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/Panc_can_ucla/")
# ```



# ```{r} 
# differentially methylated genes 
# cancer vs nor
load("hm450_panCan_files3.Rda")
mean.can <- apply(can, 1, mean, na.rm = T)
mean.nor <- apply(nor1, 1, mean, na.rm = T)
mean.difference.canvsnor <- mean.can - mean.nor
datMat<-cbind(can, nor1); dim(datMat)
rawp3 <- apply(datMat, 1, Wilcox.test, s1 = colnames(can), s2 = colnames(nor1))
adjp3 <- p.adjust(rawp3, method = "BH")
wilcox.can.nor<-matrix(NA, nrow=dim(datMat)[1], ncol=1)
row.names(wilcox.can.nor)<-row.names(datMat)
wilcox.can.nor<-as.data.frame(wilcox.can.nor)
wilcox.can.nor$p_raw<-rawp3
wilcox.can.nor$p_bh_adj<-adjp3
wilcox.can.nor$mean.can<-mean.can
wilcox.can.nor$mean.nor<-mean.nor
wilcox.can.nor$meandiff<-mean.difference.canvsnor
wilcox.can.nor$abs_meandiff<-abs(mean.difference.canvsnor)
wilcox.can.nor$V1<-NULL
can.nor.sig<-subset(wilcox.can.nor, wilcox.can.nor$p_bh_adj<=0.05); dim(can.nor.sig)
can.nor.md.p.sig<-subset(can.nor.sig, can.nor.sig$abs_meandiff>=0.2); dim(can.nor.md.p.sig) # 8611 6

# Annotate
can.nor.md.p.sig.annot<-merge(can.nor.md.p.sig, annot, by="row.names")
row.names(can.nor.md.p.sig.annot)<-can.nor.md.p.sig.annot[,1]
can.nor.md.p.sig.annot <- as.data.frame(can.nor.md.p.sig.annot)

save(datMat,wilcox.can.nor,can.nor.sig, can.nor.md.p.sig.annot, file="hm450_panCan_files4.Rda")
write.table(can.nor.md.p.sig.annot, file="can.nor.md.p.sig.annot.csv", sep=",")
# ```

# ```{r} 
#Starburst plot
load("all.probes.results.rda")
merged.md.ge <- merge(can.nor.md.p.sig.annot, all.probes.results, by.x = "nearestGeneSymbol", by.y = "gene.symbols1"); dim(merged.md.ge) # 11537   33
# some Cgs are dduplicated
merged.md.ge1 <- merged.md.ge[order(merged.md.ge[,2]),]
merged.md.ge1$dup <- duplicated(merged.md.ge1[,2])
merged.md.ge1 <- subset(merged.md.ge1 , merged.md.ge1$dup == FALSE); dim(merged.md.ge1)
rownames(merged.md.ge1)<-merged.md.ge1[,2]
merged.md.ge1 <- merged.md.ge1[ ,c(1,3:11, 27, 28, 29, 30, 31, 33)];  dim(merged.md.ge1) # 5823   16
merged.md.ge2 <- subset(merged.md.ge1, merged.md.ge1$adj.P.Val <=0.05); dim(merged.md.ge2); # 1567 16

#add a column for log 10 if md beta value is less than 0 and -log 10 if md beta is >=0

merged.md.ge2$log10adjp.me<-NA

for (i in 1: 1567)
{
  # print(cat("\n",i))
  if (merged.md.ge2[i,6] < 0)
  {
    merged.md.ge2[i,17]<-log(merged.md.ge2[i,3],10)
  }
  else if (merged.md.ge2[i,6] > 0)
  {
    merged.md.ge2[i,17]<--log(merged.md.ge2[i,3],10)
  }
}

#add a column for log 10 if expression fc value is less than 0 and -log 10 if expression fc is >=0

merged.md.ge2$log10adjp.ge<-NA

for (i in 1: 1567)
{
  # print(cat("\n",i))
  
  if (merged.md.ge2[i,16] < 0)
  {
    merged.md.ge2[i,18]<-log(merged.md.ge2[i,15],10)
  }
  else if (merged.md.ge2[i,16] > 0)
  {
    merged.md.ge2[i,18]<- -log(merged.md.ge2[i,15],10)
  }
}

merged.md.ge2$sig.logp.md.ge <- merged.md.ge2$log10adjp.me >= 1.3 &  merged.md.ge2$log10adjp.ge <= -1.3

p1 <- ggplot(merged.md.ge2,
             aes(log10adjp.me, log10adjp.ge, color = sig.logp.md.ge))

p1.1 <- p1 +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(size = 1.35) + theme_bw() +
  scale_colour_manual(values = c('FALSE' = "slateblue4",'TRUE' = "red4")) +
  #		coord_cartesian(xlim =  c(-13.0, 30.0), ylim = c(-8, 6)) +
  
  geom_vline(xintercept = c(-0.2, 0.2),
             linetype = 2, size = 0.25, color = "darkgray") +
  
  geom_hline(yintercept = c(1, -1),
             linetype = 2, size = 0.25, color = "darkgray") 
#		opts(legend.position = "none")
png("Starburst_can_nor1.3p-value.png")
p1.1
dev.off()

silenced.starburst<-subset(merged.md.ge2,merged.md.ge2$sig.logp.md.ge=="TRUE"); dim(silenced.starburst) # 555 19
save(merged.md.ge, merged.md.ge2, silenced.starburst,file="starburst_pancan.Rda")
write.table(silenced.starburst, file = "silenced.starburst.csv", sep = ",", col.names = NA)

# heatmaps of hnf4a
# hnf4a probes on 450k
dat.dm <- cbind(can, nor1)
hnf.all <- dat.dm[row.names(dat.dm)%in%row.names(annot[annot$nearestGeneSymbol == "HNF4A",]),]; dim(hnf.all)


## loaction of the probes

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
names(anno)
annot.hnf <- anno[ grep("HNF4A",anno$UCSC_RefGene_Name),]
annot.hnf <- annot.hnf[row.names(hnf.all),]
# significant probes [1] "cg06640637" "cg20848979" "cg23792485" "cg27420224"
word.list <- list(letters[1:4], letters[1:5], letters[1:2], letters[1:6])
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))

location <- strsplit(as.character(annot.hnf$UCSC_RefGene_Group),"[; ]+",fixed = FALSE)
n.obs <- sapply(location, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(location, "[", i = seq.max))
annot.hnf$location <- mat[,1]
match(row.names(hnf.all), row.names(annot.hnf))
hnf.all.can <- hnf.all[,colnames(hnf.all)%in%colnames(can)]; dim(hnf.all.can)
hnf.all.nor <- hnf.all[,colnames(hnf.all)%in%colnames(nor1)]; dim(hnf.all.nor)

# ```

```{r Visualization of differentially methylated and differentially expressed genes; Starburst plot}
p1.1

```

```{r #Heatmap}
  
  loc <- function(x) {
    if (x == "1stExon") "black"
    else if (x == "TSS1500") "black"
    else if (x == "TSS200") "black"
    else if (x == "Body") "grey"
    else if (x == "3'UTR") "grey"
  }
  col.j<-jet.colors(75)
  vlookup <- function(val, df, col){df[df[1] == val,col]}
  annot.hnf<- as.data.frame(annot.hnf)
  annot.hnf$ILMNID <-row.names(annot.hnf)
  annot.hnf <- annot.hnf[,c(16,1:15)]
  rc1<-apply(as.matrix(rownames(hnf.all.can)[1:23]),1, vlookup, annot.hnf,16)
  cc.col6<-unlist(lapply(rc1, loc))
  mat<-cbind(cc.col6,cc.col6)
  
  
  
  hv2<-heatmap.plus(as.matrix(hnf.all.can),na.rm=TRUE,scale="none",  
                    RowSideColor= mat,
                    #		ColSideColors=col.mat,
                    col=col.j,
                    zlim=c(0,1),
                    Colv=T,Rowv=T,
                    cexRow=1,cexCol=1,
                    # keep.dendro = NA,
                    main = paste("pancreatic_cancer", dim(hnf.all.can)[1],"Probes; ",dim(hnf.all.can)[2],"Samples"),
                    labRow=row.names(hnf.all.can)
  )
  
  hv2
  
  row.order <- hv2$rowInd
  #col.order <- colnames(tum.dic.pr1[ ,hv1$colInd])
  hnf.all.nor<-hnf.all.nor[row.order,]
  
  heatmap.plus(as.matrix(hnf.all.nor),na.rm=TRUE,scale="none",  
               #		RowSideColor=col.mat,
               #		ColSideColors=col.mat,
               col=col.j,
               zlim=c(0,1),
               Colv=T,Rowv=NA,
               cexRow=1,cexCol=1,
               # keep.dendro = NA,
               main = paste("Adjacent_normals", dim(hnf.all.nor)[1],"Probes; ",dim(hnf.all.nor)[2],"Samples"),
               labRow=row.names(hnf.all.nor)
  )
  
  

  
  
  