################################################################################
# function for qqplot 
################################################################################
qqplot.data <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  d <- data.frame(resids = vec)
  ggplot(d, aes(sample = resids)) + stat_qq() + geom_abline(slope = slope, 
                                                            intercept = int)
}


################################################################################
# Expressed genes: signal detected in 20% IBS or control samples
################################################################################
Expressed <- function(x) {rowSums(x[,grep("IBS", colnames(x))] > 0) >= 
dim(x)[2]/5 | rowSums(x[,grep("HC",colnames(x))] > 0) >= dim(x)[2]/5 }


################################################################################
# Fold Change: input x/y data and output foldchange values
################################################################################
foldChange <- function(x) {ifelse (x >= 0, round(2^x,2), round(-1*2^(-1*x),2))}


################################################################################
#Create a dataframe for number of DEgene
################################################################################

lmDf <- function(datMat,pData){
  TS = factor(pData,levels=levels(factor(pData)))
  design <- model.matrix(~0 + TS)
  colnames(design) <- levels(TS)
  fit <- lmFit(datMat, design)
  if (length(levels(factor(pData)))==2){
    cont.matrix <- makeContrasts(paste(paste(levels(factor(pData))[1],levels(factor(pData))[2],
                                             sep = "_"), paste(levels(factor(pData))[1],levels(factor(pData))[2],
                                                               sep = "-"), sep = "="), levels =  design)
  }
  else if (length(levels(factor(pData)))==3) {
    cont.matrix <- makeContrasts(paste(paste(levels(factor(pData))[1],levels(factor(pData))[2],
                                             sep = "_"), paste(levels(factor(pData))[1],levels(factor(pData))[2],
                                                               sep = "-"), sep = "="), ... = 
                                   paste(paste(levels(factor(pData))[1],levels(factor(pData))[3],
                                               sep = "_"), paste(levels(factor(pData))[1],levels(factor(pData))[3],
                                                                 sep = "-"), sep = "="), 
                                 paste(paste(levels(factor(pData))[2],levels(factor(pData))[3],
                                             sep = "_"), paste(levels(factor(pData))[2],levels(factor(pData))[3],
                                                               sep = "-"), sep = "="), levels =  design)
  }
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  ibs_fit2 <- eBayes(fit2)
  
  if (length(levels(factor(pData)))==2) {
    results = topTable(ibs_fit2, number=dim(datMat)[1], coef=1)
  }
  else if (length(levels(factor(pData)))==3) {
    results1 <- topTable (ibs_fit2, number=dim (datMat)[1], coef=1)
    colnames(results1) <- paste(colnames(results1), colnames(cont.matrix)[1], sep = "_")
    results2 <- topTable (ibs_fit2, number=dim(datMat)[1], coef=2)
    colnames(results2) <- paste(colnames(results2), colnames(cont.matrix)[2], sep = "_")
    results2 <- results2[row.names(results1),]
    results3 <- topTable(ibs_fit2, number=dim(datMat)[1], coef=3)
    colnames(results3) <- paste(colnames(results3), colnames(cont.matrix)[3], sep = "_")
    results3 <- results3[row.names(results2),]
    results <- cbind (results1, results2, results3) 
  }
  return(results)
}


################################################################################
# make a table of significant values
################################################################################
tableSig <- function(x)
{ a1 <- subset(x, x[,1] < 0.05)
a2 <- subset(x, x[,2] < 0.05)
a3 <- subset(x, x[,3] < 0.05)
a4 <- subset(x, x[,4] < 0.05) 
y <- rbind( round(a1,3), round(a2,3), round(a3,3), round(a4,3))
y <- as.data.frame(y)
z <- colnames(x[,1:4])
Comparison <- c(c(rep(z[1], dim(a1)[1])),rep(z[2], dim(a2)[1]), rep(z[3], dim(a3)[1]), rep(z[4], dim(a4)[1]))
y$Comparison <- Comparison
y <- y[order(y[,1],y[,2],y[,3],y[,4]),]
return(y)  
}
################################################################################
# Get the names of DE miRNAs with normalization that gives most significant 
# number and see if they overlap with any pubmed entry
################################################################################


library(RISmed)
res <- EUtilsSummary("nanoString, miRNA", type="esearch", db="pubmed", 
datetype='pdat', mindate=2000, maxdate=2015, retmax=500)
QueryCount(res)
summary(res)
t<-ArticleTitle(EUtilsGet(res))
typeof(t)
head(t,1)

QueryId(res)

records<- EUtilsGet(res)
class(records)
pubmed_data <- data.frame('Title'=ArticleTitle(records),
'Abstract'=AbstractText(records))
head(pubmed_data,1)

write.table(pubmed_data, file = "nanostring_miRNA_pubmed_abstracts.csv", 
            sep = ",")






################################################################################
# Function for correlation between two dataframes
################################################################################
cor2dfs <- function(x,y) {
  corList <- rcorr(as.matrix(x), type="spearman")
  corVal  <- corList$r
  adjP <- apply(corList$P,2,p.adjust)
  highCor <- ifelse(adjP <= 0.1, 1,0)
  corInt <- highCor[1:y, y+1:dim(x)[2]]
  corInt$sum1<-apply(corInt,1,sum)
  corInt1<-corInt[corInt$sum1>0,]
  corInt2<-apply(corInt1,2,sum)
  df2Int<-names(corInt2[corInt2>0])
  SigCor <- corInt2[,df2Int]
  return(SigCor)
}
# Usage
# df1Sel <- t(microbMat2.ibs[row.names(microbMat2.ibs)%in%
#                                  row.names(cor2dfs),]); dim(cor.microb)
# df2Sel <- t(mirna1.ibs[row.names(mirna1.ibs) %in% colnames(cor2dfs),]); 
# dim(cor.mirna)
# 
# data2<-as.data.frame(cbind(df1Sel,df2Sel));dim(data2)
################################################################################
# make ggplots scatter plot from 2 vectors and a color variable in a data frame 
################################################################################
plotAllScatter <- function(z) {
  P<-ggplot(z, aes(z[,x],z[,y])) + geom_point(size=2.5) 
  q <- P + geom_point(aes(color=as.factor(z$Group))) + xlab(colnames(z)[x]) + 
    ylab(colnames(z)[y])
  ggsave(q, filename = paste(colnames(z)[x],colnames(z)[y],".png"))
}
# usage
# for(x in 1:2) for(y in 3:4) {
#   plotAllScatter(data3)
# }
################################################################################
# make ggplots box plot from 2 vectors, (one numeric, one factor) in a dataframe 
################################################################################
plotAllBox <- function(z) {
  P <- ggplot(z, aes(as.factor(z[,x]), z[,y], fill=as.factor(z[,x]))) 
  q <- P + geom_boxplot() + xlab(colnames(z)[x]) +  ylab(colnames(z)[y])
  ggsave(q, filename = paste(colnames(z)[x],colnames(z)[y],".png"))
}
# usage
# for(y in 3:4) {
#   x = 'Group'
#   plotAllBox(data3)
# }
################################################################################
# Scatterplots with regression line with geom_smooth function
################################################################################
plotAllScatterRegr <- function(z) {
  P<-ggplot(z, aes(z[,x],z[,y])) + geom_point(size=2.5) 
  q <- P + geom_point(aes()) + xlab(colnames(z)[x]) + 
    ylab(colnames(z)[y]) + geom_smooth(aes(z[,x],z[,y]), method=lm, se=FALSE)
  ggsave(q, filename = paste(colnames(z)[x],colnames(z)[y],".png"))
}
# usage
# for(x in 1:2) for(y in 3:4) {
#   plotAllScatterRegr(data3)
# }
################################################################################
# heatmap.3 function
################################################################################
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = colorRampPalette(c("green", "black", "red")), 
                      #function(x)rev(heat.colors(x)),
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider 
            using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(colnames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors),
             las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), 
             colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), 
        axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 
                                rep(0, length(csep)), xright = csep + 0.5 + 
                                sepwidth[1], ytop = rep(ncol(x) + 1, csep), 
                              lty = 1, lwd = 1, col = sepcolor, 
                              border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, 
                              xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep)
                              - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

################################################################################
# Heatmap function; make a scaled dataframe and create the heatmap in object
################################################################################
heatmapZ <- function(x) {
  library(matlab)
  z <- t(apply(x, 1, scale))
  row.names(z)<-row.names(x)
  col.j<-jet.colors(75)
  hm1<-heatmap.3(as.matrix(z),na.rm=TRUE,scale="none", Colv=NA,Rowv=T,
                 cexRow=1,cexCol=1, main = paste("Heatmap", dim(z)[1],
                 "Probes; ",dim(z)[2],"Samples"), labRow = row.names(z), 
                 labCol = colnames(z))
  return(hm1)
}





