# 
# library("hgu95av2.db")
# ids <- getPMID(affys,"hgu95av2")
# ids <- unlist(ids,use.names=FALSE)
# ids <- unique(ids[!is.na(as.numeric(ids))])
# length(ids)
# 
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("annotate")
# library("annotate")
# x <- pubmed("nanoString")
# a <- xmlRoot(x)
# numAbst <- length(xmlChildren(a))
# numAbst


library(RISmed)
res <- EUtilsSummary("nanoString, miRNA", type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2015, retmax=500)
QueryCount(res)
summary(res)
t<-ArticleTitle(EUtilsGet(res))
typeof(t)
head(t,1)

QueryId(res)

records<- EUtilsGet(res)
class(records)
pubmed_data <- data.frame('Title'=ArticleTitle(records),'Abstract'=AbstractText(records))
head(pubmed_data,1)

write.table(pubmed_data, file = "nanostring_miRNA_pubmed_abstracts.csv", sep = ",")

