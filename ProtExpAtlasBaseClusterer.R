library(UniProt.ws)
library(data.table)
library(RCurl)

ProtExpAtlasBaseClusterer <- setClass("ProtExpAtlasBaseClusterer", 
                          slots= c(
                            experimentID = "character",
                            typeOfExp = "character"
                            ),
			  contains="EuclideanClusterer"
                          )

asNumeric <- function(x) as.numeric(x)
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.character)],   
                                                   asNumeric))


setMethod("retrieveFeatures",signature(object="ProtExpAtlasBaseClusterer"), function(object) {
  # We need first to translate proteins into genes for gene expression atlas
  keysProts<-object@proteins
  ktype<-"UNIPROTKB"
  columns<-c("ENSEMBL")
  res.dt <- data.table(select(UniProt.ws, keysProts, columns, ktype),key=c("ENSEMBL"))
  # Now we retrieve the data for the experiment, this depends on the type
  # of experiment: "baseline results" (where a matrix is downloaded
  # from https://www.ebi.ac.uk/gxa/experiments/<experimentID>.tsv ),
  # "differential expression results" 
  # and "data used for differential expression results"
  urlPrefix<-"https://www.ebi.ac.uk/gxa/experiments/"
  url<-paste(urlPrefix,experimentID,".tsv",sep="")
  d<-getURL(url)
  expData<-read.table(text=d,comment.char="#",sep="\t",header=T)
  data.table(expData[,c(-2)],key=c("Gene.ID"))->expData.dt
  expData.dt[res.dt]->expreForProts
  expreForProts[!is.na(expreForProts$Gene.ID),]->expreForProts
  expreForProts$rn<-row.names(expreForProts)
  expreForProts$UNIPROTKB[duplicated(expreForProts$UNIPROTKB)]<-paste(expreForProts$UNIPROTKB[duplicated(expreForProts$UNIPROTKB)],expreForProts$rn[duplicated(expreForProts$UNIPROTKB)],sep=".")
  as.data.frame(expreForProts)->expreForProts
  row.names(expreForProts)<-expreForProts$UNIPROTKB
  expreForProts[,!names(expreForProts) %in% c("ENSEMBL","Gene.ID","rn","UNIPROTKB")]->expreForProts
  expreForProtsDF<-as.data.frame.matrix(expreForProts)
  expreForProtsDF[expreForProtsDF=='']<-0
  object@uniqueFeaturesTable<-factorsNumeric(expreForProtsDF)
  return(object)
})

setMethod("calculateDistances",signature(object="ProtExpAtlasBaseClusterer"), function(object,minFeaturePres=0) {
  #toDrop <- names(which(colSums(object@uniqueFeaturesTable)<minFeaturePres))
  #uniqTableFeatsAboveMin <- object@uniqueFeaturesTable[!(colnames(object@uniqueFeaturesTable) %in% toDrop)]
  # here we should use jaccard or other binary distance for presence/absence, not euclidean.
  # 1 if for jaccard
  dist(x=as.matrix(object@uniqueFeaturesTable),method="euclidean")->object@proteinDist
  dist(x=t(as.matrix(object@uniqueFeaturesTable)),method="euclidean")->object@featureDist
  # According to this paper http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6421371&tag=1 complete and ward are the best options for
  # binary data
  hclust(d=object@proteinDist,method="complete")->object@clustProts
  hclust(d=object@featureDist,method="complete")->object@clustFeatures
  return(object)
})
