#' @include EuclideanClusterer.R
NULL

# the taxonomy needs to be set for Uniprot.ws (default to human)

setClass("ProtExpAtlasBaseClusterer", 
                          slots= c(
                            experimentID = "character",
                            typeOfExp = "character",
                            multiOrgExp = "logical"
                            ),
                          prototype = prototype(multiOrgExp = F),
			                    contains="EuclideanClusterer"
                          )

#' An S4 object that handles clustering of proteins based on the proteins 
#' expression in the ArrayExpress baseline experiments.
#' 
#' @param proteins List of UniProt identifiers to produce analysis for.
#' @param experimentID The experiment identifier in ArrayExpress Atlas Baseline 
#'   repository
#' @param taxonID The NCBI Taxonomy / UniProt Taxonomy identifier for the 
#'   species. Defaults to human 9606.
#' @param multiOrgExp A boolean denoting whether the ArrayExpress Atlast 
#'   Baseline experiment is multi organism or not. Default to \code{FALSE}.
#'   
#' @export
ProtExpAtlasBaseClusterer <- function(proteins, experimentID, taxonID=9606, multiOrgExp=FALSE, ...) {
    new("ProtExpAtlasBaseClusterer",proteins=proteins,experimentID=experimentID,taxonID=taxonID,
        multiOrgExp=multiOrgExp)  
}


asNumeric <- function(x) as.numeric(x)
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.character)],   
                                                   asNumeric))

taxId2Label<-data.frame(taxId=c(9606,10090),label=c("Homo sapiens","Mus musculus"))

#' @export
setMethod("retrieveFeatures",signature(object="ProtExpAtlasBaseClusterer"), function(object) {
  # We need first to translate proteins into genes for gene expression atlas
  
  if(nrow(object@annotation)==0) {
    annotateProteins(object)->object
    #keysProts<-object@proteins
    #ktype<-"UNIPROTKB"
    #columns<-c("ENSEMBL")
    #if(object@taxonID!=9606) {
    #  taxId(UniProt.ws)<-object@taxonID
    #}
    #res.dt <- data.table(select(UniProt.ws, keysProts, columns, ktype),key=c("ENSEMBL"))
  }
  res.dt = data.table(object@annotation,key=c("ENSEMBL"))
  # Now we retrieve the data for the experiment, this depends on the type
  # of experiment: "baseline results" (where a matrix is downloaded
  # from https://www.ebi.ac.uk/gxa/experiments/<experimentID>.tsv ),
  # "differential expression results" 
  # and "data used for differential expression results"
  urlPrefix<-"https://www.ebi.ac.uk/gxa/experiments/"
  urlSuffix<-""
  if(object@multiOrgExp)
    urlSuffix<-paste("?serializedFilterFactors=ORGANISM:",taxId2Label$label[taxId2Label$taxId==object@taxonID],sep="")
  # some data sets might require additional information, otherwise they are obtained
  # by default for human:
  # http://www.ebi.ac.uk/gxa/experiments/E-GEOD-30352.tsv?serializedFilterFactors=ORGANISM:Mus%20musculus
  # http://www.ebi.ac.uk/gxa/experiments/E-MTAB-599.tsv?serializedFilterFactors=ORGANISM:Mus%20musculus
  url<-paste(urlPrefix,object@experimentID,".tsv",urlSuffix,sep="")
  d<-getURL(URLencode(url))
  expData<-read.table(text=d,comment.char="#",sep="\t",header=T)
  data.table(expData[,c(-2)],key=c("Gene.ID"))->expData.dt
  expData.dt[res.dt]->expreForProts
  expreForProts[!is.na(expreForProts$Gene.ID),]->expreForProts
  expreForProts$rn<-row.names(expreForProts)
  expreForProts$UNIPROTKB[duplicated(expreForProts$UNIPROTKB)]<-paste(expreForProts$UNIPROTKB[duplicated(expreForProts$UNIPROTKB)],expreForProts$rn[duplicated(expreForProts$UNIPROTKB)],sep=".")
  as.data.frame(expreForProts)->expreForProts
  row.names(expreForProts)<-expreForProts$UNIPROTKB
  expreForProts[,!names(expreForProts) %in% c("ENSEMBL","Gene.ID","rn","UNIPROTKB","GENES")]->expreForProts
  expreForProtsDF<-as.data.frame.matrix(expreForProts)
  expreForProtsDF[expreForProtsDF=='']<-0
  object@uniqueFeaturesTable<-factorsNumeric(expreForProtsDF)
  return(object)
})
