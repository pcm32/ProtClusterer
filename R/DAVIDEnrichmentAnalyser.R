#' @include GeneralClusterer.R
NULL

DAVIDEnrichmentAnalyser <- setClass("DAVIDEnrichmentAnalyser",
				    slots= c(
					     email = "character",
					     annotCategories = "vector"	
					     )
				    )

#' \code{DAVIDEnrichmentAnalyser} is an S4 object that handles uses the DAVID 
#' site for enrichment analysis.
#' 
#' @param email The registered email to use the DAVID remote API. This email has
#'   to be registered with DAVID.
#' @param annotCategories A vector of characters with the desired categories for
#'   enrichment analysis. 
setGeneric("DAVIDEnrichmentAnalyser",function(email,annotCatergories) standardGeneric("DAVIDEnrichmentAnalyser"))
# categories are available here: http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html#approved_list

getGO_ALLCategories <- function() { c("GOTERM_BP_ALL","GOTERM_MF_ALL","GOTERM_CC_ALL") }
getPathwaysCategories <- function() { c("BBID","BIOCARTA","EC_NUMBER","KEGG_COMPOUND","KEGG_PATHWAY","KEGG_REACTION") }

getDiseaseCategories <- function() { c("GENETIC_ASSOCIATION_DB_DISEASE","OMIM_DISEASE") }

setGeneric("enrichmentAnalysisForClusters",function(object,generalClusterer,clusterNumbers,prefix,...) standardGeneric("enrichmentAnalysisForClusters"))
setMethod("enrichmentAnalysisForClusters",signature(object="DAVIDEnrichmentAnalyser",generalClusterer="GeneralClusterer",clusterNumbers="vector",prefix="character"),function(object,generalClusterer,clusterNumbers,prefix="") {
  david<-DAVIDWebService$new(email=object@email)
  #setAnnotationCategories(david, c("GOTERM_BP_ALL","GOTERM_MF_ALL","GOTERM_CC_ALL"))
  setAnnotationCategories(david, object@annotCategories)
  slicedProteins = generalClusterer@proteinGroup
  
  dir.create(file.path(getwd(), "enrichmentAnalysis"), showWarnings = FALSE)
  
  for(g in clusterNumbers) {
    group<-names(slicedProteins[slicedProteins == g])
    result<-addList(david,inputIds=group,listName=paste("group",g,sep=""),idType="UNIPROT_ACCESSION")
    getFunctionalAnnotationChart(david)->funcAnnChart
    write.table(x=funcAnnChart,file=paste("enrichmentAnalysis/enrichmentGroup_",prefix,g,".xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    funcAnnChartDAG<-1
    try(DAVIDGODag(funcAnnChart,pvalueCutoff = 0.01,"BP")->funcAnnChartDAG)
    if(!is.numeric(funcAnnChartDAG)) {  
      pdf(file=paste("enrichmentAnalysis/enrichmentDiagramBP_",prefix,g,".pdf",sep=""),paper="a4")
      plotGOTermGraph(g=goDag(funcAnnChartDAG),r=funcAnnChartDAG, max.nchar=40, node.shape="ellipse")      
      dev.off()
    }
    
    funcAnnCHartDAG<-1
    try(DAVIDGODag(funcAnnChart,pvalueCutoff = 0.01,"MF")->funcAnnChartDAG)
    if(!is.numeric(funcAnnChartDAG)) {            
      pdf(file=paste("enrichmentAnalysis/enrichmentDiagramMF_",prefix,g,".pdf",sep=""),paper="a4")
      plotGOTermGraph(g=goDag(funcAnnChartDAG),r=funcAnnChartDAG, max.nchar=40, node.shape="ellipse")      
      dev.off()
    }
    Sys.sleep(30)
    
    funcAnnCHartDAG<-1
    try(DAVIDGODag(funcAnnChart,pvalueCutoff = 0.01,"CC")->funcAnnChartDAG)
    if(!is.numeric(funcAnnChartDAG)) {      
      pdf(file=paste("enrichmentAnalysis/enrichmentDiagramCC_",prefix,g,".pdf",sep=""),paper="a4")
      plotGOTermGraph(g=goDag(funcAnnChartDAG),r=funcAnnChartDAG, max.nchar=40, node.shape="ellipse")      
      dev.off()
    }
    Sys.sleep(30)
  }
})
