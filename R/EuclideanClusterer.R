#' @include GeneralClusterer.R
EuclideanClusterer <- setClass("EuclideanClusterer", 
                                     #slots= c(),
                                     contains="GeneralClusterer"
)


#' In the euclidean implementation, a particular color mode is used.
#' 
#' @param object The ProtClusterer object, with distances calculated, and 
#'   clusters generated.
#' @param displayProt Boolean to display the protein identifier. Defaults
#'   \code{FALSE}.
#' @param displayFeat Boolean to display the feature (domain, expression tissue,
#'   etc) identifier. Defaults to \code{FALSE}.
#'   
#' @export
setMethod("plotHeatMap",signature(object="EuclideanClusterer"),
          function(object,displayProt=FALSE,displayFeat=FALSE,useProtNames=TRUE,title="") {
            protLab = c("")
            if(displayProt) { 
              if(useProtNames) {
                # first retrieve only the first name (make a new data table)
                object@annotation[,list(GENE=unlist(strsplit(GENES,split=' '))[1]),by=UNIPROTKB]->id2names
                # This part should be used by both the Euclidean and Presence absence
                protLab <- id2names[match(rownames(object@uniqueFeaturesTable),id2names$UNIPROTKB),GENE,with=T]
              } else {
                protLab = rownames(object@uniqueFeaturesTable)
              }
            }
              
              
            featLab  = if(displayFeat)  colnames(object@uniqueFeaturesTable) else c("");
            heatmap.2(as.matrix(object@uniqueFeaturesTable),
                      Rowv=object@colorProtDend,
                      Colv=as.dendrogram(object@clustFeatures),
                      col=colorpanel(99,"lightgrey",mid="purple",high="blue"),breaks=100,
                      trace="none",key=FALSE,labRow=protLab,labCol=featLab,main=title
                      )->hm
            return(hm)
          })