#' @include GeneralClusterer.R
PresenceAbsenceClusterer <- setClass("PresenceAbsenceClusterer", 
                             #slots= c(),
                             contains="GeneralClusterer"
)

#' \code{plotHeatMap} plots the heatmap for the provided ProtClusterer object. 
#' In the presence/absence implementation, a particular color mode is used.
#' 
#' @param object The ProtClusterer object, with distances calculated, and 
#'   clusters generated.
#' @param displayProt Boolean to display the protein identifier. Defaults
#'   \code{FALSE}.
#' @param displayFeat Boolean to display the feature (domain, expression tissue,
#'   etc) identifier. Defaults to \code{FALSE}.
#'   
#' @export
setMethod("plotHeatMap",signature(object="PresenceAbsenceClusterer"),
          function(object,displayProt=FALSE,displayFeat=FALSE) {
            protLab = if(displayProt) rownames(object@uniqueFeaturesTable) else c("");
            featLab  = if(displayFeat)  colnames(object@uniqueFeaturesTable) else c("");
            heatmap.2(as.matrix(object@uniqueFeaturesTable),
                      Rowv=object@colorProtDend,
                      Colv=as.dendrogram(object@clustFeatures),
                      col=colorpanel(2,"white","red"),breaks=3,trace="none",key=FALSE,labRow=protLab,labCol=featLab)->hm
            return(hm)
          })