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

#' @export
setMethod("calculateDistances",signature(object="PresenceAbsenceClusterer"), function(object,minFeaturePres=0) {
  toDrop <- names(which(colSums(object@uniqueFeaturesTable)<minFeaturePres))
  uniqTableFeatsAboveMin <- object@uniqueFeaturesTable[!(colnames(object@uniqueFeaturesTable) %in% toDrop)]
  # here we should use jaccard or other binary distance for presence/absence, not euclidean.
  # 1 if for jaccard
  dist.binary(df=as.matrix(uniqTableFeatsAboveMin),method=1)->object@proteinDist
  dist.binary(df=t(as.matrix(uniqTableFeatsAboveMin)),method=1)->object@featureDist
  # According to this paper http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6421371&tag=1 complete and ward are 
  # the best options for binary data
  hclust(d=object@proteinDist,method="complete")->object@clustProts
  hclust(d=object@featureDist,method="complete")->object@clustFeatures
  return(object)
})