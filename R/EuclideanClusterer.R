#' @include GeneralClusterer.R
EuclideanClusterer <- setClass("EuclideanClusterer", 
                                     #slots= c(),
                                     contains="GeneralClusterer"
)

#' @export
setMethod("calculateDistances",signature(object="EuclideanClusterer"), function(object,minFeaturePres=0) {
  #toDrop <- names(which(colSums(object@uniqueFeaturesTable)<minFeaturePres))
  #uniqTableFeatsAboveMin <- object@uniqueFeaturesTable[!(colnames(object@uniqueFeaturesTable) %in% toDrop)]
  # here we should use jaccard or other binary distance for presence/absence, not euclidean.
  # 1 if for jaccard
  dist(x=as.matrix(object@uniqueFeaturesTable),method="euclidean")->object@proteinDist
  dist(x=t(as.matrix(object@uniqueFeaturesTable)),method="euclidean")->object@featureDist

  hclust(d=object@proteinDist,method="complete")->object@clustProts
  hclust(d=object@featureDist,method="complete")->object@clustFeatures
  return(object)
})

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