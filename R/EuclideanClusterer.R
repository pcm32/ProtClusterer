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

setMethod("getColourModel",signature(object="EuclideanClusterer"),
          function(object) {
            return(colorpanel(99,"lightgrey",mid="purple",high="blue"))
          })
