#' @include GeneralClusterer.R
PresenceAbsenceClusterer <- setClass("PresenceAbsenceClusterer", 
                             #slots= c(),
                             contains="GeneralClusterer"
)

setMethod("getColourModel",signature(object="PresenceAbsenceClusterer"),
          function(object) {
            return(colorpanel(2,"white","red"))
          })

#' @export
setMethod("calculateDistances",signature(object="PresenceAbsenceClusterer"), function(object,minFeaturePres=0) {
  toDrop <- names(which(colSums(object@uniqueFeaturesTable)<minFeaturePres))
  uniqTableFeatsAboveMin <- object@uniqueFeaturesTable[!(colnames(object@uniqueFeaturesTable) %in% toDrop)]
  toDropProts <- names(which(rowSums(uniqTableFeatsAboveMin)==0))
  object@uniqueFeaturesTable<-uniqTableFeatsAboveMin[!(rownames(uniqTableFeatsAboveMin) %in% toDropProts),]
  # here we should use jaccard or other binary distance for presence/absence, not euclidean.
  # 1 if for jaccard
  dist.binary(df=as.matrix(object@uniqueFeaturesTable),method=1)->object@proteinDist
  dist.binary(df=t(as.matrix(object@uniqueFeaturesTable)),method=1)->object@featureDist
  # According to this paper http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6421371&tag=1 complete and ward are 
  # the best options for binary data
  hclust(d=object@proteinDist,method="complete")->object@clustProts
  hclust(d=object@featureDist,method="complete")->object@clustFeatures
  return(object)
})