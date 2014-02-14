library(dendroextras)
library(gplots)

setOldClass("hclust")
setOldClass("dendrogram")
setOldClass("dist")
setOldClass("integer")

GeneralClusterer <- setClass("GeneralClusterer", 
                          slots= c(
                            proteins = "factor", 
                            retrievalResult = "data.frame", # retrievalResult
                            uniqueFeaturesTable = "data.frame", # unique results, presence absence
                            clustProts = "hclust",
                            clustFeatures = "hclust", # cluster features
                            colorProtDend = "dendrogram",
                            proteinGroup = "numeric",
                            proteinDist = "dist",
                            featureDist = "dist"   # feature distances.
                          ),
                          prototype = list(clustProts=structure(list(),class="hclust"),
                                           clustFeatures=structure(list(),class="hclust"),
                                           colorProtDend=structure(list(),class="dendrogram"),
                                           proteinDist=structure(double(),class="dist"),
                                           featureDist=structure(double(),class="dist")
                          )
)

setGeneric("retrieveFeatures",function(object,...) standardGeneric("retrieveFeatures"))
setGeneric("calculateDistances",function(object,minFeaturePres=0,...) standardGeneric("calculateDistances"))
setGeneric("getProteinsInClusters",function(object,clusterNumbers,...) standardGeneric("getProteinsInClusters"))
setGeneric("redoWithClusters",function(object,clusterNumbers,...) standardGeneric("redoWithClusters"))
setGeneric("plotHeatMap",function(object,displayProt=F,displayFeat=F,...) standardGeneric("plotHeatMap"))
setGeneric("colourProteinClusters",function(object,k,h,groupLabels=T,...) standardGeneric("colourProteinClusters"))
setGeneric("autoColourProteinClusters",function(object,...) standardGeneric("autoColourProteinClusters"))

setMethod("getProteinsInClusters",signature(object="GeneralClusterer",clusterNumbers="vector"),function(object,clusterNumbers) {
   object@proteinGroup[object@proteinGroup %in% clusterNumbers]
})

setMethod("colourProteinClusters",signature(object="GeneralClusterer"),function(object,k,h,groupLabels=TRUE) {
   colour_clusters(d=object@clustProts,k=k,h=h,groupLabels=groupLabels)->object@colorProtDend
   slice(object@clustProts,k=k,h=h)->object@proteinGroup
   if(length(groupLabels)>1) {
     for(i in 1:length(object@proteinGroup)) {
       object@proteinGroup[i] = groupLabels[object@proteinGroup[i]]
     }
   }
   return(object)
})

setMethod("redoWithClusters",signature(object="GeneralClusterer",clusterNumbers="vector"),function(object,clusterNumbers,type=class(object)[[1]]) {
  as.factor(names(getProteinsInClusters(object,clusterNumbers)))->proteins
  #toDrop <- names(which(colSums(object@uniqueFeaturesTable[row.names(object@uniqueFeaturesTable) %in% names(proteins)])>0))
  #iproDomsUniqTableForProts <- object@uniqueFeaturesTable[!(colnames(object@uniqueFeaturesTable) %in% toDrop)]
  new(type,proteins=proteins)->npc
  # ProtClusterer(proteins=proteins)->npc
  retrieveFeatures(npc)->npc
  calculateDistances(npc)->npc
  colourProteinClusters(npc,k=length(clusterNumbers),groupLabels=clusterNumbers)->npc
  npc
})

setMethod("plotHeatMap",signature(object="GeneralClusterer"),
          function(object,displayProt=FALSE,displayFeat=FALSE) {
  protLab = if(displayProt) rownames(object@uniqueFeaturesTable) else c("");
  featLab  = if(displayFeat)  colnames(object@uniqueFeaturesTable) else c("");
  heatmap.2(as.matrix(object@uniqueFeaturesTable),
            Rowv=object@colorProtDend,
            col=colorpanel(2,"white","red"),breaks=3,trace="none",key=FALSE,labRow=protLab,labCol=featLab)->hm
  return(hm)
})
