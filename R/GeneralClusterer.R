setOldClass("hclust")
setOldClass("dendrogram")
setOldClass("dist")
setOldClass("integer")

#' An S4 class to represent the most generic type of cluster. This class is
#' never initialized as such. It is equivalent to an abstract class.
#' 
#' @slot proteins A list of UniProt identifiers as a factor for the proteins to
#'   be clustered.
#' @export
GeneralClusterer <- setClass("GeneralClusterer", 
                          slots= c(
                            proteins = "factor", 
                            retrievalResult = "data.frame", # retrievalResult
                            uniqueFeaturesTable = "data.frame", # unique results, presence absence
                            annotation = "data.frame",
                            clustProts = "hclust",
                            clustFeatures = "hclust", # cluster features
                            colorProtDend = "dendrogram",
                            proteinGroup = "numeric",
                            taxonID = "numeric",
                            proteinDist = "dist",
                            featureDist = "dist"   # feature distances.
                          ),
                          prototype = list(clustProts=structure(list(),class="hclust"),
                                           clustFeatures=structure(list(),class="hclust"),
                                           colorProtDend=structure(list(),class="dendrogram"),
                                           taxonID = 9606,
                                           proteinDist=structure(double(),class="dist"),
                                           featureDist=structure(double(),class="dist")
                          )
)

#' Generate clusters
#' 
#' \code{generateClusters} runs at once all the required method to produce data
#' to be ready to be plotted in a heatmap.
#' 
#' @param object An object which inherits from GeneralClusterer.
#' @param minFeaturePres The minimum number that a protein needs to show up in
#'   a feature to be considered for the distance calculations.
#' @param k The number of clusters (integer) to be enumerated.
#' @param h The height to cut the dendrogram at to produce the clusters. 
#'   \code{k} overrides this.
#' @param groupLabels whether to produce group/cluster labels.
#' 
#' @return a new object of the same type as the one provided with a matrix of 
#'   proteins/features filled, distances calculated and clusters colored.
#'  
#' @export
setGeneric("generateClusters",function(object,minFeaturePres=0,k,h,groupLabels=T,...) standardGeneric("generateClusters"))

#' Annotate proteins
#' 
#' \code{annotateProteins} retrieves annotation for the proteins through a web
#' service
#' 
#' @param object An object which inherits from GeneralClusterer.
#'   
#' @return a new object of the same type as the one provided, with annotations
setGeneric("annotateProteins",function(object) standardGeneric("annotateProteins"))

#' Retrieve features
#' 
#' \code{retrieveFeatures} connects to the appropiate data source and fetches 
#' whatever attributes are relevant for the clusterer (such as expression, 
#' presence of certain annotations, localizations, etc.).
#' 
#' @param object An object which inherits from GeneralClusterer.
#' @param ... parameters to be passed to the inheriting object for data 
#'   retrieval.
#'   
#' @return a new object of the same type as the one provided with a matrix of 
#'   proteins/features filled.
#' @export
setGeneric("retrieveFeatures",function(object,...) standardGeneric("retrieveFeatures"))

#' Calculate distances
#' 
#' \code{calculateDistances} uses the data retrieved in
#' \code{\link{retrieveFeatures}} to produces distances (according to the
#' object, that could be euclidian, binary, etc.) both within the proteins and
#' the features. It should be called after \code{retrieveFeatures}
#' 
#' @param object An object which inherits from GeneralClusterer.
#' @param minFeaturePres The minimum number that a protein needs to show up in
#'   a feature to be considered for the distance calculations.
#' @return a new object of the same type as the one provided with a matrix of 
#'   proteins/features filled and distances calculated.
#' @export
setGeneric("calculateDistances",function(object,minFeaturePres=0,...) standardGeneric("calculateDistances"))

#' Colour protein clusters
#' 
#' \code{colourProteinClusters} colours the protein clusters generated during 
#' the distance calculation step, according to the number of clusters \code{k} 
#' desired or the cut-off height \code{h} desired. Optionally group labels can 
#' be ommited by issuing \code{groupLabels=F} (defaults to \code{TRUE}). This 
#' method needs to be called after \code{\link{calculateDistances}}.
#' 
#' @param object An object which inherits from GeneralClusterer, for which data 
#'   has been retrieved and distances calculated
#' @param k The number of clusters (integer) to be enumerated.
#' @param h The height to cut the dendrogram at to produce the clusters. 
#'   \code{k} overrides this.
#' @param groupLabels whether to produce group/cluster labels.
#' @param ... other parameters for each clusterer.
#' @return a new object of the same type as the one provided with a matrix of 
#'   proteins/features filled, distances calculated and clusters colored.
#' @export
setGeneric("colourProteinClusters",function(object,k,h,groupLabels=T,...) standardGeneric("colourProteinClusters"))

#' Get proteins in clusters
#' 
#' \code{getProteinsInClusters} returns the set of proteins present in the 
#' clusters defined in the \code{clusterNumber} vector for the provided 
#' GeneralClusterer \code{object}. The cluster numbers are read by the user from
#' the figure produced by \code{colourProteinClusters}.
#' 
#' @param object An object which inherits from GeneralClusterer with a matrix of
#'   proteins/features filled, distances calculated and clusters colored.
#' @param clusterNumbers The clusters (defined by their numbers) for which
#'   proteins should be obtained.
#' @export
setGeneric("getProteinsInClusters",function(object,clusterNumbers,...) standardGeneric("getProteinsInClusters"))

#' Redo with Clusters
#' 
#' \code{redoWithClusters} runs the all the steps again for the proteins in the 
#' defined clusters.
#' 
#' @param object An object which inherits from GeneralClusterer with a matrix of
#'   proteins/features filled, distances calculated and clusters colored.
#' @param clusterNumbers The clusters for which the analysis should be redone.
#' @param ... additional parameters for each implementation.
#'   
#' @return a new object of the same type as the one provided with a matrix of 
#'   proteins/features filled, distances calculated and clusters colored, but
#'   based on the proteins available in the defined clusters of the object given
#'   as parameter.
#' @export
setGeneric("redoWithClusters",function(object,clusterNumbers,...) standardGeneric("redoWithClusters"))

#' Plot heatmap
#' 
#' \code{plotHeatMap} produces a heatmap of the proteins and features of the
#' object provided. Protein identifiers and feature identifiers can be shown
#' optionally.
#' 
#' @param object An object which inherits from GeneralClusterer with a matrix of
#'   proteins/features filled, distances calculated and clusters colored.
#' @param displayProt Sets whether the protein identifier should be shown in the
#'   heatmap, defaults to \code{False}.
#' @param displayFeat Sets whether the feature identifier should be shown in the
#' heatmap, defaults to \code{False}.
#' @param ... additional parameters for implementations.
#' @export
setGeneric("plotHeatMap",function(object,displayProt=F,displayFeat=F,useProtNames=TRUE,title="",...) standardGeneric("plotHeatMap"))
setGeneric("autoColourProteinClusters",function(object,...) standardGeneric("autoColourProteinClusters"))

setGeneric("getColourModel",function(object) standardGeneric("getColourModel"))

#' @export
setMethod("getProteinsInClusters",signature(object="GeneralClusterer",clusterNumbers="vector"),function(object,clusterNumbers) {
   object@proteinGroup[object@proteinGroup %in% clusterNumbers]
})

#' @export
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

#' @export
setMethod("redoWithClusters",signature(object="GeneralClusterer",clusterNumbers="vector"),function(object,clusterNumbers,type=class(object)[[1]]) {
  as.factor(names(getProteinsInClusters(object,clusterNumbers)))->proteins
  #toDrop <- names(which(colSums(object@uniqueFeaturesTable[row.names(object@uniqueFeaturesTable) %in% names(proteins)])>0))
  #iproDomsUniqTableForProts <- object@uniqueFeaturesTable[!(colnames(object@uniqueFeaturesTable) %in% toDrop)]
  new(type,proteins=proteins,taxonID=object@taxonID)->npc
  # ProtClusterer(proteins=proteins)->npc
  generateClusters(npc,k=length(clusterNumbers),groupLabels=clusterNumbers)->npc
  #retrieveFeatures(npc)->npc
  #calculateDistances(npc)->npc
  #colourProteinClusters(npc,k=length(clusterNumbers),groupLabels=clusterNumbers)->npc
  #npc
})

#' @export
setMethod("generateClusters",signature(object="GeneralClusterer"),function(object,minFeaturePres=0,k,h,groupLabels=T,...){
  annotateProteins(object)->object
  retrieveFeatures(object)->object
  calculateDistances(object,minFeaturePres=minFeaturePres)->object
  colourProteinClusters(object,groupLabels=groupLabels,k=k,h=h)->object
})

#' @export
setMethod("annotateProteins",signature(object="GeneralClusterer"),function(object) {
  if(object@taxonID!=9606) {
    taxId(UniProt.ws)<-object@taxonID
  }
  data.table(select(UniProt.ws,keys=object@proteins,columns=c("GENES"),keytype="UNIPROTKB"),key="UNIPROTKB")->prot2names
  data.table(select(UniProt.ws,keys=object@proteins,columns=c("ENSEMBL"),keytype="UNIPROTKB"),key="UNIPROTKB")->prot2ENSEMBL
  
  merge(prot2names,prot2ENSEMBL,all=TRUE)->object@annotation
  
  return(object)
})
 
#' @export
setMethod("plotHeatMap",signature(object="GeneralClusterer"),
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
            
            colourPanel<-getColourModel(object)
            
            featLab  = if(displayFeat)  colnames(object@uniqueFeaturesTable) else c("");
            heatmap.2(as.matrix(object@uniqueFeaturesTable),
                      Rowv=object@colorProtDend,
                      Colv=as.dendrogram(object@clustFeatures),
                      col=colourPanel,breaks=length(colourPanel)+1,
                      trace="none",key=FALSE,labRow=protLab,labCol=featLab,main=title,
                      margins=c(5,8),srtCol=45,cexCol=0.8,cexRow=0.8
            )->hm
            return(hm)
          })


