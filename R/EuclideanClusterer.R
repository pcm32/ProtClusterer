library(gplots)

EuclideanClusterer <- setClass("EuclideanClusterer", 
                                     #slots= c(),
                                     contains="GeneralClusterer"
)

setMethod("plotHeatMap",signature(object="EuclideanClusterer"),
          function(object,displayProt=FALSE,displayFeat=FALSE) {
            protLab = if(displayProt) rownames(object@uniqueFeaturesTable) else c("");
            featLab  = if(displayFeat)  colnames(object@uniqueFeaturesTable) else c("");
            heatmap.2(as.matrix(object@uniqueFeaturesTable),
                      Rowv=object@colorProtDend,
                      Colv=as.dendrogram(object@clustFeatures),
                      col=colorpanel(99,"lightgrey",mid="purple",high="blue"),breaks=100,
                      trace="none",key=FALSE,labRow=protLab,labCol=featLab)->hm
            return(hm)
          })