library(biomaRt)
library(dynamicTreeCut)
library(ade4)

setOldClass("hclust")
setOldClass("dendrogram")
setOldClass("dist")

ProtClusterer <- setClass("ProtClusterer", 
                          slots= c(
                            proteins = "factor", 
                            biomartDomainRes = "data.frame",
                            iproDomsUniqTable = "data.frame",
                            davidEmail = "character",
                            clustProts = "hclust",
                            clustDomains = "hclust",
                            colorProtDend = "dendrogram",
                            proteinGroup = "numeric",
                            proteinDist = "dist",
                            domainDist = "dist"
                            ),
                          prototype = list(clustProts=structure(list(),class="hclust"),
                                           clustDomains=structure(list(),class="hclust"),
                                           colorProtDend=structure(list(),class="dendrogram"),
                                           proteinDist=structure(double(),class="dist"),
                                           domainDist=structure(double(),class="dist")
                                           )
                          )

setGeneric("retrieveInterProDomains", function(object,...) standardGeneric("retrieveInterProDomains"))
setMethod("retrieveInterProDomains","ProtClusterer", function(object) {
  interproBM<-useMart("prod-intermart_1",dataset="entry")
  # listFilters(interproBM)
  filters<-c("protein_ac")
  
  attributes<-c("protein_ac","protein_name","entry_id","entry_name","entry_type"
                ,"method_id","method_name","pos_from","pos_to","fragment","protein_length")
  getBM(attributes=attributes,filters=filters,values=object@proteins,mart=interproBM)->object@biomartDomainRes
  
  object@iproDomsUniqTable<-as.data.frame.matrix(table(unique(object@biomartDomainRes[,c(1,3)])))
  return(object)
})

setGeneric("calculateDistances",function(object,...) standardGeneric("calculateDistances"))
setMethod("calculateDistances","ProtClusterer", function(object,minDomainPresence=0) {
  toDrop <- names(which(colSums(object@iproDomsUniqTable)<minDomainPresence))
  iproDomsUniqTableAboveMin <- object@iproDomsUniqTable[!(colnames(object@iproDomsUniqTable) %in% toDrop)]
  # here we should use jaccard or other binary distance for presence/absence, not euclidean.
  # 1 if for jaccard
  dist.binary(df=as.matrix(iproDomsUniqTableAboveMin),method=1)->object@proteinDist
  dist.binary(df=t(as.matrix(iproDomsUniqTableAboveMin)),method=1)->object@domainDist
  # According to this paper http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6421371&tag=1 complete and ward are the best options for
  # binary data
  hclust(d=object@proteinDist,method="complete")->object@clustProts
  hclust(d=object@domainDist,method="complete")->object@clustDomains
  return(object)
})


reorderAsSlice <- function(sliceLike,tree) {
  unique_ids=which(!duplicated(sliceLike))
  # create a numeric vector where the position in the vector indicates new group
  # number i.e. if position 1 contains group 4 then old group 4 -> new group 1
  # NB x$order contains original indices of tree elements ordered by dendrogram
  xtable=order(match(unique_ids, tree$order))
  # now map old groups to new
  newgroups=structure(match(sliceLike,xtable),.Names=names(sliceLike))
  # finally return them in dendrogram order
  newgroups[tree$order]
}

setGeneric("getProteinsInClusters",function(object,...) standardGeneric("getProteinsInClusters"))
setMethod("getProteinsInClusters","ProtClusterer",function(object,clusterNumbers) {
  object@proteinGroup[object@proteinGroup %in% clusterNumbers]
})

setGeneric("redoWithClusters",function(object,...) standardGeneric("redoWithClusters"))
setMethod("redoWithClusters","ProtClusterer",function(object,clusterNumbers) {
  as.factor(names(getProteinsInClusters(object,clusterNumbers)))->proteins
  #toDrop <- names(which(colSums(object@iproDomsUniqTable[row.names(object@iproDomsUniqTable) %in% names(proteins)])>0))
  #iproDomsUniqTableForProts <- object@iproDomsUniqTable[!(colnames(object@iproDomsUniqTable) %in% toDrop)]
  ProtClusterer(proteins=proteins,davidEmail=object@davidEmail)->npc
  retrieveInterProDomains(npc)->npc
  calculateDistances(npc)->npc
  colourProteinClusters(npc,k=length(clusterNumbers),groupLabels=clusterNumbers)->npc
  npc
})

setGeneric("plotHeatMap",function(object,...) standardGeneric("plotHeatMap"))
setMethod("plotHeatMap","ProtClusterer",
          function(object,displayProt=FALSE,displayDom=FALSE) {
  protLab = if(displayProt) rownames(object@iproDomsUniqTable) else c("");
  domLab  = if(displayDom)  colnames(object@iproDomsUniqTable) else c("");
  heatmap.2(as.matrix(object@iproDomsUniqTable),
            Rowv=object@colorProtDend,
            col=colorpanel(2,"white","red"),breaks=3,trace="none",key=FALSE,labRow=protLab,labCol=domLab)
})

setGeneric("colourProteinClusters",function(object,...) standardGeneric("colourProteinClusters"))
setMethod("colourProteinClusters","ProtClusterer", function(object,k,h,groupLabels=TRUE) {
  colour_clusters(d=object@clustProts,k=k,h=h,groupLabels=groupLabels)->object@colorProtDend
  slice(object@clustProts,k=k,h=h)->object@proteinGroup
  if(length(groupLabels)>1) {
    for(i in 1:length(object@proteinGroup)) {
      object@proteinGroup[i] = groupLabels[object@proteinGroup[i]]
    }
  }
  return(object)
})

setGeneric("enrichmentAnalysisForClusters",function(object,...) standardGeneric("enrichmentAnalysisForClusters"))
setMethod("enrichmentAnalysisForClusters","ProtClusterer",function(object,clusterNumbers,prefix="") {
  david<-DAVIDWebService$new(email=object@davidEmail)
  setAnnotationCategories(david, c("GOTERM_BP_ALL","GOTERM_MF_ALL","GOTERM_CC_ALL"))
  slicedProteins = object@proteinGroup
  
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

setGeneric("autoColourProteinClusters",function(object,...) standardGeneric("autoColourProteinClusters"))
setMethod("autoColourProteinClusters","ProtClusterer", function(object,col=rainbow,groupLabels = NULL) {
  # cutreeHybrid(dendro=object@clustProts,distM=as.matrix(object@proteinDist),deepSplit=4,useMedoids=TRUE)->dynHybTree
  cutreeDynamic(dendro=object@clustProts,distM=as.matrix(object@proteinDist))->dynHybTree
  # we need to get from the dynHybTree object and the original table to something that looks like a slice:
  # > slice(clustRows,k=3)->sliceRes
  # > sliceRes
  # A2AB22 A2AB23 A0AVI4 A1L020 A2AE51 A2AAZ4 A2AAZ5 A2AAZ6 A2AE48 A2AE50 
  # 1      1      2      3      3      3      3      3      3      3 
  # Where each element is a leaf and each index the group to where that element was assigned
  # sliceLike<-as.numeric((dynHybTree$labels + 1))  
  sliceLike<-as.numeric((dynHybTree + 1))  
  names(sliceLike)<-object@clustProts$labels
  d<-object@clustProts
  sliceLike<-reorderAsSlice(sliceLike,d)
  object@proteinGroup<-sliceLike
  
  
  
  # from dendroextras colour_clusters
  if (!inherits(d, "dendrogram") && !inherits(d, "hclust")) 
    stop("Expects a dendrogram or hclust object")
  # g = slice(d, k = k, h = h)
  g <- sliceLike
  if (inherits(d, "hclust")) 
    d = as.dendrogram(d)
  k = max(g)
  if (is.function(col)) {
    col = col(k) 
  } else if (length(col) != k) { 
    stop("Must give same number of colours as clusters")
  }
  if (!is.null(groupLabels)) {
    if (length(groupLabels) == 1) {
      if (is.function(groupLabels)) 
        groupLabels = groupLabels(seq.int(length.out = k))
      else if (is.logical(groupLabels)) {
        if (groupLabels) 
          groupLabels = seq.int(length.out = k)
        else groupLabels = NULL
      }
    }
    if (!is.null(groupLabels) && length(groupLabels) != k) 
      stop("Must give same number of group labels as clusters")
  }
  
  ## These two functions are required for colourProteinClusters
  ## and belong to dendroextras colour_clusters
  addcol <- function(n, col) {
    attr(n, "edgePar") = c(attr(n, "edgePar"), list(col = col))
    n
  }
  descendTree <- function(sd) {
    groupsinsubtree = unique(g[labels(sd)])
    if (length(groupsinsubtree) > 1) {
      for (i in seq(sd)) sd[[i]] <- descendTree(sd[[i]])
    }
    else {
      sd = dendrapply(sd, addcol, col[groupsinsubtree])
      if (!is.null(groupLabels)) {
        attr(sd, "edgetext") = groupLabels[groupsinsubtree]
        attr(sd, "edgePar") = c(attr(sd, "edgePar"), 
                                list(p.border = col[groupsinsubtree]))
      }
    }
    sd
  }

  descendTree(d)->d
  
  object@colorProtDend<-d
  return(object)
})
