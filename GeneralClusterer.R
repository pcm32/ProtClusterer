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