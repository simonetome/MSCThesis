library(WGCNA)
setwd("C:/Users/simon/Tesi_tome")
library(gprofiler2)
library(stylo)



gene.annotations <- read.csv(file = "datasets/orthologs/Gene_annotations.csv",sep=";")

perform.WGCNA <- function(dataset,title,min_size,deepsplit,cortype,type.adj,RsquaredCut){
  options(stringsAsFactors = FALSE)
  powers <- c(1:20)
  
  #==============================================================================#
  #============================ pick soft threshold =============================#
  #==============================================================================#
  
  
  sft <- pickSoftThreshold(dataset, 
                           dataIsExpr = TRUE, 
                           powerVector = powers, 
                           verbose = 5,
                           RsquaredCut = RsquaredCut,
                           )
  
  pdf(file = paste(title,"soft_threshold_picking.pdf"), width = 12, height = 9)
  par(mfrow = c(1,2))
  cex1 <-  0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red") 
  
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity")) 
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  # power decided with the scale free topology criterion
  power <- sft$powerEstimate
  
  # mask base cor 
  cor <- WGCNA::cor
  
  
  #==============================================================================#
  #============================ WGCNA steps == ==================================#
  #==============================================================================#
  
  # adjacency by soft thresholding, default is PEARSON 
  # unsigned type takes module of correlation 
  # signed type takes (cor+0.5)/2 
  # signed hybrid takes cor if > 0, 0 otherwise 
  
  adj <- adjacency(datExpr = dataset, 
                   type = type.adj,
                   power = power,
                   corFnc = cor,
                   corOptions = list(method = cortype,use="all.obs")
  )
  
  TOM <- TOMsimilarity(adj)
  dissTOM <- 1 - TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  # Module identification using dynamic tree cut
  # We like large modules, so we set the minimum module size relatively high 
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = deepsplit, 
                               pamRespectsDendro = FALSE,minClusterSize = min_size);
  dynamicColors <- labels2colors(dynamicMods)
  
  pdf(file = paste(title,"module_tree.pdf"), width = 8, height = 6);
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                      hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
  dev.off()
  
  modules <- data.frame(row.names = colnames(dataset))
  modules$Geneid <- colnames(dataset)
  modules$Label <- dynamicMods
  modules$Colors <- dynamicColors
  modules$clustering <- dynamicMods
  
  print(unique(modules$Colors))
  
  #==============================================================================#
  #============================ Eigengenes====== ================================#
  #==============================================================================#
  
  # Calculate eigengenes
  eigengenes.data <- moduleEigengenes(dataset, colors = dynamicColors)
  eigengenes <- eigengenes.data$eigengenes
  # Calculate dissimilarity of module eigengenes
  eigengenes.dissimilarity <- 1-cor(eigengenes)
  METree = hclust(as.dist(eigengenes.dissimilarity), method = "average")
  
  pdf(file = paste(title,"eigengenes_tree.pdf"), width = 8, height = 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  
  dev.off()
  
  l <- list(geneTree,sft,adj,TOM,modules,eigengenes.data)
  names(l) <- c("geneTree","sft","adj","TOM","modules","eigengenes")
  
  return(l)
  
}


perform.KMODULE <- function(dataset,wcgna.analysis,title){
  
  changing <- TRUE
  modules <- wcgna.analysis$modules
  new.modules <- modules
  greyLabel <- 0
  
  # vector that contains the mean connectivity to each module for a certain gene
  mean.connectivity <- vector(length = length(unique(modules$Label)))
  mean.connectivity[1] <- 0
  
  new.adj <- wcgna.analysis$adj
  diag(new.adj) <- 0
  
  total.changes <- 0
  nIter <- 0
  
  greyGenes <- as.list(new.modules[new.modules$Label == 0, "Geneid"])
  
  while(changing & nIter < 30){
    # changes in the current iteration
    changes <- 0
    for(i in modules$Geneid){
      if(!(i %in% greyGenes)){
        # calculate mean connectivity for gene i to each module m 
        for(m in c(1:max(unique(new.modules$Label)))){
          n.m <- nrow(new.modules[new.modules$Label == m,])
          s.m <- new.modules[new.modules$Label == m,"Geneid"]
          mean.connectivity[m+1] <- 1/n.m * (sum(new.adj[i,s.m])) 
        }
        # label of the module 
        mod.to.assign <- which.max(mean.connectivity) - 1
        if(new.modules[i,"Label"] != mod.to.assign){
          print(paste(i,new.modules[i,"Label"],"->",mod.to.assign))
          changes <- changes + 1
          total.changes <- total.changes +1
          new.modules[i,"Label"] <- mod.to.assign
        }
      }
    }
    if(changes == 0){
      changing <- FALSE
    }
    nIter <- nIter + 1
  }
  
  new.modules$Colors <- labels2colors(new.modules$Label)
  new.modules$clustering <- new.modules$Label
  
  
  
  
  pdf(file = paste(title,"kmodule_result.pdf"), width = 8, height = 6)
  plotDendroAndColors(wcgna.analysis$geneTree, cbind(modules$Colors, new.modules$Colors), 
                      c("Dynamic Tree Cut", "Kmodule"), dendroLabels = FALSE, 
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
  dev.off()
  
  # Calculate eigengenes
  eigengenes.data <- moduleEigengenes(dataset, colors = new.modules$Colors)
  eigengenes <- eigengenes.data$eigengenes
  # Calculate dissimilarity of module eigengenes
  eigengenes.dissimilarity <- 1-cor(eigengenes)
  METree = hclust(as.dist(eigengenes.dissimilarity), method = "average")
  
  pdf(file = paste(title,"eigengenes_tree_kmodule.pdf"), width = 8, height = 6)
  plot(METree, main = "Clustering of Kmodule eigengenes",
       xlab = "", sub = "")
  dev.off()
  
  l <- list(new.modules,eigengenes.data)
  names(l) <- c("modules","eigengenes")
  
  return(l)
  
  
}




convert.to.ens <- function(geneSet,species){
  result <- list()
  if(species == "mmusculus"){
    for(g in geneSet){
      x <- gene.annotations[gene.annotations$Geneid == paste("M_",g,sep=""),"Ensembl"]
      result <- append(result,x)
    }  
  }
  else{
    for(g in geneSet){
      x <- gene.annotations[gene.annotations$Geneid == paste("H_",g,sep=""),"Ensembl"]
      result <- append(result,x)
    }
  }
  return(result)
}

perform.FUNENRICHMENT <- function(modules,species){
  
  fun.analysis <- vector(mode = "list",length = length(unique(modules$Colors)))
  names(fun.analysis) <- unique(modules$Colors)
  
  for(c in unique(modules$Colors)){
    fun.analysis[[c]] <- gost(query = modules[modules$Colors == c, "Geneid"],organism = species) 
  }
  
  return(fun.analysis)
}












