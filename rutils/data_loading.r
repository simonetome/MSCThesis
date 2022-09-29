library(matrixStats)

logtpm <- function (x){
  return(log2(x + 1.))
}

zs <- function(x){
  return((x-mean(x))/(sd(x)+1e-17))
}

l2r <- function(x){
  x <- log2(x + 1.)
  return(x - mean(x))
}

split_dataset <- function (assembly){
  genes <- row.names(assembly)
  human.genes <- genes[startsWith(genes,"H")]
  murine.genes <- genes[startsWith(genes,"M")]
  assembly_new <- list()
  assembly_new[[1]] <- assembly[human.genes,]
  assembly_new[[2]] <- assembly[murine.genes,]
  return(assembly_new)
}


top_variance <- function (assembly, num){
  assembly %>% as.data.frame() -> assembly
  vars <- rowVars(as.matrix(assembly))
  variances <- as.data.frame(matrix(nrow = length(vars), ncol = 2 ))
  colnames(variances) = c("Geneid","Variance")
  
  variances$Geneid = row.names(assembly)
  variances$Variance = vars 
  
  ordered <-variances[order(variances$Variance),]
  interesting <- ordered[(nrow(ordered)-num+1):nrow(ordered),"Geneid"]
  return(assembly[interesting,])
}

load_raw <- function (path, top.num = NULL, length.filter = NULL){
  assembly <- read.csv(path)
  row.names(assembly) <- assembly$Geneid
  if(!is.null(length.filter)){
    assembly <- assembly[assembly$Length > length.filter,]
  }
  
  
  samples <- colnames(assembly)[3:length(colnames(assembly))]
  assembly <- assembly[,samples]
  
  assembly <- split_dataset(assembly)
  
  if(!is.null(top.num)){
    assembly[[1]] <- top_variance(assembly[[1]], top.num)
    assembly[[2]] <- top_variance(assembly[[2]], top.num) 
  }
  
  return(assembly)
  
 
}

load_log <- function (path, top.num = NULL, length.filter = NULL){
  assembly <- read.csv(path)
  row.names(assembly) <- assembly$Geneid
  
  if(!is.null(length.filter)){
    assembly <- assembly[assembly$Length > length.filter,]
  }
  
  samples <- colnames(assembly)[3:length(colnames(assembly))]
  assembly <- assembly[,samples]
  assembly <- (data.frame(lapply(assembly,logtpm)))
  
  assembly <- split_dataset(assembly)
  
  if(!is.null(top.num)){
    assembly[[1]] <- top_variance(assembly[[1]], top.num)
    assembly[[2]] <- top_variance(assembly[[2]], top.num) 
  }
  
  return(assembly)
  
}

load_zscored <- function (path, top.num = NULL, length.filter = NULL){
  assembly <- read.csv(path)
  row.names(assembly) <- assembly$Geneid
  
  if(!is.null(length.filter)){
    assembly <- assembly[assembly$Length > length.filter,]
  }
  
  samples <- colnames(assembly)[3:length(colnames(assembly))]
  assembly <- assembly[,samples]
  
  assembly <- split_dataset(assembly)
  
  if(!is.null(top.num)){
    assembly[[1]] <- top_variance(assembly[[1]], top.num)
    assembly[[2]] <- top_variance(assembly[[2]], top.num) 
  }
  
  assembly[[1]] <- t(apply(assembly[[1]],1,zs))
  assembly[[2]] <- t(apply(assembly[[2]],1,zs))
  
  return(assembly)
  
}

load_log2ratio <- function (path, top.num = NULL){
  assembly <- read.csv(path)
  row.names(assembly) <- assembly$Geneid
  samples <- colnames(assembly)[3:length(colnames(assembly))]
  assembly <- assembly[,samples]
  
  assembly <- split_dataset(assembly)
  
  if(!is.null(top.num)){
    assembly[[1]] <- top_variance(assembly[[1]], top.num)
    assembly[[2]] <- top_variance(assembly[[2]], top.num) 
  }
  
  assembly[[1]] <- t(apply(assembly[[1]],1,l2r))
  assembly[[2]] <- t(apply(assembly[[2]],1,l2r))
  
  return(assembly)
  
}

zscore <- function(assembly){
  assembly <- t(apply(assembly,1,zs))
  return(assembly)
}

zscore_col <- function(assembly){
  assembly <- (apply(assembly,2,zs))
  return(assembly)
}

l2r <- function(assembly){
  assembly <- t(apply(assembly,1,l2r))
  return(assembly)
}






