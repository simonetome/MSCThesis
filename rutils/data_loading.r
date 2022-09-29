library(matrixStats)

zs <- function(x){
  return((x-mean(x))/(sd(x)+1e-17))
}

l2r <- function(x){
  x <- log2(x + 1.)
  return(x - mean(x))
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

zscore_row <- function(assembly){
  assembly <- t(apply(assembly,1,zs))
  return(assembly)
}

zscore_col <- function(assembly){
  assembly <- (apply(assembly,2,zs))
  return(assembly)
}

l2r_row <- function(assembly){
  assembly <- t(apply(assembly,1,l2r))
  return(assembly)
}

# rdsFile is the path of the rds and can be either 
# "C:\Users\simon\MSCThesis\datasets\robjects\full_assembly_gt400.rds"
# "C:\Users\simon\MSCThesis\datasets\robjects\full_assembly.rds"
# the dataset is already filtered and in tpm 
# the dataset is splitted two:
# [[1]] -> human
# [[2]] -> murine
load_dataset <- function(gt,
                         zscore,
                         topNum){
  if(gt == TRUE){
    rdsFile = "C:\\Users\\simon\\MSCThesis\\datasets\\robjects\\full_assembly_gt400.rds"
  }else{
    rdsFile = "C:\\Users\\simon\\MSCThesis\\datasets\\robjects\\full_assembly.rds"
  }
  assembly = readRDS(rdsFile)
  human = assembly[[1]]
  murine = assembly[[2]]
  
  if(!is.null(topNum)){
    human = top_variance(human,topNum)
    murine = top_variance(murine,topNum)
  }
  
  if(zscore == TRUE){
    human = zscore_row(human)
    murine = zscore_row(murine)
  }
  
  assembly[[1]] = human
  assembly[[2]] = murine
  
  return(assembly)
  
}





















