setwd("C:/Users/simon/Tesi_tome")
library(factoextra)
library(NbClust)
library(clustertend)
library(stylo)
library(amap)
library(cluster)
library(ComplexHeatmap)
library(testit)

m_zscored <- read.csv("datasets/m_zscored_top1000.csv")
h_zscored <- read.csv("datasets/h_zscored_top1000.csv")

row.names(m_zscored) <- m_zscored$Geneid
row.names(h_zscored) <- h_zscored$Geneid

samples <- colnames(h_zscored)[-1]

m_zscored <- m_zscored[,samples]
h_zscored <- h_zscored[,samples]

THRESHOLD = as.numeric(0.66)


build_metagenes <- function(dataset, k, clustering){
  
  dataset$Pam <- clustering$clustering
  assert("same order", {row.names(dataset) == names(clustering$clustering)})
  
  metagenes <- matrix(ncol = 624, nrow = k)
  colnames(metagenes) <- samples
  row.names(metagenes) <- seq(1,k)
  
  for (i in 1:k){
    metagenes[i,] <- colMeans(dataset[dataset$Pam==i,samples])
  }
  
  return(metagenes)
}

compute_cor.matrix <- function(x,y){
  cor.matrix <- matrix(nrow = nrow(x), ncol = nrow(y))
  # metagenes are on rows
  cor.matrix <- cor(t(x), t(y), method = "pearson")
  return(cor.matrix)
}


# evaluate clustering 
# x clustering for murine metagenes 
# y clustering for human metagenes
eval_clustering <- function(cor.matrix,x,y,type){
  
  switch (type,
    "gtratio" = res <- ((length(cor.matrix[cor.matrix >= THRESHOLD]))/ 
                  (nrow(x) * nrow(y))),
    "gt" =  res <- length(cor.matrix[cor.matrix >= THRESHOLD]),
    "max" = res <- (max(cor.matrix)),
    "top5" = res <- (as.numeric(quantile(cor.matrix, probs = c(0.95)))),
    "top1" = res <- (as.numeric(quantile(cor.matrix, probs = c(0.99)))),
    "gt2stdratio" = res <- length(cor.matrix[cor.matrix > (2*sd(cor.matrix))])/(nrow(x) * nrow(y))
  )
  
    return(res)
  
}




# try all clustering in the range LB to UB 

LB = as.numeric(10)
NUM = as.numeric(20)

evaluation.matrix.1 <- matrix(nrow = (NUM+1), ncol = (NUM+1))
evaluation.matrix.2 <- matrix(nrow = (NUM+1), ncol = (NUM+1))
evaluation.matrix.3 <- matrix(nrow = (NUM+1), ncol = (NUM+1))
evaluation.matrix.4 <- matrix(nrow = (NUM+1), ncol = (NUM+1))
evaluation.matrix.5 <- matrix(nrow = (NUM+1), ncol = (NUM+1))

d_m <- dist.cosine(as.matrix(m_zscored))
d_h <- dist.cosine(as.matrix(h_zscored))

for(i in 0:NUM){
  for(j in 0:NUM){
    
    cat("Evaluating ",i," ",j,"\n")
    
    cluster_m <- cluster::pam(d_m, k = i + LB, diss = TRUE)
    cluster_h <- cluster::pam(d_h, k = j + LB, diss = TRUE)
    metagenes_m <- build_metagenes(m_zscored,i + LB,cluster_m)
    metagenes_h <- build_metagenes(h_zscored,j + LB,cluster_h)
    cor.matrix <- compute_cor.matrix(metagenes_m,metagenes_h)
    res.1 <- eval_clustering(cor.matrix, metagenes_m, metagenes_h,"gt2stdratio")
    res.2 <- eval_clustering(cor.matrix, metagenes_m, metagenes_h,"gtratio")
    res.3 <- eval_clustering(cor.matrix, metagenes_m, metagenes_h,"max")
    res.4 <- eval_clustering(cor.matrix, metagenes_m, metagenes_h,"top5")
    res.5 <- eval_clustering(cor.matrix, metagenes_m, metagenes_h,"top1")
    evaluation.matrix.1[i+1,j+1] <- res.1
    evaluation.matrix.2[i+1,j+1] <- res.2
    evaluation.matrix.3[i+1,j+1] <- res.3
    evaluation.matrix.4[i+1,j+1] <- res.4
    evaluation.matrix.5[i+1,j+1] <- res.5
    
  }
}

row.names(evaluation.matrix.1) <- c(LB:(NUM+LB))
colnames(evaluation.matrix.1) <- c(LB:(NUM+LB))

row.names(evaluation.matrix.2) <- c(LB:(NUM+LB))
colnames(evaluation.matrix.2) <- c(LB:(NUM+LB))

row.names(evaluation.matrix.3) <- c(LB:(NUM+LB))
colnames(evaluation.matrix.3) <- c(LB:(NUM+LB))

row.names(evaluation.matrix.4) <- c(LB:(NUM+LB))
colnames(evaluation.matrix.4) <- c(LB:(NUM+LB))

row.names(evaluation.matrix.5) <- c(LB:(NUM+LB))
colnames(evaluation.matrix.5) <- c(LB:(NUM+LB))


png("RObjects/eval/kk-validation/gt2sdratio_unilatero.png")
Heatmap(evaluation.matrix.1,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Murine metagenes",
        column_title = "Human metagenes")
dev.off()

png("RObjects/eval/kk-validation/gtratio066_unilatero.png")
Heatmap(evaluation.matrix.2,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Murine metagenes",
        column_title = "Human metagenes"
        )
dev.off()

png("RObjects/eval/kk-validation/max.png")
Heatmap(evaluation.matrix.3,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Murine metagenes",
        column_title = "Human metagenes"
)
dev.off()

png("RObjects/eval/kk-validation/top5pcnt_unilatero.png")
Heatmap(evaluation.matrix.4,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Murine metagenes",
        column_title = "Human metagenes"
)
dev.off()

png("RObjects/eval/kk-validation/top1pcnt_unilatero.png")
Heatmap(evaluation.matrix.5,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Murine metagenes",
        column_title = "Human metagenes"
)
dev.off()

title("Correlation murine metagenes - top1000 human")

