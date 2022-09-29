library(factoextra)
library(sm)

setwd("C:/Users/simon/Tesi_tome")

human <- read.csv("Datasets/h_zscored_top1000.csv")
murine <- read.csv("Datasets/m_zscored_top1000.csv")
samples <- colnames(human)[-1]
row.names(human) <- human$Geneid
row.names(murine) <- murine$Geneid
human <- human[samples]
murine <- murine[samples]


NH = as.numeric(160)
NM = as.numeric(60)

res.pca.h <- prcomp(t(human), scale = FALSE, rank. = NH)
res.pca.m <- prcomp(t(murine), scale = FALSE, rank. = NM)

# principal components 
pch <- res.pca.h$rotation
pcm <- res.pca.m$rotation

# rotate the dataset 
# t(pc) has dimension NUM.pc x 1000
# dataset has dimension 1000 x 624
# rotated dataset has dimension NUM.pc x 624
res1 <- t(pch) %*% as.matrix(human) 
res2 <- t(pcm) %*% as.matrix(murine) 

# build pairwise correlation for each dimension of the new space 
# row principal components of murine 
# col principal components of human
cor.matrix = matrix(nrow = NM , ncol = NH)
colnames(cor.matrix) <- colnames(res.pca.h$rotation)
rownames(cor.matrix) <- colnames(res.pca.m$rotation)

for(i in 1:NM){
  for(j in 1:NH){
    cor.matrix[i,j] <- cor(res2[i,], res1[j,], method = "pearson")
  }
}

# correlation matrix distribution 
d <- sm.density(x = as.vector(cor.matrix),
                xlab = "Correlation value",
                ylab = "Density")

# index in corr matrix of highly correlated principal components 
idx <- (which(data.frame(cor.matrix) > 0.5, arr.ind = TRUE)) # element 2 - 1



#2nd PCM
#1st PCH 


a <- pch[,1]
b <- pcm[,2]

res1. <- t(a) %*% as.matrix(human) 
res2. <- t(b) %*% as.matrix(murine) 

cor(res1.[1,],res2.[1,],method="pearson")

humang <- names(a[a>0.05])
murineg <- names(b[b>0.05])

# explained variance by the two components 

summary_pca <- function(x){
  eigs <- x$sdev^2
  res <- (t(data.frame(rbind(
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs),
    Cumulative = cumsum(eigs)/sum(eigs)))))
  row.names(res) <- seq(nrow(res))
  return(data.frame(res))
}

summary.h <- summary_pca(res.pca.h)
summary.m <- summary_pca(res.pca.m)

summary.h$Proportion[1]
summary.m$Proportion[2]

