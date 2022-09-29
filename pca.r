library(factoextra)
library(sm)

setwd("C:/Users/simon/Tesi_tome")

h_zscore <- read.csv("Datasets/h_zscored.csv")
m_zscore <- read.csv("Datasets/m_zscored.csv")
h_zscore_top1000 <- read.csv("Datasets/h_zscored_top1000.csv")
m_zscore_top1000 <- read.csv("Datasets/m_zscored_top1000.csv")

samples <- colnames(h_zscore)[-1]

row.names(h_zscore) <- h_zscore$Geneid
row.names(m_zscore) <- m_zscore$Geneid
row.names(h_zscore_top1000) <- h_zscore_top1000$Geneid
row.names(m_zscore_top1000) <- m_zscore_top1000$Geneid

h_zscore <- h_zscore[samples]
m_zscore <- m_zscore[samples]
h_zscore_top1000 <- h_zscore_top1000[samples]
m_zscore_top1000 <- m_zscore_top1000[samples]

# PCA utils 

summary_pca <- function(x){
  eigs <- x$sdev^2
  res <- (t(data.frame(rbind(
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs),
    Cumulative = cumsum(eigs)/sum(eigs)))))
  row.names(res) <- seq(nrow(res))
  return(data.frame(res))
}

# PCA calculation 

NUM_COMPONENTS_H = as.numeric(300)
NUM_COMPONENTS_M = as.numeric(300)
NUM_COMPONENTS_H_TOP = as.numeric(160)
NUM_COMPONENTS_M_TOP = as.numeric(60)

res.pca.h <- prcomp(t(h_zscore), scale = FALSE, rank. = NUM_COMPONENTS_H)
res.pca.m <- prcomp(t(m_zscore), scale = FALSE, rank. = NUM_COMPONENTS_M)
res.pca.h.1000 <- prcomp(t(h_zscore_top1000), scale = FALSE, rank. = NUM_COMPONENTS_H_TOP)
res.pca.m.1000 <- prcomp(t(m_zscore_top1000), scale = FALSE, rank. = NUM_COMPONENTS_M_TOP)

summary.h <- summary_pca(res.pca.h)
summary.m <- summary_pca(res.pca.m)
summary.h.1000 <- summary_pca(res.pca.h.1000)
summary.m.1000 <- summary_pca(res.pca.m.1000)

cat("Explained variance for h:",summary.h$Cumulative[NUM_COMPONENTS_H])
cat("Explained variance for m:",summary.m$Cumulative[NUM_COMPONENTS_M])
cat("Explained variance for h top 1000:",summary.h.1000$Cumulative[NUM_COMPONENTS_H_TOP])
cat("Explained variance for m top 1000:",summary.m.1000$Cumulative[NUM_COMPONENTS_M_TOP])

# principal.components <- res.pca$rotation
# transformed.data <- res.pca$x

png("RObjects/pca1.png",width=7.5,height=4,units="in",res=1200)
fviz_eig(res.pca.h, ncp = 40, geom = "line",title = "Human zscore PCA all genes")
dev.off()

png("RObjects/pca2.png",width=7.5,height=4,units="in",res=1200)
fviz_eig(res.pca.m, ncp = 40, geom = "line",title = "Murine zscore PCA all genes")
dev.off()

png("RObjects/pca3.png",width=7.5,height=4,units="in",res=1200)
fviz_eig(res.pca.h.1000, ncp = 40, geom = "line",title = "Human zscore PCA top genes")
dev.off()

png("RObjects/pca4.png",width=7.5,height=4,units="in",res=1200)
fviz_eig(res.pca.m.1000, ncp = 40, geom = "line",title = "Murine zscore PCA top genes")
dev.off()


# correlation between principal components 

pc.h <- res.pca.h$rotation
pc.m <- res.pca.m$rotation
pc.h.1000 <- res.pca.h.1000$rotation
pc.m.1000 <- res.pca.m.1000$rotation

cor1 <- matrix(nrow = NUM_COMPONENTS_H_TOP, ncol = NUM_COMPONENTS_M_TOP)
cor2 <- matrix(nrow = NUM_COMPONENTS_H_TOP, ncol = NUM_COMPONENTS_M)
cor3 <- matrix(nrow = NUM_COMPONENTS_H, ncol = NUM_COMPONENTS_M_TOP)
cor4 <- matrix(nrow = NUM_COMPONENTS_H, ncol = NUM_COMPONENTS_M)
 

# calclulate cor1 

for(i in 1:nrow(cor1)){
  for(j in 1:ncol(cor1)){
    cor1[i,j] <- cor(pc.h.1000[,i],
                     pc.m.1000[,j],
                     method = "pearson")
  }
}


png("RObjects/pca_cor.png",width=7.5,height=4,units="in",res=1200)
d <- sm.density(x = as.vector(cor1),
                xlab = "Correlation value",
                ylab = "Density")
title("Correlation between eigenvectors of top1000 h vs top1000 m")
plot(d)
dev.off()










