# single clustering robustness 

setwd("C:/Users/simon/Tesi_tome")
library(factoextra)
library(NbClust)
library(clustertend)
library(stylo)
library(amap)
library(cluster)
library(clusterSim)

m_zscored <- read.csv("datasets/m_zscored_top1000.csv")
h_zscored <- read.csv("datasets/h_zscored_top1000.csv")

row.names(m_zscored) <- m_zscored$Geneid
row.names(h_zscored) <- h_zscored$Geneid

samples <- colnames(h_zscored)[-1]

m_zscored <- m_zscored[,samples]
h_zscored <- h_zscored[,samples]

d_m <- dist.cosine(as.matrix(m_zscored))
d_h <- dist.cosine(as.matrix(h_zscored))



# eval robustness using silhouette and pearson correlation distance
silhouette_scores <- function(dataset,clustering,type){
  d <- Dist(dataset, method = type)
  return(s <- silhouette(clustering$clustering, d, diss = TRUE))
}

average_silhouette <- function(sil){
  sil.df <- as.data.frame(sil)
  means <- matrix(nrow = length(unique(sil.df$cluster)), ncol = 1)
  for(i in 1:nrow(means)){
    means[i,1] <- mean(sil.df[sil.df$cluster == i,"sil_width"])
  }
  return(mean(means[,1]))
}

db_score <- function(dataset,clustering,distance.matrix){
  return(index.DB(dataset, clustering$clustering,distance.matrix, centrotypes="medoids")$DB)
}

avg <- function(sil){
  return(mean(sil[,3]))
}

NUM <- as.numeric(50)

# evaluate total mean
res.m <- matrix(nrow = NUM,ncol = 1)
res.h <- matrix(nrow = NUM,ncol = 1)

# evaluate mean of means 

res.m.mom <- matrix(nrow = NUM,ncol = 1)
res.h.mom <- matrix(nrow = NUM,ncol = 1)


res.DB.m <- matrix(nrow = NUM,ncol = 1)
res.DB.h <- matrix(nrow = NUM,ncol = 1)


for(i in 2:NUM){
  print(i)
  cl.m <- pam(d_m,i,diss = TRUE)
  cl.h <- pam(d_h,i,diss = TRUE)
  s.m <- silhouette_scores(m_zscored,cl.m,type = "pearson")
  s.h <- silhouette_scores(h_zscored,cl.h,type = "pearson")
  res.m[i,1] <- avg(s.m)
  res.h[i,1] <- avg(s.h)
  
  res.m.mom[i,1] <- average_silhouette(s.m)
  res.h.mom[i,1] <- average_silhouette(s.h)
  
  res.DB.m[i,1] <- db_score(m_zscored,cl.m,d_m) 
  res.DB.h[i,1] <- db_score(h_zscored,cl.h,d_h)
  
}

# no calculating for single cluster 
res.m[1,] <- 0.
res.h[1,] <- 0.

res.m.mom[1,] <- 0.
res.h.mom[1,] <- 0.

res.DB.m[1,] <- 0.
res.DB.h[1,] <- 0.

png("RObjects/eval/pearson_sil_h_tot_mean.png")
plot(c(1:NUM),res.h[,1],main="Average silhouette for human clustering",
     xlab="K",
     ylab="average width",)


num <- which(res.h[,] == max(res.h[,1]))
value <- res.h[num,1]

points(num, value, col = "red", pch = 19)
text(num+1, value,as.character(num))

# add vertical line
abline(v = num, col = "blue", lty = "dashed")
?abline
dev.off()

png("RObjects/eval/pearson_sil_m_tot_mean.png")
plot(c(1:NUM),res.m[,1],main="Average silhouette for murine clustering",
     xlab="K",
     ylab="average width",)


num <- which(res.m[,1] == max(res.m[,1]))
value <- res.m[num,1]

points(num, value, col = "red", pch = 19)
text(num+1, value,as.character(num))

# add vertical line
abline(v = num, col = "blue", lty = "dashed")
?abline
dev.off()

###############################################################################

png("RObjects/eval/pearson_sil_h_mom.png")
plot(c(1:NUM),res.h.mom[,1],main="Average silhouette for human clustering - mom",
     xlab="K",
     ylab="average width",)


num <- which(res.h.mom[,] == max(res.h.mom[,1]))
value <- res.h.mom[num,1]

points(num, value, col = "red", pch = 19)
text(num+1, value,as.character(num))

# add vertical line
abline(v = num, col = "blue", lty = "dashed")
?abline
dev.off()

png("RObjects/eval/pearson_sil_m_mom.png")
plot(c(1:NUM),res.m.mom[,1],main="Average silhouette for murine clustering - mom",
     xlab="K",
     ylab="average width",)


num <- which(res.m.mom[,1] == max(res.m.mom[,1]))
value <- res.m.mom[num,1]

points(num, value, col = "red", pch = 19)
text(num+1, value,as.character(num))

# add vertical line
abline(v = num, col = "blue", lty = "dashed")
?abline
dev.off()

################################################################################

png("RObjects/eval/hDB.png")
plot(c(1:NUM),res.DB.h[,1],main="DB index for human clustering",
     xlab="K",
     ylab="average width",)


num <- which(res.DB.h[,] == min(res.DB.h[,1]))
value <- res.DB.h[num,1]

points(num, value, col = "red", pch = 19)
text(num+1, value,as.character(num))

# add vertical line
abline(v = num, col = "blue", lty = "dashed")
?abline
dev.off()

png("RObjects/eval/mDB.png")
plot(c(1:NUM),res.DB.m[,1],main="DB index for murine clustering",
     xlab="K",
     ylab="average width",)


num <- which(res.DB.m[,1] == min(res.DB.m[,1]))
value <- res.DB.m[num,1]

points(num, value, col = "red", pch = 19)
text(num+1, value,as.character(num))

# add vertical line
abline(v = num, col = "blue", lty = "dashed")
?abline
dev.off()







