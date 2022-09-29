setwd("C:/Users/simon/Tesi_tome")
library(factoextra)
library(NbClust)
library(clustertend)
library(stylo)
library(amap)
library(cluster)
library(clusterSim)
library(fpc)


m_zscored <- read.csv("datasets/m_zscored_top1000.csv")
h_zscored <- read.csv("datasets/h_zscored_top1000.csv")

row.names(m_zscored) <- m_zscored$Geneid
row.names(h_zscored) <- h_zscored$Geneid

samples <- colnames(h_zscored)[-1]

m_zscored <- m_zscored[,samples]
h_zscored <- h_zscored[,samples]



#d_m <- dist(as.matrix(m_zscored), method = "euclidean")
#d_h <- dist(as.matrix(h_zscored), method = "euclidean")

d_m <- dist.cosine(as.matrix(m_zscored))
d_h <- dist.cosine(as.matrix(h_zscored))

# interesting indexes 
# ch - Calinksi-Harabasz index -> better higher
# dunn - dunn index -> better higher 
# wb.ratio - avg.within/avg.between where -> better minimum
# avg.silwidth - average sil width -> better higher 
# clus.avg.silwidths . vector of cl avg sil width -> to avg better higher

# study cluster between 1 and 30 

res.m <- matrix(nrow = 30, ncol = 6)
res.h <- matrix(nrow = 30, ncol = 6)

colnames(res.m) <- c("ch","dunn","wb.ratio","avg.sil","avg.sil.mom","sindex")
colnames(res.h) <- c("ch","dunn","wb.ratio","avg.sil","avg.sil.mom","sindex")

res.m[1,1] <- 0. 
res.m[1,2] <- 0.
res.m[1,3] <- Inf
res.m[1,4] <- 0.
res.m[1,5] <- 0.
res.m[1,6] <- 0.
res.h[1,1] <- 0. 
res.h[1,2] <- 0.
res.h[1,3] <- Inf
res.h[1,4] <- 0.
res.h[1,5] <- 0.
res.h[1,6] <- 0.
  
for(i in 2:30){
  print(i)
  
  pam.m <- pam(d_m,i,diss = TRUE)
  pam.h <- pam(d_h,i,diss = TRUE)
  
  c.coeff <- (1000-i)/(i-1)
  
  stats.m <- cluster.stats(d = d_m, pam.m$clustering)
  stats.h <- cluster.stats(d = d_h, pam.h$clustering)
  
  ch.m <- c.coeff * (1/stats.m$wb.ratio)
  ch.h <- c.coeff * (1/stats.h$wb.ratio)
  
  res.m[i,1] <- ch.m
  res.m[i,2] <- stats.m$dunn
  res.m[i,3] <- stats.m$wb.ratio
  res.m[i,4] <- stats.m$avg.silwidth
  res.m[i,5] <- mean(stats.m$clus.avg.silwidths)
  res.m[i,6] <- stats.m$sindex
  res.h[i,1] <- ch.h
  res.h[i,2] <- stats.h$dunn
  res.h[i,3] <- stats.h$wb.ratio
  res.h[i,4] <- stats.h$avg.silwidth
  res.h[i,5] <- mean(stats.h$clus.avg.silwidths)
  res.h[i,6] <- stats.h$sindex

}

library(plotrix)


plot_metric <- function(path,metric,type,species){
  
  png(path)
  if(species == "murine"){
    res <- res.m
  }
  else{
    res <- res.h
  }

  plot(c(2:30),res[2:30,metric],main=paste(species,metric,"- better if ",type),
       xlab="K",
       ylab="Value",
       type = "b")
  
  if(type == "max"){
    num <- which(res[,metric] == max(res[,metric]))  
  }
  else{
    num <- which(res[,metric] == min(res[,metric]))
  }
  
  value <- res[num,metric]
  
  points(num, value, col = "red", pch = 19)
  text(num+1, value,as.character(num))
  
  # add vertical line
  abline(v = num, col = "blue", lty = "dashed")
  dev.off()
  
  
}


plot_metric("RObjects/eval/sil_murine.png","avg.sil","max","murine")
plot_metric("RObjects/eval/silmom_murine.png","avg.sil.mom","max","murine")
plot_metric("RObjects/eval/wbratio_murine.png","wb.ratio","min","murine")
plot_metric("RObjects/eval/ch_murine.png","ch","max","murine")
plot_metric("RObjects/eval/dunn_murine.png","dunn","max","murine")
plot_metric("RObjects/eval/sindex_murine.png","sindex","max","murine")

plot_metric("RObjects/eval/sil_human.png","avg.sil","max","human")
plot_metric("RObjects/eval/silmom_human.png","avg.sil.mom","max","human")
plot_metric("RObjects/eval/wbratio_human.png","wb.ratio","min","human")
plot_metric("RObjects/eval/ch_human.png","ch","max","human")
plot_metric("RObjects/eval/dunn_human.png","dunn","max","human")
plot_metric("RObjects/eval/sindex_human.png","sindex","max","human")







