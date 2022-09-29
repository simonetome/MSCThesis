setwd("C:/Users/simon/Tesi_tome")
source("wgcna_wrapper.r")
source("rutils/data_loading.r")
library("stringr")
library("stylo")
library(fpc)
library(cluster)
library(amap)

zscored <- load_zscored("datasets/filtered_assembly.csv", top.num = 1000)
human <- t(as.data.frame(zscored[1]))
murine <- t(as.data.frame(zscored[2]))

colnames(human) <- unlist(lapply(colnames(human),function(x){return(substr(x,3,nchar(x)))}))
colnames(murine) <- unlist(lapply(colnames(murine),function(x){return(substr(x,3,nchar(x)))}))


#==============================================================================#
#====================== Comparison to choose WGCNA parameters =================#
#==============================================================================#

library(ComplexHeatmap)

parameters.m = as.data.frame(matrix(ncol = 3, nrow = 4))
row.names(parameters.m) = seq(1,4)
colnames(parameters.m) = seq(20,30,5)

parameters.h = as.data.frame(matrix(ncol = 3, nrow = 4))
row.names(parameters.h) = seq(1,4)
colnames(parameters.h) = seq(20,30,5)


for(g in seq(20,30,5)){
  for(d in seq(1,4)){
    human.analysis <- perform.WGCNA(human,"epithelial",min_size = g, deepsplit = d,type.adj = "signed hybrid")
    murine.analysis <- perform.WGCNA(murine,"stromal",min_size = g, deepsplit = d,type.adj = "signed hybrid")
    
    parameters.h[[as.character(g)]][d] <- length(human.analysis$modules[human.analysis$modules$Colors == "grey",1])/1000
    parameters.m[[as.character(g)]][d] <- length(murine.analysis$modules[murine.analysis$modules$Colors == "grey",1])/1000
    
  }
}


?Heatmap
pdf("results/grey_ratio.pdf", width = 12, height = 6)

h1 = Heatmap(as.matrix(parameters.m),
        column_title = "Minimum module size stromal",
        row_title = "Deep split parameter",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        name = "Stromal grey module ratio")

h2 = Heatmap(as.matrix(parameters.h),
        column_title = "Minimum module size epithelial",
        row_title = "Deep split parameter",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        name = "Epithelial grey module ratio")

ht_list = h1 + h2
draw(ht_list, ht_gap = unit(1, "cm"))
dev.off()
  
  
  
#==============================================================================#
#====================== Deepsplit analysis using silwidth =====================#
#==============================================================================#


deepsplits = seq(1,4)
granularity = 30

results.m = data.frame(row.names = seq(1,4))
results.m$WGCNA = 0
results.m$PAM.euclidean = 0
results.m$PAM.cosine = 0
results.m$CLARA = 0
results.m$num = 0

results.h = data.frame(row.names = seq(1,4))
results.h$WGCNA = 0
results.h$PAM.euclidean = 0
results.h$PAM.cosine = 0
results.h$CLARA = 0
results.h$num = 0

for(d in deepsplits){
  
  human.analysis <- perform.WGCNA(human,"epithelial",min_size = granularity, deepsplit = d,type.adj = "signed hybrid")
  murine.analysis <- perform.WGCNA(murine,"stromal",min_size = granularity, deepsplit = d,type.adj = "signed hybrid")
  
  murine.analysis$modules$clustering <- murine.analysis$modules$Label + 1
  human.analysis$modules$clustering <- human.analysis$modules$Label + 1
  
  num.murine <- max(murine.analysis$modules$clustering)
  num.human <- max(human.analysis$modules$clustering)
  
  #============================================================================#
  
  d.m.eval <- Dist(t(murine), method = "pearson")
  d.h.eval <- Dist(t(human), method = "pearson")

  d.m.cosine <- dist.cosine(t(murine))
  d.h.cosine <- dist.cosine(t(human))
  
  d.m.euclidean <- dist(t(murine), method = "euclidean")
  d.h.euclidean <- dist(t(human), method = "euclidean")
  
  #============================================================================#
  
  pam.m.cosine <- pam(d.m.cosine,max(murine.analysis$modules$clustering),diss = TRUE)
  pam.h.cosine <- pam(d.h.cosine,max(human.analysis$modules$clustering),diss = TRUE)
  
  pam.m.euclidean <- pam(d.m.euclidean,max(murine.analysis$modules$clustering),diss = TRUE)
  pam.h.euclidean <- pam(d.h.euclidean,max(human.analysis$modules$clustering),diss = TRUE)
  
  clara.m.euclidean <- clara(t(murine),max(murine.analysis$modules$clustering),metric="euclidean")
  clara.h.euclidean <- clara(t(human),max(human.analysis$modules$clustering),metric="euclidean")
  
  #============================================================================#
  
  stats.m.WGCNA <- cluster.stats(d = d.m.eval, murine.analysis$modules$clustering)
  stats.h.WGCNA <- cluster.stats(d = d.h.eval, human.analysis$modules$clustering)
  
  stats.m.PAM.cosine <- cluster.stats(d = d.m.eval, pam.m.cosine$clustering)
  stats.h.PAM.cosine <- cluster.stats(d = d.h.eval, pam.h.cosine$clustering)
  
  stats.m.PAM.euclidean <- cluster.stats(d = d.m.eval, pam.m.euclidean$clustering)
  stats.h.PAM.euclidean <- cluster.stats(d = d.h.eval, pam.h.euclidean$clustering)
  
  stats.m.CLARA <- cluster.stats(d = d.m.eval, clara.m.euclidean$clustering)
  stats.h.CLARA <- cluster.stats(d = d.h.eval, clara.h.euclidean$clustering)
  
  #============================================================================#
  
  results.m[d,"WGCNA"] = mean(stats.m.WGCNA$clus.avg.silwidths)
  results.m[d,"PAM.cosine"] = mean(stats.m.PAM.cosine$clus.avg.silwidths)
  results.m[d,"PAM.euclidean"] = mean(stats.m.PAM.euclidean$clus.avg.silwidths)
  results.m[d,"CLARA"] = mean(stats.m.CLARA$clus.avg.silwidths)
  results.m[d,"num"] = num.murine
  
  results.h[d,"WGCNA"] = mean(stats.h.WGCNA$clus.avg.silwidths)
  results.h[d,"PAM.cosine"] = mean(stats.h.PAM.cosine$clus.avg.silwidths)
  results.h[d,"PAM.euclidean"] = mean(stats.h.PAM.euclidean$clus.avg.silwidths)
  results.h[d,"CLARA"] = mean(stats.h.CLARA$clus.avg.silwidths)
  results.h[d,"num"] = num.human
  
}

?png
pdf("results/clustering_comparison_sil.pdf", width = 8, height = 6)
par(mfrow = c(1,2))
plot(c(1,2,3,4),results.m$WGCNA,xaxt = "n", pch = 1, col="blue", type = "b",ylim = c(0.0,0.5),
     ylab = "Mean cluster average silhouette width", xlab = "Num modules", main ="Stromal")
lines(c(1,2,3,4),results.m$PAM.cosine,xaxt = "n", pch = 2, col="red", type = "b")
lines(c(1,2,3,4),results.m$PAM.euclidean,xaxt = "n", pch = 4, col="orange", type = "b")
lines(c(1,2,3,4),results.m$CLARA,xaxt = "n", pch = 3, col="darkgreen", type = "b")

legend(cex = 0.75,text.width = 0.8,x = 2.5,y = 0.15,legend = c("WGCNA","PAM-cosine","PAM-euclidean","CLARA-euclidean"),col=c("blue", "red","orange","darkgreen"),pch = c(1,2,4,3))

axis(1,1:4,results.m$num)


plot(c(1,2,3,4),results.h$WGCNA,xaxt = "n", pch = 1, col="blue", type = "b",ylim = c(0.0,0.175),
     ylab = "Mean cluster average silhouette width", xlab = "Num modules", main ="Epithelial")
lines(c(1,2,3,4),results.h$PAM.cosine,xaxt = "n", pch = 2, col="red", type = "b")
lines(c(1,2,3,4),results.h$PAM.euclidean,xaxt = "n", pch = 4, col="orange", type = "b")
lines(c(1,2,3,4),results.h$CLARA,xaxt = "n", pch = 3, col="darkgreen", type = "b")

legend(cex = 0.75,text.width = 0.8,x = 2.5,y = 0.05,legend = c("WGCNA","PAM-cosine","PAM-euclidean","CLARA-euclidean"),col=c("blue", "red","orange","darkgreen"),pch = c(1,2,4,3))

axis(1,1:4,results.h$num)

dev.off()



#==============================================================================#
#====================== Deepsplit analysis using dunn index ===================#
#==============================================================================#


deepsplits = seq(1,4)
granularity = 30

results.m = data.frame(row.names = seq(1,4))
results.m$WGCNA = 0
results.m$PAM.euclidean = 0
results.m$PAM.cosine = 0
results.m$CLARA = 0
results.m$num = 0

results.h = data.frame(row.names = seq(1,4))
results.h$WGCNA = 0
results.h$PAM.euclidean = 0
results.h$PAM.cosine = 0
results.h$CLARA = 0
results.h$num = 0

for(d in deepsplits){
  
  human.analysis <- perform.WGCNA(human,"epithelial",min_size = granularity, deepsplit = d,type.adj = "signed hybrid")
  murine.analysis <- perform.WGCNA(murine,"stromal",min_size = granularity, deepsplit = d,type.adj = "signed hybrid")
  
  murine.analysis$modules$clustering <- murine.analysis$modules$Label + 1
  human.analysis$modules$clustering <- human.analysis$modules$Label + 1
  
  num.murine <- max(murine.analysis$modules$clustering)
  num.human <- max(human.analysis$modules$clustering)
  
  #============================================================================#
  
  d.m.eval <- Dist(t(murine), method = "pearson")
  d.h.eval <- Dist(t(human), method = "pearson")
  
  d.m.cosine <- dist.cosine(t(murine))
  d.h.cosine <- dist.cosine(t(human))
  
  d.m.euclidean <- dist(t(murine), method = "euclidean")
  d.h.euclidean <- dist(t(human), method = "euclidean")
  
  #============================================================================#
  
  pam.m.cosine <- pam(d.m.cosine,max(murine.analysis$modules$clustering),diss = TRUE)
  pam.h.cosine <- pam(d.h.cosine,max(human.analysis$modules$clustering),diss = TRUE)
  
  pam.m.euclidean <- pam(d.m.euclidean,max(murine.analysis$modules$clustering),diss = TRUE)
  pam.h.euclidean <- pam(d.h.euclidean,max(human.analysis$modules$clustering),diss = TRUE)
  
  clara.m.euclidean <- clara(t(murine),max(murine.analysis$modules$clustering),metric="euclidean")
  clara.h.euclidean <- clara(t(human),max(human.analysis$modules$clustering),metric="euclidean")
  
  #============================================================================#
  
  stats.m.WGCNA <- cluster.stats(d = d.m.eval, murine.analysis$modules$clustering)
  stats.h.WGCNA <- cluster.stats(d = d.h.eval, human.analysis$modules$clustering)
  
  stats.m.PAM.cosine <- cluster.stats(d = d.m.eval, pam.m.cosine$clustering)
  stats.h.PAM.cosine <- cluster.stats(d = d.h.eval, pam.h.cosine$clustering)
  
  stats.m.PAM.euclidean <- cluster.stats(d = d.m.eval, pam.m.euclidean$clustering)
  stats.h.PAM.euclidean <- cluster.stats(d = d.h.eval, pam.h.euclidean$clustering)
  
  stats.m.CLARA <- cluster.stats(d = d.m.eval, clara.m.euclidean$clustering)
  stats.h.CLARA <- cluster.stats(d = d.h.eval, clara.h.euclidean$clustering)
  
  #============================================================================#
  
  results.m[d,"WGCNA"] = stats.m.WGCNA$dunn
  results.m[d,"PAM.cosine"] = stats.m.PAM.cosine$dunn
  results.m[d,"PAM.euclidean"] = stats.m.PAM.euclidean$dunn
  results.m[d,"CLARA"] = stats.m.CLARA$dunn
  results.m[d,"num"] = num.murine
  
  results.h[d,"WGCNA"] = stats.h.WGCNA$dunn
  results.h[d,"PAM.cosine"] = stats.h.PAM.cosine$dunn
  results.h[d,"PAM.euclidean"] = stats.h.PAM.euclidean$dunn
  results.h[d,"CLARA"] = stats.h.CLARA$dunn
  results.h[d,"num"] = num.human
  
}

?png
pdf("results/clustering_comparison_dunn.pdf", width = 8, height = 6)
par(mfrow = c(1,2))
plot(c(1,2,3,4),results.m$WGCNA,xaxt = "n", pch = 1, col="blue", type = "b",ylim = c(0.0,0.125),
     ylab = "Dunn index", xlab = "Num modules", main ="Stromal")
lines(c(1,2,3,4),results.m$PAM.cosine,xaxt = "n", pch = 2, col="red", type = "b")
lines(c(1,2,3,4),results.m$PAM.euclidean,xaxt = "n", pch = 4, col="orange", type = "b")
lines(c(1,2,3,4),results.m$CLARA,xaxt = "n", pch = 3, col="darkgreen", type = "b")

legend(cex = 0.75,text.width = 0.8,x = 2.5,y = 0.125,legend = c("WGCNA","PAM-cosine","PAM-euclidean","CLARA-euclidean"),col=c("blue", "red","orange","darkgreen"),pch = c(1,2,4,3))

axis(1,1:4,results.m$num)


plot(c(1,2,3,4),results.h$WGCNA,xaxt = "n", pch = 1, col="blue", type = "b",ylim = c(0.0,0.2),
     ylab = "Dunn index", xlab = "Num modules", main ="Epithelial")
lines(c(1,2,3,4),results.h$PAM.cosine,xaxt = "n", pch = 2, col="red", type = "b")
lines(c(1,2,3,4),results.h$PAM.euclidean,xaxt = "n", pch = 4, col="orange", type = "b")
lines(c(1,2,3,4),results.h$CLARA,xaxt = "n", pch = 3, col="darkgreen", type = "b")

legend(cex = 0.75,text.width = 0.8,x = 2.5,y = 0.2,legend = c("WGCNA","PAM-cosine","PAM-euclidean","CLARA-euclidean"),col=c("blue", "red","orange","darkgreen"),pch = c(1,2,4,3))

axis(1,1:4,results.h$num)

dev.off()



#==============================================================================#
#                             Alluvial Comparison                              #
#==============================================================================#
library(alluvial)
library(dplyr)

# chosen parameters are deepsplit = 3 and granularity = 20 
granularity = 30

human.analysis <- perform.WGCNA(human,"epithelial",min_size = granularity, deepsplit = 3,type.adj = "signed hybrid")
murine.analysis <- perform.WGCNA(murine,"stromal",min_size = granularity, deepsplit = 3,type.adj = "signed hybrid")

murine.analysis$modules$clustering <- murine.analysis$modules$Label + 1
human.analysis$modules$clustering <- human.analysis$modules$Label + 1

d.m.cosine <- dist.cosine(t(murine))
d.h.cosine <- dist.cosine(t(human))

pam.m.cosine <- pam(d.m.cosine,max(murine.analysis$modules$clustering),diss = TRUE)
pam.h.cosine <- pam(d.h.cosine,max(human.analysis$modules$clustering),diss = TRUE)

results.murine <- data.frame(row.names = murine.analysis$modules$Geneid)
results.murine$Geneid <- as.factor(row.names(results.murine))
results.murine$WGCNA <- murine.analysis$modules$Colors
results.murine$PAM <- pam.m.cosine$clustering

results.human <- data.frame(row.names = human.analysis$modules$Geneid)
results.human$Geneid <- as.factor(row.names(results.human))
results.human$WGCNA <- human.analysis$modules$Colors
results.human$PAM <- pam.h.cosine$clustering


results.murine %>% group_by(WGCNA, PAM) %>%
  summarise(Freq = n()) -> results.murine.2d
results.human %>% group_by(WGCNA, PAM) %>%
  summarise(Freq = n()) -> results.human.2d

?alluvial

pdf("results/alluvial_plot.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
alluvial(results.murine.2d[,1:2],
         col = results.murine.2d$WGCNA,
         freq=results.murine.2d$Freq,
         border = "grey",
         hide=results.murine.2d$Freq < 5)
mtext("Stromal", 3, line=3, font=2)

alluvial(results.human.2d[,1:2],
         col = results.human.2d$WGCNA,
         freq=results.human.2d$Freq,
         border = "grey",
         hide=results.human.2d$Freq < 5)
mtext("Epithelial", 3, line=3, font=2)

dev.off()

#==============================================================================#
#             Comparison w/ heatmaps 
#==============================================================================#

library(ComplexHeatmap)

#==============================================================================#
#=================================Stromal======================================#
#==============================================================================#
par(mfrow = c(1,2))


ha1 = rowAnnotation(pam = pam.m.cosine$clustering, 
                   col = list(pam = c( '1' = "brown",
                                       '2' = "red",
                                       '3' = "green",
                                       '4' = "yellow",
                                       '5'= "purple",
                                       '6' = "pink",
                                       '7' = "black",
                                       '8' = "greenyellow",
                                       '9' = "blue",
                                       '10' = "grey",
                                       '11' = "tan",
                                       '12' = "turquoise",
                                       '13' = "magenta")))
ht1 = Heatmap(as.matrix(t(murine)),
        column_title = "Stromal - WGCNA",
        name = "zscore",
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_title_rot = 0,
        row_gap = unit(.5,"mm"),
        row_split =murine.analysis$modules$Colors,
        right_annotation = ha1,
)


ha2 = rowAnnotation(pam = pam.h.cosine$clustering, 
                   col = list(pam = c( '1' = "brown",
                                       '2' = "red",
                                       '3' = "green",
                                       '4' = "yellow",
                                       '5'= "purple",
                                       '6' = "pink",
                                       '7' = "black",
                                       '8' = "greenyellow",
                                       '9' = "blue",
                                       '10' = "grey",
                                       '11' = "tan",
                                       '12' = "turquoise")))
ht2 = Heatmap(as.matrix(t(human)),
              column_title = "Epithelial - WGCNA",
              name = "zscore",
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              show_column_names = FALSE,
              show_row_names = FALSE,
              row_title_rot = 0,
              row_gap = unit(.5,"mm"),
              row_split =human.analysis$modules$Colors,
              right_annotation = ha2,
)

pdf("results/risultati xs/clustering comparison/stromal_heatmap.pdf")
draw(ht1)
dev.off()
pdf("results/risultati xs/clustering comparison/epithelial_heatmap.pdf")
draw(ht2)
dev.off()

