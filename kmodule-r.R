setwd("C:/Users/simon/Tesi_tome")
source("wgcna_wrapper.r")
source("rutils/data_loading.r")
library("stringr")
library("stylo")
library(fpc)
library(cluster)
library(amap)

zscored <- load_zscored("datasets/filtered_assembly.csv", top.num = 1000)
human <- as.data.frame(zscored[1])
murine <- as.data.frame(zscored[2])

?Dist
d.m.eval <- Dist(murine, method = "pearson")
d.h.eval <- Dist(human, method = "pearson")


human.WGCNA.analysis = perform.WGCNA(t(human),"1",20,4,"signed hybrid")
murine.WGCNA.analysis = perform.WGCNA(t(murine),"2",20,4,"signed hybrid")

human.KMODULE = perform.KMODULE(t(human),human.WGCNA.analysis,"3")
murine.KMODULE = perform.KMODULE(t(murine),murine.WGCNA.analysis,"4")

human.WGCNA.analysis$modules$clustering =  human.WGCNA.analysis$modules$Label + 1
murine.WGCNA.analysis$modules$clustering =  murine.WGCNA.analysis$modules$Label + 1
human.KMODULE$modules$clustering =  human.KMODULE$modules$Label + 1
murine.KMODULE$modules$clustering =  murine.KMODULE$modules$Label + 1

stat.h.wgcna = cluster.stats(d = d.h.eval, human.WGCNA.analysis$modules$clustering)
stat.m.wgcna = cluster.stats(d = d.m.eval, murine.WGCNA.analysis$modules$clustering)

stat.h.k = cluster.stats(d = d.h.eval, human.KMODULE$modules$clustering)
stat.m.k = cluster.stats(d = d.m.eval, murine.KMODULE$modules$clustering)

print(mean(stat.h.wgcna$clus.avg.silwidths))
print(mean(stat.h.k$clus.avg.silwidths))

print(mean(stat.m.wgcna$clus.avg.silwidths))
print(mean(stat.m.k$clus.avg.silwidths))


length(human.WGCNA.analysis$modules$clustering)





