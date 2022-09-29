setwd("C:/Users/simon/Tesi_tome")
source("wgcna_wrapper.r")
source("rutils/data_loading.r")
library("stringr")
library("stylo")
library(fpc)
library(cluster)
library(xlsx)

zscored <- load_zscored("datasets/filtered_assembly.csv", top.num = 1000)
human <- t(as.data.frame(zscored[1]))
murine <- t(as.data.frame(zscored[2]))

colnames(human) <- unlist(lapply(colnames(human),function(x){return(substr(x,3,nchar(x)))}))
colnames(murine) <- unlist(lapply(colnames(murine),function(x){return(substr(x,3,nchar(x)))}))


human.analysis <- perform.WGCNA(human,"results/epithelial",min_size = 20, deepsplit = 3, type.adj = "signed hybrid")
murine.analysis <- perform.WGCNA(murine,"results/stromal",min_size = 20, deepsplit = 3,type.adj = "signed hybrid")

murine.analysis$modules$clustering <- murine.analysis$modules$Label + 1
human.analysis$modules$clustering <- human.analysis$modules$Label + 1

?gost
#==============================================================================#
#                             Kmodule optimization                             #
#==============================================================================#

#human.optimization <- perform.KMODULE(human,analysis,"epithelial")
#murine.optimization <- perform.KMODULE(murine,analysis,"stromal")


#==============================================================================#
#                             Functional enrichment                            #
#==============================================================================#

  
fun.analysis.murine <- vector(mode = "list",length = length(unique(murine.analysis$modules$Colors)))
fun.analysis.human <- vector(mode = "list",length = length(unique(human.analysis$modules$Colors)))
names(fun.analysis.murine) <- unique(murine.analysis$modules$Colors)
names(fun.analysis.human) <- unique(human.analysis$modules$Colors)

for(c in unique(murine.analysis$modules$Colors)){
  fun.analysis.murine[[c]] <- 
    gost(query = murine.analysis$modules[murine.analysis$modules$Colors == c, "Geneid"],organism = "mmusculus", domain_scope = "annotated") 
}

for(c in unique(human.analysis$modules$Colors)){
  fun.analysis.human[[c]] <- 
    gost(query = human.analysis$modules[human.analysis$modules$Colors == c, "Geneid"],organism = "hsapiens", domain_scope = "annotated") 
}

counter = 0

for(c in unique(murine.analysis$modules$Colors)){
  fun.analysis.murine[[c]]$result$color = c
  fun.analysis.murine[[c]]$result$parents = NULL
  to_dump = fun.analysis.murine[[c]]$result

  if(counter == 0){
    write.table(to_dump, "results/fun_analysis_murine_WGCNA.csv", sep = ",",col.names = TRUE)
    write.xlsx(to_dump, "results/fun_analysis_murine_WGCNA_xcel.xlsx", sheetName = c)
  }else{
    write.table(to_dump, "results/fun_analysis_murine_WGCNA.csv", sep = ",",col.names = FALSE, append = TRUE)
    write.xlsx(to_dump, "results/fun_analysis_murine_WGCNA_xcel.xlsx", append = TRUE, sheetName = c)
  }
  counter = counter + 1
}

counter = 0
for(c in unique(human.analysis$modules$Colors)){
  fun.analysis.human[[c]]$result$color = c
  fun.analysis.human[[c]]$result$parents = NULL
  
  to_dump = fun.analysis.human[[c]]$result

  
  
  if(counter == 0){
    write.table(to_dump, "results/fun_analysis_human_WGCNA.csv", sep = ",", col.names = TRUE)
    write.xlsx(to_dump, "results/fun_analysis_human_WGCNA_xcel.xlsx", sheetName = c, col.names = TRUE)
  }else{
    write.table(to_dump, "results/fun_analysis_human_WGCNA.csv", sep = ",",col.names = FALSE, append = TRUE)
    write.xlsx(to_dump, "results/fun_analysis_human_WGCNA_xcel.xlsx", sheetName = c, col.names = TRUE, append = TRUE)
  }
  counter = counter + 1
}


write.table(murine.analysis$modules, "results/modules_murine_WGCNA.csv", sep = ",")
write.table(human.analysis$modules, "results/modules_human_WGCNA.csv", sep = ",")

p <- gostplot(fun.analysis.murine$blue, capped = FALSE, interactive = TRUE)
p

?




dev.off()
#==============================================================================#


?Heatmap


