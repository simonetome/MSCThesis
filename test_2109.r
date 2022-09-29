setwd("C:/Users/simon/Tesi_tome")
source("rutils/wgcna_wrapper.r")
source("rutils/data_loading.r")
library(stringr)
library(stylo)
library(dplyr)
library(ComplexHeatmap)

# This script has the aim to optimize the correlation between murine eigengenes and 
# stromal population --> separate the two genomes 
# find optimal wgcna for murine 

assembly = readRDS("datasets/robjects/full_assembly_gt400.rds")
murine = top_variance(assembly[[2]],3000)
sample.traits = readRDS("populations/candiolo_tpm_output.rds")

cuts = c(0.8,0.85,0.9,0.95)
min_sizes = c(10,15,20,25)
deepsplits = c(1,2,3,4)

# cuts = c(0.8)
# min_sizes = c(10)
# deepsplits = c(1)

eigengenes = vector(mode = "list", length = 4)
names(eigengenes) = c("cuts","min_size","deepsplits","eigengenes")
eigengenes[[1]] = vector(length = 64, mode = "numeric")
eigengenes[[2]] = vector(length = 64, mode = "numeric")
eigengenes[[3]] = vector(length = 64, mode = "numeric")
eigengenes[[4]] = vector(length = 64, mode = "list")


counter = 1

for(c in cuts){
  for(m in min_sizes){
    for(d in deepsplits){
      rsquared_cut = c
      min_clust_size = m
      deepSplit = d
      
      murine.WGCNA.analysis = perform.WGCNA(dataset = t(murine),
                                            title = "WGCNA_output/murine",
                                            min_size = min_clust_size,
                                            deepsplit = deepSplit,
                                            cortype = "spearman",
                                            type.adj = "signed hybrid",
                                            RsquaredCut = rsquared_cut)
      
      murine.eigengenes = as.data.frame(murine.WGCNA.analysis$eigengenes$eigengenes)
    
      eigengenes$cuts[counter] = c
      eigengenes$min_size[counter] = m
      eigengenes$deepsplits[counter] = d
      eigengenes$eigengenes[[counter]] = murine.eigengenes

      
      counter = counter + 1
    }
  }
}
saveRDS(eigengenes,"datasets/robjects/eigWGCNA_murine_top3000_gt400.rds")
counter = 1

eigengenes = readRDS("datasets/robjects/eigWGCNA_murine_top3000_gt400.rds")
#==============================================================================#
# Candiolo results 
#==============================================================================#

results = as.data.frame(matrix(nrow = 64, ncol = 6))
colnames(results) = c("cuts","min_size","deepsplits","Endothelial_max","Leucocyte_max","CAF_max")

for(counter in c(1:64)){
  
  c = eigengenes$cuts[counter]
  m = eigengenes$min_size[counter]
  d = eigengenes$deepsplits[counter]
  
  eig = eigengenes$eigengenes[[counter]]
  cor1 = cor(eig,sample.traits)
  results[counter,"cuts"] = c
  results[counter,"min_size"] = m
  results[counter,"deepsplits"] = d
  results[counter,"Endothelial_max"] = max(cor1[,1])
  results[counter,"Leucocyte_max"] = max(cor1[,2])
  results[counter,"CAF_max"] = max(cor1[,3])
  
}


saveRDS(results,"WGCNAcorrbenchmark/WGCNA_benchmark_candiolo_murine_check.rds")
# Parameters for WGCNA 

#==============================================================================#
# test murine based deconvolutions 
#==============================================================================#


cuts = c(0.8,0.85,0.9,0.95)
min_sizes = c(10,15,20,25)
deepsplits = c(1,2,3,4)

#mmcp_counter
#seqimmucc
#dcq
#base
sample.traits.cor = vector(mode = "list", length = 4)
names(sample.traits.cor) = c("mmcp_output.rds",
                             "seqi.llsr_output.rds",
                             "dcq_output.rds",
                             "base_output.rds")

for(s in names(sample.traits.cor)){
  sample.traits = readRDS(paste("populations/",s,sep=""))
  row.names(sample.traits) = sample.traits$cell_type
  sample.traits$cell_type = NULL
  sample.traits = t(sample.traits)
  
  
  results = as.data.frame(matrix(nrow = 64, ncol = length(colnames(sample.traits)) + 3))
  colnames(results) = c("cuts","min_size","deepsplits",colnames(sample.traits))
  
  for(counter in c(1:64)){
    c = eigengenes$cuts[counter]
    m = eigengenes$min_size[counter]
    d = eigengenes$deepsplits[counter]
    eig = eigengenes$eigengenes[[counter]]
    
    cor1 = cor(eig,sample.traits)
    results[counter,"cuts"] = c
    results[counter,"min_size"] = m
    results[counter,"deepsplits"] = d
    for(i in c(1:length(colnames(sample.traits))))
    {
      results[counter,colnames(sample.traits)[i]] = max(cor1[,i])
    }
  }
  sample.traits.cor[[s]] = results
  
}

saveRDS(sample.traits.cor,"WGCNAcorrbenchmark/WGCNA_benchmark_murine_based_check.rds")

View(sample.traits.cor$mmcp_output.rds)
View(sample.traits.cor$seqi.llsr_output.rds)
View(sample.traits.cor$dcq_output.rds)
View(sample.traits.cor$base_output.rds)


#==============================================================================#
# test human based deconvolutions 
#==============================================================================#

# cuts = c(0.8)
# min_sizes = c(10)
# deepsplits = c(1)

sample.traits.cor = vector(mode = "list", length = 7)
names(sample.traits.cor) = c("estimate_output.rds",
                             "consensus_tme_output.rds",
                             "abis_output.rds",
                             "xcell_output.rds",
                             "mcp_counter_output.rds",
                             "timer_output.rds",
                             "quantiseq_output.rds")

for(s in names(sample.traits.cor)){
  sample.traits = readRDS(paste("populations/",s,sep="")) %>% as.data.frame()
  row.names(sample.traits) = sample.traits$cell_type
  sample.traits$cell_type = NULL
  sample.traits = t(sample.traits)
  
  
  results = as.data.frame(matrix(nrow = 64, ncol = length(colnames(sample.traits)) + 3))
  colnames(results) = c("cuts","min_size","deepsplits",colnames(sample.traits))
  
  for(counter in c(1:64)){
    c = eigengenes$cuts[counter]
    m = eigengenes$min_size[counter]
    d = eigengenes$deepsplits[counter]
    eig = eigengenes$eigengenes[[counter]]
    
    cor1 = cor(eig,sample.traits)
    results[counter,"cuts"] = c
    results[counter,"min_size"] = m
    results[counter,"deepsplits"] = d
    for(i in c(1:length(colnames(sample.traits))))
    {
      results[counter,colnames(sample.traits)[i]] = max(cor1[,i])
    }
  }
  sample.traits.cor[[s]] = results
  
}

saveRDS(sample.traits.cor,"WGCNAcorrbenchmark/WGCNA_benchmark_human_based_check.rds")


#==============================================================================#

epic = (readRDS("populations/EPIC_output.rds"))$cellFractions %>% as.data.frame()
epic.no = (readRDS("populations/EPIC_output_no.rds"))$cellFractions %>% as.data.frame()


sample.traits = epic
results = as.data.frame(matrix(nrow = 64, ncol = length(colnames(sample.traits)) + 3))
colnames(results) = c("cuts","min_size","deepsplits",colnames(sample.traits))

for(counter in c(1:64)){
  c = eigengenes$cuts[counter]
  m = eigengenes$min_size[counter]
  d = eigengenes$deepsplits[counter]
  eig = eigengenes$eigengenes[[counter]]
  
  cor1 = cor(eig,sample.traits)
  results[counter,"cuts"] = c
  results[counter,"min_size"] = m
  results[counter,"deepsplits"] = d
  for(i in c(1:length(colnames(sample.traits))))
  {
    results[counter,colnames(sample.traits)[i]] = max(cor1[,i])
  }
}

saveRDS(results,"WGCNAcorrbenchmark/WGCNA_benchmark_epic.rds")

sample.traits = epic.no
results = as.data.frame(matrix(nrow = 64, ncol = length(colnames(sample.traits)) + 3))
colnames(results) = c("cuts","min_size","deepsplits",colnames(sample.traits))

for(counter in c(1:64)){
  c = eigengenes$cuts[counter]
  m = eigengenes$min_size[counter]
  d = eigengenes$deepsplits[counter]
  eig = eigengenes$eigengenes[[counter]]
  
  cor1 = cor(eig,sample.traits)
  results[counter,"cuts"] = c
  results[counter,"min_size"] = m
  results[counter,"deepsplits"] = d
  for(i in c(1:length(colnames(sample.traits))))
  {
    results[counter,colnames(sample.traits)[i]] = max(cor1[,i])
  }
}

saveRDS(results,"WGCNAcorrbenchmark/WGCNA_benchmark_epic_no.rds")




