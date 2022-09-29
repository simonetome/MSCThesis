#==============================================================================#
# Data preparing and dumping 
#==============================================================================#

source("rutils/data_loading.r")
library(Biobase)

assembly = read.csv("datasets/filtered_assembly.csv")
samples = colnames(assembly)[3:length(colnames(assembly))]
row.names(assembly) = assembly$Geneid
assembly.length = assembly[assembly$Length >= 400,]
assembly = assembly[samples]
assembly.length = assembly.length[samples]

assembly = split_dataset(assembly)
assembly.length = split_dataset(assembly.length)

saveRDS(assembly,"datasets/robjects/full_assembly.rds")
saveRDS(assembly.length,"datasets/robjects/full_assembly_gt400.rds")









