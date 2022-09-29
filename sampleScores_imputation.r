setwd("C:/Users/simon/Tesi_tome")
library("stringr")
library("stylo")
library(fpc)
library(cluster)
library(amap)
library(magrittr)
library(dplyr)
library(EPIC)
library(immunedeconv)
library(ComplexHeatmap)

tpm.assembly <- readRDS("datasets/robjects/full_assembly.rds")
human <- tpm.assembly[[1]]
murine <- tpm.assembly[[2]]

# need to convert murine genes into human genes 
gene.annotations <- read.csv(file = "datasets/orthologs/Gene_annotations.csv",sep=";")
gene.orthologs <- read.csv(file = "datasets/orthologs/Hortolog_HS_MM.txt", sep = "\t")


convert.to.h <- function(x){
  gene.orthologs %>% filter(mmusculus_GeneSymbol == x) %>% select(hsapiens_GeneSymbol) -> temp
  if(nrow(temp) > 0){
    return(substr(temp$hsapiens_GeneSymbol,3,nchar(temp$hsapiens_GeneSymbol)))
  }else{
    return("null")
  }
}

murine %>%
  mutate(Geneid = row.names(murine)) %>%
  rowwise() %>% 
  mutate(Geneid = convert.to.h(Geneid)) %>% 
  filter(Geneid != "null") -> murine.converted

murine.converted = as.data.frame(murine.converted)
row.names(murine.converted) = murine.converted$Geneid
murine.converted$Geneid = NULL

#==============================================================================#
#================================    EPIC     =================================#
#==============================================================================#
epic.out = EPIC(bulk = murine.converted, reference = "TRef",withOtherCells = TRUE)
epic.out.n = EPIC(bulk = murine.converted, reference = "TRef",withOtherCells = FALSE)
saveRDS(epic.out,"populations/EPIC_output.rds")
saveRDS(epic.out.n,"populations/EPIC_output_no.rds")


#==============================================================================#
#================================    xcell     ================================#
#==============================================================================#

xcell = immunedeconv::deconvolute(murine.converted,"xcell")
saveRDS(xcell,"populations/xcell_output.rds")

#==============================================================================#
#================================    quantiseq     ============================#
#==============================================================================#

quantiseq = immunedeconv::deconvolute(murine.converted,"quantiseq")
saveRDS(quantiseq,"populations/quantiseq_output.rds")

#==============================================================================#
#================================    mcp_counter     ==========================#
#==============================================================================#

mcp_counter = immunedeconv::deconvolute(murine.converted,"mcp_counter")
saveRDS(mcp_counter,"populations/mcp_counter_output.rds")

#==============================================================================#
#================================    abis     ==========================#
#==============================================================================#

abis = immunedeconv::deconvolute(murine.converted,"abis")
saveRDS(abis,"populations/abis_output.rds")

#==============================================================================#
#================================    estimate     ==========================#
#==============================================================================#

estimate = immunedeconv::deconvolute(murine.converted,"estimate")
saveRDS(estimate,"populations/estimate_output.rds")

#==============================================================================#
#================================    Timer       ==========================#
#==============================================================================#

indications = c(1:624)
indications[c(1:624)] = "COAD"
timer = immunedeconv::deconvolute(murine.converted,"timer",indications = indications)
saveRDS(timer,"populations/timer_output.rds")

#==============================================================================#
#================================    consensus_tme     ==========================#
#==============================================================================#

consensus_tme = deconvolute(murine.converted,"consensus_tme", indications = indications)
saveRDS(consensus_tme,"populations/consensus_tme_output.rds")


# leverage murine based
#==============================================================================#
#================================    murine     ===============================#
#==============================================================================#

# truncate prefix "M_"
row.names(murine) = substr(row.names(murine),3,nchar(row.names(murine)))

base.output = as.data.frame(deconvolute_mouse(murine,"base"))
mmcp.output = as.data.frame(deconvolute_mouse(murine,"mmcp_counter"))
seqi.llsr.output = as.data.frame(deconvolute_mouse(murine,"seqimmucc",algorithm = "LLSR"))
dcq.output = as.data.frame(deconvolute_mouse(murine,"dcq"))

saveRDS(base.output,"populations/base_output.rds")
saveRDS(mmcp.output,"populations/mmcp_output.rds")
saveRDS(seqi.llsr.output,"populations/seqi.llsr_output.rds")
saveRDS(dcq.output,"populations/dcq_output.rds")

#==============================================================================#
#================================    Candiolo sig.   ==========================#
#==============================================================================#
# not all genes in the signature are present 

library(xlsx)

signatures = read.xlsx("signatures.xlsx",1)
row.names(signatures) = signatures$Gene.Symbol
signatures$Gene.Symbol = NULL
signatures$Signature = NULL

# common genes 
sigID = row.names(signatures)[row.names(signatures) %in% row.names(murine.converted)]
# simple matrix multiplication
scores_tpm = t(as.matrix(murine.converted[sigID,])) %*% as.matrix(signatures[sigID,]) 
saveRDS(scores_tpm,"populations/candiolo_tpm_output.rds")

saveRDS(murine.converted,"datasets/robjects/murine_converted_full.rds")



#==============================================================================#
#============================ Compare results =================================#
#==============================================================================#

?apply
# zscore columns, representing cell types
h1 = Heatmap((apply(scores_tpm,2,zs)), show_column_names = FALSE, show_row_names = FALSE)
h2 = Heatmap((apply(scores_zs,2,zs)), show_column_names = FALSE, show_row_names = FALSE)
draw(h1+h2)



















