setwd("C:/Users/simon/Tesi_tome")
library(dplyr)
library(gprofiler2)
library(ComplexHeatmap)
library(gdata)
source("rutils/data_loading.r")
source("wgcna_wrapper.r")
library(xlsx)

ligand.receptor.db = read.csv("interactions/ligand_receptor_strict.csv", sep = ";")

ligand.receptor.db  %>% 
  filter(LigandM != "") %>% 
  select(ReceptorH,LigandM) %>% 
  mutate(LigandM = paste("M_",LigandM,sep=""), ReceptorH = paste("H_",ReceptorH,sep =""))-> murine.ligand

ligand.receptor.db  %>% 
  filter(ReceptorM != "") %>% 
  select(LigandH,ReceptorM) %>%
  mutate(LigandH = paste("H_",LigandH,sep=""), ReceptorM = paste("M_",ReceptorM,sep =""))-> human.ligand


zscored <- load_zscored("datasets/filtered_assembly.csv", top.num = NULL)
human <- as.data.frame(zscored[1])
murine <- as.data.frame(zscored[2])

ligand.receptor.db %>%
  mutate(LigandH = paste("H_",LigandH,sep=""),
         ReceptorH = paste("H_",ReceptorH,sep=""),
         LigandM = paste("M_",LigandM,sep=""),
         ReceptorM = paste("M_",ReceptorM,sep="")) -> ligand.receptor.db


(rbind(human,murine) -> temp)  %>% 
  mutate(Geneid = row.names(temp)) %>%
  filter(Geneid %in% human.ligand$LigandH | Geneid %in% human.ligand$ReceptorM) %>%
  select(!Geneid) -> human.ligand.exp

(rbind(human,murine) -> temp)  %>% 
  mutate(Geneid = row.names(temp)) %>%
  filter(Geneid %in% murine.ligand$LigandM | Geneid %in% murine.ligand$ReceptorH) %>%
  select(!Geneid) -> murine.ligand.exp

Heatmap(murine.ligand.exp, show_row_names = FALSE, show_column_names = FALSE)


print("HUMAN LIGAND")
human.ligand.analysis = perform.WGCNA(t(human.ligand.exp),"results/human_ligand",30,4,"unsigned")

#windows()
#human.ligand.analysis$modules %>%
#  group_by(Colors) %>%
#  summarise(num = n())-> temp
#barplot(temp$num, names.arg = temp$Colors,las = 2, horiz = TRUE, col = temp$Colors)

# how many interactions in the modules 

modules = human.ligand.analysis$modules
for(c in unique(modules$Colors)){
  
  modules %>%
    filter(Colors == c) %>% 
    select(Geneid) -> geneList
  geneList = (as.list(geneList$Geneid))
  ligand.receptor.db %>%
    filter(LigandH %in% geneList & ReceptorM %in% geneList) -> temp
  print(c)
  print(temp)
  print(nrow(temp))
}


print("MURINE LIGAND")
murine.ligand.analysis = perform.WGCNA(t(murine.ligand.exp),"results/murine_ligand",30,4,"unsigned")
modules = murine.ligand.analysis$modules
for(c in unique(modules$Colors)){
  
  modules %>%
    filter(Colors == c) %>% 
    select(Geneid) -> geneList
  
  geneList = (as.list(geneList$Geneid))
  ligand.receptor.db %>%
    filter(LigandM %in% geneList & ReceptorH %in% geneList) -> temp
  print(c)
  print(temp)
  print(nrow(temp))
}













