# wrapper for functional enrichment with gprofiler having wgcna modules as input
setwd("C:/Users/simon/MSCThesis")
library(gprofiler2)
library(xlsx)
library(dplyr)

homologs = read.csv("C:\\Users\\simon\\MSCThesis\\datasets\\orthologs/Hortolog_HS_MM.txt",sep = "\t")

map_homolog <- function(x){
  if(startsWith(x,"H")){
    return(x)
  }
  else{
    homologs %>% filter(mmusculus_GeneSymbol == x) %>% select(hsapiens_GeneSymbol) -> temp
    if(nrow(temp) > 0){
      return(temp$hsapiens_GeneSymbol)
    }else{
      return("null")
    }
  }
}


functional_enrichment <- function(modules,
                                  fileName,
                                  prefixLength,
                                  homolog_mapping){
  
  counter = 1
  for(c in c(unique(modules$Colors))){
    print(c)
    modules %>% filter(Colors == c) %>% 
      select(Geneid) %>% 
      rowwise %>% 
      mutate(Geneid_mapped = map_homolog(Geneid)) %>%
      filter(Geneid_mapped != "null") %>%
      rowwise %>%
      mutate(Geneid_mapped = substr(Geneid_mapped,prefixLength+1,nchar(Geneid_mapped))) %>%
      as.data.frame -> geneList.df
    geneList = geneList.df$Geneid_mapped
    
    analysis = gost(query = geneList, 
                    organism = "hsapiens", 
                    ordered_query = FALSE,
                    correction_method = "g_SCS")$result
    if(counter == 1){
      append= FALSE
    }else{
      append = TRUE
    }
    write.xlsx2(
      x = analysis,
      file = fileName,
      sheetName = c,
      append = append
    )
    counter = counter + 1
  }
}












unique(modules$Colors)



for(c in c(unique(modules$Colors))){
  print(c)
}









