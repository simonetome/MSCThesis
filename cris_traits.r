
cris.results = readRDS("C:/Users/simon/Tesi_tome/datasets/robjects/NTP_cpm_pdx_reference.rds")
cris.results = cris.results$result
row.names(cris.results) = cris.results$aliquot_id

samples = cris.results$aliquot_id
epic_result = (readRDS("populations/EPIC_output_no.rds"))$cellFractions %>% as.data.frame()
commonIds = samples[samples %in% row.names(epic_result)]
epic_filtered = epic_result[commonIds,]

epic_filtered$CRIS_label = cris.results[commonIds,"predict.label2"]
celltypes = colnames(epic_filtered)[colnames(epic_filtered) != "CRIS_label"]
celltypes

View(epic_filtered)

for(ct in celltypes){
  epic_filtered[c(ct,"CRIS_label")] %>% 
    ggplot(aes(x=CRIS_label, y=get(ct),fill=CRIS_label)) + 
    geom_boxplot() + 
    geom_jitter() +
    labs(y = as.character(ct)) + 
    ggtitle(paste("CRIS to cell fraction distribution",as.character(ct))) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste("figs/epic_",ct,".png",sep=""),
         width = 8,
         height = 6,)
}

#==============================================================================#
# Candiolo signatures 
#==============================================================================#

candiolo = readRDS("populations/candiolo_tpm_output.rds")
colnames(candiolo) = c("Endothelial","Leucocyte","CAF")
commonIds = samples[samples %in% row.names(candiolo)]
candiolo_filtered = candiolo[commonIds,] %>% as.data.frame()
candiolo_filtered$CRIS_label = cris.results[commonIds,"predict.label2"]
celltypes = colnames(candiolo_filtered)[colnames(candiolo_filtered) != "CRIS_label"]


for(ct in celltypes){
  candiolo_filtered[c(ct,"CRIS_label")] %>% 
    ggplot(aes(x=CRIS_label, y=get(ct),fill=CRIS_label)) + 
    geom_boxplot() + 
    geom_jitter() +
    labs(y = as.character(ct)) + 
    ggtitle(paste("CRIS to cell fraction distribution",as.character(ct))) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste("figs/candioloCRIS_",ct,".png",sep=""),
         width = 8,
         height = 6,)
}













