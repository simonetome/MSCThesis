source("rutils/wgcna_wrapper.r")
source("rutils/data_loading.r")

library(dplyr)

full.asssembly = readRDS("datasets/robjects\\full_assembly_gt400.rds")
human = top_variance(full.asssembly[[1]],3000)
murine = top_variance(full.asssembly[[2]],3000)
full = rbind(human,murine)

cuts = c(0.8,0.85,0.9)
min_clust_sizes = c(20,30,40)
deepSplits = c(1,2,3,4)


results = as.data.frame(matrix(nrow = 36, ncol = 4))
colnames(results) = c("cut","minsize","deepsplit","moduleNumber")

counter = 1

for(c in cuts){
  for(m in min_clust_sizes){
    for(d in deepSplits){
      
      wgcna.analysis = perform.WGCNA(dataset = t(full),
                                     title = "WGCNA_output/full",
                                     min_size = m,
                                     deepsplit = d,
                                     cortype = "spearman",
                                     type.adj = "signed hybrid",
                                     RsquaredCut = c)
      
      results[counter,"cut"] = c
      results[counter,"minsize"] = m
      results[counter,"deepsplit"] = d
      
      
      wgcna.analysis$modules %>% as.data.frame() -> modules
      modules$Label = modules$Label + 1
      
      balance = as.data.frame(matrix(nrow = length(unique(modules$Label)), ncol = 2))
      colnames(balance) = c("Human.percentage","Murine.percentage")
      
      for(i in c(1:length(unique(modules$Label)))){
        tmp = modules[modules$Label == i,]
        total_count = nrow(tmp)
        tmp %>% rowwise() %>% filter(startsWith(Geneid,"H")) %>% as.data.frame() -> h
        balance[i,"Human.percentage"] = nrow(h)/total_count
        balance[i,"Murine.percentage"] = 1 - balance[i,"Human.percentage"]
      }
      
      l = row.names(balance[balance["Human.percentage"] >= 0.25 & balance["Human.percentage"] <= 0.75,])
      results[counter,"moduleNumber"] = length(l)
      
      counter = counter + 1
      
    }
  }
}


View(results)
saveRDS(results,"fullWGCNA6000_nummodule_gt025_balance.rds")

c = 0.85
m = 30
d = 1

wgcna.analysis = perform.WGCNA(dataset = t(full),
                               title = "WGCNA_output/full",
                               min_size = m,
                               deepsplit = d,
                               cortype = "spearman",
                               type.adj = "signed hybrid",
                               RsquaredCut = c)

wgcna.analysis$modules %>% as.data.frame() -> modules
modules$Label = modules$Label + 1

balance = as.data.frame(matrix(nrow = length(unique(modules$Label)), ncol = 2))
colnames(balance) = c("Human.percentage","Murine.percentage")

for(i in c(1:length(unique(modules$Label)))){
  tmp = modules[modules$Label == i,]
  total_count = nrow(tmp)
  tmp %>% rowwise() %>% filter(startsWith(Geneid,"H")) %>% as.data.frame() -> h
  balance[i,"Human.percentage"] = nrow(h)/total_count
  balance[i,"Murine.percentage"] = 1 - balance[i,"Human.percentage"]
}

ggplot(balance, aes(x=Human.percentage)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity",fill="red") 
  

# extract modules where balance is at least 25 for both:

l = row.names(balance[balance["Human.percentage"] >= 0.25 & balance["Human.percentage"] <= 0.75,])
interesting.colors = modules[modules$Label %in% l,"Colors"]
interesting.colors = paste("ME",unique(interesting.colors),sep="")
eigengenes = wgcna.analysis$eigengenes$eigengenes
interesting.eigengenes = eigengenes[interesting.colors]

candiolo.results = readRDS("populations/candiolo_tpm_output.rds")
colnames(candiolo.results) = c("Endothelial","Leucocyte","CAF")



#==============================================================================#
#============================Candiolo vs EPIC no===============================#
#==============================================================================#
png("figs/WGCNAfull_vs/EPICno.png", width = 700, height = 600)
c1 = cor(candiolo.results,interesting.eigengenes,method="spearman")
h1 = Heatmap(t(c1),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Candiolo \n max:",round(max(c1),digits = 2)),
             column_title = "Candiolo")


epic.results = readRDS("populations/epic_output_no.rds")$cellFractions
c2 = cor(epic.results,interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n EPIC \n max:",round(max(c2),digits = 2)),
             column_title = "EPIC")
draw(h1+h2)
dev.off()

#==============================================================================#
#============================Candiolo vs Consensus=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/consensus.png", width = 700, height = 600)
consensus.results = readRDS("populations/consensus_tme_output.rds") %>% as.data.frame 
row.names(consensus.results) = consensus.results$cell_type
consensus.results$cell_type = NULL 

c2 = cor(t(consensus.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n consensus \n max:",round(max(c2),digits = 2)),
             column_title = "Consensus")
draw(h1+h2)
dev.off()

#==============================================================================#
#============================Candiolo vs Base=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/bases.png", width = 700, height = 600)
base.results = readRDS("populations/base_output.rds") %>% as.data.frame 
row.names(base.results) = base.results$cell_type
base.results$cell_type = NULL 

c2 = cor(t(base.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Base \n max:",round(max(c2),digits = 2)),
             column_title = "Base")
draw(h1+h2)
dev.off()
#==============================================================================#
#============================Candiolo vs Abis=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/abis.png", width = 700, height = 600)
abis.results = readRDS("populations/abis_output.rds") %>% as.data.frame 
row.names(abis.results) = abis.results$cell_type
abis.results$cell_type = NULL 

c2 = cor(t(abis.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Abis \n max:",round(max(c2),digits = 2)),
             column_title = "Abis")
draw(h1+h2)
dev.off()
#==============================================================================#
#============================Candiolo vs DCQ=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/dcq.png", width = 700, height = 600)
dcq.results = readRDS("populations/dcq_output.rds") %>% as.data.frame 
row.names(dcq.results) = dcq.results$cell_type
dcq.results$cell_type = NULL 

c2 = cor(t(dcq.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n DCQ \n max:",round(max(c2),digits = 2)),
             column_title = "DCQ")
draw(h1+h2)
dev.off()
#==============================================================================#
#============================Candiolo vs Estimate=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/estimates.png", width = 700, height = 600)
estimate.results = readRDS("populations/estimate_output.rds") %>% as.data.frame 
row.names(estimate.results) = estimate.results$cell_type
estimate.results$cell_type = NULL 

c2 = cor(t(estimate.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Estimate \n max:",round(max(c2),digits = 2)),
             column_title = "Estimate")
draw(h1+h2)
dev.off()
#==============================================================================#
#============================Candiolo vs Mcp=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/mcp.png", width = 700, height = 600)
mcp.results = readRDS("populations/mcp_counter_output.rds") %>% as.data.frame 
row.names(mcp.results) = mcp.results$cell_type
mcp.results$cell_type = NULL 

c2 = cor(t(mcp.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n MCP \n max:",round(max(c2),digits = 2)),
             column_title = "MCP")
draw(h1+h2)
dev.off()
#==============================================================================#
#============================Candiolo vs Mmcp=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/mmcp.png", width = 700, height = 600)
mmcp.results = readRDS("populations/mmcp_output.rds") %>% as.data.frame 
row.names(mmcp.results) = mmcp.results$cell_type
mmcp.results$cell_type = NULL 

c2 = cor(t(mmcp.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n MMCP \n max:",round(max(c2),digits = 2)),
             column_title = "MMCP")
draw(h1+h2)
dev.off()
#==============================================================================#
#============================Candiolo vs Quantiseq=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/quantiseq.png", width = 700, height = 600)
quantiseq.results = readRDS("populations/quantiseq_output.rds") %>% as.data.frame 
row.names(quantiseq.results) = quantiseq.results$cell_type
quantiseq.results$cell_type = NULL 

c2 = cor(t(quantiseq.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Quantiseq \n max:",round(max(c2),digits = 2)),
             column_title = "Quantiseq")
draw(h1+h2)
dev.off()

#==============================================================================#
#============================Candiolo vs SeqiMMUC=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/seqimmucc.png", width = 700, height = 600)
seqi.results = readRDS("populations/seqi.llsr_output.rds") %>% as.data.frame 
row.names(seqi.results) = seqi.results$cell_type
seqi.results$cell_type = NULL 

c2 = cor(t(seqi.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n SeqiMmucc \n max:",round(max(c2),digits = 2)),
             column_title = "SeqiMmucc")
draw(h1+h2)
dev.off()


#==============================================================================#
#============================Candiolo vs TIMER=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/timers.png", width = 700, height = 600)
timer.results = readRDS("populations/timer_output.rds") %>% as.data.frame 
row.names(timer.results) = timer.results$cell_type
timer.results$cell_type = NULL 

c2 = cor(t(timer.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Timer \n max:",round(max(c2),digits = 2)),
             column_title = "Timer")
draw(h1+h2)
dev.off()

#==============================================================================#
#============================Candiolo vs Xcell=============================#
#==============================================================================#
png("figs/WGCNAfull_vs/xcell.png", width = 700, height = 600)
xcell.results = readRDS("populations/xcell_output.rds") %>% as.data.frame 
row.names(xcell.results) = xcell.results$cell_type
xcell.results$cell_type = NULL 

c2 = cor(t(xcell.results),interesting.eigengenes,method="spearman")
h2 = Heatmap(t(c2),
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             name = paste("Spearman \n Xcell \n max:",round(max(c2),digits = 2)),
             column_title = "Xcell")
draw(h1+h2)
dev.off()












