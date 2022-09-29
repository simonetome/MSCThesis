setwd("C:/Users/simon/Tesi_tome")
source("rutils/wgcna_wrapper.r")
source("rutils/data_loading.r")

library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(tidyr)

assembly = readRDS("datasets/robjects/full_assembly_gt400.rds")
human = top_variance(assembly[[1]],3000)
murine = top_variance(assembly[[2]],3000)

sample.traits = readRDS("populations/candiolo_tpm_output.rds")
colnames(sample.traits) = c("Endothelial","Leucocyte","CAF")

# Parameters in analysis for murine
rsquared_cut = 0.85
min_clust_size = 10
deepSplit = 2
 
murine.WGCNA.analysis = perform.WGCNA(dataset = t(murine),
                                      title = "figs/murine",
                                      min_size = min_clust_size,
                                      deepsplit = deepSplit,
                                      cortype = "spearman",
                                      type.adj = "signed hybrid",
                                      RsquaredCut = rsquared_cut)
murine.eigengenes = murine.WGCNA.analysis$eigengenes$eigengenes
colnames(murine.eigengenes) = paste("M_",colnames(murine.eigengenes),sep = "")

cor1 = cor(murine.eigengenes, sample.traits, method = "spearman")
colnames(cor1) = c("Endothelium","Leucocyte","CAF")
cor1 = as.data.frame(cor1)
h1 = Heatmap(cor1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.2f", cor1[i, j]), x, y, gp = gpar(fontsize = 5))})
draw(h1)

cor.threshold = 0.6
end.eigengenes = row.names(cor1 %>% filter(Endothelium >= cor.threshold))
leu.eigengenes = row.names(cor1 %>% filter(Leucocyte >= cor.threshold))
caf.eigengenes = row.names(cor1 %>% filter(CAF >= cor.threshold))

end.eigengenes
leu.eigengenes
caf.eigengenes

interesting.eigengenes = c(end.eigengenes,leu.eigengenes,caf.eigengenes)

cor.e = cor(murine.eigengenes[,interesting.eigengenes],sample.traits)


png("figs/interestingEigengenes.png", width = 400, height = 600)
Heatmap(cor.e,
        column_title = "Correlation between sample traits and murine eigengenes",
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        use_raster = FALSE,
        name = paste("Spearman \n max:",round(max(cor.e),digits = 2)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.e[i, j]), x, y, gp = gpar(fontsize = 10))}
          )
dev.off()




#==============================================================================#
# Test int.eigengenes against human eigengenes and look for cor distribution   # 
#==============================================================================#
results = as.data.frame(matrix(nrow = 64, ncol = 9))
colnames(results) = c("Cut","Min_size","Deepsplit","Top5","Top1","Max","MaxEnd","MaxLeu","MaxCaf")
counter = 1

cuts = c(0.8,0.85,0.9,0.95)
min_sizes = c(10,15,20,25)
deepsplits = c(1,2,3,4)


for(c in cuts){
  for(m in min_sizes){
    for(d in deepsplits){
      rsquared_cut = c
      min_clust_size = m
      deepSplit = d
      
      human.analysis = perform.WGCNA(dataset = t(human),
                                      title = "WGCNA_output/human",
                                      min_size = min_clust_size,
                                      deepsplit = deepSplit,
                                      cortype = "spearman",
                                      type.adj = "signed hybrid",
                                      RsquaredCut = rsquared_cut)
      human.eigengenes = human.analysis$eigengenes$eigengenes
      
      cor.matrix = as.matrix(cor(murine.eigengenes[,interesting.eigengenes],human.eigengenes, method = "spearman"))
      cor.caf = as.matrix(cor(murine.eigengenes[,caf.eigengenes],human.eigengenes, method = "spearman"))
      cor.end = as.matrix(cor(murine.eigengenes[,end.eigengenes],human.eigengenes, method = "spearman"))
      cor.leu = as.matrix(cor(murine.eigengenes[,leu.eigengenes],human.eigengenes, method = "spearman"))
      
      top5 = (as.numeric(quantile(cor.matrix, probs = c(0.95))))
      top1 = (as.numeric(quantile(cor.matrix, probs = c(0.99))))
      max = max(cor.matrix)
      
      results[counter,"Cut"] = c
      results[counter,"Min_size"] = m
      results[counter,"Deepsplit"] = d
      results[counter,"Top5"] = top5
      results[counter,"Top1"] = top1
      results[counter,"Max"] = max
      results[counter,"MaxCaf"] = max(cor.caf)
      results[counter,"MaxEnd"] = max(cor.end)
      results[counter,"MaxLeu"] = max(cor.leu)
      
      counter = counter + 1
      
}}}


#==============================================================================#
#                                Plot results                                  #
#==============================================================================#

results[7:9] %>% 
  gather(key="Metric", value="Value") %>% 
  ggplot(aes(x=Metric, y=Value, fill=Metric)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+ 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Eigengene to eigengene correlation \n using interesting murine eigengenes \n varying human WGCNA results") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/eigengene2eigengene.png",
       width = 8,
       height = 6,)


#==============================================================================#
#  What about correlating with single genes?
#==============================================================================#

cor.single.caf = cor(murine.eigengenes[,caf.eigengenes],t(human))
cor.single.leu = cor(murine.eigengenes[,leu.eigengenes],t(human))
cor.single.end = cor(murine.eigengenes[,end.eigengenes],t(human))

# CAF
png("figs/human_genes_to_cafEigenegenes.png", width = 800, height = 400)
Heatmap(cor.single.caf,
        column_title = "CAF: human top 3000 Genes",
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        use_raster = FALSE,
        name = paste("Spearman \n max:",round(max(cor.single.caf),digits = 2)))
dev.off()
# Endothelial
png("figs/human_genes_to_endEigenegenes.png", width = 800, height = 400)
Heatmap(cor.single.end,
        column_title = "Endothelial: human top 3000 Genes",
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        use_raster = FALSE,
        name = paste("Spearman \n max:",round(max(cor.single.end),digits = 2)))
dev.off()
# Leucocyte
png("figs/human_genes_to_leuEigenegenes.png", width = 800, height = 400)
Heatmap(cor.single.leu,
        column_title = "Leucocyte: human top 3000 Genes",
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        use_raster = FALSE,
        name = paste("Spearman \n max:",round(max(cor.single.leu),digits = 2)))
dev.off()




print(interesting.eigengenes)






















































