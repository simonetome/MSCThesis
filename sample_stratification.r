setwd("C:/Users/simon/Tesi_tome")

source("wgcna_wrapper.r")
library(cluster)
library(stylo) 
library(ComplexHeatmap)
library(dplyr)
library(psych)

human = readRDS("Datasets/human_zscore_5000.rds")
murine = readRDS("Datasets/murine_zscore_2500.rds")

rsquared_cut = 0.85
min_clust_size = 20
deepSplit = 2

sample.traits = readRDS("populations/candiolo_tpm_output.rds")

murine.WGCNA.analysis = perform.WGCNA(dataset = t(murine),
                                      title = "WGCNA_output/murine",
                                      min_size = min_clust_size,
                                      deepsplit = deepSplit,
                                      cortype = "spearman",
                                      type.adj = "signed hybrid",
                                      RsquaredCut = rsquared_cut)
murine.eigengenes = murine.WGCNA.analysis$eigengenes$eigengenes
colnames(murine.eigengenes) = paste("M_",colnames(murine.eigengenes),sep = "")

cor1 = cor(murine.eigengenes, sample.traits, method = "spearman")
h1 = Heatmap(cor1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.2f", cor1[i, j]), x, y, gp = gpar(fontsize = 5))})

draw(h1)

colnames(cor1) = c("Endothelium","Leucocyte","CAF")
cor1 = as.data.frame(cor1)

# Interesting modules 

end.eigengenes = row.names(cor1 %>% filter(Endothelium >= 0.7))
leu.eigengenes = row.names(cor1 %>% filter(Leucocyte >= 0.7))
caf.eigengenes = row.names(cor1 %>% filter(CAF >= 0.7))
interesting.eigengenes = c(end.eigengenes,leu.eigengenes,caf.eigengenes)

cor = cor(murine.eigengenes[interesting.eigengenes],t(human))
Heatmap(cor)

k = 6
d_h = dist.cosine(t(human))
patient_stratification = pam(d_h, k = k)

stratification = vector(mode = "list", length = k)
for(i in c(1:k)){
  stratification[[i]] = names(patient_stratification$clustering[patient_stratification$clustering == i])
}




?points 
View(cor$r)
View(cor$p)

View(corr.test(temp1,temp2,method = "spearman"))

?text
?corr.test
