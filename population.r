setwd("C:/Users/simon/Tesi_tome")
source("wgcna_wrapper.r")
source("rutils/data_loading.r")
library("stringr")
library("stylo")
library(fpc)
library(cluster)
library(amap)
library(magrittr)
library(dplyr)
library(ComplexHeatmap)

human = readRDS("Datasets/human_zscore_1000.rds")
murine = readRDS("Datasets/murine_zscore_1000.rds")

human.WGCNA.analysis = perform.WGCNA(t(human),"WGCNA_output/human",10,4,"spearman","signed hybrid",0.95)
murine.WGCNA.analysis = perform.WGCNA(t(murine),"WGCNA_output/murine",10,4,"spearman","signed hybrid",0.95)

sample.traits = as.data.frame(readRDS("populations/candiolo_tpm_output.rds"))
#row.names(sample.traits) = sample.traits$cell_type
#sample.traits$cell_type = NULL
#sample.traits = t(sample.traits)

murine.eigengenes = murine.WGCNA.analysis$eigengenes$eigengenes
human.eigengenes = human.WGCNA.analysis$eigengenes$eigengenes
identical(row.names(murine.WGCNA.analysis$eigengenes$eigengenes),row.names(sample.traits))
colnames(murine.eigengenes) = paste("M_",colnames(murine.eigengenes),sep = "")
colnames(human.eigengenes) = paste("H_",colnames(human.eigengenes),sep = "")
cor1 = cor(murine.eigengenes, sample.traits, method = "spearman")
cor2 = cor(murine.eigengenes, human.eigengenes, method = "spearman")

h1 = Heatmap(cor1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.1f", cor1[i, j]), x, y, gp = gpar(fontsize = 5))
             })
h2 = Heatmap(cor2,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.1f", cor2[i, j]), x, y, gp = gpar(fontsize = 5))
             })
draw(h1+h2)


# correlation between traits and murine ME 
cor1 = cor(murine.eigengenes, sample.traits, method = "spearman")
cor2 = cor(murine.eigengenes, human.eigengenes, method = "spearman")

cor_threshold = 0.4
as.data.frame(cor1 > cor_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp1
intersection = cor2[row.names(temp1),]
as.data.frame(intersection > cor_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp2
as.data.frame(cor1 > cor_threshold)[row.names(temp2),]


x = murine.eigengenes$M_MEsalmon
y = human.eigengenes$H_MEblue
par(mar = c(1, 1, 1, 1))
plot(x,y, pch = 19)



cor_neg_threshold = -0.5
as.data.frame(cor1 < cor_neg_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp1
intersection = cor2[row.names(temp1),]
as.data.frame(intersection < cor_neg_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp2
as.data.frame(cor1 < cor_neg_threshold)[row.names(temp2),]


#==============================================================================#

#==============================================================================#
samples = colnames(murine)
PRX = samples[grepl("PRX",samples,fixed = TRUE)]
LMX = samples[grepl("LMX",samples,fixed = TRUE)]

human.lmx = human[LMX]
human.prx = human[PRX]

murine.lmx = murine[LMX]
murine.prx = murine[PRX]

# function(dataset,title,min_size,deepsplit,cortype,type.adj,RsquaredCut){

human.lmx.WGCNA.analysis = perform.WGCNA(dataset = t(human.lmx),
                                         title = "WGCNA_output/human LMX",
                                         min_size = 20,
                                         deepsplit = 2,
                                         cortype = "spearman",
                                         type.adj = "signed hybrid",
                                         RsquaredCut = 0.9)

murine.lmx.WGCNA.analysis =  perform.WGCNA(dataset = t(murine.lmx),
                                           title = "WGCNA_output/murine LMX",
                                           min_size = 20,
                                           deepsplit = 2,
                                           cortype = "spearman",
                                           type.adj = "signed hybrid",
                                           RsquaredCut = 0.9)

murine.lmx.ME = murine.lmx.WGCNA.analysis$eigengenes$eigengenes
human.lmx.ME = human.lmx.WGCNA.analysis$eigengenes$eigengenes

sample.traits.lmx = sample.traits[LMX,]


cor1 = cor(murine.lmx.ME, sample.traits.lmx, method = "pearson")
cor2 = cor(murine.lmx.ME, human.lmx.ME, method = "pearson")


h1 = Heatmap(cor1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.1f", cor1[i, j]), x, y, gp = gpar(fontsize = 4))
             })
h2 = Heatmap(cor2,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.1f", cor2[i, j]), x, y, gp = gpar(fontsize = 4))
             })
draw(h1+h2)


cor_threshold = 0.4
as.data.frame(cor1 > cor_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp1
intersection = cor2[row.names(temp1),]
as.data.frame(intersection > cor_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp2
as.data.frame(cor1 > cor_threshold)[row.names(temp2),]

?EPIC
#==============================================================================#

#==============================================================================#



human.prx.WGCNA.analysis = perform.WGCNA(dataset = t(human.prx),
                                         title = "WGCNA_output/human PRX",
                                         min_size = 20,
                                         deepsplit = 3,
                                         cortype = "spearman",
                                         type.adj = "signed hybrid",
                                         RsquaredCut = 0.9)

murine.prx.WGCNA.analysis =  perform.WGCNA(dataset = t(murine.prx),
                                           title = "WGCNA_output/murine PMX",
                                           min_size = 20,
                                           deepsplit = 3,
                                           cortype = "spearman",
                                           type.adj = "signed hybrid",
                                           RsquaredCut = 0.9)

murine.prx.ME = murine.prx.WGCNA.analysis$eigengenes$eigengenes
human.prx.ME = human.prx.WGCNA.analysis$eigengenes$eigengenes

sample.traits.prx = sample.traits[PRX,]

h1 = Heatmap((cor(murine.prx.ME, sample.traits.prx, method = "pearson")))
h2 = Heatmap((cor(murine.prx.ME, human.prx.ME, method = "pearson")))
draw(h1+h2)

cor1 = cor(murine.prx.ME, sample.traits.prx, method = "pearson")
cor2 = cor(murine.prx.ME, human.prx.ME, method = "pearson")

cor_threshold = 0.5
as.data.frame(cor1 > cor_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp1
intersection = cor2[row.names(temp1),]
as.data.frame(intersection > cor_threshold) %>% filter_all(any_vars(. %in% c(TRUE))) -> temp2
as.data.frame(cor1 > cor_threshold)[row.names(temp2),]

















          