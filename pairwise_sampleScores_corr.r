setwd("C:/Users/simon/Tesi_tome")
source("rutils/data_loading.r")

library(dplyr)
library(ggplot2)
library(plotly)
library(ComplexHeatmap)

#==============================================================================#
candiolo = readRDS("populations/candiolo_tpm_output.rds") %>% as.data.frame()
colnames(candiolo) = c("Endothelial","Leucocyte","CAF")
#==============================================================================#
mmcp =  readRDS("populations/mmcp_output.rds") %>% as.data.frame()
dcq =  readRDS("populations/dcq_output.rds") %>% as.data.frame()
base =  readRDS("populations/base_output.rds") %>% as.data.frame()
seqi =  readRDS("populations/seqi.llsr_output.rds") %>% as.data.frame()
#==============================================================================#
abis =  readRDS("populations/abis_output.rds") %>% as.data.frame()
xcell =  readRDS("populations/xcell_output.rds") %>% as.data.frame()
estimate =  readRDS("populations/estimate_output.rds") %>% as.data.frame()
mcp =  readRDS("populations/mcp_counter_output.rds") %>% as.data.frame()
quantiseq =  readRDS("populations/quantiseq_output.rds") %>% as.data.frame()
timer =  readRDS("populations/timer_output.rds") %>% as.data.frame()
consensus_tme =  readRDS("populations/consensus_tme_output.rds") %>% as.data.frame()
epic = (readRDS("populations/EPIC_output.rds"))$cellFractions %>% as.data.frame()
epic.no = (readRDS("populations/EPIC_output_no.rds"))$cellFractions %>% as.data.frame()

#==============================================================================#

row.names(mmcp) = mmcp$cell_type
row.names(dcq) = dcq$cell_type
row.names(base) = base$cell_type
row.names(seqi) = seqi$cell_type

row.names(abis) = abis$cell_type
row.names(xcell) = xcell$cell_type
row.names(estimate) = estimate$cell_type
row.names(mcp) = mcp$cell_type
row.names(quantiseq) = quantiseq$cell_type
row.names(timer) = timer$cell_type
row.names(consensus_tme) = consensus_tme$cell_type

#==============================================================================#

mmcp$cell_type = NULL
dcq$cell_type = NULL
base$cell_type = NULL
seqi$cell_type = NULL

abis$cell_type = NULL
xcell$cell_type = NULL
estimate$cell_type = NULL
mcp$cell_type = NULL
quantiseq$cell_type = NULL
timer$cell_type = NULL
consensus_tme$cell_type = NULL

#==============================================================================#

candiolo %>% zscore() %>% as.data.frame() -> candiolo.zs
epic %>% zscore() %>% as.data.frame() -> epic.zs
epic.no %>% zscore() %>% as.data.frame() -> epic.no.zs

mmcp %>% zscore_col() %>% as.data.frame() %>% t() -> mmcp.zs
dcq %>% zscore_col() %>% as.data.frame()  %>% t()-> dcq.zs
base %>% zscore_col() %>% as.data.frame()  %>% t()-> base.zs
seqi %>% zscore_col() %>% as.data.frame()  %>% t()-> seqi.zs

abis %>% zscore_col() %>% as.data.frame()  %>% t()-> abis.zs
xcell %>% zscore_col() %>% as.data.frame()  %>% t()-> xcell.zs
estimate %>% zscore_col() %>% as.data.frame()  %>% t()-> estimate.zs
mcp %>% zscore_col() %>% as.data.frame()  %>% t()-> mcp.zs
quantiseq %>% zscore_col() %>% as.data.frame()  %>% t()-> quantiseq.zs
timer %>% zscore_col() %>% as.data.frame()  %>% t()-> timer.zs
consensus_tme %>% zscore_col() %>% as.data.frame()  %>% t()-> consensus_tme.zs


#==============================================================================#
# zsvcore is applied sample-wise: each sample will have mean 0 and std 1 

# rows are samples
sum(rowVars(candiolo.zs %>% as.matrix)) == 624
sum(rowVars(consensus_tme.zs %>% as.matrix)) == 624

cor.candiolo = cor(candiolo.zs,candiolo.zs, method = "spearman")

cor.epic = cor(epic.zs,epic.zs, method = "spearman")
cor.epic.no = cor(epic.no.zs,epic.no.zs, method = "spearman")

cor.mmcp = cor(mmcp.zs,mmcp.zs, method = "spearman")
cor.dcq = cor(dcq.zs,dcq.zs, method = "spearman")
cor.seqi = cor(seqi.zs,seqi.zs, method = "spearman")
cor.base = cor(base.zs,base.zs, method = "spearman")

cor.abis = cor(abis.zs,abis.zs, method = "spearman")
cor.xcell = cor(xcell.zs,xcell.zs, method = "spearman")
cor.estimate = cor(estimate.zs,estimate.zs, method = "spearman")
cor.mcp = cor(mcp.zs,mcp.zs, method = "spearman")
cor.quantiseq = cor(quantiseq.zs,quantiseq.zs, method = "spearman")
cor.timer = cor(timer.zs,timer.zs, method = "spearman")
cor.consensus = cor(consensus_tme.zs,consensus_tme.zs, method = "spearman")


#==============================================================================#

png("figs/cor.candiolo.png", width = 700, height = 600)
Heatmap(cor.candiolo,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - Candiolo",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.candiolo[i, j]), x, y, gp = gpar(fontsize = 10))}
        )
dev.off()
        
#==============================================================================#
#==============================================================================#

png("figs/cor.epic.png", width = 700, height = 600)
Heatmap(cor.epic,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - EPIC all",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.epic[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()


png("figs/cor.epic.no.png", width = 700, height = 600)
Heatmap(cor.epic.no,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - EPIC no Other",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.epic.no[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

#==============================================================================#

png("figs/cor.mmcp.png", width = 700, height = 600)
Heatmap(cor.mmcp,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - MMCP",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.mmcp[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.dcq.png", width = 700, height = 600)
Heatmap(cor.dcq,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - DCQ",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.dcq[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.base.png", width = 700, height = 600)
Heatmap(cor.base,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - Base",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.base[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.seqi.png", width = 700, height = 600)
Heatmap(cor.seqi,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - SeqiMmucc",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.seqi[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

#==============================================================================#

png("figs/cor.abis.png", width = 700, height = 600)
Heatmap(cor.abis,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - Abis",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.abis[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.xcell.png", width = 1400, height = 1200)
Heatmap(cor.xcell,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - XCell",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.xcell[i, j]), x, y, gp = gpar(fontsize = 12))}
)
dev.off()

png("figs/cor.estimate.png", width = 700, height = 600)
Heatmap(cor.estimate,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - Estimate",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.estimate[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.mcp.png", width = 700, height = 600)
Heatmap(cor.mcp,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - MCP",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.mcp[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.quantiseq.png", width = 700, height = 600)
Heatmap(cor.quantiseq,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - Quantiseq",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.quantiseq[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.timer.png", width = 700, height = 600)
Heatmap(cor.timer,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - Timer",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.timer[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

png("figs/cor.consensus.png", width = 700, height = 600)
Heatmap(cor.consensus,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = "Spearman",
        column_title = "Correlation value - ConsensusTME",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", cor.consensus[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()











