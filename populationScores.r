library(ggplot2)
library(dplyr)
library(tidyr)

murine_based = readRDS("WGCNAcorrbenchmark/WGCNA_benchmark_murine_based_check.rds") 
human_based = readRDS("WGCNAcorrbenchmark/WGCNA_benchmark_human_based_check.rds")
candiolo = readRDS("WGCNAcorrbenchmark/WGCNA_benchmark_candiolo_murine_check.rds")

epic = (readRDS("WGCNAcorrbenchmark/WGCNA_benchmark_epic.rds"))
epic.no = (readRDS("WGCNAcorrbenchmark/WGCNA_benchmark_epic_no.rds"))

# WGCNA has been performed with 64 different triples of parameters (2500 murine top genes) 
# each point in the boxplot corresponds to the maximum correlation between an eigengene and 
# a sample trait in a WGCNA run 

# this has the aim to show the variance of population score correlation wrt WGCNA parameters


candiolo_wide = candiolo[4:length(colnames(candiolo))]
candiolo_wide %>% select(Endothelial_max,CAF_max,Leucocyte_max) %>%
  gather(key="CellType", value="Value") -> candiolo_long

candiolo_long %>% ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() +
  ggtitle("Candiolo signatures") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/candiolo_check.png",
       width = 8,
       height = 6,)


# MURINE BASED 
#"mmcp_output.rds"      

murine_based$mmcp_output.rds[4:length(colnames(murine_based$mmcp_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("MMCP") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/mmcp_check.png",
       width = 8,
       height = 6,)
  
#"seqi.llsr_output.rds" 
murine_based$seqi.llsr_output.rds[4:length(colnames(murine_based$seqi.llsr_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Seqi") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/seqi_check.png",
       width = 8,
       height = 6,)

#"dcq_output.rds"  
murine_based$dcq_output.rds[4:length(colnames(murine_based$dcq_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("DCQ") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/dcq_check.png",
       width = 8,
       height = 6,)


#"base_output.rds" 
murine_based$base_output.rds[4:length(colnames(murine_based$base_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Base") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/base_check.png",
       width = 8,
       height = 6,)
  
# HUMAN BASED
#"estimate_output.rds"      
human_based$estimate_output.rds[4:length(colnames(human_based$estimate_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Estimate") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/estimate_check.png",
       width = 8,
       height = 6,)


#"consensus_tme_output.rds" 
human_based$consensus_tme_output.rds[4:length(colnames(human_based$consensus_tme_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Consensus tme") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/consensus_check.png",
       width = 8,
       height = 6,)


#"abis_output.rds"       
human_based$abis_output.rds[4:length(colnames(human_based$abis_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Abis") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/abis_check.png",
       width = 8,
       height = 6,)


#"xcell_output.rds"   
human_based$xcell_output.rds[4:length(colnames(human_based$xcell_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Xcell") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/xcell_check.png",
       width = 8,
       height = 6,)


#"mcp_counter_output.rds" 
human_based$mcp_counter_output.rds[4:length(colnames(human_based$mcp_counter_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("MCP counter") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/mcp_check.png",
       width = 8,
       height = 6,)


#"timer_output.rds"         
human_based$timer_output.rds[4:length(colnames(human_based$timer_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Timer") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/timer_check.png",
       width = 8,
       height = 6,)


#"quantiseq_output.rds"  
human_based$quantiseq_output.rds[4:length(colnames(human_based$quantiseq_output.rds))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Quantiseq") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/quantiseq_check.png",
       width = 8,
       height = 6,)



#"epic"  
epic[4:length(colnames(epic))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("Epic") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/epic.png",
       width = 8,
       height = 6,)



#"epic no"  
epic.no[4:length(colnames(epic.no))] %>% 
  gather(key="CellType", value="Value") %>% 
  ggplot(aes(x=CellType, y=Value, fill=CellType)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))+
  ggtitle("epic no Other types") + theme(plot.title = element_text(hjust = 0.5))
ggsave("figs/epicno.png",
       width = 8,
       height = 6,)







