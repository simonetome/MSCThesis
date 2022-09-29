setwd("C:/Users/simon/Tesi_tome")
source("rutils/data_loading.r")

library(dplyr)
library(matrixStats)

assembly = load_raw("datasets/filtered_assembly.csv", top.num = NULL)
assembly.top5000 = load_raw("datasets/filtered_assembly.csv", top.num = 5000)

human = as.data.frame(assembly[[1]])
human.top = as.data.frame(assembly.top5000[[1]])

sum(as.numeric(row.names(human.top) %in% row.names(human)))

gene.var = rowVars(as.matrix(human))
top15Var = as.numeric(quantile(gene.var, probs = c(0.85)))

top15.names = row.names(human)[gene.var >= top15Var]
top15.names %in% row.names(top5000.h)


l1 = c("b","a")
l2 = c("a","c","d","b")
l1 %in% l2





