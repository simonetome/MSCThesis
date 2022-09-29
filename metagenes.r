setwd("C:/Users/simon/Tesi_tome")

library("ComplexHeatmap")
library("lsa")
library("circlize")
library("stylo")
library("dendsort")
library("dplyr")

# read the csv and put geneid as row names 
m_zscore <- read.csv("datasets/m_zscored_top1000.csv")
h_zscore <- read.csv("datasets/h_zscored_top1000.csv")
samples <- tail(colnames(m_zscore),-1)
rownames(m_zscore) <- m_zscore$Geneid
rownames(h_zscore) <- h_zscore$Geneid
m_zscore <- m_zscore[samples]
h_zscore <- h_zscore[samples]

# function object to hclust with cosine distance 
h_c <- function(x) hclust(dist.cosine(as.matrix(x)))

# dendrograms object for mzscore and hzscore
row_dend_m = dendsort(hclust(dist.cosine(as.matrix(m_zscore))))
col_dend_m = dendsort(hclust(dist.cosine(t(as.matrix(m_zscore)))))
row_dend_h = dendsort(hclust(dist.cosine(as.matrix(h_zscore))))
col_dend_h = dendsort(hclust(dist.cosine(t(as.matrix(h_zscore)))))

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


# MURINE

pam_m = cluster::pam(m_zscore, k = 29)

#png("RObjects/heatmap1.png",width=7.5,height=4,units="in",res=1200)
ht_opt$TITLE_PADDING = unit(c(7.5, 7.5), "points")
Heatmap(as.matrix(m_zscore),
        column_title = "Murine zscore top1000 pc genes",
        col = col_fun,
        heatmap_legend_param = list(
          title = "zscore"
        ),
        #cluster_rows = h_c,
        cluster_columns = col_dend_m,
        show_column_dend = FALSE,
        show_row_dend = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_gap = unit(.5,"mm"),
        column_title_gp = gpar(fill = "red", col = "white", border = "blue"),
        row_split = paste0("", pam_m$clustering)
        )
dev.off()
m_zscore$Pam <- pam_m$clustering

meta_genes_m <- data.frame(matrix(ncol = length(colnames(m_zscore)), nrow = 14))
colnames(meta_genes_m) <- colnames(m_zscore)
meta_genes_m$Pam <- seq(1,14)

for (i in seq(1,29)){
  meta_genes_m[meta_genes_m$Pam == i,] <- colMeans(m_zscore[m_zscore$Pam==i,])
}

#png("RObjects/heatmap2.png",width=7.5,height=4,units="in",res=1200)
# ht_opt$TITLE_PADDING = unit(c(7.5, 7.5), "points")
# Heatmap(meta_genes_m[samples],
#         column_title = "Murine metagenes",
#         show_column_names = FALSE,
#         heatmap_legend_param = list(
#           title = "avg zscore"
#         ),
#         cluster_rows = FALSE,
#         cluster_columns = col_dend_m,)

#dev.off()

# HUMAN

pam_h = cluster::pam(h_zscore, k = 15)

png("RObjects/heatmap3_15.png",width=7.5,height=4,units="in",res=1200)
ht_opt$TITLE_PADDING = unit(c(7.5, 7.5), "points")
Heatmap(as.matrix(h_zscore),
        column_title = "Human zscore top1000 pc genes",
        col = col_fun,
        heatmap_legend_param = list(
          title = "zscore"
        ),
        #cluster_rows = h_c,
        cluster_columns = col_dend_h,
        show_column_dend = FALSE,
        show_row_dend = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_gap = unit(.5,"mm"),
        column_title_gp = gpar(fill = "red", col = "white", border = "blue"),
        row_split = paste0("", pam_h$clustering)
)
dev.off()

h_zscore$Pam <- pam_h$clustering

meta_genes_h <- data.frame(matrix(ncol = length(colnames(h_zscore)), nrow = 9))
colnames(meta_genes_h) <- colnames(h_zscore)
meta_genes_h$Pam <- seq(1,9)

for (i in seq(1,9)){
  meta_genes_h[meta_genes_h$Pam == i,] <- colMeans(h_zscore[h_zscore$Pam==i,])
}

png("RObjects/heatmap4.png",width=7.5,height=4,units="in",res=1200)
ht_opt$TITLE_PADDING = unit(c(7.5, 7.5), "points")
Heatmap(meta_genes_h[samples],
        show_column_names = FALSE,
        column_title = "Human metagenes",
        heatmap_legend_param = list(
          title = "avg zscore"
        ),
        cluster_rows = FALSE,
        cluster_columns = col_dend_h,)
dev.off()

# correlation between metagenes 
num_human_metagenes = 9
num_murine_metagenes = 14
cor_bet_metagenes <- matrix(nrow = 14, ncol = 9)

for(i in 1:nrow(cor_bet_metagenes)){
  for(j in 1:ncol(cor_bet_metagenes)){
    cor_bet_metagenes[i,j] <- cor(as.numeric(meta_genes_m[i,samples]),
                                  as.numeric(meta_genes_h[j,samples]),
                                  ,method = "pearson")
  }
}

# correlation between murine metagenes vs top1000 human genes 

cor_bet_genes_1 <- matrix(nrow = num_murine_metagenes, ncol = 1000)
colnames(cor_bet_genes_1) <- row.names(h_zscore)

for(i in 1:nrow(cor_bet_genes_1)){
  for(j in 1:ncol(cor_bet_genes_1)){
    cor_bet_genes_1[i,j] <- cor(as.numeric(meta_genes_m[i,samples]),
                                as.numeric(h_zscore[j,samples]),
                                method = "pearson")
  }
}

# correlation between human metagenes vs top1000 murine genes 

cor_bet_genes_2 <- matrix(nrow = num_human_metagenes, ncol = 1000)
colnames(cor_bet_genes_2) <- row.names(m_zscore)

for(i in 1:nrow(cor_bet_genes_2)){
  for(j in 1:ncol(cor_bet_genes_2)){
    cor_bet_genes_2[i,j] <- cor(as.numeric(meta_genes_h[i,samples]),
                                as.numeric(m_zscore[j,samples]),
                                method = "pearson")
  }
}

# correlation between all top1000 genes 

# cor_bet_genes <- matrix(nrow = 1000, ncol = 1000)
# row.names(cor_bet_genes) <- row.names(m_zscore)
# colnames(cor_bet_genes) <- row.names(h_zscore)
# 
# for(i in 1:nrow(cor_bet_genes)){
#   for(j in 1:ncol(cor_bet_genes)){
#     cor_bet_genes[i,j] <- cor(as.numeric(m_zscore[i,samples]),
#                               as.numeric(h_zscore[j,samples]),
#                               method = "pearson")
#   }
# }

?saveRDS
saveRDS(cor_bet_metagenes, file = "RObjects/cor_metagenes.rds")
saveRDS(cor_bet_genes_2, file = "RObjects/cor_h_metagenes.rds")
saveRDS(cor_bet_genes_1, file = "RObjects/cor_m_metagenes.rds")

#saveRDS(cor_bet_genes, file = "RObjects/cor_genes.rds")

library("sm")
?sm.density

png("RObjects/cor1.png",width=7.5,height=4,units="in",res=1200)
d <- sm.density(as.vector(cor_bet_metagenes),
                xlab = "Correlation value",
                ylab = "Density")
plot(d)
title("Correlation between metagenes")
dev.off()

png("RObjects/cor2.png",width=7.5,height=4,units="in",res=1200)
d <- sm.density(x = as.vector(cor_bet_genes_1),
                xlab = "Correlation value",
                ylab = "Density")
plot(d)
title("Correlation murine metagenes - top1000 human")
dev.off()

png("RObjects/cor3.png",width=7.5,height=4,units="in",res=1200)
d <- sm.density(x = as.vector(cor_bet_genes_2),
                xlab = "Correlation value",
                ylab = "Density")
title("Correlation human metagenes - top1000 murine")
plot(d)
dev.off()

