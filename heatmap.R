library(ggplot2)
library(dplyr)
library(plotly)
library(ggfortify)
library(cummeRbund)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(factoextra)
library(xlsx)

setwd("~/Downloads/usc_folder/course/2020Spring/TRGN515/project")


#Load raw data
DE_gene_raw <- read.csv("gene_exp.diff", sep = "\t")
FPKM_gene_raw <- read.csv("genes.fpkm_tracking", sep = "\t")
#Drop lowly expressed genes
DE_gene <- DE_gene_raw[!(DE_gene_raw$log2.fold_change. == Inf | DE_gene_raw$log2.fold_change. == -Inf) ,]

#Find up-regulated genes
up_reg_filter <- which(DE_gene$p_value <= 0.05)
up_reg_filter <- up_reg_filter[which(DE_gene[up_reg_filter,]$log2.fold_change. >= 1.0)]
up_reg_genes <- DE_gene[up_reg_filter,]

#Find down-regulated genes
down_reg_filter <- which(DE_gene$p_value <= 0.05)
down_reg_filter <- down_reg_filter[which(DE_gene[down_reg_filter,]$log2.fold_change. <= -3)]
down_reg_genes <- DE_gene[down_reg_filter,]

#Find differentially expressed genes
sig_DE_filter <- which(DE_gene$p_value <= 0.01)
sig_DE_filter <- sig_DE_filter[which(DE_gene[sig_DE_filter,]$log2.fold_change. >= 2 | DE_gene[sig_DE_filter,]$log2.fold_change. <= -3 )]
sig_DE_genes <- DE_gene[sig_DE_filter,]
sig_genes_list <- as.character(sig_DE_genes$gene_id)

#Find DE genes in the FPKM table
FPKM_gene_sig <- FPKM_gene_raw[which(FPKM_gene_raw$gene_id %in% sig_genes_list),]
FPKM_gene_sig <- select(FPKM_gene_sig,"gene_short_name",
                           "AD4_1_FPKM","AD4_2_FPKM","AD4_3_FPKM",
                           "c_AD4_1_FPKM","c_AD4_2_FPKM","c_AD4_3_FPKM",
                           "AD5_1_FPKM","AD5_2_FPKM","AD5_3_FPKM",
                           "c_AD5_1_FPKM","c_AD5_2_FPKM","c_AD5_3_FPKM")
#Build data matrix to make a heatmap
rownames(FPKM_gene_sig) <- FPKM_gene_sig$gene_short_name
FPKM_gene_sig <- FPKM_gene_sig[,-c(1)]
colnames(FPKM_gene_sig) <- sub("_FPKM", "", colnames(FPKM_gene_sig))
FPKM_gene_sig_mat <- data.matrix(t(FPKM_gene_sig))
FPKM_gene_sig <- as.data.frame(t(FPKM_gene_sig))
#Plot heatmap
draw_colnames_45 <- function (coln, gaps, ...) {
        coord <- pheatmap:::find_coordinates(length(coln), gaps)
        x     <- coord$coord - 0.5 * coord$size
        res   <- grid::textGrob(
                coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
                vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
        )
        return(res)
}
assignInNamespace(
        x = "draw_colnames",
        value = "draw_colnames_45",
        ns = asNamespace("pheatmap")
)
pheatmap(FPKM_gene_sig_mat, scale = "column", drop_levels = TRUE,fontsize = 12, fontsize_col = 6, fontsize_row = 10)

#Build full table with annotation
sample_name <- row.names(FPKM_gene_sig_mat)
full_table <- as.data.frame(FPKM_gene_sig_mat)
full_table$sample <- sample_name
full_table$group <- c("disease","disease","disease",
                      "control","control","control",
                      "disease","disease","disease",
                      "control","control","control")
#Plot PCA
autoplot(prcomp(FPKM_gene_sig_mat), data = full_table, colour = "sample", shape = "group", size = 5)
p <- autoplot(prcomp(FPKM_gene_sig_mat), data = full_table, colour = "sample", shape = "group")
ggplotly(p = p)

#Plot dendrogram
scaled_FPKM_table <- scale(FPKM_gene_sig)
dist_cor <- get_dist(scaled_FPKM_table, method = "pearson")
hclust_ward_cor <- hclust(d = dist_cor, method = "ward.D2")
fviz_dend(hclust_ward_cor, k = 4, labels = FALSE,
          k_colors = c("#FF0000", "#2ECC71", "#8E44AD", "#2E86C1"),
          color_labels_by_k = TRUE,
          rect = TRUE,
          main = "pearson_ward_dendrogram", ylab = "")

#IPA analysis gene list
up_reg_filter_IPA <- which(DE_gene$p_value <= 0.05)
up_reg_filter_IPA <- up_reg_filter_IPA[which(DE_gene[up_reg_filter_IPA,]$log2.fold_change. >= 1.5)]
up_reg_genes_IPA <- DE_gene[up_reg_filter_IPA,]
up_genes_table <- select(up_reg_genes_IPA, "gene_id", "log2.fold_change.", "p_value", "q_value")
write.xlsx(up_genes_table, "up_genes_00515.xlsx", row.names = FALSE)


down_reg_filter_IPA <- which(DE_gene$p_value <= 0.05)
down_reg_filter_IPA <- down_reg_filter_IPA[which(DE_gene[down_reg_filter_IPA,]$log2.fold_change. <= -1.5)]
down_reg_genes_IPA <- DE_gene[down_reg_filter_IPA,]
down_genes_table <- select(down_reg_genes_IPA, "gene_id", "log2.fold_change.", "p_value", "q_value")
write.xlsx(down_genes_table, "down_genes_more.xlsx")

#log10 convertion and plot
processed_FPKM_mat <- log10(FPKM_gene_sig_mat+0.01)

autoplot(prcomp(processed_FPKM_mat), data = full_table, colour = "sample", shape = "group", size = 5)
pheatmap(processed_FPKM_mat, scale = "column", drop_levels = TRUE,fontsize = 12, fontsize_col = 6, fontsize_row = 10)

#test new gene exp file
gene_diff_ad4 <- read.csv("AD4_gene_exp.diff", sep = "\t")
gene_diff_ad4_cp <- gene_diff_ad4[!(gene_diff_ad4$log2.fold_change. == Inf | gene_diff_ad4$log2.fold_change. == -Inf) ,]

up_reg_filter_ad4 <- which(gene_diff_ad4_cp$p_value <= 0.05)
up_reg_filter_ad4 <- up_reg_filter_ad4[which(gene_diff_ad4_cp[up_reg_filter_ad4,]$log2.fold_change. >= 2)]
up_reg_genes_ad4 <- gene_diff_ad4_cp[up_reg_filter_ad4,]
up_genes_table_new <- select(up_reg_genes_ad4, "gene_id", "log2.fold_change.", "p_value", "q_value")
write.xlsx(up_genes_table_new, "up_genes_new.xlsx")

down_reg_filter_ad4 <- which(gene_diff_ad4_cp$p_value <= 0.05)
down_reg_filter_ad4 <- down_reg_filter_ad4[which(gene_diff_ad4_cp[down_reg_filter_ad4,]$log2.fold_change. <= -3)]
down_reg_genes_ad4 <- gene_diff_ad4_cp[down_reg_filter_ad4,]
down_genes_table_new <- select(down_reg_genes_ad4, "gene_id", "log2.fold_change.", "p_value", "q_value")
write.xlsx(down_genes_table_new, "down_genes_new.xlsx", row.names = FALSE)

sig_genes_ad4 <- rbind(up_reg_genes_ad4, down_reg_genes_ad4)


