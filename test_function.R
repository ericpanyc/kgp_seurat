library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(AnnotationHub)
library(limma)
library(edgeR)
library(ggplot2)
library(cowplot)
library(Matrix)
library(scales)
library(RCurl)
library(purrr)
library(metap)
library(multtest)
library(xlsx)
library(ggpubr)
library(clustree)
library(DESeq)
library(DESeq2)
library(ComplexHeatmap)
library(pheatmap)
library(reticulate)

my_colour = list(
  cells = c(cluster0="#a47d4c",cluster1="#2962ff",cluster2="#c30101",
            cluster3="#32cd32",cluster4="#7e5aa2",cluster5="#ffd7cd",
            cluster6="#03989e",cluster7="#c8b9ac",cluster8="#d61b5d",
            cluster9="#848c2f",cluster10="#bce3dd",cluster11="#383838",
            cluster12="#a2665c",cluster13="#5b56bc",cluster14="#33ccff",
            cluster15="#cc0066",cluster16="#FFC300",cluster17="#388A5C",
            cluster18="#9157A8",cluster19="#6030B9",cluster20="#CB5F9F"),
  markers = c(cluster0="#a47d4c",cluster1="#2962ff",cluster2="#c30101",
              cluster3="#32cd32",cluster4="#7e5aa2",cluster5="#ffd7cd",
              cluster6="#03989e",cluster7="#c8b9ac",cluster8="#d61b5d",
              cluster9="#848c2f",cluster10="#bce3dd",cluster11="#383838",
              cluster12="#a2665c",cluster13="#5b56bc",cluster14="#33ccff",
              cluster15="#cc0066",cluster16="#FFC300",cluster17="#388A5C",
              cluster18="#9157A8",cluster19="#6030B9",cluster20="#CB5F9F")
)

make_pheatmap <- function(seurat_object) {
  #dense_table <- seurat_object[["SCT"]]@scale.data
  #Idents(object = seurat_object) <- required_reso
  #markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  #markers <- markers[which(markers$gene %in% row.names(dense_table)),]
  pheat_data <- dense_table[markers$gene,]
  row_anno <- markers[,c("gene","cluster")] %>% data.frame()
  rownames(row_anno) <- make.names(row_anno$gene, unique = T)
  rownames(pheat_data) <- rownames(row_anno)
  row_anno <- row_anno[2]
  names(row_anno) <- "markers"
  markers_per_cluster <- lapply(unique(row_anno$markers), function(i) sum(row_anno$markers == i)) %>% flatten_int()
  current_sum <- 0
  row_gap_list <- lapply(c(1:(length(markers_per_cluster)-1)), function(i) current_sum <<- current_sum + markers_per_cluster[[i]]) %>% flatten_dbl()
  col_anno <- metadata[,c(required_reso), drop = F]
  col_anno$cell <- row.names(col_anno)
  col_anno <- col_anno[order(col_anno[,1]),]
  col_anno <- col_anno[1]
  names(col_anno) <- "cells"
  cells_per_cluster <- lapply(unique(col_anno$cells), function(i) sum(col_anno$cells == i)) %>% flatten_int()
  current_sum <- 0
  col_gap_list <- lapply(c(1:(length(cells_per_cluster)-1)), function(i) current_sum <<- current_sum + cells_per_cluster[[i]]) %>% flatten_dbl()
  row_anno$markers <- paste0("cluster", row_anno$markers) %>% as.vector()
  col_anno$cells <- paste0("cluster", col_anno$cells) %>% as.vector()
  pheat_data <- pheat_data[,rownames(col_anno)]
  pheatmap(pheat_data, annotation_row = row_anno, annotation_col = col_anno, 
           show_colnames = F, cluster_rows = F, cluster_cols = F,
           annotation_colors = my_colour, annotation_names_row = F,
           annotation_legend = T, gaps_row = row_gap_list,
           main = paste0("Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))
}

setwd("/Users/ericpan/Documents/pan/usc/kgp/seurat")
use_python("/usr/local/bin/python3")

filtered_seurat <- readRDS("data/VS4_SCT_filtered.RDS")

required_reso <- "SCT_snn_res.0.4"
Idents(object = filtered_seurat) <- required_reso

dense_table <- filtered_seurat[["SCT"]]@scale.data
metadata <- filtered_seurat@meta.data
markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
markers <- markers[which(markers$gene %in% row.names(dense_table)),]

make_pheatmap(filtered_seurat)


# dense_table <- as.matrix(GetAssayData(object = filtered_seurat, slot = "scale.data"))
# dense_table <- as.data.frame(dense_table) %>% rownames_to_column(var = "gene")
# make_pheatmap <- function(seurat_object) {
#   dense_table <- seurat_object[["SCT"]]@scale.data
#   Idents(object = seurat_object) <- required_reso
#   markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#   markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#   markers <- markers[which(markers$gene %in% row.names(dense_table)),]
#   pheat_data <- dense_table[markers$gene,]
#   row_anno <- markers[,c("gene","cluster")] %>% data.frame()
#   rownames(row_anno) <- make.names(row_anno$gene, unique = T)
#   rownames(pheat_data) <- rownames(row_anno)
#   row_anno <- row_anno[2]
#   names(row_anno) <- "markers"
#   markers_per_cluster <- lapply(unique(row_anno$markers), function(i) sum(row_anno$markers == i)) %>% flatten_int()
#   current_sum <- 0
#   row_gap_list <- lapply(c(1:(length(markers_per_cluster)-1)), function(i) current_sum <<- current_sum + markers_per_cluster[[i]]) %>% flatten_dbl()
#   col_anno <- metadata[,c(required_reso), drop = F]
#   col_anno$cell <- row.names(col_anno)
#   col_anno <- col_anno[order(col_anno[,1]),]
#   col_anno <- col_anno[1]
#   names(col_anno) <- "cells"
#   cells_per_cluster <- lapply(unique(col_anno$cells), function(i) sum(col_anno$cells == i)) %>% flatten_int()
#   current_sum <- 0
#   col_gap_list <- lapply(c(1:(length(cells_per_cluster)-1)), function(i) current_sum <<- current_sum + cells_per_cluster[[i]]) %>% flatten_dbl()
#   row_anno$markers <- paste0("cluster", row_anno$markers) %>% as.vector()
#   col_anno$cells <- paste0("cluster", col_anno$cells) %>% as.vector()
#   pheat_data <- pheat_data[,rownames(col_anno)]
#   pheatmap(pheat_data, annotation_row = row_anno, annotation_col = col_anno, 
#            show_colnames = F, cluster_rows = F, cluster_cols = F,
#            annotation_colors = my_colour, annotation_names_row = F,
#            annotation_legend = T, gaps_row = row_gap_list,
#            main = paste0("Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))
# }
