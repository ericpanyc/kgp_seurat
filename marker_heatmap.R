## Marker Heatmap
library(pheatmap)

#Suppose we have a seurat object named "filtered_seurat", which has gone through the SCTranform and UMAP clustering 

# This will extract the raw counts for all genes (not used in the following plot)
non_scale_table <- filtered_seurat[["SCT"]]@data %>% as.matrix()

# This will extract the 3000 variable features identified by SCT
dense_table <- filtered_seurat[["SCT"]]@scale.data
dense_table <- seurat_integrated[["RNA"]]@data %>% as.matrix()


# Set clustering resolution to 1
required_reso <- "SCT_snn_res.1"
Idents(object = filtered_seurat) <- required_reso
Idents(object = seurat_integrated) <- required_reso


# Find top 20 marker genes for each cluster and exclude genes not present in the dense_table matrix
VS4.markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VS4.markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

VS4_top.markers <- VS4.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
VS4_top.markers <- VS4_top.markers[which(VS4_top.markers$gene %in% row.names(dense_table)),]

# Extract expressing data of marker genes only from dense_table, this table will be used for the final heatmap plot
pheat_data <- dense_table[VS4_top.markers$gene,]

# Build annotation table for row (separate marker genes of different clusters)
row_anno <- VS4_top.markers[,c("gene","cluster")] %>% data.frame()
rownames(row_anno) <- make.names(row_anno$gene, unique = T)
rownames(pheat_data) <- rownames(row_anno)
row_anno <- row_anno[2]
names(row_anno) <- "markers"

# Set row gap list (get the number of marker genes in each cluster (accumulative))
markers_per_cluster <- sapply(unique(row_anno$markers), function(i) sum(row_anno$markers == i))
current_sum <- 0
row_gap_list <- sapply(c(1:(length(markers_per_cluster)-1)), function(i) current_sum <<- current_sum + markers_per_cluster[[i]])

# Build annotation table for column (separate cells from different clusters)
metadata <- filtered_seurat@meta.data
col_anno <- metadata[,c(required_reso), drop = F]
col_anno <- col_anno %>% arrange(SCT_snn_res.0.4)
names(col_anno) <- "cells"

# Set column gap list (get the number of cells in each cluster (accumulative))
cells_per_cluster <- sapply(unique(col_anno$cells), function(i) sum(col_anno$cells == i))
current_sum <- 0
col_gap_list <- sapply(c(1:(length(cells_per_cluster)-1)), function(i) current_sum <<- current_sum + cells_per_cluster[[i]])

# Change the legend names displayed on the heatmap 
row_anno$markers <- paste0("cluster", row_anno$markers)
col_anno$cells <- paste0("cluster", col_anno$cells)

# Arrange the pheat_data matrix by the order of col_anno table
pheat_data <- pheat_data[,rownames(col_anno)]

# Set legend color
colors <- c(cluster0="#a47d4c",cluster1="#2962ff",cluster2="#c30101",
            cluster3="#32cd32",cluster4="#7e5aa2",cluster5="#ffd7cd",
            cluster6="#03989e",cluster7="#c8b9ac",cluster8="#d61b5d",
            cluster9="#848c2f",cluster10="#bce3dd",cluster11="#383838",
            cluster12="#a2665c")
my_colour = list(
  cells = colors,
  markers = colors
)

# Plot the heatmap by pheatmap function
pheatmap(pheat_data, annotation_row = row_anno, annotation_col = col_anno, 
         show_colnames = F, cluster_rows = F, cluster_cols = F,
         annotation_colors = my_colour, annotation_names_row = F,
         annotation_legend = T, gaps_row = row_gap_list, fontsize_row = 3,
         main = paste0("VS4 Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))