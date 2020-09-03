## Marker Heatmap
library(pheatmap)
library(circlize)



#Suppose we have a seurat object named "seurat_integrated", which has gone through the SCTranform and UMAP clustering 
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)


# This will extract the raw counts for all genes (not used in the following plot)
non_scale_table <- seurat_integrated[["SCT"]]@data %>% as.matrix()

# This will extract the 3000 variable features identified by SCT
# dense_table <- seurat_integrated[["integrated"]]@scale.data
dense_table <- seurat_integrated[["RNA"]]@data %>% as.matrix()

# Set clustering resolution to 1
required_reso <- "integrated_snn_res.1"
Idents(object = seurat_integrated) <- required_reso


# Find top 20 marker genes for each cluster and exclude genes not present in the dense_table matrix

all_markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, logfc.threshold = 0.25)

#saveRDS(all_markers,"results/reso_1.0/cystic_46_all_markers.RDS")
cystic_all_markers <- readRDS("results/reso_1.0/cystic_46_all_markers.RDS")
cystic_top.markers <- cystic_all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Extract expressing data of marker genes only from dense_table, this table will be used for the final heatmap plot
cystic_cluster_data <- dense_table[cystic_top.markers$gene,]

# Build annotation table for row (separate marker genes of different clusters)
row_anno <- cystic_top.markers[,c("gene","cluster")] %>% data.frame()
rownames(row_anno) <- make.names(row_anno$gene, unique = T)
rownames(pheat_data) <- rownames(row_anno)
row_anno <- row_anno[2]
row_anno$cluster <- factor(row_anno$cluster)
names(row_anno) <- "markers"

# Set row gap list (get the number of marker genes in each cluster (accumulative))
markers_per_cluster <- sapply(unique(row_anno$markers), function(i) sum(row_anno$markers == i))
current_sum <- 0
row_gap_list <- sapply(c(1:(length(markers_per_cluster)-1)), function(i) current_sum <<- current_sum + markers_per_cluster[[i]])

# Build annotation table for column (separate cells from different clusters)
metadata <- seurat_integrated@meta.data
col_anno <- metadata[,c(required_reso), drop = F]
col_anno <- col_anno %>% arrange(integrated_snn_res.1)
col_anno$integrated_snn_res.1 <- factor(col_anno$integrated_snn_res.1)
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

# Plot the marker gene heatmap
all_genes <- rownames(dense_table)

cilia_gene <- read.csv("data/cilia_gene_table.csv",sep = ",", header = T)
cilia_gene <- cilia_gene$Genes
marker_gene_tb <- read.csv("data/marker_gene_table.csv",sep = ",", header = T,na.strings=c("","NA"))
af_up <- marker_gene_tb$AF_up[1:17]
af_down <- marker_gene_tb$AF_down[1:28]
target_gene <- c(af_up, af_down)
# target_gene <- marker_gene_tb$Wnt[1:39]


# target_gene <- c("MYC","MYCN","CCND1", "IGF1","IGF1R","IRS2", "LAD1", "CA2", "PYCR1")
# filter <- all_genes %in% target_gene
filter <- all_genes %in% target_gene
gene_table <- dense_table[filter,]
dim(gene_table)

genes <- rownames(gene_table)
gene_anno <- cbind(genes, type = vector(mode = "character", length = 42)) %>% data.frame()
gene_anno$type <- ifelse(gene_anno$genes %in% af_up, "up", "down")
# cilia_table <- dense_table[filter,rownames(col_anno)]



col_fun = colorRamp2(c(0, 3, 6), c("lightblue", "white", "red"))
# va = columnAnnotation(foo = anno_mark(at = c(0:21), labels = paste0("cluster_",c(0:21))))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
size <- c(12,12,12,11,11,10,8,8,6,6,5,4,4,4,2,2,2,2,2,1,1,1)
cystic_metadata <- seurat_integrated@meta.data

Heatmap(gene_table,
        col = col_fun,
        name="AF_gene_heatmap",
        show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = F,
        row_split = gene_anno$type, row_title_gp = gpar(font = 1, col = "black"),
        column_split = cystic_metadata$integrated_snn_res.1, column_title_gp = gpar(font = 1, col = "white"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = getPalette(22)),
                                                            labels = paste0("cluster",c(0:21)),
                                                            labels_gp = gpar(col = "black", fontsize = size))),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
        )


col_fun_2 = colorRamp2(c(0, 4.5, 9), c("lightblue", "white", "red"))
Heatmap(cystic_cluster_data,
        col = col_fun_2, name="cystic_AF_gene_heatmap", row_names_gp = gpar(fontsize = 5),
        show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = F,
        row_split = cystic_top.markers$cluster, row_title_gp = gpar(fontsize = 5, col = "black"),
        column_split = cystic_metadata$integrated_snn_res.1, column_title_gp = gpar(col = "white"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = getPalette(22)),
                                                            labels = paste0("cluster",c(0:21)),
                                                            labels_gp = gpar(col = "black", fontsize = size))),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)
