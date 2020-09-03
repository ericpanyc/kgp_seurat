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
library(SingleR)
library(scRNAseq)
library(scater)
library(DropletUtils)
library(stringr)
library(data.table)

setwd("/Users/ericpan/Documents/pan/usc/kgp/seurat")
use_python("/usr/local/bin/python3")

# example_file <- system.file("extra/TagSeqExample.tab", package="DESeq")
# data <- read.delim(example_file, header=T, row.names="gene")
# data_subset <- as.matrix(data[rowSums(data)>50000,])
# VS4.sc <- read10xCounts("./data/VS4_filtered_feature_bc_matrix")

VS4.data <- Read10X(data.dir = "./data/VS4_filtered_feature_bc_matrix")
VS4 <- CreateSeuratObject(counts = VS4.data, project = "VS4", min.features = 100)

VS4$log10GenesPerUMI <- log10(VS4$nFeature_RNA) / log10(VS4$nCount_RNA)

VS4$mitoRatio <- PercentageFeatureSet(object = VS4, pattern = "^MT-")
VS4$mitoRatio <- VS4@meta.data$mitoRatio / 100

metadata <- VS4@meta.data

metadata <- metadata %>%
  dplyr::rename(sample = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
VS4@meta.data <- metadata

filtered_seurat <- subset(x = VS4, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 400) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)



# -------------------------------------------------------------------------------------------------------
## Old Protocol 

# filtered_seurat <- NormalizeData(filtered_seurat)
# 
# filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(filtered_seurat)
# filtered_seurat <- ScaleData(filtered_seurat, features = all.genes)
# -------------------------------------------------------------------------------------------------------

# ----------------------------------
filtered_seurat <- SCTransform(filtered_seurat)

filtered_seurat <- RunPCA(filtered_seurat)
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:30)
filtered_seurat <- FindClusters(object = filtered_seurat,
                    resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# remove(counts)
# remove(filtered_counts)
# remove(nonzero)
# remove(VS4)
# remove(VS4.data)
# clustree(filtered_seurat, prefix = "SCT_snn_res.")

# saveRDS(filtered_seurat, file = "data/VS4_SCT_filtered.RDS")
filtered_seurat <- readRDS("data/VS4_SCT_filtered.RDS")
filtered_seurat@meta.data$barcode <- rownames(filtered_seurat@meta.data)
# Cell numbers
# Check the numbers of cells in each cluster
n_cells <- FetchData(filtered_seurat, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells <- n_cells %>% column_to_rownames(var="orig.ident") %>% t() %>% as.data.frame()
n_cells$V2 <- NA
n_cells$V2 <- (n_cells$SeuratProject)/(sum(n_cells$SeuratProject))
colnames(n_cells) <- c("cell_counts","percentage")

# SingleR annotation
hpca.se <- HumanPrimaryCellAtlasData()
# Create a SingleCellExperiment version of seurat object
DefaultAssay(filtered_seurat) <- "RNA"
VS4_sc <- as.SingleCellExperiment(filtered_seurat)
colnames(VS4_sc) <- colData(VS4_sc)$barcode
VS4_cluster_barcode <- lapply(levels(filtered_seurat@meta.data$SCT_snn_res.0.4), function(i) filtered_seurat@meta.data[which(filtered_seurat@meta.data$SCT_snn_res.0.4 == i),"barcode"])
pred.VS4_results <- lapply(VS4_cluster_barcode, function(i) VS4_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.main))
cluster_summary <- sapply(pred.VS4_results, function(i) table(i$labels) %>% as.data.frame() %>%
                            arrange(desc(Freq)) %>% pull(1) %>% as.character() %>% head(3)) %>% as.data.frame()
colnames(cluster_summary) <- paste0("cluster",c(0:12))
pred.VS4_results_fine <- lapply(VS4_cluster_barcode, function(i) VS4_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.fine))
cluster_summary_fine <- lapply(pred.VS4_results_fine, function(i) table(i$labels) %>% as.data.frame() %>% arrange(desc(Freq)) %>% pull(1) %>% as.character() %>% head(3)) %>% as.data.frame()
                   
cell_cluster <- as.factor(filtered_seurat@meta.data$SCT_snn_res.0.4)
pred.VS4_results_fine_cluster <- SingleR(test = VS4_sc, ref = hpca.se, labels = hpca.se$label.fine, method = "cluster", clusters = cell_cluster)



# find all markers distinguishing cluster 0 from clusters 1,3,4,6,9
cluster0_contrast_markers <- FindMarkers(filtered_seurat, ident.1 = 0, ident.2 = c(1,3,4,6,9), min.pct = 0.25, only.pos = T, logfc.threshold = 0.5)
cluster0_top_markers <- cluster0_contrast_markers %>% top_n(n = 30, wt = avg_logFC) %>% rownames_to_column(var = "gene")

## Marker Heatmap
#Suppose we have a seurat object named "filtered_seurat", which has gone through the SCTranform and UMAP clustering 
metadata <- filtered_seurat@meta.data

# This will extract the raw counts for all genes
non_scale_table <- filtered_seurat[["SCT"]]@data %>% as.matrix()

# This will extract the 3000 variable features indentified by SCT
dense_table <- filtered_seurat[["SCT"]]@scale.data

# Set clustering resolution to 0.4
required_reso <- "SCT_snn_res.0.4"
Idents(object = filtered_seurat) <- required_reso

# Find top 20 marker genes for each cluster and exclude genes not present in the dense_table matrix
VS4.markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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

# Plot the heatmap by pheatmap function
pheatmap(pheat_data, annotation_row = row_anno, annotation_col = col_anno, 
         show_colnames = F, cluster_rows = F, cluster_cols = F,
         annotation_colors = my_colour, annotation_names_row = F,
         annotation_legend = T, gaps_row = row_gap_list, fontsize_row = 3,
         main = paste0("VS4 Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))




# my_colour = list(
#   cells = c(cluster0="#a47d4c",cluster1="#2962ff",cluster2="#c30101",
#             cluster3="#32cd32",cluster4="#7e5aa2",cluster5="#ffd7cd",
#             cluster6="#03989e",cluster7="#c8b9ac",cluster8="#d61b5d",
#             cluster9="#848c2f",cluster10="#bce3dd",cluster11="#383838",
#             cluster12="#a2665c",cluster13="#5b56bc",cluster14="#33ccff",
#             cluster15="#cc0066",cluster16="#FFC300",cluster17="#388A5C",
#             cluster18="#9157A8",cluster19="#6030B9",cluster20="#CB5F9F"),
#   markers = c(cluster0="#a47d4c",cluster1="#2962ff",cluster2="#c30101",
#               cluster3="#32cd32",cluster4="#7e5aa2",cluster5="#ffd7cd",
#               cluster6="#03989e",cluster7="#c8b9ac",cluster8="#d61b5d",
#               cluster9="#848c2f",cluster10="#bce3dd",cluster11="#383838",
#               cluster12="#a2665c",cluster13="#5b56bc",cluster14="#33ccff",
#               cluster15="#cc0066",cluster16="#FFC300",cluster17="#388A5C",
#               cluster18="#9157A8",cluster19="#6030B9",cluster20="#CB5F9F")

DimPlot(filtered_seurat, label = T)

## Get annotation information
ah <- AnnotationHub()
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb, 
                     return.type = "data.frame")
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, entrezid,description)
annotations$entrezid <- as.character(annotations$entrezid)

## Integret annotations to the marker list
all.markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top20.markers <- left_join(x = top20.markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                         by = c("gene" = "gene_name"))

# Find GO pathway enrichment for each cluster (top 15 entries)
markers_by_cluster <- lapply(unique(top20.markers$cluster), function(i) top20.markers[top20.markers$cluster == i,])
markers_name <- lapply(markers_by_cluster, function(i) i[[7]])
markers_des <- lapply(markers_by_cluster, function(i) i[[8]])
markers <- lapply(markers_by_cluster, function(i) i[[9]])
cluster_go_annotation <- lapply(markers, function(i) goana(i) %>% topGO(n=15))
cluster_go_annotation <- lapply(cluster_go_annotation, function(i) i[[1]])
names(cluster_go_annotation) <- levels(top20.markers$cluster)
names(markers_name) <- levels(top20.markers$cluster)
names(markers_des) <- levels(top20.markers$cluster)

# -------------------------------------------------------------------------------------------------------

VS5.data <- Read10X(data.dir = "./data/VS5_filtered_feature_bc_matrix")
VS5 <- CreateSeuratObject(counts = VS5.data, project = "VS5", min.features = 100)

VS5$log10GenesPerUMI <- log10(VS5$nFeature_RNA) / log10(VS5$nCount_RNA)

VS5$mitoRatio <- PercentageFeatureSet(object = VS5, pattern = "^MT-")
VS5$mitoRatio <- VS5@meta.data$mitoRatio / 100

VS5_metadata <- VS5@meta.data

VS5_metadata <- VS5_metadata %>%
  dplyr::rename(sample = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
VS5@meta.data <- VS5_metadata

VS5_filtered_seurat <- subset(x = VS5, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 400) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

counts <- GetAssayData(object = VS5_filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
VS5_filtered_counts <- counts[keep_genes, ]
VS5_filtered_seurat <- CreateSeuratObject(VS5_filtered_counts, meta.data = VS5_filtered_seurat@meta.data)

# remove(counts)
# remove(VS5_filtered_counts)
# remove(nonzero)
# remove(VS5)
# remove(VS5.data)

VS5_filtered_seurat <- SCTransform(VS5_filtered_seurat)
VS5_filtered_seurat <- RunPCA(VS5_filtered_seurat)
VS5_filtered_seurat <- FindNeighbors(VS5_filtered_seurat, dims = 1:30)
VS5_filtered_seurat <- FindClusters(object = VS5_filtered_seurat,
                                resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
VS5_filtered_seurat <- RunUMAP(VS5_filtered_seurat, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
VS5_metadata <- VS5_filtered_seurat@meta.data

VS5_dense_table <- VS5_filtered_seurat[["SCT"]]@scale.data %>% as.data.frame()
Idents(object = VS5_filtered_seurat) <- required_reso

DimPlot(VS5_filtered_seurat, label = T)

VS5.markers <- FindAllMarkers(VS5_filtered_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VS5_top5.markers <- VS5.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

VS5_top30.markers <- VS5.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
VS5_top30.markers <- VS5_top30.markers[which(VS5_top30.markers$gene %in% row.names(VS5_dense_table)),]


VS5_pheat_data <- VS5_dense_table[VS5_top5.markers$gene,]

VS5_row_anno <- VS5_top5.markers[,c("gene","cluster")] %>% data.frame()
rownames(VS5_row_anno) <- make.names(VS5_row_anno$gene, unique = T)
rownames(VS5_pheat_data) <- rownames(VS5_row_anno)
VS5_row_anno <- VS5_row_anno[2]
names(VS5_row_anno) <- "markers"

VS5_markers_per_cluster <- lapply(unique(VS5_row_anno$markers), function(i) sum(VS5_row_anno$markers == i)) %>% flatten_int()
current_sum <- 0
VS5_row_gap_list <- lapply(c(1:(length(VS5_markers_per_cluster)-1)), function(i) current_sum <<- current_sum + VS5_markers_per_cluster[[i]]) %>% flatten_dbl()

VS5_col_anno <- VS5_metadata[,c(required_reso), drop = F]
VS5_col_anno$cell <- row.names(VS5_col_anno)
VS5_col_anno <- VS5_col_anno[order(VS5_col_anno[,1]),]
VS5_col_anno <- VS5_col_anno[1]
names(VS5_col_anno) <- "cells"

VS5_cells_per_cluster <- lapply(unique(VS5_col_anno$cells), function(i) sum(VS5_col_anno$cells == i)) %>% flatten_int()
current_sum <- 0
VS5_col_gap_list <- lapply(c(1:(length(VS5_cells_per_cluster)-1)), function(i) current_sum <<- current_sum + VS5_cells_per_cluster[[i]]) %>% flatten_dbl()

VS5_row_anno$markers <- paste0("cluster", VS5_row_anno$markers) %>% as.vector()
VS5_col_anno$cells <- paste0("cluster", VS5_col_anno$cells) %>% as.vector()

VS5_pheat_data <- VS5_pheat_data[,rownames(VS5_col_anno)]

pheatmap(VS5_pheat_data, annotation_row = VS5_row_anno, annotation_col = VS5_col_anno, 
         show_colnames = F, cluster_rows = F, cluster_cols = F,
         annotation_colors = my_colour, annotation_names_row = F,
         annotation_legend = T, gaps_row = VS5_row_gap_list,
         main = paste0("VS5 Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))

VS5_top20.markers <- VS5.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
VS5_top20.markers <- left_join(x = VS5_top20.markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                           by = c("gene" = "gene_name"))

# Find GO pathway enrichment for each cluster (top 15 entries)
VS5_markers_by_cluster <- lapply(unique(VS5_top20.markers$cluster), function(i) VS5_top20.markers[VS5_top20.markers$cluster == i,])
VS5_markers_name <- lapply(VS5_markers_by_cluster, function(i) i[[7]])
VS5_markers_des <- lapply(VS5_markers_by_cluster, function(i) i[[8]])
VS5_markers <- lapply(VS5_markers_by_cluster, function(i) i[[9]])
VS5_cluster_go_annotation <- lapply(VS5_markers, function(i) goana(i) %>% topGO(n=15))
VS5_cluster_go_annotation <- lapply(VS5_cluster_go_annotation, function(i) i[[1]])
names(VS5_cluster_go_annotation) <- levels(VS5_top20.markers$cluster)
names(VS5_markers_name) <- levels(VS5_top20.markers$cluster)
names(VS5_markers_des) <- levels(VS5_top20.markers$cluster)


# Feature plot
FeaturePlot(filtered_seurat, features = c("ALDH1A1"), label = T)
FeaturePlot(VS5_filtered_seurat, features = c("ALDH1A1"), label = T)

# SingleR part
hpca.se <- HumanPrimaryCellAtlasData()

hESCs <- LaMannoBrainData('human-es')

filtered_seurat@meta.data$orig.ident <- rownames(filtered_seurat@meta.data)
VS4.sc.filtered <- as.SingleCellExperiment(filtered_seurat)
colnames(VS4.sc.filtered) <- colData(VS4.sc.filtered)$orig.ident
# VS4.sc.test <- VS4.sc.filtered[,1:100]
# pred.VS4 <- SingleR(test = VS4.sc.test, ref = hpca.se, labels = hpca.se$label.main)
metadata <- filtered_seurat@meta.data
cluster_0_cells <- metadata[which(metadata$SCT_snn_res.0.4 == 0),1]

clutser_0.sc <- VS4.sc.filtered[,cluster_0_cells]

pred.VS4_cluster_0 <- SingleR(test = clutser_0.sc, ref = hpca.se, labels = hpca.se$label.main)
