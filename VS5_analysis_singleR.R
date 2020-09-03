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
VS5_filtered_seurat <- RunUMAP(VS5_filtered_seurat, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

VS5_filtered_seurat <- FindNeighbors(VS5_filtered_seurat, dims = 1:30)
VS5_filtered_seurat <- FindClusters(object = VS5_filtered_seurat,
                                    resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

required_reso <- "SCT_snn_res.0.4"
Idents(object = VS5_filtered_seurat) <- required_reso

DimPlot(VS5_filtered_seurat, label = TRUE)

# saveRDS(VS5_filtered_seurat,"data/VS5_filtered_seurat.RDS")
#VS5_seurat <- readRDS("data/VS5_filtered_seurat.RDS")
#VS4_seurat <- readRDS("data/VS4_SCT_filtered.RDS")

#clustree(VS5_seurat, prefix = "SCT_snn_res.")

#seurat_list <- list(VS4_seurat, VS5_seurat)

meta_data <- VS5_seurat@meta.data

dense_table <- VS5_seurat[["SCT"]]@scale.data

Idents(object = VS5_seurat) <- required_reso

VS5.markers <- FindAllMarkers(VS5_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VS5_top.markers <- VS5.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
VS5_top.markers <- VS5_top.markers[which(VS5_top.markers$gene %in% row.names(dense_table)),]

dense_table <- dense_table[VS5_top.markers$gene,]

row_anno <- VS5_top.markers[,c("gene","cluster")] %>% data.frame()
rownames(row_anno) <- make.names(row_anno$gene, unique = T)
rownames(dense_table) <- rownames(row_anno)
row_anno <- row_anno[2]
names(row_anno) <- "markers"

markers_per_cluster <- lapply(unique(row_anno$markers), function(i) sum(row_anno$markers == i)) %>% flatten_int()
current_sum <- 0
row_gap_list <- lapply(c(1:(length(markers_per_cluster)-1)), function(i) current_sum <<- current_sum + markers_per_cluster[[i]]) %>% flatten_dbl()

col_anno <- meta_data[,c(required_reso), drop = F]
col_anno$cell <- row.names(col_anno)
col_anno <- col_anno[order(col_anno[,1]),]
col_anno <- col_anno[1]
names(col_anno) <- "cells"

cells_per_cluster <- lapply(unique(col_anno$cells), function(i) sum(col_anno$cells == i)) %>% flatten_int()
current_sum <- 0
col_gap_list <- lapply(c(1:(length(cells_per_cluster)-1)), function(i) current_sum <<- current_sum + cells_per_cluster[[i]]) %>% flatten_dbl()

row_anno$markers <- paste0("cluster", row_anno$markers) %>% as.vector()
col_anno$cells <- paste0("cluster", col_anno$cells) %>% as.vector()

dense_table <- dense_table[,rownames(col_anno)]

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

pheatmap(dense_table, annotation_row = row_anno, annotation_col = col_anno, 
         show_colnames = F, cluster_rows = F, cluster_cols = F,
         annotation_colors = my_colour, annotation_names_row = F,
         annotation_legend = T, gaps_row = row_gap_list, fontsize_row = 5,
         main = paste0("VS5 Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))

# Cluster annotation
hpca.se <- HumanPrimaryCellAtlasData()
mid.se <- MonacoImmuneData()

VS5_sc <- as.SingleCellExperiment(VS5_seurat)
colnames(VS5_sc) <- colData(VS5_sc)$barcode
VS5_cluster_barcode <- lapply(levels(meta_data$SCT_snn_res.0.4), function(i) meta_data[which(meta_data$SCT_snn_res.0.4 == i),"barcode"])
pred.VS5_results <- lapply(VS5_cluster_barcode, function(i) VS5_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.main))
                        
VS4_seurat@meta.data$barcode <- rownames(VS4_seurat@meta.data)
VS4_sc <- as.SingleCellExperiment(VS4_seurat)
colnames(VS4_sc) <- colData(VS4_sc)$barcode
VS4_cluster_barcode <- lapply(levels(VS4_seurat@meta.data$SCT_snn_res.0.4), function(i) VS4_seurat@meta.data[which(VS4_seurat@meta.data$SCT_snn_res.0.4 == i),"barcode"])
pred.VS4_results <- lapply(VS4_cluster_barcode, function(i) VS4_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.main))



names(pred.VS4_results) <- c(0:12)
pred.VS4_results_combined <- lapply(VS4_cluster_barcode, function(i) VS4_sc[,i] %>% SingleR(ref = list(MI=mid.se, HPCA=hpca.se),labels = list(mid.se$label.main, hpca.se$label.main)))
names(pred.VS4_results_combined) <- c(0:12)



   

