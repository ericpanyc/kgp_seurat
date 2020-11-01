target_gene <- c("SOX9","PYCR1","CKAP4","STXBP6","SLIT2","NPY1R","CENPW","SCRN1","GPM6B","FOXF1")
type <- c(rep("up",4),rep("down",6))
gene_anno <- cbind(target_gene,type) %>% as.data.frame()

# target_gene_counts <- GetAssayData(naive58, slot = "counts") %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% dplyr::filter(gene %in% target_gene)
# target_gene_data <- GetAssayData(naive58[["RNA"]], slot = "data") %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% dplyr::filter(gene %in% target_gene)
# 
# target_gene_counts <- target_gene_counts %>% column_to_rownames(var="gene")
# target_gene_data <- target_gene_data %>% column_to_rownames(var="gene")
# 
# max(target_gene_counts)
# max(target_gene_data)

naive58 <- NormalizeData(naive58)
# target_gene_data_norm <- GetAssayData(naive58[["RNA"]], slot = "data") %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% dplyr::filter(gene %in% target_gene)
# target_gene_data_norm <- target_gene_data_norm %>% column_to_rownames(var="gene")
# max(target_gene_data_norm)
# 
# VlnPlot(naive58, features = target_gene, pt.size = 0)



Idents(naive58) <- seurat_integrated_naive@meta.data$integrated_snn_res.1
Idents(cystic_unintegrated) <- seurat_integrated@meta.data$integrated_snn_res.1

cluster.averages <- AverageExpression(naive58, return.seurat = TRUE)
cluster.averages <- AverageExpression(cystic_unintegrated, return.seurat = T)

ave_data <- cluster.averages@assays$RNA@scale.data
colnames(ave_data) <- paste0("cluster", colnames(ave_data))

all_genes <- rownames(ave_data)
col_fun = colorRamp2(c(-2, 1.5, 5), c("dodgerblue4", "white", "red"))


cilia_gene <- marker_gene_tb$Cilia
target_gene <- marker_gene_tb$Cilia

insulin_gene <- marker_gene_tb$Insulin[1:9]
target_gene <- marker_gene_tb$Insulin[1:9]

target_gene <- marker_gene_tb$Wnt[1:39]


target_gene <- target_gene[target_gene %in% all_genes]

ave_data_target <- ave_data[target_gene,]
max(ave_data_target)
min(ave_data_target)


## Normal heatmap 
Heatmap(ave_data_target,
        col = col_fun,
        name="naive_insulin_gene_heatmap",
        show_column_names = TRUE, column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
        cluster_rows = T, 
        cluster_columns = T,
        cluster_columns = cluster_within_group(mat, group)
        # row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)
## Heatmap with dendrogram
group = kmeans(t(ave_data_target), centers = 3)$cluster
Heatmap(ave_data_target,
        col = col_fun,
        name="cystic_AF_gene_heatmap",
        show_column_names = TRUE, column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
        cluster_rows = T, 
        cluster_columns = cluster_within_group(ave_data_target, group),
        column_title = "K-means 3 clustering"
        # row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)


# Gene level filtering 
counts <- GetAssayData(object = filtered_seurat, slot = "counts")



nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_counts <- counts[target_gene,]


filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
filtered_seurat <- NormalizeData(filtered_seurat)

filtered_seurat@meta.data$cluster <- seurat_integrated_naive@meta.data$integrated_snn_res.1
Idents(filtered_seurat) <- seurat_integrated_naive@meta.data$integrated_snn_res.1
filtered_seurat_ave <- AverageExpression(filtered_seurat, return.seurat = T)
filtered_seurat_ave_data <- filtered_seurat_ave@assays$RNA@scale.data


naive_heatmap_table <- GetAssayData(filtered_seurat, slot = "counts") %>% as.matrix()

cluster_cell_counts <- table(Idents(filtered_seurat)) %>% as.data.frame()
colnames(cluster_cell_counts) <- c("cluster","cell_counts")

meta_data <- filtered_seurat@meta.data

mtx <- NULL
for (gene in c("SOX9","PYCR1","CKAP4","STXBP6")) {
  gene_counts_cluster <- NULL
  for (cluster in cluster_cell_counts$cluster) {
    target_cluster_cells <- meta_data[meta_data$cluster==cluster,] %>% rownames()
    target_gene_counts <- naive_heatmap_table[gene, target_cluster_cells] %>% sum()
    target_gene_per_cell <- target_gene_counts / cluster_cell_counts[cluster_cell_counts$cluster==11,2]
    gene_counts_cluster <- c(gene_counts_cluster, target_gene_per_cell)
  }
  mtx <- rbind(mtx, gene_counts_cluster)
}
rownames(mtx) <- c("SOX9","PYCR1","CKAP4","STXBP6")
colnames(mtx) <- paste0("cluster", cluster_cell_counts$cluster)

mtx.rep <- mtx
# mtx.rep[mtx.rep==0] <- 0.00001
mtx.rep <- 100*mtx.rep
mtx.rep <- scale(mtx.rep)

# mtx <- log10(mtx)

gene_counts_cluster <- vector("integer")
for (cluster in cluster_cell_counts$cluster) {
  target_cluster_cells <- meta_data[meta_data$cluster==cluster,] %>% rownames()
  target_gene_counts <- naive_heatmap_table["SOX9", target_cluster_cells] %>% sum()
  target_gene_per_cell <- target_gene_counts / cluster_cell_counts[cluster_cell_counts$cluster==11,2]
}



cluster.averages <- AverageExpression(filtered_seurat, return.seurat = TRUE)
ave_data <- cluster.averages@assays$RNA@scale.data



############################################################################
# Cystic Group Heatmap
cluster_anno <- seurat_integrated@meta.data$integrated_snn_res.1
Idents(cystic_unintegrated) <- cluster_anno

cluster.averages <- AverageExpression(cystic_unintegrated, return.seurat = TRUE)
ave_data <- cluster.averages@assays$RNA@scale.data
all_genes <- row.names(all_genes <- ave_data)
filter <- all_genes %in% target_gene_full
ave_data_target <- ave_data[filter,]

genes <- rownames(ave_data_target)
gene_anno_full <- cbind(genes, type = vector(mode = "character", length = 42)) %>% data.frame()
gene_anno_full$type <- ifelse(gene_anno_full$genes %in% af_up, "up", "down")

max(ave_data_target)
min(ave_data_target)

genes_subset <- gene_anno$target_gene
ave_data_subset <- ave_data[genes_subset,]

max(ave_data_subset)
min(ave_data_subset)


col_fun = colorRamp2(c(-2, 1.5, 5), c("dodgerblue4", "white", "red"))
Heatmap(ave_data_subset,
        col = col_fun,
        name="cystic_AF_gene_heatmap",
        show_column_names = TRUE, column_names_rot = 0, column_names_centered = TRUE,
        cluster_rows = F, cluster_columns = F,
        row_split = gene_anno$type, row_title_gp = gpar(font = 1, col = "black"),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)
