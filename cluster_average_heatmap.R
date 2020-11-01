pbmc <- readRDS(file = "~/Downloads/pbmc3k_final.rds")


DefaultAssay(seurat_integrated) <- "RNA"

cluster.averages <- AverageExpression(seurat_integrated)
head(cluster.averages[["RNA"]][, 1:5])


cluster.averages <- AverageExpression(seurat_integrated, return.seurat = TRUE)
cluster.averages


DoHeatmap(cluster.averages, features = cystic_top.markers$gene, size = 3, 
          draw.lines = FALSE)


orig.levels <- levels(pbmc)
Idents(pbmc) <- gsub(pattern = " ", replacement = "_", x = Idents(pbmc))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(pbmc) <- orig.levels
cluster.averages <- AverageExpression(pbmc, return.seurat = TRUE)
cluster.averages

DoHeatmap(cluster.averages, features = unlist(TopFeatures(pbmc[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

filtered_seurat@meta.data$cluster <- seurat_integrated@meta.data$integrated_snn_res.1
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)

Idents(filtered_seurat) <- seurat_integrated@meta.data$integrated_snn_res.1
cluster.averages <- AverageExpression(filtered_seurat, return.seurat = TRUE)
cystic_average_heatmap <- DoHeatmap(cluster.averages, features = cystic_top.markers$gene, size = 3, draw.lines = FALSE, ) + NoLegend()
cystic_average_heatmap + theme(axis.text.y = element_text(size=3, angle = 0, hjust = 1.3))

ave_data <- cluster.averages@assays$RNA@scale.data
ave_metadata <- cluster.averages@meta.data

all_genes <- rownames(ave_data)
marker_gene_tb <- read.csv("data/marker_gene_table.csv",sep = ",", header = T,na.strings=c("","NA"))
af_up <- marker_gene_tb$AF_up[1:17]
af_down <- marker_gene_tb$AF_down[1:28]
target_gene_full <- c(af_up, af_down)
filter <- all_genes %in% target_gene_full
gene_table <- ave_data[filter,]
dim(gene_table)

genes <- rownames(gene_table)
gene_anno <- cbind(genes, type = vector(mode = "character", length = 42)) %>% data.frame()
gene_anno$type <- ifelse(gene_anno$genes %in% af_up, "up", "down")

col_fun = colorRamp2(c(-3, 1, 5), c("dodgerblue4", "white", "red"))
# getPalette = colorRampPalette(brewer.pal(8, "Set2"))
cystic_metadata <- seurat_integrated@meta.data

Heatmap(gene_table,
        col = col_fun,
        name="AF_gene_heatmap",
        show_column_names = TRUE,
        cluster_rows = FALSE, cluster_columns = F,
        row_split = gene_anno$type, row_title_gp = gpar(font = 1, col = "black"),
        #row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)


filtered_seurat@meta.data$cluster <- seurat_integrated_naive@meta.data$integrated_snn_res.1
Idents(filtered_seurat) <- seurat_integrated_naive@meta.data$integrated_snn_res.1

filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)

cluster.averages <- AverageExpression(filtered_seurat, return.seurat = TRUE)
cystic_average_heatmap <- DoHeatmap(cluster.averages, features = top30_markers_naive$gene, size = 3, draw.lines = FALSE, ) + NoLegend()
cystic_average_heatmap + theme(axis.text.y = element_text(size=3, angle = 0, hjust = 1.3))

ave_data <- cluster.averages@assays$RNA@scale.data
filter <- all_genes %in% top30_markers_naive$gene
gene_table <- ave_data[filter,]

gene_table <- ave_data[top30_markers_naive$gene,]

Heatmap(gene_table,
        col = col_fun,
        name="AF_gene_heatmap",
        show_column_names = TRUE,
        cluster_rows = FALSE, cluster_columns = F,
        #row_split = top30_markers_naive$cluster, row_title_gp = gpar(font = 1, col = "black"),
        #row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)
