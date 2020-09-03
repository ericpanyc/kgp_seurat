seurat_integrated_naive <- NormalizeData(seurat_integrated_naive, verbose = FALSE)
dense_table <- seurat_integrated_naive[["RNA"]]@data %>% as.matrix()

metadata <- seurat_integrated_naive@meta.data

naive_heatmap_gene <- all_markers_naive %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pheat_data <- dense_table[naive_heatmap_gene$gene,]

col_fun = colorRamp2(c(0, 5, 9), c("lightblue", "white", "red"))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
size <- c(12,12,12,11,11,10,8,8,6,6,5,4,4,4,2,2,2,2,2,1,1,1,1,1,1,1) + 2

all_genes <- rownames(dense_table)
marker_gene_tb <- read.csv("data/marker_gene_table.csv",sep = ",", header = T,na.strings=c("","NA"))

af_up <- marker_gene_tb$AF_up[1:17]
af_down <- marker_gene_tb$AF_down[1:28]

target_gene <- c(af_up, af_down)
target_gene <- marker_gene_tb$Wnt[1:39]

filter <- all_genes %in% target_gene
gene_table <- dense_table[filter,]
gent_table <- gene_table[-which(rownames(gene_table)=="GPM6B"),]
dim(gene_table)

genes <- rownames(gene_table)
gene_anno <- cbind(genes, type = vector(mode = "character", length = 42)) %>% data.frame()
gene_anno$type <- ifelse(gene_anno$genes %in% af_up, "up", "down")
gent_anno <- gene_anno[-19,]

Heatmap(pheat_data,
        col = col_fun,
        name="naive_marker_gene_heatmap",
        row_names_gp = gpar(fontsize = 3),
        show_column_names = FALSE,
        cluster_rows = F, cluster_columns = F,
        row_split = naive_heatmap_gene$cluster, row_title_gp = gpar(font = 1, col = "black"),
        column_split = metadata$integrated_snn_res.1, column_title_gp = gpar(font = 1, col = "white"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = getPalette(26)),
                                                            labels = paste0("cluster",c(0:25)),
                                                            labels_gp = gpar(col = "black", fontsize = size))),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)

col_fun_2 = colorRamp2(c(0, 2.5, 5), c("lightblue", "white", "red"))
Heatmap(gent_table,
        col = col_fun_2,
        name="naive_AF_gene_heatmap",
        row_names_gp = gpar(fontsize = 8),
        show_column_names = FALSE,
        cluster_rows = F, cluster_columns = F,
        row_split = gent_anno$type, row_title_gp = gpar(fontsize = 10, col = "black"),
        column_split = metadata$integrated_snn_res.1, column_title_gp = gpar(font = 1, col = "white"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = getPalette(26)),
                                                            labels = paste0("cluster",c(0:25)),
                                                            labels_gp = gpar(col = "black", fontsize = size))),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
)
