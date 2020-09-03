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
library(janitor)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

# Get annotation information
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

options(future.globals.maxSize = 4000 * 1024^2)
setwd("/Users/ericpan/Documents/pan/usc/kgp/seurat")
use_python("/usr/local/bin/python3")
hpca.se <- HumanPrimaryCellAtlasData()


# Read in 10x outpus
samples <- c("VS4","VS5","VS6","VS8","VS9")

for (file in samples){
  folder <- paste0(file,"_filtered_feature_bc_matrix")
  seurat_data <- Read10X(data.dir = paste0("data/", folder))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}
# Combine all samples together
merged_seurat <- merge(x = VS4, y = c(VS5,VS6,VS8,VS9), add.cell.id = samples)

# Add more information in metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Process metadata
meta_data <- merged_seurat@meta.data
meta_data$group <- ifelse(meta_data$orig.ident=="VS4","cystic",
                          ifelse(meta_data$orig.ident=="VS6","cystic",
                                 ifelse(meta_data$orig.ident=="VS9","cystic","naive")))
meta_data <- meta_data %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

merged_seurat@meta.data <- meta_data

# Cell level filtering
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Gene level filtering
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Split the merged data by disease subtype
split_seurat <- SplitObject(filtered_seurat, split.by = "group")
split_seurat <- split_seurat[c("cystic", "naive")]

# SCTransform
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Integrate data
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# saveRDS(seurat_integrated,"data/cystic_naive_merge.RDS")

########################################################
cystic_naive <- readRDS("data/cystic_naive_merge.RDS")

cystic_naive <- RunPCA(cystic_naive, verbose = FALSE)
cystic_naive <- RunUMAP(cystic_naive, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

cystic_naive <- FindNeighbors(object = cystic_naive, dims = 1:30)
cystic_naive <- FindClusters(object = cystic_naive, resolution = c(0,0.1,0.2,0.4,0.6,0.8,1.0,1.4))

meta_data <- cystic_naive@meta.data 
Idents(cystic_naive) <- "integrated_snn_res.0.6"
DimPlot(cystic_naive, split.by = "group", label = T) + NoLegend()
DimPlot(cystic_naive, label = T) + NoLegend()

# Calculate cell numbers in each cluster
n_cells <- FetchData(cystic_naive, 
                     vars = c("ident", "group")) %>%
  dplyr::count(ident, group) %>%
  tidyr::spread(ident, n) %>% column_to_rownames(var = "group")

# Feature plot
metrics <-  c("nUMI", "nGene", "mitoRatio")
FeaturePlot(cystic_naive, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = T,
            min.cutoff = 'q10',
            label = TRUE)

# Select the RNA counts slot to be the default assay
DefaultAssay(cystic_naive) <- "RNA"
cystic_naive <- NormalizeData(cystic_naive, verbose = FALSE)

# Check the expression pattern of NK cells related genes
FeaturePlot(cystic_naive, 
            reduction = "umap", 
            features = c("GNLY", "KLRC1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Get conserved cell marker gene list
get_conserved <- function(cluster){
  FindConservedMarkers(cystic_naive,
                       ident.1 = cluster,
                       grouping.var = "group",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers <- map_dfr(c(0:19), get_conserved)
# write.xlsx(conserved_markers,"results/integrated_conserved_markers.xlsx")
# Get top 10 marker genes
top30 <- conserved_markers %>% 
  mutate(avg_fc = (cystic_avg_logFC + naive_avg_logFC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 30, 
        wt = avg_fc)

# Preliminary SingleR annotation
# Create a SingleCellExperiment version of seurat object
DefaultAssay(cystic_naive) <- "RNA"
cn_sc <- as.SingleCellExperiment(cystic_naive)
colnames(cn_sc) <- colData(cn_sc)$cells

# Get all cell id in a cluster
cluster_barcode <- lapply(levels(cystic_naive@meta.data$integrated_snn_res.0.6), function(i) cystic_naive@meta.data[which(cystic_naive@meta.data$integrated_snn_res.0.6 == i),"cells"])
# Cell level SingleR annotation
pred.cell_fine_results <- lapply(cluster_barcode, function(i) cn_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.fine))
names(pred.cell_fine_results) <- c(0:19)
# saveRDS(pred.cell_fine_results,"results/pred.cell_fine_results.RDS")
# pred.cell_fine_results <- pred.cluster_fine_results
# cluster level annotation fine labels
pred.cluster_fine <- SingleR(test = cn_sc, ref = hpca.se, labels = hpca.se$label.fine, method = "cluster", clusters = cystic_naive@meta.data$integrated_snn_res.0.6)
cluster_anno <- pred.cluster_fine$labels
names(cluster_anno) <- paste0("cluster",c(0:(length(cluster_anno)-1)))
cluster_anno <- as.data.frame(cluster_anno)
# write.xlsx(cluster_anno,"results/singleR_cluster_level_anno.xlsx")

# Select the RNA counts slot to be the default assay
DefaultAssay(cystic_naive) <- "RNA"
# Normalize RNA data for visualization purposes
cystic_naive <- NormalizeData(cystic_naive, verbose = FALSE)
group <- cystic_naive@meta.data$group
cystic_naive@meta.data$type <- group
cystic_naive@meta.data[cystic_naive@meta.data=="cystic"] <- "VS4&VS6&VS9"
cystic_naive@meta.data[cystic_naive@meta.data=="naive"] <- "VS5&VS8"
cystic_naive@meta.data$group <- group

# Cilia Genes of Interest - TTC25, ARL13B, MCDIC, LARP6, RFX7, RFX4, FOXJ1
cilia_plots <- VlnPlot(cystic_naive, features = c("TTC25","ARL13B","LARP6","RFX7","RFX4","FOXJ1"), split.by = "type", group.by = "integrated_snn_res.0.6", 
                 pt.size = 0, combine = FALSE, slot = "counts")
wrap_plots(plots = cilia_plots, ncol = 1)

# Fibrosis Genes - STRAP, COL1A1, COL1A2, FK506, LARP6, RHA
fibrosis_plots <- VlnPlot(cystic_naive, features = c("STRAP","COL1A1","COL1A2","FKBP5","LARP6","DHX9"), split.by = "type", group.by = "integrated_snn_res.0.6", 
                       pt.size = 0, combine = FALSE, slot = "counts")
wrap_plots(plots = fibrosis_plots, ncol = 1)
