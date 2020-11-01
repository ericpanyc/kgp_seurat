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
library(circlize)


setwd("/Users/ericpan/Documents/pan/usc/kgp/seurat")
use_python("/usr/local/bin/python3")
hpca.se <- HumanPrimaryCellAtlasData()
options(future.globals.maxSize = 4000 * 1024^2)

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


object_list <- c("VS4","VS6")
# object_list <- c("VS4","VS5","VS6","VS8","VS9")
object_list <- paste0(object_list, "_filtered_feature_bc_matrix")

# Read in 10x outputs files
for (file in object_list){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# naive <- merge(VS5_filtered_feature_bc_matrix, y = VS8_filtered_feature_bc_matrix, add.cell.ids = c("VS5", "VS8"), project = "naive_group")
cystic <- merge(VS4_filtered_feature_bc_matrix, y = VS6_filtered_feature_bc_matrix, add.cell.ids = c("VS4", "VS6"))

# naive@meta.data$orig.ident <- sapply(naive@meta.data$orig.ident, function(i) substr(i,1,3))
cystic@meta.data$orig.ident <- sapply(cystic@meta.data$orig.ident, function(i) substr(i,1,3))

# all <- merge(x = naive, y = cystic, add.cell.ids = c("N","S"), project = "all")

# Add number of genes per UMI for each cell to metadata
cystic$log10GenesPerUMI <- log10(cystic$nFeature_RNA) / log10(cystic$nCount_RNA)

# Compute percent mito ratio
cystic$mitoRatio <- PercentageFeatureSet(object = cystic, pattern = "^MT-")
cystic$mitoRatio <- cystic@meta.data$mitoRatio / 100

metadata <- cystic@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$sample <- NA
metadata$sample <- metadata$orig.ident
# metadata$sample[which(str_detect(metadata$cells, "^N_"))] <- "naive"
# metadata$sample[which(str_detect(metadata$cells, "^S_"))] <- "cystic"

cystic@meta.data <- metadata

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Filter out low-quality cells
filtered_seurat <- subset(x = cystic, 
                          subset=
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Gene level filtering 
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
# remove(counts)
# remove(nonzero)
# remove(filtered_counts)
# remove(cystic)
# remove(VS4_filtered_feature_bc_matrix)
# remove(VS6_filtered_feature_bc_matrix)
# remove(VS9_filtered_feature_bc_matrix)
# saveRDS(filtered_seurat, "data/cystic_46_filtered_before_integration.RDS")
cystic_unintegrated <- readRDS("data/cystic_46_filtered_before_integration.RDS")

# Split the merge seurat object
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

# SCTransform on each seurat object
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Integration 2 samples
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# saveRDS(seurat_integrated, "data/VS4_VS6_integrated_seurat.RDS")
seurat_integrated <- readRDS("data/VS4_VS6_integrated_seurat.RDS")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
# UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

# Cluster
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:30)
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(0,0.1,0.4,0.6,1.0))
# clustree(seurat_integrated, prefix = "integrated_snn_res.")

# Set the clustering resolution to 1
required_reso <- "integrated_snn_res.1"
Idents(object = seurat_integrated) <- "integrated_snn_res.1"

# UMAP Plot
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 3, split.by = "sample")
# UMAP group by sample
DimPlot(seurat_integrated, label = F, group.by = "sample", pt.size = 0.2)
# UMAP
DimPlot(seurat_integrated, label = T)

# Preliminary SingleR annotation
# Create a SingleCellExperiment version of seurat object
# meta_data <- seurat_integrated@meta.data
DefaultAssay(seurat_integrated) <- "RNA"
cystic_sc <- as.SingleCellExperiment(seurat_integrated)
colnames(cystic_sc) <- colData(cystic_sc)$cells
# Get all cell id in a cluster
cluster_barcode <- lapply(levels(seurat_integrated@meta.data$integrated_snn_res.1), function(i) seurat_integrated@meta.data[which(seurat_integrated@meta.data$integrated_snn_res.1 == i),"cells"])
# cluster level annotation fine labels
pred.cluster_fine <- SingleR(test = cystic_sc, ref = hpca.se, labels = hpca.se$label.fine, method = "cluster", clusters = seurat_integrated@meta.data$integrated_snn_res.1)
cluster_anno <- pred.cluster_fine$labels
names(cluster_anno) <- paste0("cluster",c(0:(length(cluster_anno)-1)))
cluster_anno <- as.data.frame(cluster_anno)
# write.xlsx(cluster_anno,"results/reso_1.0/cystic_46_singleR_cluster_level_anno.xlsx")
# cell level annotation fine labels
pred.cell_fine <- lapply(cluster_barcode, function(i) cystic_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.fine))
names(pred.cell_fine) <- c(0:(length(pred.cell_fine)-1))
# fine_labels <- sapply(pred.cell_fine, function(i) table(i$labels) %>% as.data.frame() %>% 
#                         arrange(desc(Freq)) %>% pull(1) %>% as.character() %>% head(3)) %>% as.data.frame()
# fine_labels <- t(fine_labels) %>% as.data.frame()
# saveRDS(pred.cell_fine,"results/reso_1.0/cystic_46_singleR_cell_level_anno.RDS")

new.cluster.ids <- c("Schwann_cell","Macrophage","Macrophage","Macrophage","Macrophage","Macrophage","Macrophage",
                     "Macrophage","Macrophage","T_cell","DC","Schwann_cell","Schwann_cell","Schwann_cell","Schwann_cell",
                     "T_cell","NK","Schwann_cell","Proliferating_cell","Endothelial_cell","Macrophage","Fibroblast")
# annotation <- new.cluster.ids
# cluster <- c(0:(length(annotation)-1))
# annos <- cbind(cluster, annotation) %>% as.data.frame()
# write.xlsx(annos,"results/reso_1.0/VS4&6_cluster_annotation.xlsx", row.names = F)

names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
seurat_integrated@meta.data$anno <- Idents(seurat_integrated)
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 3)
# meta_data$anno <- Idents(seurat_integrated)

# Find cluster marker genes
DefaultAssay(seurat_integrated) <- "RNA"
all_markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)
top30_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
top30_markers <- left_join(x = top30_markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                                 by = c("gene" = "gene_name"))
top30_markers <- as.data.frame(top30_markers)
# write.xlsx(top30_markers,"results/reso_1.0/cystic_46_cluster_top30_marker.xlsx")

# GO annotation
markers_by_cluster <- lapply(unique(top30_markers$cluster), function(i) top30_markers[top30_markers$cluster == i,])
markers <- lapply(markers_by_cluster, function(i) i[[9]])
cluster_go_annotation <- lapply(markers, function(i) goana(i) %>% topGO(n=15))
cluster_go_annotation <- lapply(cluster_go_annotation, function(i) i[[1]])
names(cluster_go_annotation) <- c(0:(length(cluster_go_annotation)-1))
df_go_annotation <- as.data.frame(cluster_go_annotation)
colnames(df_go_annotation) <- paste0("cluster",c(0:(length(cluster_go_annotation)-1)))
# write.xlsx(df_go_annotation, "results/reso_1.0/cystic_46_cluster_GO_pathway.xlsx", row.names = F)


cluster9.markers <- FindMarkers(seurat_integrated, ident.1 = 9, min.pct = 0.25)
top30_cluster9.markers <- cluster9.markers %>% rownames_to_column(var = "gene") %>% top_n(n = 50, wt = avg_logFC)
write.xlsx(top30_cluster9.markers,"results/VS4_m_VS6/cluster9_top50_genes.xlsx")

cluster9.markers_compare <- FindMarkers(seurat_integrated, ident.1 = 9, ident.2 = c(3,4,8), min.pct = 0.25, only.pos = T)
top30_cluster9.markers_relative <- cluster9.markers_compare %>% arrange(desc(avg_logFC)) %>% head(n=20) %>% rownames_to_column(var = "gene")
write.xlsx(top30_cluster9.markers_relative,"results/VS4_m_VS6/cluster9_top20_relative_genes.xlsx")

FeaturePlot(seurat_integrated, features = c("LARP6"),label = T)

DimPlot(seurat_integrated, label = T) + NoLegend()
VlnPlot(seurat_integrated,features = "LARP6",pt.size = 0) + NoLegend()


cystic_reso_1_cluster12_marker <- FindMarkers(seurat_integrated, ident.1 = 12, min.pct = 0.25, only.pos = T)
top30_cluster12_marker <- cystic_reso_1_cluster12_marker %>% rownames_to_column(var = "gene") %>% top_n(n = 30, wt = avg_logFC)
 
p1 <- DimPlot(seurat_integrated, label = F, group.by = "sample", pt.size = 0.2)
p1 <- p1[[1]]
p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, pt.size = 0.2,label.size = 6)
p2 <- p2[[1]]
p3 <- FeaturePlot(seurat_integrated, 
                  reduction = "umap", 
                  features = "nUMI",
                  pt.size = 0.4, 
                  order = T,
                  min.cutoff = 'q10',
                  label = F)
p3 <- p3[[1]]


p1 <- p1+ggtitle("Patient") + border() +
  theme(legend.position = c(0.1, 0.15),
        plot.title = element_text(size=30, vjust = -7,hjust = 0.05),
        legend.text = element_text(size=20,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1

p2 <- p2+ggtitle("Cell type") + border() + 
  theme(plot.title = element_text(size=30, vjust = -7,hjust = 0.05),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p2

p3 <- p3+ggtitle("Transcript count") + border() + 
  theme(legend.position = c(0.1, 0.15),
        plot.title = element_text(size=30, vjust = -7,hjust = 0.05),
        legend.text = element_text(size=20,face="bold"))

p3

p1/p2/p3

# Marker gene plot
DefaultAssay(seurat_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

# Idents(seurat_integrated) <- "integrated_snn_res.0.4"

# TTC25, ARL13B, MCDIC, LARP6, RFX7, RFX4, FOXJ1

cilia_plot_vln_cys <- VlnPlot(seurat_integrated, features = c("TTC25","ARL13B","LARP6","RFX7","RFX4","FOXJ1"), pt.size = 0, combine = F)
cystic_ttc25 <- cilia_plot_vln_cys[[1]]+ggtitle("VS4&VS6: TTC25")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_arl13b <- cilia_plot_vln_cys[[2]]+ggtitle("VS4&VS6: ARL13B")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_larp6 <- cilia_plot_vln_cys[[3]]+ggtitle("VS4&VS6: LARP6")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_rfx7 <- cilia_plot_vln_cys[[4]]+ggtitle("VS4&VS6: RFX7")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_rfx4 <- cilia_plot_vln_cys[[5]]+ggtitle("VS4&VS6: RFX4")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_foxj1 <- cilia_plot_vln_cys[[6]]+ggtitle("VS4&VS6: FOXJ1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

fibrosis_plot_vln_cys <- VlnPlot(seurat_integrated, features = c("STRAP","COL1A1","COL1A2","FKBP3","DHX9"), pt.size = 0, combine = F)
cystic_strap <- fibrosis_plot_vln_cys[[1]]+ggtitle("VS4&VS6: STRAP")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_col1a1 <- fibrosis_plot_vln_cys[[2]]+ggtitle("VS4&VS6: COL1A1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_col1a2 <- fibrosis_plot_vln_cys[[3]]+ggtitle("VS4&VS6: COL1A2")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_fkbp3 <- fibrosis_plot_vln_cys[[4]]+ggtitle("VS4&VS6: FKBP3")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_dhx9 <- fibrosis_plot_vln_cys[[5]]+ggtitle("VS4&VS6: DHX9")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

cystic_strap/cystic_col1a1/cystic_col1a2/cystic_fkbp3/cystic_dhx9

cystic_schwannoma_marker_gene <- VlnPlot(seurat_integrated, features = c("CD74","OOEP","NRG1","SEMA5A","ADGRB3","PDE4D","C1R","CLU","LINGO1","LGI4"), pt.size = 0, combine = F)
cystic_CD74 <- cystic_schwannoma_marker_gene[[1]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_OOEP <- cystic_schwannoma_marker_gene[[2]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_NRG1 <- cystic_schwannoma_marker_gene[[3]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_SEMA5A <- cystic_schwannoma_marker_gene[[4]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_ADGRB3 <- cystic_schwannoma_marker_gene[[5]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_PDE4D <- cystic_schwannoma_marker_gene[[6]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_C1R <- cystic_schwannoma_marker_gene[[7]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_CLU <- cystic_schwannoma_marker_gene[[8]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_LINGO1 <- cystic_schwannoma_marker_gene[[9]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_LGI4 <- cystic_schwannoma_marker_gene[[10]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

(cystic_CD74|cystic_OOEP)/(cystic_NRG1|cystic_SEMA5A)/(cystic_ADGRB3|cystic_PDE4D)/(cystic_C1R|cystic_CLU)/(cystic_LINGO1|cystic_LGI4)


marker_plot <- FeaturePlot(seurat_integrated, 
                           reduction = "umap", 
                           features = c("CLDN5","COL1A1","LYZ","CD3D","MPZ","GNLY"), 
                           min.cutoff = 'q10',
                           cols = c("azure3", "firebrick2"),
                           label = F,
                           pt.size = 0.4,
                           combine = F)
marker_plot <- lapply(marker_plot, function(x){x + theme(legend.position = "none",axis.title = element_blank(),
                                                         title = element_text(size=8),
                                                         axis.text = element_blank(),
                                                         axis.ticks = element_blank(),
                                                         plot.margin = margin(0, 0, 0, 0, "cm")) + border() + ggtitle("")})

p_endo <- marker_plot[[1]]+ggtitle("Endothelial: CLDN5")+theme(plot.title = element_text(vjust = -15,hjust = 0.95))
p_fibro <- marker_plot[[2]]+ggtitle("Fibroblast: COL1A1")+theme(plot.title = element_text(vjust = -15,hjust = 0.95))
p_macro <- marker_plot[[3]]+ggtitle("Macrophage: LYZ")+theme(plot.title = element_text(vjust = -15,hjust = 0.95))
p_Tcell <- marker_plot[[4]]+ggtitle("T-cell: CD3D")+theme(plot.title = element_text(vjust = -10,hjust = 0.95))
p_sch <- marker_plot[[5]]+ggtitle("Schwann-cell: MPZ")+theme(plot.title = element_text(vjust = -10,hjust = 0.95))
p_nk <- marker_plot[[6]]+ggtitle("NK: GNLY")+theme(plot.title = element_text(vjust = -10,hjust = 0.95))

p_fibro

sig_gene_plot <- FeaturePlot(seurat_integrated, 
                             reduction = "umap", 
                             features = c("PTPRC","CD68","AIF1","PTGS2","TGFB1","IL1B",
                                          "IL6","TNF","VEGFA","ICAM1","CD70","MKI67"), 
                             min.cutoff = 'q10', 
                             label = F,
                             cols = c("azure3", "firebrick2"),
                             pt.size = 0.4,
                             combine = F)
sig_gene_list <- lapply(sig_gene_plot, function(x){x + theme(legend.position = "none",axis.title = element_blank(),
                                                             axis.text = element_blank(),
                                                             axis.ticks = element_blank(),
                                                             plot.title = element_text(vjust = -7,hjust = 0.95),
                                                             plot.margin = margin(0, 0, 0, 0, "cm")) + border()})
p_ptprc <- sig_gene_list[[1]]
p_cd68 <- sig_gene_list[[2]]
p_aif1 <- sig_gene_list[[3]]
p_ptgs2 <- sig_gene_list[[4]]
p_tgfb1 <- sig_gene_list[[5]]
p_il1b <- sig_gene_list[[6]]
p_il6 <- sig_gene_list[[7]]
p_tnf <- sig_gene_list[[8]]
p_vegfa <- sig_gene_list[[9]]
p_icam1 <- sig_gene_list[[10]]
p_cd70 <- sig_gene_list[[11]]
p_mki67 <- sig_gene_list[[12]]


p_mid_top <- wrap_plots(p_endo,p_fibro,p_macro,p_Tcell,p_sch,p_nk,ncol = 3,nrow = 2)
p_mid_mid <- wrap_plots(p_ptprc,p_cd68,p_aif1,p_ptgs2,p_tgfb1,p_il1b,ncol = 3,nrow = 2)
p_mid_bot <- wrap_plots(p_il6,p_tnf,p_vegfa,p_icam1,p_cd70,p_mki67,ncol = 3,nrow = 2)

p_mid_top/p_mid_mid/p_mid_bot

# Plot stacked bar chart
macro_sub <- c(1:6,7,8,20)
sch_sub <- c(0,11:14,17)
t_sub <- c(9,15)

meta_data <- seurat_integrated@meta.data
bar_df <- meta_data[,c("sample","integrated_snn_res.1","anno")]
macro_df <- bar_df[which(bar_df$anno=="Macrophage"),]
sample <- c(rep("VS4",9),rep("VS6",9))
cluster <- rep(macro_sub,2)
df <- cbind(sample,cluster) %>% as.data.frame()
df$percent <- NA
df$sum <- NA

value <- vector()
for (s in c("VS4","VS6")) {
  for (c in macro_sub) {
    num <- nrow(macro_df[macro_df$sample == s & macro_df$integrated_snn_res.1 == c,])
    value <- c(value, num)
  }
}
df$sum <- value

percent <- vector()
for (x in c(1:nrow(df))) {
  curr_cluster <- as.integer(df[x,"cluster"])
  val <- df[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent <- c(percent, val)
}
df$percent <- percent
df$cluster_subtype <- as.factor(rep(1:9,2))
df$anno <- rep("Macrophage",nrow(df))

fill <- c("hotpink4","limegreen")


sch_df <- bar_df[which(bar_df$anno=="Schwann_cell"),]
sample_sch <- c(rep("VS4",6),rep("VS6",6))
cluster_sch <- rep(sch_sub,2)
df_sch <- cbind(sample_sch,cluster_sch) %>% as.data.frame()
df_sch$percent <- NA
df_sch$sum <- NA

value_sch <- vector()
for (s in c("VS4","VS6")) {
  for (c in sch_sub) {
    num <- nrow(sch_df[sch_df$sample == s & sch_df$integrated_snn_res.1 == c,])
    value_sch <- c(value_sch, num)
  }
}
df_sch$sum <- value_sch

percent_sch <- vector()
for (x in c(1:nrow(df_sch))) {
  curr_cluster <- as.integer(df_sch[x,"cluster_sch"])
  val <- df_sch[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent_sch <- c(percent_sch, val)
}
df_sch$percent <- percent_sch
df_sch$cluster_subytpe <- as.factor(rep(1:6,2))
df_sch$anno <- rep("Schwann_cell",nrow(df_sch))

t_df <- bar_df[which(bar_df$anno=="T_cell"),]
sample_t <- c(rep("VS4",2),rep("VS6",2))
cluster_t <- rep(t_sub,2)
df_t <- cbind(sample_t,cluster_t) %>% as.data.frame()
df_t$percent <- NA
df_t$sum <- NA

value_t <- vector()
for (s in c("VS4","VS6")) {
  for (c in t_sub) {
    num <- nrow(t_df[t_df$sample == s & t_df$integrated_snn_res.1 == c,])
    value_t <- c(value_t, num)
  }
}
df_t$sum <- value_t

percent_t <- vector()
for (x in c(1:nrow(df_t))) {
  curr_cluster <- as.integer(df_t[x,"cluster_t"])
  val <- df_t[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent_t <- c(percent_t, val)
}
df_t$percent <- percent_t
df_t$cluster_subtype <- as.factor(rep(1:2,2))
df_t$anno <- rep("T_cell",nrow(df_t))


names(df_sch) <- names(df)
names(df_t) <- names(df)
merged_df <- rbind(df,df_sch,df_t)

p9 <- ggplot(merged_df, aes( x = cluster_subtype, y = percent, fill = sample ) ) + 
  geom_bar( stat = "identity" ) + scale_fill_manual(values=fill) +
  facet_grid(. ~ anno,
             scales="free",
             space = "free") +
  theme(panel.background = element_blank(), legend.position="right",
        axis.text.y = element_text(size=20, angle = 90),
        axis.title.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        plot.margin = margin(0, 0, 0, 0.5, "cm"))
p9


cluster_sum_macro <- sapply(cluster_barcode[(macro_sub+1)], function(x) length(x))
cluster_sum_sch <- sapply(cluster_barcode[(sch_sub+1)], function(x) length(x))
cluster_sum_t <- sapply(cluster_barcode[(t_sub+1)], function(x) length(x))

subtype_bar <- data.frame(cluster_subtype = c(1:9,1:6,1:2),
                          count = c(cluster_sum_macro,cluster_sum_sch,cluster_sum_t),
                          anno = c(rep("Macrophage",9),rep("Schwann_cell",6),rep("T_cell",2)))


my_breaks <- function(x) { if (max(x) < 4) seq(0,2.5,1) else seq(0,12.5,1) }
p10 <- ggplot(subtype_bar, aes( x = cluster_subtype, y = count)) + 
  geom_bar( stat = "identity",fill="steelblue") +
  facet_grid(. ~ anno,
             scales="free",
             space = "free") +
  scale_x_continuous(breaks = my_breaks)+
  theme(legend.position="none",
        axis.text.y = element_text(size=20, angle = 90),
        axis.title.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0.5, "cm"))

p10

box_macro <- meta_data[which(meta_data$anno=="Macrophage"),c("nGene","integrated_snn_res.1","anno")]
box_macro$integrated_snn_res.1 <- ifelse(box_macro$integrated_snn_res.1=="1","1",
                                         ifelse(box_macro$integrated_snn_res.1=="2","2",
                                         ifelse(box_macro$integrated_snn_res.1=="3","3",
                                        ifelse(box_macro$integrated_snn_res.1=="4","4",
                                        ifelse(box_macro$integrated_snn_res.1=="5","5",
                                        ifelse(box_macro$integrated_snn_res.1=="6","6",
                                        ifelse(box_macro$integrated_snn_res.1=="7","7",
                                        ifelse(box_macro$integrated_snn_res.1=="8","8","9"))))))))
box_sch <- meta_data[which(meta_data$anno=="Schwann_cell"),c("nGene","integrated_snn_res.1","anno")]
box_sch$integrated_snn_res.1 <- ifelse(box_sch$integrated_snn_res.1=="0","1",
                                       ifelse(box_sch$integrated_snn_res.1=="11","2",
                                              ifelse(box_sch$integrated_snn_res.1=="12","3",
                                                     ifelse(box_sch$integrated_snn_res.1=="13","4",
                                                            ifelse(box_sch$integrated_snn_res.1=="14","5","6")))))
                                                                  
box_t <- meta_data[which(meta_data$anno=="T_cell"),c("nGene","integrated_snn_res.1","anno")]
box_t$integrated_snn_res.1 <- ifelse(box_t$integrated_snn_res.1=="9","1","2")


box_merge <- rbind(box_macro,box_sch,box_t)
box_merge$nGene_log <- log10(box_merge$nGene)
box_merge <- dplyr::rename(box_merge, cluster_subtype = integrated_snn_res.1)
box_merge$anno_f = factor(box_merge$anno, levels=c('Macrophage','Schwann_cell','T_cell'))
box_merge$cluster_subtype <- factor(box_merge$cluster_subtype, levels = c(1:11))


p11 <- ggplot(box_merge, aes(x = cluster_subtype, y =nGene_log)) + 
  facet_grid(. ~ anno_f,scales="free",space = "free",switch = "x") +
  geom_boxplot(fill='#A4A4A4', color="darkblue",outlier.size=0.5, outlier.colour="red") + 
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, angle = 90),
        strip.text.x = element_text(size=16,face="bold"),
        strip.background.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0.5, "cm"))

p11

# p9/p10/p11

(p1/p2/p3)|(p_mid_top/p_mid_mid/p_mid_bot)|(p9/p10/p11)+plot_layout(widths = c(2,3,3,2))

#-------------------------------------------------------------------------------------------------------------
VS4_df <- seurat_integrated@meta.data[which(seurat_integrated@meta.data$sample=="VS4"),c("sample","integrated_snn_res.1","anno")]
VS6_df <- seurat_integrated@meta.data[which(seurat_integrated@meta.data$sample=="VS6"),c("sample","integrated_snn_res.1","anno")]

cell_type <- c("Macrophage","Schwann_cell","Fibroblast","Proliferating_cell","Endothelial_cell","T_cell","NK","DC","Myocyte")


VS4_cell_number <- vector("integer")
for (x in cell_type) {
  num <- nrow(VS4_df[VS4_df$anno==x,])
  VS4_cell_number <- c(VS4_cell_number,num)
}
VS6_cell_number <- vector("integer")
for (x in cell_type) {
  num <- nrow(VS6_df[VS6_df$anno==x,])
  VS6_cell_number <- c(VS6_cell_number,num)
}



VS4_percent = VS4_cell_number/nrow(VS4_df)
VS6_percent = VS6_cell_number/nrow(VS6_df)



new_df <- data.frame(cell_type = rep(cell_type,2),
                     cell_count = c(VS4_cell_number,VS6_cell_number),
                     percentage = c(VS4_percent,VS6_percent),
                     group = rep("cystic",18))







