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

object_list <- c("VS5","VS8")
object_list <- paste0(object_list, "_filtered_feature_bc_matrix")

# Read in 10x outputs files
for (file in object_list){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

naive <- merge(VS5_filtered_feature_bc_matrix, y = VS8_filtered_feature_bc_matrix, add.cell.ids = c("VS5", "VS8"), project = "naive")
# naive <- merge(VS4_filtered_feature_bc_matrix, y = c(VS6_filtered_feature_bc_matrix, VS9_filtered_feature_bc_matrix), add.cell.ids = c("VS4", "VS6", "VS9"), project = "naive")

naive@meta.data$orig.ident <- sapply(naive@meta.data$orig.ident, function(i) substr(i,1,3))
# naive@meta.data$orig.ident <- sapply(naive@meta.data$orig.ident, function(i) substr(i,1,3))

# all <- merge(x = naive, y = naive, add.cell.ids = c("N","S"), project = "all")

# Add number of genes per UMI for each cell to metadata
naive$log10GenesPerUMI <- log10(naive$nFeature_RNA) / log10(naive$nCount_RNA)

# Compute percent mito ratio
naive$mitoRatio <- PercentageFeatureSet(object = naive, pattern = "^MT-")
naive$mitoRatio <- naive@meta.data$mitoRatio / 100

metadata <- naive@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$sample <- NA
metadata$sample <- metadata$orig.ident
# metadata$sample[which(str_detect(metadata$cells, "^N_"))] <- "naive"
# metadata$sample[which(str_detect(metadata$cells, "^S_"))] <- "naive"

naive@meta.data <- metadata

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
filtered_seurat <- subset(x = naive, 
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
# remove(naive)
# remove(VS5_filtered_feature_bc_matrix)
# remove(VS8_filtered_feature_bc_matrix)
# remove(VS4_filtered_feature_bc_matrix)
# remove(VS6_filtered_feature_bc_matrix)
# remove(VS9_filtered_feature_bc_matrix)
# saveRDS(filtered_seurat, "data/naive_58_filtered_before_integration.RDS")
naive58 <- readRDS("data/naive_58_filtered_before_integration.RDS")

# Split the merge seurat object
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

# SCTransform on each seurat object
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Integration 3 samples
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# saveRDS(seurat_integrated, "data/naive_integrated_seurat.RDS")
seurat_integrated_naive <- readRDS("data/naive_integrated_seurat.RDS")

# Run PCA
seurat_integrated_naive <- RunPCA(object = seurat_integrated_naive)
# UMAP
seurat_integrated_naive <- RunUMAP(seurat_integrated_naive, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

# Cluster
# Determine the K-nearest neighbor graph
seurat_integrated_naive <- FindNeighbors(object = seurat_integrated_naive, dims = 1:30)
# Determine the clusters for various resolutions                                
seurat_integrated_naive <- FindClusters(object = seurat_integrated_naive, resolution = c(0,0.1,0.4,0.6,1.0))

# Set the clustering resolution to 0.4
required_reso <- "integrated_snn_res.1"
Idents(object = seurat_integrated_naive) <- "integrated_snn_res.1"

# cluster_20 <- FindMarkers(seurat_integrated_naive, ident.1 = 20,verbose = FALSE)



# Feature plot
metrics <-  c("nUMI", "nGene", "mitoRatio")
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

new.cluster.ids <- c("Schwann_cell","Macrophage","Schwann_cell","Schwann_cell","Macrophage","Schwann_cell","Schwann_cell",
                     "Macrophage","Macrophage","Macrophage","Schwann_cell","Fibroblast","Schwann_cell","Schwann_cell",
                     "T_cell","Macrophage","Macrophage","Schwann_cell","Schwann_cell","Macrophage","Myocyte","Endothelial_cell",
                     "Schwann_cell","Proliferating_cell","Macrophage","T_cell")
# annotation <- new.cluster.ids
# cluster <- c(0:(length(annotation)-1))
# 
# annos <- cbind(cluster, annotation) %>% as.data.frame()
# write.xlsx(annos,"results/reso_1.0/VS5&8_cluster_annotation.xlsx", row.names = F)


names(new.cluster.ids) <- levels(seurat_integrated_naive)
seurat_integrated_naive <- RenameIdents(seurat_integrated_naive, new.cluster.ids)
seurat_integrated_naive@meta.data$anno <- Idents(seurat_integrated_naive)
DimPlot(seurat_integrated_naive)
# meta_data <- seurat_integrated_naive@meta.data
# meta_data$anno <- Idents(seurat_integrated_naive)

p1 <- DimPlot(seurat_integrated_naive, label = F, group.by = "sample", pt.size = 0.2)
p1 <- p1[[1]]
p2 <- DimPlot(seurat_integrated_naive, reduction = "umap", label = T, pt.size = 0.2,label.size = 6)
p2 <- p2[[1]]
p3 <- FeaturePlot(seurat_integrated_naive, 
                  reduction = "umap", 
                  features = "nUMI",
                  pt.size = 0.4, 
                  order = T,
                  min.cutoff = 'q10',
                  label = F)
p3 <- p3[[1]]

p1 <- p1+ggtitle("Patient") + border() +
  theme(legend.position = c(0.8, 0.3),
        plot.title = element_text(size=30, vjust = -7,hjust = 0.95),
        legend.text = element_text(size=20,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1

p2 <- p2+ggtitle("Cell type") + border() + 
  theme(plot.title = element_text(size=30, vjust = -7,hjust = 0.95),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p2

p3 <- p3+ggtitle("Transcript count") + border() + 
  theme(legend.position = c(0.75, 0.3),
        plot.title = element_text(size=30, vjust = -7,hjust = 0.95),
        legend.text = element_text(size=20,face="bold"))

p3

# Marker gene plot
DefaultAssay(seurat_integrated_naive) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated_naive <- NormalizeData(seurat_integrated_naive, verbose = FALSE)

cilia_plot_vln_naive <- VlnPlot(seurat_integrated_naive, features = c("TTC25","ARL13B","LARP6","RFX7","RFX4","FOXJ1"), pt.size = 0, combine = F)
naive_ttc25 <- cilia_plot_vln_naive[[1]]+ggtitle("VS5&VS8: TTC25")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_arl13b <- cilia_plot_vln_naive[[2]]+ggtitle("VS5&VS8: ARL13B")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_larp6 <- cilia_plot_vln_naive[[3]]+ggtitle("VS5&VS8: LARP6")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_rfx7 <- cilia_plot_vln_naive[[4]]+ggtitle("VS5&VS8: RFX7")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_rfx4 <- cilia_plot_vln_naive[[5]]+ggtitle("VS5&VS8: RFX4")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_foxj1 <- cilia_plot_vln_naive[[6]]+ggtitle("VS5&VS8: FOXJ1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

fibrosis_plot_vln_naive <- VlnPlot(seurat_integrated_naive, features = c("STRAP","COL1A1","COL1A2","FKBP3","DHX9"), pt.size = 0, combine = F)
naive_strap <- fibrosis_plot_vln_naive[[1]]+ggtitle("VS5&VS8: STRAP")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_col1a1 <- fibrosis_plot_vln_naive[[2]]+ggtitle("VS5&VS8: COL1A1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_col1a2 <- fibrosis_plot_vln_naive[[3]]+ggtitle("VS5&VS86: COL1A2")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_fkbp3 <- fibrosis_plot_vln_naive[[4]]+ggtitle("VS5&VS8: FKBP3")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_dhx9 <- fibrosis_plot_vln_naive[[5]]+ggtitle("VS5&VS8: DHX9")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

# Cilia gene
(cystic_ttc25/cystic_arl13b/cystic_larp6/cystic_rfx7/cystic_rfx4/cystic_foxj1)|
  (naive_ttc25/naive_arl13b/naive_larp6/naive_rfx7/naive_rfx4/naive_foxj1)
# Fibrosis gene
(cystic_strap/cystic_col1a1/cystic_col1a2/cystic_fkbp3/cystic_larp6/cystic_dhx9)|
  (naive_strap/naive_col1a1/naive_col1a2/naive_fkbp3/naive_larp6/naive_dhx9)

# cystic_469_cilium_gene | 
#   (naive_ttc25/naive_arl13b/naive_larp6/naive_rfx7/naive_rfx4/naive_foxj1)
# 
# cystic_469_fibrosis_gene | 
#   (naive_strap/naive_col1a1/naive_col1a2/naive_fkbp3/naive_larp6/naive_dhx9)

naive_gene <- FeaturePlot(seurat_integrated_naive, 
                           reduction = "umap", 
                           features = c("CD74","OOEP",
                                        "NRG1","SEMA5A",
                                        "ADGRB3","PDE4D",
                                        "C1R","CLU",
                                        "LINGO1","LGI4"), 
                           min.cutoff = 'q8',
                           cols = c("azure3", "firebrick2"),
                           label = F,
                           pt.size = 0.4, combine = F)
naive_gene <- lapply(naive_gene, function(x){x + theme(legend.position = "none",
                                                         axis.text = element_blank(),
                                                         axis.ticks = element_blank(),
                                                         axis.title = element_blank())})
wrap_plots(naive_gene, ncol = 2)

naive_schwannoma_marker_gene <- VlnPlot(seurat_integrated_naive, features = c("CD74","OOEP","NRG1","SEMA5A","ADGRB3","PDE4D","C1R","CLU","LINGO1","LGI4"), pt.size = 0, combine = F)
saveRDS(naive_schwannoma_marker_gene,"results/reso_1.0/naive_schwannoma_marker_gene.RDS")
naive_schwannoma_marker_gene <- readRDS("results/reso_1.0/naive_schwannoma_marker_gene.RDS")
naive_CD74 <- naive_schwannoma_marker_gene[[1]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_OOEP <- naive_schwannoma_marker_gene[[2]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_NRG1 <- naive_schwannoma_marker_gene[[3]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_SEMA5A <- naive_schwannoma_marker_gene[[4]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_ADGRB3 <- naive_schwannoma_marker_gene[[5]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_PDE4D <- naive_schwannoma_marker_gene[[6]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_C1R <- naive_schwannoma_marker_gene[[7]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_CLU <- naive_schwannoma_marker_gene[[8]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_LINGO1 <- naive_schwannoma_marker_gene[[9]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_LGI4 <- naive_schwannoma_marker_gene[[10]]+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

(naive_CD74|naive_OOEP)/(naive_NRG1|naive_SEMA5A)/(naive_ADGRB3|naive_PDE4D)/(naive_C1R|naive_CLU)/(naive_LINGO1|naive_LGI4)



marker_plot <- FeaturePlot(seurat_integrated_naive, 
                           reduction = "umap", 
                           features = c("CLDN5","COL1A1","LYZ","CD3D","MPZ","MKI67"), 
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
p_Tcell <- marker_plot[[4]]+ggtitle("NK: CD3D")+theme(plot.title = element_text(vjust = -10,hjust = 0.95))
p_sch <- marker_plot[[5]]+ggtitle("Schwann-cell: MPZ")+theme(plot.title = element_text(vjust = -10,hjust = 0.95))
p_pro <- marker_plot[[6]]+ggtitle("Proliferating-cell: MKI67")+theme(plot.title = element_text(vjust = -10,hjust = 0.95))
# p_endo
# p_fibro
# p_macro
# p_Tcell
# p_sch
# p_pro

# FeaturePlot(seurat_integrated_naive, features = "MKI67",min.cutoff = 'q10',label = F)

sig_gene_plot <- FeaturePlot(seurat_integrated_naive, 
                             reduction = "umap", 
                             features = c("PTPRC","CD68","AIF1","PTGS2","TGFB1","IL1B",
                                          "IL6","TNF","VEGFA","ICAM1","NF2","APOD"), 
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
p_nf2 <- sig_gene_list[[11]]
p_apod <- sig_gene_list[[12]]

p_mid_top <- wrap_plots(p_endo,p_fibro,p_macro,p_Tcell,p_sch,p_pro,ncol = 3,nrow = 2)
p_mid_mid <- wrap_plots(p_ptprc,p_cd68,p_aif1,p_ptgs2,p_tgfb1,p_il1b,ncol = 3,nrow = 2)
p_mid_bot <- wrap_plots(p_il6,p_tnf,p_vegfa,p_icam1,p_nf2,p_apod,ncol = 3,nrow = 2)

p_mid_top/p_mid_mid/p_mid_bot
# p_mid_top
# p_mid_mid
# p_mid_bot

fill <- c("slateblue3","darkorange1")

macro_sub <- c(1,4,7,8,9,15,16,19,24)
sch_sub <- c(0,2,3,5,6,10,12,13,17,18,22)
t_sub <- c(14,25)

meta_data <- seurat_integrated_naive@meta.data
bar_df <- meta_data[,c("sample","integrated_snn_res.1","anno")]
macro_df <- bar_df[which(bar_df$anno=="Macrophage"),]
sample <- c(rep("VS5",9),rep("VS8",9))
cluster <- rep(macro_sub,2)
df <- cbind(sample,cluster) %>% as.data.frame()
df$percent <- NA
df$sum <- NA

value <- vector()
for (s in c("VS5","VS8")) {
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


sch_df <- bar_df[which(bar_df$anno=="Schwann_cell"),]
sample_sch <- c(rep("VS5",11),rep("VS8",11))
cluster_sch <- rep(sch_sub,2)
df_sch <- cbind(sample_sch,cluster_sch) %>% as.data.frame()
df_sch$percent <- NA
df_sch$sum <- NA

value_sch <- vector()
for (s in c("VS5","VS8")) {
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
df_sch$cluster_subtype <- as.factor(rep(1:11,2))
df_sch$anno <- rep("Schwann_cell",nrow(df_sch))


sch_t <- bar_df[which(bar_df$anno=="T_cell"),]
sample_t <- c(rep("VS5",2),rep("VS8",2))
cluster_t <- rep(t_sub,2)
df_t <- cbind(sample_t,cluster_t) %>% as.data.frame()
df_t$percent <- NA
df_t$sum <- NA

value_t <- vector()
for (s in c("VS5","VS8")) {
  for (c in t_sub) {
    num <- nrow(sch_t[sch_t$sample == s & sch_t$integrated_snn_res.1 == c,])
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

subtype_bar <- data.frame(cluster_subtype = c(1:9,1:11,1:2),
                          count = c(cluster_sum_macro,cluster_sum_sch,cluster_sum_t),
                          anno = c(rep("Macrophage",9),rep("Schwann_cell",11),rep("T_cell",2)))

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
                                         ifelse(box_macro$integrated_snn_res.1=="4","2",
                                          ifelse(box_macro$integrated_snn_res.1=="7","3",
                                          ifelse(box_macro$integrated_snn_res.1=="8","4",
                                          ifelse(box_macro$integrated_snn_res.1=="9","5",
                                          ifelse(box_macro$integrated_snn_res.1=="15","6",
                                          ifelse(box_macro$integrated_snn_res.1=="16","7",
                                          ifelse(box_macro$integrated_snn_res.1=="19","8","9"))))))))
                                         
box_sch <- meta_data[which(meta_data$anno=="Schwann_cell"),c("nGene","integrated_snn_res.1","anno")]
box_sch$integrated_snn_res.1 <- ifelse(box_sch$integrated_snn_res.1=="0","1",
                                       ifelse(box_sch$integrated_snn_res.1=="2","2",
                                        ifelse(box_sch$integrated_snn_res.1=="3","3",
                                        ifelse(box_sch$integrated_snn_res.1=="5","4",
                                        ifelse(box_sch$integrated_snn_res.1=="6","5",
                                        ifelse(box_sch$integrated_snn_res.1=="10","6",
                                        ifelse(box_sch$integrated_snn_res.1=="12","7",
                                        ifelse(box_sch$integrated_snn_res.1=="13","8",
                                        ifelse(box_sch$integrated_snn_res.1=="17","9",
                                        ifelse(box_sch$integrated_snn_res.1=="18","10","11"))))))))))

box_t <- meta_data[which(meta_data$anno=="T_cell"),c("nGene","integrated_snn_res.1","anno")]
box_t$integrated_snn_res.1 <- ifelse(box_t$integrated_snn_res.1=="14","1","2")

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

(p1/p2/p3)|(p_mid_top/p_mid_mid/p_mid_bot)|(p9/p10/p11)+plot_layout(widths = c(2,3,3,2))
# saveRDS(p1,"results/naive_merge/N_p1.RDS")
# saveRDS(p2,"results/naive_merge/N_p2.RDS")
# saveRDS(p3,"results/naive_merge/N_p3.RDS")
# saveRDS(p_mid_top,"results/naive_merge/N_p_mid_top.RDS")
# saveRDS(p_mid_mid,"results/naive_merge/N_p_mid_mid.RDS")
# saveRDS(p_mid_bot,"results/naive_merge/N_p_mid_bot.RDS")
# saveRDS(p9,"results/naive_merge/N_p9.RDS")
# saveRDS(p10,"results/naive_merge/N_p10.RDS")
# saveRDS(p11,"results/naive_merge/N_p11.RDS")


# Plot PCA
PCAPlot(seurat_integrated, split.by = "sample")
# UMAP Plot
DimPlot(seurat_integrated_naive, reduction = "umap", label = TRUE, label.size = 6)
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated_naive, label = TRUE, split.by = "sample")  + NoLegend()
# Clustree plot
# meta_data <- seurat_integrated@meta.data
# meta_data <- meta_data[,-c(17,14)]
# seurat_integrated@meta.data <- meta_data
clustree(seurat_integrated_naive, prefix = "integrated_snn_res.")

# Check the numbers of cells in each cluster
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

n_cells <- n_cells %>% column_to_rownames(var="orig.ident") %>% t() %>% as.data.frame()
total_cell <- sum(n_cells)
cluster_cell_sum <- apply(n_cells, 2, sum)

n_cells$VS5_percent <- NA
n_cells$VS8_percent <- NA
#n_cells$VS9_percent <- NA
n_cells$cluster_percent <- NA
n_cells$VS5_percent <- n_cells$VS5/cluster_cell_sum[1]
n_cells$VS8_percent <- n_cells$VS8/cluster_cell_sum[2]
#n_cells$VS9_percent <- n_cells$VS9/cluster_cell_sum[3]
n_cells$cluster_percent <- (n_cells$VS5+n_cells$VS8)/total_cell
cluster_cell_sum <- apply(n_cells, 2, sum)
n_cells <- rbind(n_cells, cluster_cell_sum)
# write.xlsx(n_cells, "results/naive_merge/cell_number_summary.xlsx")


# Preliminary SingleR annotation
# Create a SingleCellExperiment version of seurat object
DefaultAssay(seurat_integrated_naive) <- "RNA"
naive_sc <- as.SingleCellExperiment(seurat_integrated_naive)
colnames(naive_sc) <- colData(naive_sc)$cells
# Get all cell id in a cluster
cluster_barcode <- lapply(levels(seurat_integrated_naive@meta.data$integrated_snn_res.1), function(i) seurat_integrated_naive@meta.data[which(seurat_integrated_naive@meta.data$integrated_snn_res.1 == i),"cells"])

# cell level annotation fine labels
pred.cell_fine <- lapply(cluster_barcode, function(i) naive_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.fine))
names(pred.cell_fine) <- c(0:(length(cluster_barcode)-1))
fine_labels <- sapply(pred.cell_fine, function(i) table(i$labels) %>% as.data.frame() %>% 
                        arrange(desc(Freq)) %>% pull(1) %>% as.character() %>% head(3)) %>% as.data.frame()
fine_labels <- t(fine_labels) %>% as.data.frame()
row.names(fine_labels) <- c(0:(length(cluster_barcode)-1))
# saveRDS(pred.cell_fine, "results/reso_1.0/naive_pred.cell_fine_anno.RDS")
# pred.cell_fine <- readRDS("results/reso_1.0/naive_pred.cell_fine_anno.RDS")
# Cluster level single R annotation
pred.cluster_fine <- SingleR(test = naive_sc, ref = hpca.se, labels = hpca.se$label.fine, method = "cluster", clusters = meta_data$integrated_snn_res.1)
cluster_anno <- pred.cluster_fine$labels
names(cluster_anno) <- paste0("cluster",c(0:(length(cluster_anno)-1)))
cluster_anno <- as.data.frame(cluster_anno)
write.xlsx(cluster_anno,"results/reso_1.0/naive_cluster_level_singleR_anno_1.0.xlsx")


# Find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seurat_integrated_naive) <- "RNA"
all_markers_naive <- FindAllMarkers(object = seurat_integrated_naive, 
                              only.pos = TRUE,
                              logfc.threshold = 0.25)
# saveRDS(all_markers_naive,"results/reso_1.0/naive_58_all_markers.RDS")
top30_markers_naive <- all_markers_naive %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top30_markers_naive <- left_join(x = top30_markers_naive, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                           by = c("gene" = "gene_name"))
top30_markers_naive <- as.data.frame(top30_markers_naive)
write.xlsx(top30_markers_naive,"results/reso_1.0/naive_top30_marker.xlsx")

naive_cluster7.markers <- FindMarkers(seurat_integrated_naive, ident.1 = 7, min.pct = 0.25, only.pos = T)
top50_naive_cluster7.markers <- naive_cluster7.markers %>% rownames_to_column(var = "gene") %>% top_n(n = 50, wt = avg_logFC)
write.xlsx(top50_naive_cluster7.markers,"results/naive_merge/cluster7_top50_genes.xlsx")

naive_cluster7.markers_compare <- FindMarkers(seurat_integrated_naive, ident.1 = 7, ident.2 = c(0,3,4,5), min.pct = 0.25, only.pos = T)
top20_cluster7.markers_relative <- naive_cluster7.markers_compare %>% arrange(desc(avg_logFC)) %>% head(n=20) %>% rownames_to_column(var = "gene")
write.xlsx(top20_cluster7.markers_relative,"results/naive_merge/cluster7_top20_relative_genes.xlsx")

FeaturePlot(seurat_integrated_naive, features = c("GAP43","CCL2"))
DimPlot(seurat_integrated_naive,label = T) + NoLegend()
VlnPlot(seurat_integrated_naive,features = "LARP6",pt.size = 0) + NoLegend()



# GO annotation
markers_by_cluster <- lapply(unique(top30_markers_naive$cluster), function(i) top30_markers_naive[top30_markers_naive$cluster == i,])
markers <- lapply(markers_by_cluster, function(i) i[[9]])
cluster_go_annotation <- lapply(markers, function(i) goana(i) %>% topGO(n=15))
cluster_go_annotation <- lapply(cluster_go_annotation, function(i) i[[1]])
names(cluster_go_annotation) <- c(0:(length(cluster_go_annotation)-1))
df_go_annotation <- as.data.frame(cluster_go_annotation)
colnames(df_go_annotation) <- paste0("cluster",c(0:(length(cluster_go_annotation)-1)))
# write.xlsx(df_go_annotation, "results/reso_1.0/naive_58_cluster_GO_pathway.xlsx", row.names = F)

# View the annotation
cluster_go_annotation$'0'
cluster_go_annotation$'1'
cluster_go_annotation$'2'
cluster_go_annotation$'3'
cluster_go_annotation$'11'

# Top 10 marker gene heatmap
meta_data <- seurat_integrated@meta.data
dense_table <- seurat_integrated[["SCT"]]@scale.data
Idents(object = seurat_integrated) <- required_reso

top10_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_markers <- top10_markers[which(top10_markers$gene %in% row.names(dense_table)),]

dense_table <- dense_table[top10_markers$gene,]

row_anno <- top10_markers[,c("gene","cluster")] %>% data.frame()
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
         main = paste0("naive integrated Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))


# Gene Feature Plot
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated_naive) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated_naive <- NormalizeData(seurat_integrated_naive, verbose = FALSE)

# CD70 EPCAM gene distribution
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD70","EPCAM"), 
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD68","CD163","TNF","NOD2"), 
            min.cutoff = 'q9', 
            label = T)

# M1?
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("FCGR1A","CD86"), 
            min.cutoff = 'q9', 
            label = T)
# M2?
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD163","FCE2"), 
            min.cutoff = 'q9', 
            label = T)

# Paper ref
# Endothelial
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CLDN5","FLT1","CDH5","RAMP2"), 
            min.cutoff = 'q9', 
            label = F)

# Epithelial (not exist)
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CAPS","TMEM190","PIFO","SNTN"), 
            min.cutoff = 'q9', 
            label = F)

# Fibroblast
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("COL1A1","COL1A2","DCN","C1R"), 
            min.cutoff = 'q9', 
            label = F)

# bcell (not exist)
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD79A","IKGC","IGLC3","IGHG3"), 
            min.cutoff = 'q9', 
            label = F)

# myeloid
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD68","LYZ","MARCO","FCGR3A"), 
            min.cutoff = 'q9', 
            label = F)

# tcell
FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD3D","TRBC1","TRBC2","TRAC"), 
            min.cutoff = 'q9', 
            label = F)

FeaturePlot(seurat_integrated_naive, 
            reduction = "umap", 
            features = c("CD3D","GNLY"), 
            min.cutoff = 'q9', 
            label = F)

seurat_integrated_naive <- RenameIdents(seurat_integrated_naive, `8`="T-cell",`12`="Endothelial")


NK_naive<-VlnPlot(seurat_integrated_naive,
        features = c("ICAM1","GNLY","KLRC1"),
        pt.size = 0, combine = FALSE)
# saveRDS(NK_naive,"results/reso_1.0/NK_naive_vln.RDS")
NK_cystic<-VlnPlot(seurat_integrated,
        features = c("ICAM1","GNLY","KLRC1"),
        pt.size = 0, combine = FALSE)

naive_ICAM1 <- NK_naive[[1]]+ggtitle("Naive: ICAM1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_GNLY <- NK_naive[[2]]+ggtitle("Naive: GNLY")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
naive_KLRC1 <- NK_naive[[3]]+ggtitle("Naive: KLRC1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_ICAM1 <- NK_cystic[[1]]+ggtitle("Cystic: ICAM1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_GNLY <- NK_cystic[[2]]+ggtitle("Cystic: GNLY")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_KLRC1 <- NK_cystic[[3]]+ggtitle("Cystic: KLRC1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
(naive_ICAM1/naive_GNLY/naive_KLRC1)|(cystic_ICAM1/cystic_GNLY/cystic_KLRC1)




# -------------------------------------------------------------------------
VS5_df <- seurat_integrated_naive@meta.data[which(seurat_integrated_naive@meta.data$sample=="VS5"),c("sample","integrated_snn_res.1","anno")]
VS8_df <- seurat_integrated_naive@meta.data[which(seurat_integrated_naive@meta.data$sample=="VS8"),c("sample","integrated_snn_res.1","anno")]

cell_type <- c("Macrophage","Schwann_cell","Fibroblast","Proliferating_cell","Endothelial_cell","T_cell","NK","DC","Myocyte")


VS5_cell_number <- vector("integer")
for (x in cell_type) {
  num <- nrow(VS5_df[VS5_df$anno==x,])
  VS5_cell_number <- c(VS5_cell_number,num)
}
VS8_cell_number <- vector("integer")
for (x in cell_type) {
  num <- nrow(VS8_df[VS8_df$anno==x,])
  VS8_cell_number <- c(VS8_cell_number,num)
}

VS5_percent = VS5_cell_number/nrow(VS5_df)
VS8_percent = VS8_cell_number/nrow(VS8_df)


new_df_2 <- data.frame(cell_type = rep(cell_type,2),
                     cell_count = c(VS5_cell_number,VS8_cell_number),
                     percentage = c(VS5_percent,VS8_percent),
                     group = rep("naive",18))
# saveRDS(new_df_2,"results/reso_1.0/naive_58_cell_number_new_df_2.RDS")
# new_df_2 <- readRDS("results/reso_1.0/naive_58_cell_number_new_df_2.RDS")


total_bar_df <- rbind(new_df,new_df_2)
total_bar_df$cell_type <- factor(total_bar_df$cell_type)


# VS4&VS6 VS VS5&VS8 comparison
# total_bar_df[total_bar_df=="cystic"] <- "VS4&VS6&VS9"
total_bar_df[total_bar_df=="cystic"] <- "VS4&VS6"
total_bar_df[total_bar_df=="naive"] <- "VS5&VS8"


p_20 <- ggplot(total_bar_df, aes(cell_type, percentage, fill = group)) + scale_fill_manual(values=c('#999999','#E69F00')) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
  theme(axis.text.x = element_text(size=10)) +
  stat_compare_means(aes(group = group), method = "t.test",label = "p.format", label.y=0.80)


p_20



