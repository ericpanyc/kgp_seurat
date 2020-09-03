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

setwd("/Users/ericpan/Documents/pan/usc/kgp/seurat")
use_python("/usr/local/bin/python3")
hpca.se <- HumanPrimaryCellAtlasData()

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


object_list <- c("VS4","VS6","VS9")
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
cystic <- merge(VS4_filtered_feature_bc_matrix, y = c(VS6_filtered_feature_bc_matrix, VS9_filtered_feature_bc_matrix), add.cell.ids = c("VS4", "VS6", "VS9"), project = "cystic")

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
# saveRDS(filtered_seurat, "data/cystic_filtered_before_normalize_seurat.RDS")

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

# saveRDS(seurat_integrated, "data/cystic_integrated_seurat.RDS")
seurat_integrated <- readRDS("data/cystic_integrated_seurat.RDS")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
# UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

# Cluster
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:30)
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(0,0.1,0.4,0.6,1.0))

# Set the clustering resolution to 1
required_reso <- "integrated_snn_res.1"
Idents(object = seurat_integrated) <- "integrated_snn_res.1"

meta_data <- seurat_integrated@meta.data

# mitoRatio plot
metrics <-  c("nUMI", "nGene", "mitoRatio")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Plot PCA
PCAPlot(seurat_integrated, split.by = "sample")
# UMAP Plot
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 3)
# UMAP group by sample
DimPlot(seurat_integrated, label = F, group.by = "sample", pt.size = 0.2)
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, label = TRUE)  + NoLegend()

# Clustree plot
clustree(seurat_integrated, prefix = "integrated_snn_res.")

# Check the numbers of cells in each cluster
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

n_cells <- n_cells %>% column_to_rownames(var="orig.ident") %>% t() %>% as.data.frame()
total_cell <- sum(n_cells)
cluster_cell_sum <- apply(n_cells, 2, sum)

n_cells$VS4_percent <- NA
n_cells$VS6_percent <- NA
n_cells$VS9_percent <- NA
n_cells$cluster_percent <- NA
n_cells$VS4_percent <- n_cells$VS4/cluster_cell_sum[1]
n_cells$VS6_percent <- n_cells$VS6/cluster_cell_sum[2]
n_cells$VS9_percent <- n_cells$VS9/cluster_cell_sum[3]
n_cells$cluster_percent <- (n_cells$VS4+n_cells$VS6+n_cells$VS9)/total_cell
cluster_cell_sum <- apply(n_cells, 2, sum)
n_cells <- rbind(n_cells, cluster_cell_sum)
# write.xlsx(n_cells, "results/cystic_merge/cell_number_summary.xlsx")

# Cilium and Fibrosis gene violin plot
# Marker gene plot
DefaultAssay(seurat_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

cilia_plot_vln_cys <- VlnPlot(seurat_integrated, features = c("TTC25","ARL13B","LARP6","RFX7","RFX4","FOXJ1"), pt.size = 0, combine = F)
cystic_ttc25 <- cilia_plot_vln_cys[[1]]+ggtitle("VS4&VS6&VS9: TTC25")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_arl13b <- cilia_plot_vln_cys[[2]]+ggtitle("VS4&VS6&VS9: ARL13B")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_larp6 <- cilia_plot_vln_cys[[3]]+ggtitle("VS4&VS6&VS9: LARP6")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_rfx7 <- cilia_plot_vln_cys[[4]]+ggtitle("VS4&VS6&VS9: RFX7")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_rfx4 <- cilia_plot_vln_cys[[5]]+ggtitle("VS4&VS6&VS9: RFX4")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_foxj1 <- cilia_plot_vln_cys[[6]]+ggtitle("VS4&VS6&VS9: FOXJ1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

fibrosis_plot_vln_cys <- VlnPlot(seurat_integrated, features = c("STRAP","COL1A1","COL1A2","FKBP3","DHX9"), pt.size = 0, combine = F)
cystic_strap <- fibrosis_plot_vln_cys[[1]]+ggtitle("VS4&VS6&VS9: STRAP")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_col1a1 <- fibrosis_plot_vln_cys[[2]]+ggtitle("VS4&VS6&VS9: COL1A1")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_col1a2 <- fibrosis_plot_vln_cys[[3]]+ggtitle("VS4&VS6&VS9: COL1A2")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_fkbp3 <- fibrosis_plot_vln_cys[[4]]+ggtitle("VS4&VS6&VS9: FKBP3")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cystic_dhx9 <- fibrosis_plot_vln_cys[[5]]+ggtitle("VS4&VS6&VS9: DHX9")+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())

cystic_469_cilium_gene <- cystic_ttc25/cystic_arl13b/cystic_larp6/cystic_rfx7/cystic_rfx4/cystic_foxj1
cystic_469_fibrosis_gene <- cystic_strap/cystic_col1a1/cystic_col1a2/cystic_fkbp3/cystic_larp6/cystic_dhx9
# saveRDS(cystic_469_cilium_gene,"results/reso_1.0/cystic_469_cilium_gene_violin_plot.RDS")
# saveRDS(cystic_469_fibrosis_gene,"results/reso_1.0/cystic_469_fibrosis_gene_violin_plot.RDS")
# cystic_469_cilium_gene <- readRDS("results/reso_1.0/cystic_469_cilium_gene_violin_plot.RDS")
# cystic_469_fibrosis_gene <- readRDS("results/reso_1.0/cystic_469_fibrosis_gene_violin_plot.RDS")

# Screenshot gene
cystic_gene <- FeaturePlot(seurat_integrated, 
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
cystic_gene <- lapply(cystic_gene, function(x){x + theme(legend.position = "none",
                                                         axis.text = element_blank(),
                                                         axis.ticks = element_blank(),
                                                         axis.title = element_blank())})
wrap_plots(cystic_gene, ncol = 2)

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



# Preliminary SingleR annotation
# Create a SingleCellExperiment version of seurat object
# meta_data <- seurat_integrated@meta.data
DefaultAssay(seurat_integrated) <- "RNA"
cystic_sc <- as.SingleCellExperiment(seurat_integrated)
colnames(cystic_sc) <- colData(cystic_sc)$cells
# Get all cell id in a cluster
cluster_barcode <- lapply(levels(seurat_integrated@meta.data$integrated_snn_res.1), function(i) seurat_integrated@meta.data[which(seurat_integrated@meta.data$integrated_snn_res.1 == i),"cells"])
# cell level annotation fine labels
pred.cystic_cell_fine <- lapply(cluster_barcode, function(i) cystic_sc[,i] %>% SingleR(ref = hpca.se, labels = hpca.se$label.fine))
names(pred.cystic_cell_fine) <- c(0:(length(cluster_barcode)-1))
fine_labels <- sapply(pred.cystic_cell_fine, function(i) table(i$labels) %>% as.data.frame() %>% 
                        arrange(desc(Freq)) %>% pull(1) %>% as.character() %>% head(3)) %>% as.data.frame()
fine_labels <- t(fine_labels) %>% as.data.frame()
# saveRDS(pred.cystic_cell_fine,"results/reso_1.0/cystic_469_singleR_cell_level_result.RDS")

# cluster level annotation fine labels
pred.cluster_fine <- SingleR(test = cystic_sc, ref = hpca.se, labels = hpca.se$label.fine, method = "cluster", clusters = seurat_integrated@meta.data$integrated_snn_res.1)
cluster_anno <- pred.cluster_fine$labels
names(cluster_anno) <- paste0("cluster",c(0:(length(cluster_anno)-1)))
cluster_anno <- as.data.frame(cluster_anno)
# write.xlsx(cluster_anno,"results/reso_1.0/cystic_469_singleR_cluster_level_anno.xlsx")
# cell level annotation fine labels

## _____________________________________________________________________________________________________________________________
new.cluster.ids <- c("Schwann_cell","Macrophage","Macrophage","Macrophage","Macrophage",
                     "Schwann_cell","Macrophage","Macrophage","T_cell","Macrophage",
                     "DC","Macrophage","Schwann_cell","Macrophage","Fibroblast",
                     "Schwann_cell","Fibroblast","Schwann_cell","Macrophage","T_cell",
                     "Schwann_cell","Macrophage","NK","Macrophage",'Proliferating_cell',
                     "Endothelial_cell","Monocyte")
# annotation <- new.cluster.ids
# cluster <- c(0:26)
# 
# annos <- cbind(cluster, annotation) %>% as.data.frame()
# write.xlsx(annos,"results/reso_1.0/VS4&6&9_cluster_annotation.xlsx", row.names = F)

names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
seurat_integrated@meta.data$anno <- Idents(seurat_integrated)
DimPlot(seurat_integrated, reduction = "umap")
# meta_data$anno <- Idents(seurat_integrated)


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
  theme(legend.position = c(0.1, 0.15),
        plot.title = element_text(size=30, vjust = -7,hjust = 0.95),
        legend.text = element_text(size=20,face="bold"))

p3

p1/p2/p3

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

# p_endo


sig_gene_plot <- FeaturePlot(seurat_integrated, 
                             reduction = "umap", 
                             features = c("PTPRC","CD68","AIF1","PTGS2","TGFB1","IL1B",
                                          "IL6","TNF","VEGFA","ICAM1","NF2","MKI67"), 
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
p_mki67 <- sig_gene_list[[12]]


p_mid_top <- wrap_plots(p_endo,p_fibro,p_macro,p_Tcell,p_sch,p_nk,ncol = 3,nrow = 2)
p_mid_mid <- wrap_plots(p_ptprc,p_cd68,p_aif1,p_ptgs2,p_tgfb1,p_il1b,ncol = 3,nrow = 2)
p_mid_bot <- wrap_plots(p_il6,p_tnf,p_vegfa,p_icam1,p_nf2,p_mki67,ncol = 3,nrow = 2)

p_mid_top/p_mid_mid/p_mid_bot

#(p1|p2|p3)/(p_endo|p_fibro|p_macro|p_Tcell|p_sch|p_pro)/(p_ptprc|p_cd68|p_aif1|p_ptgs2|p_tgfb1|p_il1b)/(p_il6|p_tnf|p_vegfa|p_icam1|p_nf2|p_apod) +  plot_layout(heights = c(2,1,1,1))

# Plot stacked bar chart
meta_data <- seurat_integrated@meta.data
macro_sub <- c(1,2,3,4,6,7,9,11,13,18,21,23)
sch_sub <- c(0,5,12,15,17,20)
tcell_sub <- c(8,19)
fibro_sub <- c(14,16)
bar_df <- meta_data[,c("sample",required_reso,"anno")]
macro_df <- bar_df[which(bar_df$anno=="Macrophage"),]
sample <- c(rep("VS4",12),rep("VS6",12),rep("VS9",12))
cluster <- rep(macro_sub,3)
df <- cbind(sample,cluster) %>% as.data.frame()
df$percent <- NA
df$sum <- NA

value <- vector()
for (s in c("VS4","VS6","VS9")) {
  for (c in macro_sub) {
    num <- nrow(macro_df[macro_df$sample == s & macro_df$integrated_snn_res.1 == c,])
    value <- c(value, num)
  }
}
df$sum <- value

percent <- vector()
for (x in c(1:36)) {
  curr_cluster <- as.integer(df[x,"cluster"])
  val <- df[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent <- c(percent, val)
}

df$percent <- percent
df$cluster_subtype <- as.factor(rep(1:12,3))
df$anno <- rep("Macrophage",nrow(df))

fill <- c("hotpink4","limegreen","darkgoldenrod3","darkorchid2")

# Schwann cell count
sch_df <- bar_df[which(bar_df$anno=="Schwann_cell"),]
sample_sch <- c(rep("VS4",6),rep("VS6",6),rep("VS9",6))
cluster_sch <- rep(sch_sub,3)
df_sch <- cbind(sample_sch,cluster_sch) %>% as.data.frame()
df_sch$percent <- NA
df_sch$sum <- NA

value_sch <- vector()
for (s in c("VS4","VS6","VS9")) {
  for (c in sch_sub) {
    num <- nrow(sch_df[sch_df$sample == s & sch_df$integrated_snn_res.1 == c,])
    value_sch <- c(value_sch, num)
  }
}
df_sch$sum <- value_sch
percent_sch <- vector()
for (x in c(1:18)) {
  curr_cluster <- as.integer(df_sch[x,"cluster_sch"])
  val <- df_sch[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent_sch <- c(percent_sch, val)
}
df_sch$percent <- percent_sch
df_sch$cluster_subytpe <- as.factor(rep(1:6,3))
df_sch$anno <- rep("Schwann_cell",nrow(df_sch))

# Fibroblast cell count
fibro_df <- bar_df[which(bar_df$anno=="Fibroblast"),]
sample_fibro <- c(rep("VS4",2),rep("VS6",2),rep("VS9",2))
cluster_fibro <- rep(fibro_sub,3)
df_fibro <- cbind(sample_fibro,cluster_fibro) %>% as.data.frame()
df_fibro$percent <- NA
df_fibro$sum <- NA

value_fibro <- vector()
for (s in c("VS4","VS6","VS9")) {
  for (c in fibro_sub) {
    num <- nrow(fibro_df[fibro_df$sample == s & fibro_df$integrated_snn_res.1 == c,])
    value_fibro <- c(value_fibro, num)
  }
}
df_fibro$sum <- value_fibro

percent_fibro <- vector()
for (x in c(1:6)) {
  curr_cluster <- as.integer(df_fibro[x,"cluster_fibro"])
  val <- df_fibro[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent_fibro <- c(percent_fibro, val)
}
df_fibro$percent <- percent_fibro
df_fibro$cluster_subtype <- as.factor(rep(1:2,3))
df_fibro$anno <- rep("Fibroblast",nrow(df_fibro))

# T_cell count
tcell_df <- bar_df[which(bar_df$anno=="T_cell"),]
sample_tcell <- c(rep("VS4",2),rep("VS6",2),rep("VS9",2))
cluster_tcell <- rep(tcell_sub,3)
df_tcell <- cbind(sample_tcell,cluster_tcell) %>% as.data.frame()
df_tcell$percent <- NA
df+tcell$sum <- NA

value_tcell <- vector()
for (s in c("VS4","VS6","VS9")) {
  for (c in tcell_sub) {
    num <- nrow(tcell_df[tcell_df$sample == s & tcell_df$integrated_snn_res.1 == c,])
    value_tcell <- c(value_tcell, num)
  }
}
df_tcell$sum <- value_tcell

percent_tcell <- vector()
for (x in c(1:6)) {
  curr_cluster <- as.integer(df_tcell[x,"cluster_tcell"])
  val <- df_tcell[x,"sum"] / length(cluster_barcode[[curr_cluster+1]])
  percent_tcell <- c(percent_tcell, val)
}
df_tcell$percent <- percent_tcell
df_tcell$cluster_subtype <- as.factor(rep(1:2,3))
df_tcell$anno <- rep("T_cell",nrow(df_tcell))



# bar_fibro
names(df_sch) <- names(df)
names(df_fibro) <- names(df)
names(df_tcell) <- names(df)

merged_df <- rbind(df,df_fibro,df_sch,df_tcell)
#merged_df <- merged_df %>% dplyr::rename(cluster_subtype=cluster_in_chart)

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
cluster_sum_fibro <- sapply(cluster_barcode[(fibro_sub+1)], function(x) length(x))
cluster_sum_tcell <- sapply(cluster_barcode[(tcell_sub+1)], function(x) length(x))

subtype_bar <- data.frame(cluster_subtype = c(1:12,1:6,1:2,1:2),
                          count = c(cluster_sum_macro,cluster_sum_sch,cluster_sum_fibro,cluster_sum_tcell),
                          anno = c(rep("Macrophage",12),rep("Schwann_cell",6),rep("Fibroblast",2),rep("T_cell",2)))


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
                                           ifelse(box_macro$integrated_snn_res.1=="6","5",
                                           ifelse(box_macro$integrated_snn_res.1=="7","6",
                                           ifelse(box_macro$integrated_snn_res.1=="9","7",
                                           ifelse(box_macro$integrated_snn_res.1=="11","8",
                                           ifelse(box_macro$integrated_snn_res.1=="13","9",
                                           ifelse(box_macro$integrated_snn_res.1=="18","10",
                                           ifelse(box_macro$integrated_snn_res.1=="21","11","12")))))))))))

box_sch <- meta_data[which(meta_data$anno=="Schwann_cell"),c("nGene","integrated_snn_res.1","anno")]
box_sch$integrated_snn_res.1 <- ifelse(box_sch$integrated_snn_res.1=="0","1",
                                         ifelse(box_sch$integrated_snn_res.1=="5","2",
                                            ifelse(box_sch$integrated_snn_res.1=="12","3",
                                               ifelse(box_sch$integrated_snn_res.1=="15","4",
                                                  ifelse(box_sch$integrated_snn_res.1=="17","5","6")))))

box_fibro <- meta_data[which(meta_data$anno=="Fibroblast"),c("nGene","integrated_snn_res.1","anno")]
box_fibro$integrated_snn_res.1 <- ifelse(box_fibro$integrated_snn_res.1=="14","1","2")

box_tcell <- meta_data[which(meta_data$anno=="T_cell"),c("nGene","integrated_snn_res.1","anno")]
box_tcell$integrated_snn_res.1 <- ifelse(box_tcell$integrated_snn_res.1=="8","1","2")

box_merge <- rbind(box_macro,box_sch,box_fibro,box_tcell)
box_merge$nGene_log <- log10(box_merge$nGene)
box_merge <- dplyr::rename(box_merge, cluster_subtype = integrated_snn_res.1)
box_merge$anno_f = factor(box_merge$anno, levels=c('Fibroblast','Macrophage','Schwann_cell','T_cell'))
box_merge$cluster_subtype <- factor(box_merge$cluster_subtype, levels = c(1:12))
# box_merge <- box_merge %>%  dplyr::rename(cluster_subtype=cluster_id)
                                          

p11 <- ggplot(box_merge, aes(x = cluster_subtype, y =nGene_log)) + 
  facet_grid(. ~ anno_f,scales="free",space = "free",switch = "x") +
  geom_boxplot(fill='#A4A4A4', color="darkblue",outlier.size=0.5, outlier.colour="red") + 
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, angle = 90),
        strip.text.x = element_text(size=10,face="bold"),
        strip.background.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0.5, "cm"))

p11

p9/p10/p11

(p1/p2/p3)|(p_mid_top/p_mid_mid/p_mid_bot)|(p9/p10/p11)+plot_layout(widths = c(2,3,3,2))




# Find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seurat_integrated) <- "RNA"
all_markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)
top30_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
top30_markers <- left_join(x = top30_markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                         by = c("gene" = "gene_name"))
write.xlsx(as.data.frame(top30_markers), "results/reso_1.0/cystic_469_top30_markers.xlsx", row.names = F)

# Find conservative markers differentiate cluster 0 and 1.
cluster0vs1_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   ident.2 = 1,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25) %>% rownames_to_column(var = "gene") %>% left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),by = c("gene" = "gene_name"))
conservative_markeres <- cluster0vs1_conserved_markers[!is.infinite(cluster0vs1_conserved_markers$VS6_avg_logFC+
                                                                      cluster0vs1_conserved_markers$VS9_avg_logFC+
                                                                      cluster0vs1_conserved_markers$VS4_avg_logFC),]
top30 <- conservative_markeres %>% 
  mutate(avg_fc = (VS6_avg_logFC + VS9_avg_logFC + VS4_avg_logFC) /3) %>% 
  top_n(n = 30, wt = avg_fc)

diff_go <- goana(top30$entrezid) %>% topGO(n=15)
# # Find conservative markers
# DefaultAssay(seurat_integrated) <- "RNA"
# 
# get_conserved <- function(cluster){
#   FindConservedMarkers(seurat_integrated,
#                        ident.1 = cluster,
#                        grouping.var = "sample",
#                        only.pos = TRUE) %>%
#     rownames_to_column(var = "gene") %>%
#     left_join(y = unique(annotations[, c("gene_name", "description", "entrezid")]),
#               by = c("gene" = "gene_name")) %>%
#     cbind(cluster_id = cluster, .)
# }
# 
# # # Call the function to get the markers table
# conserved_markers <- map_dfr(c(0:14), get_conserved)

# GO annotation
markers_by_cluster <- lapply(unique(top30_markers$cluster), function(i) top30_markers[top30_markers$cluster == i,])
markers <- lapply(markers_by_cluster, function(i) i[[9]])
cluster_go_annotation <- lapply(markers, function(i) goana(i) %>% topGO(n=15))
cluster_go_annotation <- lapply(cluster_go_annotation, function(i) i[[1]])
names(cluster_go_annotation) <- c(0:14)
df_go_annotation <- as.data.frame(cluster_go_annotation)
colnames(df_go_annotation) <- paste0("cluster",c(0:14))
# write.xlsx(df_go_annotation, "results/cystic_merge/cluster_GO_pathway.xlsx", row.names = F)

# # View the annotation
# cluster_go_annotation$'0'
# cluster_go_annotation$'1'
# cluster_go_annotation$'2'
# cluster_go_annotation$'3'
# cluster_go_annotation$'11'

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
         main = paste0("Cystic integrated Cluster Marker Genes Heatmap SCT UMAP ", substring(required_reso,13,15)))

# Resolution 0.1 
Idents(seurat_integrated) <- "integrated_snn_res.0.1"
DimPlot(seurat_integrated)
macrophage <- FindMarkers(seurat_integrated,
                          ident.1 = 0,
                          ident.2 = 1)
macrophage <- macrophage %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))

macrophage <- macrophage[, c(1, 3:5,2,6:8)]
macrophage <- macrophage[order(-abs(macrophage$avg_logFC)),]


# Gene Feature Plot
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

Idents(seurat_integrated) <- "integrated_snn_res.0.4"

macrophage2 <- FindMarkers(seurat_integrated,
                          ident.1 = 0,
                          ident.2 = 3)
macrophage2 <- macrophage2 %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))

macrophage2 <- macrophage2[, c(1, 3:5,2,6:8)]

macrophage2 <- macrophage2 %>%
  dplyr::arrange(p_val_adj)


macrophage_8vs0_3 <- FindMarkers(seurat_integrated,
                           ident.1 = 8,
                           ident.2 = c(0,3))
macrophage_8vs0_3 <- macrophage_8vs0_3 %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))

macrophage_8vs0_3 <- macrophage_8vs0_3[, c(1, 3:5,2,6:8)]

macrophage_8vs0_3 <- macrophage_8vs0_3[order(-abs(macrophage_8vs0_3$avg_logFC)),]
  
macrophage_8vs0_3 <- FindMarkers(seurat_integrated,
                                 ident.1 = 8,
                                 ident.2 = c(0,3))


fibro_2vs5_9 <- FindMarkers(seurat_integrated,
                                 ident.1 = 2,
                                 ident.2 = c(5,9))
fibro_2vs5_9 <- fibro_2vs5_9 %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))
fibro_2vs5_9 <- fibro_2vs5_9[, c(1, 3:5,2,6:8)]
fibro_2vs5_9 <- fibro_2vs5_9[order(-abs(fibro_2vs5_9$avg_logFC)),]

fibro_5vs9 <- FindMarkers(seurat_integrated,
                            ident.1 = 5,
                            ident.2 = 9)
fibro_5vs9 <- fibro_5vs9 %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))
fibro_5vs9 <- fibro_5vs9[, c(1, 3:5,2,6:8)]
fibro_5vs9 <- fibro_5vs9[order(-abs(fibro_5vs9$avg_logFC)),]

fibro_10vs2_5_9 <- FindMarkers(seurat_integrated,
                          ident.1 = 10,
                          ident.2 = c(2,5,9))
fibro_10vs2_5_9 <- fibro_10vs2_5_9 %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))
fibro_10vs2_5_9 <- fibro_10vs2_5_9[, c(1, 3:5,2,6:8)]
fibro_10vs2_5_9 <- fibro_10vs2_5_9[order(-abs(fibro_10vs2_5_9$avg_logFC)),]

fibro_8vs0_3 <- FindMarkers(seurat_integrated,
                               ident.1 = 8,
                               ident.2 = c(0,3))
fibro_8vs0_3 <- fibro_8vs0_3 %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description","entrezid")]),
            by = c("gene" = "gene_name"))
fibro_8vs0_3 <- fibro_8vs0_3[, c(1, 3:5,2,6:8)]
fibro_8vs0_3 <- fibro_8vs0_3[order(-abs(fibro_8vs0_3$avg_logFC)),]

# CD70 EPCAM gene distribution
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD70","EPCAM"), 
            min.cutoff = 'q10', 
            label = TRUE)

# Non-myelinating schwann cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CDH2","L1CAM","EDNRB","EMP1"), 
            min.cutoff = 'q10', 
            label = TRUE)

# CCND1+ Schwann cells tumor cells?
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CCND1","APOE","RAB13","TMSB4X","MAP1B"), 
            min.cutoff = 'q10', 
            label = TRUE)

# myelinating Schwann
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MAG","MBP","PMP22","MPZ","PLP1"), 
            min.cutoff = 'q10', 
            label = TRUE)

# Smooth muscle and perycyte
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("DES","MYLK","ACTA2","RGS5"), 
            min.cutoff = 'q10', 
            label = TRUE)

# CAF
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("S100A4","CSPG4","ACTA2","FAP","TNC","POSTN","VIM","DES"), 
            min.cutoff = 'q10', 
            label = F)

# endoneurial mesenchymal cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("ETV1","SOX9","OSR2","MEOX1"), 
            min.cutoff = 'q10', 
            label = TRUE)

# proliferating cells – Mki67 and Top2a
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MKI67","TOP2A"), 
            min.cutoff = 'q10', 
            label = F)

# Cluster 10
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("NRXN1","EDNRB","PTN","CDH19"), 
            min.cutoff = 'q10', 
            label = TRUE)

# Macropahge anti-inflam
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD68","CD163","TNF"), 
            min.cutoff = 'q10', 
            label = F)

# test
# Cp, Il1rl1, Maoa, and Cdr2
# Clec7a, Spp1, Lpl, Axl, Ms4a4c, and Ms4a6c
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MRC1","CD68","CD163","ODC1","ITGB1","CEBPB"), 
            min.cutoff = 'q8', 
            label = F)

# Clec7a, Spp1, Lpl, Axl, Ms4a4c, and Ms4a6c
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("APOD","APOE","PLCG2","TREM2"), 
            min.cutoff = 'q9', 
            label = T)

# PCLG2 (hematopoietic cells)
# CXCL8 (M1 marker)
# CD83 DC - cluster 1 and 4
# JUN FOS M1 pathway -  cluster 3
# DCN collagen fibril assembly and tumor suppression
# APOD antioxidant and anti-inflammatory activity
# SPP1?
# MPZ nerve cells
# PDGFRA fibroblast mesenchymal
# PLEK ?
# EGR2 ?
# NCAM1 ?
# PPP1R14a schwann cells
# CLU PTX3 - cluster 2 sign of damage

# Paper ref
# Endothelial
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CLDN5","FLT1","CDH5","RAMP2"), 
            min.cutoff = 'q9', 
            label = F)

# Epithelial (not exist)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CAPS","TMEM190","PIFO","SNTN"), 
            min.cutoff = 'q9', 
            label = F)

# Fibroblast
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("COL1A1","COL1A2","DCN","C1R"), 
            min.cutoff = 'q9', 
            label = F)

# bcell (not exist)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD79A","IKGC","IGLC3","IGHG3"), 
            min.cutoff = 'q9', 
            label = F)

# myeloid
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD68","LYZ","MARCO","FCGR3A"), 
            min.cutoff = 'q9', 
            label = F)

# tcell
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD3D","GNLY","NKG7"), 
            min.cutoff = 'q9', 
            label = T)

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("APOD","PLCG2"), 
            min.cutoff = 'q9', 
            label = T)

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PTN","GAP43","WNT1","PMP22"), 
            min.cutoff = 'q9', 
            label = T)

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PECAM1","FLT1"), 
            min.cutoff = 'q9', 
            label = T)

VlnPlot(seurat_integrated,"TMEM119")

cluster10_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 10,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

DimPlot(seurat_integrated)

# ---------------------------------------------------------------
VS4_df <- seurat_integrated@meta.data[which(seurat_integrated@meta.data$sample=="VS4"),c("sample","integrated_snn_res.1","anno")]
VS6_df <- seurat_integrated@meta.data[which(seurat_integrated@meta.data$sample=="VS6"),c("sample","integrated_snn_res.1","anno")]
VS9_df <- seurat_integrated@meta.data[which(seurat_integrated@meta.data$sample=="VS9"),c("sample","integrated_snn_res.1","anno")]

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
VS9_cell_number <- vector("integer")
for (x in cell_type) {
  num <- nrow(VS9_df[VS9_df$anno==x,])
  VS9_cell_number <- c(VS9_cell_number,num)
}


VS4_percent = VS4_cell_number/nrow(VS4_df)
VS6_percent = VS6_cell_number/nrow(VS6_df)
VS9_percent = VS9_cell_number/nrow(VS9_df)


new_df <- data.frame(cell_type = rep(cell_type,3),
                     cell_count = c(VS4_cell_number,VS6_cell_number,VS9_cell_number),
                     percentage = c(VS4_percent,VS6_percent,VS9_percent),
                     group = rep("cystic",27))
# saveRDS(new_df,"results/reso_1.0/cytic_469_cell_percentage.RDS")
cystic_new_df <- readRDS("results/reso_1.0/cytic_469_cell_percentage.RDS")
