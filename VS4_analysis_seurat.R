library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(AnnotationHub)
library(limma)

setwd("/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat")


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

## Read in 10x outputs direcotry and make seurat object (barcodes, features, matrix)
VS4.data <- Read10X(data.dir = "/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat/VS4_filtered_feature_bc_matrix")
VS4 <- CreateSeuratObject(counts = VS4.data, project = "VS4", min.cells = 3, min.features = 200)

# Find mitochondrial genes
VS4[["percent.mt"]] <- PercentageFeatureSet(VS4, pattern = "^MT-")

#summary(VS4@meta.data$nFeature_RNA)

#VlnPlot(VS4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Remove outliers in the dataset
VS4 <- subset(VS4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

## Normalize
VS4 <- NormalizeData(VS4)
VS4 <- FindVariableFeatures(VS4, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(VS4), 10)

#plot1 <- VariableFeaturePlot(VS4)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

## Scale data
all.genes <- rownames(VS4)
VS4 <- ScaleData(VS4, features = all.genes)

## Find PCAs
VS4 <- RunPCA(VS4, features = VariableFeatures(object = VS4))

#DimHeatmap(VS4, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(VS4, dims = 1:15, cells = 500, balanced = TRUE)

VS4 <- FindNeighbors(VS4, dims = 1:10)
VS4 <- FindClusters(VS4, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

## Check VS4 object information
head(VS4@meta.data,10)

## Set VS4 to resolution of 0.4 (lower resolution)
Idents(object = VS4) <- "RNA_snn_res.0.4"
VS4 <- RunUMAP(VS4, dims = 1:10)

## Save object
#saveRDS(VS4, file = "/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat/VS4.rds")

## Read object
#VS4 <- readRDS("/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat/VS4.rds")

## Plot the UMAP clustering
DimPlot(VS4, reduction = "umap", label = T)

#cluster1.markers <- FindMarkers(VS4, ident.1 = 0, min.pct = 0.25)
#head(cluster1.markers, n = 20)

## Integret annotations to the marker list
VS4.markers <- FindAllMarkers(VS4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VS4.markers <- VS4.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
VS4.markers <- left_join(x = VS4.markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                         by = c("gene" = "gene_name"))
## Save marker genes of each cluster
#write.csv(VS4.markers,"/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat/VS4_marker_genes.csv",row.names = F)

# Find GO pathway enrichment for each cluster (top 15 entries)
markers_by_cluster <- lapply(unique(VS4.markers$cluster), function(i) VS4.markers[VS4.markers$cluster == i,])
markers <- lapply(markers_by_cluster, function(i) i[[9]])
cluster_go_annotation <- lapply(markers, function(i) goana(i) %>% topGO(n=15))
cluster_go_annotation <- lapply(cluster_go_annotation, function(i) i[[1]])

## Find marker gene distribution
# Injured schwann markers
FeaturePlot(VS4, features = c("SHH", "BDNF", "SEMA4F", 
                              "CCL2", "CCL3", "CLCF1", 
                              "CRLF1", "TGFB1", "UCN2"), label = T)

# mesenchymal cells
FeaturePlot(VS4, features = c("ANGPT1", "CRLF1", 
                              "CCL2", "CCL3", "SEMA4F", 
                              "INHBB", "LIF", "NGF"), label = T)
# Myofibro
FeaturePlot(VS4, features = c("ACTA2", "MCAM", "PDGFA",
                              "MYLK", "MYL9", "IL6"), label = T)
# CAF
FeaturePlot(VS4, features = c("FAP", "THY1", "PDPN",
                              "MMP2", "PDGFRA", "TGFB3"), label = T)
# Fibro
FeaturePlot(VS4, features = c("FAP", "COL1A2", "PDPN",
                              "DCN", "COL3A1", "COL6A1"), label = T)
# DC
FeaturePlot(VS4, features = c("CD83","CD80","CD40","CCR7"), label = T)

# Endo
FeaturePlot(VS4, features = c("ENG","VWF","PECAM1"), label = T)

# Macrophage
FeaturePlot(VS4, features = c("CD14","CD163","CD68",
                              "FCGR2A","CSF1R"), label = T)
# Tcell
FeaturePlot(VS4, features = c("CD2","CD3D","CD3E","CD3G"), label = T)
# Neuron
FeaturePlot(VS4, features = c("NRXN1","NLGN1","NLGN2","NLGN3","LAMB1"), label = T)


# Fibroblast
VlnPlot(VS4, features = c("ACTA2","FAP","PDGFRB","SPARC"))
# Schwann cells
VlnPlot(VS4, features = c("SOX10","S100A1","MPZ","ERBB2"))
# M2 macrophage
VlnPlot(VS4, features = c("ADGRE1","CD80","TGFB1","LGMN","CD68","CD163"))
## DC
VlnPlot(VS4, features = c("CD83","CD80","CD40","CCR7"))
## LAMB1 gene
VlnPlot(VS4, features = c("PPP1R14A","ALDH1A1","BMI1"))

# Find differentailly expressed genes between clusters
cluster1.markers <- FindMarkers(VS4, ident.1 = 1, ident.2 = 6, min.pct = 0.25) %>% rownames_to_column(var = "gene")
cluster1.markers <- left_join(x = cluster1.markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                         by = c("gene" = "gene_name"))
cluster1vs6_anno <- head(cluster1.markers$entrezid,30) %>% goana() %>% topGO(n=15)

cluster1vs9.markers <- FindMarkers(VS4, ident.1 = 1, ident.2 = 9, min.pct = 0.25) %>% rownames_to_column(var = "gene")
cluster1vs9.markers <- left_join(x = cluster1vs9.markers, y = unique(annotations[, c("gene_name", "description","entrezid")]),
                              by = c("gene" = "gene_name"))
cluster1vs9_anno <- head(cluster1vs9.markers$entrezid,30) %>% goana() %>% topGO(n=15)

top10 <- VS4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(VS4, features = top10$gene) + NoLegend()

# Find all genes expressed in a cluster
gene_cluster <- AverageExpression(object = VS4)
gene_cluster <- rownames_to_column(gene_cluster, var = "name") %>% as.data.frame()
## Export the table
# write.csv(gene_cluster,"/Users/lufeixiaoxin/Downloads/seurat_cluster_genes.csv")





