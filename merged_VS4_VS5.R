library(dplyr)
library(Seurat)
library(SeuratData)
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

setwd("/Users/ericpan/Downloads/usc/kgp/seurat")

## Load in cell-ranger outputs
## -------------------------------------------------------------------------------------
VS4.data <- Read10X(data.dir = "/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat/VS4_filtered_feature_bc_matrix")
VS4 <- CreateSeuratObject(counts = VS4.data, project = "VS4", min.features = 100)

VS5.data <- Read10X(data.dir = "/Users/lufeixiaoxin/Downloads/usc_folder/kgp/seurat/VS5_filtered_feature_bc_matrix")
VS5 <- CreateSeuratObject(counts = VS5.data, project = "VS5", min.features = 100)
## -------------------------------------------------------------------------------------

## Process the data
## -------------------------------------------------------------------------------------
merged_seurat <- merge(x = VS4, 
                       y = VS5, 
                       add.cell.id = c("VS4", "VS5"))

## Remove redundant data
#remove(VS4.data)
#remove(VS5.data)
#remove(VS4)
#remove(VS5)

head(merged_seurat@meta.data, 10)
tail(merged_seurat@meta.data, 10)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(sample = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="merged_filtered_seurat.RData")
## -------------------------------------------------------------------------------------


## Do data QC on unfiltered data
## -------------------------------------------------------------------------------------

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
## -------------------------------------------------------------------------------------

## Filter raw data
## -------------------------------------------------------------------------------------
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                                  (nGene >= 250) & 
                                  (log10GenesPerUMI > 0.80) & 
                                  (mitoRatio < 0.20))

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
## -------------------------------------------------------------------------------------

## Check data variance on cell cycle
## -------------------------------------------------------------------------------------
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
## -------------------------------------------------------------------------------------

## Remove used data
#remove(seurat_phase)
#remove(counts)
#remove(nonzero)
#remove(filtered_counts)
#remove(merged_seurat)


## SCTransform
## -------------------------------------------------------------------------------------
options(future.globals.maxSize = 4000 * 1024^5)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("VS4", "VS5")]

for (i in 1:length(split_seurat)) {
        split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
        split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
        split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

## Integrate types
## -------------------------------------------------------------------------------------
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")









