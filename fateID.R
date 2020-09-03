library(devtools)

install_github("dgrun/FateID")

library(FateID)

assay <- as.matrix(GetAssayData(seurat_integrated))
dim(assay)

anno <- seurat_integrated@meta.data$integrated_snn_res.0.4 %>% as.integer()
names(anno) <- row.names(seurat_integrated@meta.data) 
