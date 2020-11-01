library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(reticulate)

use_py("/home/yachenpan/miniconda3/envs/mypython3/bin/python3")
setwd("/home/yachenpan/single_cell/naive_group")

# Prepare file list
object_list <- c("VS5","VS8")
object_list <- paste0(object_list, "_filtered")

# Read in count matrix
for (file in object_list){
  seurat_data <- Read10X(data.dir = paste0("../data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# Merge two Seurat object
naive <- merge(x= VS5_filtered, y = VS8_filtered, add.cell.ids = c("VS5", "VS8"), project = "naive")

