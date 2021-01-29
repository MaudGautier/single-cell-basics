#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# functions.R                                                                  #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script contains the functions required to run the Seurat pipeline




# Generic -----------------------------------------------------------------

# Import libraries
import_libraries <- function() {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(patchwork))
  suppressPackageStartupMessages(library(SingleR))
}

# Create folder
init_folder <- function(folder) {
  if (!dir.exists(folder)) {
    dir.create(file.path(folder), showWarnings = FALSE, recursive = TRUE)
  }
}



# Quality check -----------------------------------------------------------

# Preprocesses the data: runs QC, normalisation and scaling
# Output: final object (with counts and data)


# Create seurat object and get mitochondrial percentage
get_seurat_object <- function(path_dataset, 
                              sample_name, 
                              min_cells = 3, 
                              min_features = 200) {
  
  # Create Seurat object
  sc.data <- Read10X(data.dir = path_dataset)
  sc.seurat <- CreateSeuratObject(counts = sc.data, 
                                  project = sample_name, 
                                  min.cells = min_cells, 
                                  min.features = min_features)
  
  # Analyse mitochondrial percentage
  sc.seurat[["percent.mt"]] <- PercentageFeatureSet(sc.seurat, pattern = "^MT-")
  
  # Return object
  return(sc.seurat)
  
}


# Plot QC metrics
plot_metrics <- function(sc.seurat, 
                         vln_plot, 
                         feature_scatter_plot, 
                         plots.on = T) {
  
  # Visualize QC metrics as a violin plot
  if (plots.on) { png(vln_plot, width=1000,height=1000,units="px") }
  print(VlnPlot(object = sc.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  if (plots.on) { dev.off() }
  
  # Visualize feature-feature relationships
  plot1 <- FeatureScatter(object = sc.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") 
  plot2 <- FeatureScatter(object = sc.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  if (plots.on) { png(feature_scatter_plot, width=2000,height=1000,units="px") }
  print(CombinePlots(plots = list(plot1,plot2)))
  if (plots.on) { dev.off() }
  
}



