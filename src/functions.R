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





# Reduce dimensionality ---------------------------------------------------

# Dimnesionality reduction:
# Scales the data, identifies main features based on highly expressed genes 
# and runs PCA to detetermine dimensionality of the dataset (nb of dimensions
# to reduce the datast to)


# Data scaling and feature selection
runSCTransform_and_PCA <- function(sc.seurat) {
  
  # Scale data with SCTransform OR ScaleData
  sc.seurat <- SCTransform(object = sc.seurat)
  
  # Run PCA
  sc.seurat <- RunPCA(object = sc.seurat)
  
  # Return
  return(sc.seurat)
  
}

# Plot variable features
plot_variable_features <- function(sc.seurat, 
                                   variable_feature_plot) {
  
  # Identify the 10 most highly variable genes
  var.features.top10 <- head(x = VariableFeatures(object = sc.seurat), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(object = sc.seurat)
  plot2 <- LabelPoints(plot = plot1, points = var.features.top10, repel = TRUE)

  if (plots.on) { png(variable_feature_plot, width = 1500, height=1000, unit="px") }
  print(plot2)
  if (plots.on) { dev.off() }
  
}

# Determine dimensionality with PCA
determine_dimensionality <- function(sc.seurat, 
                                     elbow_plot) {
  
  # Plot elbowplot to determine dimensionality
  if (plots.on) { png(elbow_plot, width = 1000, height=1000, unit="px") }
  print(ElbowPlot(object = sc.seurat, ndims=50))
  if (plots.on) { dev.off() }
  
}



# Cluster cells -----------------------------------------------------------

# Cluster and visualise UMAP
cluster_and_visualise <- function(sc.seurat, 
                                  nb_dims, 
                                  UMAP_plot, 
                                  resolution = 0.5,
                                  plots.on = TRUE) {
  
  # Run UMAP
  sc.seurat <- RunUMAP(object = sc.seurat, dims = 1:nb_dims)
  
  # Find clusters
  sc.seurat <- FindNeighbors(object = sc.seurat, dims = 1:nb_dims)
  sc.seurat <- FindClusters(object = sc.seurat, resolution = resolution)
  
  # Plot UMAP
  if (plots.on) { png(UMAP_plot, width = 800, height=800, unit="px") }
  print(DimPlot(object = sc.seurat, reduction = 'umap'))
  if (plots.on) { dev.off() }
  
  return(sc.seurat)
  
}




# Cell cycle score --------------------------------------------------------

# Get cell cycle info
get_cell_cycle_info <- function(sc.seurat) {
  sc.seurat <- Seurat::CellCycleScoring(
    object = sc.seurat,
    g2m.features = Seurat::cc.genes$g2m.genes,
    s.features = Seurat::cc.genes$s.genes
  )
  return(sc.seurat)
}

# Plot certain markers
plot_markers <- function(sc.seurat, 
                         downstream_analysis_plots, 
                         list_markers, 
                         colors = NA,
                         suffix = "",
                         reduction = "umap",
                         plots.on = T) {
  
  for (marker in list_markers) {
    print(paste0("Plotting -- marker: ", marker, "..."))
    
    
    if (!sum(c(rownames(sc.seurat@assays[["RNA"]]), rownames(sc.seurat)) == marker)) {
      # Stop if marker not in sc.seurat
      print(paste0("Gene not found:", marker))
    } else {
      # Else, create the fetaure plot
      if (plots.on) { png(file.path(downstream_analysis_plots, paste0(marker, suffix, ".png"))) }
      if (is.na(colors)) {
        print(Seurat::FeaturePlot(object = sc.seurat, marker, reduction = reduction))
      }
      else {
        print(Seurat::FeaturePlot(object = sc.seurat, marker, reduction = reduction,
                                  cols = colors))
      }
      if (plots.on) { dev.off() }
    }
  }
}




# IC scores from Zynovyev method ------------------------------------------

add_scores_IC_zynovyev <- function(sc.seurat,
                            file_IC_scores,
                            name) {

  # Ajouter IC-EwS
  read_ICs <- read.table(file = file_IC_scores, 
                         fill = TRUE,
                         header=TRUE)
  
  IC_EwS_score <- read_ICs$IC10..1
  IC_ECM_score <- read_ICs$IC30..1
  IC_G1S_score <- read_ICs$IC2.
  IC_G2M_score <- read_ICs$IC1.
  if (name == "pdx1058") { 
    names(IC_EwS_score) <- paste0(read_ICs$SAMPLE, "-1") 
    names(IC_ECM_score) <- paste0(read_ICs$SAMPLE, "-1") 
    names(IC_G1S_score) <- paste0(read_ICs$SAMPLE, "-1") 
    names(IC_G2M_score) <- paste0(read_ICs$SAMPLE, "-1") 
  } else {
    names(IC_EwS_score) <- read_ICs$SAMPLE
    names(IC_ECM_score) <- read_ICs$SAMPLE
    names(IC_G1S_score) <- read_ICs$SAMPLE
    names(IC_G2M_score) <- read_ICs$SAMPLE
  }
  sc.seurat <- AddMetaData(object = sc.seurat,
                           metadata = IC_EwS_score,
                           col.name = "Zynov_IC_EwS_score")
  sc.seurat <- AddMetaData(object = sc.seurat,
                           metadata = IC_ECM_score,
                           col.name = "Zynov_IC_ECM_score")
  sc.seurat <- AddMetaData(object = sc.seurat,
                           metadata = IC_G1S_score,
                           col.name = "Zynov_IC_G1S_score")
  sc.seurat <- AddMetaData(object = sc.seurat,
                           metadata = IC_G2M_score,
                           col.name = "Zynov_IC_G2M_score")
  
  return(sc.seurat)
  
}



