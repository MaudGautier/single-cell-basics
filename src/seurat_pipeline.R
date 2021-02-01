#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# seurat_pipeline.R                                                            #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script contains the successive steps of a basic Seurat pipeline to
# analyse single-cell datasets sequenced with the 10X technology.


# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript ./src/seurat_pipeline.R your_config_file.R




# Parse arguments ---------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]

# Deal with case where repo not given
# (in this case, assumes that it is the working directory)
if (!exists(deparse(substitute(main_folder)))) {
  main_folder <- getwd()
}



# Import variables and functions ------------------------------------------

print("Sourcing and importing libraries...")


# Source functions
source(file.path(main_folder, "src/functions.R"))

# Source configuration file
source(config_file)

# Import libraries
import_libraries()


# Define paths to output plots
#cluster_dir <- file.path(sample_dir, "clusters-markers/")
path_QC <- file.path(output_plots, "01_QC")
path_DR <- file.path(output_plots, "02_dim_reduction")
path_CL <- file.path(output_plots, "03_clustering")
path_DA <- file.path(output_plots, "04_downstream_analysis")




# Initialisation and quality check ----------------------------------------

# Create all subfolders for the sample
for (folder in c(output_plots, path_QC, path_DR, path_CL, path_DA)) {
  init_folder(folder)
}

# Get seurat object
seurat_object <- get_seurat_object(path_dataset = dataset,
                                   sample_name = name, 
                                   min_cells = 1, 
                                   min_features = 200)

# Plot metrics
plot_metrics(seurat_object, 
             vln_plot = file.path(path_QC, "vln_plot.png"), 
             feature_scatter_plot = file.path(path_QC, "feature-scatter.png"), 
             plots.on = plots.on)





# Dimensionality reduction ------------------------------------------------

# Based on plots, select cells to keep
if (all_features) { 
  seurat_object <- subset(x = seurat_object, 
                          subset = percent.mt < max_percent_mt)
} else {
  seurat_object <- subset(x = seurat_object, 
                          subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_mt)
}

# Normalise data
seurat_object_NORMALISED <- NormalizeData(object = seurat_object) 

# scTransform and detect dimensionality
seurat_object <- runSCTransform_and_PCA(seurat_object)
plot_variable_features(seurat_object, 
                       variable_feature_plot = file.path(path_DR, "variable_features.png"))
determine_dimensionality(seurat_object, 
                         elbow_plot = file.path(path_DR, "elbow-plot.png"))





# Cluster cells -----------------------------------------------------------

# Cluster cells after the SCTransform
seurat_object <- cluster_and_visualise(seurat_object, nb_dims, 
                                       UMAP_plot = file.path(path_CL, "UMAP.png"), 
                                       resolution = resolution)

# Cluster after standard normalisation
seurat_object_NORMALISED <- ScaleData(object = seurat_object_NORMALISED, features = rownames(x = seurat_object_NORMALISED))
seurat_object_NORMALISED <- FindVariableFeatures(object = seurat_object_NORMALISED)
seurat_object_NORMALISED <- RunPCA(object = seurat_object_NORMALISED)
seurat_object_NORMALISED <- cluster_and_visualise(seurat_object_NORMALISED, nb_dims, 
                                                  UMAP_plot = NA, 
                                                  resolution = resolution,
                                                  plots.on = FALSE)



# Individual markers and cell cycle ---------------------------------------


# Plot features for a list of markers
plot_markers(seurat_object, path_DA, 
             list_markers = list_markers)

# Define cell-cycle S.Score and G2M.score and plot them
seurat_object_NORMALISED <- get_cell_cycle_info(seurat_object_NORMALISED)
seurat_object <- AddMetaData(object = seurat_object,
                             metadata = seurat_object_NORMALISED$S.Score,
                             col.name = "S.Score")
seurat_object <- AddMetaData(object = seurat_object,
                             metadata = seurat_object_NORMALISED$G2M.Score,
                             col.name = "G2M.Score")
plot_markers(seurat_object, path_DA, 
             list_markers = c("S.Score", "G2M.Score", "nFeature_RNA"), 
             colors = c("blue", "yellow", "red"), 
             suffix = "_BYR")
plot_markers(seurat_object, path_DA, 
             list_markers = c("S.Score", "G2M.Score", "nFeature_RNA"))




# Save RDS ----------------------------------------------------------------

# Save RDS
if (save_rds) {
  folder_rds <- dirname(file_rds)
  init_folder(folder_rds)
  print(paste0("Seurat object will be saved at ", file_rds))
  saveRDS(object = seurat_object, file = file_rds)
}





# IC scores from Zynovyev method ------------------------------------------


# Get ICs from Zynovyev method
seurat_object <- add_scores_IC_zynovyev(seurat_object, zynov_IC_scores, name)

# Plot ICs
plot_markers(seurat_object, path_DA, 
             list_markers = c("Zynov_IC_EwS_score", "Zynov_IC_ECM_score", "Zynov_IC_G1S_score", "Zynov_IC_G2M_score")
             )


