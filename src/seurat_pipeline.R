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





# Sample-specific parameters ----------------------------------------------

# Paths specific to sample
#cluster_dir <- file.path(sample_dir, "clusters-markers/")
path_QC <- file.path(output_plots, "01_QC")
path_DR <- file.path(output_plots, "02_dim_reduction")
path_CL <- file.path(output_plots, "03_clustering")
path_DA <- file.path(output_plots, "04_downstream_analysis")



# # Path plots
# path_vln_plot <- file.path(path_QC, "vln_plot.png")
# path_feature_scatter_plot <- file.path(path_QC, "feature-scatter.png")
# path_variable_feature_plot <- file.path(path_DR, "variable_features.png")
# path_elbow_plot <- file.path(path_DR, "elbow-plot.png")
# path_UMAP_plot <- file.path(path_CL, "UMAP.png")


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


