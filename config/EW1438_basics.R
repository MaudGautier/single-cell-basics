#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# EW1438_basics.R                                                              #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This is the configuration file to analyse the EW1438 single-cell dataset
# with the basic Seurat pipeline.


# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript ./src/seurat_pipeline.R ./config/EW1438_basics.R




# Config for EW1438 dataset -----------------------------------------------


# General parameters
main_folder <- "/Users/maudgautier/Documents/github-savings/single-cell-basics"
name <- "EW1438"
dataset <- file.path(main_folder, "data", name, "filtered_feature_bc_matrix")
output_plots <- file.path(main_folder, "plots", name)


# Pipeline parameters
plots.on <- TRUE

# Parameteters to (re)define after evaluating QC plots
all_features <- TRUE
min_nFeature_RNA <- 200    # only if all_features is FALSE
max_nFeature_RNA <- 10000  # only if all_features is FALSE
max_percent_mt <- 20


# Number of dimensions to keep  (after evaluating the elbow plot)
nb_dims <- 10

# Clusterisation
resolution <- 0.5 # Play on the resolution value to modify the number of clusters

# Print individual markers
list_markers <- c("CD99", "FLI1", "PAPPA", "IL1RAP", "LOX", "ICAM1")

# Save RDS
save_rds <- TRUE
file_rds <- file.path(main_folder, "rds", paste0(name, ".rds"))



