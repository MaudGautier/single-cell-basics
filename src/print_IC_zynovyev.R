#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# print_IC_zynovyev.R                                                          #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script allows to create feature plots for ICs calculated with the method
# by A. Zynovyev available at the following github repo : 
# https://github.com/sysbio-curie/EwingSingleCellDataAnalysis



# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript ./src/print_IC_zynovyev.R your_config_file.R




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
path_IC_scores <- file.path(output_plots, "IC_scores_zynov")




# Initialisation and quality check ----------------------------------------

# Create all subfolders for the sample
for (folder in c(output_plots, path_IC_scores)) {
  init_folder(folder)
}

# Get seurat object
if (!file.exists(file_rds)) { 
  print(paste0("No file found at ", file_rds))
}
seurat_object <- readRDS(file_rds)



# IC scores from Zynovyev method ------------------------------------------

# Get ICs from Zynovyev method
seurat_object <- add_scores_IC_zynovyev(seurat_object, zynov_IC_scores, name)

# Plot ICs
plot_markers(seurat_object, path_IC_scores,
             list_markers = c("Zynov_IC_EwS_score", "Zynov_IC_ECM_score", "Zynov_IC_G1S_score", "Zynov_IC_G2M_score")
             )
