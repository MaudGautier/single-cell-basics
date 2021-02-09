#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# print_gene_signatures.R                                                      #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script contains the successive steps of a basic Seurat pipeline to
# analyse single-cell datasets sequenced with the 10X technology.
# This pipeline will run : 
# - the creation of quality check plots
# - 


# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript ./src/sprint_gene_signatures.R your_config_file.R





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
path_signatures <- file.path(output_plots, "06_gene_signatures")



# Initialisation and quality check ----------------------------------------

# Create all subfolders for the sample
for (folder in c(output_plots, path_signatures)) {
  init_folder(folder)
}

# Get seurat object
if (!file.exists(file_rds)) { 
  print(paste0("No file found at ", file_rds))
}
seurat_object <- readRDS(file_rds)



# Print gene signatures ---------------------------------------------------

plot_signature(seurat_object, file_signatures, list_signatures, path_signatures)

