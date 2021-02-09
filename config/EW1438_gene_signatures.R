#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# EW1438_gene_signatures.R                                                     #
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
# Rscript ./src/print_gene_signatures.R ./config/EW1438_gene_signatures.R






# Config for EW1438 dataset -----------------------------------------------


# General parameters
main_folder <- "/Users/maudgautier/Documents/github-savings/single-cell-basics"
name <- "EW1438"
dataset <- file.path(main_folder, "data", name, "filtered_feature_bc_matrix")
output_plots <- file.path(main_folder, "plots", name)
file_rds <- file.path(main_folder, "rds", paste0(name, ".rds"))

# List signatures
list_signatures <- c("IC10")
file_signatures <- "/Users/maudgautier/Documents/data/tmp_from_calcsub/nadege/List_IC10plus.csv"


