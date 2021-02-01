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
# Rscript ./src/print_IC_zynovyev.R ./config/EW1438_IC_zynov.R




# Config for EW1438 dataset -----------------------------------------------


# General parameters
main_folder <- "/Users/maudgautier/Documents/github-savings/single-cell-basics"
name <- "EW1438"
dataset <- file.path(main_folder, "data", name, "filtered_feature_bc_matrix")
output_plots <- file.path(main_folder, "plots", name)
file_rds <- file.path(main_folder, "rds", paste0(name, ".rds"))

# IC scores throught Zynovyev's method
zynov_IC_scores <- paste0("/Users/maudgautier/Documents/sc-projects/reanalysis-aynaud2020/github-EwingSingleCellDataAnalysis/data/",
                         name,
                         "_no_zero_col_nu2k.txt.moduleAverages")

