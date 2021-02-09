# single-cell-basics
This repo contains basic scripts to perform the first-stage analysis of single-cell 10X datasets.

## Requirements

These scripts require the installation of R, version 4.0.3.

The following libraries must be installed:

* ggplot2
* dplyr
* Seurat
* patchwork
* SingleR

These libraries can be installed as follows (from the R CLI):
```
# Base packages
install.packages("ggplot2")
install.packages('dplyr')
install.packages('Seurat')
install.packages('patchwork')

# BiocManager
install.packages("BiocManager")

# Installation with BiocManager
BiocManager::install("SingleR")
```


## Usage

**Basic seurat pipeline**
The basic seurat pipeline can be run by creating a configuration file in a similar fashion to the [``config/EW1438_basics.R`` file](https://github.com/MaudGautier/single-cell-basics/tree/main/config/EW1438_basics.R) and executing it using this command line:
```
Rscript ./src/seurat_pipeline.R ./config/EW1438_basics.R
```

**Gene signatures**
To plot gene signatures, you can create a configuration file in a similar fashion to the [``config/EW1438_gene_signatures.R`` file](https://github.com/MaudGautier/single-cell-basics/tree/main/config/EW1438_gene_signatures.R) and execute it using this command line:
```
Rscript ./src/print_gene_signatures.R ./config/EW1438_gene_signatures.R
```





