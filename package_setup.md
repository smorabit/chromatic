1. Create a new local conda environment to work with Chromatic 
    (also can serve as installation instructions)

```bash

conda create -n chromatic -c conda-forge python=3.11 mamba 

conda activate chromatic

mamba install -c conda-forge -c bioconda \
  r-base=4.4 \
  r-seurat \
  r-signac \
  r-devtools \
  r-hdf5r \
  bioconductor-tfbstools \
  bioconductor-rtracklayer \
  bioconductor-genomicranges \
  bioconductor-motifmatchr

mamba install -c bioconda macs3

```

```r 

library(devtools)
library(roxygen2)
library(pkgdown)
library(usethis)
library(dplyr)

setwd("/Users/sam.morabito/Library/CloudStorage/GoogleDrive-sam.morabito@cnag.eu/My Drive/projects/")

# Create the package!
usethis::create_package("chromatic")
devtools::create("chromatic")

setwd('chromatic')

# add functions to the .R 
# (doing this outside R with Bash)

# document functions 
devtools::document()

# install the package:
devtools::install()

usethis::use_gpl3_license()

# Set up the necessary files for pkgdown (only run this once)
usethis::use_git()
usethis::use_pkgdown_github_pages()

# build the intitial version of the site!
pkgdown::build_site()

# commit and push to Github!


```



