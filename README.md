# chromatic

`chromatic` is an R package for analyzing chromatin state heterogeneity in high-dimensional epigenomic sequencing assays, such as scATAC-seq. `chromatic` calculates cell-level *chromatin entropy* scores based on a cell's usage of different chromatin states (enhancers, promoters, silencers, etc). Chromatin entropy scores can indicate important cell populations and reveal epigenomic remodeling underlying cell state changes across many systems including neurodegeneration, cancer, and development.  

`chromatic` is under active development, and is currently an **alpha** version. This means that some features of the final package are missing, and current features are subject to change throughout development. The first major release of `chromatic` will coincide with a forthcoming publication.

## Installation

```bash
# create a new conda environment
conda create -n chromatic -c conda-forge python=3.11 mamba 
conda activate chromatic

# install dependencies using mamba
mamba install -c conda-forge -c bioconda \
  r-base=4.4 \
  r-seurat \
  r-signac \
  r-devtools \
  r-hdf5r \
  bioconductor-tfbstools \
  bioconductor-rtracklayer \
  bioconductor-genomicranges \
  bioconductor-motifmatchr \ 
  macs2
```

Now install `chromatic` within R using `devtools`.

```r 
devtools::install_github('smorabit/chromatic')
```