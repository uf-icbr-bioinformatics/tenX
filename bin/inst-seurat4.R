install.packages("BiocManager")
install.packages("remotes")

remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_github("satijalab/azimuth", "seurat5")
BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("hdf5r")
