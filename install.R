# install.R
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("remotes")
install.packages("BiocManager")

remotes::install_cran("Seurat")
BiocManager::install("MAST")
install.packages("reticulate")
