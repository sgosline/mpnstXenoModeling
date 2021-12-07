##set up R env
install.packages("rmarkdown")
install.packages("remotes")
install.packages("gridExtra")
install.packages("BiocManager")
install.packages("ggfortify")
install.packages("ggrepel")

library(BiocManager)
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("tximportData")

library(remotes)
remotes::install_github("rstudio/reticulate")
remotes::install_github("sgosline/mpnstXenoModeling")
remotes::install_github("biodataganache/leapR")
