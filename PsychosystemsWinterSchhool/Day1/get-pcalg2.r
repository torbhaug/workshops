

# use this to get the correct packages

BiocManager::install("graph","RBGL","pcalg")


# or use the following

BiocManager::install(version = "4.3")
install.packages("BiocInstaller")
source("http://bioconductor.org/biocLite.R")
biocLite("graph")