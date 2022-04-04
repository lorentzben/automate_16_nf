#! /usr/bin/env Rscript --vanilla
if(!require(renv)) {install.packages("renv",repos="http://cran.us.r-project.org")}
library(renv)
renv::init(bare=T)
renv::isolate()
renv::settings$snapshot.type("all")
renv::install("rmarkdown")
install.packages("BiocManager")
BiocManager::install(version = "3.12", ask=FALSE)
BiocManager::repositories(version="3.12")
BiocManager::install("remotes")

