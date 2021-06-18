#! /usr/bin/env Rscript --vanilla
if(!require(renv)) {install.packages("renv",repos="http://cran.us.r-project.org")}
renv::init()
install.packages("BiocManager")
BiocManager::install(version = "3.13")
if(!require(dplyr)){install.packages("dplyr",repos="http://cran.us.r-project.org")}
if(!require(tibble)) {install.packages("tibble",repos="http://cran.us.r-project.org")}
if(!require(qiime2R)) {devtools::install_github("jbisanz/qiime2R")} # current version is 0.99.20
if(!require(phyloseq)) {install.packages("phyloseq",repos="http://cran.us.r-project.org")}
if(!require(jamba)){remotes::install_github("jmw86069/jamba@0.0.6.900")}
renv::restore()
q()