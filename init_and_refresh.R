#! /usr/bin/env Rscript --vanilla
if(!require(renv)) {install.packages("renv",repos="http://cran.us.r-project.org")}
renv::init()
install.packages("BiocManager")
BiocManager::install(version = "3.13", ask=FALSE)
if(!require(rhdf5)){renv::install("bioc::rhdf5")}
if(!require(Biostrings)){renv::install("bioc::Biostrings")}
if(!require(Biobase)){renv::install("bioc::Biobase")}
if(!require(dplyr)){renv::install("dplyr")}
if(!require(tibble)) {renv::install("tibble")}
if(!require(phyloseq)) {BiocManager::install("phyloseq", version="3.13")}
if(!require(qiime2R)) {renv::install("jbisanz/qiime2R@d1ad96657ada993cf6c2841b29113a4f635c6b56")} # current version is 0.99.20
if(!require(jamba)){renv::install("jmw86069/jamba@0.0.6.900")}

#renv::restore()
q()

