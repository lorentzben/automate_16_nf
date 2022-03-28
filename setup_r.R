#!/usr/bin/env Rscript --vanilla
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
if(!require(rmarkdown)) {install.packages("rmarkdown", repos="http://cran.us.r-project.org")}
if(!require(renv)) {install.packages("renv",repos="http://cran.us.r-project.org")}
renv::init()
renv::restore(packages = "renv")
if(!require(remotes)){install.packages("remotes",repos="http://cran.us.r-project.org")}
if(!require(devtools)){install.packages("devtools",repos="http://cran.us.r-project.org")}
if(!require(jamba)){remotes::install_github("jmw86069/jamba@0.0.6.900")}
remotes::install_github("tidyverse/ggplot2@v3.3.2")
remotes::install_github("vegandevs/vegan@v2.5-7")
if(!require(ampvis2)){remotes::install_github("MadsAlbertsen/ampvis2@2.6.8")}
if(!require(ggvegan)){remotes::install_github("gavinsimpson/ggvegan@4bc6ee9945dd9229ed486409c0acab9413b8c9af")}
if(!require(ggConvexHull)){remotes::install_github("cmartin/ggConvexHull@660f4094da44dd500c3c0684b9c5c20c21ee823a")}