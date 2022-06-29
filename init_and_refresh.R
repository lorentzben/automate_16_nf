#! /usr/bin/env Rscript --vanilla
#Rscript -e ".libPaths('./r_lib/')"
.libPaths('./r_lib/')
renv::settings$use.cache(TRUE)
renv::activate()
#Rscript -e "renv::init()"
#Rscript -e "renv::install('rmarkdown')"
renv::restore(library='./renv/library/R-4.1/x86_64-pc-linux-gnu/', lockfile='./renv.lock')
q()

