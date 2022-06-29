#! /usr/bin/env Rscript --vanilla
#Rscript -e ".libPaths('./r_lib/')"
.libPaths('/renv_dev/r_lib/')
Rscript -e "renv::activate()"
#Rscript -e "renv::init()"
#Rscript -e "renv::install('rmarkdown')"
Rscript -e "renv::restore(library='./renv/library/R-4.1/x86_64-pc-linux-gnu/', lockfile='./renv.lock')"
q()

