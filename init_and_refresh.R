#! /usr/bin/env Rscript --vanilla
#Rscript -e ".libPaths('./r_lib/')"
.libPaths('./r_lib/')
Sys.setenv(R_INSTALL_STAGED = "false")
options(install.opts = "--no-test-load")
Sys.setenv(RENV_DEFAULT_R_LIBS_SITE="/renv_dev/r_lib")
Sys.setenv(RENV_PATHS_LIBRARY="/renv_dev/r_lib")
Sys.setenv(R_LIBS_SITE="/renv_dev/r_lib")
#RENV_DEFAULT_R_LIBS_SITE
#RENV_PATHS_LIBRARY
#renv::settings$use.cache(TRUE)
renv::activate()
#Rscript -e "renv::init()"
#Rscript -e "renv::install('rmarkdown')"
renv::restore(library='./renv/library/R-4.1/x86_64-pc-linux-gnu/', lockfile='./renv.lock')
q()

