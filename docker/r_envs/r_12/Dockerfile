FROM rocker/verse:4.1.0
LABEL maintainer="lorentzb"
RUN apt-get update -y && apt-get install -y  libpng-dev  libjpeg-dev  zlib1g-dev  libicu-dev  pandoc  make  libcurl4-openssl-dev  cmake  libssl-dev  git  libfreetype6-dev  libfribidi-dev  libharfbuzz-dev  libxml2-dev  libfontconfig1-dev  libtiff-dev  libgit2-dev  libglpk-dev  libgmp3-dev  imagemagick  libmagick++-dev  gsfonts && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = TRUE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages(c("renv","remotes"))'
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'
CMD ["bash"]