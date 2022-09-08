#FROM rocker/r-ver:4.1.0
#LABEL maintainer="lorentzb"
#RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
#  && apt-get install -y git-core && \
#  apt-get install -y procps && \
#  apt-get install -y pandoc && \
#  apt-get install -y nano && \
#  apt-get install -y texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra lmodern

#RUN ["install2.r", "dplyr", "forcats", "ggplot2", "purrr", "readr", "stringr", "tibble", "tidyr", "tidyverse"]
#RUN ["install2.r", "-r https://bioconductor.org/packages/3.14/bioc -r https://bioconductor.org/packages/3.14/data/annotation -r https://bioconductor.org/packages/3.14/data/experiment -r https://bioconductor.org/packages/3.14/workflows", "BiocGenerics", "phyloseq", "qiime2R"]
#WORKDIR /payload/
#CMD ["R"]

#RUN ["install2.r", "dplyr", "forcats", "ggplot2", "purrr", "readr", "remotes", "stringr", "tibble", "tidyr", "tidyverse"]
#RUN ["install2.r", "-r https://bioconductor.org/packages/3.14/bioc -r https://bioconductor.org/packages/3.14/data/annotation -r https://bioconductor.org/packages/3.14/data/experiment -r https://bioconductor.org/packages/3.14/workflows", "BiocGenerics", "phyloseq"]
#RUN ["installGithub.r", "jbisanz/qiime2R@d1ad96657ada993cf6c2841b29113a4f635c6b56"]
#WORKDIR /payload/
#COPY qiime2R/* . 
#RUN Rscript -e "install.packages('./qiime2R')"

#CMD ["R"]


FROM rocker/verse:4.1.0
LABEL maintainer="lorentzb"
RUN apt-get update -y && apt-get install -y procps libpng-dev  libjpeg-dev  zlib1g-dev  libicu-dev  pandoc  make  libcurl4-openssl-dev  libssl-dev  git  libfreetype6-dev  libfribidi-dev  libharfbuzz-dev  libxml2-dev  libfontconfig1-dev  libtiff-dev  libgit2-dev  libglpk-dev  libgmp3-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = TRUE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages(c("renv","remotes"))'
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'
CMD ["bash"]