FROM r-base:4.1.0

ENV DEBIAN_FRONTEND noninteractive
ENV INITRD No

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y zlib1g-dev && \
    apt-get install -y libudunits2-dev && \
    apt-get install -y libfontconfig1-dev && \
    apt-get install -y bash && \
    apt-get install -y libgdal-dev && \
    apt-get install -y procps && \
    apt-get install -y pandoc && \
    apt-get install -y nano && \
	apt-get install -y cmake && \
    apt-get install -y texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra lmodern &&\
	apt-get install -y libssl-dev && \
	apt-get install -y libglpk-dev && \
	apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev  libssl-dev && \
    apt-get install -y libharfbuzz-dev libfribidi-dev



CMD ["mkdir r_lib"]

#RUN Rscript -e ".libPaths('./r_lib/')"
#RUN Rscript -e "install.packages('renv')"
RUN Rscript -e "install.packages('rmarkdown')"
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "install.packages('devtools')"

RUN Rscript -e "Sys.setenv(R_INSTALL_STAGED = FALSE)"
#RUN Rscript -e "renv::init(bare=T)"

CMD ["bash"]
