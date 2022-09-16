FROM continuumio/anaconda
WORKDIR /py2_env
# Create Conda Env
COPY py2_corr_env.yml python2_env.yml
RUN conda env create -f python2_env.yml
RUN apt-get -y update --allow-releaseinfo-change && \
	apt-get -y upgrade && \
    apt-get install -y procps

# Conda command in new env
SHELL ["conda", "run", "-n", "automate_16s_nf","/bin/bash","-c"]

# Activate the environment 
RUN echo "Environment now exists"

