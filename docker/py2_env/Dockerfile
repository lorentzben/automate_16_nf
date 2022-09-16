#FROM continuumio/anaconda
#Run echo "Reconstruction of 2.0 tag"
#WORKDIR /py2_env
# Create Conda Env
#COPY python2_corr_env.yml python2_env.yml
#RUN conda env create -f python2_env.yml 

#ARG DEBIAN_FRONTEND=noninteractive

#RUN apt-get --allow-releaseinfo-change update &&  \
#    apt-get upgrade -y &&  \
#    apt-get install -y procps

# Conda command in new env
#SHELL ["conda", "run", "-n", "automate_16s_nf","/bin/bash","-c"]

# Activate the environment 
#RUN echo "Environment now exists"
#This code no longer works the biobakery channel for conda no 
#longer contains lefse:1.0.0_dev_9adc3a62460e 

FROM biobakery/lefse:1.0.0_dev_9adc3a62460e

RUN apt install -y procps

# Conda command in new env
SHELL ["conda", "run", "-n", "automate_16s_nf","/bin/bash","-c"]