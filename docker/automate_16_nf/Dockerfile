FROM qiime2/core:2020.8
WORKDIR /nextflow_env_test
# Create Conda Env
#COPY environment.yml .
#RUN conda env create -f environment.yml

# Conda command in new env
#SHELL ["conda", "run", "-n", "automate_16s_nf","/bin/bash","-c"]

# Activate the environment 
RUN qiime info
RUN python3 -m pip install Path
RUN python3 -m pip install pandas
RUN python3 -m pip install numpy
#RUN conda install -c bioconda -c conda-forge picrust2 -y 
RUN conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict
RUN pip install git+https://github.com/bokulich-lab/RESCRIPt.git
RUN apt-get --allow-releaseinfo-change update && apt-get -y upgrade && \
    apt-get clean && \
    apt-get -qq update && apt-get -y upgrade && \
    apt-get -y install procps && \ 
    apt-get -y install unzip 