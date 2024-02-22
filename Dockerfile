FROM continuumio/miniconda3:4.10.3

RUN apt-get update && apt-get install -y build-essential libarchive-dev

# Installing mamba
RUN conda install -c conda-forge mamba

# Creating Conda Envs
COPY enviroment.yml .

RUN mamba env create -f enviroment.yml
RUN echo "source activate modi-finder-web" > ~/.bashrc
ENV PATH /opt/conda/envs/modi-finder-web/bin:$PATH

COPY . /app
WORKDIR /app
RUN chmod +x /app/run_server.sh