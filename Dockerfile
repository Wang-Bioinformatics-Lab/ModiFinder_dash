FROM ubuntu:22.04

RUN apt-get update && apt-get install -y build-essential libarchive-dev wget git

# Install Mamba
ENV CONDA_DIR /opt/conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && /bin/bash ~/miniforge.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Adding to bashrc
RUN echo "export PATH=$CONDA_DIR:$PATH" >> ~/.bashrc

# Creating Conda Envs
COPY conda-env.yml .

RUN mamba env create -f conda-env.yml -n smallmol_mod_site_localization_dash

# Copying in the module and installing
#COPY ModiFinder_base /app/ModiFinder_base
#RUN /bin/bash -c "source activate smallmol_mod_site_localization_dash && pip install -e /app/ModiFinder_base"

RUN /bin/bash -c "source activate smallmol_mod_site_localization_dash && pip install modifinder==1.5.5"

COPY . /app
WORKDIR /app
RUN chmod +x /app/run_server.sh