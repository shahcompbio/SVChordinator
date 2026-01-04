# Use the Miniconda base image
FROM continuumio/miniconda3:24.7.1-0
LABEL authors="asherpreskasteinberg"

# Clone the Minda repository
RUN git clone https://github.com/shahcompbio/minda.git

# Set the working directory to Minda and create the environment
WORKDIR /minda
RUN conda env create --name minda --file environment.yml

# Configure the container to activate the environment and switch to the minda directory on start
# Activate conda environment
RUN echo "conda activate minda" >> ~/.bashrc
ENV PATH /opt/conda/envs/minda/bin:$PATH

# Install additional packages in the "minda" environment
RUN conda install -n minda -c conda-forge intervaltree -y

