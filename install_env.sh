#!/bin/bash

# install dependancies and activate SEWAGE conda environment
set -e 

conda create -n SEWAGE_env python==3.8.3 --yes
conda activate SEWAGE_env
pip install tqdm
conda install -c bioconda art --yes