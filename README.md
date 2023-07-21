# SEWAGE

## Synthetically Enriched Wastewater-like sequence data for Assessing Genomic and Environmental populations

***SEWAGE*** is entirely written in Python 3 and is tested on Python 3.8.3. As of now, the only dependencies 
are tqdm (https://github.com/tqdm/tqdm) and ART (Huang, Weichun, et al. "ART: a next-generation 
sequencing read simulator." Bioinformatics 28.4 (2012): 593-594.). However, in the future ART will 
be replaced with an our own algorithm for simulating reads.

### Installing dependencies via conda
Note: It is not necessary to use conda as long as you have **tqdm** and **art_illuina** in your $PATH
```
conda create -n SEWAGE_env python==3.8.3 --yes
conda activate SEWAGE_env
pip install pip install tqdm
conda install -c bioconda art --yes
conda install -c bioconda art --yes
```
