# SEWAGE

## Synthetically Enriched Wastewater-like sequence data for Assessing Genomic and Environmental populations

***SEWAGE*** is entirely written in Python 3 and is tested on Python 3.8.3. As of now, the only dependencies 
are tqdm<sup>1</sup> and ART<sup>2</sup>. However, in the future ART will 
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
### Citations

1. https://github.com/tqdm/tqdm
2. Huang, W., Li, L., Myers, J. R., & Marth, G. T. (2012). ART: a next-generation sequencing read simulator. Bioinformatics, 28(4), 593-594.