# SEWAGE
### Synthetically Engineered Wastewater sequence data for Assessing Genomic variants in Environmental populations 
### Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities

SEWAGE is a tool for generating reproducible sequence data representing a heterogeneous population of closely related species. Specifically, it was designed to mirror sequence data that resembles a mixed SARS-CoV-2 population derived from a wastewater sample by using targeted enrichment or tiled amplicon approaches. SEWAGE was developed to help assess the accuracy of different alignment/mappers and relative abundance calculation tools used in the National Wastewater Survaillance Systems (NWSS) SARS-CoV-2 Wastewater detection pipeline "AquaScope" (CITE GITHUB/GITLAB).

SEWAGE currently offers two main functionalities: 1) the ability to produce amplicons for each genome from a set of closely related reference genomes using a set of primers and 2) create Illumina short-read data sets from those amplicions that mimic heterogeneous populations of closely related species at various proportions. SEWAGE currently offers SARS-CoV2 [ARTIC](https://github.com/artic-network/primer-schemes) and [VarSkip](https://github.com/nebiolabs/VarSkip) primer sets for creating amplicions. However, we are actively working on allowing users to supply their own primer sets.  

## Installing dependencies via conda
SEWAGE is written in Python 3 version 3.8.3 and the only dependencies are tqdm<sup>1</sup>, NumPy<sup>2</sup> and ART<sup>3</sup>. However, ART will eventually be replaced in future version of SEWAGE with an internal algorithm for simulating long and short reads.  

*Note: It is not necessary to use conda as long as you have ***tqdm*** and ***NumPY*** installed and ***art_illuina*** in your $PATH*
```
conda create -n SEWAGE_env python==3.8.3 --yes
conda activate SEWAGE_env
pip install tqdm numpy==1.21.0
conda install -c bioconda art --yes
```
Once you have installed the conda environment you can add SEWAGE as a symlink to your bin  
```
ln -s <pathway/to/SEWAGE> <pathway/to/bin>
```

## Usage for generating amplicons

Minimal Usage:  
```
SEWAGE amplicon -f <pathway/to/fasta-files> -s <scheme>
```
Help Menu:
```

```
### SEWAGE Amplicon Output:

|Name |Type |Description |
|:----:|:----:|:-----------:|
|


For detailed information about paramters:  
```SEWAGE -h```

### Output
|Name |Type |Description |
|:----:|:----:|:-----------:|
|proportions_list.txt|text file|Two column file with col1=pathway/to/file.fasta and col2=proporiton (float)|
|<input_name>_1.fastq|fastq file|Forward reads|
|<input_name>_2.fastq|fastq file|Reverse reads|
|<input_name>_logfile.log|text file|Information about the run|

## Citations

1. https://github.com/tqdm/tqdm
2. https://github.com/numpy/numpy
2. Huang, W., Li, L., Myers, J. R., & Marth, G. T. (2012). ART: a next-generation sequencing read simulator. Bioinformatics, 28(4), 593-594.