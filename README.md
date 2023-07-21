# SEWAGE
### Synthetically Enriched Wastewater-like sequence data for Assessing Genomic and Environmental populations

***SEWAGE*** is entirely written in Python 3 and is tested on Python 3.8.3. As of now, the only dependencies 
are tqdm<sup>1</sup>, NumPy<sup>2</sup> and ART<sup>3</sup>. However, in the future ART will 
be replaced with an internal algorithm for simulating reads.

For detailed information about the tool: ```SEWAGE --details```

## Installing dependencies via conda
Note: It is not necessary to use conda as long as you have ***tqdm*** and ***NumPY*** installed for Python 3 and ***art_illuina*** in your $PATH
```
conda create -n SEWAGE_env python==3.8.3 --yes
conda activate SEWAGE_env
pip install tqdm numpy==1.21.0
conda install -c bioconda art --yes
```
Once you have installed the conda environment (if using) you can add SEWAGE as a symlink to your bin  
```
ln -s <pathway/to/SEWAGE> <pathway/to/bin>
```


## Usage
For detailed information about paramters:  
```SEWAGE -h```

Minimal Usage:  
```SEWAGE -i <input>```

Help Menu:
```
usage: SEWAGE.py [-h] [-i PATHWAY or FILE] [-o STR] [-O DIR] [-p {r,e}] [-rs INT] [-pf INT]
                 [-ss {HS10,HS20,HS25,HSXn,HSXt,MinS,MSv1,MSv3,NS50,GA1,GA2}] [-l INT] [-m INT] [-s INT] [-ir INT] [-ir2 INT] [-dr INT]
                 [-dr2 INT] [-rsA INT] [-k INT] [-nf INT] [--details]

Simulation of Environmental Wastewater sequence data for the Analysis of Genomics and Epidemiology

optional arguments:
  -h, --help            show this help message and exit

Input and Output Parameters:
  -i/--in is required

  -i PATHWAY or FILE, --in PATHWAY or FILE
                        Pathway to directory with FASTA files or a text file with a list of pathways to FASTA files. NOTE: FASTA files
                        must end in .fasta, .fa, or .fsa when a pathway is specified.
  -o STR, --out STR     Name of output directory for storage (if not specificed, default is 'SEWAGE_' + 10 random alphanumeric
                        characters).
  -O DIR, --out_pathway DIR
                        Pathway to where output directory is stored [default='.']

Proportion options [default is -p r -rs 13]:

  -p {r,e}, --proportion {r,e}
                        Generate random (r) or equal (e) proportions of reads
  -rs INT, --rndSeed INT
                        Random seed for generateing proportions [default=13]

ART parameters:
  Default parameters listed below are for simulating "perfect" reads at 150bp and can be modified as needed. All other parameters not
  listed here are in default setting or not used as defined by "art_illumia" and cannot be access via SEWAGE. Please be familiar with
  how "art_illumina" functins before modifying these parameters

  -pf INT, --pfold INT  Value to be mutiplied by proportion for use with '--fcov' from 'art_illumina' [default=1000] Example:
                        proportion*pfold=fcov or fold coverage
  -ss {HS10,HS20,HS25,HSXn,HSXt,MinS,MSv1,MSv3,NS50,GA1,GA2}, --seqSys {HS10,HS20,HS25,HSXn,HSXt,MinS,MSv1,MSv3,NS50,GA1,GA2}
                        From 'art_illumina': 'The name of Illumina sequencing system of the built-in profile used for simulation'
                        [default=HS25]. Note: chosing a differnt Illumina sequecning system may require modifying the 'art_illumina'
                        paramters beforehand
  -l INT, --len INT     From 'art_illumina': the length of reads to be simulated [default=150]
  -m INT, --mflen INT   From 'art_illumina': the mean size of DNA/RNA fragments for paired-end simulations [default=250]
  -s INT, --sdev INT    From 'art_illumina': the standard deviation of DNA/RNA fragment size for paired-end simulations [default=1]
  -ir INT, --insRate INT
                        From 'art_illumina': the first-read insertion rate [default=0]
  -ir2 INT, --insRate2 INT
                        From 'art_illumina': the second-read insertion rate [default=0]
  -dr INT, --delRate INT
                        From 'art_illumina': the first-read deletion rate [default=0]
  -dr2 INT, --delRate2 INT
                        From 'art_illumina': the second-read deletion rate [default=0]
  -rsA INT, --rndSeed_art_illumina INT
                        From 'art_illumina': the seed for random number generator [default=13]
  -k INT, --maxIndel INT
                        From 'art_illumina': the maximum total number of insertion and deletion per read [default=0]
  -nf INT, --maskN INT  From 'art_illumina': the cutoff frequency of 'N' in a window size of the read length for masking genomic
                        regions [default=0]

Tool Description:
  Detailed information about the tool

  --details

Minimal usage: ./SEWAGE -i <input>
```
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