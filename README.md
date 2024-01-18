# SEWAGE 

### Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities

SEWAGE is a tool for generating reproducible sequence data that represents a heterogeneous population of closely related species at various proportions. Specifically, it was designed to mirror sequence data that resembles a mixed SARS-CoV-2 population derived from a wastewater sample by using targeted-enrichment or tiled-amplicon PCR approaches.

SEWAGE offers ability Produces amplicons that mimic tiled-amplicon sequencing approaches for each input genome using a set of primers and create short-read or long-read data sets (i.e.,fastq files) from those amplicions that mimic a heterogeneous populations of closely related species at various proportions for producing wastewater-like sequence data.

SEWAGE currently offers SARS-CoV2 [ARTIC](https://github.com/artic-network/primer-schemes) and [VarSkip](https://github.com/nebiolabs/VarSkip) primer sets for creating amplicions.

## Dependencies
```
numpy
pandas
```

## Usage
Minimal Usage:  
```
SEWAGE -i <multi.fasta> -s <scheme>
```
Help Menu:
```
usage: SEWAGE [-h] -i STR -s STR [-afn STR] [-sd STR] [-p {r,e,d}] [-dg STR] [-dp FLOAT] [-ps INT] [-q STR]
              [-rl INT] [-cd INT] [-mr INT] [-rs INT]

Synthetically Engineered Wastewater-like sequence data for Assessing Genomic and Environmental populations

options:
  -h, --help            show this help message and exit

Input and Scheme Parameters:
  Required flags

  -i STR, --infasta STR
                        Multifasta file or single column list of pathways to fasta files
  -s STR, --scheme STR  Available primer scheme: Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"] VarSkip:
                        Long read = ["vsl1a"]; Short-read = ["vss1a", "vss2a", "vss2b"])

Amplicon Parameters:

  -afn STR, --amplicon_fasta_name STR
                        File name for amplicon data frame and fasta file [default="SEWAGE_amplicons"]
  -sd STR, --storage_dir STR
                        Directory name for amplicon data storage [default="SEWAGE_[data_time]"]

Proportion options:

  -p {r,e,d}, --proportion_model {r,e,d}
                        Generate equal (e), random (r), or dominate (d) variant of concern proportions of reads
                        [default: d]
  -dg STR, --dVOC_genome STR
                        Name of dVOC. NOTE: name of dVOC must match the defline of the reference fasta file
  -dp FLOAT, --dVOC_proporiton FLOAT
                        Proportion of dDOV [default: 0.8]
  -ps INT, --proportion_seed INT
                        Random seed number for reproducing proportions [default: 13]

Read generator options:

  -q STR, --fastq_name STR
                        Name of fastq files for F/R reads [default: SEWAGE_{R1/R2}.fastq]
  -rl INT, --read_length INT
                        Read length in bp (value should not exceed amplicon length or will workflow fail)
                        [default: 250]
  -cd INT, --coverage_depth INT
                        Total sequence depth coverage for each fastq file [default: 500]
  -mr INT, --max_reads INT
                        Total number of reads for each fastq file [default: None]. NOTE: If set,
                        --coverage_depth is ignored.
  -rs INT, --read_seed INT
                        Random seed number for reproducing reads [default: 13]

Fast Start: SEWAGE -i <multi.fasta> -s <scheme>
```
### Output:

Assuming Minimal Usage ```SEWAGE -i <multi.fasta> -s <scheme>```

From the code above, results will be stored in a directory called ```SEWAGE_{YYYYMMDD}_{HHMMSS}``` in the current working directory. You can change the name of the storage directory by using the ```--storage_dir``` flag. Several files will be created:

|File Name|Description|
|:----|:----|
|Reference_genomes.fasta|Fasta file with reference genomes|
|Proportion_Read_metadata.tsv|Proportional data used to calculate reads|
|parameters.txt|Parameters used when running ```SEWAGE```|
|SEWAGE_amplicons.fasta|Amplicons detected for all reference genomes|
|SEWAGE_amplicons_metadata.tsv|Amplicon meta data for amplified and non-amplified primers|
|SEWAGE_R1.fastq|Forward reads|
|SEWAGE_R2.fastq|Reverse reads|

With the exception Reference_genomes.fasta, Proportion_Read_metadata.tsv, and parameters.txt files,Files listed above are named with default settings and can be modified using the ```--fastq_name``` and ```--amplicon_fasta_name``` flags.

### Comments about creating amplicons:
Both forward and reverse primers for each primer set must be found in a reference sequence in order for amplification to occur. If at least one primer is missing, the amplicion will not be amplified for that primer and there will be no defline in the amplicon muliti-fasta file.  However the log file will indicate which primers did not amplify by stating the primer name followed by "No Amplification".


The ```SEWAGE_amplicon``` directory (or whatever you choose to name it) is used as input for the next command: ```SEWAGE enrich```

## Generate wastewater sequence data
Using the previously generated amplicons, we can supply the ```SEWAGE_amplicon``` directory as input with the ```--in``` flag.  There are many options when using this command; however, the only required flag is the ```--in``` and the other other flags are set with defaults which you can read below in the help menu section.  The main default setting to be aware of are the ```Proportion options```. By default, these are set to ```-p v -V 0.8 -rs 13``` which means that the proportion (```-p```) is set up with the ```v``` choice which creates a heterogeneous sequence data set at random proporitons (which can be reproduced using the random seed ```-rs``` flag) and a single reference genomes is chosen at random to be the dominant variant of concern (dVOC), representing 80% of the population (i.e. ```-V 0.8```) relative abundance.  If you do not want a single genome to represent the dVOC, use the ```-p r``` flag which creates heterogeneous sequence data at random propotions.  Using the random seed flag ```-rs``` allows for these sequence reads to be reproducible.

For those parameters available to use in ART<sup>1</sup>, we have set many of those to optimal settings so that the reads generated are not altered (i.e. no indels or mutations).  However, you have the ability to modify such paramters including but not limited to the rate of inserions and/or deletions and quality scores, in the case of creating exploratory datasets.

Minimal Usage:
```SEWAGE enrich -i pathway/to/SEWAGE_amplicon```

Help Menu:
```
usage: SEWAGE enrich [-h] [-i PATHWAY or FILE] [-o STR] [-O DIR] [-p {v,r,e}] [-V FLOAT] [-rs INT] [-pf INT]
                     [-ss {HS10,HS20,HS25,HSXn,HSXt,MinS,MSv1,MSv3,NS50,GA1,GA2}] [-l INT] [-m INT] [-s INT] [-ir INT] [-ir2 INT]
                     [-dr INT] [-dr2 INT] [-rsA INT] [-k INT] [-nf INT] [-qL INT] [-qU INT] [-qs INT] [-qs2 INT]

optional arguments:
  -h, --help            show this help message and exit

Input and Output Parameters:
  -i/--in is required

  -i PATHWAY or FILE, --in PATHWAY or FILE
                        Pathway to directory with FASTA files or a text file with a list of pathways to FASTA files. NOTE: FASTA files
                        must end in .fasta, .fa, or .fsa when a pathway is specified.
  -o STR, --out STR     Name of output directory for storage [default='SEWAGE_enrich'].
  -O DIR, --out_pathway DIR
                        Pathway to where output directory is stored [default='.']

Proportion options [default is -p v -V 0.8 -rs 13]:

  -p {v,r,e}, --proportion {v,r,e}
                        Generate a dVOC (v), random (r), or equal (e) proportion read set
  -V FLOAT, -dVOC FLOAT
                        Proporiton of dVOC [default=0.8]. NOTE: End proportion value might vary slighlty. Valuse >= 1 are converted to
                        0.99
  -rs INT, --rndSeed INT
                        Random seed used for generating dVOC and random proportions [default=13]

ART parameters:
  Default parameters listed below are for simulating "perfect" reads at 150bp and can be modified as needed. All other parameters not
  listed below are either not in use or in default setting as defined by "art_illumia" and cannot be access via SEWAGE. Please be
  familiar with how "art_illumina" functins before modifying these parameters

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
  -qL INT, --minQ INT   From 'art_illumina': the minimum base quality score [default=28]
  -qU INT, --maxQ INT   From 'art_illumina': the maximum base quality score [default=40]
  -qs INT, --qShift INT
                        From 'art_illumina': the amount to shift every first-read quality score by
  -qs2 INT, --qShift2 INT
                        From 'art_illumina': the amount to shift every second-read quality score by
```

### Output data from ```SEWAGE enrich```

The output data will be stored in the current working directory with the default name ```SEWAGE_enrich``` if the ```--out``` flag is not specified.  These files and directories include:
|Name |Type |Description |
|:----:|:----:|:-----------:|
|proportions_list.txt|text file|Two column file with col1=pathway/to/file.fasta and col2=proporiton (float)|
|SEWAGE_enrich_1.fastq|fastq file|Forward reads for all regerence genomes|
|SEWAGE_enrich_2.fastq|fastq file|Reverse reads for all regerence genomes|
|SEWAGE_enrich_logfile.log|text file|Information about the run|
|00.RAWREADS|directory|Raw reads storage for each reference genome|
|01.LOGS|directory|Information about the "art_illumina" sequencing run|

### Comments about creating sequence data
As of v0.1.0, amplicon coverage is generated equally 

## Citations

1. Huang, W., Li, L., Myers, J. R., & Marth, G. T. (2012). ART: a next-generation sequencing read simulator. Bioinformatics, 28(4), 593-594.
2. https://github.com/tqdm/tqdm
3. https://github.com/numpy/numpy
