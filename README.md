# SEWAGE 

### Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities

SEWAGE is a tool for generating reproducible sequence data that represents a heterogeneous population of closely related species at various proportions. Specifically, it was designed to mirror sequence data that resembles a mixed SARS-CoV-2 population derived from a wastewater sample by using targeted-enrichment or tiled-amplicon PCR approaches.

SEWAGE currently offers two main functionalities for producing wastewater-like sequence data: 
1. The ability to produce amplicons that mimic tiled-amplicon sequencing approaches for each input genome using a set of primers.
2. Create Illumina short-read data sets from those amplicions that mimic a heterogeneous populations of closely related species at various proportions.

### Comments about SEWAGE:
1. SEWAGE currently offers only SARS-CoV2 [ARTIC](https://github.com/artic-network/primer-schemes) and [VarSkip](https://github.com/nebiolabs/VarSkip) primer sets for creating amplicions. However, we are actively working on allowing users to supply their own primer sets.
2. As of now, SEWAGE can only create Illumia short-reads dur to its reliability on the tool ART<sup>1</sup>.  We are actively working creating an internal alrogithm for generating short and long read data. Updates will be made when available.

## Installing dependencies via conda
SEWAGE is written in Python 3 version 3.8.3 and the only dependencies are tqdm<sup>2</sup>, NumPy<sup>3</sup> and ART<sup>1</sup>. However, ART will eventually be replaced in future versions of SEWAGE with an internal algorithm for simulating long and short reads.  

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

## Generate amplicons from a set of reference genomes
The first step is to generate amplicons from reference genomes of any size.  Reference genomes must be store as fasta files and each reference fasta file must only include a single reference genome. The reference fasta files should be stored in a directory together with the extention ".fasta", ".fa", or ".fsa".  Other files may be stored in the directory with the reference files.

There are two required flags: The ```--fasta``` flag is used to tell SEWAGE where the fasta files are stored and the ```--scheme``` flag sets the primers to use (stored in the ```scheme``` directory). SEWAGE currently offers only SARS-CoV2 [ARTIC](https://github.com/artic-network/primer-schemes) and [VarSkip](https://github.com/nebiolabs/VarSkip) primer sets for creating amplicions. However, we are actively working on allowing users to supply their own primer sets to use with other species. The optional flags ```--output``` and ```--pathway``` set the name of the output (default="SEWAGE_amplicons") and the pathway for those amplicons to be stored (default='.'), respectivly.

Minimal Usage:  
```
SEWAGE amplicon -f <pathway/to/fasta-files> -s <scheme>
```
Help Menu:
```
usage: SEWAGE amplicon [-h] -f FASTA -s SCHEME [-o STR] [-p PATHWAY]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Single or multi-fasta reference. Multi-fasta files should be unique genomes for each defline/sequence. Output
                        will produce "${defline}_amplicon.fasta" files for each genome
  -s SCHEME, --scheme SCHEME
                        Primer scheme: (Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"], VarSkip = ["vsl1a", "vss1a", "vss2a",
                        "vss2b"])
  -o STR, --output STR  Output Prefix name for fasta file [default="SEWAGE_amplicons.fasta"]
  -p PATHWAY, --pathway PATHWAY
                        Pathway to storgage directory [default="."]
```
### SEWAGE Amplicon Output:

From the code above, the results would be stored in a directory called ```SEWAGE_amplicon``` by default in the current working directory.  For each reference genome, a ```_amplicons.fasta``` and ```_amplicons.log``` file is created.  The multi-fasta file contains the amplicons generated and the log is a tab seperated file with:  

|primer-name|reference-name|start-position|end-potision|amplicon-length (bp)|
|:----:|:----:|:-----------:|:-----------:|:-----------:|

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
