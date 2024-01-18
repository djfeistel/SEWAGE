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

**Reference_genomes.fasta**, **Proportion_Read_metadata.tsv**, and **parameters.txt** file names cannot be modified. The fasta, meta data, and fastq files listed above are named with default settings and can be modified using the ```--fastq_name``` and ```--amplicon_fasta_name``` flags.

### Comments about creating amplicons:
For each primer set, both the forward and reverse primer must be found in a reference sequence in order for amplification to occur, i.e., a single nucleotide mismatch between the reference sequence and a primer sequence results in a non-amplification event. Thus, if at least one primer is missing, the amplicion will not be 'amplified' for that primer and there will be no defline in the amplicon muliti-fasta file.  However the meta data file will indicate which primers did not amplify by stating the primer name followed by "No Amplification" in the **amplicon_sequence** column.