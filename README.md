# SEWAGE 

### <u>S</u>ynthetically <u>E</u>ngineered <u>W</u>astewater sequence data for <u>A</u>ssessing <u>G</u>enomic <u>E</u>ntities

SEWAGE is a tool for generating reproducible sequence data representing a heterogeneous population of closely related species at various proportions. Specifically, it was designed to mirror wastewater sequencing data that uses a tiled-amplicon PCR approaches resulting in a mixed SARS-CoV-2 population derived from a wastewater sample.

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
usage: SEWAGE [-h] -i STR -s STR [-n STR] [-sd STR] [-p {r,e,d}] [-dg STR] [-dp FLOAT] [-ps INT] [-q STR]
              [-rl INT] [-auto] [-cd INT] [-mr INT] [-rs INT]

Synthetically Engineered Wastewater-like sequence data for Assessing Genomic and Environmental populations

optional arguments:
  -h, --help            show this help message and exit

Input and Scheme Parameters:
  Required flags

  -i STR, --infasta STR
                        Multifasta file or single column list of pathways to fasta files
  -s STR, --scheme STR  Available primer scheme: Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"] VarSkip:
                        Long read = ["vsl1a"]; Short-read = ["vss1a", "vss2a", "vss2b"])

Naming Parameters:

  -n STR, --file_prefix_name STR
                        File name prefix for all amplicon and meta data generated [default="SEWAGE_"]
  -sd STR, --storage_dir STR
                        Directory name for data storage [default="SEWAGE_[data_time]"]

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
  -auto, --auto_read_length_detection
                        Automatically set the read length to the maximum possible based on amplicon length
                        (i.e., half max amplicon + 50bp). [default: False]
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
|SEWAGE_Reference_genomes.fasta|Fasta file with reference genomes|
|SEWAGE_parameters.txt|Parameters used when running ```SEWAGE```|
|SEWAGE_Proportion_Read_metadata.tsv|Proportional data used to calculate reads|
|SEWAGE_amplicons.fasta|Amplicons detected for all reference genomes|
|SEWAGE_amplicons_metadata.tsv|Amplicon meta data for amplified and non-amplified primers|
|SEWAGE_R1.fastq|Forward reads|
|SEWAGE_R2.fastq|Reverse reads|

All files listed above are named with default settings and can be modified using the either the ```--fastq_name``` flag to change the prefix for the fastq files or the ```--file_prefix_name``` flag to change the prefix to all other files listed above.  

### Comments about generating amplicons:
For each primer set, both the forward and reverse primer must be found in a reference sequence in order for amplification to occur, i.e., a single nucleotide mismatch between the reference sequence and a primer sequence results in a non-amplification event. Thus, if at least one primer is missing, the amplicion will not be 'amplified' for that primer and there will be no defline in the amplicon muliti-fasta file.  However the meta data file will indicate which primers did not amplify by stating the primer name followed by "No Amplification" in the **amplicon_sequence** column.

### Comments about generating proportions
By default, the ```--proportion_model``` flag is set to ```r``` or random. For a more refined proporiton, users can set ```--proportion_model d``` which assigns a single genome at random as the dominant variant of concern (dVOC) at a 0.8 default proporiton, but the proportion can be modified using the ```--dVOC_proporiton``` flag. Users can also use the ```--dVOC_genome``` flag to assign which genome is the dVOC. By setting the ```--proportion_seed``` flag, users can reproduce proportions.

### Comments about generating reads
Default parameters when using the ```minimal usage``` command have been optimized for the **ARTIC V5.3.2** primer set and will result in 250bp F/R reads at approximatly 500X total depth of coverage between all reference seequences supplied. This approach works well when generateing short-read paired-end data on default settings. However, when using differnet schemes or adjusting other paramters, users should refer to the **SEWAGE_Proportion_Read_metadata.tsv** file for the total number of read generated per ampicon and the length of reads before running any downstream analyses and adjust accordingly to the experiment. Users can either modify the length of reads with the ```--read_length``` flag to mimic different sequencing platforms (e.g. 75bp, 150bp, 250bp, etc...). Note that depending on the scheme used, the default 250bp might need to be adjusted as only the ends of the amplicons are 'sequenced' and are nto broken up into smaller fragments. Users can also use the ```--auto_read_length_detection``` to guarantee reads that overlap. This works by taking half the length of the largest amplicon deteted and divided it by two and then add 50pb. Note that using the ```--auto_read_length_detection``` is not intended to mimic real-world sequening platforms and is intended to be used in experimental cases.  

### Additional Notes
**NOTE:** Reads are assigned the highest Q-score possible.  
**NOTE:** When using the long-read primer scheme for Varskips "vsl1a", default paramters should be modified and the resulting data should be manually checked before any downstream analysis with the fastq files as this has not been fully worked out yet. 