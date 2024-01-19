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
usage: SEWAGE [-h] -i STR -s STR [-n STR] [-o STR] [-t] [-p {r,e,d}] [-dg STR] [-dp FLOAT] [-ps INT] [-rl INT]
              [-fl INT] [-cd INT] [-rs INT]

Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities

options:
  -h, --help            show this help message and exit

Input and Scheme Parameters:
  Required flags

  -i STR, --infasta STR
                        Multifasta file or single column list of pathways to fasta files
  -s STR, --scheme STR  Available primer scheme: Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"] VarSkip:
                        Long read = ["vsl1a"]; Short-read = ["vss1a", "vss2a", "vss2b"])

Output naming:

  -n STR, --file_prefix_name STR
                        File name prefix for all generated data [default="SEWAGE"]
  -o STR, --storage_dir STR
                        Directory name for data storage [default="SEWAGE_Workspace"]
  -t, --time_stamp      Append date and time stamp to the storage directory [default: False]

Proportion options:

  -p {r,e,d}, --proportion_model {r,e,d}
                        Generate equal (e), random (r), or dominate (d) variant of concern proportions of reads
                        [default: r]
  -dg STR, --dVOC_genome STR
                        Name of dVOC. NOTE: name of dVOC must match the defline of the reference fasta file
  -dp FLOAT, --dVOC_proporiton FLOAT
                        Proportion of dDOV [default: 0.8]
  -ps INT, --proportion_seed INT
                        Random seed number for reproducing proportions [default: 13]

Read generator options:

  -rl INT, --read_length INT
                        Read length in bp (value should not exceed amplicon length or will workflow fail)
                        [default: 150]
  -fl INT, --frag_length INT
                        Fragment length in bp (value should not exceed amplicon length or will workflow fail)
                        [default: 300]
  -cd INT, --coverage_depth INT
                        Total sequence depth coverage for each fastq file [default: 500]
  -rs INT, --read_seed INT
                        Random seed number for reproducing reads [default: 13]

Minimal Usage: SEWAGE -i <multi.fasta> -s <scheme>
```
### Output:

Assuming Minimal Usage ```SEWAGE -i <multi.fasta> -s <scheme>```

From the code above, results will be stored in a directory called ```SEWAGE_Workspace``` in the current working directory. You can change the name of the storage directory by using the ```--storage_dir``` flag. You can also append a date/time stamp with the ```--time_stamp``` flag. The files created are:

|File Name|Description|
|:----|:----|
|SEWAGE_Reference_genomes.fasta|Fasta file with reference genomes|
|SEWAGE_parameters.txt|Parameters used when running ```SEWAGE```|
|SEWAGE_metadata.tsv|Metadata related for calculating reads|
|SEWAGE_amplicons.fasta|Detected amplicons for all reference genomes|
|SEWAGE_R1.fastq|Forward reads|
|SEWAGE_R2.fastq|Reverse reads|

All files listed above are named with default settings and can be modified using the either the ```--file_prefix_name``` flag to change the prefix for all other files listed above.  

### Comments about generating amplicons:
For each primer set, both the forward and reverse primer must be found in a reference sequence in order for amplification to occur; i.e., a single nucleotide mismatch between the reference sequence and a primer sequence results in a non-amplification event. Thus, if at least one primer is missing, the amplicion will not be 'amplified' for that primer and there will be no defline in the ```SEWAGE_amplicons.fasta``` file.  However the ```SEWAGE_metadata.tsv``` file will indicate which primers did not amplify by stating the primer name followed by "No Amplification" in the **amplicon_sequence** column.

### Comments about generating proportions
By default, the ```--proportion_model``` flag is set to ```r``` or random. For a more refined proporiton, users can set ```--proportion_model d``` which assigns a single genome at random as the dominant variant of concern (dVOC) at a 0.8 default proporiton, and the proportion can be modified using the ```--dVOC_proporiton``` flag. Users can also use the ```--dVOC_genome``` flag to assign which genome is the dVOC (full defline in original reference fasta as of now). By setting the ```--proportion_seed``` flag, users can reproduce proportions.

### Comments about generating reads
Default parameters when using the ```minimal usage``` command have been optimized to result in 150bp F/R reads form 300bp fragments at a total depth of coverage of approximatly 500X between all reference sequences supplied (NOTE: estimated depth of coverage will result in higher coverage due to overlaping regions of tiled amplicons when using defult settings).  
When parameters are adjusted, users are recommended to manually values in the **SEWAGE_metadata.tsv** file for the proporitons, total number of fragments generated per ampicon, and the length of reads before running any downstream analyses and adjust accordingly to the experiment.  
Users can modify the length of reads with the ```--read_length``` flag and length of fragments with the ```--frag_length``` flag to mimic different sequencing platforms (e.g. 75bp, 150bp, 250bp, etc...). Note that depending on the scheme used, the default 250bp might need to be adjusted as only the ends of the amplicons are 'sequenced' and are nto broken up into smaller fragments.   

### Additional Notes
**NOTE:** Reads are assigned the highest Q-score possible for now.  
**NOTE:** Only "short-read" fastq data is available as of now. I am activly working on adding long read output data.

### Future Additions:
1. Fragments pulled from normal distribution with std for variable read lengths 
2. Quality scores assigned to bases for high adn low quality
3. INDELS
4. Long Read 
5. Conda package
