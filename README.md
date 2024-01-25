# SEWAGE 

### Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities

SEWAGE is a tool for generating reproducible *in silico* sequence data representing a heterogeneous population of closely related species at various proportions. Specifically, it was designed to mirror SARS-CoV-2 sequencing data from a wastewater sample representing a heterogeneous population using a tiled-amplicon PCR approaches.

SEWAGE currently provides two main functionalities:
 1. Generate amplicons from a set of primers using a tiled-amplicon approach for each SARS-CoV-2 genome.  
 2. Produce short-read, paired-end shotgun sequencing data formatted as fastq files, which mimics a heterogeneous population at varying proportions, similar to sequence data derived from wastewater samples.   

SEWAGE currently offers SARS-CoV2 [ARTIC](https://github.com/artic-network/primer-schemes) and [VarSkip](https://github.com/nebiolabs/VarSkip) primer sets for creating amplicions.

### SEWAGE caveats
1. no variation in read length (yet)
2. no INDELs (yet)
3. estimated depth of coverage > input coverage depth bc of tiled amplicon approach. primers are not removed from amplicions (but the option will be available)
4. Amplicon dropout only occurs if one or both primers are not found in reference sequence (but option for regional drop out will soon be available)
5. Short-read paired end only (for now)
## Dependencies
```
numpy
pandas
```

## Usage
```
Minimal Usage: SEWAGE -i <multi.fasta> -s <scheme>
```
Using the ```Minimal Usage``` command is the simplest way to generate amplicons and short-read paired-end *in silico* data. If you prefer to run this and move one, refer to the **Output** section below for more detail.  
If you prefer to have more control over the output data, here are a few helpful commands explaining what most flags perfrom.  
```
SEWAGE --file_prefix_name <prefix_name> --storage_dir <sotrage_dir_name> --time_stamp -i <multi.fasta> -s <scheme>
```
The ```--file_prefix_name``` flag will attached a prefix to the begining of all files generated in the main workflow. Use ```--storage_dir``` to give a name or a pathway to where the data is stored in the current working directory (if the directory exists then the code terminates). The ```--time_stamp``` flag will append the date and time as **YYYMMDD_HHMMSS** to the end of the directoy.
Help Menu:
```
usage: SEWAGE [-h] -i STR -s STR [-n STR] [-o STR] [-t] [-p {r,e,d}] [-dg STR] [-dp FLOAT] [-ps INT] [-rl INT] [-fl INT] [-cd INT] [-rs INT]
              [-m INT INT] [-sd INT] [-bp INT] [-red INT]

Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities

options:
  -h, --help            show this help message and exit

Input and Scheme Parameters:
  Required flags

  -i STR, --infasta STR
                        Multifasta file or single column list of pathways to fasta files
  -s STR, --scheme STR  Available primer scheme: Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"] VarSkip: Long read = ["vsl1a"]; Short-read =
                        ["vss1a", "vss2a", "vss2b"])

Output naming:

  -n STR, --file_prefix_name STR
                        File name prefix for all generated data [default="SEWAGE"]
  -o STR, --storage_dir STR
                        Directory name for data storage [default="SEWAGE_Workspace"]
  -t, --time_stamp      Append date and time stamp to the storage directory [default: False]

Proportion options:

  -p {r,e,d}, --proportion_model {r,e,d}
                        Generate equal (e), random (r), or dominate (d) variant of concern proportions of reads [default: r]
  -dg STR, --dVOC_genome STR
                        Name of dVOC. NOTE: name of dVOC must match the defline of the reference fasta file
  -dp FLOAT, --dVOC_proporiton FLOAT
                        Proportion of dDOV [default: 0.8]
  -ps INT, --proportion_seed INT
                        Random seed number for reproducing proportions [default: 13]

Read generator options:

  -rl INT, --read_length INT
                        Read length in bp (value should not exceed amplicon length or will workflow fail) [default: 150]
  -fl INT, --frag_length INT
                        Fragment length in bp (value should not exceed amplicon length or will workflow fail) [default: 300]
  -cd INT, --coverage_depth INT
                        Total sequence depth coverage for each fastq file [default: 500]
  -rs INT, --read_seed INT
                        Random seed number for reproducing reads [default: 13]

Read generator options:

  -m INT INT, --min_max_q INT INT
                        Minimum and Maximum mean range quality scores [default: 30 40]
  -sd INT, --std_dev_q INT
                        Standard deviation for quality scores [default: 3]
  -bp INT, --startup_effect_bp INT
                        Number of bp's at the begining of a read that have reduced quality [default: 10]
  -red INT, --startup_effect_q_reduction INT
                        Quality score reduction for '-bp' (maximum mean minus '-red') [default: 4]

Minimal Usage: SEWAGE -i <multi.fasta> -s <scheme>
```


## Output

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
Default parameters when using the ```minimal usage``` command have been optimized to result in 150bp F/R reads form 300bp fragments at a total depth of coverage of approximatly 500X between all reference sequences supplied (NOTE: estimated depth of coverage will result in higher coverage due to overlaping regions of tiled amplicons when using defult settings). Furthermore, there are no qualitye scores other than 40 assotiated with bases (as of now). 
When parameters are adjusted, users are recommended to manually inspect values in the **SEWAGE_metadata.tsv** file for the proporitons, total number of fragments generated per ampicon, and the length of reads, and any other *in silico* data before running any downstream analyses (adjust accordingly to the experiment).  
Users can modify the length of reads with vi the ```--read_length``` flag and length of fragments with the ```--frag_length``` flag to mimic different sequencing platforms (e.g. 75bp, 150bp, 250bp, etc...).

### Additional Notes
**NOTE:** Reads are assigned the highest Q-score possible for now.  
**NOTE:** Only "short-read" fastq data is available as of now. I am activly working on adding long read output data.

### Future Additions:
1. Fragments pulled from normal distribution with std for variable read lengths 
2. Quality scores assigned to bases for high and low quality
3. INDELS addedd 
4. Long Read 
5. Conda package
6. Possibly adding adapters to ends of smaller sequences?
