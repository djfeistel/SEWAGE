# SEWAGE 

### <u>S</u>ynthetically <u>E</u>ngineered <u>W</u>astewater sequence data for <u>A</u>ssessing <u>G</u>enomic <u>E</u>ntities

SEWAGE is a tool for generating reproducible sequence data representing a heterogeneous population of closely related species at various proportions. Specifically, it was designed to mirror wastewater sequencing data that uses a tiled-amplicon PCR approaches resulting in a mixed SARS-CoV-2 population derived from a wastewater sample.

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
6. deflines renamed to genome_1, genome_2 for dVOCname? but deflines stay the same in results files
## Dependencies
```
numpy
pandas
```
## Workflow
1. Load reference fasta file
2. Generate amplicons
3. Assign proporitons
4. Calculate number of reads per amplicon per reference genomes based on coverage
5. Write fastq R1 and R2 files
6. Save all metadata assotiated with generatign amplicons and reads.
## Usage
```
Minimal Usage: SEWAGE -i <multi.fasta> -s <scheme>
```
Using the ```Minimal Usage``` command is the simplest way to generate amplicons and short-read paired-end *in silico* data. If you prefer to run this and move one, refer to the **Output** section below for more detail.  

If you prefer to have more control over the output data, here are a few helpful commands explaining what most flags perfrom.  
### Output options
```
SEWAGE --file_prefix_name <prefix_name> --storage_dir <sotrage_dir_name> --time_stamp -i <multi.fasta> -s <scheme>
```
The ```--file_prefix_name``` flag will attached a prefix to the begining of all files generated in the main workflow. Use ```--storage_dir``` to give a name or a pathway to where the data is stored in the current working directory (if the directory exists then the code terminates). The ```--time_stamp``` flag will append the date and time as **YYYMMDD_HHMMSS** to the end of the directoy.  

### Proportion options
```
SEWAGE  --proportion_model {r,e,d} --dVOC_genome <defline> --dVOC_proporiton <FLOAT> --proportion_seed <INT> -i <multi.fasta> -s <scheme>
```
When using the ```--proportion_model``` flag, there are three options: **random** or ```r``` will randomly assign proportions to the reference genomes (this model is the default), **equal** or ```e``` will assign equal proprtions across all refreence gneomes, and **dominante variant of concern (dVOC)** or ```d``` will randomly choose one reference genome and assign a 0.8 proporiton as the default. When the model is set to ```d```, you can modify the default value with ```--dVOC_proporiton```. Depending on how mnay reference genomes used, the actual proportion may differet slightly from the assigned proporion. To assign a reference genome with the dVOC proporiton, use the ```--dVOC_genome``` flag. Note that, as of now, you must supply the complete defline found in the reference fasta file in order to assing it as the dVOC. Set the ```--proportion_seed``` flag to reproduce results (default 13).  

### Read generator options
```
SEWAGE  --read_length <INT> --frag_length <INT> --coverage_depth <INT> --read_seed INT -i <multi.fasta> -s <scheme>
```
User can modify the default values for ```--read_length``` and ```--frag_length```. Note that there is no variation in read/fragment length in this algorithm as of now i.e., all reads/fragments will be the size supplied. Increasing/decreasing the ```--coverage_depth``` will produce more or less reads. Set the ```--read_seed``` flag to reproduce results (default 13).  

### Quality score options
```
SEWAGE  --min_max_q < INT INT> --std_dev_q <INT> --startup_effect_bp <INT> --startup_effect_q_reduction <INT> -i <multi.fasta> -s <scheme>
```
Setting the ```--min_max_q``` with two integers between 0-40 will give the minimum and maximum Phred qualitiy scores (Q). The ```--std_dev_q``` flag sets the standard deviation at which those values are randomly selected to produced variation in the Q-scores. Modifying the ```--startup_effect_bp``` flag sets the number of base pairs at the begining of a read that have slighly decreased quality relative to the rest of the reads due to the initial cycle effect, and the ```--startup_effect_q_reduction``` changed the magnatude of that reduction (default 4; MAX(Q) - 4 = startup_effect_bp). Note: these parematers do not output data the align perfectly with real data, though it shows similarites.  Even though Q-scores are assigned, no mutaitons or indels are introduced into any reads. Thus, the Q-scores do not directly effect the *in silico* data.  

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

Output options:

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

Quality score options:

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
|SEWAGE_amplicons.tsv|Metadata related for generating amplicons|
|SEWAGE_amplicons.fasta|Detected amplicons for all reference genomes|
|SEWAGE_reads.tsv|Metadata related for calculating reads|
|SEWAGE_R1.fastq|Forward reads|
|SEWAGE_R2.fastq|Reverse reads|

All files listed above are named with default settings and can be modified using the either the ```--file_prefix_name``` flag to change the prefix for all other files listed above. Refer to **Usage** seciton from detailed inforamtoin on flag usage.  

## Discussion 
### Comments about generating amplicons:
For each primer pair in a set, both the forward and reverse primer must be found in a reference sequence in order for amplification to occur; i.e., a single nucleotide mismatch between the reference sequence and a primer sequence results in a non-amplification event. That is to say, if one primer cannot 'anneal', then that amplicion will not amplified for that primer pair and there will be no defline in the ```SEWAGE_amplicons.fasta``` file.  However, the ```SEWAGE_amplicon.tsv``` file will indicate which primers did not produce amplicons by stating the primer name followed by "No Amplification" in the **amplicon_sequence** column. Generating amplcions is perfrom similarly to the [in_silico_PCR.pl](https://github.com/egonozer/in_silico_pcr) tool.  

### Comments about generating reads
Default parameters when using the ```minimal usage``` command have been optimized to result in 150bp F/R reads form 300bp fragments at a total depth of coverage of approximatly 500X between all reference sequences supplied (NOTE: estimated depth of coverage will result in higher coverage due to overlaping regions of tiled amplicons when using defult settings). When **Read generator options** are adjusted, users are recommended to manually inspect values in the **SEWAGE_metadata.tsv** file for the proporitons, total number of fragments generated per ampicon, the length of reads, etc... and any other *in silico* data before running any downstream analyses.  

### Additional Notes 
**NOTE:** Only "short-read" fastq data is available as of now. I am activly working on adding long read output data. There is probaly some way to modify the paramters to achieve amplicon long reads using the **vsl1a** scheme. If so please let me know!  

### SEWAGE is always a WIP (Work In Progress)

### Future Additions:
1. Variation in fragement and read length
3. Mutaitons and indel additions 
4. Long Read simulations
5. Conda package
6. Remove primers from amplicons prior generating reads  

## LICENSE
SEWAGE
Copyright (C) 2024 Dorian J. Feistel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
