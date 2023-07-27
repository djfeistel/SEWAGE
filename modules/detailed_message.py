def details():
    message = '''
    SEWAGE: Synthetically Enriched Wastewater sequence data for Assessing Genomic variants of Epidemiological surveilance

    SEWAGE is a tool for generating reproducible sequence data representing a heterogeneous population of 
    closely related species. Specifically, it was designed to mirror sequence data that resembles a mixed 
    SARS-CoV-2 population derived from a wastewater sample by using targeted enrichment or tiled amplicon 
    approaches. 

    SEWAGE currently offers two main functionalities: 1) the ability to produce amplicons for each genome 
    from a set of closely related reference genomes using a set of primers and 2) create Illumina short-read 
    data sets from those amplicions that mimic heterogeneous populations of closely related species at 
    various proportions.  SEWAGE currently only offers SARS-CoV2 ARTIC 
    (https://github.com/artic-network/primer-schemes) and VarSkip (https://github.com/nebiolabs/VarSkip) 
    primer sets for creating amplicions. However, we are currently working on allowing users to supply 
    their own uniqe primer sets for generating amplicions. 

    SEWAGE was developed out of necessity to create reproducible SARS-CoV-2 mixed population in silico sequence data to
    assess the accurcy of alignment/mapper tools as well as relative abundance calculation 
    tools for the National Wastewater Survaillance Systems (NWSS) SARS-CoV-2 Wastewater detection 
    pipeline titled "AquaScope" (CIT GITHUB/GITLAB).

    The idea for SEWAGE is to create enriched wastewater-like sequence data that are "ideal", 
    meaning that the simulated reads do not deviated from the genomes or amplicons that are supplied 
    as input. The reason we say "enriched wastewater-like" is due to how SARS-CoV-2 has been sequenced
    from wastewater (WW) samples. Briefly, SARS-CoV-2 is typically (but not exclusivly) amplified from 
    a WW sample by tiled amplicon or target enrichment sequencing, and those amplicons are then 
    sequenced using some desiered platform.  

    Currently, SEWAGE's default settings create Illumina 150bp forward (F) and reverse (R) reads and 
    heavily rely on the bioinformatics read simulation tool "art_illumina" (Huang, Weichun, et al. 
    "ART: a next-generation sequencing read simulator." Bioinformatics 28.4 (2012): 593-594.) to create 
    enriched wastewater-like sequence data. However, in the future we plan to update the tool with our 
    own simulation algorithm so that users do not need to install outside tools. Furthermore, we plan 
    to update SEWAGE with other functionallity which (should) include functions like amplicon and long 
    read sequencing and addition of barcodes to sequences (for more 'realworld' sequencing approachs).
    Additionally, we hope to create a conda package.

    Some notes about using SEWAGE:
    Leaving the majority of the paramters on their defult setting should result in reads that are 
    "perfect" (i.e. no indels or mutations), matching the parent sequence .  If you want to introduce 
    more variation in to the F/R reads, please familiarize yourself with the "art_illumina" tools 
    before manipulating the parameters until we add in out internal sequencing algorithm.

    Dependencies:
    python 3.8.3
    tqdm (pip install tqdm)
    ART (conda install -c bioconda art)

    Minimal usage: SEWAGE -i <input>

    SEWAGE -h for help menu
    
    https://github.com/djfeistel/SEWAGE
    '''
    return print(message)