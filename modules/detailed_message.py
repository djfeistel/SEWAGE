def details():
    message = '''
    SEWAGE: Synthetically Enriched Wastewater-like sequence data for Assessing Genomic and Environmental populations

    SEWAGE was originally created out of a need to create reproducable in silico sequence data to
    assess the accurcy of different alignment/mapper tools as well as relative abundance calculation 
    tools for the National Wastewater Survaillance Systems (NWSS) SARS-CoV-2 Wastewater detection 
    pipeline titled "AquaScope" (CIT GITHUB/GITLAB).

    The general idea for SEWAGE is to create enriched wastewater-like sequence data that are "ideal", 
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
    "perfect" i.e. matching the parent sequence (no indels or mutations).  If you want to introduce 
    more variation in to the F/R reads, please familiorize yourself with the "art_illumina" tools 
    before manipulating the parameters.

    Dependencies:
    python 3.8.3
    tqdm (pip install tqdm)
    ART (conda install -c bioconda art)

    Minimal usage: SEWAGE -i <input>

    SEWAGE -h for help menu
    
    '''
    return print(message)