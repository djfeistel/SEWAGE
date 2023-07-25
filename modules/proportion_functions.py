import numpy as np

def equalProportions(fasta_pathway_list:list)->list:
    equal_proportions = [1/len(fasta_pathway_list)] * len(fasta_pathway_list)
    fasta_proportions = zip(fasta_pathway_list, equal_proportions)
    print(f"Equal proportion is ~{round(1/len(fasta_pathway_list), 3)} for each FASTA file and is summed to {sum(equal_proportions)}")
    return list(fasta_proportions)

def randomProportions(fasta_pathway_list:list, random_seed:int=13)->list:
    '''return a list of float values that sum to 1 which are use as proportions'''
    num_proportions = len(fasta_pathway_list)
    np.random.seed(random_seed)
    random_proportions = np.random.rand(num_proportions)
    # Normalize the proportions to sum up to 1
    normalized_proportions =  random_proportions / np.sum(random_proportions)
    fasta_proportions = zip(fasta_pathway_list, normalized_proportions)
    #print(f"Proportion summed to {np.sum(normalized_proportions)}", file=sys.stderr)
    return list(fasta_proportions)