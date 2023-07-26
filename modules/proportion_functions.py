import numpy as np
import random

def equalProportions(fasta_pathway_list:list)->list:
    equal_proportions = [1/len(fasta_pathway_list)] * len(fasta_pathway_list)
    fasta_proportions = zip(fasta_pathway_list, equal_proportions)
    print(f"Equal proportion is ~{round(1/len(fasta_pathway_list), 3)} for each FASTA file and is summed to {sum(equal_proportions)}")
    return list(fasta_proportions)

def randomProportions(fasta_pathway_list: list, random_seed: int = 13) -> list:
    '''
    return a list of float values that sum to 1 which are use as proportions
    no value value can be less thant 0.001
    '''
    num_proportions = len(fasta_pathway_list)
    np.random.seed(random_seed)

    random_proportions = np.random.rand(num_proportions)

    normalized_proportions = random_proportions / np.sum(random_proportions)
    
    idx = normalized_proportions < 0.001

    normalized_proportions[idx] += 0.001
    
    normalized_proportions /= np.sum(normalized_proportions)

    fasta_proportions = zip(fasta_pathway_list, normalized_proportions)

    return list(fasta_proportions)


def dominantVariantProportions(fasta_pathway_list:list, random_seed:int=13, dVOC:float=0.8):
    # Generate a list of random numbers between 0 and 1
    np.random.seed(random_seed)
    random_list = [random.uniform(0.001, 1) for _ in range(len(fasta_pathway_list))]
    
    # Choose one element from the random list and set it to 0.8
    random_index = random.randint(0, len(random_list) - 1)
    random_list[random_index] = dVOC
    
    # Adjust the other elements to maintain the sum of the list as 1
    total_sum = sum(random_list)
    remaining_sum = 1.0 - dVOC
    scaling_factor = remaining_sum / (total_sum - dVOC)
    for i in range(len(random_list)):
        if i != random_index:
            random_list[i] *= scaling_factor
    random.shuffle(random_list)
    fasta_proportions = list(zip(fasta_pathway_list, random_list))
    return fasta_proportions
