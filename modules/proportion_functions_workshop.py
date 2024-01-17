import numpy as np
import random
import sys
from fasta_funcitons_class_workshop import load_fasta_workflow_test


'''
calculate various proportions values for reference genomes for determining relative abundance in fastq files
'''
class GenomeProporitons():

    def __init__(self,
                reference_fasta_dict: dict
                ):
        self.reference_fasta_dict = reference_fasta_dict

    def equalProportions(self)->list:
        reference_genomes_keys = self.reference_fasta_dict.keys()
        equal_proportions = [1/len(reference_genomes_keys)] * len(reference_genomes_keys)
        return dict(zip(reference_genomes_keys, equal_proportions))

    def randomProportions(self, random_seed: int=13) -> list:
        '''
        return a list of float values that sum to 1 which are use as proportions
        no value value can be less thant 0.001
        '''
        np.random.seed(random_seed)
        # shuffle genomes and create list
        reference_genomes_keys = list(self.reference_fasta_dict.keys())
        np.random.shuffle(reference_genomes_keys)
        # count genomes
        reference_genomes_count = len(reference_genomes_keys)
        # generate array of random proportions
        random_proportions = np.random.rand(reference_genomes_count)
        # normalize proportions 
        normalized_proportions = random_proportions / np.sum(random_proportions)
        # index proporions less than 0.01
        idx = normalized_proportions < 0.01
        # add 0.01 to those values < 0.01
        normalized_proportions[idx] += 0.01
        # renormalize proportions
        normalized_proportions /= np.sum(normalized_proportions)
        # create tuple
        return dict(zip(reference_genomes_keys, normalized_proportions))

    def dominantVariantProportions_og(self, dVOC_genome:str=None, dVOC_proporiton:float=0.8, random_seed:int=13):
        # Generate a list of random numbers between 0 and 1
        np.random.seed(random_seed)

        reference_genomes_keys = list(self.reference_fasta_dict.keys())
        
        np.random.shuffle(reference_genomes_keys)
        
        random_list = [random.uniform(0.001, 1) for _ in range(len(reference_genomes_keys))]
        
        
        # Choose one element (i.e., defline) from the random list and set it to 0.8
        random_index = random.randint(0, len(random_list) - 1)
        random_list[random_index] = dVOC_proporiton

        # Adjust the other elements to maintain the sum of the list as 1
        total_sum = sum(random_list)
        remaining_sum = 1.0 - dVOC_proporiton
        scaling_factor = remaining_sum / (total_sum - dVOC_proporiton)
        for i in range(len(random_list)):
            if i != random_index:
                random_list[i] *= scaling_factor
            for i in range(len(random_list)):
                if random_list[i] < 0.001:
                    random_list[i] += 0.001

        # Normalize the list to ensure the sum of the values equals 1 again
        total_sum = sum(random_list)
        random_list = [val / total_sum for val in random_list]

        random.shuffle(random_list)
        return dict(zip(reference_genomes_keys, random_list))

    def dominantVariantProportions(self, dVOC_genome: str = None, dVOC_proporiton: float = 0.8, random_seed: int = 13):
        # Set random seed for reproducibility
        np.random.seed(random_seed)
        random.seed(random_seed)

        reference_genomes_keys = list(self.reference_fasta_dict.keys())
        np.random.shuffle(reference_genomes_keys)

        # Check if dVOC_genome is provided and in reference_genomes_keys
        if dVOC_genome and dVOC_genome in reference_genomes_keys:
            # Assign the dVOC_proportion to the dVOC_genome
            proportions = {genome: 0 for genome in reference_genomes_keys}
            proportions[dVOC_genome] = dVOC_proporiton

            # Calculate the remaining proportion to distribute among other genomes
            remaining_sum = 1.0 - dVOC_proporiton
            other_genomes = [genome for genome in reference_genomes_keys if genome != dVOC_genome]

            # Generate random proportions for other genomes
            random_proportions = [random.uniform(0.001, 1) for _ in other_genomes]
            scaling_factor = remaining_sum / sum(random_proportions)
            scaled_proportions = [prop * scaling_factor for prop in random_proportions]

            # Ensure no proportion is less than 0.001
            for i in range(len(scaled_proportions)):
                if scaled_proportions[i] < 0.001:
                    scaled_proportions[i] = 0.001

            # Update the proportions dictionary with scaled values for other genomes
            for genome, prop in zip(other_genomes, scaled_proportions):
                proportions[genome] = prop

            # Normalize the proportions to sum up to 1
            total_sum = sum(proportions.values())
            proportions = {genome: prop / total_sum for genome, prop in proportions.items()}

        else:
            # If dVOC_genome is not provided or not in the list, apply the original logic
            random_list = [random.uniform(0.001, 1) for _ in range(len(reference_genomes_keys))]
            random_index = random.randint(0, len(random_list) - 1)
            random_list[random_index] = dVOC_proporiton
            total_sum = sum(random_list)
            remaining_sum = 1.0 - dVOC_proporiton
            scaling_factor = remaining_sum / (total_sum - dVOC_proporiton)
            for i in range(len(random_list)):
                if i != random_index:
                    random_list[i] *= scaling_factor
                for i in range(len(random_list)):
                    if random_list[i] < 0.001:
                        random_list[i] += 0.001
            total_sum = sum(random_list)
            random_list = [val / total_sum for val in random_list]
            random.shuffle(random_list)
            proportions = dict(zip(reference_genomes_keys, random_list))

        return proportions


def proportion_workflow_test():
    reference_fasta_file = load_fasta_workflow_test(sys.argv[1])
    genome_props = GenomeProporitons(reference_fasta_dict=reference_fasta_file)
    dVOC_props = genome_props.dominantVariantProportions()
    print(dVOC_props)
    

if __name__ == "__main__":
    proportion_workflow_test()
