#!/usr/bin/env python3
import numpy as np
import random
import sys
from scipy.stats import truncnorm
import matplotlib.pyplot as plt



def get_truncated_normal(mean=0, sd=4, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def skewed_normal_distribution(read_size, mean:int=40, lower_limit:int=20, upper_limit:int=40, sd:int=5):
    
    
    # Create a truncated normal distribution
    tr_norm = get_truncated_normal(mean=mean, sd=sd, low=lower_limit, upp=upper_limit)

    # Generate values and round to get discrete numbers
    values = tr_norm.rvs(read_size)
    discrete_values = np.round(values)

    return discrete_values

def map_values_to_chars(array_2d):
    illumina_qscore_dict = {
    0: '!',
    1: '"',
    2: '#',
    3: '$',
    4: '%',
    5: '&',
    6: "'",
    7: '(',
    8: ')',
    9: '*',
    10: '+',
    11: ',',
    12: '-',
    13: '.',
    14: '/',
    15: '0',
    16: '1',
    17: '2',
    18: '3',
    19: '4',
    20: '5',
    21: '6',
    22: '7',
    23: '8',
    24: '9',
    25: ':',
    26: ';',
    27: '<',
    28: '=',
    29: '>',
    30: '?',
    31: '@',
    32: 'A',
    33: 'B',
    34: 'C',
    35: 'D',
    36: 'E',
    37: 'F',
    38: 'G',
    39: 'H',
    40: 'I'
}
    # Determine the maximum key value for the size of the lookup array
    max_index = max(illumina_qscore_dict.keys())
    
    # Create a lookup array with a size one greater than the max index
    # '<U1' assumes single-character strings; adjust if needed
    lookup_array = np.empty(max_index + 1, dtype='<U1')
    
    # Populate the lookup array based on the dictionary
    for key, value in illumina_qscore_dict.items():
        lookup_array[key] = value

    # Use the lookup array to convert the 2D array
    return lookup_array[array_2d]

def quailty_scores_generator_array(
        num_reads,
        read_length,
        max_q:int = 40,
        min_q:int = 20,
        begining_q:int=4,
        begining_bp:int=10
):
    '''workshoping this one
    currently creates a linear decrease in qvalues'''

    array_2d = np.tile(np.linspace(max_q ,min_q, read_length), (num_reads, 1))
    flux_2d = np.random.normal(0,5, (num_reads, read_length))
    
    final_array = array_2d + flux_2d
    final_array = np.clip(final_array, 0, 40) # qscores are between 0-40
    final_array = np.round(final_array)
    final_array[:, :begining_bp] -= begining_q
    return final_array
    # means = np.mean(final_array, axis=0)
    # std_devs = np.std(final_array, axis=0)

    # # Number of columns in the data
    # num_columns = final_array.shape[1]

    # # Plotting
    # plt.figure(figsize=(12, 6))

    # # Adjust the x values to match the number of columns
    # x_values = np.arange(1, num_columns + 1)

    # # Plotting the means and standard deviations
    # plt.errorbar(x=x_values, y=means, yerr=std_devs, fmt='o', ecolor='r', capsize=5)
    # plt.title("Mean and Standard Deviation of Each Column")
    # plt.xlabel("Column Number")
    # plt.ylabel("Mean Value")
    # plt.grid(True)

    # plt.show()
    


########



def generate_quality_scores(read_length, std_dev=5):
    # Starting quality score (near the maximum of 40)
    start_quality = 40
    end_quality = 20
    # Generate a gradual decrease in quality towards the end of the read
    quality_scores = np.linspace(start_quality, end_quality, read_length, endpoint=True)

    # Add random noise to introduce variation, controlled by std_dev
    random_noise = np.random.normal(0, std_dev, read_length)
    quality_scores = quality_scores + random_noise

    # Ensure quality scores are within the bounds (0 to 40)
    quality_scores = np.clip(quality_scores, 0, 40)

    # Round the scores to the nearest integer and convert to a list
    quality_scores = np.round(quality_scores).astype(int).tolist()

    return quality_scores

    

if __name__ == "__main__":

    scores_numbers = quailty_scores_generator_array()
    scores_q = map_values_to_chars(array=scores_numbers, mapping=illumina_qscore_dict)
    print(scores_q)
    