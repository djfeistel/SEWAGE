#!/usr/bin/env python3
import numpy as np

mean = 1000
size = 10

exponential_numbers = np.random.exponential(scale=1/mean, size=size)
poisson_numbers = np.random.poisson(mean, size=size)
poisson_props = poisson_numbers / sum(poisson_numbers)
print(poisson_props, sum(poisson_props))