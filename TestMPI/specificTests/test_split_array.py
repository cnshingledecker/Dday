from customFunctions import split_array
import numpy as np
import itertools

for j in range(1, 10):
    print("Splitting into " + str(j) + " chunks.")
    for i in range(500,1000):
        data = [list(np.linspace(1.7,2.7,i)), list(np.linspace(1,2,1)), list(np.linspace(0.3,0.35,1))]
        data = itertools.product(*data)
        data = [list(fitting_factor_combination) for fitting_factor_combination in data]
        print("The number of big chunks is " + str(j) +", and length of array is " + str(len(data)))
        data = split_array(data, j)