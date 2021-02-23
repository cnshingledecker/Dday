from customFunctions import split_array
import numpy as np
import itertools

for j in range(1, 10): # J is the number of chunks the array is split into using split_array
    print("\n\n\nSplitting into " + str(j) + " chunks.")
    for i in range(500,1000,102): # I is the number of fitting factor combinations in the array
        data = [list(np.linspace(1.7,2.7,i)), list(np.linspace(1,2,1)), list(np.linspace(0.3,0.35,1))]
        data = itertools.product(*data)
        data = [list(fitting_factor_combination) for fitting_factor_combination in data]
        print("\nThe number of big chunks is " + str(j) +", and length of array is " + str(len(data)))
        total_elements_from_chunks = 0
        data = split_array(data, j)
        for m in range(0, len(data)):
            print("Chunk " + str(m) + " has length " + str(len(data[m])))
            total_elements_from_chunks += len(data[m])
        print("For a total of " + str(total_elements_from_chunks) + " elements among all the chunks.")