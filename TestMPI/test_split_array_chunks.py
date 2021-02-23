from broadcastSplit.customFunctions import *
import numpy as np
import itertools

for j in range(1, 5):
    print("Splitting into " + str(j) + " chunks.")
    for i in range(61,79,2):
        print("\nNumber of elements: " + str(i))
        data = [list(np.linspace(1.7,2.7,i)), list(np.linspace(1,2,1)), list(np.linspace(0.3,0.35,1))]
        data = itertools.product(*data)
        data = [list(fitting_factor_combination) for fitting_factor_combination in data]
        data = split_array(data, j)
        data = split_array_chunks(data,15)