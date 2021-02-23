from customFunctions import *
import itertools
import numpy as np

for j in range(1, 5):
    print("\n\n\nSplitting into " + str(j) + " chunks.")
    total_elements_from_chunks = 0
    for i in range(61,79,2):
        print("Number of elements: " + str(i) + "\n")
        data = [list(np.linspace(1.7,2.7,i)), list(np.linspace(1,2,1)), list(np.linspace(0.3,0.35,1))]
        data = itertools.product(*data)
        data = [list(fitting_factor_combination) for fitting_factor_combination in data]
        data = split_array(data, j)
        print("After running split_array(data, " + str(j) + "):")
        for k in range(0, len(data)):
            print("Chunk " + str(k) + " has " + str(len(data[k])) + " elements.")
            total_elements_from_chunks += len(data[k])
        print("For a total of " + str(total_elements_from_chunks) + " among all the chunks.")
        data = split_array_chunks(data,15)
        print("\nAfter running split_array_chunks(data,15):\n")
        for k in range(0, len(data)):
            print("Printing lengths of mini-chunks in chunk " + str(k) + ", which has " + str(len(data[k])) + " elements:")
            for m in range(0, len(data[k])):
                print("Mini-chunk " + str(m) + " has " + str(len(data[k][m])) + " elements.")