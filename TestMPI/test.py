from customFunctions import split_array
import itertools
import numpy as np

data = [list(np.linspace(1.7,2.7,2)), list(np.linspace(1,2,3)), list(np.linspace(0.3,0.35,13))]
data = itertools.product(*data)

# i = [j for j in range(1, 79)]
# print("Len of i is " + str(len(i)))
j = split_array([list(j) for j in data], 4)
for item in j:
    print("\nLength is " + str(len(item)))
    print(str(item))