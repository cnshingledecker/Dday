import numpy as np
import itertools
import math
from customFunctions import split_array

data = [list(np.linspace(1.7,2.7,4)), list(np.linspace(1,2,4)), list(np.linspace(0.3,0.35,4))]
data = itertools.product(*data)
data = [list(j) for j in np.array_split([list(k) for k in data],4)]
data = split_array(data, 15)
i = -1
for node_chunks in data:
    j = -1
    i += 1
    for mini_chunk in node_chunks:
        j += 1
        print("\nFor node " + str(i) + ", chunk " + str(j) + " is , with length " + str(len(data[i][j])))