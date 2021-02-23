# This file tests: The time taken to create the combinations of fitting factors using the itertools product function 
#                  (with 2 different methods), and splitting up these combinations using the numpy function array_split 

import itertools, time
import numpy as np

a = list(np.linspace(1.7,2.7,10))
b = list(np.linspace(1,2,1))
c = list(np.linspace(0.3,0.35,2))

ts1 = time.time() * 1000
data = [a,b,c]
data = itertools.product(*data)
l = [list(j) for j in np.array_split([list(j) for j in data],6)]
time.sleep(0.5)
te1 = time.time() * 1000
print("Way 1 (itertools.product(*data)) took " + str(te1 - ts1) + " milliseconds.")

ts2 = (time.time() * 1000)
data = itertools.product(*[a,b,c]) # This way takes longer as the number of arrays in data increases, but is comparable in terms of time for smaller test sets
l = [list(j) for j in np.array_split([list(j) for j in data],6)]
time.sleep(0.5)
te2 = (time.time() * 1000)
print("Way 2 (itertools.product(*[a,b,c])) took " + str(te2 - ts2) + " milliseconds.")