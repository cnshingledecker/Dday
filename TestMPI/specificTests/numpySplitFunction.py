import numpy as np
import itertools
import time

a = list(np.linspace(1.7,2.7,10))
b = list(np.linspace(1,2,1))
c = list(np.linspace(0.3,0.35,2))

ts1 = time.time() * 1000
data = [a,b,c]
data = itertools.product(*data)
l = [list(j) for j in np.array_split([list(j) for j in data],6)]
time.sleep(0.5)
te1 = time.time() * 1000
print("Way 1 took " + str(te1 - ts1) + " milliseconds.")

ts2 = (time.time() * 1000)
# data = [a,b,c]
data = itertools.product(*[a,b,c]) # This way takes longer as the number of arrays in data increases, but is comparable in terms of time for smaller test sets
l = [list(j) for j in np.array_split([list(j) for j in data],6)]
time.sleep(0.5)
te2 = (time.time() * 1000)
print("Way 2 took " + str(te2 - ts2) + " milliseconds.")