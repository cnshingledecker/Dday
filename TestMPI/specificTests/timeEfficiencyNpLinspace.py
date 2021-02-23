# This file tests: Whether creating numpy linspaces and then converting them to lists or combining the two operations is faster

import numpy as np
import time

# Between way 1 and 2, doesn't appear to be much time differences, even for a large number of arrays to split

ts = time.time() * 1000
a = np.linspace(1.7,2.7,10)
b = np.linspace(1,2,30)
c = np.linspace(0.3,0.35,100000)
a = list(a)
b = list(b)
c = list(c)
time.sleep(0.5)
te = time.time() * 1000
print("Way 1 (creating the numpy linspaces and then converting them to lists) took " + str(te-ts) + " milliseconds.")

ts = time.time() * 1000
a = list(np.linspace(1.7,2.7,10))
b = list(np.linspace(1,2,1))
c = list(np.linspace(0.3,0.35,1))
time.sleep(0.5)
te = time.time() * 1000
print("Way 2 (creating the numpy linspaces and converting them to lists inside a conversion of them to lists) took " + str(te-ts) + " milliseconds.")