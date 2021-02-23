# This file tests: The amount of time it takes to generate a list of 3 numpy.linspaces using 2 different methods

import numpy as np
import time

# No more than 10 to 12 milliseconds' difference between the two ways, 
#     except when one the lists (not the ones with 10 and 1 element, respectively) is longer than 100 thousand elements; this increases the time
#     A runtime increase of 80 milliseconds was seen with the third list had 1 million elements

ts = time.time() * 1000
a = list(np.linspace(1.7,2.7,10))
b = list(np.linspace(1,2,11))
c = list(np.linspace(0.3,0.35,1000000))
j = [a,b,c]
time.sleep(0.5)
te = time.time() * 1000
print("Way 1 (creating the numpy linspaces and then creating a list of them) took " + str(te-ts) + " milliseconds.")

ts = time.time() * 1000
j = [list(np.linspace(0.3,0.35,1)), list(np.linspace(1,2,1)), list(np.linspace(1.7,2.7,10))]
time.sleep(0.5)
te = time.time() * 1000
print("Way 2 (creating the numpy linspaces inside the creation of the list) took " + str(te-ts) + " milliseconds.")