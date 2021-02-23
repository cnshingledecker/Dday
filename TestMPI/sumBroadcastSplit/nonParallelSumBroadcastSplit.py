# This file does: Creates 3 numpy.linspaces and their Cartesian product, adds the sum of each combination to a sum variable, and returns that sum variable

import numpy as np, itertools

def sum_num(a,b,c):
    sum_total = 0
    data = [list(np.linspace(1.7,2.7,a)), list(np.linspace(1,2,b)), list(np.linspace(0.3,0.35,c))] # Creates the numpy linspaces to be usedin making the cartesian product (all possible fitting factor combinations)
    data = itertools.product(*data) # Creates the cartesian product of all possible fitting factor combinations
    data = [list(fitting_factor_combination) for fitting_factor_combination in data] # Converts the result of itertools.product into a list (completely array-indexible)
    
    for k in data:
        sum_total += sum(k)

    return sum_total

