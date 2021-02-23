import numpy as np, itertools

def sum_num(a,b,c):
    sum_total = 0
    data = [list(np.linspace(1.7,2.7,a)), list(np.linspace(1,2,b)), list(np.linspace(0.3,0.35,c))]
    data = itertools.product(*data)
    data = [list(fitting_factor_combination) for fitting_factor_combination in data]

    for k in data:
        sum_total += sum(k)

    return sum_total

