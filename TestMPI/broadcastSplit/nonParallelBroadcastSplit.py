import numpy as np
import csv, itertools, time, os
from customFunctions import *


num_cores = 4
data = [list(np.linspace(1.7,2.7,23)), list(np.linspace(1,2,15)), list(np.linspace(0.3,0.35,1))]
data = itertools.product(*data)
data = [list(fitting_factor_combination) for fitting_factor_combination in data]
fitting_factors_and_least_rmsd = [1e80,0,0,0]
reactions = ["Reaction 1", "Reaction 2", "Reaction 3"]
results = open("results_file",'w')

matrix_dimension = 0
with open('num_matrices_iteration.csv', 'r') as file2:
    reader = csv.reader(file2, delimiter=",")
    for line in reader:
        matrix_dimension = int(line[0])

for k in data:
    for i in range(100):
        j = np.random.rand(matrix_dimension,matrix_dimension)
    fake_performance_measure = math.sqrt(sum([fitting_factor**2 for fitting_factor in k]))
    if (fake_performance_measure < fitting_factors_and_least_rmsd[0]):
        fitting_factors_and_least_rmsd[0] = fake_performance_measure
        fitting_factors_and_least_rmsd[1] = k[0]
        fitting_factors_and_least_rmsd[2] = k[1]
        fitting_factors_and_least_rmsd[3] = k[2]
    output_string  = ""
    for i in range(0, len(reactions)): 
        output_string = output_string + str(k[i]) + "".join(" "*(23 - len(str(k[i])))) + reactions[i] + " delta values \n"
    output_string += str(fake_performance_measure) + "".join(" "*(23 - len(str(fake_performance_measure)))) + "RMSD" + "\n\n"
    results.write(output_string)

results.close()

least_rmsd_index = 0
for i in range(1, len(data)):
    if data[i][0] < data[least_rmsd_index][0]:
        least_rmsd_index = i

print("The smallest performance metric value and associated fitting factors: ")
print(fitting_factors_and_least_rmsd)
print(data[least_rmsd_index])
