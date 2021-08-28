# This file does: For 1-4 cores, it runs nonParallelSumBroadcastSplit and sumBroadcastSplit (the parallel version),
#                 having them each sum up certain numbers determined by the number of elements in linspaces (which come from the sum_linspace_vars.csv file)

import csv, os, time
import numpy as np
from nonParallelSumBroadcastSplit import sum_num

numpy_vars = []
with open('sum_linspace_vars.csv', 'r') as file2:
    reader = csv.reader(file2, delimiter=",")
    for line in reader: # Add all number of elements in the numpy linspaces to a list
        numpy_vars.append([int(line[0]), int(line[1]), int(line[2])])
    file2.close()

for j in range(1,5): # J is the number of cores used in the parallel version
    sums_non_parallel = []
    sums_parallel = []

    with open('sum_num_cores.csv', 'w') as num_cores_file: # Writes the number of cores to be used by and in running the parallel sumBroadcastSplit script
        num_cores_file.write(str(j))
        num_cores_file.close()

    for numpy_var_set in numpy_vars: 
        with open('sum_iter_vars.csv', 'w') as write_file: # Write numpy linspace parameters to a file that the parallel sumBroadcastSplit script reads from and uses
            write_file.write(str(numpy_var_set[0]) + "," + str(numpy_var_set[1]) + "," + str(numpy_var_set[2]))
            write_file.close()
        sums_non_parallel.append(sum_num(numpy_var_set[0], numpy_var_set[1], numpy_var_set[2])) # Find sum with non-parallel method and add it to a list
        os.system("mpiexec -n " + str(j) + " python3 sumBroadcastSplit.py")
        with open('sum_parallel.csv', 'r') as file2: # Read in sum result (from parallel version) and add it to a list
            reader = csv.reader(file2, delimiter=",")
            for line in reader: # There is only 1 line in the file
                sums_parallel.append((line[0]))
            file2.close()

    print("\n\n\nResults for running the code with " + str(j) + " cores:\n")

    if(len(sums_non_parallel) != len(sums_parallel)):
        print("An error occurred with the parallel sum calculator: \n")

    for i in range(0, len(sums_non_parallel)):
        print("For the np.linspace parameters " + str(numpy_vars[i]) + " "*(48-len("For the np.linspace parameters " + str(numpy_vars[i]))) + ", the difference between the parallel and non-parallel sums was " + str(np.float64(sums_parallel[i])) + " "*(18-len((str(np.float64(sums_parallel[i]))))) + " - " + str(sums_non_parallel[i]) + " "*(20-len(str(sums_non_parallel[i]))) +   " = " + str(np.float64(sums_parallel[i]) - sums_non_parallel[i]))