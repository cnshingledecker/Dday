import csv, os, time
import numpy as np
from nonParallelSumBroadcastSplit import sum_num

numpy_vars = []
with open('sum_linspace_vars.csv', 'r') as file2:
    reader = csv.reader(file2, delimiter=",")
    for line in reader:
        numpy_vars.append([int(line[0]), int(line[1]), int(line[2])])
    file2.close()

for j in range(1,5):
    sums_non_parallel = []
    sums_parallel = []

    with open('sum_num_processors.csv', 'w') as num_processors_file:
        num_processors_file.write(str(j))
        num_processors_file.close()

    for numpy_var_set in numpy_vars:
        with open('sum_iter_vars.csv', 'w') as write_file:
            write_file.write(str(numpy_var_set[0]) + "," + str(numpy_var_set[1]) + "," + str(numpy_var_set[2]))
            write_file.close()
        sums_non_parallel.append(sum_num(numpy_var_set[0], numpy_var_set[1], numpy_var_set[2]))
        os.system("mpiexec -n " + str(j) + " python3 sumBroadcastSplit.py")
        with open('sum_parallel.csv', 'r') as file2:
            reader = csv.reader(file2, delimiter=",")
            for line in reader:
                sums_parallel.append((line[0]))
            file2.close()

    print("\n\n\nResults for running the code with " + str(j) + " processors:\n")

    if(len(sums_non_parallel) != len(sums_parallel)):
        print("An error occurred with the parallel sum calculator: \n")

    for i in range(0, len(sums_non_parallel)):
        # print("J: " + str(j) + " I: " + str(i))
        print("For the np.linspace parameters " + str(numpy_vars[i]) + " "*(44-len("For the np.linspace parameters " + str(numpy_vars[i]))) + ", the difference between the parallel and non-parallel sums was " + str(np.float64(sums_parallel[i])) + " "*(18-len((str(np.float64(sums_parallel[i]))))) + " - " + str(sums_non_parallel[i]) + " "*(20-len(str(sums_non_parallel[i]))) +   " = " + str(np.float64(sums_parallel[i]) - sums_non_parallel[i]))