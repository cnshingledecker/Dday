import csv, os, time

matrix_dimensions = []
with open('num_matrices.csv', 'r') as file2:
    reader = csv.reader(file2, delimiter=",")
    matrix_dimensions = [int(line[0]) for line in reader]
    file2.close()

j = []

for matrix_dimension in matrix_dimensions:
    with open('num_matrices_iteration.csv', 'w') as write_file:
        write_file.write(str(matrix_dimension))
        write_file.close()

    ts1 = time.time() * 1000
    os.system("mpiexec -n 4 python3 broadcastSplit.py")
    te1 = time.time() * 1000

    ts2 = time.time() * 1000
    os.system("python3 nonParallelBroadcastSplit.py")
    te2 = time.time() * 1000
    j.append("Matrix dim: " + str(matrix_dimension) + " "*(6-len(str(matrix_dimension))) + " Time: Non-par: " + (str(te2-ts2)) + " ms" + " "*(40-len(" Time: Non-par: " + (str(te2-ts2)) + " ms")) + " Par: " + str(te1-ts1) + " ms   Runtime ratio p/np: " + str((te1-ts1)/(te2-ts2)))

for s in j:
    print(s)