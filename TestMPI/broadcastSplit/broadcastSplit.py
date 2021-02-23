import numpy as np
import csv, itertools, math, os, random, time

from ../customFunctions import split_array, split_array_chunks
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

matrix_dimension = 0
with open('num_matrices_iteration.csv', 'r') as file2:
    reader = csv.reader(file2, delimiter=",")
    for line in reader:
        matrix_dimension = int(line[0])

if rank == 0:
    num_processors = 4
    files_to_copy_to_new_dir = ["stlCardinals.txt", "alarmClock.dat", "drawer.csv"]
    data = [list(np.linspace(1.7,2.7,23)), list(np.linspace(1,2,15)), list(np.linspace(0.3,0.35,1))]
    data = itertools.product(*data)
    data = [list(fitting_factor_combination) for fitting_factor_combination in data]
    data = split_array(data, num_processors)
    data = split_array_chunks(data, 15)
    for i in range(0, num_processors):
        new_dir_name = "files_processor" + str(i)
        os.system("mkdir " + new_dir_name)
        for file_name in files_to_copy_to_new_dir:
            os.system("cp " + file_name + " " + new_dir_name)
        os.system("touch " + new_dir_name + "/results_file")
    for i in range(0, len(data)):
        comm.send(len(data[i]), dest=i)
    mini_chunk_to_send = [0 for i in range(num_processors)]
    chunks_left_to_send = [True for i in range(num_processors)]
    while(sum(chunks_left_to_send) > 0):
        for j in range(0, len(data)):
            if mini_chunk_to_send[j] < len(data[j]):
                comm.send(data[j][mini_chunk_to_send[j]], dest=j)
                mini_chunk_to_send[j] += 1
            else:
                chunks_left_to_send[j] = False

if rank >= 0:
    random.seed()
    num_mini_chunks_to_recv = comm.recv(source=0)
    processor_fitting_factor_combinations = []
    fitting_factors_and_least_rmsd = [1e80, 0,0,0]
    reactions = ["Reaction 1", "Reaction 2", "Reaction 3"]

    print("rank is " + str(rank) + " and trying to open " + "files_processor" + str(rank) + "/results_file")
    results = open("files_processor" + str(rank) + "/results_file",'w')
    num_items_recv = 0
    for i in range(0, num_mini_chunks_to_recv):
        data = comm.recv(source=0)
        num_items_recv += len(data)
        processor_fitting_factor_combinations.append(data)
    num_samples_processed = 0
    for mini_chunk in processor_fitting_factor_combinations:
        for combination in mini_chunk:
            for i in range(100):
                j = np.random.rand(matrix_dimension,matrix_dimension)
            # Do work of inserting new parameters; running the run.sh file for this processor (includes monaco), and calculation of a performance measure  
            fake_performance_measure = math.sqrt(sum([fitting_factor**2 for fitting_factor in combination]))
            if (fake_performance_measure < fitting_factors_and_least_rmsd[0]):
                fitting_factors_and_least_rmsd[0] = fake_performance_measure
                fitting_factors_and_least_rmsd[1] = combination[0]
                fitting_factors_and_least_rmsd[2] = combination[1]
                fitting_factors_and_least_rmsd[3] = combination[2]
            output_string  = ""
            for i in range(0, len(reactions)): 
                output_string = output_string + str(combination[i]) + "".join(" "*(23 - len(str(combination[i])))) + reactions[i] + " delta values \n"
            output_string += str(fake_performance_measure) + "".join(" "*(23 - len(str(fake_performance_measure)))) + "RMSD" + "\n\n"
            results.write(output_string)
    results.close()

least_rmsds_and_fitting_factors = comm.gather(fitting_factors_and_least_rmsd, root=0)

if rank == 0:
    least_rmsd_index = 0
    for i in range(1, len(least_rmsds_and_fitting_factors)):
        if least_rmsds_and_fitting_factors[i][0] < least_rmsds_and_fitting_factors[least_rmsd_index][0]:
            least_rmsd_index = i
    # os.system("rm -r files_processor*")

    print("The smallest performance metric value and associated fitting factors: ")
    print(least_rmsds_and_fitting_factors[least_rmsd_index])
    # for array in least_rmsds_and_fitting_factors:
    #     fake_performance_measure = 0
    #     for i in range(1,4):
    #         fake_performance_measure += array[i]**2
    #     assert array[0] == math.sqrt(fake_performance_measure)