# This file tests: 

import numpy as np
import csv, itertools, math, os, random, time
from customFunctions import split_array, split_array_chunks
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

matrix_dimension = 0
with open('num_matrices_iteration.csv', 'r') as file2: # This csv file holds the N dimension (A NxN matrix will be created further down in the program)
    reader = csv.reader(file2, delimiter=",")
    for line in reader: # There will be only 1 line
        matrix_dimension = int(line[0])

if rank == 0:
    num_processors = 4
    data = [list(np.linspace(1.7,2.7,23)), list(np.linspace(1,2,15)), list(np.linspace(0.3,0.35,1))] # Creates the numpy linspaces to be usedin making the cartesian product (all possible fitting factor combinations)
    data = itertools.product(*data) # Creates the cartesian product of all possible fitting factor combinations
    data = [list(fitting_factor_combination) for fitting_factor_combination in data] # Converts the result of itertools.product into a list (completely array-indexible)
    data = split_array(data, num_processors) # Splits the array into 'num_processors' chunks
    data = split_array_chunks(data, 15) # Splits each of the array chunks into mini-chunks of up to size 15 (creates as many with size 15 as possible)
    
    files_to_copy_to_new_dir = ["stlCardinals.txt", "alarmClock.dat", "drawer.csv"] # Note; this can be used with copying files into the directory when writing the parallel version of the generalized baragiola optimization

    for i in range(0, num_processors): # Create directory for the files for each processor, copy into it the files specified in the above array, and create the results file in each array
        new_dir_name = "files_processor" + str(i)
        os.system("mkdir " + new_dir_name)
        for file_name in files_to_copy_to_new_dir:
            os.system("cp " + file_name + " " + new_dir_name)
        os.system("touch " + new_dir_name + "/results_file")
    for i in range(0, len(data)): # Tell each processor how many mini-chunks it is going to be receiving
        comm.send(len(data[i]), dest=i)
    mini_chunk_to_send = [0 for i in range(num_processors)] # Tells comm.send which mini-chunk to send (tells the index) to a certain processor
    chunks_left_to_send = [True for i in range(num_processors)] # List of boolean variables that says whether there are mini-chunks left to send to each processor
    while(sum(chunks_left_to_send) > 0):
        for j in range(0, len(data)):
            if mini_chunk_to_send[j] < len(data[j]):
                comm.send(data[j][mini_chunk_to_send[j]], dest=j)
                mini_chunk_to_send[j] += 1
            else:
                chunks_left_to_send[j] = False

if rank >= 0:
    random.seed()
    num_mini_chunks_to_recv = comm.recv(source=0) # Receive the number of mini-chunks it is going to receive
    processor_fitting_factor_combinations = []
    fitting_factors_and_least_rmsd = [1e80, 0,0,0]
    reactions = ["Reaction 1", "Reaction 2", "Reaction 3"] # Used in writing output string of fitting factors and the fake_performance_measure to output file

    results = open("files_processor" + str(rank) + "/results_file",'w')
    num_items_recv = 0
    for i in range(0, num_mini_chunks_to_recv): # Receive all the mini-chunks
        data = comm.recv(source=0)
        num_items_recv += len(data)
        processor_fitting_factor_combinations.append(data)
    num_samples_processed = 0
    for mini_chunk in processor_fitting_factor_combinations:
        num_samples_processed += len(mini_chunk)
        for combination in mini_chunk:
            for i in range(100): # This is to simulate the work done by the processor (and the time it takes) when running monaco
                j = np.random.rand(matrix_dimension,matrix_dimension)
            # For the parallel generalized baragiola optimization, will do work of inserting new parameters; running the run.sh file for this processor (includes monaco), and calculation of a performance measure  
            fake_performance_measure = math.sqrt(sum([fitting_factor**2 for fitting_factor in combination]))
            if (fake_performance_measure < fitting_factors_and_least_rmsd[0]):
                fitting_factors_and_least_rmsd[0] = fake_performance_measure
                fitting_factors_and_least_rmsd[1] = combination[0]
                fitting_factors_and_least_rmsd[2] = combination[1]
                fitting_factors_and_least_rmsd[3] = combination[2]
            output_string  = ""
            for i in range(0, len(reactions)): # Write to the output file for this processor (write reactions and associated fitting factors and the fake performance metric)
                output_string = output_string + str(combination[i]) + "".join(" "*(23 - len(str(combination[i])))) + reactions[i] + " delta values \n"
            output_string += str(fake_performance_measure) + "".join(" "*(23 - len(str(fake_performance_measure)))) + "RMSD" + "\n\n"
            results.write(output_string)
        # print("The number of samples processed is " + str(num_samples_processed))
    results.close()

least_rmsds_and_fitting_factors = comm.gather(fitting_factors_and_least_rmsd, root=0) # Gather the least fake performance metric value and associated fitting factors from each processor back to the root processor (0)

if rank == 0:
    least_rmsd_index = 0
    for i in range(1, len(least_rmsds_and_fitting_factors)): # Loop through the least fake performance metric value and associated fitting factors 
                                                             # from each processor and find the ones with the least value for the fake performance metrix
        if least_rmsds_and_fitting_factors[i][0] < least_rmsds_and_fitting_factors[least_rmsd_index][0]:
            least_rmsd_index = i
    # os.system("rm -r files_processor*")

    print("The smallest performance metric value and associated fitting factors: ")
    print(least_rmsds_and_fitting_factors[least_rmsd_index])
    
    # Optional test to verify that the fake performance metric value is the square root of the sum of the squares of its associated fitting factors
    # for array in least_rmsds_and_fitting_factors:
    #     fake_performance_measure = 0
    #     for i in range(1,4):
    #         fake_performance_measure += array[i]**2
    #     assert array[0] == math.sqrt(fake_performance_measure)