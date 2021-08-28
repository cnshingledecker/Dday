# This file does: Runs a parallelized version of nonParallelSumBroadcastSplit
#                 It splits a list of fitting factor combinations among 'num_cores' cores, has them
#                 sum up all of the fitting factor combinations they have, returns the sum to the root core (0),
#                 the root core adds up the sums from each of the cores and writes this sum to a file

import csv, itertools, math, os, random, time
import numpy as np

from customFunctions import split_array, split_array_chunks
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    num_cores = 0
    with open('sum_num_cores.csv', 'r') as num_cores_file:
        reader = csv.reader(num_cores_file, delimiter=",")
        for line in reader:
            num_cores = int(line[0])
        num_cores_file.close()
    num1 = 0
    num2 = 0
    num3 = 0
    with open('sum_iter_vars.csv', 'r') as csv_file: # This file holds 1 line, which contains the number of elements in each of the 3 numpy.linspaces created below 
        reader = csv.reader(csv_file)
        for line in reader:
            num1 = int(line[0])
            num2 = int(line[1])
            num3 = int(line[2])
        csv_file.close()
    data = [list(np.linspace(1.7,2.7,num1)), list(np.linspace(1,2,num2)), list(np.linspace(0.3,0.35,num3))] # Creates the numpy linspaces to be usedin making the cartesian product (all possible fitting factor combinations)
    data = itertools.product(*data) # Creates the cartesian product of all possible fitting factor combinations
    data = [list(fitting_factor_combination) for fitting_factor_combination in data] # Converts the result of itertools.product into a list (completely array-indexible)
    data = split_array(data, num_cores) # Splits the array into 'num_cores' chunks
    data = split_array_chunks(data, 15) # Splits each of the array chunks into mini-chunks of up to size 15 (creates as many with size 15 as possible)
    for i in range(0, len(data)): # Tell each core how many mini-chunks it is going to be receiving
        comm.send(len(data[i]), dest=i)
    mini_chunk_to_send = [0 for i in range(num_cores)] # Tells comm.send which mini-chunk to send (tells the index) to a certain core
    chunks_left_to_send = [True for i in range(num_cores)] # List of boolean variables that says whether there are mini-chunks left to send to each core
    while(sum(chunks_left_to_send) > 0):
        for j in range(0, len(data)):
            if mini_chunk_to_send[j] < len(data[j]):
                comm.send(data[j][mini_chunk_to_send[j]], dest=j)
                mini_chunk_to_send[j] += 1
            else:
                chunks_left_to_send[j] = False

if rank >= 0:
    num_mini_chunks_to_recv = comm.recv(source=0) # Receive the number of mini-chunks it is going to receive
    core_fitting_factor_combinations = []
    sum_total = 0
    for i in range(0, num_mini_chunks_to_recv): # Receive all the mini-chunks
        data = comm.recv(source=0)
        core_fitting_factor_combinations.append(data)
    num_samples_processed = 0
    for mini_chunk in core_fitting_factor_combinations:
        num_samples_processed += len(mini_chunk)
        for combination in mini_chunk:
            sum_total += sum(combination)
    # print("Sum for rank " + str(rank) + " is " + str(sum_total))

sums = comm.gather(sum_total, root=0) # Gather sum from each core to the root core (core 0)

if rank == 0:
    with open('sum_parallel.csv', 'w') as file_write: # Write to a file the sum of the resultant sums from each of the cores
        file_write.write(str(sum(sums)))
        file_write.close()