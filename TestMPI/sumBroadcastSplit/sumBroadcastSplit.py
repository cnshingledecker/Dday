import numpy as np
import csv, itertools, math, os, random, time

from customFunctions import split_array, split_array_chunks
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    num_processors = 0
    with open('sum_num_processors.csv', 'r') as num_processors_file:
        reader = csv.reader(num_processors_file, delimiter=",")
        for line in reader:
            num_processors = int(line[0])
        num_processors_file.close()
    num1 = 0
    num2 = 0
    num3 = 0
    with open('sum_iter_vars.csv', 'r') as csv_file:
        reader = csv.reader(csv_file)
        for line in reader:
            num1 = int(line[0])
            num2 = int(line[1])
            num3 = int(line[2])
        csv_file.close()
    data = [list(np.linspace(1.7,2.7,num1)), list(np.linspace(1,2,num2)), list(np.linspace(0.3,0.35,num3))]
    data = itertools.product(*data)
    data = [list(fitting_factor_combination) for fitting_factor_combination in data]
    data = split_array(data, num_processors)
    data = split_array_chunks(data, 15)
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
    num_mini_chunks_to_recv = comm.recv(source=0)
    processor_fitting_factor_combinations = []
    sum_total = 0
    for i in range(0, num_mini_chunks_to_recv):
        data = comm.recv(source=0)
        processor_fitting_factor_combinations.append(data)
    num_samples_processed = 0
    for mini_chunk in processor_fitting_factor_combinations:
        for combination in mini_chunk:
            sum_total += sum(combination)
    # print("Sum for rank " + str(rank) + " is " + str(sum_total))

sums = comm.gather(sum_total, root=0)

if rank == 0:
    with open('sum_parallel.csv', 'w') as file_write:
        file_write.write(str(sum(sums)))
        file_write.close()