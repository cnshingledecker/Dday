import numpy as np
import csv, itertools, math, os, random, time
from exportable_custom_functions import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
reactions = []
all_vector_args = [] # Holds a series of arrays (one for each of the linspaces that is created for a reaction)

with open('reaction_fitting_factor_linspace_args/reaction_fitting_factor_vector_arguments.csv', newline='') as vector_creation_args_csv:  # Read in the parameters from the csv file for the creation of the linspaces (for each fitting factor to be varied)
    reader = csv.reader(vector_creation_args_csv, delimiter=',')      
    for i in range(0, 3): # Skips the first 3 lines of the csv file (lines which are comments)
        fields = next(reader)
    for row in reader:  # Each row is a set of arguments to be used to create the linspace for fitting factors for a reaction
        args_for_single_vector = []  # A vector to hold the set of arguments to be used to create the linspace for the fitting factors for a reaction
        for argument in row:
            if is_float(argument):
                args_for_single_vector.append(float(argument)) # Adds the arguments to the previously created vector (3 lines above)
            else:  # It is text (it is not an index specifying a delta value to choose; because it is text, it must be the reaction(s) for which the fitting factors are being varied using the linspace created using some of the values in this row)
                reactions.append(argument)
        all_vector_args.append(args_for_single_vector)

if rank == 0:
    fitting_factors = [] # Holds arrays containing the fitting_factors for the reactions (each element of the array is a list with the various fitting factors for a reaction)
    num_processors = 4
    for vector_args in all_vector_args:
        single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
        single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
        fitting_factors.append(single_vector_fitting_factors)     
    
    all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of fitting factor values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates an array of n arrays (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n arrays has k elements (where k is the number of fitting factors (num_modified_fitting_factors))
    all_fitting_factor_combinations = [list(fitting_factor_combination) for fitting_factor_combination in all_fitting_factor_combinations]
    all_fitting_factor_combinations = split_array(all_fitting_factor_combinations, num_processors) # Splits the array into 'num_processors' chunks
    all_fitting_factor_combinations = split_array_chunks(all_fitting_factor_combinations, 15) # Splits each of the array chunks into mini-chunks of up to size 15 (creates as many with size 15 as possible)
                                                                                              # Note: Sending a mini-chunk on my (Daniel's) machine does not arrive at the destination (the program just sits and the mini-chunk
                                                                                              #       never gets to its destination). Feel free to change this if desired and if your machine can handle a smaller or bigger mini-chunk size.
    files_to_copy_to_new_dir = ["TestREADME.md"] 

    for i in range(0, num_processors): # Create directory for the files for each processor, copy into it the files specified in the above array, and create the results file in each array
        new_dir_name = "baragiola_files_processor" + str(i)
        os.system("mkdir " + new_dir_name)
        os.system("touch " + new_dir_name + "/results_file")
        for file_name in files_to_copy_to_new_dir:
            os.system("cp " + file_name + " " + new_dir_name)
    for i in range(0, len(all_fitting_factor_combinations)): # Tell each processor how many mini-chunks it is going to be receiving
        comm.send(len(all_fitting_factor_combinations[i]), dest=i)
    mini_chunk_to_send = [0 for i in range(num_processors)] # Tells comm.send which mini-chunk to send (tells the index) to a certain processor
    chunks_left_to_send = [True for i in range(num_processors)] # List of boolean variables that says whether there are mini-chunks left to send to each processor
    while(sum(chunks_left_to_send) > 0):
        for j in range(0, len(all_fitting_factor_combinations)):
            if mini_chunk_to_send[j] < len(all_fitting_factor_combinations[j]):
                comm.send(all_fitting_factor_combinations[j][mini_chunk_to_send[j]], dest=j)
                mini_chunk_to_send[j] += 1
            else:
                chunks_left_to_send[j] = False

if rank >= 0:
    random.seed()
    num_mini_chunks_to_recv = comm.recv(source=0) # Receive the number of mini-chunks it is going to receive
    processor_fitting_factor_combinations = []
#     fitting_factors_and_least_rmsd = [1e80, 0,0,0]
#     reactions = ["Reaction 1", "Reaction 2", "Reaction 3"] # Used in writing output string of fitting factors and the fake_performance_measure to output file

#     results = open("files_processor" + str(rank) + "/results_file",'w')
#     num_items_recv = 0
    for i in range(0, num_mini_chunks_to_recv): # Receive all the mini-chunks
        data = comm.recv(source=0)
#         num_items_recv += len(data)
        processor_fitting_factor_combinations.append(data)
#     num_samples_processed = 0
    for mini_chunk in processor_fitting_factor_combinations:
        pass
#         num_samples_processed += len(mini_chunk)
#         for combination in mini_chunk:
#             for i in range(100): # This is to simulate the work done by the processor (and the time it takes) when running monaco
#                 j = np.random.rand(256,256)
#             # For the parallel generalized baragiola optimization, will do work of inserting new parameters; running the run.sh file for this processor (includes monaco), and calculation of a performance measure  
#             fake_performance_measure = math.sqrt(sum([fitting_factor**2 for fitting_factor in combination]))
#             if (fake_performance_measure < fitting_factors_and_least_rmsd[0]):
#                 fitting_factors_and_least_rmsd[0] = fake_performance_measure
#                 fitting_factors_and_least_rmsd[1] = combination[0]
#                 fitting_factors_and_least_rmsd[2] = combination[1]
#                 fitting_factors_and_least_rmsd[3] = combination[2]
#             output_string  = ""
#             for i in range(0, len(reactions)): # Write to the output file for this processor (write reactions and associated fitting factors and the fake performance metric)
#                 output_string = output_string + str(combination[i]) + "".join(" "*(23 - len(str(combination[i])))) + reactions[i] + " delta values \n"
#             output_string += str(fake_performance_measure) + "".join(" "*(23 - len(str(fake_performance_measure)))) + "RMSD" + "\n\n"
#             results.write(output_string)
#         # print("The number of samples processed is " + str(num_samples_processed))
#     results.close()

# least_rmsds_and_fitting_factors = comm.gather(fitting_factors_and_least_rmsd, root=0) # Gather the least fake performance metric value and associated fitting factors from each processor back to the root processor (0)

# if rank == 0:
#     least_rmsd_index = 0
#     for i in range(1, len(least_rmsds_and_fitting_factors)): # Loop through the least fake performance metric value and associated fitting factors 
#                                                              # from each processor and find the ones with the least value for the fake performance metrix
#         if least_rmsds_and_fitting_factors[i][0] < least_rmsds_and_fitting_factors[least_rmsd_index][0]:
#             least_rmsd_index = i
#     # os.system("rm -r files_processor*")

#     print("The smallest performance metric value and associated fitting factors: ")
#     print(least_rmsds_and_fitting_factors[least_rmsd_index])
    
#     # Optional test to verify that the fake performance metric value is the square root of the sum of the squares of its associated fitting factors
#     # for array in least_rmsds_and_fitting_factors:
#     #     fake_performance_measure = 0
#     #     for i in range(1,4):
#     #         fake_performance_measure += array[i]**2
#     #     assert array[0] == math.sqrt(fake_performance_measure)