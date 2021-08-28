import numpy as np
import csv, itertools, math, os, time
from exportable_custom_functions import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

reactions = [] # Will hold the reactions for which the fitting factors are being modified
all_vector_args = [] # Holds a series of lists (with each list containing the values in each of the linspaces that is created for a reaction)
base_dir_name = "baragiola_files_core" # The partial name for each directory of the files for a core

# Note: the below code is ran on every core because each core needs the reaction, and reading it on each core means the data from the file doesn't have to be sent to each core
with open('reaction_fitting_factor_linspace_args/reaction_fitting_factor_vector_arguments.csv', newline='') as vector_creation_args_csv:  # Read in the parameters from the csv file for the creation of the linspaces (for each fitting factor to be varied)
    reader = csv.reader(vector_creation_args_csv, delimiter=',')      
    for i in range(0, 4): # Skips the first 4 lines of the csv file (lines which are comments)
        fields = next(reader)
    
    i = 0
    for row in reader:  # Each row is a set of arguments to be used to create the linspace for fitting factors for a reaction
        i += 1
        args_for_single_linspace = []  # A vector to hold the set of arguments to be used to create the linspace for the fitting factors for a reaction
        for argument in row:
            if is_float(argument): # It is one of the arguments to be used to make a numpy linspace
                args_for_single_linspace.append(float(argument)) # Adds the arguments to the previously created vector (3 lines above)
            else:  # It is text (it is not an index specifying a delta value to choose; because it is text, it must be the reaction(s) for which the fitting factors are being varied using the linspace created using some of the values in this row)
                reactions.append(argument) # Note that for each fitting factor combination; fitting_factor_combination[i] is the fitting factor associated with reaction[i]
        all_vector_args.append(args_for_single_linspace)

# Notes: The core is the one that handles the generation and distribution of fitting factors and the collection of data.
#        In this if statement, a directory is created for each core with the files necessary for it to run monaco as well as find the rmsd and write to output files,
#        the fitting factor combinations are generated and distributed to each core (including itself (core 0)),         
if rank == 0: 
    fitting_factors = [] # Holds lists (there will be 3; each list is a numpy linspace (converted to a list) that was created using the arguments read in from the csv input file above)
    num_cores = 4
    for vector_args in all_vector_args: # Create numpy linspace out using the parameters in vector_args (read from an input file)
        single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
        single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
        fitting_factors.append(single_vector_fitting_factors)  
    
    all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of fitting factor values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates a list of n lists (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n lists has k elements (where k is the number of fitting factors (num_modified_fitting_factors))
    all_fitting_factor_combinations = [list(fitting_factor_combination) for fitting_factor_combination in all_fitting_factor_combinations]
    all_fitting_factor_combinations = split_list(all_fitting_factor_combinations, num_cores) # Splits the list into 'num_cores' chunks
    all_fitting_factor_combinations = split_list_chunks(all_fitting_factor_combinations, 15) # Splits each of the list chunks into mini-chunks of up to size 15 (creates as many with size 15 as possible)
                                                                                              # Note: Sending a mini-chunk on my (Daniel's) machine does not arrive at the destination (the program just sits and the mini-chunk
                                                                                              #       never gets to its destination). Feel free to change this if desired and if your machine can handle a smaller or bigger mini-chunk size.
    # These file are copied into a directory for each core
    files_to_copy_to_new_dir = ["clean.sh", "run.sh", "rd_eff.txt", "radiolysis.dat", "photo_processes.dat", "parameter_inputs_template.dat",
                                "network.dat", "monaco", "model.inp", "Lee_ea_17.txt", "init_surf_ab.inp", "init_gas_ab.inp", 
                                "init_bulk_ab.inp", "enthalpias.txt", "class_2_suprathermal.dat"] 

    for i in range(0, num_cores): # Create directory for the files for each core, copy into it the files specified in the above array, copy the experimental_data directory into it, and create the results file in each array
        new_dir_name = base_dir_name + str(i)
        os.system("mkdir " + new_dir_name)
        os.system("touch " + new_dir_name + "/output_file_results") # Create file that each core will write to (will write fitting factors and associated rmsds)
        for file_name in files_to_copy_to_new_dir:
            os.system("cp " + file_name + " " + new_dir_name)
        os.system("cp -R experimental_data " + new_dir_name)
    
    for i in range(0, len(all_fitting_factor_combinations)): # Tell each core how many mini-chunks it is going to be receiving
        comm.send(len(all_fitting_factor_combinations[i]), dest=i)
    
    # Send all of the mini_-chunks to the different cores
    mini_chunk_to_send = [0 for i in range(num_cores)] # Tells comm.send which mini-chunk to send (tells the index) to a certain core
    chunks_left_to_send = [True for i in range(num_cores)] # List of boolean variables that says whether there are mini-chunks left to send to each core
    while(sum(chunks_left_to_send) > 0):
        for j in range(0, len(all_fitting_factor_combinations)):
            if mini_chunk_to_send[j] < len(all_fitting_factor_combinations[j]): # If the possible index (of a mini-chunk to send) is a valid one (less than the number of mini-chunks)
                comm.send(all_fitting_factor_combinations[j][mini_chunk_to_send[j]], dest=j)
                mini_chunk_to_send[j] += 1 # Increment the index for the next run of the loop
            else:
                chunks_left_to_send[j] = False

if rank >= 0:
    new_dir_name = base_dir_name + str(rank) # The name of the directory conaining the files for this core
    num_mini_chunks_to_recv = comm.recv(source=0) # Receive the number of mini-chunks it is going to receive
    core_fitting_factor_combinations = []
    fitting_factors_and_least_rmsd = [1e80, 0,0,0] # Holds the least rmsd and the fitting factors that produced it; fitting_factors_and_least_rmsd[0] is the rmsd, indices 1-3 are the fitting factors that led to that rmsd value

    results = open(new_dir_name + "/output_file_results",'w')
    for i in range(0, num_mini_chunks_to_recv): # Receive all the mini-chunks
        data = comm.recv(source=0)
        core_fitting_factor_combinations.append(data)
    for mini_chunk in core_fitting_factor_combinations:
        for fitting_factor_combination in mini_chunk:
            infile = open(new_dir_name + "/parameter_inputs_template.dat",'r')
            outfile = open(new_dir_name + "/photo_processes_2.dat",'w')
            for line in infile:
                line_as_list = line.split()  # Convert the DAT file line into a list
                possible_fitting_factor_index = line_as_list[len(line_as_list) - 1] # See below comment for the meaning of this variable
                if(is_int(possible_fitting_factor_index)): # If the last value of the line is an integer, which means it is one of the reactions for which we are modifying the fitting factors
                    if(float(line_as_list[len(line_as_list) - 3]) > 0):
                        new_fitting_factor_val = fitting_factor_combination[int(possible_fitting_factor_index)] # Convert the possible fitting factor index to an int and get the correct fitting factor (where fitting_factor_index tells the computer which fitting factor to get)
                        new_fitting_factor_val = np.format_float_scientific(new_fitting_factor_val, precision=2,unique=False)  # Convert the new fitting factor value into a string of a number in scientific notation rounded to 2 places after the decimal point
                        line = line[0:106] + new_fitting_factor_val + line[114:len(line)] # In the line, replace the old fitting factor with the new value             
                outfile.write(line)
            infile.close()
            outfile.close()
            print("Running model...core " + str(rank))
            os.system('cd ' + new_dir_name + '; ./run.sh') # Includes a run of monaco (note that what these commands do is temporarily dipping down into the directory of the files for this core and running run.sh); after runnign these commands, the current directory is the same as it was before these commands were run
            print("Finding RMSD...")  # RMSD is root-mean square deviation

            num_experimental_data_points = 0
            with open(new_dir_name + '/experimental_data/experimental_o3.csv') as csv_file: # Experimental data
                deviations = [] # The deviations for each model value from the experimental data value (at the closest time)
                csv_reader = csv.reader(csv_file, delimiter=',')
                csv_model_data = open(new_dir_name + '/csv/total_ice_o3.csv')
                csv_model_data_reader = csv.reader(csv_model_data, delimiter=',')
                csv_model_data_list = list(csv_model_data_reader)
                for row in csv_reader: # This assumes that row[0] is the time, row[1] is the y-value (I'm not sure what this is),
                    closest_model_values = csv_model_data_list[find_nearest_index(row[0], 0, csv_model_data_list)]
                    deviation = float(closest_model_values[1])*1e7 - float(row[1]) # Deviation of the model value from the actual (experimental) value
                    
                    # the deviation of the model from the y-value is allowed to be up to 10% away from the y-value
                    if 0.9 * float(row[1]) <= deviation <= 1.1 * float(row[1]): # 0.9 * float(row[1]) is the allowed_lower_deviation, 1.1 * float(row[1]) is the allowed_upper_deviation
                        deviation = 0
                    deviations.append(deviation)
                    num_experimental_data_points += 1
                sum = 0
                for value in deviations:
                    sum += (value**2)
                rmsd = (sum / (num_experimental_data_points - 2))**0.5   # Formula for RMSD

                # Should I create a boolean to only write to the output string if there was data in the experimental data csv file?
                
                # Create a string to hold the rmsd along with the fitting factor value for each reaction set (the fitting factor values combination)
                output_string = ""
                for i in range(0, len(reactions)): 
                    output_string = output_string + str(fitting_factor_combination[i]) + "".join(" "*(23 - len(str(fitting_factor_combination[i])))) + reactions[i] + " delta values \n"
                output_string += str(rmsd) + "".join(" "*(23 - len(str(rmsd)))) + "RMSD" + "\n\n"
                results.write(output_string)

                if (rmsd < fitting_factors_and_least_rmsd[0]): # fitting_factors_and_least_rmsd[0] is the least rmsd; if the new rmsd is less than it, store the new rmsd and the fitting factors that produced it
                    fitting_factors_and_least_rmsd[0] = rmsd
                    fitting_factors_and_least_rmsd[1] = fitting_factor_combination[0]
                    fitting_factors_and_least_rmsd[2] = fitting_factor_combination[1]
                    fitting_factors_and_least_rmsd[3] = fitting_factor_combination[2]
                csv_file.close()
    results.close()

least_rmsds_and_fitting_factors = comm.gather(fitting_factors_and_least_rmsd, root=0) # Gather the least fake performance metric value and associated fitting factors from each core back to the root core (0)

if rank == 0:
    least_rmsd_index = 0
    for i in range(1, len(least_rmsds_and_fitting_factors)): # Loop through the least fake performance metric value and associated fitting factors 
                                                             # from each core and find the ones with the least value for the fake performance metrix
        if least_rmsds_and_fitting_factors[i][0] < least_rmsds_and_fitting_factors[least_rmsd_index][0]:
            least_rmsd_index = i
    # print("The smallest performance metric value and fitting factors that produced it: ")
    # print(least_rmsds_and_fitting_factors[least_rmsd_index])
    results_file = open("results_generalized_parallel", 'w')
    output_string = ""
    # print(least_rmsds_and_fitting_factors[least_rmsd_index])
    for i in range(1, len(reactions) + 1): 
        output_string = output_string + str(least_rmsds_and_fitting_factors[least_rmsd_index][i]) + "".join(" "*(23 - len(str(least_rmsds_and_fitting_factors[least_rmsd_index][i])))) + reactions[i-1] + " delta values \n"
    output_string += str(least_rmsds_and_fitting_factors[least_rmsd_index][0]) + "".join(" "*(23 - len(str(least_rmsds_and_fitting_factors[least_rmsd_index][0])))) + "RMSD" + "\n\n"
    results_file.write(output_string)
    results_file.close()