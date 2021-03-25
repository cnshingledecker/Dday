import numpy as np
import csv, itertools, math, os, time
from exportable_custom_functions import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

reactions = [] # Will hold the reactions for which the fitting factors are being modified
all_vector_args = [] # Holds a series of lists (with each list containing the values in each of the linspaces that is created for a reaction)
base_dir_name = "baragiola_files_processor" # The partial name for each directory of the files for a processor

# Note: the below code is ran on every processor because each processor needs the reaction, and reading it on each processor means the data from the file doesn't have to be sent to each processor
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
        if(i == 1):
            num_fitting_factor_combinations = 0
            with open("experimental_data/num_fitting_factor_combinations.csv", "r") as num_file:
                for line in num_file:
                    num_fitting_factor_combinations = int(line)
            args_for_single_linspace[2] = float(num_fitting_factor_combinations)
        else:
            args_for_single_linspace[2] = 1
        all_vector_args.append(args_for_single_linspace)

# Notes: Processor is the one that handles the generation and distribution of fitting factors and the collection of data.
#        In this if statement, a directory is created for each processor with the files necessary for it to run monaco as well as find the rmsd and write to output files,
#        the fitting factor combinations are generated and distributed to each processor (including itself (processor 0)),         
if rank == 0: 
    fitting_factors = [] # Holds lists (there will be 3; each list is a numpy linspace (converted to a list) that was created using the arguments read in from the csv input file above)
    num_processors = 4
    for vector_args in all_vector_args: # Create numpy linspace out using the parameters in vector_args (read from an input file)
        single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
        single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
        fitting_factors.append(single_vector_fitting_factors)  
    
    all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of fitting factor values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates a list of n lists (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n lists has k elements (where k is the number of fitting factors (num_modified_fitting_factors))
    all_fitting_factor_combinations = [list(fitting_factor_combination) for fitting_factor_combination in all_fitting_factor_combinations]
    all_fitting_factor_combinations = split_list(all_fitting_factor_combinations, num_processors) # Splits the list into 'num_processors' chunks
    all_fitting_factor_combinations = split_list_chunks(all_fitting_factor_combinations, 15) # Splits each of the list chunks into mini-chunks of up to size 15 (creates as many with size 15 as possible)
                                                                                              # Note: Sending a mini-chunk on my (Daniel's) machine does not arrive at the destination (the program just sits and the mini-chunk
                                                                                              #       never gets to its destination). Feel free to change this if desired and if your machine can handle a smaller or bigger mini-chunk size.
    # These file are copied into a directory for each processor
    files_to_copy_to_new_dir = ["clean.sh", "run.sh", "save_results.mod", "run_dvode.mod", "read_rate06.mod", "read_model.mod", 
                                "rd_eff.txt", "radiolysis.dat", "photo_processes.dat", "parameter_inputs_template.dat",
                                "network.dat", "monaco", "model.inp", "mod_save_results.f90", "mod_run_dvode.f90", 
                                "mod_read_rate06.f90", "mod_read_model.f90", "mod_global_variables.f90", "mod_global_functions.f90",
                                "mod_calculate_rates.f90", "Lee_ea_17.txt", "init_surf_ab.inp", "init_gas_ab.inp", "init_bulk_ab.inp",
                                "global_variables.mod", "global_functions.mod", "enthalpias.txt", "dvode_f90_m.mod",
                                "dvode_f90_m.f90", "class_2_suprathermal.dat", "chem_rate06_dvode.f90", "calculate_rates.mod"] 

    for i in range(0, num_processors): # Create directory for the files for each processor, copy into it the files specified in the above array, copy the experimental_data directory into it, and create the results file in each array
        new_dir_name = base_dir_name + str(i)
        os.system("mkdir " + new_dir_name)
        os.system("touch " + new_dir_name + "/output_file_results") # Create file that each processor will write to (will write fitting factors and associated rmsds)
        for file_name in files_to_copy_to_new_dir:
            os.system("cp " + file_name + " " + new_dir_name)
        os.system("cp -R experimental_data " + new_dir_name)
    
    for i in range(0, len(all_fitting_factor_combinations)): # Tell each processor how many mini-chunks it is going to be receiving
        comm.send(len(all_fitting_factor_combinations[i]), dest=i)
    
    # Send all of the mini_-chunks to the different processors
    mini_chunk_to_send = [0 for i in range(num_processors)] # Tells comm.send which mini-chunk to send (tells the index) to a certain processor
    chunks_left_to_send = [True for i in range(num_processors)] # List of boolean variables that says whether there are mini-chunks left to send to each processor
    while(sum(chunks_left_to_send) > 0):
        for j in range(0, len(all_fitting_factor_combinations)):
            if mini_chunk_to_send[j] < len(all_fitting_factor_combinations[j]): # If the possible index (of a mini-chunk to send) is a valid one (less than the number of mini-chunks)
                comm.send(all_fitting_factor_combinations[j][mini_chunk_to_send[j]], dest=j)
                mini_chunk_to_send[j] += 1 # Increment the index for the next run of the loop
            else:
                chunks_left_to_send[j] = False

if rank >= 0:
    new_dir_name = base_dir_name + str(rank) # The name of the directory conaining the files for this processor
    num_mini_chunks_to_recv = comm.recv(source=0) # Receive the number of mini-chunks it is going to receive
    processor_fitting_factor_combinations = []
    fitting_factors_and_least_rmsd = [1e80, 0,0,0] # Holds the least rmsd and the fitting factors that produced it; fitting_factors_and_least_rmsd[0] is the rmsd, indices 1-3 are the fitting factors that led to that rmsd value

    results = open(new_dir_name + "/output_file_results",'w')
    for i in range(0, num_mini_chunks_to_recv): # Receive all the mini-chunks
        data = comm.recv(source=0)
        processor_fitting_factor_combinations.append(data)
    for mini_chunk in processor_fitting_factor_combinations:
        for fitting_factor_combination in mini_chunk:
            infile = open(new_dir_name + "/parameter_inputs_template.dat",'r')
            outfile = open(new_dir_name + "/photo_processes_2.dat",'w')
            for line in infile:
                line_as_list = line.split()  # Convert the DAT file line into a list
                # We need to find the first index in the line (as a list) where there is a fitting factor value (it is the first number besides the first element of the list, hence why i (below) is set to 1)
                index_first_float_value_in_line = 1  # The first entry in the line is a float, but we are looking for the first fitting factor (the 2nd float value in the line); starting at 1 therefore makes the implementation of the below while loop easier
                # The while loop below skips the first value because it is an integer and would thus be counted as a float. 
                while(index_first_float_value_in_line < len(line_as_list) and not(is_float(line_as_list[index_first_float_value_in_line]))):  # Increments the index value until the first float value (fitting factor value) is found
                    index_first_float_value_in_line += 1
                possible_fitting_factor_index = line_as_list[len(line_as_list) - 1] # See below comment for the meaning of this variable
                if(is_int(possible_fitting_factor_index)): # If the last value of the line is an integer, which means it is one of the reactions for which we are modifying the fitting factors
                    new_fitting_factor_val = fitting_factor_combination[int(possible_fitting_factor_index)] # Convert the possible fitting factor index to an int and get the correct fitting factor (where fitting_factor_index tells the computer which fitting factor to get)
                    new_fitting_factor_val = np.format_float_scientific(new_fitting_factor_val, precision=2,unique=False)  # Convert the new fitting factor value into a string of a number in scientific notation rounded to 2 places after the decimal point
                    line = line[0:106] + new_fitting_factor_val + line[114:len(line)] # In the line, replace the old fitting factor with the new value             
                outfile.write(line)
            infile.close()
            outfile.close()
            print("Running model...processor " + str(rank))
            os.system('cd ' + new_dir_name + '; ./run.sh') # Includes a run of monaco (note that what these commands do is temporarily dipping down into the directory of the files for this processor and running run.sh); after runnign these commands, the current directory is the same as it was before these commands were run
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

least_rmsds_and_fitting_factors = comm.gather(fitting_factors_and_least_rmsd, root=0) # Gather the least fake performance metric value and associated fitting factors from each processor back to the root processor (0)

if rank == 0:
    least_rmsd_index = 0
    for i in range(1, len(least_rmsds_and_fitting_factors)): # Loop through the least fake performance metric value and associated fitting factors 
                                                             # from each processor and find the ones with the least value for the fake performance metrix
        if least_rmsds_and_fitting_factors[i][0] < least_rmsds_and_fitting_factors[least_rmsd_index][0]:
            least_rmsd_index = i
    # print("The smallest performance metric value and fitting factors that produced it: ")
    # print(least_rmsds_and_fitting_factors[least_rmsd_index])
    results_file = open("results_generalized_parallel", 'w')
    output_string = ""
    for i in range(1, len(reactions)): 
        output_string = output_string + str(least_rmsds_and_fitting_factors[least_rmsd_index][i]) + "".join(" "*(23 - len(str(least_rmsds_and_fitting_factors[least_rmsd_index][i])))) + reactions[i] + " delta values \n"
    output_string += str(least_rmsds_and_fitting_factors[least_rmsd_index][0]) + "".join(" "*(23 - len(str(least_rmsds_and_fitting_factors[least_rmsd_index][0])))) + "RMSD" + "\n\n"
    results_file.write(output_string)
    results_file.close()