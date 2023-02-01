import numpy as np
import csv, itertools, math, os, time
from exportable_custom_functions import *
from mpi4py import MPI

startTime = time.time()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

reactions = [] # Will hold the reactions for which the fitting factors are being modified
all_vector_args = [] # Holds a series of lists (with each list containing the values in each of the linspaces that is created for a reaction)
base_dir_name = "baragiola_files_core" # The partial name for each directory of the files for a core

to_modify_modelInp_values = True # Set this to True if you want to modify model.inp values using the below array 
reset_modelInp = True  # If this is true, model.inp will be reset to default values 
                        #     (specified in modelCopy.inp,; model.inp will be overwritten with the contents of this file)

experimental_data = setup_experimental_data() # The experimental data we compare the model to
initialO2 = 5.7E22
num_delta_values = 3

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
            argument = argument.strip() # Not working for some reason
            if is_float(argument): # It is one of the arguments to be used to make a numpy linspace
                args_for_single_linspace.append(float(argument)) # Adds the arguments to the previously created vector (3 lines above)
            else:  # It is text (it is not an index specifying a delta value to choose; because it is text, it must be the reaction(s) for which the fitting factors are being varied using the linspace created using some of the values in this row)
                reactions.append(argument) # Note that for each fitting factor combination; fitting_factor_combination[i] is the fitting factor associated with reaction[i]
        all_vector_args.append(args_for_single_linspace)
 
modified_lines_to_modify_modelInp = [] # Needed by the core with rank 0 but we create it here so we only process the file once
lines_to_modify_modelInp = get_data_to_modify_modelInp()
# If the user wants to modify model.inp values, the user should insert lists of the form [lineNumber, variableVal, variableName] into this list,
#     where lineNumber is an integer that is the 0-indexed line number of the line in model.inp of the variable the user wants to change,
#     variableVal is of a real number data type (the value that the user wants to set the variable in lineNumber to), and
#     and variableName is a string that is the exact case of the variable name in the line in model.inp

if(to_modify_modelInp_values):
    for line in lines_to_modify_modelInp:
        line_linspace = np.linspace(line[0], line[1], int(line[2])) # Range of values we want to try for that line's value
        line_linspace_list = list(line_linspace)
        line_linspace_list_power10 = []
        for value in line_linspace_list:
            line_linspace_list_power10.append(10**value)
        modified_lines_to_modify_modelInp.append(line_linspace_list_power10)

# Notes: The core is the one that handles the generation and distribution of fitting factors and the collection of data.
#        In this if statement, a directory is created for each core with the files necessary for it to run monaco as well as find the rmsd and write to output files,
#        the fitting factor combinations are generated and distributed to each core (including itself (core 0)),         
if rank == 0:
    if(reset_modelInp == True):
        os.system("cat modelCopy.inp > model.inp") # Reset model.inp to the default values in modelCopy.inp only if the user wants to (if reset_modelInp is true)
    
    fitting_factors = [] # Holds lists (each list is a numpy linspace (converted to a list) that was created using the arguments read in from the csv input file above) 
                         #     and lists for modifying model.inp values  
    num_processors = 4 # IMPORTANT: Need to adjust if running on a different number of processors
    for vector_args in all_vector_args: # Create numpy linspace out using the parameters in vector_args (read from an input file)
        single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
        single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
        single_vector_fitting_factors_power10 = []
        for value in single_vector_fitting_factors:
            single_vector_fitting_factors_power10.append(10**value)
        fitting_factors.append(single_vector_fitting_factors_power10)

    if(to_modify_modelInp_values == True):
        for list_of_values in modified_lines_to_modify_modelInp:
            fitting_factors.append(list_of_values)
    
    all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of fitting factor values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates a list of n lists (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n lists has k elements (where k is the number of fitting factors (num_modified_fitting_factors) plus the number of model.inp values we want to modify)

    all_fitting_factor_combinations = [list(fitting_factor_combination) for fitting_factor_combination in all_fitting_factor_combinations]
    all_fitting_factor_combinations = split_list(all_fitting_factor_combinations, num_processors) # Splits the list into 'num_processors' chunks
    all_fitting_factor_combinations = split_list_chunks(all_fitting_factor_combinations, 15) # Splits each of the list chunks into mini-chunks of up to size 15 (creates as many with size 15 as possible)
                                                                                              # Note: Sending a mini-chunk on my (Daniel's) machine does not arrive at the destination (the program just sits and the mini-chunk
                                                                                              #       never gets to its destination). Feel free to change this if desired and if your machine can handle a smaller or bigger mini-chunk size.
    # These file are copied into a directory for each core
    files_to_copy_to_new_dir = ["clean.sh", "run.sh", "rd_eff.txt", "radiolysis.dat", "photo_processes.dat", "parameter_inputs_template.dat",
                                "network.dat", "monaco", "model.inp", "Lee_ea_17.txt", "init_surf_ab.inp", "init_gas_ab.inp", 
                                "init_bulk_ab.inp", "enthalpias.txt", "class_2_suprathermal.dat"] 

    for i in range(0, num_processors): # Create directory for the files for each core, copy into it the files specified in the above array, copy the experimental_data directory into it, and create the results file in each array
        new_dir_name = base_dir_name + str(i)
        os.system("mkdir " + new_dir_name)
        os.system("touch " + new_dir_name + "/output_file_results") # Create file that each core will write to (will write fitting factors and associated rmsds)
        for file_name in files_to_copy_to_new_dir:
            os.system("cp " + file_name + " " + new_dir_name)
        os.system("cp -R experimental_data " + new_dir_name)
    
    for i in range(0, len(all_fitting_factor_combinations)): # Tell each core how many mini-chunks it is going to be receiving
        comm.send(len(all_fitting_factor_combinations[i]), dest=i)
    
    # Send all of the mini_-chunks to the different cores
    mini_chunk_to_send = [0 for i in range(num_processors)] # Tells comm.send which mini-chunk to send (tells the index) to a certain core
    chunks_left_to_send = [True for i in range(num_processors)] # List of boolean variables that says whether there are mini-chunks left to send to each core
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

    # Holds the least rmsd and the fitting factors that produced it; fitting_factors_and_least_rmsd[0] is the rmsd, 
    #    indices 1 to len(fitting_factors_and_least_rmsd) (inclusive) are the fitting factors that led to that rmsd value
    fitting_factors_and_least_rmsd = [0] * (len(reactions) + 1 + len(modified_lines_to_modify_modelInp)) # there needs to be 1 spot for the rmsd and len(reactions) spots for the fitting factor for each set of reactions, plus a spot for each of the model.inp lines to modify
    fitting_factors_and_least_rmsd[0] = 1e80 # Initialize the RMSD to a high value so during testing a lower RMSD will likely be found and be saved in the array (along with the fitting factors that produced it)

    results = open(new_dir_name + "/output_file_results",'w')
    for i in range(0, num_mini_chunks_to_recv): # Receive all the mini-chunks
        data = comm.recv(source=0)
        core_fitting_factor_combinations.append(data)
    for mini_chunk in core_fitting_factor_combinations:
        for fitting_factor_combination in mini_chunk:
            if(to_modify_modelInp_values == True):
                lines_to_modify_modelInp_local = []
                fitting_factor_combination_modelInp_index = num_delta_values # Because the fitting factors come first
                for line in lines_to_modify_modelInp:
                    lines_to_modify_modelInp_local.append([line[3], fitting_factor_combination[fitting_factor_combination_modelInp_index], line[4]])
                    fitting_factor_combination_modelInp_index += 1
                modify_modelInp_values(lines_to_modify_modelInp_local, new_dir_name)

            infile = open(new_dir_name + "/parameter_inputs_template.dat",'r')
            outfile = open(new_dir_name + "/photo_processes.dat",'w')
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
            # Run model, deal with files, and silence output
            os.system('cd ' + new_dir_name + '; ./run.sh > /dev/null') # Includes a run of monaco (note that what these commands do is temporarily dipping down into the directory of the files for this core and running run.sh); after runnign these commands, the current directory is the same as it was before these commands were run

            print("Finding RMSD...")  # RMSD is root-mean square deviation

            # Calculate the RMSD, write it and the parameters that produced it to an output file, 
            #     compare it to the best one found so far, and if it's less, store it and the parameters that produced it
            try: # If the model finished running, the CSV file will exist
                num_experimental_data_points = 0
                deviations = [] # The deviations for each model value from the experimental data value (at the closest time)

                csv_model_data = open(new_dir_name + '/csv/bO3.csv')
                csv_model_data_reader = csv.reader(csv_model_data, delimiter=',')
                throwAway = next(csv_model_data_reader)
                csv_model_data_list = list(csv_model_data_reader)
                for i in range(0, len(experimental_data['expX'])): 
                    experimentalY = experimental_data["expY"][i]
                    closest_model_values = csv_model_data_list[find_nearest_index(experimental_data["expX"][i], 0, csv_model_data_list)]

                    modelY = (float(closest_model_values[1]) / initialO2) * 100

                    deviation = modelY - float(experimentalY) # Deviation of the model value from the actual (experimental) value

                    # Daniel Lopez-Sanders: Not sure why this was in here; it didn't make sense so I commented it out
                    # the deviation of the model from the y-value is allowed to be up to 10% away from the y-value
                    # if 0.9 * float(experimentalY) <= deviation <= 1.1 * float(experimentalY): # 0.9 * float(experimentalY) is the allowed_lower_deviation, 1.1 * float(experimentalY) is the allowed_upper_deviation
                    #     deviation = 0
                    deviations.append(deviation)
                    num_experimental_data_points += 1
                sum = 0

                for value in deviations:
                    sum += (value**2)

                # NOTE: RMSD's of parallel and serial scripts were off for one run 18.8ish vs 8ish). Not an immediate significant cause for concern.
                rmsd = (sum / (num_experimental_data_points - 2))**0.5   # Formula for RMSD

                # Create a string to hold the rmsd along with the fitting factor value for each reaction set (the fitting factor values combination)
                output_string = ""
                for i in range(0, len(reactions)): 
                    fitting_factor_combination_formatted = np.format_float_scientific(fitting_factor_combination[i], precision=20,unique=False)
                    output_string = output_string + str(fitting_factor_combination_formatted) + "".join(" "*(30 - len(str(fitting_factor_combination_formatted)))) + reactions[i] + " delta values \n"
                for i in range(0, len(modified_lines_to_modify_modelInp)):
                    model_Inp_value_formatted = np.format_float_scientific(fitting_factor_combination[i + num_delta_values], precision=20,unique=False)
                    output_string = output_string + str(model_Inp_value_formatted) + "".join(" "*(30 - len(str(model_Inp_value_formatted)))) + lines_to_modify_modelInp[i][4] + " model.inp value\n"
                rmsd_formatted = np.format_float_scientific(rmsd, precision=20,unique=False)
                output_string += str(rmsd_formatted) + "".join(" "*(30 - len(str(rmsd_formatted)))) + "RMSD" + "\n\n"
                results.write(output_string)

                if (rmsd < fitting_factors_and_least_rmsd[0]): # fitting_factors_and_least_rmsd[0] is the least rmsd; if the new rmsd is less than it, store the new rmsd and the fitting factors that produced it
                    fitting_factors_and_least_rmsd[0] = rmsd
                    for i in range(1, len(reactions) + 1):
                        fitting_factors_and_least_rmsd[i] = fitting_factor_combination[i-1]
                    for i in range(1, len(modified_lines_to_modify_modelInp) + 1):
                        fitting_factors_and_least_rmsd[i + num_delta_values] = fitting_factor_combination[i + num_delta_values - 1]
            except FileNotFoundError:
                pass # Do nothing because the CSV file does not exist
    results.close()
    print(fitting_factors_and_least_rmsd)

fitting_factors_and_least_rmsd = comm.gather(fitting_factors_and_least_rmsd, root=0) # Gather the least performance metric value and associated fitting factors from each core back to the root core (0)

if rank == 0:
    least_rmsd_index = 0
    for i in range(1, len(fitting_factors_and_least_rmsd)): # Loop through the least fake performance metric value and associated fitting factors 
                                                             # from each core and find the ones with the least value for the fake performance metrix
        if fitting_factors_and_least_rmsd[i][0] < fitting_factors_and_least_rmsd[least_rmsd_index][0]:
            least_rmsd_index = i
    print("The smallest performance metric value and fitting factors that produced it: ")
    print(fitting_factors_and_least_rmsd[least_rmsd_index])
    results_file = open("resultsGeneralizedParallel", 'w')

    # Create a string to hold the rmsd along with the fitting factor value for each reaction set (the fitting factor values combination)
    output_string = ""
    for i in range(0, len(reactions)): 
       fitting_factor_combination_formatted = np.format_float_scientific(fitting_factors_and_least_rmsd[least_rmsd_index][i+1], precision=20,unique=False)
       output_string = output_string + str(fitting_factor_combination_formatted) + "".join(" "*(30 - len(str(fitting_factor_combination_formatted)))) + reactions[i] + " delta values \n"
    for i in range(1, len(modified_lines_to_modify_modelInp) + 1): # Because the RMSD is at the beginning of the list, 
                                                                   #     so we have to start 1 after we would if it wasn't at the beginning of the list
       model_Inp_value_formatted = np.format_float_scientific(fitting_factors_and_least_rmsd[least_rmsd_index][i + num_delta_values], precision=20,unique=False)
       output_string = output_string + str(model_Inp_value_formatted) + "".join(" "*(30 - len(str(model_Inp_value_formatted)))) + lines_to_modify_modelInp[i - 1][4] + " model.inp value\n"
    rmsd_formatted = np.format_float_scientific(fitting_factors_and_least_rmsd[least_rmsd_index][0], precision=20,unique=False)
    output_string += str(rmsd_formatted) + "".join(" "*(30 - len(str(rmsd_formatted)))) + "RMSD" + "\n\n"
    results_file.write(output_string)
    results_file.close()
    endTime = time.time()
    timeTaken = endTime - startTime
    print("Time taken: " + str(timeTaken / 60) + " minutes.")

    # Put in best fitting factor combination found (before generating a plot for it)
    new_dir_name = "."
    fitting_factor_combination = [0]*(len(fitting_factors_and_least_rmsd[least_rmsd_index])-1)
    for i in range(1, len(fitting_factors_and_least_rmsd[least_rmsd_index])):
        fitting_factor_combination[i-1] = fitting_factors_and_least_rmsd[least_rmsd_index][i]

    print("fitting_factor_combination is " + str(fitting_factor_combination))

    infile = open(new_dir_name + "/parameter_inputs_template.dat",'r')
    outfile = open(new_dir_name + "/photo_processes.dat",'w')
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

    if(to_modify_modelInp_values == True):
        lines_to_modify_modelInp_local = []
        fitting_factor_combination_modelInp_index = num_delta_values # Because the fitting factors come first
        for line in lines_to_modify_modelInp:
            lines_to_modify_modelInp_local.append([line[3], fitting_factor_combination[fitting_factor_combination_modelInp_index], line[4]])
            fitting_factor_combination_modelInp_index += 1
        modify_modelInp_values(lines_to_modify_modelInp_local, ".")

    print("Running model with best fit parameters...")
    os.system('./run.sh > /dev/null') # Run model, deal with files, and silence output

    # Create the plot
    os.system("python3 plotting.py")
