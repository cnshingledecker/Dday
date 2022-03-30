#! /usr/bin/python

import os
import csv
import numpy as np
import itertools
import time
from exportable_custom_functions import find_nearest_index,is_float, is_int, modify_modelInp_values

startTime = time.time()

results = open("resultsFile_2", 'w')

reset_model_inp = False # If this is true, model.inp will be reset to default values 
                        #     (specified in modelCopy.inp,; model.inp will be overwritten with the contents of this file)
if(reset_model_inp == True):
        os.system("cat modelCopy.inp > model.inp") # Reset model.inp to the default values in modelCopy.inp only if the user wants to (if reset_model_inp is true)

to_modify_modelInp_values = False # Set this to True if you want to modify model.inp values using the below array 
lines_to_modify_modelInp = [] 
                                # If the user wants to modify model.inp values, the user should insert lists of the form [lineNumber, variableVal, variableName] into this list,
                                #     where lineNumber is an integer that is the 0-indexed line number of the line in model.inp of the variable the user wants to change,
                                #     variableVal is of a real number data type (the value that the user wants to set the variable in lineNumber to), and
                                #     and variableName is a string that is the exact case of the variable name in the line in model.inp
if(to_modify_modelInp_values):
    modify_modelInp_values(lines_to_modify_modelInp)

reactions = [] # Will hold the reactions for which the fitting factors are being modified
all_vector_args = [] # Holds a series of arrays (one for each of the linspaces that is created for a reaction)
num_modified_fitting_factors = 0
fitting_factors = [] # Holds arrays containing the fitting_factors for the reactions (each element of the array is a list with the various fitting factors for a reaction)

# Note: the below code is ran on every core because each core needs the reaction, and reading it on each core means the data from the file doesn't have to be sent to each core
with open('reaction_fitting_factor_linspace_args/reaction_fitting_factor_vector_arguments.csv', newline='') as vector_creation_args_csv:  # Read in the parameters from the csv file for the creation of the linspaces (for each fitting factor to be varied)
    reader = csv.reader(vector_creation_args_csv, delimiter=',')      
    for i in range(0, 4): # Skips the first 4 lines of the csv file (lines which are comments)
        fields = next(reader)
    
    i = 0
    for row in reader:  # Each row is a set of arguments to be used to create the linspace for fitting factors for a reaction
        i += 1
        num_modified_fitting_factors += 1
        args_for_single_linspace = []  # A vector to hold the set of arguments to be used to create the linspace for the fitting factors for a reaction
        for argument in row:
            if is_float(argument): # It is one of the arguments to be used to make a numpy linspace
                args_for_single_linspace.append(float(argument)) # Adds the arguments to the previously created vector (3 lines above)
            else:  # It is text (it is not an index specifying a delta value to choose; because it is text, it must be the reaction(s) for which the fitting factors are being varied using the linspace created using some of the values in this row)
                reactions.append(argument) # Note that for each fitting factor combination; fitting_factor_combination[i] is the fitting factor associated with reaction[i]
        all_vector_args.append(args_for_single_linspace)

# Holds the least rmsd and the fitting factors that produced it; fitting_factors_and_least_rmsd[0] is the rmsd, 
#    indices 1 to len(fitting_factors_and_least_rmsd) (inclusive) are the fitting factors that led to that rmsd value
fitting_factors_and_least_rmsd = [0] * (len(reactions) + 1) # there needs to be 1 spot for the rmsd and len(reactions) spots for the fitting factor for each set of reactions 
fitting_factors_and_least_rmsd[0] = 1e80 # Initialize the RMSD to a high value so during testing a lower RMSD will likely be found and be saved in the array (along with the fitting factors that produced it)

for vector_args in all_vector_args:
    single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
    single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
    fitting_factors.append(single_vector_fitting_factors) 

all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of fitting factor values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates a list of n lists (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n lists has k elements (where k is the number of fitting factors (num_modified_fitting_factors))

for fitting_factor_combination in all_fitting_factor_combinations:  # fitting_factor_combinations[n] is a value from the (n + 1)th np.linspace created 
    infile = open("parameter_inputs_template.dat",'r')
    outfile = open("photo_processes.dat",'w')
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
    print("Running model...")
    os.system('./runGeneralization.sh') # Executes a version of run.sh that runs the monaco executable for this generalization (a monaco that uses certain fortran files for this generalization)
    print("Finding RMSD...")  # RMSD is root-mean square deviation

    num_experimental_data_points = 0
    with open('experimental_data/experimental_o3.csv') as csv_file: # Experimental data
        deviations = [] # The deviations for each model value from the experimental data value (at the closest time)
        csv_reader = csv.reader(csv_file, delimiter=',')
        csv_model_data = open('csv/total_ice_O3.csv')
        csv_model_data_reader = csv.reader(csv_model_data, delimiter=',')
        csv_model_data_list = list(csv_model_data_reader)
        csv_model_data.close()
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
            for i in range(1, len(fitting_factors_and_least_rmsd)):
                fitting_factors_and_least_rmsd[i] = fitting_factor_combination[i-1]
        csv_file.close()
results.close()

print("The smallest performance metric value and fitting factors that produced it: ")
print(fitting_factors_and_least_rmsd)
results_file = open("resultsRMSDGeneralization", 'w')
output_string = ""
for i in range(0, len(reactions)): 
    output_string = output_string + str(fitting_factors_and_least_rmsd[i+1]) + "".join(" "*(23 - len(str(fitting_factors_and_least_rmsd[i+1])))) + reactions[i] + " delta values \n"
output_string += str(fitting_factors_and_least_rmsd[0]) + "".join(" "*(23 - len(str(fitting_factors_and_least_rmsd[0])))) + "RMSD" + "\n\n"
results_file.write(output_string)
results_file.close()

endTime = time.time()
timeTaken = endTime - startTime
print("Time taken: " + str(timeTaken / 60) + " minutes.")

# Put in best fitting factor combination found (before generating a plot for it)

infile = open("parameter_inputs_template.dat",'r')
outfile = open("photo_processes.dat",'w')
for line in infile:
    line_as_list = line.split()  # Convert the DAT file line into a list
    possible_fitting_factor_index = line_as_list[len(line_as_list) - 1] # See below comment for the meaning of this variable
    if(is_int(possible_fitting_factor_index)): # If the last value of the line is an integer, which means it is one of the reactions for which we are modifying the fitting factors
        if(float(line_as_list[len(line_as_list) - 3]) > 0):
            new_fitting_factor_val = fitting_factors_and_least_rmsd[int(possible_fitting_factor_index) + 1] # Convert the possible fitting factor index to an int and get the correct fitting factor (where fitting_factor_index tells the computer which fitting factor to get)
            new_fitting_factor_val = np.format_float_scientific(new_fitting_factor_val, precision=2,unique=False)  # Convert the new fitting factor value into a string of a number in scientific notation rounded to 2 places after the decimal point
            line = line[0:106] + new_fitting_factor_val + line[114:len(line)] # In the line, replace the old fitting factor with the new value             
    outfile.write(line)
infile.close()
outfile.close()

print("Running model with best fit parameters...")
os.system('./runGeneralization.sh')

# Do stuff to create plot
os.system("python3 dataFrameCreation.py")
os.system("python3 plotting.py")

print("Done!")