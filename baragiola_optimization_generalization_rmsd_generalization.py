#! /usr/bin/python

import os
import csv
import numpy as np
import itertools
from exportable_custom_functions import find_nearest_index,is_float, is_int

results = open("resultsFile_2", 'w')

reactions = [] # Will hold the reactions for which the fitting factors are being modified
all_vector_args = [] # Holds a series of arrays (one for each of the linspaces that is created for a reaction)
num_modified_fitting_factors = 0
fitting_factors = [] # Holds arrays containing the fitting_factors for the reactions (each element of the array is a list with the various fitting factors for a reaction)

with open('reaction_fitting_factor_linspace_args/reaction_fitting_factor_vector_arguments.csv', newline='') as vector_creation_args_csv:  # Read in the parameters from the csv file for the creation of the linspaces (for each fitting factor to be varied)
    reader = csv.reader(vector_creation_args_csv, delimiter=',')      
    for i in range(0, 4): # Skips the first 4 lines of the csv file (lines which are comments)
        fields = next(reader)
    for row in reader:  # Each row is a set of arguments to be used to create the linspace for fitting factors for a reaction
        num_modified_fitting_factors += 1
        args_for_single_vector = []  # A vector to hold the set of arguments to be used to create the linspace for the fitting factors for a reaction
        for argument in row:
            if is_float(argument):
                args_for_single_vector.append(float(argument)) # Adds the arguments to the previously created vector (3 lines above)
            else:  # It is text (it is not an index specifying a delta value to choose; because it is text, it must be the reaction(s) for which the fitting factors are being varied using the linspace created using some of the values in this row)
                reactions.append(argument)
        all_vector_args.append(args_for_single_vector) 

for vector_args in all_vector_args:
    single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
    single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
    fitting_factors.append(single_vector_fitting_factors) 

all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of fitting factor values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates a list of n lists (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n lists has k elements (where k is the number of fitting factors (num_modified_fitting_factors))

for fitting_factor_combination in all_fitting_factor_combinations:  # fitting_factor_combinations[n] is a value from the (n + 1)th np.linspace created 
    infile = open("parameter_inputs_template.dat",'r')
    outfile = open("photo_processes_2.dat",'w')
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
    print("Running model...")
    os.system('./runGeneralization.sh') # Executes a version of run.sh that runs the monaco executable for this generalization (a monaco that uses certain fortran files for this generalization)
    print("Finding RMSD...")  # RMSD is root-mean square deviation

    num_experimental_data_points = 0
    with open('experimental_data/experimental_o3.csv') as csv_file: # Experimental data
        deviations = [] # The deviations for each model value from the experimental data value (at the closest time)
        csv_reader = csv.reader(csv_file, delimiter=',')
        csv_model_data = open('csv/total_ice_o3.csv')
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
        csv_file.close()
results.close()
print("Done!")