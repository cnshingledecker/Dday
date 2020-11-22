#! /usr/bin/python

import os
import csv
import numpy as np
import itertools

def isInt(val):
    try:
        int(val)
        return True
    except ValueError:
        return False

def isFloat(val):
    try:
        float(val)
        return True
    except ValueError:
        return False

def find_nearest_index(value_to_find, pos, in_order_list):
    index = 0
    smallest_index = 0
    smallest_diff = 1e20
    for value in in_order_list:
        diff = abs(value_to_find - value[pos])
        if diff < smallest_diff:
            smallest_index = index
            smallest_diff = diff
        elif diff > smallest_diff:
            return smallest_index
        index += 1
    return smallest_index

results = open("results_2", 'w')

reactions = []

all_vector_args = [] # Holds a series of arrays (one for each of the linspaces that is cretaed for a reaction)

with open('reaction_fitting_factor_linspace_args/reaction_fitting_factor_vector_arguments.csv', newline='') as vector_creation_args_csv:  # Read in the parameters from the csv file for the creation of the linspaces (for each fitting factor to be varied)
    reader = csv.reader(vector_creation_args_csv, delimiter=',')      
    fields = next(reader) # Skips first line of the csv file (the line contains headers and info about what each line means)
    fields = next(reader) # Skips second line of the csv file (the line describes how each line in the csv file corresponds to the values in each fitting_factor_combination in all_fitting_factor_combinations)
    fields = next(reader) # Skips the third line of the csv file (a contnuation of the second line)
    for row in reader:  # Each row is a set of arguments to be used to create the linspace for fitting factors for a reaction
        args_for_single_vector = []  # A vector to hold the set of arguments to be used to create the linspace for the fitting factors for a reaction
        for argument in row:
            if isFloat(argument):
                args_for_single_vector.append(float(argument)) # Adds the arguments to the previously created vector (2 lines above)
            else:  # It is text (it must be the reaction(s) for which the fitting factors are being varied using the linspace created using some of the values in this row)
                reactions.append(argument)
        all_vector_args.append(args_for_single_vector) 

fitting_factors = [] # Holds arrays containing the fitting_factors for the reactions (each element of the array is a list with the various fitting factors for a reaction)

for vector_args in all_vector_args:
    single_vector_fitting_factors = np.linspace(vector_args[0], vector_args[1], num=int(vector_args[2])) # Creates numpy linspace of the fitting factors (for the reaction) using arguments peeviously retrieved from the csv file
    single_vector_fitting_factors = list(single_vector_fitting_factors) # Converts the numpy linspace to a list
    fitting_factors.append(single_vector_fitting_factors) 

all_fitting_factor_combinations = itertools.product(*fitting_factors)    # Creates all possible combinations of delta values from the previously created lists (created in the 'for vector args in all_vector_args' for loop)
                                                                         #    Note that this creates an array of n arrays (where n is the product of the number of values for each list of fitting factors), and each of these 
                                                                         #    n arrays has k elements (where k is the number of fitting factors)

for fitting_factor_combination in all_fitting_factor_combinations:  # fitting_factor_combinations[n] is a value from the (n + 1)th np.linspace created 
    infile = open("parameter_inputs_template.dat",'r')
    outfile = open("photo_processes.dat",'w')
    for line in infile:
        line_as_list = line.split()  # Convert the DAT file line into an array
        # Need to find the first index in the line (as a list) where there is a delta value (it is the first number besides the first element of the list, hence why i (below) is set to 1)
        index_first_float_value_in_line = 1  # Skips the first value because it is an integer and would thus be counted as a float. 
        while(index_first_float_value_in_line < len(line_as_list) and not(isFloat(line_as_list[index_first_float_value_in_line]))):
            index_first_float_value_in_line += 1
        possible_fitting_factor_index = line_as_list[len(line_as_list) - 1] # See below comment for the meaning of this variable
        if(isInt(possible_fitting_factor_index)): # If the last value of the line is an integer, which means it is one of the reactions for which we are modifying the fitting factors
            fitting_factor_index = int(possible_fitting_factor_index)            
            new_fitting_factor_val = fitting_factor_combination[fitting_factor_index] # Get the correct fitting factor (where fitting_factor_index tells the computer which fitting factor to get)
            if(new_fitting_factor_val < 1.00):
                # Note that in baragiola_optimization.py, new_fitting_factor_value was multiplied by 10 here; in this code, I commented out the line to produce the same results as baragiola_optimization.py--I'm not sure why
                new_fitting_factor_val = round(new_fitting_factor_val, 2)
            new_fitting_factor_val = np.format_float_scientific(new_fitting_factor_val, precision=2,unique=False)  # Convert the new fitting factor value into a string of a number in scientific notation rounded to 2 places after the decimal point
            line = line[0:106] + new_fitting_factor_val + line[114:len(line)]
        outfile.write(line)
    infile.close()
    outfile.close()
    print("Running model...")
    os.system('./run.sh')
    print("Finding RMSD...")

    with open('experimental_data_file.csv') as csv_file:
        deviations = []
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader: # This assumes that row[0] is the time, row[1] is the y-value (I'm not sure what this is), 
                            #     if the length of the row is 4, row[2] is the specified lower bound that the deviation of the model from the y-value 
                            #     is allowed to be and row[3] is the specified upper bound for the deviation. If the row length is 3, row[2] is 
                            #     assumed to be a tolerance percentage (where the deviation (previously described) is allowed to be less than 
                            #     that percentage of the y-value away from the y-value)
            csv_model_data = open('name_of_model_data_file.csv')
            csv_model_data_reader = csv.reader(csv_model_data, delimiter=',')
            csv_model_data_list = list(csv_model_data_reader)
            closest_model_values = csv_model_data_list[find_nearest_index(row[0], 0, csv_model_data_list)]
            deviation = closest_model_values[1] - row[1]
            # Use a default tolerance value of 10% away from the y-value (the deviation of the model from the y-value is allowed to be that much away from the y-value)
            allowed_lower_deviation = 0.9 * row[1]
            allowed_upper_deviation = 1.1 * row[1]
            if allowed_lower_deviation < deviation < allowed_upper_deviation:
                deviation = 0
            deviations.append(deviation)
        sum = 0
        for value in deviations:
            sum += (value**2)
        rmsd = (sum / 16)**0.5
        output_string = "" # Create a string to hold the rmsd along with the delta value for each reaction set (the delta values combination)
        i = 0
        for reaction in reactions:  
            output_string = output_string + " \n" + reaction + " delta values: " + str(fitting_factor_combination[i])
            i += 1
        output_string += "\nRMSD: " + str(rmsd) + "\n\n"
        results.write(output_string)
results.close()
print("Done!")