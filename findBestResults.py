import os
import numpy as np
from exportable_custom_functions import is_int, modify_modelInp_values, get_data_to_modify_modelInp

fileNames = ["baragiola_files_core0/output_file_results", "baragiola_files_core1/output_file_results",
             "baragiola_files_core2/output_file_results", "baragiola_files_core3/output_file_results"]

bestRMSD = 1e80
lineValues = ["", "", "", "", "", "", ""]

for fileName in fileNames:
    with open(fileName, 'r') as file:
        lines = file.readlines()
        linesReversed = []
        lines = reversed(lines)
        for value in lines:
            linesReversed.append(value)
        for i in range(len(linesReversed)):
            currentLine = linesReversed[i]
            if currentLine.find("RMSD") != -1:
                currentLineSplit = currentLine.split(" ")
                if(float(currentLineSplit[0]) < bestRMSD):
                    bestRMSD = float(currentLineSplit[0])
                    lineValues[0] = currentLine
                    j = i + 6
                    k = 0
                    while(i != j):
                        k += 1
                        i += 1
                        lineValues[k] = linesReversed[i]

lineValuesReversed = []
lineValues = reversed(lineValues)
for value in lineValues:
    lineValuesReversed.append(value)

with open('resultsGeneralizedParallel', 'w') as resultsFile:
    for line in lineValuesReversed:
        resultsFile.write(line)

fitting_factor_combination = [0, 0, 0]
for i in range(3):
    fitting_factor_combination[i] = float(lineValuesReversed[i].split(" ")[0])

modelInpValues = [[45, 0, "TRIAL_NU"], [46, 0, "ION_NU"], [47, 0, "IONION_NU"]]
for i in range(3):
    modelInpValues[i][1] = float(lineValuesReversed[i + 3].split(" ")[0])

infile = open("./parameter_inputs_template.dat",'r')
outfile = open("./photo_processes.dat",'w')
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

modify_modelInp_values(modelInpValues, ".")

print("Best fit parameters: ")
for line in lineValuesReversed:
    print(line)

print("\n\n\nRunning model with best fit parameters...")
os.system('./run.sh > /dev/null') # Run model, deal with files, and silence output

# Create the plot
os.system("python3 plotting.py")