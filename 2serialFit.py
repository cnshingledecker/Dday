import pandas as pd
import numpy as np
from itertools import product
import re
import subprocess
from exportable_custom_functions import split_list, split_list_chunks, find_nearest_index,is_float, is_int, modify_modelInp_values, get_data_to_modify_modelInp, setup_experimental_data, format_data_with_spaces


###Switches###
varyModel = True
varyPhotoProcesses = True
Interactive = False

###Key Parameters###
initialO2 = 5.7E22 # O2 molecules cm^-3
flux = 2.2E14 #particles cm^-2 s^-1

# Initialize Data Structures
model_inp_ranges = pd.DataFrame()
photo_processes_ranges = pd.DataFrame()
print("DataFrames initialize")

expData = setup_experimental_data() # The experimental data we compare the model to
print(expData)
if Interactive:
    input("Press Enter to continue...")

# Read model.inp variables to adjust and confirm min, max, and step
if varyModel:
    model_inp_file = "./model_inp_values/model_ranges.dat"
    model_inp_ranges = pd.read_csv(model_inp_file)
    print(model_inp_ranges)
    print(model_inp_ranges.shape)
    print(model_inp_ranges.loc[0,"toModify"])
    print(model_inp_ranges.index)
    print("Model inp ranges set")
    if Interactive:
        input("Press Enter to continue...")

# Read photo_processes.dat variables to adjust and confirm min, max, and step
if varyPhotoProcesses:
    photo_processes_file = "./photo_processes_values/photo_processes_ranges.dat"
    photo_processes_ranges = pd.read_csv(photo_processes_file)
    print(photo_processes_ranges)
    print(photo_processes_ranges.shape)
    print(photo_processes_ranges.loc[0,"toModify"])
    print(photo_processes_ranges.index)
    print("Photo processes ranges set")
    if Interactive:
        input("Press Enter to continue...")

# Generate combinations
variables = pd.concat([model_inp_ranges, photo_processes_ranges], ignore_index=True)
print(variables.loc[0])
if Interactive:
    input("Press Enter to continue...")

# Generate ranges for each parameter in log-space
print("Generating ranges")
ranges = {
    row['toModify']: np.logspace(
        np.log10(row['minValue']),
        np.log10(row['maxValue']),
        int(row['Steps'])
    )
    for _, row in variables.iterrows()
}

# Generate all combinations of the ranges
print("Determining combinations")
combinations = list(product(*ranges.values()))

# Create the result DataFrame
parameterSets = pd.DataFrame(combinations, columns=ranges.keys())
RMSDvals = [0.0] * len(parameterSets.index)
print(parameterSets.loc[0])
print("There are ",len(parameterSets.index)+1," parameter sets")
if Interactive:
    input("Press Enter to continue...")

# Create a data frame for fitness values
fitness_results = pd.DataFrame(columns=["RMSD"], index=parameterSets.index)

# Modify model.inp
print("Now looping through parameter sets")
for index in parameterSets.index:
    print("***************************************************************************")
    print("Parameter Set",index," of ",parameterSets.shape[0])
    print(parameterSets.loc[index])
    for column in parameterSets.columns:
        print(column)
        # Open model.inp to read and modify its content
        with open("model.inp", "r") as file:
            lines = file.readlines()  # Read all lines into a list
        with open("model.inp", "w") as file:
            for line in lines:
                if re.search(rf"\b{column}\b", line):  # Using raw f-string (rf-string) for dynamic regex
                    print(f"Old Line: {line}")
                    # New float value to replace
                    new_value = parameterSets.loc[index,column]  # Replace with your desired float
                    # Format the new value as scientific notation with the same format (e.g., 1.0000E+15)
                    formatted_value = f"{new_value:.4E}"
                    # Regular expression to match the exponential number
                    pattern = r"=\s+[-+]?\d+\.\d+E[-+]?\d+"
                    # Replace the matched value with the new formatted value
                    line = re.sub(pattern, f"= {formatted_value}", line)
                    print(f"New Line: {line}")
                    file.write(line)
                else:
                    file.write(line)

        # Open photo_processes.dat and modify
        with open("photo_processes.dat", "r") as file:
            lines = file.readlines()  # Read all lines into a list
        with open("photo_processes.dat", "w") as file:
            for line in lines:
                if re.search(rf"\b{column}\b", line):  # Using raw f-string (rf-string) for dynamic regex
                    print(f"Old Line: {line}")
                    # New float value to replace
                    new_value = parameterSets.loc[index,column]  # Replace with your desired float
                    # Format the new value as scientific notation with the same format (e.g., 1.0000E+15)
                    formatted_value = f"{new_value:.2E}"
                    # Regular expression to match the exponential number
                    pattern = r"([-+]?\d+\.\d+E[-+]?\d+)\s*$"
                    # Replace the matched value with the new formatted value
#                    line = re.sub(pattern, f"= {formatted_value}", line)
                    start_col = 107  # Starting column (1-based index)
                    end_col = 118  # Ending column (1-based index)
                    line = line[:start_col - 1] + formatted_value + line[end_col:] + "\n"
                    print(f"New Line: {line}")
                    file.write(line)
                else:
                    file.write(line)
    # Now run monaco
    subprocess.run("./monaco")

    # Import bO3.csv and read into df
    calcData = pd.read_csv("csv/bO3.csv", header=1, names=["Fluence", "Abundance"])
    calcData["Abundance"] = (calcData["Abundance"] / initialO2) * 100.0

    # Find the closest calculated values to experiment and calculate deviation
    deviations = [0.0] * len(calcData.index)
    for ind in expData.index:
        # Find the index of the nearest entry
        scalar = expData.loc[ind, "expX"]
        closest_index = (calcData["Fluence"] - scalar).abs().idxmin()
        closest_value = calcData.loc[closest_index, "Abundance"]
        deviations[ind] = (closest_value - expData.loc[ind, "expY"]) ** 2

    RMSDvals[index] = (sum(deviations) / len(expData)) ** 0.5  # Formula for RMSD
    if Interactive:
        print("RMSD = ", RMSDvals[index])
        input("Press Enter to continue...")

# Now determine which parameter sets yielded the lowest RMSD
bestRMSD = min(RMSDvals)
bestRMSDindex = RMSDvals.index(bestRMSD)

print("Best RMSD = ", bestRMSD," for parameter set")
print(parameterSets.loc[bestRMSDindex])
parameterSets.to_csv("RMSD_vals.out")


