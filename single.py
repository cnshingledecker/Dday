import pandas as pd
import numpy as np
from itertools import product
import re
import subprocess
from exportable_custom_functions import split_list, split_list_chunks, find_nearest_index,is_float, is_int, modify_modelInp_values, get_data_to_modify_modelInp, setup_experimental_data, format_data_with_spaces


###Switches###
Interactive = True

###Key Parameters###
initialO2 = 5.7E22 # O2 molecules cm^-3
flux = 2.2E14 #particles cm^-2 s^-1


expData = setup_experimental_data() # The experimental data we compare the model to
print(expData)
if Interactive:
    input("Press Enter to continue...")


RMSDvals = [0.0]


# Modify model.inp
# Open model.inp to read and modify its content
with open("model.inp", "r") as file:
  lines = file.readlines()  # Read all lines into a list
  # If interactive, print relavant parameters

# Open photo_processes.dat and modify
with open("photo_processes.dat", "r") as file:
  lines = file.readlines()  # Read all lines into a list
  # If interactive, print relavant parameters

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

RMSDvals = (sum(deviations) / len(expData)) ** 0.5  # Formula for RMSD
print("RMSD = ", RMSDvals)
