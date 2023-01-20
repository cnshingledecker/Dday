# Run DataFrameCreation.py before this

import csv
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from exportable_custom_functions import setup_experimental_data

# Get ozone from file bO3.csv, do comparison of data based on Dr. Shingledecker email. We only need to compare this data, not anything else 
#                                                                                      (look at the below code to make sure this is the case).

#Set custom default colors
custom_cycler = (cycler(color=['#098ec3', #blue
                               '#ec8013', #orange
                               '#53b234', #green
                               '#df2020', #red
                               '#8463e9', #indigo
                               '#ccac66', #sand
                               '#ff8680', #pink
                               '#c170db', #violet
                               '#acd926', #yellowgreen
                               '#00CED1', #cyan
                               '#ec4913', #redorange
                               '#805e00' #brown
                              ]) 
                )

# Change style defaults
# plt.rc('text', usetex=True) # Allows LaTeX formatting with r'Formatted Text'
# plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']}) # Changes font to default LaTeX fount
plt.rcParams['axes.prop_cycle'] = custom_cycler # Sets plot color cycler to the cycle defined above

# Define filepath
# wo_ions_version = 'wo_ions/old_output' # Mullikin results
w_ions_version = 'csv/bO3.csv'
w_ions_modified_version = "bO3_for_dataframe.csv"

# Format csv file correctly with correct column names; we do this by writing the results to another file
with open(w_ions_version, newline='') as bO3_csv:  # Read in the lines from the file
   reader = csv.reader(bO3_csv, delimiter=',')      
   fields = next(reader) # Skip the line with the bad header
   with open(w_ions_modified_version, 'w') as bO3_for_dataframe_csv:
      bO3_csv_writer = csv.writer(bO3_for_dataframe_csv)
      bO3_csv_writer.writerow(["Fluence", "bO3"])
      try:
         while(True):
            row = next(reader)
            if(float(row[0]) != 0 and float(row[1]) != 0): # Gets rid of a value that sometimes occurs that causes weird plotting results for bO3.csv
               bO3_csv_writer.writerow(row)
      except StopIteration:
         pass # We are done processing lines; this is the exception (expected) which is thrown when we are out of lines to write to the file. 

# Load in data
# woions_df = pd.read_pickle(wo_ions_version + '/pickle_dataframes/csv_dataframe.pkl')
wions_df = pd.read_csv(w_ions_modified_version)

os.remove("./bO3_for_dataframe.csv") # Because we no longer need the file


# Isolate data of interest
# model_data_woions = woions_df[['Fluence', 'total_O3', 'total_O2']]
model_data_wions = wions_df[['Fluence', 'bO3']]

initialO2 = 5.7E22

exp_data = setup_experimental_data()

# model_data_woions["Fluence"] = model_data_woions["Fluence"]

model_data_wions["Fluence"] = model_data_wions["Fluence"]
model_data_wions["bO3"] = (model_data_wions["bO3"] / initialO2) * 100

# Create figure

# Add data to figure
ax = model_data_wions.plot(x = 'Fluence', y = 'bO3', label = 'Model with ions')
# model_data_woions.plot(ax = ax, x = 'Fluence', y = 'bO3', label = 'Without ions: Mullikin et al.')
exp_data.plot.scatter(ax = ax, x = 'expX', y = 'expY', marker="x", c='dimgrey', label = 'Experiment: Gerakines et al.') #  c='#53b234'

# Add plot elements
plt.legend()
plt.yscale("log")
plt.xscale("log")

plt.xlim([8e14,1e18])
plt.ylim([0.7,30])

plt.xlabel(r'Fluence $\left( \frac{particles}{cm^2} \right)$')
plt.ylabel(r'$[O_3]/[O_2]_{initial} \times 100 \%$')
plt.title(r'$O_3$ Production During $O_2$ Irradiation', fontsize=16)

# Show plot
# plt.savefig(w_ions_version + '/figures/f2.png')
# plt.savefig(w_ions_version + '/figures/f2.pgf')
# plt.savefig(w_ions_version + '/figures/f2.eps')
# plt.show() # Because not displaying it in a ipynb file
plt.savefig('plot.jpg', bbox_inches='tight')