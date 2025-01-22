import plotext as plt
import pandas as pd
import math

# Parameters
initialO2        = 5.699E22
flux             = 2.2E14
O3bandConversion = 24.0/1.6182984

# Data from Gerakines et al.
O3t = [4.709604, 13.496563, 50.07993, 146.55681, 465.09692, 2581.1648]
O3 = [0.055385794, 0.15238622, 0.52023816, 1.062651, 1.4836915, 1.6182984]

# Get data in correct form
O3t = [x*flux for x in O3t]
O3 = [x*O3bandConversion for x in O3]

# Read the CSV file
file_path = 'csv/bO3.csv'  # Replace with your file path
data = pd.read_csv(file_path,header=1)

# Display the first few rows
print(data.head())

# Access individual columns
x = data.iloc[:, 0]  # First column
y = data.iloc[:, 1]  # Second column
print(x.shape)

# Get y data in proper units
y = [(el/initialO2)*100.0 for el in y]

# Convert x-axes to log
x = [math.log10(el) for el in x]
O3t = [math.log10(el) for el in O3t]

# Plotext 
plt.xlabel("log10(Fluence)")
plt.ylabel("Percent O3")
plt.scatter(x,y,label="Calculated")
plt.scatter(O3t,O3, label="Experimental",color="red")
plt.xlim(O3t[0],O3t[-1])
plt.show()

print(O3t[0],O3t[-1])