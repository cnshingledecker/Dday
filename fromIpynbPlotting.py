# %% [markdown]
# ### How it works:
# 
# Data is cleaned elsewhere and saved using ```pandas.DataFrame.to_pickle('file_name.pkl')``` and are loaded in here using ```pandas.DataFrame.read_pickle('file_name.pkl')```. Values are then scaled and plotted as needed. 
# 
# Universal plot settings, such as font and the color cycler, are set in 2nd cell below this.
# 
# I will want to use .pgf to import the plots into LaTeX: https://riptutorial.com/matplotlib/example/10066/saving-and-exporting-plots-that-use-tex

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from exportable_custom_functions import setup_experimental_data
from baragiola_file_and_data_functions import getFlux, getInitialO2, process_model_data

onlyPaperPlots = True

# %%
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
plt.rc('text', usetex=True) # Allows LaTeX formatting with r'Formatted Text'
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']}) # Changes font to default LaTeX fount
plt.rcParams['axes.prop_cycle'] = custom_cycler # Sets plot color cycler to the cycle defined above

# %% [markdown]
# ## Define filepath
# 
# This is the filepath from the current directory to the directory where your data is stored. This directory should contain the subdirectories: csv, analytics, pickle_dataframes, paper_figures

# %%
# Define filepath
wo_ions_version = 'wo_ions/old_output' # Mullikin results
w_ions_version = '.'

# %% [markdown]
# ## 2test_fig2.pro: 
# 
# This code is for the model and network used in my thesis
# 
# See the Jupyter notebook file for the code if needed

# %% [markdown]
# ## figure_2_DO_NOT_MODIFY.pro
# 
# This code is for the model and network published in Mullikin et al.
# 
# See the Jupyter notebook file for the code if needed

# %% [markdown]
# ## Plotting comparison to Gerakines et al.
# 
# This cell recreated the plotting scripts from 2test_fig2.pro (New code) and figure_2_DO_NOT_MODIFY.pro (published data) to compare the fit of the current run to the previously published fit and the experimental data

# %%
# Load in data
woions_df = pd.read_pickle(wo_ions_version + '/pickle_dataframes/csv_dataframe.pkl')
wions_df = pd.read_pickle(w_ions_version + '/pickle_dataframes/csv_dataframe.pkl')

# Set constants: Variables to be taken from model.inp
# ICE_AREA = 1.0e-20 # cm2
# ICE_THICK = 1.0e-4 # cm
# ICE_DENSITY = 5.7e22 # molecules/cm3 Called RHO_ICE in model.inp
# PHI_EXP = 1.0e15

initialO2 = getInitialO2() # FLUX WAS TWICE THE VALUE USED FOR THE OTHER PAPER PLOTS. NOTE THIS IN THE EMAIL TO DR. SHINGLEDECKER.
flux = getFlux()

# Calculate initial O2 
# initialO2 = ICE_DENSITY * ICE_AREA * ICE_THICK # # of initial O2 molecules
# initial_oxygen = 2 * ICE_DENSITY * ICE_AREA * ICE_THICK # # of initial O atoms (in O2 molecules)

# Isolate data of interest
model_data_woions = woions_df[['Fluence', 'total_O3', 'total_O2']]
model_data_wions = wions_df[['Fluence', 'bO3']]

# Scale modifications
# exp_data["expY"] = exp_data["expY"] * 15 # 15 is an artifact of the digitization of the gerakines data. I can confirm with the plot in the paper that the dta is plotted correctly.
# exp_data["expX"] = exp_data["expX"] * flux

# Note: flux from the Jupyter notebook file is 2.2e14, but the one used for the model runs for the Ion ice paper 
#     (which is being used below) is 2.33e14. 
#     MAKE SURE TO ASK IF IT'S OKAY IN THE EMAIL TO DR. SHINGLEDECKER.
exp_data = setup_experimental_data()

model_data_woions["Fluence"] = model_data_woions["Fluence"]
#model_data_woions["bO3"] = 3 * (model_data_woions["bO3"] / initialO2) * 100
model_data_woions["bO3"] = (model_data_woions['total_O3'] / model_data_woions.at[0,'total_O2']) * 100

model_data_wions["Fluence"] = model_data_wions["Fluence"]
model_data_wions["bO3"] = process_model_data(model_data_wions["bO3"]) # See email from Kristen for the reason the bO3 was multiplied by 3

#model_data.head()
# Create figure

# Add data to figure
ax = model_data_wions.plot(x = 'Fluence', y = 'bO3', label = 'Model with ions')
model_data_woions.plot(ax = ax, x = 'Fluence', y = 'bO3', label = 'Without ions: Mullikin et al.')
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

# Save plot (of the abundance of O3 relative to the initial O2 vs. fluence)
plt.savefig(w_ions_version + '/paper_figures/O3AbundanceVsFluence.jpg')


# %% [markdown]
# ## Percent Ion
# 
# Plots the percentage of ionic species in the model with respect to fluence
#     This plot wasn't needed for the ion ice paper, so the code was removed (still in the )


# %% [markdown]
# ## Ion Abundances
# 
# Plots the percentage of ionic species in the model with respect to fluence

# %%
# Load in data
csv_df = pd.read_pickle(w_ions_version + '/pickle_dataframes/csv_dataframe.pkl')
initial_O2 = getInitialO2()

# Create figure
#fig = plt.figure(figsize=[12.8, 9.6])

# Set legend labels list
plot_ion_list = ['Ion volume density', 'be-', 'bO+', 'bO-', 'bO2+', 'bO2-', 'bO3+', 'bO3-']
plot_ion_labels = [r'Ion Sum', r'$e^-$', r'$O^+$', r'$O^-$', r'$O_2^+$', r'$O_2^-$', r'$O_3^+$', r'$O_3^-$']

# Scale modifications
i=0
while i < len(plot_ion_list):
    #csv_df[plot_ion_list[i]] = csv_df[plot_ion_list[i]] / csv_df['Total volume density']
    csv_df[plot_ion_list[i]] = process_model_data(csv_df[plot_ion_list[i]])
    i += 1

# Add the data
#ax = csv_df.plot(x = 'Fluence', y = 'Percent Ion', legend = 'All ions')
csv_df.plot(x = 'Fluence', y = plot_ion_list[:], label = plot_ion_labels[:])

# Add plot elements
plt.legend()
#plt.yscale("log")
plt.xscale("log")

# plt.xlim([1e14,1e17])
# #plt.ylim([1e-5,1e1])

plt.xlabel(r'Fluence $\left( \frac{particles}{cm^2} \right)$')
plt.ylabel(r'$[A]/[O_2]_{initial} \times 100 \%$')
plt.title(r'Ion Abundances in the Bulk over the Model Time',fontsize=16)

# Show plot
plt.savefig(w_ions_version + '/paper_figures/IonAbundance.jpg')

# %% [markdown]
# ## Reaction Contribution
# 
# Reaction labels and reaction types must be confirmed after new runs! (Likely only after modifying input files, but don't assume.) Because the code assigns reaction numbers, there is no guarentee that these labels are current. 

# %%
# Labels for the numberes reactions

# Highest reaction number + 1 in the spc_list[i]_rxn_dataframe.pkl
last_rxn = 186

Rxn_label = pd.DataFrame({'Rxn': list(range(0,last_rxn,1)),
                          'Label': [(r'Add reaction ' + str(i)) for i in range(last_rxn)],
                          'Type': [0 for i in range(last_rxn)]
                          })

# Rxn_label.at[0,'Label'] = r'$ \rightarrow $' #
# Rxn_label.at[0,'Type'] = 0 #

#Types:
# 0 = Neutral-Neutral reactions
# 1 = ion-neutral reactions
# 2 = ion recombination reactions
# 3 = photoprocesses
# 4 = suprathermal/quenching

Rxn_label.at[56,'Label'] = r'$2 O \rightarrow O_2$' #56
Rxn_label.at[56,'Type'] = 0 #56
Rxn_label.at[58,'Label'] = r'$O + e^- \rightarrow O^-$' #58
Rxn_label.at[58,'Type'] = 1 #58
Rxn_label.at[59,'Label'] = r'$O + O^- \rightarrow O_2^-$' #59
Rxn_label.at[59,'Type'] = 1 #59
Rxn_label.at[60,'Label'] = r'$O + O_2^- \rightarrow O_3^-$' #60
Rxn_label.at[60,'Type'] = 1 #
Rxn_label.at[61,'Label'] = r'$O^+ + e^- \rightarrow O^*$' #61
Rxn_label.at[61,'Type'] = 1 #
Rxn_label.at[62,'Label'] = r'$O^+ + O^- \rightarrow O_2^*$' #62
Rxn_label.at[62,'Type'] = 2 #
Rxn_label.at[63,'Label'] = r'$O^+ + O^- \rightarrow 2 O^*$' #63
Rxn_label.at[63,'Type'] = 2 #
Rxn_label.at[64,'Label'] = r'$O^+ + O_2 \rightarrow O_2^+ + O$' #64
Rxn_label.at[64,'Type'] = 1 #
Rxn_label.at[65,'Label'] = r'$O^+ + O_2^- \rightarrow O_3^*$' #65
Rxn_label.at[65,'Type'] = 2 #
Rxn_label.at[66,'Label'] = r'$O^+ + O_2^- \rightarrow O^* + O_2^*$' #66
Rxn_label.at[66,'Type'] = 2 #
Rxn_label.at[67,'Label'] = r'$O^+ + O_3^- \rightarrow O^* + O_3^*$' #67
Rxn_label.at[67,'Type'] = 2 #
Rxn_label.at[68,'Label'] = r'$O^+ + O_3^- \rightarrow 2 O_2^*$' #68
Rxn_label.at[68,'Type'] = 2 #
Rxn_label.at[69,'Label'] = r'$O_2 + e^- \rightarrow O_2^-$' #69
Rxn_label.at[69,'Type'] = 1 #
Rxn_label.at[70,'Label'] = r'$O_2 + O^- \rightarrow O_3^-$' #70
Rxn_label.at[70,'Type'] = 1 #
Rxn_label.at[71,'Label'] = r'$O_2^+ + e^- \rightarrow O_2^*$' #71
Rxn_label.at[71,'Type'] = 2 #
Rxn_label.at[72,'Label'] = r'$O_2^+ + e^- \rightarrow 2 O^*$' #72
Rxn_label.at[72,'Type'] = 2 #
Rxn_label.at[73,'Label'] = r'$O_2^+ + O^- \rightarrow O_3^*$' #73
Rxn_label.at[73,'Type'] = 2 #
Rxn_label.at[74,'Label'] = r'$O_2^+ + O^- \rightarrow O_2^* + O^*$' #74
Rxn_label.at[74,'Type'] = 2 #
Rxn_label.at[75,'Label'] = r'$O_2^+ + O_2 \rightarrow O_3^+ + O$' #75
Rxn_label.at[75,'Type'] = 1 #
Rxn_label.at[76,'Label'] = r'$O_2^+ + O_2^- \rightarrow O_3^* + O^*$' #76
Rxn_label.at[76,'Type'] = 2 #
Rxn_label.at[77,'Label'] = r'$O_2^+ + O_2^- \rightarrow 2 O_2^*$' #77
Rxn_label.at[77,'Type'] = 2 #
Rxn_label.at[78,'Label'] = r'$O_2^+ + O_3^- \rightarrow O_2^* + O_3^*$' #78
Rxn_label.at[78,'Type'] = 2 #
Rxn_label.at[79,'Label'] = r'$O_2^+ + O_3^- \rightarrow 2 O_2^* + O^*$' #79
Rxn_label.at[79,'Type'] = 2 #
Rxn_label.at[80,'Label'] = r'$O_2^+ + O_3^- \rightarrow O_3^* + 2 O^*$' #80
Rxn_label.at[80,'Type'] = 2 #
Rxn_label.at[81,'Label'] = r'$O_3 + e^- \rightarrow O_3^-$' #81
Rxn_label.at[81,'Type'] = 1 #
Rxn_label.at[82,'Label'] = r'$O_3 + O^- \rightarrow O_2 + O_2^-$' #82
Rxn_label.at[82,'Type'] = 1 #
Rxn_label.at[83,'Label'] = r'$O_3 + O_2^- \rightarrow O_2 + O_3^-$' #83
Rxn_label.at[83,'Type'] = 1 #
Rxn_label.at[84,'Label'] = r'$O_3 + O_3^+ \rightarrow 2 O_2 + O_2^+$' #84
Rxn_label.at[84,'Type'] = 1 #
Rxn_label.at[85,'Label'] = r'$O_3^- + O^- \rightarrow 2 O_2^-$' #85
Rxn_label.at[85,'Type'] = 2 #
Rxn_label.at[86,'Label'] = r'$O_3^+ + e^- \rightarrow O_3^*$' #86
Rxn_label.at[86,'Type'] = 2 #
Rxn_label.at[87,'Label'] = r'$O_3^+ + e^- \rightarrow O_2^* + O^*$' #87
Rxn_label.at[87,'Type'] = 2 #
Rxn_label.at[88,'Label'] = r'$O_3^+ + O^- \rightarrow O^* + O_3^*$' #88
Rxn_label.at[88,'Type'] = 2 #
Rxn_label.at[89,'Label'] = r'$O_3^+ + O^- \rightarrow 2 O_2^*$' #89
Rxn_label.at[89,'Type'] = 2 #
Rxn_label.at[90,'Label'] = r'$O_3^+ + O_2 \rightarrow O_2^+ + O_3$' #90
Rxn_label.at[90,'Type'] = 1 #
Rxn_label.at[91,'Label'] = r'$O_3^+ + O_2 \rightarrow O_3^+ + O_2$' #91
Rxn_label.at[91,'Type'] = 1 #
Rxn_label.at[92,'Label'] = r'$O_3^+ + O_2 \rightarrow 2 O_2 + O^+$' #92
Rxn_label.at[92,'Type'] = 1 #
Rxn_label.at[93,'Label'] = r'$O_3^+ + O_2^- \rightarrow O_3^* + O_2^*$' #93
Rxn_label.at[93,'Type'] = 2 #
Rxn_label.at[94,'Label'] = r'$O_3^+ + O_2^- \rightarrow 2 O_2^* + O^*$' #94
Rxn_label.at[94,'Type'] = 2 #
Rxn_label.at[95,'Label'] = r'$O_3^+ + O_2^- \rightarrow O_3^* + 2 O^*$' #95
Rxn_label.at[95,'Type'] = 2 #
Rxn_label.at[96,'Label'] = r'$O_3^+ + O_3^- \rightarrow 2 O_3^*$' #96
Rxn_label.at[96,'Type'] = 2 #
Rxn_label.at[97,'Label'] = r'$O_3^+ + O_3^- \rightarrow O_3^* + O_2^* + O^*$' #97
Rxn_label.at[97,'Type'] = 2 #
Rxn_label.at[98,'Label'] = r'$O_3^+ + O_3^- \rightarrow 3 O_2^*$' #98
Rxn_label.at[98,'Type'] = 2 #
Rxn_label.at[101,'Label'] = r'$O^* + O \rightarrow O_2$' #101
Rxn_label.at[101,'Type'] = 4 #
Rxn_label.at[102,'Label'] = r'$O_3^* + O_3 \rightarrow 3 O_2$' #102
Rxn_label.at[102,'Type'] = 4 #
Rxn_label.at[117,'Label'] = r'$O^* + O \rightarrow 2 O$' #117
Rxn_label.at[117,'Type'] = 4 #
Rxn_label.at[118,'Label'] = r'$O^* + O \rightarrow O_2^*$' #118
Rxn_label.at[118,'Type'] = 4 #
Rxn_label.at[121,'Label'] = r'$O^* + O_2 \rightarrow O_3^*$' #121
Rxn_label.at[121,'Type'] = 4 #
Rxn_label.at[123,'Label'] = r'$O^* + O_2 \rightarrow O + O_2$' #123
Rxn_label.at[123,'Type'] = 4 #
Rxn_label.at[125,'Label'] = r'$O + O_3 \rightarrow 2 O_2$' #125
Rxn_label.at[125,'Type'] = 0 #
Rxn_label.at[127,'Label'] = r'$O^* + O_3 \rightarrow 2 O_2$' #127
Rxn_label.at[127,'Type'] = 4 #
Rxn_label.at[129,'Label'] = r'$O_2^* + O \rightarrow O_3^*$' #129
Rxn_label.at[129,'Type'] = 4 #
Rxn_label.at[131,'Label'] = r'$O_2^* + O \rightarrow O_2 + O$' #131
Rxn_label.at[131,'Type'] = 4 #
Rxn_label.at[133,'Label'] = r'$O_2^* + O_2 \rightarrow 2 O_2$' #133
Rxn_label.at[133,'Type'] = 4 #
Rxn_label.at[135,'Label'] = r'$O_2^* + O_3 \rightarrow 2 O_2 + O$' #135
Rxn_label.at[135,'Type'] = 4 #
Rxn_label.at[137,'Label'] = r'$2 O^* \rightarrow 2 O$' #137
Rxn_label.at[137,'Type'] = 4 #
Rxn_label.at[138,'Label'] = r'$O^* + O \rightarrow O_2^*$' #138
Rxn_label.at[138,'Type'] = 4 #
Rxn_label.at[143,'Label'] = r'$O^* + O_2^* \rightarrow O + O_2$' #143
Rxn_label.at[143,'Type'] = 4 #
Rxn_label.at[145,'Label'] = r'$O^* + O_3^* \rightarrow 2 O_2$' #145
Rxn_label.at[145,'Type'] = 4 #
Rxn_label.at[147,'Label'] = r'$2 O_2^* \rightarrow 2 O_2$' #147
Rxn_label.at[147,'Type'] = 4 #
Rxn_label.at[149,'Label'] = r'$O_2^* + O_3^* \rightarrow 2 O_2 + O$' #149
Rxn_label.at[149,'Type'] = 4 #
Rxn_label.at[151,'Label'] = r'$2 O_3^* \rightarrow 3 O_2$' #151
Rxn_label.at[151,'Type'] = 4 #
Rxn_label.at[153,'Label'] = r'$O + O_3^* \rightarrow 2 O_2$' #153
Rxn_label.at[153,'Type'] = 4 #
Rxn_label.at[155,'Label'] = r'$O_2 + O_3^* \rightarrow 2 O_2 + O$' #155
Rxn_label.at[155,'Type'] = 4 #
Rxn_label.at[156,'Label'] = r'$O + O_2 \rightarrow O_3^*$' #156
Rxn_label.at[156,'Type'] = 0 #
Rxn_label.at[158,'Label'] = r'$O^* \rightarrow O$' #158
Rxn_label.at[158,'Type'] = 4 #
Rxn_label.at[159,'Label'] = r'$O_2^* \rightarrow O_2$' #159
Rxn_label.at[159,'Type'] = 4 #
Rxn_label.at[160,'Label'] = r'$O_3^* \rightarrow O_3$' #160
Rxn_label.at[160,'Type'] = 4 #
Rxn_label.at[169,'Label'] = r'$O_2 \rightarrow 2 O$' #169
Rxn_label.at[169,'Type'] = 3 #
Rxn_label.at[171,'Label'] = r'$O_2 \rightarrow O^+ + O^-$' #171
Rxn_label.at[171,'Type'] = 3 #
Rxn_label.at[173,'Label'] = r'$O_2 \rightarrow O_2^*$' #173
Rxn_label.at[173,'Type'] = 3 #
Rxn_label.at[175,'Label'] = r'$O_2 \rightarrow 2 O^*$' #175
Rxn_label.at[175,'Type'] = 3 #
Rxn_label.at[177,'Label'] = r'$O_2 \rightarrow O_2^+ + e^-$' #177
Rxn_label.at[177,'Type'] = 3 #
Rxn_label.at[179,'Label'] = r'$O_3 \rightarrow O^+ + O_2^-$' #179
Rxn_label.at[179,'Type'] = 3 #
Rxn_label.at[181,'Label'] = r'$O_3 \rightarrow O_2 + O$' #181
Rxn_label.at[181,'Type'] = 3 #
Rxn_label.at[183,'Label'] = r'$O_3 \rightarrow O_2^+ + O^-$' #183
Rxn_label.at[183,'Type'] = 3 #
Rxn_label.at[185,'Label'] = r'$O_3 \rightarrow O_3^*$' #185
Rxn_label.at[185,'Type'] = 3 #

#Rxn_label

# %%
# Load in data
spc_list = ['be-', 'bO', 'bO-', 'bO+', 'bO2', 'bO2-', 'bO2+', 'bO3', 'bO3-', 'bO3+']
spc_label = [r'$e^-$', r'$O$', r'$O^-$', r'$O^+$', r'$O_2$', r'$O_2^-$', r'$O_2^+$', r'$O_3$', r'$O_3^-$', r'$O_3^+$']
flux = getFlux()

i=0
while i < len(spc_list) :
    rxn_data = pd.read_pickle(w_ions_version + '/pickle_dataframes/' + spc_list[i]+'_rxn_dataframe.pkl')
    # Get list of reaction numbers
    rxn_list = rxn_data.columns.values.tolist()
    rxn_list.pop(0)
    
    rxn_data["Fluence"] = rxn_data["Fluence"] * flux
    
    # Create figure
    #fig = plt.figure(figsize=[12.8, 9.6])
    
    # Add the data
    j=0
    while j < len(rxn_list) :
        if Rxn_label.loc[rxn_list[j],'Type'] == 1:
            style = '--'
        elif Rxn_label.loc[rxn_list[j],'Type'] == 2:
            style = '-.'
        elif Rxn_label.loc[rxn_list[j],'Type'] == 3:
            style = ':'
        elif Rxn_label.loc[rxn_list[j],'Type'] == 4:
            style = (0, (3, 2, 1, 2, 1, 2, 1, 2))
        else:
            style = '-'
        
        if j == 0: # Plot each reaction with the appropriate label
            ax = rxn_data.plot(x = 'Fluence', \
                               y = rxn_list[j], \
                               linestyle = style, \
                               label = Rxn_label.loc[rxn_list[j],'Label'])
        else:
            rxn_data.plot(ax = ax, \
                          x = 'Fluence', \
                          y = rxn_list[j], \
                          linestyle = style, \
                          label = Rxn_label.loc[rxn_list[j],'Label'])
        j += 1
        
    if spc_list[i] == 'bO2' :
        #plt.legend(loc='upper left', bbox_to_anchor=(1, 1)) # To the right of the plot
        plt.legend(ncol= 3, loc='upper center', bbox_to_anchor=(.5, -.15)) # Below the plot
    
    if spc_list[i] == 'bO3' :
        plt.legend(loc='best', bbox_to_anchor=(0, 0.3, 0.5, 0.5))
        
#     if len(rxn_list) > 10 : 
#         plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
#     else:
#         plt.legend(loc='best')

    # Add plot elements
    #plt.yscale("log")
    plt.xscale("log")
    
    #plt.axhline(y=0, color='k', linewidth = .5)
    #plt.axvline(x=1.5e15, color='k', linewidth = .5)
    plt.grid(True, axis='y')
       
#     plt.xlim([1e10,1e18])
#     #plt.ylim([-50, 50])
    
    #plt.legend(bbox_to_anchor=(1, 1))
    plt.xlabel(r'Fluence $\left( \frac{particles}{cm^2} \right)$')
    plt.ylabel(r'Percent Contribution of Reaction $X$ to Rate' 
               '\n'
               r'$\left(\frac{\mathrm{d}X}{\mathrm{d}t} \middle/ \sum^N \left| \frac{\mathrm{d}X_\mathrm{n}}{\mathrm{d}t}\right| \right) \times 100\%$', 
               multialignment='center')
    plt.title(r'Reaction Contributions to Rate of Formation of ' + spc_label[i], fontsize=16)
    
    # Save plot
    plt.savefig(w_ions_version + '/paper_figures/' + spc_list[i] + '_rxn_contribution.jpg', bbox_inches='tight')
    
    i += 1

# %%
# Creates a legend of the line styles (I am sure there is a better way of doing this)

from matplotlib.lines import Line2D

fig, ax = plt.subplots()
# 0 = Neutral-Neutral reactions
# 1 = ion-neutral reactions
# 2 = ion recombination reactions
# 3 = photoprocesses
# 4 = suprathermal/quenching

legend_elements = [Line2D([0], [0], color='k', ls='-', label= 'neutral-neutral reactions'),
                   Line2D([0], [0], color='k', ls='--', label= 'ion-neutral reactions'),
                   Line2D([0], [0], color='k', ls='-.', label= 'ion recombination reactions'),
                   Line2D([0], [0], color='k', ls=':', label= 'photoprocesses'),
                   Line2D([0], [0], color='k', ls=(0, (3, 2, 1, 2, 1, 2, 1, 2)), label= 'suprathermal/quenching reactions')]

legend = ax.legend(handles=legend_elements, loc='center', framealpha=1, frameon=True)

def export_legend(legend, filename=w_ions_version + '/paper_figures/legend.pgf', expand=[-5,-5,5,5]):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)

# %%



