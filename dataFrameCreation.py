from multiprocessing.sharedctypes import Value
import pandas as pd
import numpy as np
import csv
import os

species_list = ['be-', 'bO', 'bO2', 'bO2-', 'bO2*', 'bO2+', 'bO3', 'bO3-', 'bO3*', 'bO3+', 'bO-', 'bO*', 'bO+', \
                'e-', 'G0', 'G-', 'ge-', \
                'gO', 'gO2', 'gO2-', 'gO2*', 'gO2+', 'gO3', 'gO3-', 'gO3*', 'gO3+', 'gO-', 'gO*', 'gO+', \
                'O', 'O2', 'O2-', 'O2+', 'O3', 'O3-', 'O3+', 'O-', 'O+', \
                'total_ice_O', 'total_ice_O2', 'total_ice_O3']
bulk_list = ['be-', 'bO', 'bO2', 'bO2-', 'bO2*', 'bO2+', 'bO3', 'bO3-', 'bO3*', 'bO3+', 'bO-', 'bO*', 'bO+']
bulk_list_woions = ['bO','bO*', 'bO2', 'bO2*', 'bO3', 'bO3*']
all_list_woions = ['bO', 'bO2', 'bO3', 'gO', 'gO2', 'gO3']
bulk_list_2 = ['be-', 'be-*', 'bO-', 'bO-*', 'bO', 'bO+', 'bO+*', 'bO2-','bO2-*', 'bO2', 'bO2+', 'bO2+*', 'bO2*', 'bO3-', 'bO3-*', 'bO3', "bO3+", 'bO3+*', 'bO3*', 'bO*'] #The suprathermal species disappeared from my folder
bulk_ion_list = ['be-', 'be-*', 'bO-', 'bO-*', 'bO+', 'bO+*', 'bO2-','bO2-*', 'bO2+', 'bO2+*', 'bO3-', 'bO3-*', "bO3+", 'bO3+*']
ion_list = ['bO+', 'bO2+', 'bO3+'] # for testing code

version = '.'

## For new code outputs

# Change to the list you want the data for
wk_lst = bulk_list_2

# Rewrite csv file in order to modify the header, must be done to then read into a dataframe. I used code from  https://stackoverflow.com/questions/16306819/python-edit-csv-headers 
i = 0
while i < len(wk_lst) :
    inputFileName = version + "/csv/" + wk_lst[i] + ".csv"
    outputFileName = os.path.splitext(inputFileName)[0] + "_modified.csv"
    
    with open(inputFileName, newline='') as inFile, open(outputFileName, 'w', newline='') as outfile:
        r = csv.reader(inFile)
        w = csv.writer(outfile)
        
        next(r, None)  # skip the first row from the reader, the old header
        # write new header
        w.writerow(['Fluence', wk_lst[i]])
        
        # copy the rest
        for row in r:
            w.writerow(row)
        
    i += 1

# read and merge data
i=0

while i < len(wk_lst) :
    csv = pd.read_csv(version + "/csv/" + wk_lst[i] + "_modified.csv")
    if i == 0:
        merged_data = csv
    else:
        merged_data = merged_data.merge(csv, on=["Fluence"])
    i += 1
    
# Some data analysis
merged_data['Ion volume density'] = merged_data[bulk_ion_list].sum(axis=1)
merged_data['Total volume density'] = merged_data[bulk_list_2].sum(axis=1)
merged_data['Percent Ion'] = merged_data['Ion volume density']/merged_data['Total volume density']



# Save dataframe
merged_data.to_pickle(version + '/pickle_dataframes/csv_dataframe.pkl')
print(merged_data.head(n=25))

# Defining values needed to read the analytics files
# with ions
ana_list = ['bO', 'bO2', 'bO2*', 'bO3', 'bO3*', 'bO*', 'gO', 'gO2', 'gO2*', 'gO3', 'gO3*', 'gO*']
num_rxn = [21, 36, 29, 14, 24, 25, 18, 33, 10, 16, 7, 11]

# #w/o ions
# ana_list = ['bO', 'bO2', 'bO3']
# num_rxn = [19, 24, 8]


# Defining the fixed width columns
flu_col_lbl = ['NA', 'Fluence']
flu_specs = [(0,22), (22,-1)]
column_label = ['Index', 'Rxn', 'R1', 'R2', 'P1', 'P2', 'P3', 'D1', 'D2', 'D3', 'D4', 'D5']
col_specs = [(0, 4), (4, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50,80), (80, 91), (91, 98), (98, 103), (103, 116), (116, -1)]

# A comand that returns true if the row is not a fluence row
def logic1(index):
    if (index+1) % flu_row == 0:
       return False
    return True
# A comand that returns true if the row is a fluence row
def logic2(index):
    if (index+1) % flu_row == 0:
       return True
    return False

i=0
while i < len(ana_list) :
    print("\n")
    print(i)
    print(ana_list[i])
    print(len(ana_list))
    filename = version + "/analytics_" + ana_list[i]
    flu_row = num_rxn[i] + 1
    
    # Read only the fluence rows
    temp_df = pd.read_fwf(filename, \
                          skiprows= lambda x: logic1(x), \
                          names=flu_col_lbl, \
                          colspecs=flu_specs \
                         )
    # Duplicate each row by number of reactions
    temp_df_2 = pd.concat([temp_df] * num_rxn[i], ignore_index = True)
    # Sort so all x repeats are sequestial
    flu_df = temp_df_2.sort_values(by=['Fluence'], ignore_index = True)
    
    # Read only the non-fluence rows 
    df_rxn = pd.read_fwf(filename, \
                         skiprows = lambda x: logic2(x), \
                         colspecs=col_specs, \
                         names=column_label \
                        )
    # Merge the dataframes and delete unnecessary columns
    merged_df = flu_df.join(df_rxn)
    del merged_df['NA']
    del merged_df['Index']
    del merged_df['R1']
    del merged_df['R2']
    del merged_df['P1']
    del merged_df['P2']
    del merged_df['P3']
    del merged_df['D1']
    del merged_df['D3']
    del merged_df['D4']
    del merged_df['D5']
    
    j=0
    plot_rxn_list = []
    
    while j < num_rxn[i]:
        # Get reaction number from dataframe
        rxn_num = merged_df['Rxn'].values[j] 
        plot_rxn_list.append(rxn_num)
        # Create temporary df of only one reaction
        temp_df = merged_df[merged_df['Rxn'] == rxn_num]
        # Rename D2 to the reaction number
        temp2_df = temp_df.rename(columns={"D2": rxn_num})
        # Delete rxn column
        del temp2_df['Rxn']
        # Merge temp dataframe into reaction dataframe
        if j == 0: 
            rxn_data = temp2_df
        else:
            rxn_data = rxn_data.merge(temp2_df, on=["Fluence"])
        
        j += 1
        
    #Drop all reactions that contribute less than 5% to the rate at all times
    k=0
    
    while k < len(plot_rxn_list):
        try:
            print(rxn_data[plot_rxn_list[k]])
            if ((abs(rxn_data[plot_rxn_list[k]]) < 5).all()) :
                #print('Reaction', plot_rxn_list[k], 'never contributes more than 5% to the rate.')
                rxn_data.drop(plot_rxn_list[k] , inplace=True, axis=1)
    #         else :
    #             # Use the else bit if you want only the low contributing reactions
    #             rxn_data.drop(plot_rxn_list[k] , inplace=True, axis=1)
        except TypeError: # Error; value of rxn_data[plot_rxn_list[k]] is not a number so it can't be added. 
            pass
        k += 1
        
    #print(plot_rxn_list)
    rxn_data.to_pickle(version + '/pickle_dataframes/' + ana_list[i]+'_rxn_dataframe.pkl')
    print(ana_list[i])
    print(rxn_data.head(n=5))
    i += 1
    
#rxn_data
#print(plot_rxn_list)