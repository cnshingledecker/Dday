#! /usr/bin/python

import os
import csv
import subprocess
import numpy as np

results = open("resultsFile_3", 'w')

vectord = np.linspace(1.7, 2.7, num=10)
vectore = np.linspace(0, 2, num=3)
vectorf = np.linspace(0.25, 0.35, num=3)

# first for loop to modify fitting factor of O2->2O photoexc reaction
for d in range(0,10):
    dvald = vectord[d]
    #second for loop is to modify fitting factor of O2->2O* photoion reaction
    for e in range(1,2):
        dvale = vectore[e]
        #third for loop is to modify fitting factor for O2->O2*, O3->O2+O, and O3->O3* reactions
        for f in range(1,2):
            dvalf = vectorf[f]
            # fixed 'no such file' for parameter_inputs_template by adding .dat which is hidden in Folders
            infile = open("parameter_inputs_template.dat",'r')
            outfile = open("photo_processes_3.dat",'w')
            # we know ahead of time that d is 4TH and 7TH lines of file, e is 3RD and 6TH lines of file, and f is the 5TH, 8TH, 10TH, 11TH, 13TH, 14TH lines of the file so we are working under that assumption
            lines = infile.readlines()
            linecounter = 0
            for line in lines:
            #THIS IS ONLY FOR BARAGIOLA EXPERIMENT (due to cross sections); might automate switching between experiments once this works
                # (Daniel Lopez-Sanders comment) line identified by reactant (2nd thing in line) and products (4th (and sometimes 5th) thing in line)
                linecounter = linecounter + 1
                if (linecounter == 3): # e
                    if dvale >=1.00:
                        outfile.write(" 1467 gO2       PHOION              gO*       gO*                                     1.00E+00  0.00E-19  "+str(dvale)+"0E+00\n")
                    else:
                        dvalemod = dvale*10 
                        dvalemod = round(dvalemod, 2)
                        outfile.write(" 1467 gO2       PHOION              gO*       gO*                                     1.00E+00  0.00E-19  "+str(dvalemod)+"E-01\n")
                elif (linecounter == 4):   # d
                    outfile.write(" 1468 gO2       PHOEXC              gO        gO                                      1.00E+00  2.10E-19  "+str(dvald)+"0E+00\n")
                elif (linecounter == 5): # f first reaction
                    if dvalf >= 1.00:
                        outfile.write(" 1469 gO2       PHOEXC              gO2*                                              1.00E+00  2.10E-19  "+str(dvalf)+"E+00\n")
                    else:
                        dvalfmod = dvalf*10
                        dvalfmod = round(dvalfmod, 2)
                        outfile.write(" 1469 gO2       PHOEXC              gO2*                                              1.00E+00  2.10E-19  "+str(dvalfmod)+"E-01\n")
                elif (linecounter == 6): # e
                    if dvale >=1.00:
                        outfile.write(" 1467 bO2       PHOION              bO*       bO*                                     1.00E+00  0.00E-19  "+str(dvale)+"0E+00\n")
                    else:
                        dvalemod = dvale*10
                        dvalemod = round(dvalemod, 2)
                        outfile.write(" 1467 bO2       PHOION              bO*       bO*                                     1.00E+00  0.00E-19  "+str(dvalemod)+"E-01\n")
                elif (linecounter == 7):   # d
                    outfile.write(" 1468 bO2       PHOEXC              bO        bO                                      1.00E+00  2.10E-19  "+str(dvald)+"0E+00\n")
                elif (linecounter == 8): # f first reaction
                    if dvalf >= 1.00:
                        outfile.write(" 1469 bO2       PHOEXC              bO2*                                              1.00E+00  2.10E-19  "+str(dvalf)+"E+00\n")
                    else:
                        dvalfmod = dvalf*10
                        dvalfmod = round(dvalfmod, 2)
                        outfile.write(" 1469 bO2       PHOEXC              bO2*                                              1.00E+00  2.10E-19  "+str(dvalfmod)+"E-01\n")
                elif (linecounter == 10): # f second reaction
                    if dvalf >= 1.00:
                        outfile.write(" 1474 gO3       PHOEXC              gO2       gO                                      1.00E+00  2.15E-18  "+str(dvalf)+"E+00\n")
                    else:
                        dvalfmod = dvalf*10
                        dvalfmod = round(dvalfmod, 2)
                        outfile.write(" 1474 gO3       PHOEXC              gO2       gO                                      1.00E+00  2.15E-18  "+str(dvalfmod)+"E-01\n")
                elif (linecounter == 11): # f third reaction
                    if dvalf >= 1.00:
                        outfile.write(" 1474 gO3       PHOEXC              gO3*                                              1.00E+00  2.15E-18  "+str(dvalf)+"E+00\n")
                    else:
                        dvalfmod = dvalf*10
                        dvalfmod = round(dvalfmod, 2)
                        outfile.write(" 1474 gO3       PHOEXC              gO3*                                              1.00E+00  2.15E-18  "+str(dvalfmod)+"E-01\n")
                elif (linecounter == 13): # f second reaction
                    if dvalf >= 1.00:
                        outfile.write(" 1474 bO3       PHOEXC              bO2       bO                                      1.00E+00  2.15E-18  "+str(dvalf)+"E+00\n")
                    else:
                        dvalfmod = dvalf*10
                        dvalfmod = round(dvalfmod, 2)
                        outfile.write(" 1474 bO3       PHOEXC              bO2       bO                                      1.00E+00  2.15E-18  "+str(dvalfmod)+"E-01\n")
                elif (linecounter == 14): # f third reaction
                    if dvalf >= 1.00:
                        outfile.write(" 1474 bO3       PHOEXC              bO3*                                              1.00E+00  2.15E-18  "+str(dvalf)+"E+00\n")
                    else:
                        dvalfmod = dvalf*10
                        dvalfmod = round(dvalfmod, 2)
                        outfile.write(" 1474 bO3       PHOEXC              bO3*                                              1.00E+00  2.15E-18  "+str(dvalfmod)+"E-01\n")
                else:
                    outfile.write(line)
            outfile.close()
            infile.close()
            print("Running model...")
            # make a system call to model
            # fixed syntax error due to using os.system() with shell file
            # fixed OSError: [Errno 8] Exec format error by adding #!/bin/sh to top of run.sh file
            os.system('./runTest.sh')
            # open the output file created when the model runs
            # i got much of the following about reading csv from the internet so it may not be correct
            print("Finding RMSD...")
            # fixed [Errno 2] No such file or directory: 'total_ice_o3.csv' by directing to the csv folder
            with open('csv/total_ice_o3.csv') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                linecount = 0
                # find O3 abundance data corresponding to experimental data points
                for row in csv_reader:
                    # is there a way to be more concise? also not sure that it can handle all the decimal places?
                    # write O3 abundances to another file created OUTSIDE set of for loops
                    linecount = linecount + 1
                    if (linecount == 73):
                        # fixed IndexError: list index out of range because forgot indexing from 0; all the row[2] changed to row[1]
                        # fixed TypeError: unsupported operand type for - by defining a variable for the model data and changing the string to float type
                        model73 = float(row[1])*1e7
                        data73 = 1.236090330782555e15
                        dev73 = (model73 - data73)
                        if -4.56087068782555e14 <= dev73 <= 4.469206770810233e14:
                            dev73 = 0
                    elif (linecount == 75):
                        model75 = float(row[1])*1e7
                        data75 = 1.8348405368189518e15
                        dev75 = (model75 - data75)
                        if -3.882995488189518e14 <= dev75 <= 3.927024151806072e14:
                            dev75 = 0
                    elif (linecount == 77):
                        model77 = float(row[1])*1e7
                        data77 = 2.379956600877962e15
                        dev77 = (model77 - data77)
                        if -6.52098520877962e14 <= dev77 <= 6.279259171651384e14:
                            dev77 = 0
                    elif (linecount == 78):
                        model78 = float(row[1])*1e7
                        data78 = 3.057423846387185e15
                        dev78 = (model78 - data78)
                        if -5.31220309387185e14 <= dev78 <= 4.72993814151646e14:
                            dev78 = 0
                    elif (linecount == 80):
                        model80 = float(row[1])*1e7
                        data80 = 3.743016635679190e15
                        dev80 = (model80 - data80)
                        if -8.75247903679190e14 <= dev80 <= 8.370264676223665e14:
                            dev80 = 0
                    elif (linecount == 81):
                        model81 = float(row[1])*1e7
                        data81 = 5.094633778182841e15
                        dev81 = (model81 - data81)
                        if -1.072526096182841e15 <= dev81 <= 9.06880857922256e14:
                            dev81 = 0
                    elif (linecount == 82):
                        model82 = float(row[1])*1e7
                        data82 = 5.143953939483878e15
                        dev82 = (model82 - data82)
                        if -8.97956760483878e14 <= dev82 <= 7.97777845558957e14:
                            dev82 = 0
                    elif (linecount == 84):
                        model84 = float(row[1])*1e7
                        data84 = 6.867838653131900e15
                        dev84 = (model84 - data84)
                        if -1.128324726131900 <= dev84 <= 1.075443694110916e15:
                            dev84 = 0
                    elif (linecount == 85):
                        model85 = float(row[1])*1e7
                        data85 = 8.571451918411129e15
                        dev85 = (model85 - data85)
                        if -1.525168381411129e15 <= dev85 <= 1.328935142710727e15:
                            dev85 = 0
                    elif (linecount == 87):
                        model87 = float(row[1])*1e7
                        data87 = 8.407873955784528e15
                        dev87 = (model87 - data87)
                        if -2.442385366784528e15 <= dev87 <= 2.425970816169801e15:
                            dev87 = 0
                    elif (linecount == 88):
                        model88 = float(row[1])*1e7
                        data88 = 8.407873955784528e15
                        dev88 = (model88 - data88)
                        if -1.474429304784528e15 <= dev88 <= 1.296254164528092e15:
                            dev88 = 0
                    elif (linecount == 90):
                        model90 = float(row[1])*1e7
                        data90 = 1.0697657836811518e16
                        dev90 = (model90 - data90)
                        if -1.910124407811518e15 <= dev90 <= 1.642034959856675e15:
                            dev90 = 0
                    elif (linecount == 91):
                        model91 = float(row[1])*1e7
                        data91 = 1.0392891699183728e16
                        dev91 = (model91 - data91)
                        if -1.892178521183728e15 <= dev91 <= 1.462421817165012e15:
                            dev91 = 0
                    elif (linecount == 93):
                        model93 = float(row[1])*1e7
                        data93 = 1.0697657836811518e16
                        dev93 = (model93 - data93)
                        if -1.971830682811518e15 <= dev93 <= 1.519115653156401e15:
                            dev93 = 0
                    elif (linecount == 95):
                        model95 = float(row[1])*1e7
                        data95 = 1.1893613120972284e16
                        dev95 = (model95 - data95)
                        if -2.165568556972284e15 <= dev95 <= 1.745389369874400e15:
                            dev95 = 0
                    elif (linecount == 96):
                        model96 = float(row[1])*1e7
                        data96 = 1.1779577318070568e16
                        dev96 = (model96 - data96)
                        if -2.109441554070568e15 <= dev96 <= 1.723563060628167e15:
                            dev96 = 0
                    elif (linecount == 98):
                        model98 = float(row[1])*1e7
                        data98 = 1.1334263513008370e16
                        dev98 = (model98 - data98)
                        if -2.196769977008370e15 <= dev98 <= 1.638827899078336e15:
                            dev98 = 0
                    # fixed error of last elif originally being else
                    elif (linecount == 100):
                        model100 = float(row[1])*1e7
                        data100 = 9.529705899902194e15
                        dev100 = (model100 - data100)
                        if -1.791785213902194e15 <= dev100 <= 1.413143814022726e15:
                            dev100 = 0
                # now calculate RMSD from deviations and write 1 value into results file. if i write the deviation values inside the if statements can i still do a calculation with them outside the if statements?
                rmsd = ((dev73**2 + dev75**2 + dev77**2 + dev78**2 + dev80**2 + dev81**2 + dev82**2 + dev84**2 + dev85**2 + dev87**2 + dev88**2 + dev90**2 + dev91**2 + dev93**2 + dev95**2 + dev96**2 + dev98**2 + dev100**2)/16)**0.5
                # i can't tell if the code is obeying the dev = 0 between errors because no set came back as 0..?
                results.write(str(dvald) + "".join(" "*(23 - len(str(dvald)))) +"O2-->2O delta value\n" + str(dvale) + "".join(" "*(23 - len(str(dvale)))) + "O2-->2O* delta value\n" + str(dvalf) + "".join(" "*(23 - len(str(dvalf)))) + "O2-->O2*, O3-->O2+O, and O3-->O3* delta value\n" + str(rmsd) + "".join(" "*(23 - len(str(rmsd)))) +"RMSD" + "\n\n")
                # will it do two enters if i put two \n's like that or does it ignore the second one. also is it ok to put \n at the beginning of a quote like i did
results.close()
print("Done!")
# would like to have it read through the results file and pull the smallest RMSD value and associated parameters but i think i would have to write to a csv in order to do that so i'll leave that part manual for now