# MONACO CODE
written initially by Anton Vasyunin

This code is written in the Fortran 90 language, which means that one can have a line with more than 80 characters in it, and that `MAGIC` isn't real unless explicitly stated. 

## Compiling the code

To compile this code run `make` _**three**_ times. If you encounter any errors, or if the compilation doesn't finish after the third attempt, then it's time to pay attention to the compiler messages. 

## Running the code

Once complete, the code will produce an excecutable called `monaco`, which can be run using

`./monaco`

Note: when running the code, you must call it from the directory that contains the input files, most of which end in ".in" or ".dat"

## Output

The code produces a number of output files, the most relevant of which are contained in the directories, `ab` and `csv`. 

`ab` contains, for each species in the model, a line which shows the time (in yr) and the abundance relative to hydrogen. These files are whitespace separated.

`csv` contains the same information as `ab`, but is formated, as the name suggest, in a more ASCII reader-friendly csv format.

# Above is the readme from the 'main' branch on 11/21

### **Note: This code is written on a Linux Platform with Python 3.8.5 installed.**

# Note on How to Compare the Model and Experimental Data:

1. In the reaction_fitting_factors directory, locate the csv file and in the 4th, 5th, and 6th lines, put arguments to be used for the creation of a numpy linspace (the values contained in it will be the possibilities of the fitting factors for a set of reactions (where each of the 4th, 5th and 6th lines are for a set of reactions)). Also, in the line, put the set of reactions that the fitting factor is modified for. Ensure that the generalized optimization python script (currently baragiola_generalized_serial) contains the paths to these files in the places where they are opened (**Should I specify line numbers in the current generalization file? I don't think so, because they could change**).

2. In the experimental_data directory, place csv files containing (in each line) the time and (CHECK IF THIS IS WHAT THE DATA IS) the abundance alone (or is it as a percentage of something else).

3. Ensure that the original optimization python script (currently baragiola_optimization_test.py) and the generalized optimization python script (currently baragiola_generalized_serial.py) are writing to separate photo_processes and results files (conventionally, the files are numbered for consistency (e.g. baragiola_optimization_test.py writing to photo_processes_2.dat and results_2, and baragiola_generalized_serial writing to photo_processes_3.dat and results_3)). In compareOptimizations.py, ensure the correct file names and paths are referenced for the outputted files made by the original and generalized python scripts (**I think I should be more specific; I have no idea how to do that right now**).

4. In the generalized optimization python script, ensure that the experimental_data file you put in the experimental_data file is the one being read **this seems like it needs to be more specific, but I'm not sure how**.

5. Provide a parameter_inputs_template.dat file with reactions (**I'm not sure why reactions not being processed in baragiola_optimization_test.py are included here**), using the example on the main branch as a guide. As you can see in that file on the main branch, some lines have an integer in the last column be an integer. A line has an integer in that column if it is a reaction for which the fitting factor is being modified. A line having an integer n in this column means that this reaction is part of the group in the (n + 3)th line of the csv file you put in the reaction_fitting_factors directory. This is crucial for the correct function of the generalized optimization python script.

6. Run compareOptimizations.py to compare the runtimes and results of the original and generalized optimization python scripts. This file will print to the command line the following things: these runtimes, as well as any instances where the results and photo_processes files produced by the original and generalized python scripts are not the same, do not contain the same values, or (**is there something else I'm forgetting or should I just be more general?**).
 
##### Is this next one ok? 

7. (Optional--Is it?) Run clean.sh to get rid of some files produced by running monaco and compareOptimizations.py (**how much more specific should we get on which of these files are produced by monaco and which are produced by running compareOptimizations.py).

