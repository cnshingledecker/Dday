# Generation of Models Predicting the Content of Ice in Space Upon Irradiation

***

### Project Background

The general goal of this specific project is, and has been in previous related work, to develop better models to predict the content of certain ices in space upon irradiation. To accomplish this goal, a program called monaco has been used to run models given certain input files and parameters. Later, Baragiola's work was built upon by a student at Wellesley College, Ella Mullikin, in the 2021 paper [A New Method for Simulating Photoprocesses in Astrochemical Models](https://arxiv.org/abs/2101.01209). Specifically, Mullikin's work can be seen in the file baragiola_optimization.py, which has a comment that has been added but is otherwise not modified. This code does a good job of evaulating a certain model using a performance metric known as Root Mean Square Deviation (RMSD). This code was customized for a specific experiment, and only one model simulation (using certain fitting factors for reactions) was ran at a time.

The James Webb Space Telescope is scheduled to be launched in October of 2021, and will be collecting data using a wide variety of instruments. **(check) However, certain experiments, including certain ones that look at the content of ice in space, can only occur while the telescope uses up a limited supply of a substance used for heat that these experiments cannot be conducted without (check). Thus, increasing the number of models that can be developed while this is the case would be incredible, as these models can be used in the JWST mission (check).**

Thus, this current project was intiated to both generalize the previously written baragiola script (so different inputs can be used to develop a wide variety of models) and to decrease the runtime of the code (using High-Performance Computing techniques; namely, the Message Passing Interface (MPI) is used in the parallelization of the previously-written baragiola script).

***

### Project Information (Generalization and Parallelization)
Note that monaco reads from certain input files (including reading fitting factors). Thus, when writing a different baragiola_optimization where input files were different that the original, the monaco script that a new version used had to use the files that were different because of differences from the original input files (such as when certain fitting factors were written to files, and other differences).

#### Development Machine and Environment Specifications
- Machine specifications:
    * Development Machine (Link to specs): [Dell Inspiron 7586](https://www.dell.com/support/manuals/en-us/inspiron-15-7586-2-in-1-laptop/inspiron_7586_setupandspecs/specifications-of-inspiron-7586?guid=guid-7c9f07ce-626e-44ca-be3a-a1fb036413f9&lang=en-us)
        - General Information: Intel i7 8th Generation processor (4 cores, speed: Up to 1.60/1.80 GHz (Boost - 3.90/4.60 GHz), 512 GB SSD, 16 GB RAM)
- Environment Specifications
    * Code ran using the Ubuntu (version 20.04) terminal
- Package/Language specifications
    * Python: Developed on version 3.8.5
    * Packages installed (on Ubuntu; for a full list, see installedPackagesUbuntu.txt)
        - openmpi-bin (possibly 'make all install' ran after running the initial command to install openmpi-bin)
        - pip3 (version 20.0.2)
        - Installed using pip3: mpi4py


#### Description of 2 Main Parts of the Project: (1) Generalization and (2) Parallelization

##### Generalization

How to run this: (how you would run a python script). Example "python3 baragiola_optimization_generalization_rmsd_generalization.py".

As previously described, one of the goals of this project is to write the code in such a way that it can be ran with various input files to produce different models.

Generally, what the original Python script (baragiola_optimization.py) did was as follows:

1. Generate 3 numpy linspaces and access them using 3 nested for loops so on each iteration of the innermost for loop, it is as if you are accessing combinations Cartesian combinations of their values (each combination contains fitting factors that will be used by monaco). **Note that each numpy linspace contains the fitting factors associated with a certain reaction; thus, in each fitting_factor_combination (combination), combination[i] is associated with a certain reaction.**
2. For each fitting factor combination:
    1. Replace certain values in a file containing various reactions and fitting factors (only certain reactions had certain fitting factors replaced); this file is photo_processes.dat.
    2. Run monaco using this photo_processes.dat file (and other files).
    3. From the results of the model (produced by monaco), compare the model predictions (from this iteration of the model) with experimental data, calculate the Root Mean Square Deviation (of the model data from the experimental data), and write to a file the RMSD and the fitting factors that led to it.

Note that in photo_processes.dat is where certain reactions and the fitting factors used are located (the third fitting factor (from left to right), where the rightmost 3 columns are the fitting factors) is the only one being modified.

The fitting factors that produced the least RMSD are thus the ones that produced the most accurate model.



The following changes were made in the generalization the script (see baragiola_optimization_generalization_rmsd_generalization) (note: this list may not be all changes):

1. Place the arguments for the 3 numpy linspaces (and the associated reaction) in a csv file, to be read in when the program runs.
    * This data is in the file 'reaction_fitting_factor_vector_arguments.csv', which is in the directory 'reaction_fitting_factor_linspace_args'. This is done because otherwise, the file is deleted when clean.sh is run.
2. Instead of accessing a combination of fitting factors in the innermost loop of 3 nested for loops, generate a Cartesian combination of the linspaces (after converting them to lists), and run the code using only 1 for loop that uses 1 fitting factor combination per iteration of the for loop. **The rest of the script is run almost the same inside this for loop.**
3. Note that the experimental file that the experimental data file that the baragiola optimization reads from is a csv file in the 'experimental_data' directory. Make sure the script reads from this file (change the file name that the python script reads from).
4. The output of the program (fitting factors, their reactions, and rmsd that they produced, all in a string written to a file) was composed in a different way:
    1. It was written in an aligned way that is easier to read when writing tests to compare the results to the original baragiola optimization.
    2. The reaction names, having been added to an array when being read in from the file 'reaction_fitting_factor_vector_arguments.csv' (referenced above in (1)), are added to this output string with their associated fitting factors using a for loop.
5. Generating the number of data points (used in calculating the Root Mean Square Deviation) by incrementing a counter variable in a for loop.
6. Having it so the deviation was set to 0 if the deviation of the model value from the experimental value be up to 10% of the experimental value (this was a change from the baragiola_optimization where the lower and upper bounds for the deviation were hard-coded).
7. A version of photo_processes.dat was created (photo_processes_2.dat) that if it has an integer in the last column for each row, that integer indicates the reaction for which the third fitting factor (from left to right) is being modified **(note: this integer (n) must correspond to the linspace created using the (n+4)th line of the file 'reaction_fitting_factor_vector_arguments.csv').**
8. A different version of monaco was created that runs using the copies of the input files that are specific to the generalization.
9. Also, certain FORTRAN files had to be changed to read from the files that the generalized python script writes to.

#### Parallelization

How to run this: mpiexec -n (num_processors) (how you would run any other python file). Example: "mpiexec -n 4 python3 baragiola_generalized_parallel.py". Make sure the value you use for num_processors in the command line has the same value as the same-named variable in baragiola_generalized_parallel.py, and if you run compareOptimizations.py, which runs the original, generalized, and parallel/generalized scripts, you must change num_processors in the command line argument that runs the generalized/parallel baragiola script. 

If you want to run this with different input and output files, you need to overwrite the contents of the following files:
1. photo_processes.dat, with the same conditions specified in (7) under the 'Generalization' header above.
2. A new experimental data csv file, put in the directory 'experimental_data' (note that the python script must reference this file).
You must also change the numpy linspacae arguments in the file 'reaction_fitting_factor_vector_arguments.csv' in the 'reaction_fitting_factor_linspace_args' directory, and change the reactions in each of those lines.

Why write a parallel version of the generalized script? The parallelized version has one significant benefit: multiple models can be ran at once, which can decrease runtime. However, the parallelized script has some drawbacks and possible risks:

Drawbacks:
- There is an overhead time associated with dividing up the fitting factor combinations, sending them to the different processors, and each processor sending the lowest RMSD (and the fitting factors that produced it) back to the root node. For example, if the number of models to be run (the number of fitting factor combinations) is low enough, the improvement in runtime between the original version and the parallel version can be small, nonexistent, or the runtime could even increase.

Possible Risks:
- Loss of precision in sending and receiving data, although significant examples of this have not been observed.

Benefits:
- Possible Decreased Overall Runtime
- Each processor can find the lowest RMSD that was produced from the fitting factors it was sent and send back to the root processor that RMSD and the fitting factors that produced it. This simplifies the process of finding the fitting factors that produce the least RMSD.

Below is an overview of the changes from the parallelized version of the script to the generalized version: 
1. For each processor, a directory is created, and a file is created in it where the fitting factors and RMSD that they produced is written for each fitting factor combination. The files that monaco uses (and many others; don't copy files if they don't need to be (for monaco)).
2. Once the fitting factor combinations are created (a list of them), that list is split into n chunks (where n is the number of processors used).
3. Once step 2 is completed, each chunk (there is 1 for each processor) is split up into mini-chunks, starting with size 15 (the last mini-chunk created is of size chunk-size mod 15). This function was written because on the machine this code was developed on, the processors did not receive chunks of size 16 or greater. The mini-chunk size is specified in the function call, and the for loops that send the fitting factor combinations to each processor are set up to send these mini-chunks (the for loops will send all the mini-chunks, regardless of the mini-chunk size).
3. Steps 8 and 9 from the generalization additions above are not in place for this version, because for each processor, the files used with the original monaco are copied to the directory for each processor. Therefore, there is no need to change the files that monaco reads from.
4. Instead of all of the processors writing to 1 common results file, each processor writes to its own results file.
5. Each processor finds the lowest RMSD and the fitting factors that produced it, and sends that back to the root processor (processor 0). The root processor finds the lowest RMSD (and fittign factors that produced it) from among those returned and writes the one with the lowest RMSD and writes it to a results file.
6. **Note: the [OpenMPI library](https://www.open-mpi.org/) and the python package [mpi4py](https://mpi4py.readthedocs.io/en/stable/index.html) (version 3.0.3) were used in this parallelization. Click on the links to learn more about them.**  

#### Comparison Script

The python script compareOptimizations.py runs the original, generalized and (parallel and generalized) versions of the baragiola script, one after another, in that order, prints out the runtimes for each, and compares the results files for the original and generalized versions. Note that this will only work if all 3 scripts have (a) the same numpy linspace arguments and thus are dealing with the same reactions, (b) the same experimental data files, and (c) the same photo_processes.dat files. 