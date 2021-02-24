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

# Above is the readme from the 'main' branch on 11/21/2020

Notes for running the parallel script (baragiola_generalized_parallel.py):

0. Note that this application was developed with the Windows Subsystem for Linux (WSL) using Ubuntu, version 20.04.1 LTS. Also, using Ubuntu, the following command line arguments were used to install packages: 
- pip3 install mpi4py
- sudo apt install openmpi-bin
If the above doesn't work, reach out for support.

Also, this code works for sure with Python version 3.8.5.

1. The original (non-parallel, non-generalized) version of the script is baragiola_optimization.py. There is a generalized version of the script (baragiola_optimization_generalization_rmsd_generalization.py) which produces the same results as the original version of the script. The parallel version of the script is a base MPI program that had the generalized script integrated into it (the parallel version produces the same results using the same code as the generalized, non-parallel version).

2. The original version of the script is built to run with certain files, written-in deviation values, etc. The generalized version of the script (and thus the parallel version as well) are built to handle a wide variety of input files, with the following files needing to be changed (that the baragiola script directly interacts with) if you desire to have the program generate different results:

- (File) photo_processes_2.dat
- (File in Subdirectory) experimental_data: put the experimental csv file in this directory
- (File) reaction_fitting_factor_linspace_args/reaction_fitting_factor_vector_arguments.csv
- (File) Make the csv_model_data file in the parallel baragiola script reference the correct model csv file.

3. Once you are sure the program is reading from and writing to the correct files, you can run the model. The file baragiola_generalized_parallel.py can be ran using the command "mpiexec -n [num_procesors_used] [command to run the file]"; an example is "mpiexec -n 4 python3 baragiola_generalized_parallel.py". Be sure to change the num_processors variable in baragiola_generalized_parallel.py to match num_processors_used.