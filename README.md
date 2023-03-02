![Example Results Plot](/good_trial_runs/num8_plot.jpg "Example Results Plot")

# Generation of Models Predicting the Content of Ice in Space Upon Irradiation

***

### Table of Contents
1. [Project Introduction](#project-introduction-and-motivation)
2. [Technology and Software Used](#technology-and-software-used)
3. [Project Information and Features](#project-information-and-features)
4. [Code Examples](#code-examples)
5. [Installation](#installation)
6. [Making Changes](#making-changes)
7. [Tests](#tests)
8. [Credits](#credits)


### Project Introduction and Motivation

The general goal of this specific project is, and has been in previous related work, to develop better models to predict the content of certain ices in space upon irradiation. To accomplish this goal, a program called monaco has been used to run models given certain input files and parameters. Later, Baragiola's work was built upon by a student at Wellesley College, Ella Mullikin, in the 2021 paper [A New Method for Simulating Photoprocesses in Astrochemical Models](https://arxiv.org/abs/2101.01209). Specifically, Mullikin's work can be seen in the file baragiola_optimization.py, which has a comment that has been added but is otherwise not modified. This code does a good job of evaulating a certain model using a performance metric known as Root Mean Square Deviation (RMSD). This code was customized for a specific experiment, and only one model simulation (using certain fitting factors for reactions) was ran at a time.

The James Webb Space Telescope was launched in December of 2021, and is collecting data using a wide variety of instruments. However, certain experiments, including certain ones that look at the content of ice in space, can only occur while the telescope uses up a limited supply of a substance used for heat that these experiments cannot be conducted without. Thus, increasing the number of models that can be developed while this is the case would be incredible, as these models can be used in the JWST mission.

Thus, this current project was intiated to both generalize the previously written baragiola script (so different inputs can be used to develop a wide variety of models) and to decrease the runtime of the code (using High-Performance Computing techniques; namely, the Message Passing Interface (MPI) is used in the parallelization of the previously-written baragiola script). Specifically, [Open MPI](https://www.open-mpi.org/), an open-source MPI implementation, is used. 

### Technology and Software Used
- Machine specifications:
    * Development Machine (Link to specs): [Dell Inspiron 7586](https://www.dell.com/support/manuals/en-us/inspiron-15-7586-2-in-1-laptop/inspiron_7586_setupandspecs/specifications-of-inspiron-7586?guid=guid-7c9f07ce-626e-44ca-be3a-a1fb036413f9&lang=en-us)
        - General Information: Intel i7 8th Generation core (4 cores, speed: Up to 1.60/1.80 GHz (Boost - 3.90/4.60 GHz), 512 GB SSD, 16 GB RAM)
- Environment Specifications
    * Code ran using the Ubuntu (version 20.04) terminal (Ubuntu 20.04.1 LTS)
    * **Note: The code (for this project) has developed to run in an environment that accepts Linux shell commands.**
- Package/Language specifications
    * Python: Developed on version 3.8.10
    * Packages installed:
        - On Ubuntu
            - openmpi-bin/focal,now 4.0.3-0ubuntu1 amd64 [installed] 
        - Python:
            - mpi4py (version 3.1.4)
            - numpy (version 1.24.2)

Any used Python packages not installed should be put in this section, as well as in [`requirements.txt`](/requirements.txt).

### Project Information and Features
Note that monaco reads from certain input files (including reading fitting factors). Thus, when writing a different baragiola_optimization where input files were different that the original, the monaco script that a new version used had to use the files that were different because of differences from the original input files (such as when certain fitting factors were written to files, and other differences).

There are three versions of this generalized script:
1. Serial, where only one model can be run at a time.
2. Parallel, where up to the number of cores on the machine can be running at a time.
3. Parallel where `TRIAL_NU` and `ION_NU`, both `model.inp` parameters, are forced to be the same value.


**Note: the README on the 'main' branch may be helpful to read.**

#### Description of 2 Main Parts of the Project: (1) Generalization and (2) Parallelization

#### Note for Parallel and Generalized Scripts
Please note that at the end of the successful model runs that it was given to run, both the parallel and serial scripts run the Python file `plotting.py`, which generates a plot abundance vs. fluence using the best fit parameters the Python script found (with the least error). This happens after a run of monaco with these best fit parameters, which happens after these best fit parameters are written to the `photo_processes` file that monaco reads from. If youdon't want plots generated, simply comment out the lines inserting the best fit parameters into `photo_processes.dat` and the lines calling the two aforementioend Python scripts. These lines are at the bottom of both the serial and parallel Python scripts.

##### Generalization

How to run this: (how you would run a python script). Example `python3 baragiola_generalized_serial.py`.

As previously described, one of the goals of this project is to write the code in such a way that it can be ran with various input files to produce different models.

Generally, what the original Python script (baragiola_optimization.py) did was as follows:

1. Generate 3 numpy linspaces and access them using 3 nested for loops so on each iteration of the innermost for loop, it is as if you are accessing combinations Cartesian combinations of their values (each combination contains fitting factors that will be used by monaco). **Note that each numpy linspace contains the fitting factors associated with a certain reaction; thus, in each fitting_factor_combination (combination), combination[i] is associated with a certain reaction.**
2. For each fitting factor combination:
    1. Replace certain values in a file containing various reactions and fitting factors (only certain reactions had certain fitting factors replaced); this file is photo_processes.dat.
    2. Run monaco using this photo_processes.dat file (and other files).
    3. From the results of the model (produced by monaco), compare the model predictions (from this iteration of the model) with experimental data, calculate the Root Mean Square Deviation (RMSD) of the model data from the experimental data, and write to a file the RMSD and the fitting factors that led to it.

Note that in photo_processes.dat is where certain reactions and the fitting factors used are located (the third fitting factor (from left to right), where the rightmost 3 columns are the fitting factors) is the only one being modified.

The fitting factors that produced the least RMSD are thus the ones that produced the most accurate model.



The following changes were made in the generalization of the script (see baragiola_generalized_serial) (note: this list may not be all changes):

1. Place the arguments for the numpy linspaces (and the associated reaction) in a csv file, to be read in when the program runs.
    1. This data is in the file `reaction_fitting_factor_vector_arguments.csv`, which is in the directory `reaction_fitting_factor_linspace_args`. This is done because otherwise, the file is deleted when clean.sh is run. The format for this is in that file; do this for each of the reaction rates for which you want to try different values.
2. Place the arguments for creating numpy linspaces for `model.inp` values that are desired to be changed in the file [`model_inp_values/model_inp_values.csv`](/model_inp_values/model_inp_values.csv). As with the above file, the format for this is in that file; do this for each of the `model.inp` values for which you want to try different values.
2. Instead of accessing a combination of fitting factors in the innermost loop of 3 nested for loops, generate a Cartesian combination of the linspaces (after converting them to lists), and run the code using only 1 for loop that uses 1 fitting factor combination per iteration of the for loop. **The rest of the script is run the same inside this for loop.**
3. Note that the experimental file that the experimental data file that the baragiola optimization reads from is in a custom function in `exportable_custom_functions.py`.
4. The output of the program (fitting factors, their reactions, the model.inp values, and RMSD that they produced, all in a string written to a file) was composed in a different way:
    1. It was written in an aligned way that is easier to read when writing tests to compare the results to the original baragiola optimization.
    2. The reaction names, having been added to an array when being read in from the file `reaction_fitting_factor_vector_arguments.csv` (referenced above in (1)), are added to this output string with their associated fitting factors using a for loop.
5. Generating the number of data points (used in calculating the Root Mean Square Deviation) by incrementing a counter variable in a for loop.
6. A version of photo_processes.dat was created (`photo_processes_2.dat`) that if it has an integer in the last column for each row, that integer indicates the reaction for which the third fitting factor (from left to right) is being modified **(note: this integer (n) must correspond to the linspace created using the (n+4)th line of the file 'reaction_fitting_factor_vector_arguments.csv').**
7. A different version of monaco was created that runs using the copies of the input files that are specific to the generalization.
8. Also, certain FORTRAN files had to be changed to read from the files that the generalized python script writes to.
9. Wrote lowest RMSD (and the fitting factors and `model.inp` values that produced it) to a results file.

#### Parallelization

How to run this: `mpiexec -n (num_cores)` followed by how you would run any other Python file. An example is as follows: `mpiexec -n 4 python3 baragiola_generalized_parallel.py`. Make sure the value you use for num_cores in the command line has the same value as the same-named variable in `baragiola_generalized_parallel.py`. 

Why write a parallel version of the generalized script? The parallelized version has one significant benefit: multiple models can be ran at once, which can decrease runtime. However, the parallelized script has some drawbacks and possible risks:

Drawbacks:
- There is an overhead time associated with dividing up the fitting factor combinations, sending them to the different cores, and each core sending the lowest RMSD (and the fitting factors that produced it) back to the root node. For example, if the number of models to be run (the number of fitting factor combinations) is low enough, the improvement in runtime between the original version and the parallel version can be small, nonexistent, or the runtime could even increase.

Possible Risks:
- Loss of precision in sending and receiving data, although significant examples of this have not been observed.

Benefits:
- Possible Decreased Overall Runtime
- Each core can find the lowest RMSD that was produced from the fitting factors it was sent and send back to the root core that RMSD and the fitting factors that produced it. This simplifies the process of finding the fitting factors that produce the least RMSD.

Below is an overview of the changes from the parallelized version of the script to the generalized version: 
1. For each core, a directory is created, and a file is created in it where the fitting factors and RMSD that they produced is written for each fitting factor combination. The files that monaco uses (and many others; don't copy files if they don't need to be (for monaco)).
2. Once the fitting factor combinations are created (a list of them), that list is split into n chunks (where n is the number of cores used).
3. Once step 2 is completed, each chunk (there is 1 for each core) is split up into mini-chunks, starting with size 15 (the last mini-chunk created is of size chunk-size mod 15). This function was written because on the machine this code was developed on, the cores did not receive chunks of size 16 or greater. The mini-chunk size is specified in the function call, and the for loops that send the fitting factor combinations to each core are set up to send these mini-chunks (the for loops will send all the mini-chunks, regardless of the mini-chunk size).
3. Steps 8 and 9 from the generalization additions above are not in place for this version, because for each core, the files used with the original monaco are copied to the directory for each core. Therefore, there is no need to change the files that monaco reads from.
4. Instead of all of the cores writing to 1 common results file, each core writes to its own results file.
5. Each core finds the lowest RMSD and the fitting factors that produced it, and sends that back to the root core (core 0). The root core finds the lowest RMSD (and fitting factors that produced it) from among those returned and writes the one with the lowest RMSD and writes it to a results file.
6. **Note: the [OpenMPI library](https://www.open-mpi.org/) and the python package [mpi4py](https://mpi4py.readthedocs.io/en/stable/index.html) (version 3.0.3) were used in this parallelization. Click on the links to learn more about them.**  

The following was done for both the generalized and parallel scripts: a custom function (in exportable_custom_functions.py) was written to allow the user to change values of variables in model.inp. If you wish to reset the values of `model.inp` for any version of this script, set the variable `reset_modelInp` to `True`. If not, set it to `False`.

### Code Examples
The Baragiola scripts are as follows:
    - [Original -- `baragiola_optmization.py`](/baragiola_optimization.py)
    - [Generalized Serial -- `baragiola_generalized_serial.py`](/baragiola_generalized_serial.py)
    - [Generalized Parallel - `baragiola_generalized_parallel.py`](/baragiola_generalized_parallel.py)
    - [Generalized Parallel with Keeping Trial Nu and Ion Nu the same - `baragiola_generalized_parallel_trial_nu_ion_nu_same.py`](/baragiola_generalized_parallel_trial_nu_ion_nu_same.py)

    The parallel scripts contain examples of running models in parallel, and all scripts contain examples of running the models and calculating the RMSD for a model run.

### Installation
To install the software for this project, make sure to clone the code off of GitHub into a Linux environment (Ubuntu was used for the development of this software). Make sure Python 3 is installed (Python 3.8.10 was used to develop this project, so at least that is safe). Make sure to install the packages listed in the [Technology and Software Used](#technology-and-software-used) section of this documentation, which for the Python ones can be done with the [`requirements.txt`](/requirements.txt) file.

### Making Changes

When switching the network, make sure to change the files mentioned below as follows. You may need to change them as necessary throughout your work. 
- [`baragiola_file_and_data_functions.py`](/baragiola_file_and_data_functions.py)
    - Make sure to change the names of the output files for the different versions of the generalized script as desired.
    - Set the minimum field width (when writing `model.inp` values, fitting factors, and the `RMSD` to a file) in the body of that function.
    - Rewrite the definition of the `process_model_data` function so it formats the string data outputted by the model (the dependent variable) so it is the same format as the experimental data, which is in the `exportable_custom_functions.py` file.
- [`model_inp_values/model_inp_values.csv`](/model_inp_values/model_inp_values.csv)
    - If you wish to modify values in `model.inp` to optimize them (by trying out ranges of values, similar to modifying the delta values for the reaction rates), edit the file `model_inp_values/model_inp_values.csv` according to the rules in there, and make sure to set the variable `to_modify_modelInp_values` to `True` if you want to test out the different `model.inp` values specified in the file. If you don't want the contents of that file to be read (whether it contains data for modifying lines or not), set `to_modify_modelInp_values` to `False`. 
- [`reaction_fitting_factor_linspace_args/reaction_fitting_factor_linspace_args.csv`](/reaction_fitting_factor_linspace_args/reaction_fitting_factor_linspace_args.csv)
    - Make sure to put in one line per reaction rate for which you want to try a range of values. The rules for how to format this are in that file.
- [`exportable_custom_functions.py`](/exportable_custom_functions.py)
    - You must change the `setup_experimental_data` function--put your experimental data in there; keep the names `expX` (independent variable) and `expY` (dependent variable).
- [`model.inp`](/model.inp)
    - Make sure to set the initial values as desired.
- [`plotting.py`](/plotting.py)
    - Make sure to change the `w_ions_version` and `w_ions_modified_version` declaration lines to reference the correct output file for the model data. Also, make sure to change any other lines that reference the output file species name, and don't forget to change the axis titles, axis ranges, and the plot title as necessary.

### Tests
To test this code, you should make sure that the outputted CSV files are generating the data you think they should, the models are running with minimal errors, and (if you change it) that the code to modify `model.inp` and fitting factor values is working as expected.

### Credits
Anton Vasyunin's development of `monaco` and Christopher Shingledecker's work on it are invaluable to this project. Also, Ella Mullikin's writing of the original Baragiola script [`baragiola_optimization.py`](/baragiola_optimization.py) was invaluable to the development of this script. Finally, the advising of especially Dr. Christopher Shingledecker, as well as Dr. Juan Carlos Araque, have helped Daniel Lopez-Sanders work on this project.