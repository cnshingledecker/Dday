# Generation of Models Predicting the Content of Ice Upon Irradiation
***

### General Project Buildup
- What is the goal of this? To use monaco (a program that gives model predictions using certain input files) to come up with good models to predict the content of interstellar ice upon irradiation (over time)
- Work done by Baragiola
- Work done by Ella Mullikin (Wellesley College) and others on the paper (cite paper)
- Note how the previous work was excellent, but in order to be able to more easily come up with different and more complex models and to do so in less time (especially for the JWST), certain changes had to be made. **Intro project goals of generalization and speed-up of code using parallelization techniques (HPC techniques).**

### This project
- Hardware, environments, and packages/languages/software used in project
- Restate goals of generalizations and parallelization.
- Writing a generalized version of this script.
    * List specific changes (generalizations) that were made 
    * Note that generalization now includes reading from input files; lsit files, explain what they do
    * Note that certain fortran files were adapted to deal with different files (explain better) due to the changing of input files (specially photo_processes_2.dat)
    * Note writing of version of the original script that wrote its output the same as the generalized version, which allowed comparison of the results of both files (compareOptimizations.py)
- Writing a parallel version of the generalized version of the script
    * Why write a parallel version? benefits and drawbacks; mention overhead
    * Give overview of changes
    * Explain what OpenMPI is and why it was used
    * Mention use of mpi4py
    * Describe why certain functions were written and chnages werre made when writing this version.
- Time comparisions between the different scripts

# README file from main branch
Include contents of the README file from the main branch for reference.