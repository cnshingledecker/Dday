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
