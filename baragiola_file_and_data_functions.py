import math, time

def getInitialO2():
    return 5.7E22

def getFlux():
    return 2.33e14

def modelCSVFileName():
    return "./csv/bO3.csv"

def serialProgressShellScriptGeneration(numFittingFactorCombinations):
    numFittingFactorCombinationsDivided = math.floor(numFittingFactorCombinations / 100)
    progressShellScriptString = f"""#!/usr/bin/bash

resultsFilepathSerial=$(pwd)
resultsFilepathFormatted="$resultsFilepathSerial/{serialAllOutputFileName()}"

totalNumLines=`wc --lines < "$resultsFilepathFormatted"`

linesPerSimulation=8
approxNumSimulationsDone=$(expr $totalNumLines / $linesPerSimulation)

# The operand on the left is the number of simulations being run
approxNumSimulationsLeft=$(expr {numFittingFactorCombinations} - $approxNumSimulationsDone)

# The operand on the right is for the number of simulations being run divided by 100, if it is not an integer, the nearest integer to it
percentOfSimulationsDone=$(expr $approxNumSimulationsDone / {numFittingFactorCombinationsDivided})

echo "Approximate Number of Simulations Done: $approxNumSimulationsDone of {numFittingFactorCombinations}"
echo "Percent of Simulations Done: $percentOfSimulationsDone" 

startTime={int(time.time())}

echo "Start Time: " 
TZ=America/Chicago date -d @$startTime +"%Y-%m-%d %I:%M %P" """
    
    with open("baragiolaProgress.sh", "w") as progressShellScriptFile:
        progressShellScriptFile.write(progressShellScriptString)

# original is a boolean parameter--supply True for the original parallel script, and False for the Model Nu Trial Nu Same parallel script
def parallelProgressShellScriptGeneration(numFittingFactorCombinations, original):  
    numFittingFactorCombinationsDivided = math.floor(numFittingFactorCombinations / 100)
    parallelScriptOutputFileName = parallelAllOutputForCoreFileName() if original == True else parallelTrialNuIonNuSameAllOutputForCoreFileName()

    progressShellScriptString = f"""#!/usr/bin/bash

resultsFilepathParallel=$(pwd)

core0filepath="$resultsFilepathParallel/baragiola_files_core0/{parallelScriptOutputFileName}"
core1filepath="$resultsFilepathParallel/baragiola_files_core1/{parallelScriptOutputFileName}"
core2filepath="$resultsFilepathParallel/baragiola_files_core2/{parallelScriptOutputFileName}"
core3filepath="$resultsFilepathParallel/baragiola_files_core3/{parallelScriptOutputFileName}"

core0numLines=`wc --lines < "$core0filepath"`
core1numLines=`wc --lines < "$core1filepath"`
core2numLines=`wc --lines < "$core2filepath"`
core3numLines=`wc --lines < "$core3filepath"`

resultsFilepath=$(pwd)
resultsFilepathFormatted="$resultsFilepath/{serialAllOutputFileName()}"

totalNumLines=$(expr $core0numLines + $core1numLines + $core2numLines + $core3numLines)

linesPerSimulation=8
approxNumSimulationsDone=$(expr $totalNumLines / $linesPerSimulation)

# The operand on the left is the number of simulations being run
approxNumSimulationsLeft=$(expr {numFittingFactorCombinations} - $approxNumSimulationsDone)

# The operand on the right is for the number of simulations being run divided by 100, if it is not an integer, the nearest integer to it
percentOfSimulationsDone=$(expr $approxNumSimulationsDone / {numFittingFactorCombinationsDivided})

echo "Approximate Number of Simulations Done: $approxNumSimulationsDone of {numFittingFactorCombinations}"
echo "Percent of Simulations Done: $percentOfSimulationsDone" 

startTime={int(time.time())}

echo "Start Time: " 
TZ=America/Chicago date -d @$startTime +"%Y-%m-%d %I:%M %P" """
    
    with open("baragiolaProgress.sh", "w") as progressShellScriptFile:
        progressShellScriptFile.write(progressShellScriptString)

# For the parallel scripts
def num_processors_to_use():
    return 4

def min_field_width():
    return 30

def serialAllOutputFileName():
    return "resultsSerialAll.txt"

def serialBestResultsOutputFileName():
    return "resultsSerialBest.txt"

def parallelAllOutputForCoreFileName():
    return "resultsParallelAll.txt"

def parallelBestResultsFileName():
    return "resultsParallelBest.txt"

def parallelTrialNuIonNuSameAllOutputForCoreFileName():
    return "resultsParallelTrialNuIonNuSameAll.txt"

def parallelTrialNuIonNuSameBestResultsFileName():
    return "resultsParallelTrialNuIonNuSameBest.txt"

def process_model_data(data):
    return (data / getInitialO2()) * 100

def projectTimeZone():
    return "America/Chicago"