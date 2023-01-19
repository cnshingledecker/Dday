#!/usr/bin/bash

core0filepath="/mnt/c/Users/danie/Desktop/Discovery Day/2020-2021 School Year/Actual Project/Dday/baragiola_files_core0/output_file_results"
core1filepath="/mnt/c/Users/danie/Desktop/Discovery Day/2020-2021 School Year/Actual Project/Dday/baragiola_files_core1/output_file_results"
core2filepath="/mnt/c/Users/danie/Desktop/Discovery Day/2020-2021 School Year/Actual Project/Dday/baragiola_files_core2/output_file_results"
core3filepath="/mnt/c/Users/danie/Desktop/Discovery Day/2020-2021 School Year/Actual Project/Dday/baragiola_files_core3/output_file_results"

core0numLines=`wc --lines < $core0filepath`
core1numLines=`wc --lines < $core1filepath`
core2numLines=`wc --lines < $core2filepath`
core3numLines=`wc --lines < $core3filepath`

totalNumLines=`expr $core0numLines + $core1numLines + $core2numLines + $core3numLines`
linesPerSimulation=8
approxNumSimulationsDone=`expr $totalNumLines / $linesPerSimulation`

# The operand on the left is the number of simulations being run
approxNumSimulationsLeft=`expr 42875 - $approxNumSimulationsDone`

# The operand on the right is for the number of simulations being run divided by 100, if it is not an integer, the nearest integer higher than it 
percentOfSimulationsDone=`expr $approxNumSimulationsDone / 429`

beginTime=1674104850 # The beginning time of the simulation in seconds since the Unix epoch
currentTime=$EPOCHSECONDS
timeElapsed=`expr $currentTime - $beginTime`
approxNumSimulationsDoneModified=`expr $approxNumSimulationsDone \* 1000`
millisecondsPerSimulation=`expr $approxNumSimulationsDoneModified / $timeElapsed`
timeForRemainingSimulations=`expr $millisecondsPerSimulation \* $approxNumSimulationsLeft`
timeForRemainingSimulationsModified=`expr $timeForRemainingSimulations / 1000`
finishTime=`expr $currentTime + $timeForRemainingSimulationsModified`
echo "Approximate Number of Simulations Done: $approxNumSimulationsDone"
echo "Percent of Simulations Done: $percentOfSimulationsDone"

# Note: time is not working perfectly yet
echo "Approximate Finish Time: " 
TZ=America/Chicago date -d @${finishTime} +"%Y-%m-%d %I:%M %P"
