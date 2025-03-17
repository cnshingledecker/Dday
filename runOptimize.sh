#!/bin/bash

# Script to run optimize.py with specified parameters

# Define input parameters and their example values
POP_SIZE=100          # Population size (default: 50, increased for more exploration)
GENERATIONS=50        # Number of generations (default: 30, increased for better convergence)
MUTATION_RATE=0.15    # Mutation rate (default: 0.2, slightly reduced for stability)
CROSSOVER_RATE=0.85   # Crossover rate (default: 0.8, slightly increased for diversity)
TOURNAMENT_SIZE=5     # Tournament size (default: 3, increased for better selection pressure)
PROCESSES=4           # Number of processes (default: os.cpu_count(), set to a modest value)
OUTPUT_FILE="optimization_results.csv"  # Output file (default: genetic_algorithm_results.csv, customized name)

# Print the parameters for reference
echo "Running optimize.py with the following parameters:"
echo "-----------------------------------------------"
echo "Population Size: $POP_SIZE"
echo "Generations: $GENERATIONS"
echo "Mutation Rate: $MUTATION_RATE"
echo "Crossover Rate: $CROSSOVER_RATE"
echo "Tournament Size: $TOURNAMENT_SIZE"
echo "Number of Processes: $PROCESSES"
echo "Output File: $OUTPUT_FILE"
echo "-----------------------------------------------"

# Call the Python script with the specified arguments
python3 optimize.py \
    --pop-size "$POP_SIZE" \
    --generations "$GENERATIONS" \
    --mutation-rate "$MUTATION_RATE" \
    --crossover-rate "$CROSSOVER_RATE" \
    --tournament-size "$TOURNAMENT_SIZE" \
    --processes "$PROCESSES" \
    --output "$OUTPUT_FILE"

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "optimize.py completed successfully."
else
    echo "Error: optimize.py failed to run."
fi