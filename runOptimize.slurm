#!/bin/bash
#SBATCH -J genetic_opt         # Job name
#SBATCH -o genetic_opt_%j.out  # Output file name (%j expands to jobID)
#SBATCH -e genetic_opt_%j.err  # Error file name (%j expands to jobID)
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 48                  # Number of tasks (should match PROCESSES in runOptimize.sh)
#SBATCH -p skx                 # Queue/partition name
#SBATCH -t 24:00:00            # Run time (hh:mm:ss)
#SBATCH -A TG-PHY240246        # Replace with your project allocation account
#SBATCH --mail-user=shingledeckercn@vmi.edu  # Email notification
#SBATCH --mail-type=ALL        # Type of email notification (BEGIN, END, FAIL, ALL)

# Load required modules (modify as needed for Stampede3)
module load python3

# Display job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"

# Go to the directory where your scripts are located
cd $SLURM_SUBMIT_DIR

# Make sure the runOptimize.sh script is executable
chmod +x runOptimize.sh

# Run the optimization script
./runOptimize.sh

# Display job completion information
echo "Job completed at: $(date)"