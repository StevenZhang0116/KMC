#!/bin/bash
#SBATCH --job-name=my_cpp_job       # Job name
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1        # Number of tasks (processes) per node
#SBATCH --cpus-per-task=1          # Number of CPU cores per task (process)
#SBATCH --mem=128G                  # Memory per node (e.g., 4 GB)
#SBATCH --time=30:00:00           # Time limit (hh:mm:ss)
#SBATCH --output=output.log       # Output log file
#SBATCH --error=error.log         # Error log file

# Load necessary modules (if needed)
module load cmake
module load gcc
module -q purge
module -q load openmpi
module list
module load python

lscpu
nvidia-smi

# Run your C++ executable
./test_main
