#!/bin/bash#SBATCH --job-name=PE_equilibration
#SBATCH -o out.o%j
#SBATCH -e err.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
module load intel/18.0.2 impi/18.0.2 lammps/9Jan20
ibrun lmp -in lmp.in
