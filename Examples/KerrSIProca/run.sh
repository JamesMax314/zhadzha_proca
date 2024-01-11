#!/bin/sh

#SBATCH -J GRDZhadzha #Give it something meaningful.
#SBATCH -o standard_output_file.%J.out
#SBATCH -e standard_error_file.%J.err
#SBATCH -p cosma8 #or some other partition, e.g. cosma, cosma6, etc.
#SBATCH -A dp202 #e.g. dp004
#SBATCH -o sbatch.out
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=james.l.maxwell@durham.ac.uk

source ~/chombo/Cosma-Durham-modules.sh
source ~/chombo/paths.sh 

# Run the program
mpirun -np 4 ./Main_ProcaField3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.COSMA8.Intel2022.ex params.txt 
