#!/bin/bash
#SBATCH --job-name=Mono_JS
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-socket=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=160:00:00
#SBATCH --account=jsampath
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
echo "Date   = $(date)"
echo "host   = $(hostname -s)"
echo "Directory = $(pwd)"
ml purge
ml ngc-lammps/10Feb2021
INPUT=in_nve.poly
mpirun -np 1 lmp -in $INPUT -pk kokkos -k on g 1 -sf kk
#mpirun lmp -in $INPUT -pk kokkos neigh half newton off comm no -k on g 2 -sf kk #gpu +KOKKOS
