#!/usr/bin/bash
#SBATCH --job-name Ar_conv_ip_1.000
#SBATCH --output Ar_conv_ip_1.000.o%j
#SBATCH --partition=jila
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --nice=1.0
#SBATCH --mem=16G

source ~/.bashrc

REPO_DIR="/data/becker/begh0305/Research/bspline_tdse"
hostname
pwd

mpirun -np $SLURM_NTASKS $REPO_DIR/bspline_tdse.out tise >> results.log
mpirun -np $SLURM_NTASKS $REPO_DIR/bspline_tdse.out >> results.log
