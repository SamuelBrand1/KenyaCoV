#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=n1-standard-32
#SBATCH --time=10:00:00
#SBATCH --job-name=job_p08
#BATCH --workdir=/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/
#SBATCH --output=/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening_CTC/simFiles_5/fixed_K_1000tests/job_p08.o
#SBATCH --error=/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening_CTC/simFiles_5/fixed_K_1000tests/job_p08.o

module purge
module load julia

cd "/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/"
julia /home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening_CTC/simFiles_5/fixed_K_1000tests/p08_MSCT_CTC.jl
