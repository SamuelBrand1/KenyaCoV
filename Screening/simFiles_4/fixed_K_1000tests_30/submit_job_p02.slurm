#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --partition=n1-standard-32
#SBATCH --time=10:00:00
#SBATCH --job-name=job_p02
#SBATCH --mail-user=rabia.aziza@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#BATCH --workdir=/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/
#SBATCH --output=/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_1000tests_30/job_p02.o
#SBATCH --error=/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_1000tests_30/job_p02.o

module purge
module load julia

cd "/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/"
julia /home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/distribute.jl
julia /home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_1000tests_30/p02_SympS.jl
