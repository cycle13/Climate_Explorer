#!/bin/bash -l
#SBATCH --mem=30G
#SBATCH --ntasks=1
#SBATCH --output=/data/users/hadkw/WORKING_HADISDH/slurm_logs/MergeERArh_pentad.txt
#SBATCH --time=300
#SBATCH --qos=normal
module load scitools/default-current
python MergeAggRegridERA5.py --var rh --freq pentad
