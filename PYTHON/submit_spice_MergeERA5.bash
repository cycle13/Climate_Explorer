#!/bin/bash
# set -x
# ************************************************************************
# Simple script to create and submit SPICE jobs to Slurm.
#
# Author: hadkw
# Date: 16 March 2020
#
# Does each year combination as a separate script.  
#    Polls the queue to not overload the system or quotas (Thanks to JJK for this)
#
# ************************************************************************
# START

#for var in q rh t e td tw dpd;
for var in t e td tw dpd;
do
    # check for number of jobs running at the moment
    n_jobs=`squeue -l | grep hadkw | wc -l`

    # if more than 300 wait and try again in 2 mins
    while [ $n_jobs -gt 15 ];
    do
        echo `date` "  SPICE queue for user hadkw maxed out - sleeping 2mins"
        sleep 2m
        n_jobs=`squeue -l | grep hadkw | wc -l`
    done

    res='pentad'

    # separate file for each job
    spice_script=spice_mergeera5_${var}_${res}.bash
        
    echo "#!/bin/bash -l" > ${spice_script}
    echo "#SBATCH --mem=30G" >> ${spice_script}
    echo "#SBATCH --ntasks=1" >> ${spice_script}
    echo "#SBATCH --output=/data/users/hadkw/WORKING_HADISDH/slurm_logs/mergeera5_${var}_${res}.txt" >> ${spice_script}
    echo "#SBATCH --time=300" >> ${spice_script}
    echo "#SBATCH --qos=normal" >> ${spice_script}

    echo module load scitools/default-current
    echo "python MergeAggRegridERA5.py --var ${var} --freq ${res}" >> ${spice_script}

    res='monthly'

    # separate file for each job
    spice_script=spice_mergeera5_${var}_${res}.bash
        
    echo "#!/bin/bash -l" > ${spice_script}
    echo "#SBATCH --mem=30G" >> ${spice_script}
    echo "#SBATCH --ntasks=1" >> ${spice_script}
    echo "#SBATCH --output=/data/users/hadkw/WORKING_HADISDH/slurm_logs/mergeera5_${var}_${res}.txt" >> ${spice_script}
    echo "#SBATCH --time=300" >> ${spice_script}
    echo "#SBATCH --qos=normal" >> ${spice_script}

    echo module load scitools/default-current
    echo "python MergeAggRegridERA5.py --var ${var} --freq ${res}" >> ${spice_script}
        
    sbatch ${spice_script}       

    echo "Submitted ${var}"

done

# remove all job files on request
#rm -i ${cwd}/spice_hadisdh_grid_*

#  END
#************************************************************************
