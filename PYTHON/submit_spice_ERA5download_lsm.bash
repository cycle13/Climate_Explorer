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

year=1979
end=1979

echo "Running all months between ${year} and ${end} inclusive"

while [ $year -le $end ];
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

    # once sufficient space in the queue
    # submit the next twelve months to run
#    for month in 01 02 03 04 05 06 07 08 09 10 11 12;
##    for month in 08 10;
#    do
        # separate file for each job
        spice_script=spice_era5get_${year}.bash
        
        echo "#!/bin/bash -l" > ${spice_script}
        echo "#SBATCH --mem=30G" >> ${spice_script}
        echo "#SBATCH --ntasks=1" >> ${spice_script}
        echo "#SBATCH --output=/data/users/hadkw/WORKING_HADISDH/slurm_logs/era5_lsm_get_${year}.txt" >> ${spice_script}
        echo "#SBATCH --time=300" >> ${spice_script}
        echo "#SBATCH --qos=normal" >> ${spice_script}

        echo module load scitools/default-current
#        echo "python get_era5_lsm.py --start ${year} --end ${year}" >> ${spice_script}
        echo "python get_era5_lsm.py --start ${year} --end ${year} --remove" >> ${spice_script}
        
        sbatch ${spice_script}       
#    done

    echo "Submitted ${year}"

    let year=${year}+1
    
done

# remove all job files on request
#rm -i ${cwd}/spice_hadisdh_grid_*

#  END
#************************************************************************
