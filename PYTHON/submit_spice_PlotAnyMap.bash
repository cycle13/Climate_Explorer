#!/bin/bash
# set -x
# ************************************************************************
# Simple script to create and submit SPICE jobs to Slurm.
#
# Author: hadkw
# Date: 4 Feb 2021
#
#
# ************************************************************************
# START

year=1973
end=2020

echo "Running all year averages between ${year} and ${end} inclusive"

while [ $year -le $end ];
do

    for var in q rh t td tw e dpd;
#    for month in 01 02 03 04 05 06 07 08 09 10 11 12;
##    for month in 08 10;
    do
        # separate file for each job
        spice_script=spice_PlotAnyMap_${year}.bash
        
        echo "#!/bin/bash -l" > ${spice_script}
        echo "#SBATCH --mem=1G" >> ${spice_script}
        echo "#SBATCH --ntasks=1" >> ${spice_script}
        echo "#SBATCH --output=/scratch/hadkw/slurm_logs/PlotAnyMap_${year}_${var}.txt" >> ${spice_script}
        echo "#SBATCH --time=30" >> ${spice_script}
        echo "#SBATCH --qos=normal" >> ${spice_script}

        echo "module load scitools/default-current" >> ${spice_script}
        echo "export MPLBACKEND='Agg'" >> ${spice_script}
	echo "python PlotAnyMap.py --var ${var} --mchoice 0 --mchoiceend 11 --mmult False --ychoice ${year} --ychoiceend ${end} --ymult True" >> ${spice_script}
        echo "unset MPLBACKEND" >> ${spice_script}        
        
        sbatch ${spice_script}       

        echo "Submitted ${var}"

    done

    echo "Submitted ${year}"

    let year=${year}+1
    
done

# remove all job files on request
#rm -i ${cwd}/spice_hadisdh_grid_*

#  END
#************************************************************************
