#!/bin/bash
# set -x
# ************************************************************************
# Simple script to create and submit SPICE jobs to Slurm.
#
# Author: hadkw
# Date:3 April 2020
#
# Does each year combination as a separate script.  
#    Polls the queue to not overload the system or quotas (Thanks to JJK for this)
#
# ************************************************************************
# START

year=1973
end=2020

echo "Running trends between ${year} and ${end} inclusive"

#while [ $year -le $end ];
for var in q rh t td e tw dpd;
#for var in q rh t;
#for var in td;
#for var in t;
do

	typee='LAND'
	##typee='MARINE'
	#typee='MARINESHIP'
	##typee='BLEND'
	#typee='BLENDSHIP'
	#typee='ERA5'
	#typee='ERA5LAND'
	#typee='ERA5LANDMASK'
	#typee='ERA5MARINE'
	#typee='ERA5MARINEMASK'

#	typee='IDPHA'
#	if [ $var = 'dpd' ]
#	then   
#	    typee='PHA'
#	    
#	fi
#	
#	if [ $var = 'td' ]
#	then
#	    typee='PHADPD'
#	    
#	fi

        echo $var $typee
	
	# separate file for each job
        spice_script=spice_PlotTrendMap_${var}.bash
        
        echo "#!/bin/bash -l" > ${spice_script}
        echo "#SBATCH --mem=5G" >> ${spice_script}
        echo "#SBATCH --ntasks=1" >> ${spice_script}
        echo "#SBATCH --output=/scratch/hadkw/slurm_logs/PlotTrendMap_${year}${end}_${var}_${typee}.txt" >> ${spice_script}
        echo "#SBATCH --time=20" >> ${spice_script}
        echo "#SBATCH --qos=normal" >> ${spice_script}

        echo "module load scitools/default-current" >> ${spice_script}
        echo "export MPLBACKEND='Agg'" >> ${spice_script}
        echo "python PlotTrendMap.py --var ${var} --typee ${typee} --year1 ${year} --year2 ${end}" >> ${spice_script}
        echo "unset MPLBACKEND" >> ${spice_script}        
        
        sbatch ${spice_script}       

    echo "Submitted ${var}"

done

# remove all job files on request
#rm -i ${cwd}/spice_hadisdh_grid_*

#  END
#************************************************************************
