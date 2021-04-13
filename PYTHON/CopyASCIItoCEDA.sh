#********************************************
# Script to copy ascii versions to CEDA names
#********************************************

# Configuration file to read in UPDATE INFORMATION WITHIN THIS FILE!!!
source ../../HADISDH_BUILD/F1_HadISDHBuildConfig.txt

echo $StartYear # This should be 1973
echo $EndYear # This should be 2019 or later
echo $Version # This should be v420_2019f or later
echo $VersionDots # This should be 4.2.0.2019f or later

WorkDir="UPDATE"$EndYear
MyDir="/scratch/hadkw/"

for var in q RH T Td e Tw DPD;

do 

       if [ "$var" == "q" ];
       then   
	   cedavar='huss'
	   
       elif [ "$var" == "RH" ];
       then
	   cedavar='hurs'
	   
       elif [ "$var" == "T" ];
       then
	   cedavar='tas'
	   
       elif [ "$var" == "Td" ];
       then
	   cedavar='tds'
	   
       elif [ "$var" == "Tw" ];
       then
	   cedavar='tws'
	   
       elif [ "$var" == "e" ];
       then
	   cedavar='vps'
	   
       elif [ "$var" == "DPD" ];
       then
	   cedavar='dpds'
	   
       fi


    cp $MyDir$WorkDir/STATISTICS/GRIDS/HadISDH.land${var}.${VersionDots}_FLATgridHOM5by5_anoms8110_actuals.dat $MyDir$WorkDir/STATISTICS/GRIDS/${cedavar}-land_HadISDH_HadOBS_${VersionDash}_19730101-${EndYear}1231_actual.dat
    cp $MyDir$WorkDir/STATISTICS/GRIDS/HadISDH.land${var}.${VersionDots}_FLATgridHOM5by5_anoms8110_anomalies.dat $MyDir$WorkDir/STATISTICS/GRIDS/${cedavar}-land_HadISDH_HadOBS_${VersionDash}_19730101-${EndYear}1231_anomaly8110.dat
    cp $MyDir$WorkDir/STATISTICS/GRIDS/HadISDH.land${var}.${VersionDots}_FLATgridHOM5by5_anoms8110_2sig_uncertainty.dat $MyDir$WorkDir/STATISTICS/GRIDS/${cedavar}-land_HadISDH_HadOBS_${VersionDash}_19730101-${EndYear}1231_uncertainty2sig.dat

done 
