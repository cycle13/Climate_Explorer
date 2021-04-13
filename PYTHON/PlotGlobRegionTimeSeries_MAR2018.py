#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 22 April 2014
# Last update: 5 May 2019
# Location: /data/users/hadkw/WORKING_HADISDH/UPDATE2018/PROGS/PYTHON/	# this will probably change
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Plots regional average timeseries for all variables with uncertainty estimates
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Inbuilt:
# import matplotlib.pyplot as plt
# import numpy as np
# import numpy.ma as ma
# import sys, os
# import scipy.stats
# import struct
# import cartopy.crs as ccrs
# import cartopy.feature as cpf
# import datetime as dt
# from matplotlib.dates import date2num,num2date
# from scipy.io import netcdf
# import matplotlib.colors as mc
# import matplotlib.cm as mpl_cm
# import pdb
#
# Other:
# from RandomsRanges import LetterRange
# from LinearTrends import MedianPairwise
# from LinearTrends import OLS_AR1Corr
# 
# -----------------------
# DATA
# -----------------------
# directory for timeseries with uncertainties - these will have been created by PROGS/HADISDH_BUILD/hadisdh_error_calculations.py:
# /data/users/hadkw/WORKING_HADISDH/UPDATE2018/STATISTICS/TIMESERIES/:
# 	HadISDH.<domain><var>.<version>_<region>_ts_<timeres>_anoms8110_<mon><year>.dat
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# Set up the editables as desired: land or marine, region, version, thenmon, thenyear, nowmon, nowyear, start year, end year
# >module load scitools/default-current
# >python PlotGlobRegionTimeSeries_MAR2018.py
# 
# -----------------------
# OUTPUT
# -----------------------
# directory for output images:
# /data/users/hadkw/WORKING_HADISDH/UPDATE2018/IMAGES/TIMESERIES/
# 		Plot<domain><region>TimeSeries.<version>_<timeres>_anoms8110_<nowmon><nowyear>.png/eps
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
#
# Version 2 (4 May 2020)
# ---------
#  
# Enhancements
# Can now work with OLS_AR1Corr
#  
# Changes
#  
# Bug fixes
#
# Version 1 (5 May 2019)
# ---------
#  
# Enhancements
#  
# Changes
#  
# Bug fixes
#
#  
# -----------------------
# OTHER INFORMATION
# -----------------------
#
#************************************************************************
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import scipy.stats
import struct
import os.path
import math
#from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
#from netCDF4 import Dataset
from scipy.io import netcdf
from scipy.stats.stats import pearsonr   
from RandomsRanges import LetterRange
from LinearTrends import MedianPairwise
from LinearTrends import OLS_AR1Corr
import pdb

# Set up initial run choices
TrendChoice = 'OLS' # OLS or MPW
RegChoice = 'Trop' # 'Glob','Nhem','Trop','Shem' for global, n hemi, trop, s hemi
timetype='monthly'	#'monthly', 'annual'
nparams=7
param=list(['t','tw','td','q','e','rh','dpd'])	# tw, q, e, rh, t, td, dpd
param2=list(['T','Tw','Td','q','e','RH','DPD'])	# Tw, q, e, RH, T, Td, DPD
param3=list(['T','T$_{w}$','T$_{d}$','q','e','RH','DPD'])	# Tw, q, e, RH, T, Td, DPD
unitees=list(['$^{o}$C','$^{o}$C','$^{o}$C','g kg$^{-1}$','hPa','%rh','$^{o}$C'])
#homogtype=list(['IDPHA','IDPHA','PHADPD','IDPHA','IDPHA','IDPHA','PHA'])	# 'IDPHA','PHA','PHADPD'

Domain = 'land' # land or marine, marineSHIP, blend, blendSHIP
#Domain = 'marineSHIP' # land or marine, marineSHIP, blend, blendSHIP
#Domain = 'blendSHIP' # land or marine, marineSHIP, blend, blendSHIP

lversion='4.3.0.2020f'
mversion='1.1.0.2020f'
bversion='1.1.0.2020f'

styr=1973
edyr=2020
nyrs=(edyr-styr)+1
nmons=(nyrs)*12
climst=1981
climed=2010
stcl=climst-styr
edcl=climed-styr
if (timetype == 'annual'):
    ntims = nyrs
elif (timetype == 'monthly'):
    ntims = nmons

# Set up directories and files

PLOTDIR='/scratch/hadkw/UPDATE'+str(edyr)+'/IMAGES/TIMESERIES/'
DATADIR='/scratch/hadkw/UPDATE'+str(edyr)+'/STATISTICS/TIMESERIES/'
#PLOTDIR='/data/users/hadkw/WORKING_HADISDH/UPDATE'+str(edyr)+'/IMAGES/TIMESERIES/'
#DATADIR='/data/users/hadkw/WORKING_HADISDH/UPDATE'+str(edyr)+'/STATISTICS/TIMESERIES/'


#edyr=2018
#nyrs=(edyr-styr)+1
#nmons=(nyrs)*12
#if (timetype == 'annual'):
#    ntims = nyrs
#elif (timetype == 'monthly'):
#    ntims = nmons

IfType='.dat'	#'.nc'
if (Domain == 'land'):
    version = lversion
    INHFILEST='HadISDH.land'
elif (Domain == 'marine') | (Domain == 'marineSHIP'):
    version = mversion
    INHFILEST='HadISDH.marine'    
    if (Domain == 'marineSHIP'):
        version = version+'SHIP'    
elif (Domain == 'blend') | (Domain == 'blendSHIP'):
    version = bversion
    INHFILEST='HadISDH.blend'    
    if (Domain == 'blendSHIP'):
        version = version+'SHIP'    

#INHFILEED='5by5_'+thenmon+thenyear+'_areaTS_19732013'
if (RegChoice == 'Glob'):
    if timetype == 'monthly':
        INHFILEED='.'+version+'_global_ts_monthly_anoms8110.dat'
    else:
        INHFILEED='.'+version+'_global_ts_annual_anoms8110.dat'
elif (RegChoice == 'Nhem'):
    if timetype == 'monthly':
        INHFILEED='.'+version+'_northern_hemisphere_ts_monthly_anoms8110.dat'
    else:
        INHFILEED='.'+version+'_northern_hemisphere_ts_annual_anoms8110.dat'
elif (RegChoice == 'Shem'):
    if timetype == 'monthly':
        INHFILEED='.'+version+'_southern_hemisphere_ts_monthly_anoms8110.dat'
    else:
        INHFILEED='.'+version+'_southern_hemisphere_ts_annual_anoms8110.dat'
elif (RegChoice == 'Trop'):
    if timetype == 'monthly':
        INHFILEED='.'+version+'_tropics_ts_monthly_anoms8110.dat'
    else:
        INHFILEED='.'+version+'_tropics_ts_annual_anoms8110.dat'
OUTPLOT='Plot'+Domain+RegChoice+'TimeSeries.'+version+'_'+timetype+'_anoms8110_'+TrendChoice

# Set up variables
mdi = -1e30

varH=[]	# nvars(rows) by 4 regions, by mons masked array
uncsHtot=[] # nvars(rows) by 4 regions, by mons masked array
uncsHcov=[] # nvars(rows) by 4 regions, by mons masked array
uncsHsamp=[] # nvars(rows) by 4 regions, by mons masked array
uncsHstat=[] # nvars(rows) by 4 regions, by mons masked array


#************************************************************************
# Subroutines
#************************************************************************
# READDATA
def ReadData(FileName,typee,delimee,skipee):
    ''' Use numpy genfromtxt reading to read in all rows from a complex array '''
    ''' Need to specify format as it is complex '''
    ''' outputs an array of tuples that in turn need to be subscripted by their names defaults f0...f8 '''
    return np.genfromtxt(FileName, dtype=typee,delimiter=delimee,skip_footer=skipee) # ReadData

#************************************************************************
# PlotNiceTimeSeries
def PlotNiceTimeSeries(TheFile,TheHvars,TheHuncsC,TheHuncsSp,TheHuncsSt,
                       TheUnitees,TheMCount,TheYCount,TheTimeType,
		       TheStYr,TheEdYr,TheMDI,TheParams,TheReg,TheTrendChoice):
    ''' Plot a panel for each element of TheHvars '''
    ''' Add Coverage, Sampling and Station uncertainty ranges '''
    ''' Add lines for any extra estimates (TheVars) and HadISDH MASKED versions '''
    ''' Save as png and eps '''
    ''' TheHvars is a multi-row array: rows for vars, columns for months '''
    ''' Ditto TheHuncs C=coverage, Sp=sampling, St=station '''
    ''' TheUnitees is the units name list '''
    ''' TheMCount - number of months, TheStYr/EdYr - start and end years '''
    ''' TheMDI - missing data indicator for masking '''
    ''' TheColls - dictionary of colours for each dataset '''
     
    # SOrt out region title
    if TheReg == 'Glob':
        RegString = 'Globe (70$^{o}$S to 70$^{o}$N)'
    elif TheReg == 'Nhem':
        RegString = 'N. Hemi (20$^{o}$N to 70$^{o}$N)'
    elif TheReg == 'Trop':
        RegString = 'Tropics (20$^{o}$S to 20$^{o}$N)'
    elif TheReg == 'Shem':
        RegString = 'S. Hemi (70$^{o}$S to 20$^{o}$S)'
    
    
    # for clearer plotting
    TheMCount=TheMCount+1
    TheYCount=TheYCount+1    
    
    # set up number of panels and number of lines
    nplots=len(TheParams[:])
    print('PLOT NUMBERS: ',nplots)
#    nlines=[]
#    for n in range(nplots):
#        print(n,TheLablees[n][:])
#        nlines.append(len(TheLablees[n][:]))
	
    Letteree=[]
    Letteree=LetterRange(0,nplots)
    
    # set up x axes
    if TheTimeType == 'monthly':
        TheMonths=[]
        yr=TheStYr
        mon=1
        for m in range(TheMCount):
            TheMonths.append(dt.date(yr,mon,1))
            mon=mon+1
            if mon == 13:
                mon=1
                yr=yr+1   
        TheMonths=np.array(TheMonths)	    
    else:
        TheMonths=[]
        yr=TheStYr
        mon=1
        for y in range(TheYCount):
            TheMonths.append(dt.date(yr,mon,1))
            yr=yr+1   
        TheMonths=np.array(TheMonths)	    
    
    xtitlee='Years'
            
    # set up dimensions and plot - this is a 3 column nvar rows plot
    # Panel 1
    xpos=[]
    ypos=[]
    xfat=[]
    ytall=[]
    totalyspace=0.90	# start 0.08 end 0.98
    totalxspace=0.85	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.12)
        ypos.append(0.98-((n+1)*(totalyspace/nplots)))
        xfat.append(totalxspace)
        ytall.append(totalyspace/nplots)
    
    f,axarr=plt.subplots(nplots,figsize=(7,12),sharex=False)	#6,18
    
    for pp in range(nplots):
        print('Plot: ',pp,TheParams[pp])
#        print(TheHvars[pp,0:10])
        print(TheHuncsC[pp,0:10])
        print(TheHuncsSp[pp,0:10])
        print(TheHuncsSt[pp,0:10])
	#axarr[pp].set_size(14)
        axarr[pp].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        if TheTimeType == 'monthly':
            axarr[pp].set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
        else:
            axarr[pp].set_xlim([TheMonths[0],TheMonths[TheYCount-1]])
        if pp < nplots-1:
            axarr[pp].set_xticklabels([])    
        if timetype == 'monthly':
            miny=min(np.reshape(TheHvars[pp,:],(TheMCount-1)))-max(np.reshape(TheHuncsC[pp,:],(TheMCount-1)))  
            maxy=max(np.reshape(TheHvars[pp,:],(TheMCount-1)))+max(np.reshape(TheHuncsC[pp,:],(TheMCount-1)))  
        else:
            miny=min(np.reshape(TheHvars[pp,:],(TheYCount-1)))-max(np.reshape(TheHuncsC[pp,:],(TheYCount-1)))  
            maxy=max(np.reshape(TheHvars[pp,:],(TheYCount-1)))+max(np.reshape(TheHuncsC[pp,:],(TheYCount-1)))  
        axarr[pp].set_ylim([(math.floor(10.*miny))/10.,(math.ceil(10.*maxy))/10.])
        if len(TheHuncsC[pp,:]) > 0:
            axarr[pp].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,:]+TheHuncsC[pp,:],TheHvars[pp,:]-TheHuncsC[pp,:],
                           facecolor='Gold',edgecolor='none')
            axarr[pp].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,:]+TheHuncsSp[pp,:],TheHvars[pp,:]-TheHuncsSp[pp,:],
                           facecolor='OrangeRed',edgecolor='none')
            axarr[pp].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,:]+TheHuncsSt[pp,:],TheHvars[pp,:]-TheHuncsSt[pp,:],
                           facecolor='DeepSkyBlue',edgecolor='none')
			   
        if timetype == 'monthly':
            axarr[pp].plot(TheMonths[0:len(TheMonths)-1],TheHvars[pp,:],
                c='black',linewidth=0.25)
        else:
            axarr[pp].plot(TheMonths[0:len(TheMonths)-1],TheHvars[pp,:],
                c='black',linewidth=0.5)
	
        axarr[pp].annotate(Letteree[pp]+') '+TheParams[pp],xy=(0.03,0.86),
                  xycoords='axes fraction',size=12)

# get linear trend and annotate (does this work with masked arrays?)
        lintrend=[0.,0.,0.] # median, 5th adn 95th percentile rate of change per time step
        if (TheTrendChoice == 'MPW'):
            lintrend=MedianPairwise(TheHvars[pp,:],TheMDI,lintrend)
        elif (TheTrendChoice == 'OLS'):
            slopes=OLS_AR1Corr(TheHvars[pp,:],TheMDI,0.9)
            lintrend[0] = slopes[0]
            lintrend[1] = slopes[1]
            lintrend[2] = slopes[2]
	    	
        #pdb.set_trace()
	
        if timetype == 'monthly':
            linstr="%5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (lintrend[0]*120,lintrend[1]*120,lintrend[2]*120,TheUnitees[pp])
        else:
            linstr="%5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (lintrend[0]*10,lintrend[1]*10,lintrend[2]*10,TheUnitees[pp])

        axarr[pp].annotate(linstr,xy=(0.5,0.08),xycoords='axes fraction',size=12,ha='center')
	 
        axarr[pp].set_ylabel(TheUnitees[pp],fontsize=12)
        if TheTimeType == 'monthly':
            axarr[pp].hlines(0,TheMonths[0],TheMonths[TheMCount-1],
                color='black',linewidth=0.5)
        else:
            axarr[pp].hlines(0,TheMonths[0],TheMonths[TheYCount-1],
                color='black',linewidth=0.5)
	
    axarr[0].annotate(RegString,xy=(0.5,0.84),
             xycoords='axes fraction',size=16,ha='center')
    axarr[nplots-1].set_xlabel(xtitlee,fontsize=12)
         
# Figure Watermark and Labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    

    #plt.show()
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")


     
    return #PlotNiceDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# set up loops to read in all time series

if timetype == 'monthly':
    varH=np.zeros((nparams,nmons))
    uncsHcov=np.zeros((nparams,nmons))
    uncsHsamp=np.zeros((nparams,nmons))
    uncsHstat=np.zeros((nparams,nmons))
    uncsHtot=np.zeros((nparams,nmons))
else:
    varH=np.zeros((nparams,nyrs))
    uncsHcov=np.zeros((nparams,nyrs))
    uncsHsamp=np.zeros((nparams,nyrs))
    uncsHstat=np.zeros((nparams,nyrs))
    uncsHtot=np.zeros((nparams,nyrs))

varH[:,:]=mdi
uncsHcov[:,:]=mdi
uncsHsamp[:,:]=mdi
uncsHstat[:,:]=mdi
uncsHtot[:,:]=mdi

for nv in range(nparams):
    tmpvar=[]
    tmpvarUcov=[]
    tmpvarUsamp=[]
    tmpvarUstat=[]
    tmpvarUtot=[]
    print('Reading in: ',param2[nv])
    FILIN=INHFILEED
# read in HadISDH time series
    if IfType == '.nc':
        MyNCFile=DATADIR+INHFILEST+param2[nv]+'.'+version+'_FLATgridHOM5by5'+FILIN+'.nc'
        f=netcdf.netcdf_file(MyNCFile,'r')
        if param[nv]=='q': 
            var=f.variables['glob_q_anoms']
        elif param[nv]=='e':
            var=f.variables['glob_e_anoms']
        elif param[nv]=='rh':
            var=f.variables['glob_RH_anoms']
        elif param[nv]=='t':
            var=f.variables['glob_T_anoms']
        elif param[nv]=='tw':
            var=f.variables['glob_Tw_anoms']
        elif param[nv]=='td':
            var=f.variables['glob_Td_anoms']
        elif param[nv]=='dpd':
            var=f.variables['glob_DPD_anoms']	    
        f.close()
        tmpvar=np.array(var.data)
    else:	# its a text file
        MyDatFile=DATADIR+INHFILEST+param2[nv]+FILIN
        MyTypes=("|S10","float","float","float","float","float")
        MyDelimiters=[10,11,11,11,11,11]
        MySkips=0
        RawData=ReadData(MyDatFile,MyTypes,MyDelimiters,MySkips)
        tmpvar=np.array(RawData['f1'])
        tmpvarUcov=np.array(RawData['f3'])
        tmpvarUsamp=np.array(RawData['f2'])
        tmpvarUstat=np.array(RawData['f4'])
        tmpvarUtot=np.array(RawData['f5'])
    
    # rezero HadISDH to non 1981-2010 anomalies if not 1981-2010 - ASSUME NO MISSING DATA!!!
    if (climst != 1981) | (climed != 2010):
        print('Renorming...')
        if timetype == 'monthly':
            tmpvar=np.reshape(tmpvar,(nyrs,12))
            for mm in range(12):
                subarr=tmpvar[:,mm]
                climarr=subarr[stcl:edcl]
                subarr[:]=subarr[:]-np.mean(climarr)
                tmpvar[:,mm]=subarr[:]
            varH[nv,:]=np.reshape(tmpvar,(1,nmons))
        else:
            climarr=tmpvar[stcl:edcl]
            tmpvar[:]=tmpvar[:]-np.mean(climarr)
            varH[nv,nr,:]=np.reshape(tmpvar,(1,nyrs))
       
    # Still need to populate main array   
    else:         
        varH[nv,:] = tmpvar[0:ntims]
	            
    if len(tmpvarUcov) > 0:
        uncsHcov[nv,:]=tmpvarUcov[0:ntims]
        uncsHsamp[nv,:]=tmpvarUsamp[0:ntims]
        uncsHstat[nv,:]=tmpvarUstat[0:ntims]
        uncsHtot[nv,:]=tmpvarUtot[0:ntims]
	
# convert to masked arrays and mask out missing data
print('Masking')
varH=np.ma.masked_array(varH)
varH[varH <= mdi]=np.ma.masked

uncsHcov=np.ma.masked_array(uncsHcov)
uncsHcov[uncsHcov <= mdi]=np.ma.masked

uncsHsamp=np.ma.masked_array(uncsHsamp)
uncsHsamp[uncsHsamp <= mdi]=np.ma.masked

uncsHstat=np.ma.masked_array(uncsHstat)
uncsHstat[uncsHstat <= mdi]=np.ma.masked

uncsHtot=np.ma.masked_array(uncsHtot)
uncsHtot[uncsHtot <= mdi]=np.ma.masked

# sort out in quadrature quantities for uncs where 
#	uncsHcov is total combined in quadrature
#       uncsHsamp is uncsHstat+uncsHsamp quadrature contributions
#	uncsHstat is uncsHstat quadrature contribution
print('Sorting out Uncs...')
for nv in range(nparams):
    RatsSamp=[]
    RatsStat=[]
    RatsSamp=(uncsHsamp[nv,:]**2)/((uncsHcov[nv,:]**2)+(uncsHsamp[nv,:]**2)+(uncsHstat[nv,:]**2))
    RatsStat=(uncsHstat[nv,:]**2)/((uncsHcov[nv,:]**2)+(uncsHsamp[nv,:]**2)+(uncsHstat[nv,:]**2))
    print(len(RatsSamp),len(RatsStat))
    uncsHcov[nv,:]=uncsHtot[nv,:]
    uncsHsamp[nv,:]=(uncsHtot[nv,:]*RatsSamp[:])+(uncsHtot[nv,:]*RatsStat[:])
    uncsHstat[nv,:]=uncsHtot[nv,:]*RatsStat[:]
    

# call plotter
print('Plotting...')
MyFile=PLOTDIR+OUTPLOT
PlotNiceTimeSeries(MyFile,varH,uncsHcov,uncsHsamp,uncsHstat,
                       unitees,nmons,nyrs,timetype,styr,edyr,mdi,
        	       param3,RegChoice,TrendChoice)
		
#    stop()

print("And, we are done!")

