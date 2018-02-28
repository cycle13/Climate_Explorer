#!/usr/local/sci/bin/python

#***************************************
# 1 April 2014 KMW - v1
# Plots station locations for each variable including removed sub/super sats  
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotRegionsTimeSeries_MAY2014.py
#
# REQUIRES
# RandomsRanges.py
# LinearTrends.py
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import scipy.stats
import struct
import os.path
import math
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
#from netCDF4 import Dataset
from scipy.io import netcdf
from scipy.stats.stats import pearsonr   
from RandomsRanges import LetterRange
from LinearTrends import MedianPairwise

# Set up initial run choices
timetype='annual'	#'monthly', 'annual'
nparams=7
param=list(['t','tw','td','q','e','rh','dpd'])	# tw, q, e, rh, t, td, dpd
param2=list(['T','Tw','Td','q','e','RH','DPD'])	# Tw, q, e, RH, T, Td, DPD
param3=list(['T','T$_{w}$','T$_{d}$','q','e','RH','DPD'])	# Tw, q, e, RH, T, Td, DPD
unitees=list(['$^{o}$C','$^{o}$C','$^{o}$C','g kg$^{-1}$','hPa','%rh','$^{o}$C'])
homogtype=list(['IDPHA','IDPHA','PHADPD','IDPHA','IDPHA','IDPHA','PHA'])	# 'IDPHA','PHA','PHADPD'
#param=list(['q','e','rh','tw','td','t','dpd'])	# tw, q, e, rh, t, td, dpd
#param2=list(['q','e','RH','Tw','Td','T','DPD'])	# Tw, q, e, RH, T, Td, DPD
#param3=list(['q','e','RH','T$_{w}$','T$_{d}$','T','DPD'])	# Tw, q, e, RH, T, Td, DPD
#unitees=list(['g kg$^{-1}$','hPa','%rh','$^{o}$C','$^{o}$C','$^{o}$C','$^{o}$C'])
#homogtype=list(['IDPHA','IDPHA','IDPHA','IDPHA','PHADPD','IDPHA','PHA'])	# 'IDPHA','PHA','PHADPD'

nowmon='MAY'
nowyear='2014'
thenmon='MAY'
thenyear='2014'
version='2.0.0.2013p'
styr=1973
edyr=2013
nyrs=(edyr-styr)+1
nmons=(nyrs)*12
climst=1981
climed=2010
stcl=climst-styr
edcl=climed-styr

# Set up directories and files

PLOTDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/'
DATADIR='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/'

IfType='.dat'	#'.nc'
INHFILEST='HadISDH.land'
#INHFILEED='5by5_'+thenmon+thenyear+'_areaTS_19732013'
if timetype == 'monthly':
    INHFILEEDG='.'+version+'_global_ts_monthly_'+thenmon+thenyear+'.dat'
else:
    INHFILEEDG='.'+version+'_global_ts_annual_'+thenmon+thenyear+'.dat'
if timetype == 'monthly':
    INHFILEEDN='.'+version+'_nhemi_ts_monthly_'+thenmon+thenyear+'.dat'
else:
    INHFILEEDN='.'+version+'_nhemi_ts_annual_'+thenmon+thenyear+'.dat'
if timetype == 'monthly':
    INHFILEEDS='.'+version+'_shemi_ts_monthly_'+thenmon+thenyear+'.dat'
else:
    INHFILEEDS='.'+version+'_shemi_ts_annual_'+thenmon+thenyear+'.dat'
if timetype == 'monthly':
    INHFILEEDT='.'+version+'_trop_ts_monthly_'+thenmon+thenyear+'.dat'
else:
    INHFILEEDT='.'+version+'_trop_ts_annual_'+thenmon+thenyear+'.dat'
OUTPLOT='PlotRegionsTimeSeries.'+version+'_'+timetype+'_'+nowmon+nowyear

# Set up variables
mdi=-1e30

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
		       TheStYr,TheEdYr,TheMDI,TheParams):
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
    Letteree=LetterRange(0,nplots*3)
    
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
    totalxspace=0.29	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.10)
	ypos.append(0.98-((n+1)*(totalyspace/nplots)))
	xfat.append(totalxspace)
	ytall.append(totalyspace/nplots)
    
    f,axarr=plt.subplots(21,figsize=(15,12),sharex=False)	#6,18
    
    for pp in range(nplots):
        print('Plot: ',pp,TheParams[pp])
#        print(TheHvars[pp,0,0:10])
        print(TheHuncsC[pp,0,0:10])
        print(TheHuncsSp[pp,0,0:10])
        print(TheHuncsSt[pp,0,0:10])
	#axarr[pp].set_size(14)
	axarr[pp].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        if TheTimeType == 'monthly':
	    axarr[pp].set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
	else:
	    axarr[pp].set_xlim([TheMonths[0],TheMonths[TheYCount-1]])
	if pp < nplots-1:
	    axarr[pp].set_xticklabels([])    
	miny=min(np.reshape(TheHvars[pp,:,:],(TheYCount-1)*3))-max(np.reshape(TheHuncsC[pp,:,:],(TheYCount-1)*3))  
	maxy=max(np.reshape(TheHvars[pp,:,:],(TheYCount-1)*3))+max(np.reshape(TheHuncsC[pp,:,:],(TheYCount-1)*3))  
	axarr[pp].set_ylim([(math.floor(10.*miny))/10.,(math.ceil(10.*maxy))/10.])
	if len(TheHuncsC[pp,0,:]) > 0:
	    axarr[pp].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,0,:]+TheHuncsC[pp,0,:],TheHvars[pp,0,:]-TheHuncsC[pp,0,:],
	                   facecolor='Gold',edgecolor='none')
	    axarr[pp].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,0,:]+TheHuncsSp[pp,0,:],TheHvars[pp,0,:]-TheHuncsSp[pp,0,:],
	                   facecolor='OrangeRed',edgecolor='none')
	    axarr[pp].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,0,:]+TheHuncsSt[pp,0,:],TheHvars[pp,0,:]-TheHuncsSt[pp,0,:],
	                   facecolor='DeepSkyBlue',edgecolor='none')
			   
	if timetype == 'monthly':
	    axarr[pp].plot(TheMonths[0:len(TheMonths)-1],TheHvars[pp,0,:],
	        c='black',linewidth=0.25)
	else:
	    axarr[pp].plot(TheMonths[0:len(TheMonths)-1],TheHvars[pp,0,:],
	        c='black',linewidth=0.5)
	
	axarr[pp].annotate(Letteree[pp]+') '+TheParams[pp],xy=(0.03,0.86),
	          xycoords='axes fraction',size=12)

# get linear trend and annotate (does this work with masked arrays?)
        lintrend=[0.,0.,0.] # median, 5th adn 95th percentile rate of change per time step
        lintrend=MedianPairwise(TheHvars[pp,0,:],TheMDI,lintrend)
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
	
    axarr[0].annotate('N. Hemi (20$^{o}$N to 70$^{0}$N)',xy=(0.5,0.84),
             xycoords='axes fraction',size=16,ha='center')
    axarr[nplots-1].set_xlabel(xtitlee,fontsize=12)
     
    # Panel 2
    xpos=[]
    ypos=[]
    xfat=[]
    ytall=[]
    totalyspace=0.90	# start 0.08 end 0.98
    totalxspace=0.29	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.39)
	ypos.append(0.98-((n+1)*(totalyspace/nplots)))
	xfat.append(totalxspace)
	ytall.append(totalyspace/nplots)
    
    
    for pp in range(nplots):
        print('Plot: ',pp+nplots,TheParams[pp])
	axarr[pp+nplots].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        if TheTimeType == 'monthly':
	    axarr[pp+nplots].set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
	else:
	    axarr[pp+nplots].set_xlim([TheMonths[0],TheMonths[TheYCount-1]])
	if pp < nplots-1:
	    axarr[pp+nplots].set_xticklabels([])    
	miny=min(np.reshape(TheHvars[pp,:,:],(TheYCount-1)*3))-max(np.reshape(TheHuncsC[pp,:,:],(TheYCount-1)*3))  
	maxy=max(np.reshape(TheHvars[pp,:,:],(TheYCount-1)*3))+max(np.reshape(TheHuncsC[pp,:,:],(TheYCount-1)*3))  
	axarr[pp+nplots].set_ylim([(math.floor(10.*miny))/10.,(math.ceil(10.*maxy))/10.])

	axarr[pp+nplots].set_yticklabels([])		    
	if len(TheHuncsC[pp,1,:]) > 0:
	    axarr[pp+nplots].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,1,:]+TheHuncsC[pp,1,:],TheHvars[pp,1,:]-TheHuncsC[pp,1,:],
	                   facecolor='Gold',edgecolor='none')
	    axarr[pp+nplots].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,1,:]+TheHuncsSp[pp,1,:],TheHvars[pp,1,:]-TheHuncsSp[pp,1,:],
	                   facecolor='OrangeRed',edgecolor='none')
	    axarr[pp+nplots].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,1,:]+TheHuncsSt[pp,1,:],TheHvars[pp,1,:]-TheHuncsSt[pp,1,:],
	                   facecolor='DeepSkyBlue',edgecolor='none')
			   
	if timetype == 'monthly':
	    axarr[pp+nplots].plot(TheMonths[0:len(TheMonths)-1],TheHvars[pp,1,:],
	        c='black',linewidth=0.25)
	else:
	    axarr[pp+nplots].plot(TheMonths[0:len(TheMonths)-1],TheHvars[pp,1,:],
	        c='black',linewidth=0.5)
	
	axarr[pp+nplots].annotate(Letteree[pp+nplots]+') '+TheParams[pp],xy=(0.03,0.86),
	                 xycoords='axes fraction',size=12)

# get linear trend and annotate (does this work with masked arrays?)
        lintrend=[0.,0.,0.] # median, 5th adn 95th percentile rate of change per time step
        lintrend=MedianPairwise(TheHvars[pp,1,:],TheMDI,lintrend)
        if timetype == 'monthly':
	    linstr="%5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (lintrend[0]*120,lintrend[1]*120,lintrend[2]*120,TheUnitees[pp])
        else:
	    linstr="%5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (lintrend[0]*10,lintrend[1]*10,lintrend[2]*10,TheUnitees[pp])
	axarr[pp+nplots].annotate(linstr,xy=(0.5,0.08),xycoords='axes fraction',size=12,ha='center')

#        axarr[pp+nplots].set_ylabel(TheUnitees[pp],fontsize=10)
        if TheTimeType == 'monthly':
	    axarr[pp+nplots].hlines(0,TheMonths[0],TheMonths[TheMCount-1],
	        color='black',linewidth=0.5)
        else:
	    axarr[pp+nplots].hlines(0,TheMonths[0],TheMonths[TheYCount-1],
	        color='black',linewidth=0.5)

    axarr[7].annotate('S. Hemi (70$^{o}$S to 20$^{0}$S)',xy=(0.5,0.84),
             xycoords='axes fraction',size=16,ha='center')	
    axarr[(nplots*2)-1].set_xlabel(xtitlee,fontsize=12)

    # Panel 3
    xpos=[]
    ypos=[]
    xfat=[]
    ytall=[]
    totalyspace=0.90	# start 0.08 end 0.98
    totalxspace=0.29	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.68)
	ypos.append(0.98-((n+1)*(totalyspace/nplots)))
	xfat.append(totalxspace)
	ytall.append(totalyspace/nplots)
    
    
    for pp in range(nplots):
        print('Plot: ',pp+nplots+nplots,TheParams[pp])
	axarr[pp+nplots+nplots].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        if TheTimeType == 'monthly':
	    axarr[pp+nplots+nplots].set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
	else:
	    axarr[pp+nplots+nplots].set_xlim([TheMonths[0],TheMonths[TheYCount-1]])
	if pp < nplots-1:
	    axarr[pp+nplots+nplots].set_xticklabels([])    
	miny=min(np.reshape(TheHvars[pp,:,:],(TheYCount-1)*3))-max(np.reshape(TheHuncsC[pp,:,:],(TheYCount-1)*3))  
	maxy=max(np.reshape(TheHvars[pp,:,:],(TheYCount-1)*3))+max(np.reshape(TheHuncsC[pp,:,:],(TheYCount-1)*3))  
	axarr[pp+nplots+nplots].set_ylim([(math.floor(10.*miny))/10.,(math.ceil(10.*maxy))/10.])

	axarr[pp+nplots+nplots].set_yticklabels([])		    
	if len(TheHuncsC[pp,2,:]) > 0:
	    axarr[pp+nplots+nplots].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,2,:]+TheHuncsC[pp,2,:],TheHvars[pp,2,:]-TheHuncsC[pp,2,:],
	                   facecolor='Gold',edgecolor='none')
	    axarr[pp+nplots+nplots].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,2,:]+TheHuncsSp[pp,2,:],TheHvars[pp,2,:]-TheHuncsSp[pp,2,:],
	                   facecolor='OrangeRed',edgecolor='none')
	    axarr[pp+nplots+nplots].fill_between(TheMonths[0:len(TheMonths)-1],TheHvars[pp,2,:]+TheHuncsSt[pp,2,:],TheHvars[pp,2,:]-TheHuncsSt[pp,2,:],
	                   facecolor='DeepSkyBlue',edgecolor='none')
			   
	if timetype == 'monthly':
	    axarr[pp+nplots+nplots].plot(TheMonths[0:len(TheMonths)-1],
	        TheHvars[pp,2,:],c='black',linewidth=0.25)
	else:
	    axarr[pp+nplots+nplots].plot(TheMonths[0:len(TheMonths)-1],
	        TheHvars[pp,2,:],c='black',linewidth=0.5)
	
	axarr[pp+nplots+nplots].annotate(Letteree[pp+nplots+nplots]+') '+TheParams[pp],
	          xy=(0.03,0.86),xycoords='axes fraction',size=12)

# get linear trend and annotate (does this work with masked arrays?)
        lintrend=[0.,0.,0.] # median, 5th adn 95th percentile rate of change per time step
        lintrend=MedianPairwise(TheHvars[pp,2,:],TheMDI,lintrend)
        if timetype == 'monthly':
	    linstr="%5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (lintrend[0]*120,lintrend[1]*120,lintrend[2]*120,TheUnitees[pp])
        else:
	    linstr="%5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (lintrend[0]*10,lintrend[1]*10,lintrend[2]*10,TheUnitees[pp])

	axarr[pp+nplots+nplots].annotate(linstr,xy=(0.5,0.08),xycoords='axes fraction',
	        size=12,ha='center')

#        axarr[pp+nplots].set_ylabel(TheUnitees[pp],fontsize=10)
        if TheTimeType == 'monthly':
	    axarr[pp+nplots+nplots].hlines(0,TheMonths[0],
	        TheMonths[TheMCount-1],color='black',linewidth=0.5)
        else:
	    axarr[pp+nplots+nplots].hlines(0,TheMonths[0],
	        TheMonths[TheYCount-1],color='black',linewidth=0.5)
	
    axarr[14].annotate('Tropics (20$^{o}$S to 20$^{0}$N)',xy=(0.5,0.84),
             xycoords='axes fraction',size=16,ha='center')
    axarr[(nplots*3)-1].set_xlabel(xtitlee,fontsize=12)
    
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
    varH=np.zeros((nparams,3,nmons))
    uncsHcov=np.zeros((nparams,3,nmons))
    uncsHsamp=np.zeros((nparams,3,nmons))
    uncsHstat=np.zeros((nparams,3,nmons))
    uncsHtot=np.zeros((nparams,3,nmons))
else:
    varH=np.zeros((nparams,3,nyrs))
    uncsHcov=np.zeros((nparams,3,nyrs))
    uncsHsamp=np.zeros((nparams,3,nyrs))
    uncsHstat=np.zeros((nparams,3,nyrs))
    uncsHtot=np.zeros((nparams,3,nyrs))

varH[:,:,:]=mdi
uncsHcov[:,:,:]=mdi
uncsHsamp[:,:,:]=mdi
uncsHstat[:,:,:]=mdi
uncsHtot[:,:,:]=mdi

for nr in range(3):
    for nv in range(nparams):
        tmpvar=[]
        tmpvarUcov=[]
        tmpvarUsamp=[]
        tmpvarUstat=[]
        tmpvarUtot=[]
        print('Reading in: ',param2[nv])
        if nr == 0:
	    FILIN=INHFILEEDN
        if nr == 1:
	    FILIN=INHFILEEDS
        if nr == 2:
	    FILIN=INHFILEEDT
# read in HadISDH time series
        if IfType == '.nc':
            MyNCFile=DATADIR+INHFILEST+param2[nv]+'.'+version+'_FLATgrid'+homogtype[nv]+FILIN+'.nc'
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
            MyDelimiters=[10,10,10,10,10,10]
	    MySkips=1
            RawData=ReadData(MyDatFile,MyTypes,MyDelimiters,MySkips)
            tmpvar=np.array(RawData['f1'])
            tmpvarUcov=np.array(RawData['f3'])
            tmpvarUsamp=np.array(RawData['f2'])
            tmpvarUstat=np.array(RawData['f4'])
            tmpvarUtot=np.array(RawData['f5'])
    
    # rezero HadISDH to 1981-2010 anomalies - ASSUME NO MISSING DATA!!!
        print('Renorming...')
        if timetype == 'monthly':
            tmpvar=np.reshape(tmpvar,(nyrs,12))
            for mm in range(12):
                subarr=tmpvar[:,mm]
	        climarr=subarr[stcl:edcl]
	        subarr[:]=subarr[:]-np.mean(climarr)
	        tmpvar[:,mm]=subarr[:]
            varH[nv,nr,:]=np.reshape(tmpvar,(1,nmons))
        else:
            climarr=tmpvar[stcl:edcl]
	    tmpvar[:]=tmpvar[:]-np.mean(climarr)
            varH[nv,nr,:]=np.reshape(tmpvar,(1,nyrs))
        
        if len(tmpvarUcov) > 0:
            uncsHcov[nv,nr,:]=tmpvarUcov
            uncsHsamp[nv,nr,:]=tmpvarUsamp
            uncsHstat[nv,nr,:]=tmpvarUstat
            uncsHtot[nv,nr,:]=tmpvarUtot
	
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
for nr in range(3):
    for nv in range(nparams):
        RatsSamp=[]
        RatsStat=[]
        RatsSamp=(uncsHsamp[nv,nr,:]**2)/((uncsHcov[nv,nr,:]**2)+(uncsHsamp[nv,nr,:]**2)+(uncsHstat[nv,nr,:]**2))
        RatsStat=(uncsHstat[nv,nr,:]**2)/((uncsHcov[nv,nr,:]**2)+(uncsHsamp[nv,nr,:]**2)+(uncsHstat[nv,nr,:]**2))
        print(len(RatsSamp),len(RatsStat))
        uncsHcov[nv,nr,:]=uncsHtot[nv,nr,:]
        uncsHsamp[nv,nr,:]=(uncsHtot[nv,nr,:]*RatsSamp[:])+(uncsHtot[nv,nr,:]*RatsStat[:])
        uncsHstat[nv,nr,:]=uncsHtot[nv,nr,:]*RatsStat[:]
    

# call plotter
print('Plotting...')
MyFile=PLOTDIR+OUTPLOT
PlotNiceTimeSeries(MyFile,varH,uncsHcov,uncsHsamp,uncsHstat,
                       unitees,nmons,nyrs,timetype,styr,edyr,mdi,
		       param3)
		
#    stop()

print("And, we are done!")

