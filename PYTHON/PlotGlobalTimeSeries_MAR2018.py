#!/usr/local/sci/bin/python

#***************************************
# 7 March 2018 KMW - v2
# Plots global time series from ASCII
# Can plot against other datasets
# Requires uncertainty values for HadISDH
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotGlobalTimeSeries_MAR2018.py
#
# REQUIRES
# 
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
from RandomsRanges import LetterRange

# Set up initial run choices
timetype='annual'	#'monthly', 'annual'
nparams=7
param=list(['q','e','rh','tw','td','t','dpd'])	# tw, q, e, rh, t, td, dpd
param2=list(['q','e','RH','Tw','Td','T','DPD'])	# Tw, q, e, RH, T, Td, DPD
unitees=list(['g/kg','hPa','%rh','degrees C','degrees C','degrees C','degrees C'])
homogtype=list(['IDPHA','IDPHA','IDPHA','IDPHA','PHADPD','IDPHA','PHA'])	# 'IDPHA','PHA','PHADPD'
sourceslist=list(['ERAINTERIM','CRUTS3_21','CRUTEM4','GHCNM'])
#others=list([['ERAINTERIM'],
#             ['ERAINTERIM','CRUTS3_21'],
#	     ['ERAINTERIM'],
#	     ['ERAINTERIM'],
#	     ['ERAINTERIM'],
#	     ['ERAINTERIM','CRUTEM4','GHCNM'],
#	     ['ERAINTERIM']])           
others=list([[],
             [],
	     [],
	     [],
	     [],
	     [],
	     []])	    
nowmon='MAR'
nowyear='2018'
thenmon='JAN'
thenyear='2018'
version='4.0.0.2017f'
styr=1973
edyr=2017
nyrs=(edyr-styr)+1
nmons=(nyrs)*12
climst=1981
climed=2010
stcl=climst-styr
edcl=climed-styr

# Set up directories and files

PLOTDIR='/data/local/hadkw/HADCRUH2/UPDATE'+str(edyr)+'/IMAGES/TIMESERIES/'
DATADIR='/data/local/hadkw/HADCRUH2/UPDATE'+str(edyr)+'/STATISTICS/TIMESERIES/'

IfType='.dat'	#'.nc'
INHFILEST='HadISDH.land'
#INHFILEED='5by5_'+thenmon+thenyear+'_areaTS_19732013'
if timetype == 'monthly':
    INHFILEED='.'+version+'_global_ts_monthly_anoms8110_'+thenmon+thenyear+'.dat'
else:
    INHFILEED='.'+version+'_global_ts_annual_anoms8110_'+thenmon+thenyear+'.dat'
INOTHFULL='_areaTS_1973'+str(edyr)
INOTHMASK='_HadISDHMASKarea'+str(edyr)
OUTPLOT='PlotNiceTimeSeries.'+version+'_'+timetype+'_anoms8110_'+nowmon+nowyear

# Set up variables
mdi=-1e30

varH=[]	# nvars(rows) by mons masked array
uncsHtot=[] # nvars(rows) by mons masked array
uncsHcov=[] # nvars(rows) by mons masked array
uncsHsamp=[] # nvars(rows) by mons masked array
uncsHstat=[] # nvars(rows) by mons masked array
othervarsFULL=[] # nvars(rows) by nothers (max) by mons masked array
othervarsMASK=[] # nvars(rows) by nothers (max) by mons masked array (HadISDH mask too)


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
                       TheLablees,TheUnitees,TheVars,TheMASKVars,
		       TheMCount,TheYCount,TheTimeType,TheStYr,TheEdYr,TheMDI,TheColls,TheParams):
    ''' Plot a panel for each element of TheHvars '''
    ''' Add Coverage, Sampling and Station uncertainty ranges '''
    ''' Add lines for any extra estimates (TheVars) and HadISDH MASKED versions '''
    ''' Save as png and eps '''
    ''' TheHvars is a multi-row array: rows for vars, columns for months '''
    ''' Ditto TheHuncs C=coverage, Sp=sampling, St=station '''
    ''' TheLablees is the name list for all other vars '''
    ''' TheUnitees is the units name list '''
    ''' TheVars is a multirow array if there is 1+ var available - or [] '''
    ''' TheMASKVars - ditto above but masked to HadISDH coverage '''
    ''' TheMCount - number of months, TheStYr/EdYr - start and end years '''
    ''' TheMDI - missing data indicator for masking '''
    ''' TheColls - dictionary of colours for each dataset '''
        
    
    # set up number of panels and number of lines
    nplots=len(TheParams[:])
    print('PLOT NUMBERS: ',nplots)
    nlines=[]
    for n in range(nplots):
        print(n,TheLablees[n][:])
        nlines.append(len(TheLablees[n][:]))
	
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
            
    # set up dimensions and plot
    xpos=[]
    ypos=[]
    xfat=[]
    ytall=[]
    totalyspace=0.90	# start 0.08 end 0.98
    totalxspace=0.84	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.14)
	ypos.append(0.98-((n+1)*(totalyspace/nplots)))
	xfat.append(totalxspace)
	ytall.append(totalyspace/nplots)
    
#    plt.clf()
#    fig = plt.figure(1,figsize=(8,12))
#    plt.axes([0.15,0.1,0.8,0.80])
    f,axarr=plt.subplots(7,figsize=(6,12),sharex=True)	#6,18
    
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
	axarr[pp].set_ylim([math.floor(min(TheHvars[pp,:]-TheHuncsC[pp,:])),
	                    math.ceil(max(TheHvars[pp,:]+TheHuncsC[pp,:]))])
	if len(TheHuncsC[pp,:]) > 0:
	    axarr[pp].fill_between(TheMonths,TheHvars[pp,:]+TheHuncsC[pp,:],TheHvars[pp,:]-TheHuncsC[pp,:],
	                   facecolor='LightGray',edgecolor='none')
	    axarr[pp].fill_between(TheMonths,TheHvars[pp,:]+TheHuncsSp[pp,:],TheHvars[pp,:]-TheHuncsSp[pp,:],
	                   facecolor='LightSlateGray',edgecolor='none')
	    axarr[pp].fill_between(TheMonths,TheHvars[pp,:]+TheHuncsSt[pp,:],TheHvars[pp,:]-TheHuncsSt[pp,:],
	                   facecolor='LightSlateGray',edgecolor='none')
			   
	if timetype == 'monthly':
	    axarr[pp].plot(TheMonths,TheHvars[pp,:],c='black',linewidth=0.5)
	else:
	    axarr[pp].plot(TheMonths,TheHvars[pp,:],c='black',linewidth=2)
	
	axarr[pp].annotate(Letteree[pp]+') '+TheParams[pp],xy=(0.03,0.9),xycoords='axes fraction',size=10)
	axarr[pp].annotate('HadISDH',xy=(0.14,0.9),xycoords='axes fraction',color='black',size=10)
	for ll in range(nlines[pp]):	# no problem if 0
            print('Other: ',ll,TheLablees[pp][ll])
#	    axarr[pp].plot(TheMonths,TheVars[pp,ll,:],c=TheColls[TheLablees[pp][ll]],linewidth=0.5)
#	    axarr[pp].plot(TheMonths,TheMASKVars[pp,ll,:],c=TheColls[TheLablees[pp][ll]],linestyle='dotted',linewidth=0.5)
	    if timetype == 'monthly':
	        axarr[pp].plot(TheMonths,TheMASKVars[pp,ll,:],c=TheColls[TheLablees[pp][ll]],linewidth=0.5)
	    else:
	        axarr[pp].plot(TheMonths,TheMASKVars[pp,ll,:],c=TheColls[TheLablees[pp][ll]],linewidth=2)
	    
	    axarr[pp].annotate(TheLablees[pp][ll],xy=(0.14,0.82-(ll*0.08)),xycoords='axes fraction',
	               color=TheColls[TheLablees[pp][ll]],size=10)
        axarr[pp].set_ylabel(TheUnitees[pp],fontsize=10)
        if TheTimeType == 'monthly':
	    axarr[pp].hlines(0,TheMonths[0],TheMonths[TheMCount-1],color='black')
        else:
	    axarr[pp].hlines(0,TheMonths[0],TheMonths[TheYCount-1],color='black')
	
    axarr[nplots-1].set_xlabel(xtitlee,fontsize=10)
     
    
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
    othervarsFULL=np.zeros((nparams,3,nmons))
    othervarsMASK=np.zeros((nparams,3,nmons))
else:
    varH=np.zeros((nparams,nyrs))
    uncsHcov=np.zeros((nparams,nyrs))
    uncsHsamp=np.zeros((nparams,nyrs))
    uncsHstat=np.zeros((nparams,nyrs))
    uncsHtot=np.zeros((nparams,nyrs))
    othervarsFULL=np.zeros((nparams,3,nyrs))
    othervarsMASK=np.zeros((nparams,3,nyrs))

varH[:,:]=mdi
uncsHcov[:,:]=mdi
uncsHsamp[:,:]=mdi
uncsHstat[:,:]=mdi
uncsHtot[:,:]=mdi
othervarsFULL[:,:,:]=mdi
othervarsMASK[:,:,:]=mdi

for nv in range(nparams):
    tmpvar=[]
    tmpvarUcov=[]
    tmpvarUsamp=[]
    tmpvarUstat=[]
    tmpvarUtot=[]
    print('Reading in: ',param2[nv])
# read in HadISDH time series
    if IfType == '.nc':
        MyNCFile=DATADIR+INHFILEST+param2[nv]+'.'+version+'_FLATgrid'+homogtype[nv]+INHFILEED+'.nc'
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
	tmpvar=np.array(var.data)
        f.close()
    else:	# its a text file
        MyDatFile=DATADIR+INHFILEST+param2[nv]+INHFILEED
        MyTypes=("|S10","float","float","float","float","float")
        MyDelimiters=[10,10,10,10,10,10]
	MySkips=1
        RawData=ReadData(MyDatFile,MyTypes,MyDelimiters,MySkips)
        tmpvar=np.array(RawData['f1'])
        tmpvarUcov=np.array(RawData['f3'])
        tmpvarUsamp=np.array(RawData['f2'])
        tmpvarUstat=np.array(RawData['f4'])
        tmpvarUtot=np.array(RawData['f5'])
    
# rezero HadISDH to non 1981-2010 anomalies if necessary - ASSUME NO MISSING DATA!!!
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
            varH[nv,:]=np.reshape(tmpvar,(1,nyrs))
        
        if len(tmpvarUcov) > 0:
            uncsHcov[nv,:]=tmpvarUcov
            uncsHsamp[nv,:]=tmpvarUsamp
            uncsHstat[nv,:]=tmpvarUstat
            uncsHtot[nv,:]=tmpvarUtot
# GOT TO HERE	
# read in all time series
    for no in range(len(others[nv][:])):
        print('Reading in others: ',no,others[nv][no])
	MyNCFile=DATADIR+others[nv][no]+'_'+param2[nv]+INOTHFULL+'.nc'
        f=netcdf.netcdf_file(MyNCFile,'r')
        var=f.variables['glob_anoms']
	newvar=np.array(var.data)
        f.close()
	if timetype == 'annual':
	    newvar=np.reshape(newvar,(nyrs,12))
	    for yy in range(nyrs):
	        if newvar[yy,0] > mdi:
		    othervarsFULL[nv,no,yy]=np.mean(newvar[yy,:])
	else:
            othervarsFULL[nv,no,:]=np.reshape(newvar,(1,nmons))

        MyNCFile=DATADIR+others[nv][no]+'_'+param2[nv]+INOTHMASK+'.nc'
        f=netcdf.netcdf_file(MyNCFile,'r')
        var=f.variables['glob_anoms']
	newvar=np.array(var.data)
        f.close()
	if timetype == 'annual':
	    newvar=np.reshape(newvar,(nyrs,12))
	    for yy in range(nyrs):
		if newvar[yy,0] > mdi:
	            othervarsMASK[nv,no,yy]=np.mean(newvar[yy,:])
	else:
            othervarsMASK[nv,no,:]=np.reshape(newvar,(1,nmons))

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

othervarsFULL=np.ma.masked_array(othervarsFULL)
othervarsFULL[othervarsFULL <= mdi]=np.ma.masked

othervarsMASK=np.ma.masked_array(othervarsMASK)
othervarsMASK[othervarsMASK <= mdi]=np.ma.masked

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
    
# set up colour dictionary - so that each dataset has an associated colour
diccols={}
diccols[sourceslist[0]]='Red'
diccols[sourceslist[1]]='MediumBlue'
diccols[sourceslist[2]]='DarkOrange'
diccols[sourceslist[3]]='MediumSlateBlue'

# call plotter
print('Plotting...')
MyFile=PLOTDIR+OUTPLOT
PlotNiceTimeSeries(MyFile,varH,uncsHcov,uncsHsamp,uncsHstat,
                       others,unitees,othervarsFULL,othervarsMASK,
		       nmons,nyrs,timetype,styr,edyr,mdi,
		       diccols,param2)
		
#    stop()

print("And, we are done!")

