#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 8th October 2015
# Last update: 8th October 2015
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	# this will probably change
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Reads in monthly mean anomaly regional average time series for q, T and RH from HadISDH
# Can plot monthly or annual data
# Can plot one region or all four
# Plots a T q scatter with each year as the point (or MONYY for monthly)
# Colours the points by simultaneous RH value
# Plots RH colour bar to the right
# Adds Tq and TRH correlation to plot
# Adds Tq and TRH slope to plot
#
# NO MISSING DATA IN TIME SERIES!!!!
# 
# <references to related published material, e.g. that describes data set>
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
# from mpl_toolkits.basemap import Basemap
# import datetime as dt
# from matplotlib.dates import date2num,num2date
# from scipy.io import netcdf
# import matplotlib.colors as mc
# import matplotlib.cm as mpl_cm
# import pdb
#
# Other:
# ReadNetCDFTS - infile function to read in netCDF timeseries, written by Kate Willett
# PlotScatter - infile function to plot, written by Kate Willett
# 
# -----------------------
# DATA
# -----------------------
# directory for regional timeseries:
# /data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/
# files currently worked on:
# Specific humidity:
# HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc
# Relative humidity:
# HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc
# Temperature:
# HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# Select 'TimeRes'  to be 'M' or 'Y' for month or year
# Ensure correct file paths and files
# Ensure start year (styr) and end year (edyr) are correct
# Select 'Region' to be 'A', 'G','N','T' or 'S' for All, Globe, NHemi, Tropics, SHemi
#
# run:
# python2.7 PlotTqRhScatter_OCT2015.py
# 
# -----------------------
# OUTPUT
# -----------------------
# directory for output images:
# /data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/ANALYSIS/
# Output image file: (nowmon+nowyear= e.g., OCT2015):
# ScatterTqRH_HadISDH.landq.2.0.1.2014p_'+nowmon+nowyear+
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 8 October 2015
# ---------
#  
# Enhancements
#  
# Changes
#  
# Bug fixes
#  
# -----------------------
# OTHER INFORMATION
# -----------------------
#
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys, os
import scipy.stats
import struct
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from scipy.io import netcdf
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import pdb #stop: pdb.set_trace(), start: c

# Set up initial run choices
TimeRes='Y'	# M=month, Y=year	
Region='A'	# A=All, G=Globe, N=NHemi, T=Tropics, S=SHemi
homogtype='IDPHA'	# 'IDPHA','PHA','PHADPD'
nowmon='OCT'
nowyear='2015'

styr=1973
edyr=2014
nyrs=(edyr-styr)+1
nmons=(nyrs)*12
if (TimeRes == 'Y'):
    ntims=nyrs
else:
    ntims=nmons    
YrStr=np.array(range(styr,edyr+1),dtype=str)
YrStr=np.array(([i[2:5] for i in YrStr])) # now a string array of the last two digits

# Set up directories and files
INDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/ANALYSIS/'

In_q='HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc'
In_RH='HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc'
In_T='HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc'

OutPlot='ScatterTqRH_HadISDH.2.0.1.2014p_'+Region+'_'+nowmon+nowyear
 
# Set up variables
q_arr=0	#set once file read in
T_arr=0		#set once file read in
RH_arr=0		#set once file read in

#************************************************************************
# Subroutines
#************************************************************************
# READNETCDFTS
def ReadNetCDFTS(FileName,ReadInfo,TheData):
    ''' Open the NetCDF File
        Get the data 
	FileName: stroing containing filepath/name
	TheData: an empty 2D array big enough for 1 or 4 regions worth of data
	ReadInfo: list of 1 or 4 strings of variable name/s for the globe, N Hemi, Tropics and S.Hemi '''

    ncf=netcdf.netcdf_file(FileName,'r')

    # ncf.variables this lists the variable names
    for loo in range(len(ReadInfo)):
        print(loo)
        var=ncf.variables[ReadInfo[loo]]
        TheData[loo,:]=np.array(var.data)

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData # ReadNetCDFTS

#************************************************************************
# PlotScatter
def PlotScatter(TheFile,TheYrStr,Thentims,Theq_arr,TheRH_arr,TheT_arr,TheReg):
    ''' Plot Tq scatter with colours related to RH'''
    ''' Points are either the last two years YY or MONYY '''
    ''' Save as png and eps '''
    ''' TheFile - the filepath and filename for the image '''
    ''' TheYrStr - a string array of the last two digits for years NYrs long '''
    ''' Thentims - an integer for the number of points to be plotted '''
    ''' Theq_arr - the specific humidity data (can be monthly or yearly ''' 
    ''' TheRH_arr - the relative humidity data (can be monthly or yearly ''' 
    ''' TheT_arr - the temperature data (can be monthly or yearly ''' 
        
    # set up the points if monthly - concatenating MONYY
    if (ntims > len(TheYrStr)):
        MONLABS=np.array(('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'))
	Pointees=[i+j for i in TheYrStr for j in MONLABS]
    else:
        Pointees=TheYrStr
		
    # Load colours and set up bounds
    cmap=plt.get_cmap('BrBG') # BrownBlueGreen
        
    cmaplist=[cmap(i) for i in range(cmap.N)]
    for loo in range((cmap.N/2)-30,(cmap.N/2)+30):
        cmaplist.remove(cmaplist[(cmap.N/2)-30]) # remove the very pale colours in the middle
#    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle
#
## remove the darkest and lightest (white and black) - and reverse
#    for loo in range(40):
#        cmaplist.remove(cmaplist[0])
##    cmaplist.reverse()
##    for loo in range(10):
##        cmaplist.remove(cmaplist[0])
##    cmaplist.reverse()

    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    # work out best max and min values for colourbar and try to make them 'nice
    # must be an odd number of steps
    # for less than 0.1 must be 14 or fewer steps
    vmax=np.int(np.ceil(np.max(abs(TheRH_arr))*10))/10.    
    vmin=-vmax
    if (vmax <= 0.3):
        nsteps=np.int((vmax-vmin)/0.05)+1
    elif (vmax <= 0.5):
        vmax=np.ceil(np.max(abs(TheRH_arr))/0.06)*0.06  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.06)+1
    elif (vmax <= 1.0):
        vmax=np.ceil(np.max(abs(TheRH_arr))/0.1)*0.1  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.1)+1
    elif (vmax <= 1.4):
        vmax=np.ceil(np.max(abs(TheRH_arr))/0.2)*0.2  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.2)+1
    elif (vmax > 1.4):
        vmax=np.ceil(np.max(abs(TheRH_arr))/0.3)*0.3  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.3)+1
	#    pdb.set_trace() # stop here and play
    
    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)
        
    ytitlee='Specific Humidity Anomalies (g kg$^{-1}$)'
    xtitlee='Temperature Anomalies ($^{o}$C)'
    titlees=['Globe 70$^{o}$S to 70$^{o}$N','N. Hemisphere 20$^{o}$N to 70$^{o}$N','Tropics 20$^{o}$S to 20$^{o}$N','S. Hemisphere 70$^{o}$S to 20$^{o}$S']

    # set up max and min of q and T for axes - keep same for all regions
    qmax=np.ceil(np.max(abs(Theq_arr))/0.1)*0.1  
    qmin=-qmax
    tmax=np.ceil(np.max(abs(TheT_arr))/0.1)*0.1  
    tmin=-tmax
                
    # set up plot - are we working with one region or four?
    if (TheReg != 'A'):

        # Single plot scenario
        fig = plt.figure(1,figsize=(8,8))
        plt.clf()	# needs to be called after figure!!! (create the figure, then clear the plot space)

        ax1=plt.axes([0.1,0.1,0.75,0.8]) # left, bottom, width, height

        ax1.set_xlim([tmin,tmax])   
        ax1.set_ylim([qmin,qmax])   

        # make blank plot with zero lines on
	ax1.plot(np.zeros(100),np.linspace(qmin,qmax,100),color='black',linewidth=2)
	ax1.plot(np.linspace(tmin,tmax,100),np.zeros(100),color='black',linewidth=2)

        # plot black dots for the goods
        for vv in range(Thentims):
	    #print(vv,TheT_arr[0,vv],Theq_arr[0,vv],TheRH_arr[0,vv],r"$ {} $".format(Pointees[vv]))
	    scats=ax1.scatter(TheT_arr[0,vv],Theq_arr[0,vv],c=TheRH_arr[0,vv],marker=r"$ {} $".format(Pointees[vv]),s=250,cmap=cmap,norm=norm, edgecolors='none' ) # s=1

        ax1.set_xlabel(xtitlee,size=14)
        ax1.set_ylabel(ytitlee,size=14)
        ax1.tick_params(axis='both', which='major', labelsize=14)
        
	cbax=fig.add_axes([0.86,0.1,0.03,0.8])
        cb=plt.colorbar(scats,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=14) 
        plt.figtext(0.97,0.5,'RH Anomalies (%rh)',size=14,ha='center',rotation='vertical',va='center')

    # add watermark and plot labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)

#    plt.figtext(0.02,0.96,TheLetter,size=18)
        if (TheReg == 'G'):
	    PointTitle=0
        if (TheReg == 'N'):
	    PointTitle=1
        if (TheReg == 'T'):
	    PointTitle=2
        if (TheReg == 'S'):
	    PointTitle=3
        ax1.set_title(titlees[PointTitle],size=18)

    else:

        # Four plot scenario
        fig,ax=plt.subplots(4,figsize=(8,8))	#6,18
        plt.clf()	# needs to be called after figure!!! (create the figure, then clear the plot space)
        TheLetter=['a)','b)','c)','d)']
	
	xstart=[0.1,0.48,0.1,0.48]
	xwide=0.36
	ystart=[0.54,0.54,0.08,0.08]
	ytall=0.36

	for pp in range(4):
            ax[pp]=plt.axes([xstart[pp],ystart[pp],xwide,ytall]) # left, bottom, width, height

            ax[pp].set_xlim([tmin,tmax])   
            ax[pp].set_ylim([qmin,qmax])   

            # make blank plot with zero lines on
	    ax[pp].plot(np.zeros(100),np.linspace(qmin,qmax,100),color='black',linewidth=2)
	    ax[pp].plot(np.linspace(tmin,tmax,100),np.zeros(100),color='black',linewidth=2)

            # plot black dots for the goods
            for vv in range(Thentims):
	        scats=ax[pp].scatter(TheT_arr[pp,vv],Theq_arr[pp,vv],c=TheRH_arr[pp,vv],marker=r"$ {} $".format(Pointees[vv]),s=200,cmap=cmap,norm=norm, edgecolors='none' ) # s=1

            if (pp == 2) | (pp == 3):
	        ax[pp].set_xlabel(xtitlee,size=12)
            if (pp == 0) | (pp == 2):
                ax[pp].set_ylabel(ytitlee,size=12)
	    if (pp == 0) | (pp == 1):	
		ax[pp].xaxis.set_ticklabels([])
	    if (pp == 1) | (pp == 3):	
		ax[pp].yaxis.set_ticklabels([])
            ax[pp].tick_params(axis='both', which='major', labelsize=12)

            plt.figtext((xstart[pp]+0.02),ystart[pp]+ytall-0.05,TheLetter[pp],size=14)
        
            ax[pp].set_title(titlees[pp],size=14)
        
	cbax=fig.add_axes([0.86,0.1,0.03,0.8])
        cb=plt.colorbar(scats,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=12) 
        plt.figtext(0.97,0.5,'RH Anomalies (%rh)',size=12,ha='center',rotation='vertical',va='center')

    # add watermark and plot labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


    #plt.show()
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")

#    raw_input("stop")	# REALLY USEFUL TO INTERACT WITHIN SUBROUTINE ctrl C
    # plt.ion()
    # plt.show() can then zoom and save
     
    return #PlotNiceDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# Read in region data for each variable
if (Region == 'A'):
   nReg=4
else:
    nReg=1   
tmpq_arr=np.empty((nReg,nmons))
tmpRH_arr=np.empty((nReg,nmons))
tmpT_arr=np.empty((nReg,nmons))
q_arr=np.empty((nReg,ntims))
RH_arr=np.empty((nReg,ntims))
T_arr=np.empty((nReg,ntims))

MyFile=INDIR+In_q
if (Region == 'A'):
    ReadInfo=['glob_q_anoms','nhem_q_anoms','trop_q_anoms','shem_q_anoms']
elif (Region == 'G'):
    ReadInfo=['glob_q_anoms']
elif (Region == 'N'):
    ReadInfo=['nhem_q_anoms']
elif (Region == 'T'):
    ReadInfo=['trop_q_anoms']
elif (Region == 'S'):
    ReadInfo=['shem_q_anoms']
tmpq_arr=ReadNetCDFTS(MyFile,ReadInfo,tmpq_arr)

MyFile=INDIR+In_RH
if (Region == 'A'):
    ReadInfo=['glob_RH_anoms','nhem_RH_anoms','trop_RH_anoms','shem_RH_anoms']
elif (Region == 'G'):
    ReadInfo=['glob_RH_anoms']
elif (Region == 'N'):
    ReadInfo=['nhem_RH_anoms']
elif (Region == 'T'):
    ReadInfo=['trop_RH_anoms']
elif (Region == 'S'):
    ReadInfo=['shem_RH_anoms']
tmpRH_arr=ReadNetCDFTS(MyFile,ReadInfo,tmpRH_arr)

MyFile=INDIR+In_T
if (Region == 'A'):
    ReadInfo=['glob_T_anoms','nhem_T_anoms','trop_T_anoms','shem_T_anoms']
elif (Region == 'G'):
    ReadInfo=['glob_T_anoms']
elif (Region == 'N'):
    ReadInfo=['nhem_T_anoms']
elif (Region == 'T'):
    ReadInfo=['trop_T_anoms']
elif (Region == 'S'):
    ReadInfo=['shem_T_anoms']
tmpT_arr=ReadNetCDFTS(MyFile,ReadInfo,tmpT_arr)

# If annual - convert monthly mean anomalies to annual mean anomalies
# THERE SHOULD BE NO MISSING DATA IN THESE!!!!
for rr in range(nReg):
    if (TimeRes == 'Y'):
#        pdb.set_trace()
	q_arr[rr,:]=np.mean(np.reshape(tmpq_arr[rr,:],(ntims,12)),axis=1)
        RH_arr[rr,:]=np.mean(np.reshape(tmpRH_arr[rr,:],(ntims,12)),axis=1)
        T_arr[rr,:]=np.mean(np.reshape(tmpT_arr[rr,:],(ntims,12)),axis=1)
    else:
        q_arr[rr,:]=tmpq_arr[rr,:]
        RH_arr[rr,:]=tmpRH_arr[rr,:]
        T_arr[rr,:]=tmpT_arr[rr,:]
    
# Plot the scatter
MyFile=OUTDIR+OutPlot
PlotScatter(MyFile,YrStr,ntims,q_arr,RH_arr,T_arr,Region)
		
#    stop()

print("And, we are done!")

