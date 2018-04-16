#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 8th October 2015
# Last update: 16th April 2018
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
# import scipy.stats as ss # for pearsonr
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
# Version 3 16 April 2018
# ---------
#  
# Enhancements
# Updated editable info so fewer edits are required to run for the most recent version/year
#  
# Changes
#  
# Bug fixes
#
# Version 2 9 August 2016
# ---------
#  
# Enhancements
# Can also plot T vs RH coloured by q anomaly
#  
# Changes
#  
# Bug fixes
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
import scipy.stats as ss # for pearsonr
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
thenmon='MAR'
thenyear='2018'
nowmon='APR'
nowyear='2018'
version='4.0.0.2017f'

styr=1973
edyr=2017
nyrs=(edyr-styr)+1
nmons=(nyrs)*12
if (TimeRes == 'Y'):
    ntims=nyrs
else:
    ntims=nmons    
YrStr=np.array(range(styr,edyr+1),dtype=str)
YrStr=np.array(([i[2:5] for i in YrStr])) # now a string array of the last two digits

# Set up directories and files
INDIR='/data/local/hadkw/HADCRUH2/UPDATE'+str(edyr)+'/STATISTICS/TIMESERIES/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE'+str(edyr)+'/IMAGES/ANALYSIS/'

In_q='HadISDH.landq.'+version+'_FLATgridIDPHA5by5_anoms8110_'+thenmon+thenyear+'_areaTS_1973'+str(edyr)+'.nc'
In_RH='HadISDH.landRH.'+version+'_FLATgridIDPHA5by5_anoms8110_'+thenmon+thenyear+'_areaTS_1973'+str(edyr)+'.nc'
In_T='HadISDH.landT.'+version+'_FLATgridIDPHA5by5_anoms8110_'+thenmon+thenyear+'_areaTS_1973'+str(edyr)+'.nc'

OutPlotTq='ScatterTqbyRH_HadISDH.'+version+'_'+Region+'_'+nowmon+nowyear
OutPlotTRH='ScatterTRHbyq_HadISDH.'+version+'_'+Region+'_'+nowmon+nowyear
 
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
# MakeUpSteps
def MakeUpSteps(TheArray,stepsies=9):
    ''' Given a max and min, make up NICE step sizes for a 9 element colourbar '''
    ''' Currently works with a minimum range of 0.2 and a maximum or 3.0 '''
    ''' Can only deal with symmetric ranges '''
    ''' READS: 	TheArray - an array of data '''
    '''		stepsies (OPTIONAL) - number of colours in colourbar - default 9 is NICE '''
    ''' RETURNS: vmin - minimum threshold of range '''
    '''		 vmax - maximum threshold of range '''
    '''		 bounds - stepsies linear increments through the range from vmin to vmax '''
    ''' 	 strcounds - strings of the bounds for labelling the colourbar '''

    vmax=np.int(np.ceil(np.max(abs(TheArray))*10))/10.    
    vmin=-vmax
    nsteps = stepsies
    if (vmax <= 0.2):
        vmax = 0.2
	vmin = -0.2
    if (vmax <= 0.3):
        vmax = 0.32
	vmin = -0.32
    elif (vmax <= 0.4):
        vmax = 0.4
	vmin = -0.4
    elif (vmax <= 0.6):
        vmax = 0.6
	vmin = -0.6
    elif (vmax <= 0.8):
        vmax = 0.8
	vmin = -0.8
    elif (vmax <= 1.0):
        vmax = 1.0
	vmin = -1.0
    elif (vmax <= 1.2):
        vmax = 1.2
	vmin = -1.2
    elif (vmax <= 1.6):
        vmax = 1.6
	vmin = -1.6
    elif (vmax <= 2.0):
        vmax = 2.0
	vmin = -2.0
    elif (vmax <= 3.0):
        vmax = 3.0
	vmin = -3.0
	#    pdb.set_trace() # stop here and play

    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
        
    return vmax,vmin,strbounds,bounds

#************************************************************************
# PlotScatter
def PlotScatter(TheFileTq,TheFileTRH,TheYrStr,Thentims,Theq_arr,TheRH_arr,TheT_arr,TheReg):
    ''' Plot Tq scatter with colours related to RH'''
    ''' Plot TRH scatter with colours related to q'''
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
    
    # FIRST MAKE UP THE TqbyRH plot

    # Call MakeUpSteps routine to get a NICE set of colourbar indices
    vmin,vmax,strbounds,bounds=MakeUpSteps(TheRH_arr)
    
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
	
	# plot 1:1 line dashed
	ax1.plot(np.linspace(-5,5,100),np.linspace(-5,5,100),color='black',linewidth=2,linestyle='dashed')

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

        # Get correlation and slope of scatter and add to plot
	#pcorr = ss.pearsonr(TheT_arr[0,:],Theq_arr[0,:]) # element 0 = pearson correlation coefficient, element 1 = two-tailed p-value 
	linvals = ss.linregress(TheT_arr[0,:],Theq_arr[0,:]) # 0 = slope, 1 = intercept, 2 = r-value, 3 = two-tailed p-value, 4 = sterr
        plt.figtext(0.05,0.96,'r = '+"{:3.2f}".format(linvals[2]),size=12)
        plt.figtext(0.05,0.9,'m = '+"{:3.2f}".format(linvals[0]),size=12)
        plt.figtext(0.05,0.84,'p = '+"{:1.2f}".format(linvals[3]),size=12)


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

	    # plot 1:1 line dashed
	    ax[pp].plot(np.linspace(-5,5,100),np.linspace(-5,5,100),color='black',linewidth=2,linestyle='dashed')

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

            # Get correlation and slope of scatter and add to plot
	    #pcorr = ss.pearsonr(TheT_arr[pp,:],Theq_arr[pp,:]) # element 0 = pearson correlation coefficient, element 1 = two-tailed p-value 
	    linvals = ss.linregress(TheT_arr[pp,:],Theq_arr[pp,:]) # 0 = slope, 1 = intercept, 2 = r-value, 3 = two-tailed p-value, 4 = sterr
            plt.figtext((xstart[pp]+0.05),ystart[pp]+ytall-0.05,'r = '+"{:3.2f}".format(linvals[2]),size=12)
            plt.figtext((xstart[pp]+0.05),ystart[pp]+ytall-0.07,'m = '+"{:3.2f}".format(linvals[0]),size=12)
            plt.figtext((xstart[pp]+0.05),ystart[pp]+ytall-0.09,'p = '+"{:1.2f}".format(linvals[3]),size=12)
        
	cbax=fig.add_axes([0.86,0.1,0.03,0.8])
        cb=plt.colorbar(scats,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=12) 
        plt.figtext(0.97,0.5,'RH Anomalies (%rh)',size=12,ha='center',rotation='vertical',va='center')

    # add watermark and plot labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


    #plt.show()
    plt.savefig(TheFileTq+".eps")
    plt.savefig(TheFileTq+".png")

#    raw_input("stop")	# REALLY USEFUL TO INTERACT WITHIN SUBROUTINE ctrl C
    # plt.ion()
    # plt.show() can then zoom and save

    #***********************************
    # SECOND MAKE UP THE TRHbyq plot

    # Call MakeUpSteps routine to get a NICE set of colourbar indices
    vmin,vmax,strbounds,bounds=MakeUpSteps(Theq_arr)
    
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)
        
    ytitlee='Relative Humidity Anomalies (%rh)'
    xtitlee='Temperature Anomalies ($^{o}$C)'
    titlees=['Globe 70$^{o}$S to 70$^{o}$N','N. Hemisphere 20$^{o}$N to 70$^{o}$N','Tropics 20$^{o}$S to 20$^{o}$N','S. Hemisphere 70$^{o}$S to 20$^{o}$S']

    # set up max and min of RH and T for axes - keep same for all regions
    rhmax=np.ceil(np.max(abs(TheRH_arr))/0.1)*0.1  
    rhmin=-rhmax
    tmax=np.ceil(np.max(abs(TheT_arr))/0.1)*0.1  
    tmin=-tmax
                
    # set up plot - are we working with one region or four?
    if (TheReg != 'A'):

        # Single plot scenario
        fig = plt.figure(1,figsize=(8,8))
        plt.clf()	# needs to be called after figure!!! (create the figure, then clear the plot space)

        ax1=plt.axes([0.1,0.1,0.75,0.8]) # left, bottom, width, height

        ax1.set_xlim([tmin,tmax])   
        ax1.set_ylim([rhmin,rhmax])   

        # make blank plot with zero lines on
	ax1.plot(np.zeros(100),np.linspace(rhmin,rhmax,100),color='black',linewidth=2)
	ax1.plot(np.linspace(tmin,tmax,100),np.zeros(100),color='black',linewidth=2)

	# plot 1:1 line dashed
	ax1.plot(np.linspace(-5,5,100),np.linspace(-5,5,100),color='black',linewidth=2,linestyle='dashed')

        # plot YEAR LABELS for the goods
        for vv in range(Thentims):
	    #print(vv,TheT_arr[0,vv],Theq_arr[0,vv],TheRH_arr[0,vv],r"$ {} $".format(Pointees[vv]))
	    scats=ax1.scatter(TheT_arr[0,vv],TheRH_arr[0,vv],c=Theq_arr[0,vv],marker=r"$ {} $".format(Pointees[vv]),s=250,cmap=cmap,norm=norm, edgecolors='none' ) # s=1

        ax1.set_xlabel(xtitlee,size=14)
        ax1.set_ylabel(ytitlee,size=14)
        ax1.tick_params(axis='both', which='major', labelsize=14)
        
	cbax=fig.add_axes([0.86,0.1,0.03,0.8])
        cb=plt.colorbar(scats,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=14) 
        plt.figtext(0.97,0.5,'q Anomalies (g kg$^{-1}$)',size=14,ha='center',rotation='vertical',va='center')

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

        # Get correlation and slope of scatter and add to plot
	#pcorr = ss.pearsonr(TheT_arr[0,:],TheRH_arr[0,:]) # element 0 = pearson correlation coefficient, element 1 = two-tailed p-value 
	linvals = ss.linregress(TheT_arr[0,:],TheRH_arr[0,:]) # 0 = slope, 1 = intercept, 2 = r-value, 3 = two-tailed p-value, 4 = sterr
        plt.figtext(0.05,0.96,'r = '+"{:3.2f}".format(linvals[2]),size=12)
        plt.figtext(0.05,0.9,'m = '+"{:3.2f}".format(linvals[0]),size=12)
        plt.figtext(0.05,0.84,'p = '+"{:1.2f}".format(linvals[3]),size=12)

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
            ax[pp].set_ylim([rhmin,rhmax])   

            # make blank plot with zero lines on
	    ax[pp].plot(np.zeros(100),np.linspace(rhmin,rhmax,100),color='black',linewidth=2)
	    ax[pp].plot(np.linspace(tmin,tmax,100),np.zeros(100),color='black',linewidth=2)

	    # plot 1:1 line dashed
	    ax[pp].plot(np.linspace(-5,5,100),np.linspace(-5,5,100),color='black',linewidth=2,linestyle='dashed')

            # plot black dots for the goods
            for vv in range(Thentims):
	        scats=ax[pp].scatter(TheT_arr[pp,vv],TheRH_arr[pp,vv],c=Theq_arr[pp,vv],marker=r"$ {} $".format(Pointees[vv]),s=200,cmap=cmap,norm=norm, edgecolors='none' ) # s=1

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

            # Get correlation and slope of scatter and add to plot
	    #pcorr = ss.pearsonr(TheT_arr[pp,:],TheRH_arr[pp,:]) # element 0 = pearson correlation coefficient, element 1 = two-tailed p-value 
	    linvals = ss.linregress(TheT_arr[pp,:],TheRH_arr[pp,:]) # 0 = slope, 1 = intercept, 2 = r-value, 3 = two-tailed p-value, 4 = sterr
            plt.figtext((xstart[pp]+0.05),ystart[pp]+ytall-0.05,'r = '+"{:3.2f}".format(linvals[2]),size=12)
            plt.figtext((xstart[pp]+0.05),ystart[pp]+ytall-0.07,'m = '+"{:3.2f}".format(linvals[0]),size=12)
            plt.figtext((xstart[pp]+0.05),ystart[pp]+ytall-0.09,'p = '+"{:1.2f}".format(linvals[3]),size=12)
        
	cbax=fig.add_axes([0.86,0.1,0.03,0.8])
        cb=plt.colorbar(scats,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=12) 
        plt.figtext(0.97,0.5,'q Anomalies (g kg$^{-1}$)',size=12,ha='center',rotation='vertical',va='center')

    # add watermark and plot labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


    #plt.show()
    plt.savefig(TheFileTRH+".eps")
    plt.savefig(TheFileTRH+".png")

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
MyFileTq=OUTDIR+OutPlotTq
MyFileTRH=OUTDIR+OutPlotTRH
PlotScatter(MyFileTq,MyFileTRH,YrStr,ntims,q_arr,RH_arr,T_arr,Region)
		
#    stop()

print("And, we are done!")

