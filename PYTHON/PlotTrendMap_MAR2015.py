
#***************************************
# 22nd April 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# 22 Aprl 2015 KMW - v1
# Reads in decadal trend map
# Works out the mean trend for each latitude (simple average)  
# Plots the trend for each gridbox, grixbox boundary IF trend is sig at 95%
# Also adds a vertical latitude by trend figure - scatter of gridbox trends for each latitude and latitude average
# Fills 0-1 with dark grey to represent actual land fraction in terms of gridboxs at that latitude band
# Overlays light grey to represent fraction of land gridboxes actually observed

# now also counts % of land gridboxes present for globe (70S to 70N), N Hemi (20N to 70N), Tropics (20S to 20N) and S Hemi (70S to 20S)
# This isn't perfect because the land/sea mask may not contain all small islands. It does contain some at least. CRUTEM has a lot of these
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotTrendMap_MAR2015.py
#
# REQUIRES
# 

#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 22 April 2015
# Last update: 5 October 2015
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	# this will probably change
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Reads in decadal trend map netCDF file (any 5 by 5 degree grid)
# Works out the mean trend for each latitude (simple average)  
# Plots the trend for each gridbox, grixbox boundary IF trend is sig at 95% (STILL NEED TO SORT OUT HOW THIS IS CALCULATED!!!)
# Also adds a vertical latitude by trend figure - scatter of gridbox trends for each latitude and latitude average
# Fills 0-1 with dark grey to represent actual land fraction in terms of gridboxes at that latitude band
# Overlays light grey to represent fraction of land gridboxes actually observed as % of dark grey bar
# Counts % of land gridboxes present for globe (70S to 70N), N Hemi (20N to 70N), Tropics (20S to 20N) and S Hemi (70S to 20S)
# This isn't perfect because the land/sea mask may not contain all small islands. It does contain some at least. CRUTEM has a lot of these
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
# ReadNetCDFGrid - infile function to read in netCDF grid, written by Kate Willett
# PlotTrendMap - infile function to plot the map, written by Kate Willett
# 
# -----------------------
# DATA
# -----------------------
# directory for trendmaps:
# /data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TRENDS/
# directory for land cover map:
# /data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/
# land cover file used to calculate % land represented:
# HadCRUT.4.3.0.0.land_fraction.nc
# files currently worked on:
# HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014.nc
# HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014.nc
# HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_MPtrends_19732014.nc
# HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014.nc
# BERKELEY_T_5by519762005clim_anoms_19732014_MPtrends_19732014.nc
# GISS_T_5by519762005clim_anoms_19732014_MPtrends_19732014.nc
# CRUTEM.4.3.0.0.anomalies_MPtrends_19732014.nc
# GHCNM_18802014_MPtrends_19732014.nc'
# HadISDH.lande.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014.nc
# HadISDH.landTw.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014.nc
# HadISDH.landTd.2.0.1.2014p_FLATgridPHADPD5by5_JAN2015_MPtrends_19732014.nc
# HadISDH.landDPD.2.0.1.2014p_FLATgridPHA5by5_JAN2015_MPtrends_19732014.nc
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# Select/add 'candidate' input file (or add a new one)
# Check filenames/directory paths are correct
# Select/add 'OUTPLOT' output file
# Select/add 'Unit' variable unit
# Select/add 'Letty' desired letters to apply as labels to plots
# Select/add 'Namey' Data product name / plot title or make a blank
#
# run:
# python2.7 PlotTrendMap_MAR2015.py
# 
# -----------------------
# OUTPUT
# -----------------------
# directory for output images:
# /data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/
# Output image files (nowmon+nowyear= e.g., OCT2015):
# TrendMap_HadISDH.landq.2.0.1.2014p_'+nowmon+nowyear+
# TrendMap_HadISDH.landRH.2.0.1.2014p_'+nowmon+nowyear
# TrendMap_HadISDH.landT.2.0.1.2014p_RAW_'+nowmon+nowyear
# TrendMap_HadISDH.landT.2.0.1.2014p_'+nowmon+nowyear
# TrendMap_BERKELEY_'+nowmon+nowyear
# TrendMap_GISS_'+nowmon+nowyear
# TrendMap_CRUTEM4.3.0.0_'+nowmon+nowyear
# TrendMap_GHCNM3_'+nowmon+nowyear
# TrendMap_HadISDH.lande.2.0.1.2014p_'+nowmon+nowyear+
# TrendMap_HadISDH.landTw.2.0.1.2014p_'+nowmon+nowyear+
# TrendMap_HadISDH.landTd.2.0.1.2014p_'+nowmon+nowyear+
# TrendMap_HadISDH.landDPD.2.0.1.2014p_'+nowmon+nowyear+
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 2 (5 October 2015)
# ---------
#  
# Enhancements
# Set up is now done in bundles - easier to make sure you have all the correct file paths and info
# Read in info for netCDF files is now part of set up rather than clumsily coded into the ReadNetCDFGrid function
# which means you don't have to edit functions for each different file
# Now sorts out nice colour bar index depending on data
#  
# Changes
#  
# Bug fixes
#
# 
# Version 1 22 April 2015)
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
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

# Generic things:******************************************

# Missing data
mdi=-1e30 # may set up as masked arrays later

# Output date stamp
nowmon='OCT'
nowyear='2015'

# Set up directories and files
INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TRENDS/'
#INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'	# may be needed for non-HadISDH variables
INDIRO='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/'

# Land cover file:
#incover='new_coverpercentjul08'
incover='HadCRUT.4.3.0.0.land_fraction'

# Variables
CandTrends=[]
CompTrends=[]
RatTrends=[]
LatList=[]
LonList=[]
Letty=['a)','b)']


#**********************************************************
# Run choice bundle for input/output files, units, names, letters, read in varnames, colourmap
# CHOOSE/ADD A BUNDLE!!!

#candidate='HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.landq.2.0.1.2014p_'+nowmon+nowyear
#Unit='g kg$^{-1}$'  #'degrees C'
#Namey='HadISDH.landq.2.0.1.2014p decadal trends'
#nlats=36	       #set once file read in
#nlons=72	       #set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['q_MPtrend','q_MP5th','q_MP95th'])
#ColourMapChoice=('BrBG','noflip')

#candidate='HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.landRH.2.0.1.2014p_'+nowmon+nowyear
#Unit='%rh'  
#Namey='HadISDH.landRH.2.0.1.2014p decadal trends'
#nlats=36	      #set once file read in
#nlons=72	      #set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['RH_MPtrend','RH_MP5th','RH_MP95th'])
#ColourMapChoice=('BrBG','noflip')

#candidate='HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.landT.2.0.1.2014p_RAW_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='HadISDH.landT.2.0.1.2014p RAW decadal trends'
#nlats=36		#set once file read in
#nlons=72		#set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['T_MPtrend','T_MP5th','T_MP95th'])
#ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'

#candidate='HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.landT.2.0.1.2014p_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='HadISDH.landT.2.0.1.2014p decadal trends'
#nlats=36	       #set once file read in
#nlons=72	       #set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['T_MPtrend','T_MP5th','T_MP95th'])
#ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'

#candidate='BERKELEY_T_5by519762005clim_anoms_19732014_MPtrends_19732014'
#OUTPLOT='TrendMap_BERKELEY_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='BERKELEY decadal trends'
#nlats=36		#set once file read in
#nlons=72		#set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['T_MPtrend','T_MP5th','T_MP95th'])
#ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'

#candidate='GISS_T_5by519762005clim_anoms_19732014_MPtrends_19732014'
#OUTPLOT='TrendMap_GISS_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='GISS decadal trends'
#nlats=36		#set once file read in
#nlons=72		#set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['T_MPtrend','T_MP5th','T_MP95th'])
#ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'

#candidate='CRUTEM.4.3.0.0.anomalies_MPtrends_19732014'
#OUTPLOT='TrendMap_CRUTEM4.3.0.0_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='CRUTEM4.3.0.0 decadal trends'
#nlats=36		#set once file read in
#nlons=72		#set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['T_MPtrend','T_MP5th','T_MP95th'])
#ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'

#candidate='GHCNM_18802014_MPtrends_19732014'
#OUTPLOT='TrendMap_GHCNM3_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='GHCNM3 decadal trends'
#nlats=36		#set once file read in
#nlons=72		#set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['T_MPtrend','T_MP5th','T_MP95th'])
#ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'

#candidate='HadISDH.lande.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.lande.2.0.1.2014p_'+nowmon+nowyear
#Unit='hPa'  #'degrees C'
#Namey='HadISDH.lande.2.0.1.2014p decadal trends'
#nlats=36	       #set once file read in
#nlons=72	       #set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['e_MPtrend','e_MP5th','e_MP95th'])
#ColourMapChoice=('BrBG','noflip')

#candidate='HadISDH.landTw.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.landTw.2.0.1.2014p_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='HadISDH.landTw.2.0.1.2014p decadal trends'
#nlats=36	      #set once file read in
#nlons=72	      #set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['Tw_MPtrend','Tw_MP5th','Tw_MP95th'])
#ColourMapChoice=('BrBG','noflip')

#candidate='HadISDH.landTd.2.0.1.2014p_FLATgridPHADPD5by5_JAN2015_MPtrends_19732014'
#OUTPLOT='TrendMap_HadISDH.landTd.2.0.1.2014p_'+nowmon+nowyear
#Unit='$^{o}$C'  #'degrees C'
#Namey='HadISDH.landTd.2.0.1.2014p decadal trends'
#nlats=36	      #set once file read in
#nlons=72	      #set once file read in
#LatInfo=list(['latitude',nlats,-87.5])
#LonInfo=list(['longitude',nlons,-177.5])
#ReadInfo=list(['Td_MPtrend','Td_MP5th','Td_MP95th'])
#ColourMapChoice=('BrBG','noflip')

candidate='HadISDH.landDPD.2.0.1.2014p_FLATgridPHA5by5_JAN2015_MPtrends_19732014'
OUTPLOT='TrendMap_HadISDH.landDPD.2.0.1.2014p_'+nowmon+nowyear
Unit='$^{o}$C'  #'degrees C'
Namey='HadISDH.landDPD.2.0.1.2014p decadal trends'
nlats=36	      #set once file read in
nlons=72	      #set once file read in
LatInfo=list(['latitude',nlats,-87.5])
LonInfo=list(['longitude',nlons,-177.5])
ReadInfo=list(['DPD_MPtrend','DPD_MP5th','DPD_MP95th'])
ColourMapChoice=('BrBG','flip')

#************************************************************************
# Subroutines
#************************************************************************
# READNETCDFGRID
def ReadNetCDFGrid(FileName,ReadInfo,LatInfo,LonInfo):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data 
	FileName: stroing containing filepath/name
	ReadInfo: list of strings of variable names for the trend, 5th percentile, 95th percentile
	LatInfo: list of a string for the latitude variable name, an integer for number of lats, a float for the start latitude
	LonInfo: list of a string for the longitude variable name, an integer for number of lons, a float for the start longitude '''

    ncf=netcdf.netcdf_file(FileName,'r')

    # ncf.variables this lists the variable names
    #var=f.variables['latitude']
    #TheLatList=var.data
    # lats currently screwy so make up
    if (LatInfo[2] < 0):
        TheLatList=np.arange(LatInfo[2], LatInfo[2]+180.,(180./LatInfo[1]))
    else:
        TheLatList=np.arange(LatInfo[2], LatInfo[2]-180.,-(180./LatInfo[1]))    
    #var=f.variables['longitude']
    #TheLonList=var.data
    # lons currently screwy so make up
    TheLonList=np.arange(-177.5,182.5,5.)
    if (LonInfo[2] < 10):
        TheLonList=np.arange(LonInfo[2], LonInfo[2]+360.,(360./LonInfo[1]))
    else:
        TheLonList=np.arange(LonInfo[2], LonInfo[2]-360.,-(360./LonInfo[1]))    

    var=ncf.variables[ReadInfo[0]]
    TheData=np.array(var.data)
    var=ncf.variables[ReadInfo[1]]
    TheLower=np.array(var.data)
    var=ncf.variables[ReadInfo[2]]
    TheUpper=np.array(var.data)

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheUpper,TheLower,TheLatList,TheLonList # ReadNetCDFGrid

#************************************************************************
# PlotTrendMap
def PlotTrendMap(TheFile,LandCover,TheLatList,TheLonList,TheCandData,TheCandUppers,TheCandLowers,
                    TheUnitee,TheLetter,TheNamee,ColourMapChoice):
    ''' Create a masked array of trends
        Create a masked array of significant trends
        Plot trends on map, bounded if significant
        Add vertical latitude/trend scatter with average trend overlaid
        Save as eps and png '''

    # Missing data
    mdi=-1e30 # may set up as masked arrays later
      
    # Create the masked array of trends
    MSKTheCandData=ma.masked_where(TheCandData == mdi,TheCandData)
    
    # Create the masked array of significant trends
    SigTrends=TheCandData.copy()
    #SigTrends.fill(mdi)
    SigTrends[np.where(((TheCandData > 0.) & (TheCandLowers <= 0.)) | ((TheCandData <= 0.) & (TheCandUppers >= 0.)))]=mdi
    MSKSigTrends=ma.masked_where(SigTrends == mdi,SigTrends)
    #stop
    # Get Latitude Average Trend
    LatTrend=np.zeros(len(TheLatList))
    LatTrend.fill(mdi)
    # For each latitude find average trend
    for ltt in range(len(TheLatList)):
        gots=np.where(TheCandData[ltt,:] != mdi)[0]
	if (len(gots) > 0):
	    LatTrend[ltt]=np.mean(TheCandData[ltt,np.where(TheCandData[ltt,:] != mdi)])

    # make 2d arrays of lats and lons
    # nudge -2.5 degrees to make them south/west gridbox corners, not centres
    # add extra row/column to bound the data
    ArrLons,ArrLats=np.meshgrid(TheLonList,TheLatList)
    LngArrLons,LngArrLats=np.meshgrid(np.append(TheLonList-2.5,180.),np.append(TheLatList-2.5,90.))
    
    # set up plot
    plt.clf()
    fig=plt.figure(figsize=(10,5))
    #plt.figure(1,figsize=(10,3))
    plt1=plt.axes([0.01,0.05,0.64,0.9]) # left, bottom, width, height
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))

    # make up a blue to red (reverse) colour map
    cmap=plt.get_cmap(ColourMapChoice[0])
    
#    fullcols=list()
#    [fullcols.append(bluelist[i]) for i in range(len(bluelist))]
#    [fullcols.append(redlist[i]) for i in range(len(redlist))]
    
    cmaplist=[cmap(i) for i in range(cmap.N)]
    for loo in range((cmap.N/2)-10,(cmap.N/2)+10):
        cmaplist.remove(cmaplist[(cmap.N/2)-10]) # remove the very pale colours in the middle
    if (ColourMapChoice[1] == 'flip'):	# then reverse the colours
        cmaplist.reverse()
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    
    # work out best max and min values for colourbar and try to make them 'nice
    # must be an odd number of steps
    # for less than 0.1 must be 14 or fewer steps
    vmax=np.int(np.ceil(np.max(abs(MSKTheCandData))*10))/10.    
    vmin=-vmax
    if (vmax <= 0.3):
        nsteps=np.int((vmax-vmin)/0.05)+1
    elif (vmax <= 0.5):
        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.06)*0.06  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.06)+1
    elif (vmax <= 1.0):
        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.1)*0.1  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.1)+1
    elif (vmax <= 1.4):
        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.2)*0.2  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.2)+1
    elif (vmax > 1.4):
        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.3)*0.3  
        vmin=-vmax
        nsteps=np.int((vmax-vmin)/0.3)+1
	#    pdb.set_trace() # stop here and play
    
    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCandData,cmap=cmap,norm=norm,latlon='TRUE')

    grids=m.pcolor(LngArrLons,LngArrLats,MSKSigTrends,cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5,latlon='TRUE')

    cbax=fig.add_axes([0.03,0.08,0.59,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=10) 
    plt.figtext(0.31,0.01,'Anomalies ('+TheUnitee+')',size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.9,TheLetter[0],size=16)
    plt.figtext(0.5,0.95,TheNamee,size=16,ha='center')
    
#    print(CountGoods,CountLargeNegs,CountSmallNegs,CountSmallPos,CountLargePos)
#    pctLNs=round((float(CountLargeNegs)/CountGoods)*100.,1)
#    pctSNs=round((float(CountSmallNegs)/CountGoods)*100.,1)
#    pctVSs=round((float(CountVSmalls)/CountGoods)*100.,1)
#    pctSPs=round((float(CountSmallPos)/CountGoods)*100.,1)
#    pctLPs=round((float(CountLargePos)/CountGoods)*100.,1)
    
#    ax1=plt.axes([0.01,0.01,0.64,0.09],frameon=False) # map only
#    ax1.set_ylim(0,1)
#    ax1.set_xlim(0,1)
#    ax1.axes.get_yaxis().set_visible(False)
#    ax1.axes.get_xaxis().set_visible(False)
#    plt.annotate(str(pctLNs)+"% ratio < -1",xy=(0.07,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.05,0.8,markersize=10,marker='s',color='Firebrick')
#    plt.annotate(str(pctSNs)+"% ratio -1 to -0.05",xy=(0.07,0.1),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.05,0.2,markersize=10,marker='s',color='LightCoral')
#    plt.annotate(str(pctVSs)+"% small trend",xy=(0.40,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.38,0.8,markersize=10,marker='s',color='Khaki')
#    plt.annotate(str(pctSPs)+"% ratio 0.05 to 1",xy=(0.73,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.71,0.8,markersize=10,marker='s',color='LightSkyBlue')
#    plt.annotate(str(pctLPs)+"% ratio > 1",xy=(0.73,0.1),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.71,0.2,markersize=10,marker='s',color='DodgerBlue')

    # Attempt to plot the latitude/ratio scatter
    ax2=plt.axes([0.73,0.12,0.24,0.76]) # map only
    ax2.set_ylim(-90,90)
    #ax2.set_ylabel('Latitude')
    ax2.set_yticks(list([-60,-30,0,30,60]))
    ax2.set_yticklabels(list(['-60$^{o}$N','-30$^{o}$S','Equator','30$^{o}$N','60$^{o}$N'])) #,rotation=90)
    ax2.set_xlabel('Decadal Trend ('+TheUnitee+')')
    minx=np.min(TheCandData[np.where(TheCandData != mdi)])-0.1
    maxx=np.max(TheCandData[np.where(TheCandData != mdi)])+0.1
    ax2.set_xlim(minx,maxx)
    ax2.tick_params(axis='both',labelsize=10)
    # background fill the scatter plot with % land present (light grey) and % land with data (dark gre        
    for ltt in range(len(TheLatList)):
        # First merge LandCover with Data coverage so that at least some small islands are counted initially
	LandCover[ltt,np.where(TheCandData[ltt,:] > mdi)[0]]=100.
        landcount=len(np.where(LandCover[ltt,:] > 0)[0])
	if (landcount > 0):
	    pctland=float(landcount)/72.
	else:
	    pctland=0
	landpresent=len(np.where(TheCandData[ltt,:] > mdi)[0])
#	if (landpresent > landcount):	# small islands not picked up in landcover file so have to add in those that have been observed - not ideal
#	    landcount=landpresent
#	    pctland=float(landcount)/72.
	if (pctland > 0) & (landpresent > 0):
	    pctdata=landpresent/float(landcount)
        else:
	    pctdata=0
	print(TheLatList[ltt],pctland,pctdata)
#	stop()
	if (pctland > 0): # fill row with appropriate amount of light grey
	    newxmin=0
	    newxmax=maxx*pctland
	    plt.fill_between(np.array([newxmin,newxmax]),TheLatList[ltt]-2.5,TheLatList[ltt]+2.5,facecolor='DimGrey',edgecolor='0.0',linewidth=0.0)
	if (pctdata > 0): # fill row with appropriate amount of dark grey
	    newxmin=0
	    newxmax=(maxx*pctland)*pctdata
	    plt.fill_between(np.array([newxmin,newxmax]),TheLatList[ltt]-2.5,TheLatList[ltt]+2.5,facecolor='Silver',edgecolor='0.0',linewidth=0.0)
    #ax.set_xmargin(0.01)
    plt.plot(np.zeros(2),np.array((-90,90)),color='black')
    plt.scatter(MSKTheCandData,ArrLats,c=MSKTheCandData,marker='o',cmap=cmap,norm=norm,edgecolor='0.0',linewidth=0.0)
    plt.plot(LatTrend[np.where(LatTrend != mdi)],TheLatList[np.where(LatTrend != mdi)],color='black',linewidth=2)

    plt.figtext(0.7,0.9,TheLetter[1],size=16)
    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotTrendRatMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in trend maps
MyFile=INDIRC+candidate+'.nc'
CandData,CandUpper,CandLower,LatList,LonList=ReadNetCDFGrid(MyFile,ReadInfo,LatInfo,LonInfo)

ncf=netcdf.netcdf_file(INDIRO+incover+'.nc','r')
#var=ncf.variables['pct_land']
var=ncf.variables['land_area_fraction']
PctLand=np.array(var.data)
#PctLand=np.transpose(PctLand)
PctLand=np.flipud(PctLand)
PctLand=PctLand[0,:,:]
ncf.close()

# pass to plotter
MyFile=OUTDIR+OUTPLOT
PlotTrendMap(MyFile,PctLand,LatList,LonList,CandData,CandUpper, CandLower,
                        Unit,Letty,Namey,ColourMapChoice)

# output pct of land boxes represented per region
for ltt in range(len(LatList)):
    # First merge LandCover with Data coverage so that at least some small islands are counted initially
    PctLand[ltt,np.where(CandData[ltt,:] > mdi)[0]]=100.
        
GlobCount=len(np.where(PctLand[4:31,:] > 0)[0])
AACount=len(np.where(PctLand[0:3,:] > 0)[0])
SHCount=len(np.where(PctLand[4:13,:] > 0)[0])
TropCount=len(np.where(PctLand[14:21,:] > 0)[0])
NHCount=len(np.where(PctLand[22:31,:] > 0)[0])
ACount=len(np.where(PctLand[32:35,:] > 0)[0])

ActGlobCount=len(np.where(CandData[4:31,:] != mdi)[0])
ActAACount=len(np.where(CandData[0:3,:] != mdi)[0])
ActSHCount=len(np.where(CandData[4:13,:] != mdi)[0])
ActTropCount=len(np.where(CandData[14:21,:] != mdi)[0])
ActNHCount=len(np.where(CandData[22:31,:] != mdi)[0])
ActACount=len(np.where(CandData[32:35,:] != mdi)[0])
	
print('GLOB: ',72*28,GlobCount,ActGlobCount,float(GlobCount)/(72*28),float(ActGlobCount)/GlobCount)	
print('ARCT: ',72*4,ACount,ActACount,float(ACount)/(72*4),float(ActACount)/ACount)	
print(' HEM: ',72*10,NHCount,ActNHCount,float(NHCount)/(72*10),float(ActNHCount)/NHCount)	
print('TROP: ',72*8,TropCount,ActTropCount,float(TropCount)/(72*8),float(ActTropCount)/TropCount)	
print('S HEM: ',72*10,SHCount,ActSHCount,float(SHCount)/(72*10),float(ActSHCount)/SHCount)	
print('AARCT: ',72*4,AACount,ActAACount,float(AACount)/(72*4),float(ActAACount)/AACount)	
	
		
#    stop()

print("And, we are done!")

