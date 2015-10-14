#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 22 June 2015
# Last update: 13 October 2015
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This reads in a netCDF file of a gridded dataset and plots a map
# It can renormalise data to the desired climatology period (and deal with actual monthly means)
# It can plot: 
#	a single month, 
#	an average of months within a year (or adjacent for DJF) up to annual - set minimum data presence
#	an average of single months across a period of years (climatology) - set minimum data presence
# 	an average of several months across a period of years (climatology) up to annual - set minimum data presence
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Inbuilt: (may not all be required actually)
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
# import pdb # pdb.set_trace() or c
# 
# -----------------------
# DATA
# -----------------------
# The code requires a monthly resolution gridded dataset (anomalies or monthly means) in netCDF format
# Currently it works with:
# candidate='HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
# candidate='HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
# candidate='HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_cf'
# candidate='HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
# candidate='BERKELEY_T_5by519762005clim_anoms_19732014'
# candidate='GISS_T_5by519762005clim_anoms_19732014'
# candidate='CRUTEM.4.3.0.0.anomalies'
# candidate='GHCNM_18802014'
#
# Current directories are:
# INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
# #INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# Set up your choices: 
#     PctDataPresent = 0.6	# can be any value between 0 and 1
#     StYr = 1973 # Start Year of dataset (assumes January)
#     EdYr = 2014 # End Year of dataset (assumes December)
#     # Are you reading in actuals or anomalies? 'anoms' for anomalies, 'abs' for actuals
#     vartype = 'anoms' 
#     # Select your month of choice, or a range for an average
#     # 0...11 represent Jan..Dec, [2,4] for Mar-Apr-May average, [0,11] for annual average, [11,1] for Dec-Jan-Feb average
#     # For month ranges than span 11 to 0, December will be taken from the first year of ChooseYr - will NOT work for last year!
#     ChooseMon = [11] 
#     # Select your year of choice, or a range for an average
#     # 1973...2014 for individual years, [1973,1982] for decadal average etc
#     ChooseYr = [2014] 
#     # Choose your start year of climatology: 0 if not relevant, not important if DoReClim = False
#     ChooseClimSt = 1976 
#     # Choose your start year of climatology: 0 if relevant, not important if DoReClim = False
#     ChooseClimEd = 2005 
#     # Does the data need renormalising to the chosen climatology? True if need to renormalise, False if not
#     DoReClim = False 
#     # Are we plotting anomalies or absolute values? 'actual' for actual values, 'anomaly' for anomaly values
#     PlotType = 'actual' 
#
# Set up correct file paths:
#     # Input Directories:
#     INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
#     #INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
#     # Output Directories:
#     OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'
#
# Set up correct files and info bundles:
#     # Files and bundle example:
#     # Input file:
#     candidate = 'HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#     # Output file:
#     OUTPLOT = 'Map_HadISDH.landq.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonStr+YrStr
#     # Variable name used in netCDF file, based on vartype
#     NameOfVar = 'q_'+vartype
#     # Long name of variable for plot
#     VarName = 'Specific Humidity'
#     # Unit for variable - used for plot
#     Unit = 'g kg$^{-1}$'  
#     # Plot label
#     Namey = 'HadISDH.landq.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#     # Latitude info: variable name, number of latitudes, start latitude
#     LatInfo = list(['latitude',36,-87.5]) 
#     # Longitude info: variable name, number of longitudes, start longitude
#     LonInfo = list(['longitude',72,-177.5]) 
#     # Choice of colourmap for plot and whether it should be flipped: 'flip','noflip'
#     ColourMapChoice = ('BrBG','noflip')
#     # Optional hardwire info for plot - can be 0,0,0,' ' - use if wanting to make identical colour bars
#     PlotInfo = [0,28,14,' '] # vmin,vmax,nsteps,letter_for_plot_label 
#
# python2.7 PlotAnyMap_JUN2015.py
#
# -----------------------
# OUTPUT
# -----------------------
# An .eps and .png map
# OUTDIR = '/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'
# Example output filename built from run choices:
# OUTPLOT = 'Map_HadISDH.landq.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 (13th October 2015)
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
import pdb # pdb.set_trace() or c

# Set in stone things:
MonArr = (['January','February','March','April','May','June','July','August','September','October','November','December'])

# Missing data
mdi = -1e30 # may set up as masked arrays later

# Set up initial run choices
# What is the minimum percent of data present for each gridbox value? This can be any value between >0 and 1
PctDataPresent = 0.6	

# Start Year of dataset (assumes January)
StYr = 1973 		

# End Year of dataset (assumes December)
EdYr = 2014 		

# Are you reading in actuals or anomalies? 'anoms' for anomalies, 'abs' for actuals
vartype = 'abs' 

# Select your month of choice, or a range for an average
# 0...11 represent Jan..Dec, [2,4] for Mar-Apr-May average, [0,11] for annual average, [11,1] for Dec-Jan-Feb average
# For month ranges than span 11 to 0, December will be taken from the first year of ChooseYr - will NOT work for last year!
ChooseMon = [11,1] 

# Select your year of choice, or a range for an average
# 1973...2014 for individual years, [1973,1982] for decadal average etc
ChooseYr = [1976,2005] 

# Choose your start year of climatology: 0 if not relevant, not important if DoReClim = False
ChooseClimSt = 1976 

# Choose your start year of climatology: 0 if relevant, not important if DoReClim = False
ChooseClimEd = 2005 

# Does the data need renormalising to the chosen climatology? True if need to renormalise, False if not
DoReClim = False 

# Are we plotting anomalies or absolute values? 'actual' for actual values, 'anomaly' for anomaly values
PlotType = 'actual' 

# From these choices now work out output file string and plot title label
if (len(ChooseMon) > 1):
    MonStr = MonArr[ChooseMon[0]]+'-'+MonArr[ChooseMon[1]]
else:
    MonStr = MonArr[ChooseMon[0]]
if (len(ChooseYr) > 1):
    YrStr = str(ChooseYr[0])+'-'+str(ChooseYr[1])
else:
    YrStr = str(ChooseMon[0])

# Set up directories:
INDIRC = '/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
#INDIRC = '/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
OUTDIR = '/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'

# Set up files and variable bundles:
#candidate = 'HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landq.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 'q_'+vartype
#VarName = 'Specific Humidity'
#Unit = 'g kg$^{-1}$'  
#Namey = 'HadISDH.landq.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('BrBG','noflip')
#PlotInfo = [0,28,15,' '] # vmin,vmax,nsteps+1,letter_for_plot_label 

candidate = 'HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
OUTPLOT = 'Map_HadISDH.landRH.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
NameOfVar = 'rh_'+vartype
VarName = 'Relative Humidity'
Unit ='%rh'  
Namey = 'HadISDH.landRH.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
ColourMapChoice = ('BrBG','noflip')
PlotInfo = [1,99,15,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landT.2.0.1.2014p_RAW_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 't_'+vartype
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'HadISDH.landT.2.0.1.2014p RAW '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landT.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 't_'+vartype
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'HadISDH.landT.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'BERKELEY_T_5by519762005clim_anoms_19732014'
#OUTPLOT = 'Map_BERKELEY_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 't_'+vartype
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'BERKELEY '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'GISS_T_5by519762005clim_anoms_19732014'
#OUTPLOT = 'Map_GISS_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 't_'+vartype
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'GISS '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'CRUTEM.4.3.0.0.anomalies'
#OUTPLOT = 'Map_CRUTEM4.3.0.0_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 't_'+vartype
#VarName = 'Temperature'
#Unit = '$^{o}$C' 
#Namey = 'CRUTEM4.3.0.0 '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'GHCNM_18802014'
#OUTPLOT = 'Map_GHCNM3_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = 't_'+vartype
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'GHCNM3 '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = list(['latitude',36,-87.5]) # variable name, number of latitudes, start latitude
#LonInfo = list(['longitude',72,-177.5]) # variable name, number of longitudes, start longitude
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

# Set up other variables
NYrs = (EdYr+1)-StYr
NMons = NYrs*12
if (len(ChooseMon) > 1) | (len(ChooseYr) > 1): # Only need to point to relevant years, this works if BOTH are > 1
    StTimePointer = (ChooseYr[0]-StYr)
    if (len(ChooseYr) > 1):
        EdTimePointer = (ChooseYr[1]-StYr)
    else:
        EdTimePointer = StTimePointer
else:
    StTimePointer = (ChooseYr[0]-StYr)*12+ChooseMon[0]
    EdTimePointer = StTimePointer
TmpCandFields = np.reshape(np.repeat(mdi,NMons*LonInfo[1]*LatInfo[1]),(NMons,LatInfo[1],LonInfo[1])) # entire field of data read in
CandFields = np.reshape(np.repeat(mdi,NMons*LonInfo[1]*LatInfo[1]),(NMons,LatInfo[1],LonInfo[1])) # entire field of data read in
CandData = np.reshape(np.repeat(mdi,LonInfo[1]*LatInfo[1]),(LatInfo[1],LonInfo[1])) # chosen month/average of months to plot
LatList = [] # list of latitude gridbox centres
LonList = [] # list of longitude gridbox centres
    
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
	ReadInfo: a string variable names for the grid to read in
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
    if (LonInfo[2] < 10):
        TheLonList=np.arange(LonInfo[2], LonInfo[2]+360.,(360./LonInfo[1]))
    else:
        TheLonList=np.arange(LonInfo[2], LonInfo[2]-360.,-(360./LonInfo[1]))    

    var=ncf.variables[ReadInfo]
    TheData=np.array(var.data) # times, lats, lons

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # ReadNetCDFGrid

#************************************************************************
# Get Anomalies
def GetAnomalies(TheData,TheMDI,TheStYr,TheEdYr,TheClimSt,TheClimEd,TheMDITol,TheAnoms):
    ''' This code takes in a 3D gridded field of monthly mean/anomalies 
        and anomalises/renormalises (climate anomalies) to a new climatology period
	which is supplied along side start and end years of the data.
	It assumes all data go from Jan to Dec '''
    ''' TheData: A 3D numpy array (months, latitude, longitude) of monthly means or anomalies '''
    ''' TheMDI: The missing data indicator '''
    ''' TheStYr: The start year of the provided data '''
    ''' TheEdYr: The end year of the provided data '''
    ''' TheClimSt: The start year of the new climatology period '''
    ''' TheClimEd: The end year of the new climatology period '''
    ''' TheMDITol: The percentage of data required for a gridbox climatology to be calculated '''
    ''' TheAnoms: a 3D array identical in shape to TheData for the output, filled with missing data indicator ''' 	
    
    # Work out how many years of data in total and for the climatology period
    TheYCount = (TheEdYr - TheStYr) + 1
    TheCCount = (TheClimEd - TheClimSt) + 1
    
    # Loop through each latitude and longitude to work on a gridbox at a time
    for lnn in range(len(TheData[0,0,:])):
        for ltt in range(len(TheData[0,:,0])):
	    subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(TheYCount,12)))
	    newsubarr = np.empty_like(subarr)
	    newsubarr.fill(TheMDI)
            # loop through each month
	    for mm in range(12):
	        subclimarr = np.copy(subarr[(TheClimSt-TheStYr):(TheClimEd-TheStYr)+1,mm])	# always need +1 for working with subranges
		subsubarr = np.copy(subarr[:,mm])
		subclimarr[subclimarr == TheMDI] = np.nan
		subsubarr[subsubarr == TheMDI] = np.nan
		#print(np.float(len(subclimarr[np.isfinite(subclimarr)]))/np.float(TheCCount))
		if (np.float(len(subclimarr[np.isfinite(subclimarr)]))/np.float(TheCCount) >= TheMDITol):
		    subsubarr[np.isfinite(subsubarr)] = subsubarr[np.isfinite(subsubarr)] - np.nanmean(subclimarr)
		    subsubarr[np.isnan(subsubarr)] = TheMDI
		    newsubarr[:,mm] = np.copy(subsubarr)
            TheAnoms[:,ltt,lnn] = np.copy(np.reshape(newsubarr,TheYCount*12))
	    
    return TheAnoms # GetAnomalies
    
#************************************************************************
# PlotAnyMap
def PlotAnyMap(TheFile,TheLatList,TheLonList,TheCandData,TheUnitee,TheNamee,TheCmap,TheTypee,TheVar,TheMDI,ThePlotInfo):
    ''' Create a masked array of variable
        Plot on map
        Save as eps and png '''
      
    # Create the masked array of trends
    MSKTheCandData=ma.masked_where(TheCandData == mdi,TheCandData)
    
    # make 2d arrays of lats and lons
    # nudge -2.5 degrees to make them south/west gridbox corners, not centres
    # add extra row/column to bound the data
    ArrLons,ArrLats=np.meshgrid(TheLonList,TheLatList)
    HalfLon=(TheLonList[1]-TheLonList[0])/2.
    HalfLat=(TheLatList[1]-TheLatList[0])/2.
    LngArrLons,LngArrLats=np.meshgrid(np.append(TheLonList-HalfLon,180.),np.append(TheLatList-HalfLat,90.))
    
    # set up plot
    plt.clf()
    fig=plt.figure(figsize=(10,6))
    plt1=plt.axes([0.05,0.12,0.9,0.8]) # left, bottom, width, height
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))

    # make up a blue to red (reverse) colour map
    cmap=plt.get_cmap(TheCmap[0])     
    
    cmaplist=[cmap(i) for i in range(cmap.N)]
    for loo in range((cmap.N/2)-20,(cmap.N/2)+20):
        cmaplist.remove(cmaplist[(cmap.N/2)-20]) # remove the very pale colours in the middle
    if (TheCmap[1] == 'flip'):	# then reverse the colours
        cmaplist.reverse()
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)

    # Is there isn't a hardwired vmin and vmax then make one that is nice
    if (ThePlotInfo[2] == 0):
        # work out best max and min values for colourbar and try to make them 'nice
        # must be an odd number of steps
        # for less than 0.1 must be 14 or fewer steps
        if (TheTypee == 'anomaly'):
            vmax=np.int(np.ceil(np.max(abs(MSKTheCandData))*10))/10.    
            vmin=-vmax
            print(vmin,vmax)
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
        else:
            vmax=np.ceil(np.max(abs(MSKTheCandData))) 
            vmin=np.floor(np.min(abs(MSKTheCandData)))
	    vrange=vmax-vmin
            print(vmin,vmax,vrange)
            if (vmax-vmin <  14):
                nsteps=np.int((vmax-vmin))+1
            elif (vmax <= 21):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(1.5*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(1.5*7) 
                nsteps=np.int((vmax-vmin)/1.5)+1
            elif (vmax <= 30):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(2*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(2*7) 
                nsteps=np.int((vmax-vmin)/2)+1
            elif (vmax <= 35):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(2.5*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(2.5*7) 
                nsteps=np.int((vmax-vmin)/2.5)+1
            elif (vmax <= 42):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(3*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(3*7) 
                nsteps=np.int((vmax-vmin)/3)+1
            elif (vmax <= 56):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(4*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(4*7) 
                nsteps=np.int((vmax-vmin)/4)+1
            elif (vmax <= 70):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(5*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(5*7) 
                nsteps=np.int((vmax-vmin)/5)+1
            elif (vmax <= 84):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(6*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(6*7) 
                nsteps=np.int((vmax-vmin)/6)+1
            elif (vmax <= 98):
                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(7*7) 
                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(7*7) 
                nsteps=np.int((vmax-vmin)/7)+1
	    #    pdb.set_trace() # stop here and play
    else:
        vmin = ThePlotInfo[0]
	vmax = ThePlotInfo[1]
	nsteps = ThePlotInfo[2]
    print(vmin,vmax,nsteps)
	
    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCandData,cmap=cmap,norm=norm,latlon='TRUE')

    cbax=fig.add_axes([0.05,0.08,0.9,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=12) 
    plt.figtext(0.5,0.01,TheVar+' ('+TheUnitee+')',size=14,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.5,0.95,TheNamee,size=16,ha='center')
    # If a letter is supplied then print on plot
    if (ThePlotInfo[3] != ' '):
        plt.figtext(0.1,0.92,ThePlotInfo[3],size=16,ha='center')            
    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotAnyMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in trend maps
MyFile=INDIRC+candidate+'.nc'
#print(NameOfVar,LatInfo,LonInfo)
TmpCandFields,LatList,LonList = ReadNetCDFGrid(MyFile,NameOfVar,LatInfo,LonInfo)

##########################################################################
## TESTING CODE ##########################################################
##########################################################################
## create a data array with an identical field for each month within year but increments annually
#TmpCandFields[:,:,:] =  np.reshape(np.array(np.repeat(range(NYrs),12*LatInfo[1]*LonInfo[1]),dtype=float),(NMons,LatInfo[1],LonInfo[1]))
#
#if (LatInfo[2] < 0):
#    LatList=np.arange(LatInfo[2], LatInfo[2]+180.,(180./LatInfo[1]))
#else:
#    LatList=np.arange(LatInfo[2], LatInfo[2]-180.,-(180./LatInfo[1]))    
#if (LonInfo[2] < 10):
#    LonList=np.arange(LonInfo[2], LonInfo[2]+360.,(360./LonInfo[1]))
#else:
#    LonList=np.arange(LonInfo[2], LonInfo[2]-360.,-(360./LonInfo[1]))    
#
## Testing GetAnomalies by choosing to anomalise relative to 1981-2010
#ChooseMon = [0,11] 
## Select your year of choice, or a range for an average
## 1973...2014 for individual years, [1973,1982] for decadal average etc
#ChooseYr = [1981,2010] 
## Choose your start year of climatology: 0 if not relevant, not important if DoReClim = False
#ChooseClimSt = 1981 
## Choose your start year of climatology: 0 if relevant, not important if DoReClim = False
#ChooseClimEd = 2010 
## Does the data need renormalising to the chosen climatology? True if need to renormalise, False if not
#DoReClim = True 
## Are we plotting anomalies or absolute values? 'actual' for actual values, 'anomaly' for anomaly values
#PlotType = 'anomaly' 
## From these choices now work out output file string and plot title label
#if (len(ChooseMon) > 1):
#    MonStr = MonArr[ChooseMon[0]]+'-'+MonArr[ChooseMon[1]]
#else:
#    MonStr = MonArr[ChooseMon[0]]
#if (len(ChooseYr) > 1):
#    YrStr = str(ChooseYr[0])+'-'+str(ChooseYr[1])
#else:
#    YrStr = str(ChooseMon[0])
#if (len(ChooseMon) > 1) | (len(ChooseYr) > 1): # Only need to point to relevant years, this works if BOTH are > 1
#    StTimePointer = (ChooseYr[0]-StYr)
#    if (len(ChooseYr) > 1):
#        EdTimePointer = (ChooseYr[1]-StYr)
#    else:
#        EdTimePointer = StTimePointer
#else:
#    StTimePointer = (ChooseYr[0]-StYr)*12+ChooseMon[0]
#    EdTimePointer = StTimePointer
##########################################################################    

# Do you need to renormalise (or create anomalies from actuals?
if (DoReClim):
    CandFields = GetAnomalies(TmpCandFields,mdi,StYr,EdYr,ChooseClimSt,ChooseClimEd,PctDataPresent,CandFields)
else:
    CandFields = TmpCandFields    

##########################################################################
## TESTING CODE ##########################################################
##########################################################################
## Check if GetAnomalies works - tested for 1981-2010
## Make sure DoReClim = True 
#moo=np.array(range(42),dtype=float)
#print(np.mean(moo[(ChooseClimSt-StYr):(ChooseClimEd-StYr)+1])) # need +1 for ranges
## 22.5
#print(moo-np.mean(moo[(ChooseClimSt-StYr):(ChooseClimEd-StYr)+1])) # need +1 for ranges
##array([-22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5
## -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5
##   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5
##  13.5  14.5  15.5  16.5  17.5  18.5]
#
#foo=np.reshape(CandFields[:,0,0],(NYrs,12))
#print(foo[:,0])
##array([-22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5
##       -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5
##         1.5	2.5   3.5   4.5   5.5	6.5   7.5   8.5   9.5  10.5  11.5  12.5
##        13.5  14.5  15.5  16.5  17.5  18.5]
#
#pdb.set_trace()#
#
## GetAnomalies works!
##########################################################################

# Extract chosen month
# Easy case of single month from single year
if (len(ChooseMon) == 1) & (len(ChooseYr) == 1):
    print("One Month One Year")
    CandData = CandFields[StTimePointer,:,:]
# Easy case of single month from multiple years - requires ?% presence for each gridbox
elif (len(ChooseMon) == 1) & (len(ChooseYr) > 1):
    print("One Month X Year")
    for lnn in range(LonInfo[1]):
        for ltt in range(LatInfo[1]):
	    subarr = np.copy(np.reshape(CandFields[:,ltt,lnn],(NYrs,12))) # fills columns first
            subsubarr = np.copy(subarr[StTimePointer:EdTimePointer+1,ChooseMon[0]]) # need +1 for ranges
	    if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= PctDataPresent):
	        CandData[ltt,lnn] = np.nanmean(subsubarr)
elif (len(ChooseMon) > 1) & (len(ChooseYr) == 1):
    print("X Month One Year")
    if (ChooseMon[1] > ChooseMon[0]):	# simple run of months within a year
        for lnn in range(LonInfo[1]):
            for ltt in range(LatInfo[1]):
	        subarr = np.copy(np.reshape(CandFields[:,ltt,lnn],(NYrs,12))) 
	        subsubarr = np.copy(subarr[StTimePointer,ChooseMon[0]:ChooseMon[1]+1]) # need +1 for ranges
		subsubarr[subsubarr == mdi] = np.nan # set mdi to NaN
	        if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= PctDataPresent):
	            CandData[ltt,lnn] = np.nanmean(subsubarr)
    else: # more complex as need to pull out from two years
        for lnn in range(LonInfo[1]):
            for ltt in range(LatInfo[1]):
	        subarr = np.copy(np.reshape(CandFields[:,ltt,lnn],(NYrs,12)))
	        subsubarr = np.copy(subarr[StTimePointer,ChooseMon[0]:12]) # need +1 for ranges
		subsubarr = np.append(subsubarr,subarr[StTimePointer+1,0:ChooseMon[1]+1]) # need +1 for ranges
		subsubarr[subsubarr == mdi] = np.nan # set mdi to NaN
	        if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= PctDataPresent):
	            CandData[ltt,lnn] = np.nanmean(subsubarr)
else: # now we're dealing with seasonal/annual average climatology
    print("X Month X Year")
    if (ChooseMon[1] > ChooseMon[0]):	# simple run of months and run of years
        for lnn in range(LonInfo[1]):
            for ltt in range(LatInfo[1]):
	        subarr = np.copy(np.reshape(CandFields[:,ltt,lnn],(NYrs,12)))
	        subsubarr = np.copy(subarr[StTimePointer:EdTimePointer+1,ChooseMon[0]:ChooseMon[1]+1]) # need +1 for ranges
		subsubarr[subsubarr == mdi] = np.nan # set mdi to NaN
	        if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= PctDataPresent):
	            CandData[ltt,lnn] = np.nanmean(subsubarr)
    else: # more complex as need to pull out month runs across years 
        if (EdTimePointer < EdYr): # then we can go to the next year to get the extra months
            ExtraPointer=EdTimePointer+1
        else:
	    ExtraPointer=EdTimePointer
	for lnn in range(LonInfo[1]):
            for ltt in range(LatInfo[1]):
	        subarr = np.copy(np.reshape(CandFields[:,ltt,lnn],(NYrs,12)))
	        subsubarr = np.copy(subarr[StTimePointer:EdTimePointer+1,ChooseMon[0]:12,]) # need +1 for ranges
		subsubarr = np.append(subsubarr,subarr[StTimePointer+1:ExtraPointer,0:ChooseMon[1]+1])
		subsubarr[subsubarr == mdi] = np.nan # set mdi to NaN
	        if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= PctDataPresent):
	            CandData[ltt,lnn] = np.nanmean(subsubarr)

##########################################################################
## TESTING CODE ##########################################################
##########################################################################
## Check the selection output works on actual values - all should be ltt,lnn arrays of identical numbers
## Make sure DoReClim = False, PlotType = 'actual' 
#print(CandData[0:6,0:6])
#pdb.set_trace()
## One month, one year: tested for June, 1980 = 7
## This works!
## One month, multiple years: tested October, 2000-2010 = mean of 27:37 = 32
## This works!
## Multiple months, one year: tested MAM, 1991 = mean of [18,18,18] = 18, tested DJF, 1992 = mean of [19,20,20] = 19.66666.
## This works for both!
## Multiple months, multiple years: tested SON, 1973-1982 = mean of 0:9,0:9,0:9 = 4.5, tested JAN-DEC, 1981-2010 mean of 8:37, 30 times = 22.5
## This works for both!
##########################################################################

# pass to plotter
MyFile=OUTDIR+OUTPLOT
PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,ColourMapChoice,PlotType,VarName,mdi,PlotInfo)
		
#    stop()

print("And, we are done!")

