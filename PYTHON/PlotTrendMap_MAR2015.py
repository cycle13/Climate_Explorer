#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 22 April 2015
# Last update: 21 Jul 2020
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	# this will probably change
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Reads in decadal trend map netCDF file (any 5 by 5 degree grid)
# Works out the mean trend for each latitude (simple average)  
# Plots the trend for each gridbox, grixbox boundary OR + IF trend is sig (5-95% Conf Int of same sign as trend)
# Also adds a vertical latitude by trend figure - scatter of gridbox trends for each latitude and latitude average
# Fills 0-1 with dark grey to represent actual land or ocean fraction in terms of gridboxes at that latitude band
# Overlays light grey to represent fraction of land or ocean gridboxes actually observed as % of dark grey bar
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
# import cartopy.crs as ccrs
# import cartopy.feature as cpf
## from mpl_toolkits.basemap import Basemap
# import datetime as dt
# from matplotlib.dates import date2num,num2date
# from scipy.io import netcdf
# import matplotlib.colors as mc
# import matplotlib.cm as mpl_cm
# import pdb
#
# Other:
# ReadNetCDFGrid - infile function to read in netCDF grid, written by Kate Willett
# ReadNetCDFGrid4 - infile function to read in netCDF4 grid, written by Kate Willett
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
## python2.7 PlotTrendMap_MAR2015.py
# This is to use python 3
# >module load scitools/default-current
# >python PlotTrendMap_MAR2014.py
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
# Version 6 (21 Jul 2020)
# ---------
#  
# Enhancements
# Now doesn't use Bundles so should be easier to use - need to set up typee to determine land, marine, blend
# It means it doesn't work with other datasets so will need tweaking a little for non-HadISDH
#  
# Changes
#  
# Bug fixes
#
#
# Version 5 (25 Jan 2020)
# ---------
#  
# Enhancements
# Now works with NetCDF4 and OLS trends
# HOwever, this really needs reworking as its very tedious
#  
# Changes
#  
# Bug fixes
#
#
#
# Version 4 (2 April 2019)
# ---------
#  
# Enhancements
# Now can plot significance as a + (for IPCC) if desired
# Also the latitude distribution is optional
#  
# Changes
# This is now python 3 - so uses cartopy rather than basemap
#  
# Bug fixes
#
# Version 3 (1 February 2016)
# ---------
#  
# Enhancements
# Now reads in data in bundles with only one line to comment in/out.
# Now has an extra level for sorting out plotting intervals for the map - prettier
#  
# Changes
#  
# Bug fixes
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
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cpf
import datetime as dt
from matplotlib.dates import date2num,num2date
from scipy.io import netcdf
import netCDF4 as nc4
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

# Generic things:******************************************

# Do you want boundary boxes or pluses?
TrendSig = 'BB' # 'BB' for boundary boxes, '+' for plus sign

# Do you want the latitudinal distribution plotting?
# NOTE: This will affect the filename
LatDist = True # True for having the latitudinal distribution, False for not

# Missing data
mdi=-1e30 # may set up as masked arrays later

# Input date stamp
thenmon='JAN'
thenyear='2020'

# Working version
lversion = '4.2.0.2019f' # land
mversion = '1.0.0.2019f' # marine
bversion = '1.0.0.2019f' # blend

# Working Directory
WorkingDir = 'UPDATE2019'

# Output date stamp
nowmon='JAN'
nowyear='2020'

# Trend selection
sttrend = '1973'
edtrend = '1999'
trendchoice = sttrend+edtrend

# Variable
#Var = 'q'
#Var = 'rh'
#Var = 'e'
#Var = 't'
#Var = 'td'
#Var = 'tw'
Var = 'dpd'

VarDict = dict([('q',['q','g kg$^{-1}$','Specific Humidity',('BrBG','noflip'),dict([('MinVal',-0.3),('MaxVal',0.3),('StepVal',9.),('LetterVal',['a)','b)'])])]),
                ('rh',['RH','%rh','Relative Humidity',('BrBG','noflip'),dict([('MinVal',-3.),('MaxVal',3.),('StepVal',9.),('LetterVal',['a)','b)'])])]),
		('e',['e','hPa','Vapour Pressure',('BrBG','noflip'),dict([('MinVal',-0.5),('MaxVal',0.5),('StepVal',9.),('LetterVal',['a)','b)'])])]),
		('t',['T','$^{o}$C','Air Temperature',('coolwarm','noflip'),dict([('MinVal',-0.8),('MaxVal',0.8),('StepVal',9.),('LetterVal',['a)','b)'])])]),
		('tw',['Tw','$^{o}$C','Wetbulb Temperature',('BrBG','noflip'),dict([('MinVal',-0.6),('MaxVal',0.6),('StepVal',9.),('LetterVal',['a)','b)'])])]),
		('td',['Td','$^{o}$C','Dew Point Temperature',('BrBG','noflip'),dict([('MinVal',-0.8),('MaxVal',0.8),('StepVal',9.),('LetterVal',['a)','b)'])])]),
		('dpd',['DPD','$^{o}$C','Dew Point Depression',('BrBG','flip'),dict([('MinVal',-0.5),('MaxVal',0.5),('StepVal',9.),('LetterVal',['a)','b)'])])])])

if (Var == 'dpd'):
    DatTyp = 'PHA'
elif (Var == 'td'):
    DatTyp = 'PHADPD'
else:
    DatTyp = 'IDPHA'

# Do you want to make up the best colour ranges or use set ones?
#ColourRange = dict([('MinVal',0.),('MaxVal',0.),('StepVal',0.),('LetterVal',['a)','b)'])]) # the blank version so that the code finds its own
ColourRange = VarDict[Var][4]

# Type
#typee = 'LAND' # 'LAND','RAW','OTHER'
typee = 'MARINESHIP' # 'MARINE','MARINESHIP'
#typee = 'BLENDSHIP' # 'BLEND','BLENDSHIP'

# Domain
if (typee == 'MARINE') | (typee == 'MARINESHIP'):
    domain = 'marine'
    IsLand = False # True for land, False for marine, None for blend
    version = mversion
    if (typee == 'MARINESHIP'):
        ExpBit = 'BClocalSHIP5by5both'
    else:
        ExpBit = 'BClocal5by5both'
elif (typee == 'BLEND') | (typee == 'BLENDSHIP'):
    domain = 'blend'
    IsLand = None # True for land, False for marine, None for blend
    version = bversion
    if (typee == 'BLENDSHIP'):
        ExpBit = 'FLATgrid'+DatTyp+'BClocalSHIPboth5by5'
    else:
        ExpBit = 'FLATgrid'+DatTyp+'BClocalSHIPboth5by5'
else:
    domain = 'land'
    IsLand = True # True for land, False for marine, None for blend
    version = lversion
    ExpBit = 'FLATgrid'+DatTyp+'5by5'
    
# TrendType
TrendType = 'OLS' # 'OLS' or 'MP'

# Set up directories and files - now for /data/users/hadkw/WORKING_HADISDH/
INDIRC='/data/users/hadkw/WORKING_HADISDH/'+WorkingDir+'/STATISTICS/TRENDS/'
INDIRO='/data/users/hadkw/WORKING_HADISDH/'+WorkingDir+'/OTHERDATA/'
OUTDIR='/data/users/hadkw/WORKING_HADISDH/'+WorkingDir+'/IMAGES/OTHER/'

# Land cover file:
#incover='new_coverpercentjul08'
incover='HadCRUT.4.3.0.0.land_fraction'

# Variables
CandTrends=[]
CompTrends=[]
RatTrends=[]
LatList=[]
LonList=[]

candidate = 'HadISDH.'+domain+VarDict[Var][0]+'.'+version+'_'+ExpBit+'_anoms8110_'+thenmon+thenyear+'_'+TrendType+'trends_'+trendchoice
OUTPLOT = 'TrendMap'+TrendType+'_HadISDH.'+domain+VarDict[Var][0]+'.'+version+'_'+ExpBit+'_'+trendchoice
Namey = 'HadISDH.'+domain+VarDict[Var][0]+'.'+version+' '+TrendType+' decadal trends '+trendchoice

Unit = VarDict[Var][1]
nlats = 36	       #set once file read in
nlons = 72	       #set once file read in
LatInfo = list(['latitude',nlats,-87.5])
LonInfo = list(['longitude',nlons,-177.5])
ReadInfo = list([VarDict[Var][0]+'_trend',VarDict[Var][0]+'_lowCI',VarDict[Var][0]+'_upperCI'])
ColourMapChoice = VarDict[Var][3]


#**********************************************************
# Run choice bundle for input/output files, units, names, letters, read in varnames, colourmap

# CHOOSE/ADD A DICTIONARY BUNDLE!!!
MyBundle = ' '

#MyBundle = 'Berkeley'
#MyBundle = 'GISS'
#MyBundle = 'CRUTEM'
#MyBundle = 'GHCNM'

if (MyBundle == 'Berkeley'):
    candidate='BERKELEY_T_5by519762005clim_anoms_19732015_MPtrends_19732016'
    OUTPLOT='TrendMap_BERKELEY_'+nowmon+nowyear+'_'+trendchoice
    Unit='$^{o}$ C'  #'degrees C'
    Namey='BERKELEY decadal trends'
    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
    LatInfo=list(['latitude',nlats,-87.5])
    LonInfo=list(['longitude',nlons,-177.5])
    ReadInfo=list(['T_trend','T_lowCI','T_upperCI'])
    ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
    IsLand = True # True for land, False for marine, None for blend

if (MyBundle == 'GISS'):
    candidate='GISS_T_5by519762005clim_anoms_19732015_MPtrends_19732016'
    OUTPLOT='TrendMap_GISS_'+nowmon+nowyear+'_'+trendchoice
    Unit='$^{o}$ C'  #'degrees C'
    Namey='GISS decadal trends'
    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
    LatInfo=list(['latitude',nlats,-87.5])
    LonInfo=list(['longitude',nlons,-177.5])
    ReadInfo=list(['T_trend','T_lowCI','T_upperCI'])
    ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
    IsLand = True # True for land, False for marine, None for blend

if (MyBundle == 'CRUTEM'):
    candidate='CRUTEM.4.4.0.0.anomalies_MPtrends_19982014'
    OUTPLOT='TrendMap_CRUTEM4.4.0.0_'+nowmon+nowyear+'_'+trendchoice
    Unit='$^{o}$ C'  #'degrees C'
    Namey='CRUTEM4.4.0.0 decadal trends'
    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
    LatInfo=list(['latitude',nlats,-87.5])
    LonInfo=list(['longitude',nlons,-177.5])
    ReadInfo=list(['T_trend','T_lowCI','T_upperCI'])
    ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
    IsLand = True # True for land, False for marine, None for blend

if (MyBundle == 'GHCNM'):
    candidate='GHCNM_18802014_MPtrends_19732016'
    OUTPLOT='TrendMap_GHCNM3_'+nowmon+nowyear+'_'+trendchoice
    Unit='$^{o}$ C'  #'degrees C'
    Namey='GHCNM3 decadal trends'
    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
    LatInfo=list(['latitude',nlats,-87.5])
    LonInfo=list(['longitude',nlons,-177.5])
    ReadInfo=list(['T_trend','T_lowCI','T_upperCI'])
    ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
    IsLand = True # True for land, False for marine, None for blend


# If LatDist is False then add 'nolat' to file name
if not(LatDist):
    OUTPLOT = OUTPLOT+'_nolat'

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
# READNETCDFGRID4
def ReadNetCDFGrid4(FileName,ReadInfo,LatInfo,LonInfo):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data 
	FileName: stroing containing filepath/name
	ReadInfo: list of strings of variable names for the trend, 5th percentile, 95th percentile
	LatInfo: list of a string for the latitude variable name, an integer for number of lats, a float for the start latitude
	LonInfo: list of a string for the longitude variable name, an integer for number of lons, a float for the start longitude '''

    ncf=nc4.Dataset(FileName,'r')

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
    TheData=np.copy(var[:])
    var=ncf.variables[ReadInfo[1]]
    TheLower=np.copy(var[:])
    var=ncf.variables[ReadInfo[2]]
    TheUpper=np.copy(var[:])

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheUpper,TheLower,TheLatList,TheLonList # ReadNetCDFGrid

#************************************************************************
# PlotTrendMap
def PlotTrendMap(TheFile,LandCover,TheLatList,TheLonList,TheCandData,TheCandUppers,TheCandLowers,
                    TheUnitee,TheNamee,ColourMapChoice,IsLand,TrendSig,LatDist,TheColourRange):
    ''' Create a masked array of trends
        Create a masked array of significant trends
        Plot trends on map, bounded/+ if significant
        Add vertical latitude/trend scatter with average trend overlaid if LatDist == True
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
    
    #pdb.set_trace()
    
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
    
    #pdb.set_trace()
    
    # set up plot
    plt.clf()
    # Size depends on whether we're plotting LatDist
    if (LatDist):
        fig=plt.figure(figsize=(10,5))
        #plt.figure(1,figsize=(10,3))
        plt1=plt.axes([0.01,0.05,0.64,0.9],projection=ccrs.Robinson()) # left, bottom, width, height
# For basemap
#        plt1=plt.axes([0.01,0.05,0.64,0.9]) # left, bottom, width, height
    else:
        fig=plt.figure(figsize=(7,5))
        plt1=plt.axes([0.01,0.05,0.98,0.9],projection=ccrs.Robinson()) # left, bottom, width, height
# For basemap
#        plt1=plt.axes([0.01,0.05,0.98,0.9]) # left, bottom, width, height
    
    
#    # plot map without continents and coastlines
#    m = Basemap(projection='kav7',lon_0=0)
#    # draw map boundary, transparent
#    m.drawmapboundary()
#    m.drawcoastlines()
#    # draw paralells and medians, no labels
#    m.drawparallels(np.arange(-90,90.,30.))
#    m.drawmeridians(np.arange(-180,180.,60.))

# NOW USING CARTOPY
    # plot map without continents and coastlines
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    #ax = plt.axes(projection=ccrs.Robinson())
    #if (GlobSwitch == 1): # this is a global projection
    #    ax.set_global()
    #    ax.set_extent([-180.,180.,-90.0,90.0],crs=ccrs.PlateCarree())
    #
    #else: # this is regional
    #    ax.set_extent(RegionExtents) # THIS WILL NEED TESTING AT SOME POINT
	
    plt1.coastlines()
    #plt1.set_boundary # not sure what this does? maybe useful for plotting regions?
    plt1.gridlines(draw_labels=False) # probably a way to specify these exactly if we want something different to default
    # This line background fills the land with light grey
    #plt1.add_feature(cpf.LAND, zorder = 0, facecolor = "0.9", edgecolor = "k") # may or may not need this
    
    ext = plt1.get_extent()
    # End of CARTOPY

    # make up a blue to red (reverse) colour map
    cmap=plt.get_cmap(ColourMapChoice[0])
    
#    fullcols=list()
#    [fullcols.append(bluelist[i]) for i in range(len(bluelist))]
#    [fullcols.append(redlist[i]) for i in range(len(redlist))]
    
    cmaplist=[cmap(i) for i in range(cmap.N)]
#    pdb.set_trace()
    for loo in range(int((cmap.N/2)-10),int((cmap.N/2)+10)):
        cmaplist.remove(cmaplist[int((cmap.N/2)-10)]) # remove the very pale colours in the middle
    if (ColourMapChoice[1] == 'flip'):	# then reverse the colours
        cmaplist.reverse()
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    if (TheColourRange['StepVal'] == 0):
        # work out best max and min values for colourbar and try to make them 'nice
        # must be an odd number of steps
        # for less than 0.1 must be 14 or fewer steps
        vmax=np.int(np.ceil(np.max(abs(MSKTheCandData))*10))/10.    
        vmin=-vmax
        nsteps = 9
        if (vmax <= 0.2):
            vmax = 0.2
            vmin = -0.2
        if (vmax <= 0.3):
#        nsteps=np.int((vmax-vmin)/0.05)+1
            vmax = 0.32
            vmin = -0.32
        elif (vmax <= 0.4):
#        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.06)*0.06  
#        vmin=-vmax
#        nsteps=np.int((vmax-vmin)/0.06)+1
            vmax = 0.4
            vmin = -0.4
        elif (vmax <= 0.6):
#        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.08)*0.08  
#        vmin=-vmax
#        nsteps=np.int((vmax-vmin)/0.08)+1
            vmax = 0.6
            vmin = -0.6
        elif (vmax <= 0.8):
            vmax = 0.8
            vmin = -0.8
        elif (vmax <= 1.0):
#        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.1)*0.1  
#        vmin=-vmax
#        nsteps=np.int((vmax-vmin)/0.1)+1
            vmax = 1.0
            vmin = -1.0
        elif (vmax <= 1.2):
#        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.2)*0.2  
#        vmin=-vmax
#        nsteps=np.int((vmax-vmin)/0.2)+1
            vmax = 1.2
            vmin = -1.2
        elif (vmax <= 1.6):
#        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.2)*0.2  
#        vmin=-vmax
#        nsteps=np.int((vmax-vmin)/0.2)+1
            vmax = 1.6
            vmin = -1.6
        elif (vmax <= 2.0):
#        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.3)*0.3  
#        vmin=-vmax
#        nsteps=np.int((vmax-vmin)/0.5)+1
            vmax = 2.0
            vmin = -2.0
        elif (vmax <= 3.0):
            vmax = 3.0
            vmin = -3.0
	#    pdb.set_trace() # stop here and play

    else:
        nsteps = TheColourRange['StepVal']
        vmin = TheColourRange['MinVal']
        vmax = TheColourRange['MaxVal']

    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=plt1.pcolor(LngArrLons,LngArrLats,MSKTheCandData,transform = ccrs.PlateCarree(),cmap=cmap,norm=norm)
# Was basemap
#    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCandData,cmap=cmap,norm=norm,latlon='TRUE')
    # If We want boundary boxes then use pcolor, else, overplot +
    if (TrendSig == 'BB'):
        grids=plt1.pcolor(LngArrLons,LngArrLats,MSKSigTrends,transform = ccrs.PlateCarree(),cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5)
# Was basemap
#        grids=m.pcolor(LngArrLons,LngArrLats,MSKSigTrends,cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5,latlon='TRUE')

    elif (TrendSig == '+'):
        GotSigs = np.where(MSKSigTrends > mdi)
#        pdb.set_trace()
        for i in range(len(GotSigs[0])):
            plt.text(ArrLons[GotSigs[0][i],GotSigs[1][i]],ArrLats[GotSigs[0][i],GotSigs[1][i]],'+',va='center', ha='center',transform = ccrs.PlateCarree(),)
	
#	plt.scatter(ArrLons[np.where(MSKSigTrends > mdi)],ArrLats[np.where(MSKSigTrends > mdi)],c='black',marker='+',norm=norm,edgecolor='0.0',linewidth=0.0)
    
    
    # If LatDist then positions will differ
    if (LatDist):
        cbax=fig.add_axes([0.03,0.10,0.59,0.03])
        cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=10) 
        plt.figtext(0.31,0.01,'Trend ('+TheUnitee+' decade$^{-1}$ )',size=12,ha='center')
    else:
        cbax=fig.add_axes([0.03,0.09,0.94,0.03])
        cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=10) 
        plt.figtext(0.5,0.01,'Trend ('+TheUnitee+' decade$^{-1}$ )',size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.9,TheColourRange['LetterVal'][0],size=16)
# Switch this off if you don't want the title on the plot!!!
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

    # If LatDist then plot the latitude scatter
    if (LatDist):
        # Attempt to plot the latitude/ratio scatter
        ax2=plt.axes([0.73,0.12,0.24,0.76]) # map only
        ax2.set_ylim(-90,90)
        #ax2.set_ylabel('Latitude')
        ax2.set_yticks(list([-60,-30,0,30,60]))
        ax2.set_yticklabels(list(['-60$^{o}$N','-30$^{o}$S','Equator','30$^{o}$N','60$^{o}$N'])) #,rotation=90)
        ax2.set_xlabel('Trend ('+TheUnitee+' decade$^{-1}$ )')
        minx=np.min(TheCandData[np.where(TheCandData != mdi)])-0.1
        maxx=np.max(TheCandData[np.where(TheCandData != mdi)])+0.1
        ax2.set_xlim(minx,maxx)
        ax2.tick_params(axis='both',labelsize=10)
        # background fill the scatter plot with % land /ocean present (light grey) and % land /ocean with data (dark gre        
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

        plt.figtext(0.7,0.9,TheColourRange['LetterVal'][1],size=16)
    
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
CandData,CandUpper,CandLower,LatList,LonList=ReadNetCDFGrid4(MyFile,ReadInfo,LatInfo,LonInfo)
#CandData,CandUpper,CandLower,LatList,LonList=ReadNetCDFGrid(MyFile,ReadInfo,LatInfo,LonInfo)

ncf=netcdf.netcdf_file(INDIRO+incover+'.nc','r')
#var=ncf.variables['pct_land']
var=ncf.variables['land_area_fraction']
PctLand=np.array(var.data)
#PctLand=np.transpose(PctLand)
PctLand=np.flipud(PctLand)
PctLand=PctLand[0,:,:]
# If its marine data then swap to % ocean
if (IsLand == False):
    PctLand = 1. - PctLand
if (IsLand == None):
    PctLand[:,:] = 1.
ncf.close()

# pass to plotter
MyFile=OUTDIR+OUTPLOT
PlotTrendMap(MyFile,PctLand,LatList,LonList,CandData,CandUpper, CandLower,
                        Unit,Namey,ColourMapChoice,IsLand,TrendSig,LatDist,ColourRange)

# output pct of land /ocean boxes represented per region
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

