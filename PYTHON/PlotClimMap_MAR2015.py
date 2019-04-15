#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 22 April 2015
# Last update: 2 Apr 2019
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	# this will probably change
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Reads in monthyl netCDF file (any 5 by 5 degree grid)
# Works out the climatology (1981-2010 default) for each gridbox
# Can do monthly, seasonal or annual  
# Plots the climatology for each gridbox# 
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
# import datetime as dt
# from matplotlib.dates import date2num,num2date
## from scipy.io import netcdf
# import matplotlib.colors as mc
# import matplotlib.cm as mpl_cm
# import pdb
#
# Other:
# from ReadNetCDF import GetGrid4 - function to read in netCDF grid, written by Kate Willett
# from ReadNetCDF import GetGrid - function to read in netCDF grid, written by Kate Willett
# PlotClimMap - infile function to plot the map, written by Kate Willett
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
# HadISDH.landq.4.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014.nc
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
## python2.7 PlotClimMap_MAR2015.py
# This is to use python 3
# >module load scitools/default-current
# >python PlotClimMap_MAR2015.py
# 
# -----------------------
# OUTPUT
# -----------------------
# directory for output images:
# /data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/
# Output image files (nowmon+nowyear= e.g., OCT2015):
# ClimMap_<period>_19812010_HadISDH.landq.2.0.1.2014p_'+nowmon+nowyear+
# ClimMap_HadISDH.landRH.2.0.1.2014p_'+nowmon+nowyear
# ClimMap_HadISDH.landT.2.0.1.2014p_RAW_'+nowmon+nowyear
# ClimMap_HadISDH.landT.2.0.1.2014p_'+nowmon+nowyear
# ClimMap_BERKELEY_'+nowmon+nowyear
# ClimMap_GISS_'+nowmon+nowyear
# ClimMap_CRUTEM4.3.0.0_'+nowmon+nowyear
# ClimMap_GHCNM3_'+nowmon+nowyear
# ClimMap_HadISDH.lande.2.0.1.2014p_'+nowmon+nowyear+
# ClimMap_HadISDH.landTw.2.0.1.2014p_'+nowmon+nowyear+
# ClimMap_HadISDH.landTd.2.0.1.2014p_'+nowmon+nowyear+
# ClimMap_HadISDH.landDPD.2.0.1.2014p_'+nowmon+nowyear+
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 1 (5 April 2019)
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
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys, os
import scipy.stats
import struct
import cartopy.crs as ccrs
import cartopy.feature as cpf
import datetime as dt
from matplotlib.dates import date2num,num2date
#from scipy.io import netcdf
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

from ReadNetCDF import GetGrid4
from ReadNetCDF import GetGrid

# Generic things:******************************************

# What variable?
MyVar = 'td' # q, rh, e, t, td, tw, dpd
MyVarBig = 'Td' # q, RH, e, T, Td, Tw, DPD
# letter to annotate figure with
Letty=['d)']

# What Dataset
MyBundle = 'marine' # 'marine','land','blend','ERA-Interim' 

# What period do you want?
MyPeriod = 'ANN' # 'ANN','JAN','FEB','MAR','APR etc or 'MAM','JJA','SON','DJF'
PeriodDict = dict([('ANN',[0,1,2,3,4,5,6,7,8,9,10,11]),
                   ('JAN',[0]),
                   ('FEB',[1]),
                   ('MAR',[2]),
                   ('APR',[3]),
                   ('MAY',[4]),
                   ('JUN',[5]),
                   ('JUL',[6]),
                   ('AUG',[7]),
                   ('SEP',[8]),
                   ('OCT',[9]),
                   ('NOV',[10]),
                   ('DEC',[11]),
                   ('MAM',[2,3,4]),
                   ('JJA',[5,6,7]),
                   ('SON',[8,9,10]),
                   ('DJF',[11,0,1]),])

# Missing data
mdi=-1e30 # may set up as masked arrays later

# Input date stamp
thenmon='FEB'
thenyear='2019'

# Working version
version = '1.0.0.2018f'   # marine = '1.0.0.2018f' land = '4.1.0.2018f' ERA-I = ''

# Working Directory
WorkingDir = 'UPDATE2018'

# Output date stamp
nowmon='APR'
nowyear='2019'

# Climatology period
stcl = 1981
edcl = 2010
climchoice = str(stcl) + str(edcl)
styr = 1973
edyr = 2018

# Set up directories and files - now for /data/users/hadkw/WORKING_HADISDH/
INDIR = '/data/users/hadkw/WORKING_HADISDH/'+WorkingDir+'/STATISTICS/GRIDS/'
# INDIR FOR OTHER DATA
INDIRO = '/data/users/hadkw/WORKING_HADISDH/'+WorkingDir+'/OTHERDATA/'	# may be needed for non-HadISDH variables
# INDIR FOR MARINE DATA IN OLDER PROCESSING STEPS
#INDIR = '/project/hadobs2/hadisdh/marine/' # this could be otherdata, ICOADS.3.0.0/GRIDSOBSclim2NBC/ etc
OUTDIR = '/data/users/hadkw/WORKING_HADISDH/'+WorkingDir+'/IMAGES/OTHER/'
#OUTDIR = '/data/users/hadkw/WORKING_HADISDH/MARINE/IMAGES/CLIMMAPS/'

# Land cover file:
#incover='new_coverpercentjul08'
incover='HadCRUT.4.3.0.0.land_fraction'

# Variables
MonthlyData = []
ClimData = []
LatList = []
LonList = []


#**********************************************************
# Run choice bundle for input/output files, units, names, letters, read in varnames, colourmap
if (MyBundle == 'land'):
    if (MyVar == 'dpd'):
        candidate='HadISDH.land'+MyVarBig+'.'+version+'_FLATgridPHA5by5_anoms8110_'+thenmon+thenyear+'_cf'
    elif (MyVar == 'td'):
        candidate='HadISDH.land'+MyVarBig+'.'+version+'_FLATgridPHADPD5by5_anoms8110_'+thenmon+thenyear+'_cf'
    else:
        candidate='HadISDH.land'+MyVarBig+'.'+version+'_FLATgridIDPHA5by5_anoms8110_'+thenmon+thenyear+'_cf'

    OUTPLOT='ClimMap_'+MyPeriod+'_'+climchoice+'_HadISDH.land'+MyVarBig+'.'+version+'_'+nowmon+nowyear

    if (MyVar == 'q'):
        Unit='g kg$^{-1}$'  #'degrees C'
    elif (MyVar == 'e'):
        Unit='hPa'  #'degrees C'
    elif (MyVar == 'rh'):
        Unit='%rh'  #'degrees C'
    else:
        Unit='$^{o}$ C'  #'degrees C'

    #Namey='HadISDH.landq.'+version+' decadal trends '+trendchoice
    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
#    LatInfo=list(['latitude',nlats,-87.5])
#    LonInfo=list(['longitude',nlons,-177.5])
    LatInfo=list(['latitude'])
    LonInfo=list(['longitude'])
    ReadInfo=list([MyVar+'_clims'])
    if (MyVar == 't'):
        ColourMapChoice=('coolwarm','noflip')
    elif (MyVar == 'dpd'):
        ColourMapChoice=('BrBG','flip')
    else:    
        ColourMapChoice=('BrBG','noflip')   
    IsLand = True # True for land, False for marine, None for blend

if (MyBundle == 'marine'):

# Bias corrected version SHIP	
    if (MyVar == 'dpd'):
        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_BClocalSHIP5by5both_anoms8110_'+thenmon+thenyear+'_cf'
    elif (MyVar == 'td'):
        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_BClocalSHIP5by5both_anoms8110_'+thenmon+thenyear+'_cf'
    else:
        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_BClocalSHIP5by5both_anoms8110_'+thenmon+thenyear+'_cf'
    OUTPLOT='ClimMap_'+MyPeriod+'_'+climchoice+'_HadISDH.marine'+MyVarBig+'.'+version+'_BClocalSHIP5by5both_'+nowmon+nowyear

## bias corrected version ALL
#    if (MyVar == 'dpd'):
#        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_BClocal5by5both_anoms8110_'+thenmon+thenyear+'_cf'
#    elif (MyVar == 'td'):
#        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_BClocal5by5both_anoms8110_'+thenmon+thenyear+'_cf'
#    else:
#        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_BClocal5by5both_anoms8110_'+thenmon+thenyear+'_cf'
#    OUTPLOT='ClimMap_'+MyPeriod+'_'+climchoice+'_HadISDH.marine'+MyVarBig+'.'+version+'_BClocal5by5both_'+nowmon+nowyear


## QC'd, no bias corrected version SHIP
#    if (MyVar == 'dpd'):
#        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_NBCSHIP5by5both_anoms8110_'+thenmon+thenyear+'_cf'
#    elif (MyVar == 'td'):
#        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_NBCSHIP5by5both_anoms8110_'+thenmon+thenyear+'_cf'
#    else:
#        candidate='HadISDH.marine'+MyVarBig+'.'+version+'_NBCSHIP5by5both_anoms8110_'+thenmon+thenyear+'_cf'
#    OUTPLOT='ClimMap_'+MyPeriod+'_'+climchoice+'_HadISDH.marine'+MyVarBig+'.'+version+'_NBCSHIP5by5both_'+nowmon+nowyear

    if (MyVar == 'q'):
        Unit='g kg$^{-1}$'  #'degrees C'
    elif (MyVar == 'e'):
        Unit='hPa'  #'degrees C'
    elif (MyVar == 'rh'):
        Unit='%rh'  #'degrees C'
    else:
        Unit='$^{o}$ C'  #'degrees C'

    #Namey='HadISDH.landq.'+version+' decadal trends '+trendchoice

    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
#    LatInfo=list(['latitude',nlats,-87.5])
#    LonInfo=list(['longitude',nlons,-177.5])
    LatInfo=list(['latitude'])
    LonInfo=list(['longitude'])
    ReadInfo=list([MyVar+'_clims'])
    if (MyVar == 't'):
        ColourMapChoice=('coolwarm','noflip')
    elif (MyVar == 'dpd'):
        ColourMapChoice=('BrBG','flip')
    else:    
        ColourMapChoice=('BrBG','noflip')   
    IsLand = False # True for land, False for marine, None for blend


if (MyBundle == 'BLEND_HadISDH.landq'):
    candidate='BLEND_HadISDH.landq.'+version+'.marineq.QC0.0.0_'+thenmon+thenyear+'_MPtrends_'+trendchoice
    OUTPLOT='TrendMap_BLEND_HadISDH.landq.'+version+'_'+nowmon+nowyear+'_'+trendchoice
    Unit='g kg$^{-1}$'  #'degrees C'
    Namey='BLEND_HadISDH.landq.'+version+' decadal trends'
    nlats=36	       #set once file read in
    nlons=72	       #set once file read in
#    LatInfo=list(['latitude',nlats,-87.5])
#    LonInfo=list(['longitude',nlons,-177.5])
    LatInfo=list(['latitude'])
    LonInfo=list(['longitude'])
    ReadInfo=list(['q_MPtrend','q_MP5th','q_MP95th'])
    ColourMapChoice=('BrBG','noflip')
    IsLand = None # True for land, False for marine, None for blend


#************************************************************************
# Subroutines
#************************************************************************
# PlotClimMap
def PlotClimMap(TheFile,LandCover,TheLatList,TheLonList,TheCandData,ThePeriod,
                    TheUnitee,TheLetter,ColourMapChoice,IsLand):
    ''' Create a masked array of clims
        Plot clims on map
        Save as eps and png '''

    # Missing data
    mdi=-1e30 # may set up as masked arrays later
      
    # Create the masked array of trends
    MSKTheCandData=ma.masked_where(TheCandData == mdi,TheCandData)    
    #pdb.set_trace()
    
    # make 2d arrays of lats and lons
    # nudge -2.5 degrees to make them south/west gridbox corners, not centres
    # add extra row/column to bound the data
    ArrLons,ArrLats=np.meshgrid(TheLonList,TheLatList)
    if (TheLatList[0] > TheLatList[1]):
        LngArrLons,LngArrLats=np.meshgrid(np.append(TheLonList-2.5,180.),np.append(TheLatList+2.5,-90.))    
    else:
        LngArrLons,LngArrLats=np.meshgrid(np.append(TheLonList-2.5,180.),np.append(TheLatList-2.5,90.))    
    
    #pdb.set_trace()
    
    # set up plot
    plt.clf()
    fig=plt.figure(figsize=(7,5))
    plt1=plt.axes([0.01,0.05,0.98,0.9],projection=ccrs.Robinson()) # left, bottom, width, height
# For basemap
#    plt1=plt.axes([0.01,0.05,0.98,0.9]) # left, bottom, width, height
    
    
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
    
    
    # work out best max and min values for colourbar and try to make them 'nice
    # must be an odd number of steps
    # for less than 0.1 must be 14 or fewer steps
    vmax=np.int(np.ceil(np.max(MSKTheCandData)))    
    vmin=np.int(np.floor(np.min(MSKTheCandData))) 
    nsteps = 9
    
    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=plt1.pcolor(LngArrLons,LngArrLats,MSKTheCandData,transform = ccrs.PlateCarree(),cmap=cmap,norm=norm)
# Was basemap
#    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCandData,cmap=cmap,norm=norm,latlon='TRUE')
       
    cbax=fig.add_axes([0.03,0.09,0.94,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=10) 
    plt.figtext(0.5,0.01,ThePeriod+' Climatology ('+TheUnitee+')',size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.9,TheLetter[0],size=16)
# Switch this off if you don't want the title on the plot!!!
#    plt.figtext(0.5,0.95,TheNamee,size=16,ha='center')
    
#    print(CountGoods,CountLargeNegs,CountSmallNegs,CountSmallPos,CountLargePos)
#    pctLNs=round((float(CountLargeNegs)/CountGoods)*100.,1)
#    pctSNs=round((float(CountSmallNegs)/CountGoods)*100.,1)
#    pctVSs=round((float(CountVSmalls)/CountGoods)*100.,1)
#    pctSPs=round((float(CountSmallPos)/CountGoods)*100.,1)
#    pctLPs=round((float(CountLargePos)/CountGoods)*100.,1)
    
    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotClimMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in climatology fields
MyFile=INDIR+candidate+'.nc'
AllCandData,LatList,LonList=GetGrid4(MyFile,ReadInfo,LatInfo,LonInfo)
PctLand,LLatList,LLonList=GetGrid(INDIRO+incover+'.nc',['land_area_fraction'],['latitude'],['longitude'])

#ncf=netcdf.netcdf_file(INDIRO+incover+'.nc','r')
##var=ncf.variables['pct_land']
#var=ncf.variables['land_area_fraction']
#PctLand=np.array(var.data)
##PctLand=np.transpose(PctLand)
PctLand=np.flipud(PctLand)
PctLand=PctLand[0,:,:]
# If its marine data then swap to % ocean
if (IsLand == False):
    PctLand = 1. - PctLand
if (IsLand == None):
    PctLand[:,:] = 1.
#ncf.close()

# Get the climatology period of choice - month, season, annual?
CandData = np.empty((nlats,nlons),dtype = float)
CandData.fill(mdi)

#pdb.set_trace()

# Calculate the mean of climatological months present
# MISSING DATA: 11 out of 12 for annual, 2 out of 3 for seasons, 1 out of 1 for individual months? Should be very cautious about missing a month?
if (len(PeriodDict[MyPeriod]) == 12):
    MissingThresh = 9
elif (len(PeriodDict[MyPeriod]) == 3):
    MissingThresh = 2
elif (len(PeriodDict[MyPeriod]) == 1):
    MissingThresh = 1
    
for ltt in range(nlats):
    for lnn in range(nlons):
        TmpData = AllCandData[PeriodDict[MyPeriod],ltt,lnn]
        if (len(np.where(TmpData > mdi)[0]) >= MissingThresh):	
            CandData[ltt,lnn] = np.mean(TmpData[np.where(TmpData > mdi)])
	    

# pass to plotter
MyFile=OUTDIR+OUTPLOT
PlotClimMap(MyFile,PctLand,LatList,LonList,CandData,MyPeriod,Unit,Letty,ColourMapChoice,IsLand)
		
#    stop()

print("And, we are done!")

