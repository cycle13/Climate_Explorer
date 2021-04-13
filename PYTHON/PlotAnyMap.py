#!/usr/local/sci/bin/python
# PYTHON23
# 
# Author: Kate Willett
# Created: 22 June 2015
# Last update: 15 Feb 2021
# Location: /home/h04/hadkw/HadISDH_Code/CLIMEXPLORER/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This reads in a netCDF file of a gridded dataset and plots a map of either a single field or a processed field (average etc)
# It can renormalise data to the desired climatology period (and deal with actual monthly means)
# It can plot: 
#	a single month, 
#	an average of months within a year (or adjacent for DJF) up to annual - set minimum data presence
#	an average of single months across a period of years (climatology) - set minimum data presence
# 	an average of several months across a period of years (climatology) up to annual - set minimum data presence
# It also saves a netCDF of the  plot for replotting if desired.
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
# import cartopy.crs as ccrs # NO MORE BASEMAP
# import cartopy.feature as cpf
# import datetime as dt
# from matplotlib.dates import date2num,num2date
# from scipy.io import netcdf
# import matplotlib.colors as mc
# import matplotlib.cm as mpl_cm
# import pdb # pdb.set_trace() or c
#
# Kate's:
# from ReadNetCDF import GetGrid - written by Kate Willett, reads in any netCDF grid, can cope with multiple fields
# from ReformatGrids import GetAnomalies - written by Kate Willett, converts gridded fields to anomalies or renormalises to a different period
# import SelectSlice - written by Kate Willett, takes a single month or average of months/years for each gridbox
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
#     StYr = 1973 # Start Year of dataset (assumes start in January)
#     EdYr = 2014 # End Year of dataset
#     EdMon = 10 # End Month of dataset (1 to 12) 		
#     # Are you reading in actuals or anomalies? 'anoms' for anomalies, 'abs' for actuals
#     vartype = 'anoms' 
#     # Choose your start year of climatology: 0 if not relevant, not important if DoReClim = False
#     ChooseClimSt = 1976 
#     # Choose your start year of climatology: 0 if relevant, not important if DoReClim = False
#     ChooseClimEd = 2005 
#     # Does the data need renormalising to the chosen climatology? True if need to renormalise, False if not
#     DoReClim = False 
#     # Are we plotting anomalies or absolute values? 'actual' for actual values, 'anomaly' for anomaly values
#     PlotType = 'actual' 
#     # Are we saving the plot as a netCDF too?
#     SaveData = True

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
#     # LIST Variable name used in netCDF file, based on vartype
#     NameOfVar = ['q_'+vartype]
#     # Long name of variable for plot
#     VarName = 'Specific Humidity'
#     # Unit for variable - used for plot
#     Unit = 'g kg$^{-1}$'  
#     # Plot label
#     Namey = 'HadISDH.landq.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#     # Latitude info: LIST of variable name OR number of latitudes, start latitude
#     LatInfo = ['latitude'] # [36,-87.5] 
#     # Longitude info: LIST of variable name OR number of longitudes, start longitude
#     LonInfo = ['longitude'] # [72,-177.5]) 
#     # Choice of colourmap for plot and whether it should be flipped: 'flip','noflip'
#     ColourMapChoice = ('BrBG','noflip')
#     # Optional hardwire info for plot - can be 0,0,0,' ' - use if wanting to make identical colour bars
#     PlotInfo = [0,28,14,' '] # vmin,vmax,nsteps,letter_for_plot_label 

# Then from the command line select var, month (range) and year (range):
# --var q (rh, t, td, tw, e, dpd
# --mchoice 0 (0 to 11) 
# --mchoiceend '' (0 to 11 if you want to plot an average over months or each individual month in a range)
# --mmult False (True for every month in range, False for the average)
# --ychoice 2020 (a single year) 
# --ychoiceend '' (a single year if you want to plot an average over years or each individual year in a range)
# --ymult False (True for every year in range, False for the average)
#
# > module load scitools/default-current
# > python PlotAnyMap_JUN2015.py --var q --mchoice 0 --mchoiceend '' --mmult False --ychoice 2020 --ychoiceend '' --ymult False
#
# Or use ./submit_spice_PlotAnyMap.bash
#
# -----------------------
# OUTPUT
# -----------------------
# An .eps and .png map
# OUTDIR = '/scratch/hadkw/UPDATE2020/IMAGES/MAPS/'
# OUTDIR = '/data/users/hadkw/WORKING_HADISDH/UPDATE2020/IMAGES/MAPS/'
# Example output filename built from run choices:
# OUTPLOT = 'Map_HadISDH.landq.4.3.0.2020f_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#
# A .nc
# OUTDIR = '/scratch/hadkw/UPDATE2020/STATISTICS/GRIDS/'
# OUTDIR = '/data/users/hadkw/WORKING_HADISDH/UPDATE2020/STATISTICS/GRIDS/'
# Example output filename built from run choices:
# OUTPLOT = 'Map_HadISDH.landq.4.3.0.2020f_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 3 (15th Feb 2015)
# ---------
#  
# Enhancements
# Now call var month and year choices from command line so that it can run on spice in batch mode
# Use with 'submit_apice_PlotAnyMap.bash'
# Now also has SaveData=True option to output a netCDF file - replaces function of old IDL code 
#  
# Changes
#  
# Bug fixes
#  
#
# Version 2 (20th July 2020)
# ---------
#  
# Enhancements
# Now Python 3 rather than 2.7
#  
# Changes
#  
# Bug fixes
#  
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
import sys, os, getopt
import scipy.stats
import struct
import cartopy.crs as ccrs # NO MORE BASEMAP
import cartopy.feature as cpf
import datetime as dt
from matplotlib.dates import date2num,num2date
#from scipy.io import netcdf
from netCDF4 import Dataset
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import pdb # pdb.set_trace() or c

from ReadNetCDF import GetGrid4
from ReformatGrids import GetAnomalies
from SelectSlice import SelectSlice

# Set in stone things:
MonArr = ['January','February','March','April','May','June','July','August','September','October','November','December']
VarDict = dict([('q',['q','g kg$^{-1}$','Specific Humidity',('BrBG','noflip'),dict([('MinVal',-2.),('MaxVal',2.),('StepVal',9.),('LetterVal','')])]),
                ('rh',['RH','%rh','Relative Humidity',('BrBG','noflip'),dict([('MinVal',-12.),('MaxVal',12.),('StepVal',9.),('LetterVal','')])]),
		('e',['e','hPa','Vapour Pressure',('BrBG','noflip'),dict([('MinVal',-2.),('MaxVal',2.),('StepVal',9.),('LetterVal','')])]),
		('t',['T','$^{o}$C','Air Temperature',('coolwarm','noflip'),dict([('MinVal',-2.5),('MaxVal',2.5),('StepVal',9.),('LetterVal','')])]),
		('tw',['Tw','$^{o}$C','Wetbulb Temperature',('BrBG','noflip'),dict([('MinVal',-2.),('MaxVal',2.),('StepVal',9.),('LetterVal','')])]),
		('td',['Td','$^{o}$C','Dew Point Temperature',('BrBG','noflip'),dict([('MinVal',-2.),('MaxVal',2.),('StepVal',9.),('LetterVal','')])]),
		('dpd',['DPD','$^{o}$C','Dew Point Depression',('BrBG','flip'),dict([('MinVal',-3.),('MaxVal',3.),('StepVal',9.),('LetterVal','')])])])

# EDITABLES!!!

# Domain if HadISDH
Domain = 'land' # 'land', 'marine',' blend'
#Domain = 'blend' # 'land', 'marine',' blend'
#Domain = 'marine' # 'land', 'marine',' blend'

# Version
lversion = '4.3.1.2020f'
lExpType = 'FLATgridHOM5by5'
mversion = '1.1.0.2020f'
mExpType = 'BClocalSHIP5by5both'
bversion = '1.1.1.2020f'
bExpType = 'FLATgridHOMBClocalSHIPboth5by5'
    
# Working Year (e.g. for HadISDH this would correspond to UPDATE<YYYY>
workingyear = '2020'

# Missing data
mdi = -1e30 # may set up as masked arrays later

# Set up initial run choices
# What is the minimum percent of data present for each gridbox value? This can be any value between >0 and 1
PctDataPresent = 0.6	

# Start Year of dataset (assumes January)
StYr = 1973 		

# End Year of dataset
EdYr = 2020 		

# End Month of dataset (assumes it starts in January) - 1 to 12
EdMon = 12 		

# Are you reading in actuals or anomalies? 'anoms' or 'anomalies' for anomalies, 'abs' for actuals
vartype = 'anoms' # Check PlotType = 'actuals' if we want to plot actuals
#vartype = 'abs' # Check PlotType = 'actuals' if we want to plot actuals

# Choose your start year of climatology: 0 if not relevant, not important if DoReClim = False
ChooseClimSt = 1981 

# Choose your start year of climatology: 0 if relevant, not important if DoReClim = False
ChooseClimEd = 2010 

# Does the data need renormalising to the chosen climatology? True if need to renormalise, False if not
DoReClim = False 

# Are we plotting anomalies or absolute values? 'actual' for actual values, 'anomaly' for anomaly values
PlotType = 'anomaly' # Check RangeDict is set to 0,0,0,'?' or it will default to anomalies ranges that have been hard wired.
#PlotType = 'actual' # Check RangeDict is set to 0,0,0,'?' or it will default to anomalies ranges that have been hard wired.

# Are we saving the data for te plot?
SaveData = True # True for save to netCDF, False for plot only

# Set up forced vmin,vmax,nsteps and plotlabel letter if needed or leave as default
RangeDict = dict([('MinVal',0.),('MaxVal',0.),('StepVal',0.),('LetterVal','c)')]) # default = 0., 0. ,0., '' 
#Fix vals ,
#RangeDict = VarDict[Var][4]

# Set up directories:
INDIRC = '/scratch/hadkw/UPDATE'+workingyear+'/STATISTICS/GRIDS/'
#INDIRC = '/data/users/hadkw/WORKING_HADISDH/UPDATE'+workingyear+'/STATISTICS/GRIDS/'
#INDIRC = '/data/local/hadkw/HADCRUH2/UPDATE'+workingyear+'/OTHERDATA/'
#OUTDIRP = '/data/users/hadkw/WORKING_HADISDH/UPDATE'+workingyear+'/IMAGES/MAPS/'
OUTDIRP = '/scratch/hadkw/UPDATE'+workingyear+'/IMAGES/MAPS/'
#OUTDIRD = '/data/users/hadkw/WORKING_HADISDH/UPDATE'+workingyear+'/STATISTICS/GRIDS/'
OUTDIRD = '/scratch/hadkw/UPDATE'+workingyear+'/STATISTICS/GRIDS/'

#DataSetShort = 'BERKELEY_T_'
#DataSet = 'BERKELEY_T_'
#candidate = 'BERKELEY_T_5by519762005clim_anoms_19732014'
#NameOfVar = ['t_'+vartype]
#LatInfo = ['latitude']
#LonInfo = ['longitude']

#DataSetShort = 'BERKELEY_T_'
#DataSet = 'BERKELEY_T_'
#candidate = 'GISS_T_5by519762005clim_anoms_19732014'
#NameOfVar = ['t_'+vartype]
#LatInfo = ['latitude']
#LonInfo = ['longitude']

#DataSetShort = 'BERKELEY_T_'
#DataSet = 'BERKELEY_T_'
#candidate = 'CRUTEM.4.3.0.0.anomalies'
#NameOfVar = ['t_'+vartype]
#LatInfo = ['latitude']
#LonInfo = ['longitude']

#DataSetShort = 'BERKELEY_T_'
#DataSet = 'BERKELEY_T_'
#candidate = 'HadCRUT.4.4.0.0.median'
#NameOfVar = ['temperature_anomaly']
#LatInfo = ['latitude']
#LonInfo = ['longitude']

#DataSet = 'BERKELEY_T_'
#candidate = 'GHCNM_18802014'
#NameOfVar = ['t_'+vartype]
#LatInfo = ['latitude']
#LonInfo = ['longitude']

#************************************************************************
# Subroutines
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
    plt1=plt.axes([0.05,0.12,0.9,0.8],projection=ccrs.Robinson()) # left, bottom, width, height
    
    plt1.coastlines()
    #plt1.set_boundary # not sure what this does? maybe useful for plotting regions?
    plt1.gridlines(draw_labels=False) # probably a way to specify these exactly if we want something different to default
    # This line background fills the land with light grey
    #plt1.add_feature(cpf.LAND, zorder = 0, facecolor = "0.9", edgecolor = "k") # may or may not need this
    
    ext = plt1.get_extent()
    # End of CARTOPY

    # make up a blue to red (reverse) colour map
    cmap=plt.get_cmap(TheCmap[0])     
    
    cmaplist=[cmap(i) for i in range(cmap.N)]
    #pdb.set_trace()
    for loo in range(np.int((cmap.N/2)-20),np.int((cmap.N/2)+20)):
        cmaplist.remove(cmaplist[np.int((cmap.N/2)-20)]) # remove the very pale colours in the middle
    if (TheCmap[1] == 'flip'):	# then reverse the colours
        cmaplist.reverse()
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)

    # Is there isn't a hardwired vmin and vmax then make one that is nice
    # nsteps always 9 to avoid toomany colours
    nsteps = 9
    if (ThePlotInfo[2] == 0):
        # work out best max and min values for colourbar and try to make them 'nice
        # must be an odd number of steps
        # for less than 0.1 must be 14 or fewer steps
        if (TheTypee == 'anomaly'):
            vmax=np.int(np.ceil(np.max(abs(MSKTheCandData))*10))/10.    
            vmin=-vmax
            print(vmin,vmax)
#            if (vmax <= 0.3):
#                nsteps=np.int((vmax-vmin)/0.05)+1
#            elif (vmax <= 0.5):
#                vmax=np.ceil(np.max(abs(MSKTheCandData))/0.06)*0.06  
#                vmin=-vmax
#                nsteps=np.int((vmax-vmin)/0.06)+1
#            elif (vmax <= 1.0):
#                vmax=np.ceil(np.max(abs(MSKTheCandData))/0.1)*0.1  
#                vmin=-vmax
#                nsteps=np.int((vmax-vmin)/0.1)+1
#            elif (vmax <= 1.4):
#                vmax=np.ceil(np.max(abs(MSKTheCandData))/0.2)*0.2  
#                vmin=-vmax
#                nsteps=np.int((vmax-vmin)/0.2)+1
#            elif (vmax > 1.4):
#                vmax=np.ceil(np.max(abs(MSKTheCandData))/0.3)*0.3  
#                vmin=-vmax
#                nsteps=np.int((vmax-vmin)/0.3)+1
	    #    pdb.set_trace() # stop here and play
        else:
#            vmax=np.ceil(np.max(abs(MSKTheCandData))) 
#            vmin=np.floor(np.min(abs(MSKTheCandData)))
            vmax=np.ceil(np.max(MSKTheCandData)) 
            vmin=np.floor(np.min(MSKTheCandData))
            vrange=vmax-vmin
            print(vmin,vmax,vrange)
#            if (vmax-vmin <  14):
#                nsteps=np.int((vmax-vmin))+1
#            elif (vmax <= 21):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(1.5*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(1.5*7) 
#                nsteps=np.int((vmax-vmin)/1.5)+1
#            elif (vmax <= 30):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(2*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(2*7) 
#                nsteps=np.int((vmax-vmin)/2)+1
#            elif (vmax <= 35):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(2.5*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(2.5*7) 
#                nsteps=np.int((vmax-vmin)/2.5)+1
#            elif (vmax <= 42):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(3*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(3*7) 
#                nsteps=np.int((vmax-vmin)/3)+1
#            elif (vmax <= 56):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(4*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(4*7) 
#                nsteps=np.int((vmax-vmin)/4)+1
#            elif (vmax <= 70):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(5*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(5*7) 
#                nsteps=np.int((vmax-vmin)/5)+1
#            elif (vmax <= 84):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(6*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(6*7) 
#                nsteps=np.int((vmax-vmin)/6)+1
#            elif (vmax <= 98):
#                vmax=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))+(7*7) 
#                vmin=(np.floor(np.min(abs(MSKTheCandData)))+(vrange/2.))-(7*7) 
#                nsteps=np.int((vmax-vmin)/7)+1
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
    grids = plt1.pcolor(LngArrLons,LngArrLats,MSKTheCandData,transform = ccrs.PlateCarree(),cmap=cmap,norm=norm)
#    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCandData,cmap=cmap,norm=norm,latlon='TRUE')

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
    plt.clf()
     
    return #PlotAnyMap

#**********************************************************************
# WriteNetCDF
def WriteNetCDF(Filename,TheGrids,TheMDI,TheVar,TheLongName,TheUnit,TheLats,TheLons):
    '''
    This function writes out a NetCDF 4 file
    
    INPUTS:
    Filename - string file name
    TheGrids[:,:] - 2D array of decadal average trends
    TheMDI - the missind data value
    TheVar - string short name of var q
    TheLongName - strong long name of variable
    TheUnit - string unit of variable
    TheLats[:] - vector of latitudes from -90 to 90
    TheLons[:] - vector of longitudes from -180 to 180
    OUTPUTS:
    None
    
    '''
    
    # No need to convert float data using given scale_factor and add_offset to integers - done within writing program (packV = (V-offset)/scale
    # Not sure what this does to float precision though...

    # Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
    ncfw = Dataset(Filename,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first

    # Set up the dimension names and quantities
    ncfw.createDimension('latitude',len(TheLats))
    ncfw.createDimension('longitude',len(TheLons))
    
    # Go through each dimension and set up the variable and attributes for that dimension if needed
    MyVarLt = ncfw.createVariable('latitude','f4',('latitude',))
    MyVarLt.standard_name = 'latitude'
    MyVarLt.long_name = 'latitude'
    MyVarLt.units = 'degrees_north'
#    MyVarLt.valid_min = np.min(TheLats)
#    MyVarLt.valid_max = np.max(TheLats)
    MyVarLt.point_spacing = 'even'
    MyVarLt.axis = 'X'
    MyVarLt[:] = TheLats

    MyVarLn = ncfw.createVariable('longitude','f4',('longitude',))
    MyVarLn.standard_name = 'longitude'
    MyVarLn.long_name = 'longitude'
    MyVarLn.units = 'degrees east'
#    MyVarLn.valid_min = np.min(TheLons)
#    MyVarLn.valid_max = np.max(TheLons)
    MyVarLn.point_spacing = 'even'
    MyVarLn.axis = 'X'
    MyVarLn[:] = TheLons

    # Go through each variable and set up the variable attributes
    # I've added zlib=True so that the file is in compressed form
    # I've added least_significant_digit=4 because we do not need to store information beyone 4 significant figures.
    MyVar = ncfw.createVariable('anomalies','f4',('latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
#    MyVar.standard_name = TheStandardName
    MyVar.long_name = TheLongName
    MyVar.units = TheUnit
#    MyVar.valid_min = np.min(TheTrendGrids)
#    MyVar.valid_max = np.max(TheTrendGrids)
#    MyVar.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVar[:] = TheGrids[:,:]  	

    ncfw.close()

    return
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
def main(argv):
    # INPUT PARAMETERS AS STRINGS!!!!
    var = 'q'	    # 'q','rh','e','td','tw','t','dpd'
    mchoice = 0 # 0 to 11 this is either the month you want plotted or start of the month range to plot (can do 11-1)
    mchoiceend = '' # '' for no range of months (single) or end month of rande 0-11 (can do 11-1)
    mmult = False # False for an average, True for plotting each individual month of range
    ychoice = 1973 # 0this is either the year you want plotted or start of the year range to plot - for 11-1 need year of 11
    ychoiceend = '' # '' for no range of years (single) or end year of range
    ymult = False # False for an average, True for plotting each individual year of range
    
    try:
        opts, args = getopt.getopt(argv, "hi:",
	                           ["var=","mchoice=","mchoiceend=","mmult=","ychoice=","ychoiceend=","ymult="])
    except getopt.GetoptError:
        print('Usage (as strings) PlotAnyMap.py --var <q> --mchoice <0> --mchoiceend <11> --mmult False --ychoice <1973> --ychoiceend <''> --ymult True')
        sys.exit(2)

    for opt, arg in opts:
        if opt == "--var":
            try:
                Var = arg
            except:
                sys.exit("Failed: var not a string")
        elif opt == "--mchoice":
            try:
                Mchoice = arg
            except:
                sys.exit("Failed: mchoice not an integer")
        elif opt == "--mchoiceend":
            try:
                Mchoiceend = arg
            except:
                sys.exit("Failed: mchoiceend not an integer")
        elif opt == "--mmult":
            try:
                Mmult = arg
            except:
                sys.exit("Failed: mmult not a boolean")
        elif opt == "--ychoice":
            try:
                Ychoice = arg
            except:
                sys.exit("Failed: ychoice not an integer")
        elif opt == "--ychoiceend":
            try:
                Ychoiceend = arg
            except:
                sys.exit("Failed: ychoiceend not an integer")
        elif opt == "--ymult":
            try:
                Ymult = arg
            except:
                sys.exit("Failed: ymult not a boolean")

#    assert year1 != -999 and year2 != -999, "Year not specified."

    print(Var,Mchoice,Mchoiceend,Mmult,Ychoice,Ychoiceend,Ymult)
    # Set argparse so that key arguments can be passed to this program externally rather than being set above

    # Version
    if (Domain == 'land'):
        Version = lversion
        ExpType = lExpType
    elif (Domain == 'marine'):    
        Version = mversion
        ExpType = mExpType
    elif (Domain == 'blend'):
        Version = bversion
        ExpType = bExpType

    # Select your month of choice, or a range for an average, or a range for looping through
    # 0...11 represent Jan..Dec, [2,4] for Mar-Apr-May average, [0,11] for annual average, [11,1] for Dec-Jan-Feb average
    # For month ranges than span 11 to 0, December will be taken from the first year of ChooseYr - will NOT work for last year!
    if (Mchoiceend != ''):
        ChooseMon = [int(Mchoice),int(Mchoiceend)]
    else:
        ChooseMon = [int(Mchoice)]
    if (Mmult == 'True'):
        PlotMonMultiple = True # False to plot an average over range, True to plot each individual month in range
    else:
        PlotMonMultiple = False # False to plot an average over range, True to plot each individual month in range
    
    # Select your year of choice, or a range for an average, or a range for looping through
    # 1973...2014 for individual years, [1973,1982] for decadal average etc
    if (Ychoiceend != ''):
        ChooseYr = [int(Ychoice),int(Ychoiceend)]
    else:
        ChooseYr = [int(Ychoice)]
    if (Ymult == 'True'):
        PlotYrMultiple = True # False to plot an average over range, True to plot each individual month in range
    else:
        PlotYrMultiple = False # False to plot an average over range, True to plot each individual month in range
        
    # Set up files and variable bundles:
    # For HadISDH the next bits are automatically populated or uncomment a bundle if its a different dataset
    DataSetShort = 'HadISDH.'+Domain+VarDict[Var][0]+'.'+Version
    DataSet = DataSetShort+'_'+ExpType
    candidate = DataSet+'_anoms8110'
    NameOfVar = [Var+'_'+vartype]
    LatInfo = ['latitude'] # variable name or number of latitudes, start latitude
    LonInfo = ['longitude']  # variable name or number of longitudes, start longitude

    # Set up other variables
    VarName = VarDict[Var][2]
    Unit = VarDict[Var][1]  
    ColourMapChoice = VarDict[Var][3]
    PlotInfo = [RangeDict['MinVal'], RangeDict['MaxVal'], RangeDict['StepVal'], RangeDict['LetterVal']]
    NYrs = (EdYr+1)-StYr
    NMons = NYrs*12-(12-EdMon)

    # read in trend maps
    MyFile = INDIRC+candidate+'.nc'
    #print(NameOfVar,LatInfo,LonInfo)
    TmpCandFields,LatList,LonList = GetGrid4(MyFile,NameOfVar,LatInfo,LonInfo)

    # If the data do not end in December then pad the file with missing data
    if (EdMon < 12):
        #pdb.set_trace()
        TmpCandFields = np.concatenate((TmpCandFields,np.reshape(np.repeat(mdi,((12-EdMon)*len(LonList)*len(LatList))),(12-EdMon,len(LatList),len(LonList)))))
        #pdb.set_trace()

    # force mdi to equal mdi without floating point errors
    print(len(TmpCandFields[TmpCandFields < mdi]))
    TmpCandFields[TmpCandFields < mdi] = mdi
    print(len(TmpCandFields[TmpCandFields < mdi]))
    #pdb.set_trace()

    # Do you need to renormalise (or create anomalies from actuals?
    if (DoReClim):
        CandFields = GetAnomalies(TmpCandFields,StYr,EdYr,ChooseClimSt,ChooseClimEd,mdi,PctDataPresent)
    else:
        CandFields = TmpCandFields    

    # Now work on plots either singular or multiple
    # If we're looping through years then start loop
    if (PlotYrMultiple):

        for yy in range(ChooseYr[0],ChooseYr[1]+1): # needs extra to include last month within the range

            YrStr = str(yy)
	    
	    # If we're looping through months then start loop
            if (PlotMonMultiple):
	
                for mm in range(ChooseMon[0],ChooseMon[1]+1): # needs extra to include last month within the range
	
                    MonStr = MonArr[mm]

                    # Extract chosen month
                    CandData = SelectSlice(CandFields,StYr,EdYr,[mm],[yy],mdi,PctDataPresent)

                    # pass to plotter
                    OUTPLOT = 'Map_'+DataSet+'_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
                    Namey = DataSetShort+' '+PlotType+' '+MonStr+' '+YrStr
                    MyFile = OUTDIRP+OUTPLOT
                    PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,ColourMapChoice,PlotType,VarName,mdi,PlotInfo)
                    if (SaveData):
                        WriteNetCDF(OUTDIRD+OUTPLOT+'.nc',CandData,mdi,Var,VarName,Unit,LatList,LonList)

	    # If we're not then work on the individual or average
            else:
       
                if (len(ChooseMon) == 1):
               
                    MonStr = MonArr[ChooseMon[0]]
    
                else:
               
                    MonStr = MonArr[ChooseMon[0]]+'-'+MonArr[ChooseMon[1]]
		
                # Extract chosen month
                CandData = SelectSlice(CandFields,StYr,EdYr,ChooseMon,[yy],mdi,PctDataPresent)

                # pass to plotter
                OUTPLOT = 'Map_'+DataSet+'_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
                Namey = DataSetShort+' '+PlotType+' '+MonStr+' '+YrStr
                MyFile = OUTDIRP+OUTPLOT
                PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,ColourMapChoice,PlotType,VarName,mdi,PlotInfo)
                if (SaveData):
                    WriteNetCDF(OUTDIRD+OUTPLOT+'.nc',CandData,mdi,Var,VarName,Unit,LatList,LonList)
    
    
    # If we're not looping through years then check month multiples
    else:
        
	# If we're looping through months then start loop
        if (PlotMonMultiple):
	
            for mm in range(ChooseMon[0],ChooseMon[1]+1): # needs extra to include last month within the range
	
                MonStr = MonArr[mm]

                if (len(ChooseYr) == 1):
    
                    YrStr = str(ChooseYr[0])
   
                else:
    
                    YrStr = str(ChooseYr[0])+'-'+str(ChooseYr[1])

                # Extract chosen month
                CandData = SelectSlice(CandFields,StYr,EdYr,[mm],ChooseYr,mdi,PctDataPresent)

                # pass to plotter
                OUTPLOT = 'Map_'+DataSet+'_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
                Namey = DataSetShort+' '+PlotType+' '+MonStr+' '+YrStr
                MyFile = OUTDIRP+OUTPLOT
                PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,ColourMapChoice,PlotType,VarName,mdi,PlotInfo)
                if (SaveData):
                    WriteNetCDF(OUTDIRD+OUTPLOT+'.nc',CandData,mdi,Var,VarName,Unit,LatList,LonList)
		
	# If we're not then work on the individual or average
        else:
       
            if (len(ChooseMon) == 1):
            
                MonStr = MonArr[ChooseMon[0]]
    
            else:
               
                MonStr = MonArr[ChooseMon[0]]+'-'+MonArr[ChooseMon[1]]
		
            if (len(ChooseYr) == 1):
    
                YrStr = str(ChooseYr[0])
   
            else:
    
                YrStr = str(ChooseYr[0])+'-'+str(ChooseYr[1])

            # Extract chosen month
            CandData = SelectSlice(CandFields,StYr,EdYr,ChooseMon,ChooseYr,mdi,PctDataPresent)
            #pdb.set_trace()
            # pass to plotter
            OUTPLOT = 'Map_'+DataSet+'_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#            Namey = DataSetShort+' '+PlotType+' '+MonStr+' '+YrStr
            Namey = ''
            MyFile = OUTDIRP+OUTPLOT
            PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,ColourMapChoice,PlotType,VarName,mdi,PlotInfo)
            if (SaveData):
                WriteNetCDF(OUTDIRD+OUTPLOT+'.nc',CandData,mdi,Var,VarName,Unit,LatList,LonList)

    #pdb.set_trace()
    
    print("And, we are done!")

if __name__ == '__main__':
    
    main(sys.argv[1:])

