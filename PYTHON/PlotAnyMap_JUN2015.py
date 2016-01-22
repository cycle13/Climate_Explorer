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

from ReadNetCDF import GetGrid
from ReformatGrids import GetAnomalies
from SelectSlice import SelectSlice

# Set in stone things:
MonArr = ['January','February','March','April','May','June','July','August','September','October','November','December']

# Missing data
mdi = -1e30 # may set up as masked arrays later

# Set up initial run choices
# What is the minimum percent of data present for each gridbox value? This can be any value between >0 and 1
PctDataPresent = 0.6	

# Start Year of dataset (assumes January)
StYr = 1850 		

# End Year of dataset
EdYr = 2015 		

# End Month of dataset (assumes it starts in January) - 1 to 12
EdMon = 10 		

# Are you reading in actuals or anomalies? 'anoms' or 'anomalies' for anomalies, 'abs' for actuals
vartype = 'anoms' 

# Select your month of choice, or a range for an average
# 0...11 represent Jan..Dec, [2,4] for Mar-Apr-May average, [0,11] for annual average, [11,1] for Dec-Jan-Feb average
# For month ranges than span 11 to 0, December will be taken from the first year of ChooseYr - will NOT work for last year!
ChooseMon = [0,9] 

# Select your year of choice, or a range for an average
# 1973...2014 for individual years, [1973,1982] for decadal average etc
ChooseYr = [2015,2015] 

# Choose your start year of climatology: 0 if not relevant, not important if DoReClim = False
ChooseClimSt = 1981 

# Choose your start year of climatology: 0 if relevant, not important if DoReClim = False
ChooseClimEd = 2010 

# Does the data need renormalising to the chosen climatology? True if need to renormalise, False if not
DoReClim = True 

# Are we plotting anomalies or absolute values? 'actual' for actual values, 'anomaly' for anomaly values
PlotType = 'anomaly' 

# From these choices now work out output file string and plot title label
if (len(ChooseMon) > 1):
    MonStr = MonArr[ChooseMon[0]]+'-'+MonArr[ChooseMon[1]]
else:
    MonStr = MonArr[ChooseMon[0]]
if (len(ChooseYr) > 1):
    YrStr = str(ChooseYr[0])+'-'+str(ChooseYr[1])
else:
    YrStr = str(ChooseYr[0])

# Set up directories:
#INDIRC = '/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
INDIRC = '/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
OUTDIR = '/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'

# Set up files and variable bundles:
#candidate = 'HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landq.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['q_'+vartype]
#VarName = 'Specific Humidity'
#Unit = 'g kg$^{-1}$'  
#Namey = 'HadISDH.landq.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude'] # variable name or number of latitudes, start latitude
#LonInfo = ['longitude']  # variable name or number of longitudes, start longitude
#ColourMapChoice = ('BrBG','noflip')
#PlotInfo = [0,28,15,' '] # vmin,vmax,nsteps+1,letter_for_plot_label 

#candidate = 'HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landRH.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['rh_'+vartype]
#VarName = 'Relative Humidity'
#Unit ='%rh'  
#Namey = 'HadISDH.landRH.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('BrBG','noflip')
#PlotInfo = [1,99,15,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landT.2.0.1.2014p_RAW_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['t_'+vartype]
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'HadISDH.landT.2.0.1.2014p RAW '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#OUTPLOT = 'Map_HadISDH.landT.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['t_'+vartype]
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'HadISDH.landT.2.0.1.2014p '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'BERKELEY_T_5by519762005clim_anoms_19732014'
#OUTPLOT = 'Map_BERKELEY_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['t_'+vartype]
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'BERKELEY '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'GISS_T_5by519762005clim_anoms_19732014'
#OUTPLOT = 'Map_GISS_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['t_'+vartype]
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'GISS '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'CRUTEM.4.3.0.0.anomalies'
#OUTPLOT = 'Map_CRUTEM4.3.0.0_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['t_'+vartype]
#VarName = 'Temperature'
#Unit = '$^{o}$C' 
#Namey = 'CRUTEM4.3.0.0 '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps,letter_for_plot_label 

candidate = 'HadCRUT.4.4.0.0.median'
OUTPLOT = 'Map_HadCRUT4.4.0.0_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
NameOfVar = ['temperature_anomaly']
VarName = 'Temperature'
Unit = '$^{o}$C' 
Namey = 'HadCRUT.4.4.0.0 '+PlotType+' '+MonStr+' '+YrStr
LatInfo = ['latitude']
LonInfo = ['longitude']
ColourMapChoice = ('coolwarm','noflip')
PlotInfo = [-3.,3.,13,' '] # vmin,vmax,nsteps,letter_for_plot_label 

#candidate = 'GHCNM_18802014'
#OUTPLOT = 'Map_GHCNM3_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+PlotType+MonStr+YrStr
#NameOfVar = ['t_'+vartype]
#VarName = 'Temperature'
#Unit = '$^{o}$C'  
#Namey = 'GHCNM3 '+PlotType+' '+MonStr+' '+YrStr
#LatInfo = ['latitude']
#LonInfo = ['longitude']
#ColourMapChoice = ('coolwarm','noflip')
#PlotInfo = [0,0,0,' '] # vmin,vmax,nsteps(+1),letter_for_plot_label 

# Set up other variables
NYrs = (EdYr+1)-StYr
NMons = NYrs*12-(12-EdMon)
    
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
if (__name__ == "__main__"):
    
    # Set argparse so that key arguments can be passed to this program externally rather than being set above
    
    # read in trend maps
    MyFile=INDIRC+candidate+'.nc'
    #print(NameOfVar,LatInfo,LonInfo)
    TmpCandFields,LatList,LonList = GetGrid(MyFile,NameOfVar,LatInfo,LonInfo)

    # If the data do not end in January then pad the file with missing data
    if (EdMon < 12):
        #pdb.set_trace()
	TmpCandFields=np.concatenate((TmpCandFields,np.reshape(np.repeat(mdi,((12-EdMon)*len(LonList)*len(LatList))),(12-EdMon,len(LatList),len(LonList)))))
        #pdb.set_trace()

    # force mdi to equal mdi without floating point erros
    print(len(TmpCandFields[TmpCandFields < mdi]))
    TmpCandFields[TmpCandFields < mdi]=mdi
    print(len(TmpCandFields[TmpCandFields < mdi]))
    #pdb.set_trace()

    # Do you need to renormalise (or create anomalies from actuals?
    if (DoReClim):
        CandFields = GetAnomalies(TmpCandFields,StYr,EdYr,ChooseClimSt,ChooseClimEd,mdi,PctDataPresent)
    else:
        CandFields = TmpCandFields    

    # Extract chosen month
    CandData = SelectSlice(CandFields,StYr,EdYr,ChooseMon,ChooseYr,mdi,PctDataPresent)

    # pass to plotter
    MyFile=OUTDIR+OUTPLOT
    PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,ColourMapChoice,PlotType,VarName,mdi,PlotInfo)
		
    # stop()
    #pdb.set_trace()
    
    print("And, we are done!")

