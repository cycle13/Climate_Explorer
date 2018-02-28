#!/usr/local/sci/bin/python

#***************************************
# 22 March 2014 KMW - v1
# Plots station locations for each variable including removed sub/super sats  
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotStationLocs_MAR2014.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import sys, os
import scipy.stats
import struct
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from netCDF4 import Dataset
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LATITUDE_FORMATTER
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER


# RESTART VALUE
Restarter='------'				#'------'		#'681040'

# Set up initial run choices
param='q'	# tw, q, e, rh, t, td, dpd
param2='q'	# Tw, q, e, RH, T, Td, DPD
homogtype='IDPHA'	# 'IDPHA','PHA','PHADPD'
nowmon='OCT'
nowyear='2014'
thenmon='JAN'
thenyear='2014'
version='2.0.0.2013p'

# Set up directories and files
INDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/LISTS_DOCS/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/'

#INLIST='Posthomog'+homogtype+param+'_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INLIST='USA_post2005dropoff_OCT2014'
OUTPLOT='StationLocations_'+INLIST+nowmon+nowyear

PlotTitle='Post-2005 drop off stations in the USA'
MapBoundary=[-170.,20.,-60.,75.0]

# Set up variables
ngoods=0	#set once file read in
nsats=0		#set once file read in
nsubs=0		#set once file read in

StWMOs=[]
StLats=[]
StLons=[]

#************************************************************************
# Subroutines
#************************************************************************
# READDATA
def ReadData(FileName,typee,delimee):
    ''' Use numpy genfromtxt reading to read in all rows from a complex array '''
    ''' Need to specify format as it is complex '''
    ''' outputs an array of tuples that in turn need to be subscripted by their names defaults f0...f8 '''
    return np.genfromtxt(FileName, dtype=typee,delimiter=delimee) # ReadData

#************************************************************************
# PlotNiceMap
def PlotNiceMap(TheFile,TheGLats,TheGLons,TheGCount,TheTitle,TheBoundaries):
    ''' Plot all good lats and lons TheGLats, TheGLons '''
    ''' If there are sats (TheSTCount > 0) plot all sat lats TheSTLats, TheSTLons '''
    ''' If there are subzeros (TheSBCount > 0) plot all sat lats TheSBLats, TheSBLons '''
    ''' Label plot with totals '''
    ''' Save as png and eps '''
      
    # set up plot
 
    fig=plt.figure(1,figsize=(8,5))
    plt.clf()
    
    ax=plt.axes([0.08,0.05,0.84,0.88],projection=ccrs.PlateCarree(central_longitude=0),) # LEFT BOTTOM WIDTH HEIGH
    ax.set_extent([TheBoundaries[0],TheBoundaries[2],TheBoundaries[1],TheBoundaries[3]]) 
    ax.coastlines()

    gl = ax.gridlines(draw_labels=True)
    #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabels_top = False
    gl.xlabel_style = {'size': 16,'color': 'black'} # 'weight': 'bold'
    gl.ylabel_style = {'size': 16,'color': 'black'}
    
    ax.scatter(np.array(TheGLons),np.array(TheGLats),c='R',s=50,transform=ccrs.Geodetic())        
#    ax.plot(TheGLons,TheGLats,color='Red',linestyle=' ',markersize=10,transform=ccrs.Geodetic())

#    plt.tight_layout()
        
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
#    plt.figtext(0.05,0.9,TheLetter,size=16)
    plt.figtext(0.5,0.92,TheTitle,size=16,ha='center')

#    raw_input("stop")	# REALLY USEFUL TO INTERACT WITHIN SUBROUTINE ctrl C
    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotNiceMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list
MyTypes=("|S6","|S5","float","float","float","|S4","|S30","|S7")
MyDelimiters=[6,5,8,10,7,4,30,2]
MyFile=INDIR+INLIST+'.txt'
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
StWMOs=np.array(RawData['f0'])
StLats=np.array(RawData['f2'])
StLons=np.array(RawData['f3'])
nstats=len(StWMOs)


# pass to plotter
    
MyFile=OUTDIR+OUTPLOT
PlotNiceMap(MyFile,StLats,StLons,nstats,PlotTitle,MapBoundary)
		
#    stop()

print("And, we are done!")

