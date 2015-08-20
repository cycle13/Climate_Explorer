#!/usr/local/sci/bin/python

#***************************************
# 22nd June 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# 22 June 2015 KMW - v1
# This reads in a netCDF file of a gridded dataset
# The month of choice is plotted as an anomaly map
# Data can be renormalised to the desired climatology period

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotAnyMap_JUN2015.py
#
# REQUIRES
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

# Set up initial run choices
ChooseMon=11	# 0 to 11 for Jan to Dec
ChooseYr=2014
ChooseClimSt=1976
ChooseClimEd=2005
DoReClim=False # True if need to renormalise, False if not

MonArr=(['January','February','March','April','May','June','July','August','September','October','November','December'])

candidate='HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#candidate='HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#candidate='HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_cf'
#candidate='HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#candidate='BERKELEY_T_5by519762005clim_anoms_19732014'
#candidate='GISS_T_5by519762005clim_anoms_19732014'
#candidate='CRUTEM.4.3.0.0.anomalies'
#candidate='GHCNM_18802014'

NameOfVar='q_anoms'
#NameOfVar='rh_anoms'
#NameOfVar='q_anoms'
#NameOfVar='q_anoms'
#NameOfVar='q_anoms'

nowmon='JUN'
nowyear='2015'

# Set up directories and files
INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
#INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'

OUTPLOT='Map_HadISDH.landq.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_HadISDH.landRH.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_HadISDH.landT.2.0.1.2014p_RAW_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_HadISDH.landT.2.0.1.2014p_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_BERKELEY_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_GISS_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_CRUTEM4.3.0.0_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)
#OUTPLOT='Map_GHCNM3_clim'+str(ChooseClimSt)+str(ChooseClimEd)+'_'+MonArr[ChooseMon]+str(ChooseYr)

# Set up variables
nlats=36		#set once file read in
nlons=72		#set once file read in
CandFields=[]
CandData=[]
LatList=[]
LonList=[]
StYr=1973
EdYr=2014
NYrs=(EdYr+1)-StYr
NMons=NYrs*12
TimePointer=(ChooseYr-StYr)*12+ChooseMon

Unit='g kg$^{-1}$'  #'degrees C'
#Unit='%rh'  #'degrees C'
#Unit='$^{o}$C'  #'degrees C'

Namey='HadISDH.landq.2.0.1.2014p '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='HadISDH.landRH.2.0.1.2014p '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='HadISDH.landT.2.0.1.2014p RAW '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='HadISDH.landT.2.0.1.2014p '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='BERKELEY '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='GISS '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='CRUTEM4.3.0.0 '+MonArr[ChooseMon]+' '+str(ChooseYr)
#Namey='GHCNM3 '+MonArr[ChooseMon]+' '+str(ChooseYr)

# Missing data
mdi=-1e30 # may set up as masked arrays later

# Specific Humidity
vmin=-1.6 # Specific Humidity
vmax=1.6
nsteps=17
## Relative Humidity
#vmin=-10. 
#vmax=10.
#nsteps=21


#************************************************************************
# Subroutines
#************************************************************************
# READNETCDFGRID
def ReadNetCDFGrid(FileName,VarName):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data '''

    ncf=netcdf.netcdf_file(FileName,'r')
    # ncf.variables this lists the variable names
    #var=f.variables['latitude']
    #TheLatList=var.data
    # lats currently screwy so make up
    TheLatList=np.arange(-87.5,92.5,5.)
    #var=f.variables['longitude']
    #TheLonList=var.data
    # lons currently screwy so make up
    TheLonList=np.arange(-177.5,182.5,5.)
    var=ncf.variables[VarName]
    TheData=np.array(var.data)

#    # Maybe I've done something wrong but its reading it transposed
    #print(np.shape(TheData))
    #stop
    #TheData=np.transpose(TheData)
    ncf.close()

#    tmpvar=np.array(np.transpose(var.data[:]))
    
    return TheData,TheLatList,TheLonList # ReadNetCDFGrid

#************************************************************************
# PlotAnyMap
def PlotAnyMap(TheFile,TheLatList,TheLonList,TheCandData,TheUnitee,TheNamee,vmin,vmax,nsteps):
    ''' Create a masked array of variable
        Plot on map
        Save as eps and png '''

    # Missing data
    mdi=-1e30 # may set up as masked arrays later
      
    # Create the masked array of trends
    MSKTheCandData=ma.masked_where(TheCandData == mdi,TheCandData)
    
    # make 2d arrays of lats and lons
    # nudge -2.5 degrees to make them south/west gridbox corners, not centres
    # add extra row/column to bound the data
    ArrLons,ArrLats=np.meshgrid(TheLonList,TheLatList)
    LngArrLons,LngArrLats=np.meshgrid(np.append(TheLonList-2.5,180.),np.append(TheLatList-2.5,90.))
    
    # set up plot
    plt.clf()
    fig=plt.figure(figsize=(6,5))
    plt1=plt.axes([0.01,0.01,0.98,0.98]) # left, bottom, width, height
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))

    # make up a blue to red (reverse) colour map
    cmap=plt.get_cmap('BrBG')      # HUMIDITY
#    cmap=plt.get_cmap('coolwarm') # TEMPERATURE
#    cmap=plt.get_cmap('bwr')
    
    cmaplist=[cmap(i) for i in range(cmap.N)]
    for loo in range((cmap.N/2)-20,(cmap.N/2)+20):
        cmaplist.remove(cmaplist[(cmap.N/2)-20]) # remove the very pale colours in the middle
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCandData,cmap=cmap,norm=norm,latlon='TRUE')

    cbax=fig.add_axes([0.05,0.08,0.9,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=10) 
    plt.figtext(0.5,0.01,'Anomalies ('+TheUnitee+')',size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.5,0.95,TheNamee,size=16,ha='center')
        
    
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
CandFields,LatList,LonList=ReadNetCDFGrid(MyFile,NameOfVar)

# Do you need to renormalise?
#if (DoReClim):
    #MAKE UP CODE HERE

# Extract chosen month
CandData=CandFields[TimePointer,:,:]

# pass to plotter
MyFile=OUTDIR+OUTPLOT
PlotAnyMap(MyFile,LatList,LonList,CandData,Unit,Namey,vmin,vmax,nsteps)
		
#    stop()

print("And, we are done!")

