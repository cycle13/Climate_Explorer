#!/usr/local/sci/bin/python

#***************************************
# 7th May 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# 7 May 2015 KMW - v1
# Reads in monthly mean anomalies for candidate (q and RH)
# Reads in the global mean time series for the candidate (q and RH)
# Reads in the global mean time series (not masked to HadISDH) e.g. SOI, PDO, CRUTEM4, HadSST3, HadCRUT4
# Correlates the global mean comparison time series with each candidate gridbox
# Works out the mean trend for each latitude (simple average)  
# Plots the correlation for each gridbox on a latitude/corr scatter with overlayed latitudinal average
# Also adds a vertical latitude by trend figure - scatter of gridbox trends for each latitude and latitude average
# Fills 0-1 with dark grey to represent actual land fraction in terms of gridboxs at that latitude band
# Overlays light grey to represent fraction of land gridboxes actually observed
# Places a time series underneath the map and latitude scatter for the candidate global mean and comparison time series

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotCorrMap_MAY2015.py
#
# REQUIRES
# RandomsRanges.py
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys, os
import scipy.stats
import struct
from mpl_toolkits.basemap import Basemap
from matplotlib.dates import date2num,num2date
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import datetime as dt
from scipy.io import netcdf
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

CompDiff='False' # 'True' means look at land-ocean T
DoCorr='False' # If 'True' then do correlation else do regression (need some spread on that - could bound boxes if slope > 2sigma of spread in residuals)

incover='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/new_coverpercentjul08'

nowmon='MAY'
nowyear='2015'

# Set up initial run choices for q
candidateGB='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
candidateTS='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014'

#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landq.2.0.1.2014p_CRUTEM4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landq.2.0.1.2014p_HadCRUT4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landq.2.0.1.2014p_HadSST3110_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landq.2.0.1.2014p_PDO_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landq.2.0.1.2014p_SOI_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landq.2.0.1.2014p_CRUTEMHadSST_'+nowmon+nowyear

#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landq.2.0.1.2014p_CRUTEM4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landq.2.0.1.2014p_HadCRUT4300_'+nowmon+nowyear
OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landq.2.0.1.2014p_HadSST3110_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landq.2.0.1.2014p_PDO_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landq.2.0.1.2014p_SOI_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landq.2.0.1.2014p_CRUTEMHadSST_'+nowmon+nowyear

#Namey=['HadISDH','land surface q','CRUTEM','land surface air T']
#Namey=['HadISDH','land surface q','HadCRUT','global surface T']
Namey=['HadISDH','land surface q','HadSST','ocean surface T']
#Namey=['HadISDH (LOWESS detrended)','land surface q','PDO (lag 7)', 'Pacific Decadal Oscillation']
#Namey=['HadISDH (LOWESS detrended)','land surface q','SOI (reversed, lag 5)','El Nino Southern Oscillation']
#Namey=['HadISDH','land surface q','CRUTEM-HadSST','land-ocean difference']

CandGridName='q_anoms' # 'q_anoms','rh_anoms'		
CandTSName='glob_q_anoms' # 'glob_q_anoms','glob_RH_anoms'		
CompTSName='glob_T_anoms' # 'q_anoms','rh_anoms'		
UnitCand='Anomalies (g kg$^{-1}$)'
#MapUnit='Correlation'
MapUnit='g kg$^{-1}$ $^{o}$C$^{-1}$'
#MapUnit='g kg$^{-1}$ SU$^{-1}$'

# set up initial run choices for RH
#candidateGB='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#candidateTS='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014'

#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landRH.2.0.1.2014p_CRUTEM4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landRH.2.0.1.2014p_HadCRUT4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landRH.2.0.1.2014p_HadSST3110_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landRH.2.0.1.2014p_PDO_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landRH.2.0.1.2014p_SOI_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CorrMapTS_HadISDH.landRH.2.0.1.2014p_CRUTEMHadSST_'+nowmon+nowyear

#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landRH.2.0.1.2014p_CRUTEM4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landRH.2.0.1.2014p_HadCRUT4300_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landRH.2.0.1.2014p_HadSST3110_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landRH.2.0.1.2014p_PDO_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landRH.2.0.1.2014p_SOI_'+nowmon+nowyear
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/RegressMap_HadISDH.landRH.2.0.1.2014p_CRUTEMHadSST_'+nowmon+nowyear

#Namey=['HadISDH','land surface RH','CRUTEM','land surface air T']
#Namey=['HadISDH','land surface RH','HadCRUT','global surface T']
#Namey=['HadISDH','land surface RH','HadSST','ocean surface T']
#Namey=['HadISDH (LOWESS detrended)','land surface RH','PDO (lag 7)', 'Pacific Decadal Oscillation']
#Namey=['HadISDH (LOWESS detrended)','land surface RH','SOI (reversed, lag 5)','El Nino Southern Oscillation']
#Namey=['HadISDH','land surface RH','CRUTEM-HadSST','land-ocean difference']

#CandGridName='rh_anoms' # 'q_anoms','rh_anoms'		
#CandTSName='glob_RH_anoms' # 'glob_q_anoms','glob_RH_anoms'		
#CompTSName='glob_T_anoms' # 'q_anoms','rh_anoms'		
#UnitCand='Anomalies (%rh)'  
#MapUnit='Correlation'
#MapUnit='%rh $^{o}$C$^{-1}$'
#MapUnit='%rh SU$^{-1}$'

#--------------------------
# Set up variables
#comparisonTS='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/CRUTEM.4.3.0.0.anomalies_areaTS_18502014'
#comparisonTS='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadCRUT.4.3.0.0.median_areaTS_18502014'
comparisonTS='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadSST.3.1.1.0.median_areaTS_18502014'
#comparisonTS='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/PDO_JISAO_19002014'
#comparisonTS='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/SOI_18762014_BOM'
comparisonTS2='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadSST.3.1.1.0.median_areaTS_18502014'

nlats=36		#set once file read in
nlons=72
Styr=1973
Edyr=2014
nyrs=Styr-(Edyr+1)
nmons=nyrs*12

CompStyr=1850 #1850, 1850, 1850, 1900, 1876
CompEdyr=2014
Compnyrs=CompStyr-(CompEdyr+1)
Compnmons=Compnyrs*12

ClimSt=1976	# 1976 or 0 if PDO or SOI
ClimEd=2005     # 2005 or 0 if PDO or SOI	

UnitComp='Anomalies ($^{o}$C)'  
#UnitComp='Difference ($^{o}$C)'  
#UnitComp='Standard Units'  

# Arrays but not setting these up yet as the shape will then be set
# Could set up to required shape?		
#CandGrid=[] # array for candidate grids
#CorrGrid=[] # array for correlation of candidate grids with comparison time series
#CandTS=[]   # array for candidate time series
#CompTS=[]   # array for comparison time series
#LatList=[]
#LonList=[]

mdi=-1e30

Letty=['a)','b)','c)']

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
    #stop()
    #TheData=var
    #    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # ReadNetCDFGrid

#************************************************************************
# READNETCDFTS
def ReadNetCDFTS(FileName,VarName):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data '''

    ncf=netcdf.netcdf_file(FileName,'r')
    var=ncf.variables[VarName]
    #TheData=var
    TheData=np.array(var.data)
#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData # ReadNetCDFTS

#************************************************************************
# READASCIITS
def ReadASCIITS(FileName,Columee,Typee,Delimee, Skippee):
    ''' Use numpy genfromtxt reading to read in all rows from a complex array '''
    ''' Typee: Need to specify format as it is complex'''
    ''' Delimee: Need to specify widths to read '''
    ''' Skippee: the number of rows to skip '''
    ''' outputs an array of tuples that in turn need to be subscripted by their names defaults f0...f8 '''
    if (len(Delimee) > 1):
        return np.genfromtxt(FileName,usecols=Columee,dtype=Typee,delimiter=Delimee,skip_header=Skippee) # ReadASCIITS
    else:
        return np.genfromtxt(FileName,usecols=Columee,dtype=Typee,skip_header=Skippee) # ReadASCIITS
#************************************************************************
# ExtractTimes
def ExtractTimes(TheData,OldSt,OldEd,NewSt,NewEd,ClimSt,ClimEd,TheMDI):
    ''' Tests for shape of array '''
    ''' If array is 2d (timeseries) - just extract new times '''
    ''' If array is 3d (grids) - extract new times for each lon/lat '''
    
    pointst=(NewSt-OldSt)*12
    pointed=((OldEd+1)-OldSt)*12
    sizee=np.shape(TheData)
    #stop()
    if (len(sizee) < 2):
        # Its a time series
	TheData=TheData[pointst:pointed+1]
    else:
        # Its a set of grids
	TheData=TheData[pointst:pointed+1,:,:] # times, lats (rows), columns
 
    if (ClimSt != 0) & (ClimEd != 0):
        # renormalise (need to sort this out for grids?
	#stop()
	TheData=np.reshape(TheData,(len(TheData)/12,12))
	for mm in range(12):
	    subarr=TheData[:,mm]
	    subclim=subarr[ClimSt-NewSt:(ClimEd+1)-NewSt]
	    climgots=np.where(subclim > TheMDI)[0]
	    gots=np.where(subarr > TheMDI)[0]
	    if (len(climgots) > 15):
	        subarr[gots]=subarr[gots]-np.mean(subclim[climgots])
	    TheData[:,mm]=subarr
	TheData=np.reshape(TheData,np.size(TheData))
	    	
    return TheData
#************************************************************************
# DeTrend
def DeTrend(TheCandData,TheMDI):
    '''' Fit a lowess to the data and remove the curve '''
    ''' default smoothing looks ok on test '''
    ''' Works on grids as well as vectors '''
    
    sizee=np.shape(TheCandData)
    if (len(sizee) < 2):
        gots=np.where(TheCandData > TheMDI)[0]
        los=lowess(TheCandData[gots],range(len(TheCandData[gots])))[:,1]
        TheCandData[gots]=TheCandData[gots]-los
    else:
        for ltt in range(len(TheCandData[0,:,0])):
	    for lnn in range(len(TheCandData[0,0,:])):
                gots=np.where(TheCandData[:,ltt,lnn] > TheMDI)[0]
		if (len(gots) > 12):
                    los=lowess(TheCandData[gots,ltt,lnn],range(len(TheCandData[gots,ltt,lnn])))[:,1]
                    TheCandData[gots,ltt,lnn]=TheCandData[gots,ltt,lnn]-los	 
		else:
		    TheCandData[:,ltt,lnn]=TheMDI           
    
    return TheCandData # DETREND
#***********************************************************************
# CorrWindow
def CorrWindow(TheCandData,TheCompTS,TheMDI):
    ''' For each gridbox that has at least 15 years of data present: '''
    ''' Temporally match with candidate time series'''
    ''' Perform correlation and save value to grid box '''
    ''' Save correlation value '''
    
    TheData=np.empty_like(TheCandData[0,:,:]) # check dimensions
    TheData.fill(TheMDI)
    #stop()
    for ltt in range(len(TheData[:,0])):
        for lnn in range(len(TheData[0,:])):
	    gots=np.where((TheCandData[:,ltt,lnn] > TheMDI) & (TheCompTS > TheMDI))[0]
	    if (len(gots) > 180): 	# 15 years * 12 months
                TheData[ltt,lnn]=np.corrcoef(TheCandData[gots,ltt,lnn],TheCompTS[gots])[0,1]
 
    return TheData
#************************************************************************
# RegressWindow
def RegressWindow(TheCandData,TheCompTS,TheMDI):
    ''' For each gridbox that has at least 15 years of data present: '''
    ''' Temporally match with candidate time series'''
    ''' Perform regression and save value to grid box '''
    ''' Save regression value and 2sigma of residuals (standard error) value '''

    TheData=np.empty_like(TheCandData[0,:,:]) # check dimensions
    TheData.fill(TheMDI)
    #stop()
    for ltt in range(len(TheData[:,0])):
        for lnn in range(len(TheData[0,:])):
	    gots=np.where((TheCandData[:,ltt,lnn] > TheMDI) & (TheCompTS > TheMDI))[0]
	    if (len(gots) > 180): 	# 15 years * 12 months
#                TheData[ltt,lnn],intercept,r_value,p_value,std_err=stats.linregress(TheCandData[gots,ltt,lnn],TheCompTS[gots])
                TheData[ltt,lnn],intercept,r_value,p_value,std_err=stats.linregress(TheCompTS[gots],TheCandData[gots,ltt,lnn])
		print('Slope: '+'{0:5.2f}'.format(TheData[ltt,lnn]),' Intercept: '+'{0:5.2f}'.format(intercept),' R value: '+'{0:5.2f}'.format(r_value),' P value: '+'{0:5.2f}'.format(p_value),' Std Err: '+'{0:5.2f}'.format(std_err))

 
    return TheData
#************************************************************************
# PlotCorrMapTS
def PlotCorrMapTS(TheFile,LandCover,TheLatList,TheLonList,TheCorrData,TheCandTS,TheCompTS,
                    TheUniteeCand,TheUniteeComp,TheLetter,TheNamee,TheMDI,YrStart,CorrNotRegress,TheMapUnit):
    ''' Create a masked array of the correlations grids '''
    ''' Plot corrs on map '''
    ''' Add vertical latitude/corr scatter with average corr overlaid '''
    ''' Add time series of CandTS and CompTS with overall correlation annotated '''
    ''' Save as eps and png '''

    # Create the masked array of corrss
    MSKTheCorrData=ma.masked_where(TheCorrData == TheMDI,TheCorrData)
    MSKTheCandTS=ma.masked_where(TheCandTS == TheMDI,TheCandTS)
    MSKTheCompTS=ma.masked_where(TheCompTS == TheMDI,TheCompTS)
    
    # Get Latitude Average Trend
    LatCorr=np.zeros(len(TheLatList))
    LatCorr.fill(TheMDI)
    # For each latitude find average correlation (not sure this is strictly sound stats - better to correlate the lat average?'
    for ltt in range(len(TheLatList)):
        gots=np.where(TheCorrData[ltt,:] != TheMDI)[0]
	if (len(gots) > 0):
	    LatCorr[ltt]=np.mean(TheCorrData[ltt,np.where(TheCorrData[ltt,:] != TheMDI)])

    # make 2d arrays of lats and lons
    # nudge -2.5 degrees to make them south/west gridbox corners, not centres
    # add extra row/column to bound the data
    ArrLons,ArrLats=np.meshgrid(TheLonList,TheLatList)
    LngArrLons,LngArrLats=np.meshgrid(np.append(TheLonList-2.5,180.),np.append(TheLatList-2.5,90.))
    
    # set up plot
    plt.clf()
    fig=plt.figure(figsize=(10,8))
    #plt.figure(1,figsize=(10,3))
    plt1=plt.axes([0.01,0.3,0.64,0.65]) # left, bottom, width, height
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))

    # make up a blue to red (reverse) colour map
    #cmaps=[plt.cm.Blues_r,plt.cm.Reds]
    #bluelist=[plt.cm.Blues_r(i) for i in range(plt.cm.Blues_r.N)]
    #redlist=[plt.cm.Reds(i) for i in range(plt.cm.Reds.N)]
#    cmap=plt.get_cmap('BrBG')
#    cmap=plt.get_cmap('coolwarm')
#    cmap=plt.get_cmap('bwr')
#    cmap=plt.get_cmap('RdGy')
    cmap=plt.get_cmap('RdBu')
        
    cmaplist=[cmap(i) for i in range(cmap.N)]
    cmaplist.reverse()
    for loo in range((cmap.N/2)-30,(cmap.N/2)+30):
        cmaplist.remove(cmaplist[(cmap.N/2)-30]) # remove the very pale colours in the middle
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    if (CorrNotRegress == 'True'):
        vmin=-0.8
        vmax=0.8
        nsteps=17
    else:
        if np.max(MSKTheCorrData) > 1.6:
	    vmin=-3
            vmax=3
            nsteps=21
	else:
	    vmin=-1.2
            vmax=1.2
            nsteps=13
	  
    
    bounds=np.linspace(vmin,vmax,nsteps)
    strbounds=["%4.1f" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCorrData,cmap=cmap,norm=norm,latlon='TRUE')

#    grids=m.pcolor(LngArrLons,LngArrLats,MSKSigTrends,cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5,latlon='TRUE')

    cbax=fig.add_axes([0.03,0.33,0.59,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=9) 
    plt.figtext(0.31,0.28,TheMapUnit,size=12,ha='center')
    plt.figtext(0.05,0.9,TheLetter[0],size=16)

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.5,0.95,TheNamee[0]+' vs '+TheNamee[2],size=16,ha='center')
    
    # Attempt to plot the latitude/ratio scatter
    ax2=plt.axes([0.73,0.393,0.25,0.46]) # map only
    ax2.set_ylim(-90,90)
    #ax2.set_ylabel('Latitude')
    ax2.set_yticks(list([-60,-30,0,30,60]))
    ax2.set_yticklabels(list(['-60$^{o}$N','-30$^{o}$S','Equator','30$^{o}$N','60$^{o}$N'])) #,rotation=90)
    ax2.set_xlabel(TheMapUnit)
    minx=np.min(TheCorrData[np.where(TheCorrData != TheMDI)])-0.1
    maxx=np.max(TheCorrData[np.where(TheCorrData != TheMDI)])+0.1
    ax2.set_xlim(minx,maxx)
    ax2.tick_params(axis='both',labelsize=10)
    # background fill the scatter plot with % land present (light grey) and % land with data (dark gre        
    for ltt in range(len(TheLatList)):
        # First merge LandCover with Data coverage so that at least some small islands are counted initially
	LandCover[ltt,np.where(TheCorrData[ltt,:] > TheMDI)[0]]=100.
        landcount=len(np.where(LandCover[ltt,:] > 0)[0])
	if (landcount > 0):
	    pctland=float(landcount)/72.
	else:
	    pctland=0
	landpresent=len(np.where(TheCorrData[ltt,:] > TheMDI)[0])
#	if (landpresent > landcount):	# small islands not picked up in landcover file so have to add in those that have been observed - not ideal
#	    landcount=landpresent
#	    pctland=float(landcount)/72.
	if (pctland > 0) & (landpresent > 0):
	    pctdata=landpresent/float(landcount)
        else:
	    pctdata=0
#	print(TheLatList[ltt],pctland,pctdata)
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
    plt.scatter(MSKTheCorrData,ArrLats,c=MSKTheCorrData,marker='o',cmap=cmap,norm=norm,edgecolor='0.0',linewidth=0.0)
    plt.plot(LatCorr[np.where(LatCorr != TheMDI)],TheLatList[np.where(LatCorr != TheMDI)],color='black',linewidth=2)
    plt.figtext(0.75,0.9,TheLetter[1],size=16)

    # Attempt to plot the timeseries at the bottom
    ax3=plt.axes([0.1,0.05,0.8,0.20]) # map only
    ax3.set_ylabel(TheUniteeCand)
    
    # set up x axes
    TheMonths=[]
    yr=YrStart
    mon=1
    for m in range(len(TheCandTS)):
        TheMonths.append(dt.date(yr,mon,1))
	mon=mon+1
	if mon == 13:
	    mon=1
	    yr=yr+1   
    TheMonths=np.array(TheMonths)
    ax3.set_xlim([TheMonths[0],TheMonths[len(TheCandTS)-1]])
    ax3.set_xlabel('Years')
    ax3.tick_params(axis='both',labelsize=10)
   
    ax3.plot(TheMonths[0:len(TheMonths)],MSKTheCandTS,c='DodgerBlue',linewidth=2)

    ax4=ax3.twinx()

    # set up the y axis - make symetrical around zero and see if you need a seperate RHS axis?
    if (TheUniteeComp == 'Standard Units'):
        miny=np.floor((np.min(TheCandTS[np.where(TheCandTS != TheMDI)])-0.1)/0.5)*0.5
        maxy=np.ceil((np.max(TheCandTS[np.where(TheCandTS != TheMDI)])+0.1)/0.5)*0.5
	if (abs(miny) > maxy):
	    maxy=abs(miny)
	else:
	    miny=-(maxy)
        ax3.set_ylim(miny,maxy)
        miny=np.floor((np.min(TheCompTS[np.where(TheCompTS != TheMDI)])-1)/5)*5
        maxy=np.ceil((np.max(TheCompTS[np.where(TheCompTS != TheMDI)])+1)/5)*5
	if (abs(miny) > maxy):
	    maxy=abs(miny)
	else:
	    miny=-(maxy)
        ax4.set_ylim(miny,maxy)
    else:
        miny=np.floor((np.min(np.array((TheCandTS[np.where(TheCandTS != TheMDI)],TheCompTS[np.where(TheCompTS != TheMDI)])))-0.1)/0.5)*0.5
        maxy=np.ceil((np.max(np.array((TheCandTS[np.where(TheCandTS != TheMDI)],TheCompTS[np.where(TheCompTS != TheMDI)])))+0.1)/0.5)*0.5
	if (abs(miny) > maxy):
	    maxy=abs(miny)
	else:
	    miny=-(maxy)
        ax3.set_ylim(miny,maxy)


    # You can do this:
    #ax4.set_yticks(range(10))
    #ax4.set_yticklabels(range(10),fontsize=16)
    # Or this:
    ax4.set_ylabel(TheUniteeComp)
    ax4.set_ylim(miny,maxy)
    ax4.tick_params(axis='both',which='major',labelsize=10)

    ax4.plot(TheMonths[0:len(TheMonths)],MSKTheCompTS,c='Firebrick',linewidth=2,alpha=0.5)
    ax4.plot(TheMonths[0:len(TheMonths)],np.zeros(len(TheMonths)),c='Black',linewidth=1)
    
    ax3.annotate(TheNamee[0],xy=(0.03,0.9),xycoords='axes fraction',size=12,ha='left',color='DodgerBlue')
    ax3.annotate(TheNamee[2],xy=(0.03,0.8),xycoords='axes fraction',size=12,ha='left',color='Firebrick')
    gots=np.where((TheCandTS > TheMDI) & (TheCompTS > TheMDI))[0]
    if (CorrNotRegress == 'True'):
        CorrVal=np.corrcoef(TheCandTS[gots],TheCompTS[gots])[0,1]
        linstr='r = '+'{0:5.2f}'.format(CorrVal)
    else:
#        CorrVal,intercept,r_value,p_value,std_err=stats.linregress(TheCandTS[gots],TheCompTS[gots])
        CorrVal,intercept,r_value,p_value,std_err=stats.linregress(TheCompTS[gots],TheCandTS[gots])
	print('Slope: '+'{0:5.2f}'.format(CorrVal))
	print('Intercept: '+'{0:5.2f}'.format(intercept))
	print('R value: '+'{0:5.2f}'.format(r_value))
	print('P value: '+'{0:5.2f}'.format(p_value))
	print('Std Err: '+'{0:5.2f}'.format(std_err))
        linstr='slope = '+'{0:5.2f}'.format(CorrVal)+' +/- '+'{0:5.2f}'.format(std_err)+' '+TheMapUnit+' ('+'{0:4.2f}'.format(p_value)+')'
   
    
    ax3.annotate(linstr,xy=(0.03,0.7),xycoords='axes fraction',size=12,ha='left',color='Black')
    plt.figtext(0.01,0.20,TheLetter[2],size=16)

    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotTrendRatMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in grids

# read in the percentage land cover dataset
ncf=netcdf.netcdf_file(incover+'.nc','r')
var=ncf.variables['pct_land']
PctLand=np.array(var.data)
#PctLand=np.transpose(PctLand)
PctLand=np.flipud(PctLand)
ncf.close()

# read in the main dataset (q, RH etc)
MyFile=candidateGB+'.nc'
CandGrid,LatList,LonList=ReadNetCDFGrid(MyFile,CandGridName)
print('Read in Cand Grid')
#stop()

# read in the global average time series of the main dataset (this is only essential for this code which plots it at the bottom)
MyFile=candidateTS+'.nc'
CandTS=ReadNetCDFTS(MyFile,CandTSName)
print('Read in Cand TS')

# Test the name of the comparisonTS (time series to compare e.g., ENSO SOI) to see if it is a netCDF file - I'm looking at the 61st character
if (comparisonTS[60] == 'C') | (comparisonTS[60] == 'H'):
    MyFile=comparisonTS+'.nc'
    CompTS=ReadNetCDFTS(MyFile,CompTSName)
    if  (CompDiff):
        MyFile=comparisonTS2+'.nc'
        CompTS2=ReadNetCDFTS(MyFile,CompTSName)
	CompTS=CompTS-CompTS2
    print('Read in Comp TS')
else: # if its not a netCDF file then it is a text file so read in the time series from the text file.
    MyFile=comparisonTS+'.txt'
    MyTypes="float" # if all elements read in are same type then output is a normal array, if not it is a structured array with a tuple for each row ANNOYING
    MyColumns=(1,2,3,4,5,6,7,8,9,10,11,12) # ignore the year (first column)
    MyDelimiters=[]
    #MyTypes=np.append("int",12*("float",))
    if comparisonTS[48] == 'P':
        MySkip=40
    else:
        MySkip=9
    RawData=ReadASCIITS(MyFile,MyColumns,MyTypes,MyDelimiters,MySkip)
    # need to rearrange this to make a month time series
    CompTS=np.reshape(RawData,np.size(RawData)) # puts each year of months consecutively in one vector

# do some lagging to get the best fit if comparisonTS is the SOI or PDO - already tested this and decided on lags 5 for SOI and 7 for PDO
# this rolls the data forward so CandTS Jan 1973 Compares with CompTS Dec 1972
# The comparisonTS is a long time series so we'll need to chop it to match the time period of the candidate data
# However, rolling the data needs to be done before chopping up so that it doesn't move Dec 2014 to Jan 1973
# for the ENSO SOI:
if comparisonTS[48] == 'S':
    CompTS=np.roll(-CompTS,5)
# for the PDO
if comparisonTS[48] == 'P':
    CompTS=np.roll(CompTS,7)

# cut grids and ts down to desired time points and renormalise if necessary
#CandGrid=ExtractTimes(CandGrid,CompStyr,CompEdyr,Styr,Edyr,ClimSt,ClimEd,mdi)    
#CandTS=ExtractTimes(CandTS,CompStyr,CompEdyr,Styr,Edyr,ClimSt,ClimEd,mdi) 
#print('Extracted CandTS')   
CompTS=ExtractTimes(CompTS,CompStyr,CompEdyr,Styr,Edyr,ClimSt,ClimEd,mdi)    
print('Extracted CompTS')   

# detrend the Candidate data if the comparisonTS is SOI and PDO
if (comparisonTS[48] == 'S') | (comparisonTS[48] == 'P'):
    # first check what the correlation is without detrending.
    gots=np.where((CandTS > mdi) & (CompTS > mdi))[0]
    print('Non-detrended correlation: ',np.corrcoef(CandTS[gots],CompTS[gots])[0,1])
    CandTS=DeTrend(CandTS,mdi)
    CandGrid=DeTrend(CandGrid,mdi)

# get the correlations for each gridbox if DoCorr is true
if (DoCorr =='True'):
    CorrGrid=CorrWindow(CandGrid,CompTS,mdi)
    print('Got Corrs')
# get the regressions for each gridbox if DoCorr is true
else:
    CorrGrid=RegressWindow(CandGrid,CompTS,mdi)
    print('Got Regresses')

# pass to plotter
PlotCorrMapTS(OUTPLOT,PctLand,LatList,LonList,CorrGrid,CandTS,CompTS,
                        UnitCand,UnitComp,Letty,Namey,mdi,Styr,DoCorr,MapUnit)
		
#    stop()

print("And, we are done!")

