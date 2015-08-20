#!/usr/local/sci/bin/python

#***************************************
# 17th March 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# 4 March 2015 KMW - v1
# Reads in decadal trend map from two different sources
# Works out the absolute ratio of trends of each gridbox  
# Identifies the sign of the ratio
# Grades the gridbox by the underlying candidate trend
# Plots the ratio for each gridbox, faded relative to the magnitude of the trend:
# white=0, bright=max
# Also adds a vertical latitude by ratio figure - ratio of average trend at that latitude
# Lists the 10 grid boxes with the 'worst' ratios - largest negative

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotTrendRatioMap_MAR2015.py
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

# Set up initial run choices
#candidate='HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
candidate='CRUTEM.4.3.0.0.anomalies_MPtrends_19732014'
#comparree='CRUTEM.4.3.0.0.anomalies_MPtrends_19732014'
comparree='GHCNM_18802014_MPtrends_19732014'
#covermap='new_coverpercentjul08.nc'
coverage='HadCRUT.4.3.0.0.land_fraction'
nowmon='MAR'
nowyear='2015'

# Set up directories and files
#INDIR='/DATA/'
#OUTDIR='/FIGURES/'
INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TRENDS/'
INDIRO='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TRENDS/'
#INDIRC='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
#INDIRO='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/'

#OUTPLOT='TrendRatioMap_HadISDHCRUTEM4_'+nowmon+nowyear
#OUTPLOT='TrendRatioMap_HadISDHGHCNM_'+nowmon+nowyear
OUTPLOT='TrendRatioMap_CRUTEM4GHCNM_'+nowmon+nowyear

# Set up variables
nlats=36		#set once file read in
nlons=72		#set once file read in
CandTrends=[]
CompTrends=[]
RatTrends=[]
LatList=[]
LonList=[]

Unit='$^{o}$C'  #'degrees C'
Letty=' '
#Namey='HadISDH vs CRUTEM4'
#Namey='HadISDH vs GHCNM'
Namey='CRUTEM4 vs GHCNM'

#************************************************************************
# Subroutines
#************************************************************************
# READNETCDFGRID
def ReadNetCDFGrid(FileName):
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
    var=ncf.variables['T_MPtrend']
    TheData=np.array(var.data)
#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # ReadNetCDFGrid

#************************************************************************
# PlotTrendRatMap
def PlotTrendRatMap(TheFile,TheLatList,TheLonList,TheCandData,TheCompData,
                    TheUnitee,TheLetter,TheNamee):
    ''' Create a ratio field
        Create a light/dark field (eventually just scale transparency)
        Create a vector of lats for each longitude
        Create a vector of all gridboxes loop lats, loop lons
        Plot ratios on map, scaled light/dark
        Add vertical latitude/ratio scatter with ratio of average trend overlaid
        Save as eps and png '''

    # Missing data
    mdi=-1e30 # may set up as masked arrays later
      
    # Create the ratio fields
    # Set up empty arrays
    RatTrends=np.zeros((len(TheLatList),len(TheLonList)))
    RatTrends[:,:]=mdi
    LatTrendCand=np.zeros(len(TheLatList))
    LatTrendCand[:]=mdi
    LatTrendComp=LatTrendCand.copy()
    LatRats=LatTrendCand.copy()
    
    # Find non-missing pairs of grids and get ratios of trends
    gots=np.where((TheCandData > mdi) & (TheCandData != 0.0) & (TheCompData > mdi) & (TheCompData != 0.0))
    RatTrends[gots]=TheCandData[gots]/TheCompData[gots]
#    stop

    # For each latitude find average trend for each source, then get ratio
    for ltt in range(len(TheLatList)):
        gots=np.where(TheCandData[ltt,:] != mdi)[0]
	if (len(gots) > 0):
	    LatTrendCand[ltt]=np.mean(TheCandData[ltt,np.where(TheCandData[ltt,:] != mdi)])
        gots=np.where(TheCompData[ltt,:] != mdi)[0]
	if (len(gots) > 0):
	    LatTrendComp[ltt]=np.mean(TheCompData[ltt,np.where(TheCompData[ltt,:] != mdi)])
        if (LatTrendComp[ltt] != mdi) and (LatTrendCand[ltt] != mdi):
	    LatRats[ltt]=LatTrendCand[ltt]/LatTrendComp[ltt]

    # make 2d arrays of lats and lons
    # these need to be shifted south/west so that they are bottom left corners, not gridbox centres
    # these should have one extra on the end to bound the data so that the last row/column is plotted
    ShtArrLons,ShtArrLats=np.meshgrid(TheLonList,TheLatList)
    ArrLons,ArrLats=np.meshgrid(np.append(TheLonList-2.5,180),np.append(TheLatList-2.5,90))
    
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

    # Firebrick = Ratio < -1, Cand Trend > abs(0.05) or not ARBITRARY
    # LightCoral = Ratio > -1 and < 0, Cand Trend > abs(0.05) or not
    # Silver -0.05 to 0.05
    # LightSkyBlue = Ratio > 0 and < 1
    # DodgerBlue = Ratio > 1
    
    cmap = mc.ListedColormap(['DodgerBlue','LightSkyBlue','Khaki','LightCoral','Firebrick'])

    # Set up the colourmap
    RatTrendCOL=np.zeros((len(TheLatList),len(TheLonList)))
    #RatTrendCOL=np.array(RatTrends.copy(),dtype='i2')

    # Now large positives
    RatTrendCOL[np.where(RatTrends > 1.)]=1
    # Now small positives
    RatTrendCOL[np.where((RatTrends > 0.01) & (RatTrends <= 1.))]=2
    # Now very smalls
    RatTrendCOL[np.where((TheCandData > -0.01) & (TheCandData <= 0.01) & (RatTrends != mdi))]=3
    # Now small negatives
    RatTrendCOL[np.where((RatTrends > -1.) & (RatTrends <= -0.01))]=4
    # Now large negatives
    RatTrendCOL[np.where((RatTrends < -1.) & (RatTrends != mdi))]=5

#    # Now nest in darker colours where The Candidate trends are teeny tiny values
#
#    # Now large positives
#    RatTrendCOL[np.where((RatTrendCOL == 2) & (abs(TheCandData) < 0.05))]=1
#    # Now small positives
#    RatTrendCOL[np.where((RatTrendCOL == 4) & (abs(TheCandData) < 0.05))]=3
#    # Now small negatives
#    RatTrendCOL[np.where((RatTrendCOL == 6) & (abs(TheCandData) < 0.05))]=5
#    # Now large negatives
#    RatTrendCOL[np.where((RatTrendCOL == 8) & (abs(TheCandData) < 0.05))]=7

    # so last find mdis to mask#
    MSKRatTrendCOL=ma.masked_where(RatTrends == mdi,RatTrendCOL)
    
    # count the LN, SN, SP and LPs
    #CountGoods=len(RatTrends[np.where(RatTrends > mdi)])
    CountGoods=len(np.where(RatTrends > mdi)[0])
    CountLargeNegs=len(np.where(RatTrendCOL == 5)[0])
    CountSmallNegs=len(np.where(RatTrendCOL == 4)[0])
    CountVSmalls=len(np.where(RatTrendCOL == 3)[0])
    CountVSmallsP=len(np.where((RatTrendCOL == 3) & (RatTrends > 0.))[0])
    CountVSmallsN=len(np.where((RatTrendCOL == 3) & (RatTrends <= 0.))[0])
    CountSmallPos=len(np.where(RatTrendCOL == 2)[0])
    CountLargePos=len(np.where(RatTrendCOL == 1)[0])
    
    m.pcolormesh(ArrLons,ArrLats,MSKRatTrendCOL,cmap=cmap,vmin=1,vmax=5,latlon='TRUE')
    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
#    plt.figtext(0.05,0.9,TheLetter,size=16)
    plt.figtext(0.5,0.95,TheNamee,size=16,ha='center')
    
    print('TOTALS: ',CountGoods,CountLargeNegs,CountSmallNegs,CountVSmalls,CountVSmallsP,CountVSmallsN,CountSmallPos,CountLargePos)
    pctLNs=round((float(CountLargeNegs)/CountGoods)*100.,1)
    pctSNs=round((float(CountSmallNegs)/CountGoods)*100.,1)
    pctVSs=round((float(CountVSmalls)/CountGoods)*100.,1)
    pctVSsP=round((float(CountVSmallsP)/CountGoods)*100.,1)
    pctVSsN=round((float(CountVSmallsN)/CountGoods)*100.,1)
    pctSPs=round((float(CountSmallPos)/CountGoods)*100.,1)
    pctLPs=round((float(CountLargePos)/CountGoods)*100.,1)
    print('PCTS: ',pctLNs,pctSNs,pctVSs,pctVSsP,pctVSsN,pctSPs,pctLPs)
    print('MAX/MIN ratios: ',np.min(RatTrends[np.where(RatTrends > mdi)]),np.max(RatTrends[np.where(RatTrends > mdi)]))
    
    ax1=plt.axes([0.01,0.01,0.64,0.09],frameon=False) # map only
    ax1.set_ylim(0,1)
    ax1.set_xlim(0,1)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.annotate(str(pctLNs)+"% ratio < -1",xy=(0.07,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.05,0.8,markersize=10,marker='s',color='Firebrick')
    plt.annotate(str(pctSNs)+"% ratio -1 to -0.05",xy=(0.07,0.1),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.05,0.2,markersize=10,marker='s',color='LightCoral')
    plt.annotate(str(pctVSs)+"% small trend",xy=(0.40,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.38,0.8,markersize=10,marker='s',color='Khaki')
    plt.annotate(str(pctSPs)+"% ratio 0.05 to 1",xy=(0.73,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.71,0.8,markersize=10,marker='s',color='LightSkyBlue')
    plt.annotate(str(pctLPs)+"% ratio > 1",xy=(0.73,0.1),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.71,0.2,markersize=10,marker='s',color='DodgerBlue')

    # Attempt to plot the latitude/ratio scatter
    ax2=plt.axes([0.73,0.12,0.24,0.76]) # map only
    ax2.set_ylim(-90,90)
   # ax2.set_ylabel('Latitude')
    ax2.set_xlabel('Ratio')
    ax2.set_xlim(np.min(RatTrends[np.where(RatTrends != mdi)])-1,np.max(RatTrends[np.where(MSKRatTrendCOL != 3)])+1)

    ax2.set_yticks(list([-60,-30,0,30,60]))
    ax2.set_yticklabels(list(['-60$^{o}$N','-30$^{o}$S','Equator','30$^{o}$N','60$^{o}$N'])) #,rotation=90)
    ax2.tick_params(axis='both',labelsize=10)
    #ax.set_xmargin(0.01)
    plt.plot(np.zeros(2),np.array((-90,90)),color='black')
    plt.scatter(RatTrends[np.where(MSKRatTrendCOL == 1)],ShtArrLats[np.where(MSKRatTrendCOL == 1)],marker='o',color='DodgerBlue')
    plt.scatter(RatTrends[np.where(MSKRatTrendCOL == 2)],ShtArrLats[np.where(MSKRatTrendCOL == 2)],marker='o',color='LightSkyBlue')
    plt.scatter(RatTrends[np.where(MSKRatTrendCOL == 4)],ShtArrLats[np.where(MSKRatTrendCOL == 4)],marker='o',color='LightCoral')
    plt.scatter(RatTrends[np.where(MSKRatTrendCOL == 5)],ShtArrLats[np.where(MSKRatTrendCOL == 5)],marker='o',color='Firebrick')
    plt.plot(LatRats[np.where(LatRats != mdi)],TheLatList[np.where(LatRats != mdi)],color='black')
        
    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    # Output list of all 'worst' ratios (largest negatives to smallest negatives)
    bads=np.where(MSKRatTrendCOL >= 4)
    print("Number of BAD ratios: ",len(bads[0]))
    print("Ratios GBCentre Longs GBCentre Lats")
    SubRats=RatTrends[bads]
    SubLons=ArrLons[bads]
    SubLats=ArrLats[bads]
    SortingHat=np.argsort(SubRats)
#    stop
    for loo in range(len(bads[0])):
        print("%6.2f %8.3f %8.3f" % (SubRats[SortingHat[loo]],SubLons[SortingHat[loo]],SubLats[SortingHat[loo]]))
     
#    stop
     
    return #PlotTrendRatMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in trend maps
MyFile=INDIRC+candidate+'.nc'
CandData,LatList,LonList=ReadNetCDFGrid(MyFile)
MyFile=INDIRO+comparree+'.nc'
CompData,LatList,LonList=ReadNetCDFGrid(MyFile)

# pass to plotter
MyFile=OUTDIR+OUTPLOT
PlotTrendRatMap(MyFile,LatList,LonList,CandData,CompData,
                        Unit,Letty,Namey)
		
#    stop()

print("And, we are done!")

