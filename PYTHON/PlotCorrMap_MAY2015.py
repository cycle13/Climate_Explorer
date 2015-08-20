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
# python2.7 PlotCorrMap_MAY2015.py
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

Typee='Humidity' #'Temperature', 'HUmidity' - for telling the plotter which colours/scales etc.

candidate='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#candidate='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies'
#comparree='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies'
#comparree='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/GHCNM_18802014'
comparree='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'
#comparree='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landRH.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf'

#incover='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/new_coverpercentjul08'
incover='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/HadCRUT.4.3.0.0.land_fraction'

nowmon='JUN'
nowyear='2015'

# Set up directories and files
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/'

#OUTPLOT='CorrMap_HadISDHCRUTEM4_'+nowmon+nowyear
#OUTPLOT='CorrMap_HadISDHGHCNM_'+nowmon+nowyear
#OUTPLOT='CorrMap_CRUTEM4GHCNM_'+nowmon+nowyear
OUTPLOT='CorrMap_HadISDHT_q_'+nowmon+nowyear
#OUTPLOT='CorrMap_HadISDHT_RH_'+nowmon+nowyear

# Set up variables
nlats=36		#set once file read in
nlons=72		#set once file read in
CandTrends=[]
CompTrends=[]
RatTrends=[]
LatList=[]
LonList=[]

Unit='correlation'  #'degrees C'
#Namey='HadISDH vs CRUTEM4'
#Namey='HadISDH vs GHCNM'
#Namey='CRUTEM4 vs GHCNM'
Namey='HadISDH.landT vs HadISDH.landq'
#Namey='HadISDH.landT vs HadISDH.landRH'

CandGridName='t_anoms' # 'q_anoms','rh_anoms'		
#CandGridName='temperature_anomaly' # 'q_anoms','rh_anoms'		
#CompGridName='temperature_anomaly' # 'q_anoms','rh_anoms'		
CompGridName='q_anoms' # 'q_anoms','rh_anoms'		
#CompGridName='rh_anoms' # 'q_anoms','rh_anoms'		

nlats=36		#set once file read in
nlons=72
Styr=1973
Edyr=2014
nyrs=Styr-(Edyr+1)
nmons=nyrs*12

CandStyr=1973 #1850, 1880
#CandStyr=1850 #1850, 1880
CompStyr=1973 #1850, 1880, 1973
CompEdyr=2014
Compnyrs=CompStyr-(CompEdyr+1)
Compnmons=Compnyrs*12

ClimSt=1976	# 1976 or 0 if PDO or SOI
ClimEd=2005     # 2005 or 0 if PDO or SOI	

# Arrays but not setting these up yet as the shape will then be set
# Could set up to required shape?		
#CandGrid=[] # array for candidate grids
#CorrGrid=[] # array for correlation of candidate grids with comparison time series
#CandTS=[]   # array for candidate time series
#CompTS=[]   # array for comparison time series
#LatList=[]
#LonList=[]

mdi=-1e30

Letty=['a)','b)']

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
        
	NewData=np.empty_like(TheData)
	NewData.fill(TheMDI)
    
        if (ClimSt != 0) & (ClimEd != 0):
            # renormalise (need to sort this out for grids?
	    #stop()
	    NewData=np.reshape(TheData,(len(TheData)/12,12))
	    for mm in range(12):
	        subarr=NewData[:,mm]
	        subclim=subarr[ClimSt-NewSt:(ClimEd+1)-NewSt]
	        climgots=np.where(subclim > TheMDI)[0]
	        gots=np.where(subarr > TheMDI)[0]
	        if (len(climgots) > 15):
	            subarr[gots]=subarr[gots]-np.mean(subclim[climgots])
	        NewData[:,mm]=subarr
	    NewData=np.reshape(TheData,np.size(TheData))
    else:
        # Its a set of grids
	TheData=TheData[pointst:pointed+1,:,:] # times, lats (rows), columns
 
        NewData=np.empty_like(TheData)
	NewData.fill(TheMDI)
    
        if (ClimSt != 0) & (ClimEd != 0):
            # renormalise (need to sort this out for grids?
	    #stop()
	    for ltt in range(len(TheData[0,:,0])):
	        for lnn in range(len(TheData[0,0,:])):
	            ShortData=np.reshape(TheData[:,ltt,lnn],(len(TheData[:,ltt,lnn])/12,12))
	            for mm in range(12):
	                subarr=ShortData[:,mm]
	                subclim=subarr[ClimSt-NewSt:(ClimEd+1)-NewSt]
	                climgots=np.where(subclim > TheMDI)[0]
	                gots=np.where(subarr > TheMDI)[0]
	                if (len(climgots) > 15):
	                    subarr[gots]=subarr[gots]-np.mean(subclim[climgots])
	                ShortData[:,mm]=subarr
	            NewData[:,ltt,lnn]=np.reshape(ShortData,np.size(ShortData))
	    	
    return NewData
#************************************************************************
# CorrWindow
def CorrWindow(TheCandData,TheCompData,TheMDI):
    ''' For each gridbox that has at least 15 years of data present: '''
    ''' Temporally match with candidate time series'''
    ''' Perform correlation and save value to grid box '''
    ''' Save correlation value '''
    
    TheData=np.empty_like(TheCandData[0,:,:]) # check dimensions
    TheData.fill(TheMDI)
    #stop()
    for ltt in range(len(TheData[:,0])):
        for lnn in range(len(TheData[0,:])):
	    gots=np.where((TheCandData[:,ltt,lnn] > TheMDI) & (TheCompData[:,ltt,lnn] > TheMDI))[0]
	    if (len(gots) > 180): 	# 15 years * 12 months
                TheData[ltt,lnn]=np.corrcoef(TheCandData[gots,ltt,lnn],TheCompData[gots,ltt,lnn])[0,1]
 
    return TheData
#************************************************************************
# PlotCorrMapTS
def PlotCorrMapTS(TheFile,LandCover,TheLatList,TheLonList,TheCorrData,
                    TheUnitee,TheLetter,TheNamee,TheMDI,YrStart,VarType):
    ''' Create a masked array of the correlations grids '''
    ''' Plot corrs on map '''
    ''' Add vertical latitude/corr scatter with average corr overlaid '''
    ''' Add time series of CandTS and CompTS with overall correlation annotated '''
    ''' Save as eps and png '''

    # Create the masked array of corrss
    MSKTheCorrData=ma.masked_where(TheCorrData == TheMDI,TheCorrData)
    
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
    #cmaps=[plt.cm.Blues_r,plt.cm.Reds]
    #bluelist=[plt.cm.Blues_r(i) for i in range(plt.cm.Blues_r.N)]
    #redlist=[plt.cm.Reds(i) for i in range(plt.cm.Reds.N)]
#    cmap=plt.get_cmap('BrBG')
#    cmap=plt.get_cmap('coolwarm')
#    cmap=plt.get_cmap('bwr')
#    cmap=plt.get_cmap('RdGy')
#    cmap=plt.get_cmap('RdBu')
#    cmap=plt.get_cmap('hot')

    if (VarType == 'Temperature'):
        # Colour for T T
        cmap=plt.get_cmap('YlOrRd')
    elif (VarType == 'Humidity'):    
	# Colour for T /RH
        cmap=plt.get_cmap('bwr')
        
    cmaplist=[cmap(i) for i in range(cmap.N)]

    if (VarType == 'Temperature'):
        # Only for T T
        for loo in range(30):
            cmaplist.remove(cmaplist[0])
        cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    #vmin=0.5
    #vmax=1
    #nsteps=11	  
    #bounds=np.round(1-np.logspace(0.7,-1,num=15,endpoint='True')/10.,3)
    if (VarType == 'Temperature'):
        # bounds for T T
        bounds=np.array([0,0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.99,0.995,0.999,1])
    elif (VarType == 'Humidity'):    
	# bounds for T q/RH
        bounds=np.array([-1,-0.9,-0.8,-0.7,-0.6,-0.5,0.0,0.5,0.6,0.7,0.8,0.9,1])

    #bounds=np.logspace(vmin,vmax,nsteps)
#    strbounds=["%4.1f" % i for i in bounds]
    strbounds=["%4.3g" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    grids=m.pcolor(LngArrLons,LngArrLats,MSKTheCorrData,cmap=cmap,norm=norm,latlon='TRUE')


    # extract the poorly correlating stations - reshape to a single array, index by lowest to highest, create a new zero array where 1:10 only have values
    # replot with boundaries.
    sortedlowhigh=np.argsort(np.reshape(MSKTheCorrData,np.size(MSKTheCorrData)))	# Does this work on masks?
    #stop
    LowCorrs=np.zeros(np.size(TheCorrData))
    TheCorrData=np.reshape(TheCorrData,np.size(TheCorrData))
    LowCorrs[sortedlowhigh[0:10]]=TheCorrData[sortedlowhigh[0:10]]
    LowCorrs=np.reshape(LowCorrs,(len(TheLatList),len(TheLonList)))
    TheCorrData=np.reshape(TheCorrData,(len(TheLatList),len(TheLonList)))
    MSKLowCorrs=ma.masked_where(LowCorrs == 0,LowCorrs)
    
    grids=m.pcolor(LngArrLons,LngArrLats,MSKLowCorrs,cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5,latlon='TRUE')

    # Output list of all 'worst' correlations
    bads=np.where(MSKLowCorrs > 0)
    print("Number of Low Corrs: ",len(bads[0]))
    print("Corrs GBCentre Longs GBCentre Lats")
    SubCorrs=LowCorrs[bads]
    SubLons=ArrLons[bads]
    SubLats=ArrLats[bads]
    SortingHat=np.argsort(SubCorrs)
#    stop
    for loo in range(len(bads[0])):
        print("%6.2f %8.3f %8.3f" % (SubCorrs[SortingHat[loo]],SubLons[SortingHat[loo]],SubLats[SortingHat[loo]]))
     

    cbax=fig.add_axes([0.03,0.08,0.59,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=9) 
    cb.ax.set_xticklabels(strbounds)
    plt.figtext(0.31,0.01,TheUnitee,size=12,ha='center')
    plt.figtext(0.05,0.9,TheLetter[0],size=16)

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.5,0.95,TheNamee,size=16,ha='center')
    
    # Attempt to plot the latitude/ratio scatter
    ax2=plt.axes([0.73,0.12,0.24,0.76]) # map only
    ax2.set_ylim(-90,90)
    #ax2.set_ylabel('Latitude')
    ax2.set_yticks(list([-60,-30,0,30,60]))
    ax2.set_yticklabels(list(['-60$^{o}$N','-30$^{o}$S','Equator','30$^{o}$N','60$^{o}$N'])) #,rotation=90)
    ax2.set_xlabel(TheUnitee)
#    minx=np.min(TheCorrData[np.where(TheCorrData != TheMDI)])-0.1
#    maxx=np.max(TheCorrData[np.where(TheCorrData != TheMDI)])+0.1
    if (VarType == 'Temperature'):
        # min for T T
        minx=0.
    elif (VarType == 'Humidity'):    
	# min for T q/RH
        minx=-1.
    maxx=1.
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
    plt.figtext(0.7,0.9,TheLetter[1],size=16)
    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotTrendRatMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in grids

ncf=netcdf.netcdf_file(incover+'.nc','r')
#var=ncf.variables['pct_land']
var=ncf.variables['land_area_fraction']
PctLand=np.array(var.data)
#PctLand=np.transpose(PctLand)
PctLand=np.flipud(PctLand)
PctLand=PctLand[0,:,:]
ncf.close()

#stop

MyFile=candidate+'.nc'
CandGrid,LatList,LonList=ReadNetCDFGrid(MyFile,CandGridName)
print('Read in Cand Grid')
MyFile=comparree+'.nc'
CompGrid,CLatList,CLonList=ReadNetCDFGrid(MyFile,CompGridName)
print('Read in Comp Grid')
#stop()

# cut grids and ts down to desired time points and renormalise if necessary
CandGrid=ExtractTimes(CandGrid,CandStyr,CompEdyr,Styr,Edyr,ClimSt,ClimEd,mdi)    
CompGrid=ExtractTimes(CompGrid,CompStyr,CompEdyr,Styr,Edyr,ClimSt,ClimEd,mdi)    
print('Extracted CompGRID')   


# get the correlations for each gridbox if DoCorr is true
CorrGrid=CorrWindow(CandGrid,CompGrid,mdi)
print('Got Corrs')

print(Namey)
print('Total Data Present: ',len(np.where(CorrGrid > mdi)[0]),len(np.where(CorrGrid > mdi)[0])/float(len(np.where(PctLand > 0.)[0])))
print('Total Greater than 0.9: ',len(np.where(CorrGrid > 0.9)[0]),len(np.where(CorrGrid > 0.9)[0])/float(len(np.where(CorrGrid > mdi)[0])))
print('Total Greater than 0.8: ',len(np.where(CorrGrid > 0.8)[0]),len(np.where(CorrGrid > 0.8)[0])/float(len(np.where(CorrGrid > mdi)[0])))
print('Total Greater than 0.7: ',len(np.where(CorrGrid > 0.7)[0]),len(np.where(CorrGrid > 0.7)[0])/float(len(np.where(CorrGrid > mdi)[0])))
print('Total Greater than 0.6: ',len(np.where(CorrGrid > 0.6)[0]),len(np.where(CorrGrid > 0.6)[0])/float(len(np.where(CorrGrid > mdi)[0])))
#stop

# pass to plotter
PlotCorrMapTS(OUTDIR+OUTPLOT,PctLand,LatList,LonList,CorrGrid,
                        Unit,Letty,Namey,mdi,Styr,Typee)
		
#    stop()

print("And, we are done!")

