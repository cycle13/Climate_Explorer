#!/usr/local/sci/bin/python

#***************************************
# 23 April 2015 KMW - v1
# Plots a time series for any gridbox
# Option to overplot the same gridbox from another source too 
# 
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotGBTimeSeries_APR2015.py
#
# REQUIRES
# RandomsRanges.py
# LinearTrends.py
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy.ma as ma
import numpy as np
import sys, os
import scipy.stats
import struct
import os.path
import math
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
#from netCDF4 import Dataset
from scipy.io import netcdf
from scipy.stats.stats import pearsonr   
from RandomsRanges import LetterRange
from LinearTrends import MedianPairwise

# Set up initial run choices
timetype='monthly'	#'monthly', 'annual'
fitlin=True # False - fit a linear trend?
addcor=True # False - correlate the data?

# What are the units?
unitees='$^{o}$C'  #'deg C'

# Number of sources to plot
nsources=3 	# this can be one to three

# Source Names for Title
#Namey='HadISDH.2.0.1.2014p vs CRUTEM4.3.0.0'
#Namey='HadISDH.2.0.1.2014p vs GHCNM3'
Namey='CRUTEM4.3.0.0 vs GHCNM3'
SourceNames=['HadISDH.2.0.1.2014p','CRUTEM4.3.0.0','GHCNM3','HadISDH.2.0.1.2014p']
##SourceNames=['HadISDH.2.0.1.2014p','GHCNM3','CRUTEM4.3.0.0']
##SourceNames=['CRUTEM4.3.0.0','GHCNM3','HadISDH.2.0.1.2014p']
ShortSourceNames=['HadISDH','CRUTEM','GHCNM','HadISDH adj-raw']
##ShortSourceNames=['HadISDH','GHCNM','CRUTEM','HadISDH adj-raw']
##ShortSourceNames=['CRUTEM','GHCNM','HadISDH','HadISDH adj-raw']
ColorArr=['DimGrey','Red','DodgerBlue','Gold']	# HadISDH (adj), CRUTEM, GHCNM, adj-raw		   
##ColorArr=['DimGrey','DodgerBlue','Red','Gold']	# HadISDH (adj), GHCNM, CRUTEM, adj-raw		   
##ColorArr=['Red','DodgerBlue','DimGrey','Gold']	# CRUTEM, GHCNM,HadISDH (adj),  adj-raw		   

# Infiles
INFIL1='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf.nc'
##INFIL1='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies.nc'
INFIL2='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies.nc'
##INFIL2='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/GHCNM_18802014.nc'
INFIL3='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/GHCNM_18802014.nc'
##INFIL3='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies.nc'
##INFIL3='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf.nc'
INFILRAW='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_cf.nc'
infilee=list([INFIL1,INFIL2,INFIL3,INFILRAW])

# What gridbox to pull out (use long (-180 to 180) and lat (-90 to 90) centres?
# If this is greater than one element then make multiple time series
# HadISDH vs CRUTEM
#nGBs=13
#LongPoints=list([-57.500,157.500,167.500,-67.500,-67.500,-132.500,147.500,-62.500,-82.500,122.500,147.500,112.500,112.500])
#LatPoints=list([7.500,7.500,-47.500,-52.500,-17.500,-22.500,-22.500,-42.500,-7.500,-27.500,-17.500,-67.500,-27.500])
#OtherInfo=list([-21.98,-3.17,-2.35,-0.77,-0.72,-0.57,-0.45,-0.41,-0.35,-0.18,-0.16,-0.09,-0.05])	# could be trends or correlations or ratios
#nGBs=10
#LongPoints=list([107.5,102.5,-67.5,37.5,-72.5,-67.5,122.5,-62.5,-57.5,-177.5])
#LatPoints=list([-7.5,-2.5,7.5,2.5,-17.5,-2.5,7.5,2.5,7.5,-17.5])
#OtherInfo=list([0.29,0.31,0.36,0.37,0.38,0.39,0.42,0.42,0.46,0.47])	# could be trends or correlations or ratios

# HadISDH vs GHCNM
#nGBs=11
#LongPoints=list([37.500,-172.500,177.500,167.500,-167.500,147.500,-67.500,-157.500,112.500,112.500,147.500])
#LatPoints=list([-67.500,57.500,-12.500,-47.500,62.500,-22.500,-52.500,17.500,-67.500,-27.500,-17.500])
#OtherInfo=list([-4.74,-3.15,-2.59,-2.12,-1.82,-0.56,-0.36,-0.29,-0.08,-0.08,-0.07])
#nGBs=10
#LongPoints=list([17.5,-72.5,-67.5,-77.5,17.5,37.5,-37.5,12.5,-37.5,82.5])
#LatPoints=list([-27.5,-17.5,7.5,2.5,2.5,-2.5,-7.5,-2.5,-12.5,22.5])
#OtherInfo=list([0.28,0.46,0.48,0.48,0.48,0.53,0.55,0.58,0.60,0.61])

# CRUTEM4 vs GHCNM
nGBs=8
LongPoints=list([177.500,-167.500,-82.500,-132.500,-172.500,157.500,127.500,-57.500])
LatPoints=list([-12.500,62.500,-7.500,-22.500,57.500,7.500,-17.500,7.500])
OtherInfo=list([-4.18,-1.73,-1.37,-0.56,-0.42,-0.26,-0.25,-0.06])
#nGBs=10
#LongPoints=list([-77.5,17.5,137.5,37.5,177.5,-132.5,-62.5,157.5,-67.5,-77.5])
#LatPoints=list([2.5,-27.5,7.5,-2.5,-12.5,-22.5,-2.5,7.5,7.5,17.5])
#OtherInfo=list([0.45,0.56,0.61,0.64,0.65,0.65,0.65,0.65,0.66,0.67])

# Outfiles
#OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/GBTimeSeriesHadISDH.landT.2.0.1.2014pCRUTEM4.3.0.0_' # add long box, lat box and .eps, .png
#OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/GBTimeSeriesHadISDH.landT.2.0.1.2014pGHCNM3_' # add long box, lat box and .eps, .png
OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/GBTimeSeriesCRUTEM4.3.0.0GHCNM3_' # add long box, lat box and .eps, .png

#OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/GBTimeSeriesCORHadISDH.landT.2.0.1.2014pCRUTEM4.3.0.0_' # add long box, lat box and .eps, .png
#OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/GBTimeSeriesCORHadISDH.landT.2.0.1.2014pGHCNM3_' # add long box, lat box and .eps, .png
#OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/GBTimeSeriesCORCRUTEM4.3.0.0GHCNM3_' # add long box, lat box and .eps, .png

# Set up time info
# HadISDH vs CRUTEM
styr=list([1973,1973,1850,1880,1973])	# actual start year, start year of source 1, start year of source 2 etc.
edyr=list([2014,2014,2014,2014,2014])	# actual end year, end year of source 1, end year of source 2 etc.
# HadISDH vs GHCNM
##styr=list([1973,1973,1880,1850,1973])	   # actual start year, start year of source 1, start year of source 2 etc.
##edyr=list([2014,2014,2014,2014,2014])	   # actual end year, end year of source 1, end year of source 2 etc.
# CRUTEM4 vs GHCNM
##styr=list([1973,1850,1880,1973,1973])	   # actual start year, start year of source 1, start year of source 2 etc.
##edyr=list([2014,2014,2014,2014,2014])	   # actual end year, end year of source 1, end year of source 2 etc.
nyrs=(edyr[0]-styr[0])+1
nmons=(nyrs)*12
climst=1981
climed=2010
stcl=climst-styr[0]
edcl=climed-styr[0]

# Set up variables
mdi=-1e30

GBts=[]	# nvars(rows) by 4 regions, by mons masked array
GBtrends=[] # nvars(rows) by 4 regions, by mons masked array

#************************************************************************
# Subroutines
#************************************************************************
# ExtractGB
def ExtractGB(FileName,LongPoint,LatPoint,GBarr,YrStart,YrEnd,NameSources,ThisSource):
    ''' Read in netCDF grid and pull out specific gridbox '''
    f=netcdf.netcdf_file(FileName,'r')
    if (NameSources[ThisSource] == 'HadISDH.2.0.1.2014p'):
        var=f.variables['t_anoms']
    else:
        var=f.variables['temperature_anomaly']
    tmpvar=np.array(np.transpose(var.data[:]))
    f.close()
    
    lnn=np.int(np.floor((LongPoint-(-180))/5.))
    ltt=np.int(np.floor((LatPoint-(-90))/5.))
    st=(YrStart[0]-YrStart[ThisSource+1])*12
    ed=((YrEnd[0]+1)-YrStart[ThisSource+1])*12
    
    print(lnn,ltt,st,ed)
    
    GBarr=tmpvar[lnn,ltt,st:ed+1]

    return GBarr # ExtractGB

#************************************************************************
# CALCANOMS
def CalcAnoms(GBarr,StClim,EdClim,TheMDI):
    ''' Calculate climatology for each month if 50% of months present '''
    ''' Get climate anomalies for all months '''
    GBarr=np.reshape(GBarr,(len(GBarr)/12,12))
    for mm in range(12):
        subarr=GBarr[:,mm]
	climarr=GBarr[StClim:EdClim+1,mm]
	climgots=np.where(climarr != TheMDI)[0]
	gots=np.where(subarr != TheMDI)[0]
	if len(climgots) >= 15:
	    subarr[gots]=subarr[gots]-np.mean(climarr[climgots])
	else:
	    subarr[:]=TheMDI
	GBarr[:,mm]=subarr
    GBarr=np.reshape(GBarr,np.size(GBarr))	
	
    return GBarr # CalcAnoms

#************************************************************************
# AVERAGEANNUALS
def AverageAnnuals(GBarr,TheMDI):
    ''' Calculate the annual mean anomalies where at least 75 % of months are present '''
    AnnArr=np.empty(len(GBarr)/12)
    AnnArr.fill(TheMDI)
    GBarr=np.reshape(GBarr,(len(GBarr)/12,12))
    for yy in range(len(AnnArr)):
        subarr=GBarr[yy,:]
	gots=np.where(subarr != TheMDI)[0]
	if len(gots) >= 9:
	    AnnArr[yy]=np.mean(subarr[gots])
	    
    return AnnArr # AverageAnnuals

#************************************************************************
# FITTREND
def FitTrend(GBarr,TypeTime,TheMDI):
    ''' Use MedianPairwise to fit a trend '''
    ''' return the 5th, median and 95th percentile slopes '''
    ''' If annual then multiply by 10, if monthly then multiply by 120 '''
    TrendStats=np.empty(3)
    TrendStats=MedianPairwise(GBarr,TheMDI,TrendStats)
    
    if TypeTime == 'monthly':
        TrendStats=np.array(TrendStats)*120.
    else:
        TrendStats=np.array(TrendStats)*10.
	
    return TrendStats # FitTrend

#************************************************************************
# GETCORR
def GetCorr(ALLarr,TypeTime,TheMDI,TheSCount):
    ''' Get simple correlation over months present '''
    CorrCount=sum(range(TheSCount)) # number of unique correlating pairs
    CorrStats=np.empty((CorrCount))
    point1=0	# source number
    point2=1 	# source number
    for cloo in range(CorrCount):
        print('Source points: ',point1,point2)
	tmp1=ALLarr[point1,:]
	tmp2=ALLarr[point2,:]
	gots=np.where((tmp1 > TheMDI) & (tmp2 > TheMDI))[0]
	if (len(gots) > 60):
            CorrStats[cloo]=np.corrcoef(tmp1[gots],tmp2[gots])[0,1]
        else:
	    CorrStats[cloo]=0.0
	point2=point2+1
	if (point2 == TheSCount):	# loop through
	    point1=point1+1
	    point2=1+point1
	
    return CorrStats # GetCorr

#************************************************************************
# PlotTimeSeries
def PlotTimeSeries(TheFile,allGB,PointLong,PointLat,trendGB,corrGB,ExtraInfo,TheUnits,
                   MonCount,YrCount,TypeTime,YrStart,YrEnd,TheMDI,SourceCount,Titlee,NameSource,ShortNameSource):
    ''' Plot a single panel and overlay the time series '''
    ''' Add a legend with colours and sources '''
    ''' Annotate with trends and otherinfo '''
    ''' Save as png and eps '''
     
    # for clearer plotting
    MonCount=MonCount+1
    YrCount=YrCount+1    
    
    # Mask the time series
    MSKallGB=ma.masked_where(allGB <= -999.,allGB)	# a fudge here as floating point precision issues in MDI
    
    # set up x axes
    if TypeTime == 'monthly':
        TheMonths=[]
        yr=YrStart[0]
        mon=1
        for m in range(MonCount):
            TheMonths.append(dt.date(yr,mon,1))
	    mon=mon+1
	    if mon == 13:
	        mon=1
	        yr=yr+1   
        TheMonths=np.array(TheMonths)	    
    else:
        TheMonths=[]
        yr=YrStart[0]
        mon=1
        for y in range(YrCount):
            TheMonths.append(dt.date(yr,mon,1))
	    yr=yr+1   
        TheMonths=np.array(TheMonths)	    
    
    xtitlee='Years'
    ytitlee='Anomalies ('+TheUnits+')'
            
    # set up dimensions and plot - this is a 3 column nvar rows plot
    f,ax=plt.subplots(2,figsize=(8,6))	#6,18
    ax[0].set_position([0.1,0.32,0.85,0.58])
    if TypeTime == 'monthly':
	ax[0].set_xlim([TheMonths[0],TheMonths[MonCount-1]])
    else:
	ax[0].set_xlim([TheMonths[0],TheMonths[YrCount-1]])
    #ax.set_xticklabels([])    
    miny=(math.floor(10.*np.min(MSKallGB)))/10.  
    maxy=(math.ceil(10.*np.max(MSKallGB)))/10. 
    ax[0].set_ylim([miny,maxy])
    
    cols=['DimGrey','Red','DodgerBlue','DarkOrange']			   
    if TypeTime == 'monthly':
        for sn in range(SourceCount):
	    ax[0].plot(TheMonths[0:len(TheMonths)-1],MSKallGB[sn,:],c=cols[sn],linewidth=1)
    else:
        for sn in range(SourceCount):
	    ax[0].plot(TheMonths[0:len(TheMonths)-1],MSKallGB[sn,:],c=cols[sn],linewidth=2)
	
    
    ax[0].set_title(Titlee+' Lon: '+'{0:8.3f}'.format(PointLong)+' Lat: '+'{0:8.3f}'.format(PointLat)+' Ratio: '+'{0:6.2f}'.format(ExtraInfo),
                 fontdict=None,loc='center',size=12)
#    ax[0].set_title(Titlee+' Lon: '+'{0:8.3f}'.format(PointLong)+' Lat: '+'{0:8.3f}'.format(PointLat)+' Correlation: '+'{0:6.2f}'.format(ExtraInfo),
#                 fontdict=None,loc='center',size=12)

# get linear trend and annotate (does this work with masked arrays?)
    if trendGB[0,0] != TheMDI:
	scaletrend=([0.94,0.88,0.82])
	for sn in range(SourceCount):
	    linstr='{0:5.2f}'.format(trendGB[sn,0])+' ('+'{0:5.2f}'.format(trendGB[sn,1])+'  to '+'{0:5.2f}'.format(trendGB[sn,2])+') '+'{:s}'.format(TheUnits)+' decade$^{-1}$'
	    #linstr="%19s %5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (NameSource[sn],trendGB[sn,0],trendGB[sn,1],trendGB[sn,2],TheUnits)
	    ax[0].annotate(ShortNameSource[sn],xy=(0.02,scaletrend[sn]),xycoords='axes fraction',size=12,ha='left',color=cols[sn])
	    ax[0].annotate(linstr,xy=(0.4,scaletrend[sn]),xycoords='axes fraction',size=12,ha='left',color=cols[sn])


# get correlation and annotate (does this work with masked arrays?)
    if corrGB[0] != TheMDI:
	scaletrend=([0.02,0.08,0.14,0.2,0.26,0.32])
	point1=0
	point2=1
	for sn in range(len(corrGB[:])):
	    linstr='r = '+'{0:5.2f}'.format(corrGB[sn])
	    #linstr="%19s %5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (NameSource[sn],trendGB[sn,0],trendGB[sn,1],trendGB[sn,2],TheUnits)
	    ax[0].annotate(ShortNameSource[point1]+','+ShortNameSource[point2],xy=(0.6,scaletrend[sn]),xycoords='axes fraction',size=10,ha='left',color='Black')
	    ax[0].annotate(linstr,xy=(0.85,scaletrend[sn]),xycoords='axes fraction',size=10,ha='left',color='Black')
	    point2=point2+1
	    if (point2 == SourceCount):	# loop through
	        point1=point1+1
	        point2=1+point1

    if TypeTime == 'monthly':
        ax[0].hlines(0,TheMonths[0],TheMonths[MonCount-1],
	        color='black',linewidth=0.5)
    else:
	ax[0].hlines(0,TheMonths[0],TheMonths[YrCount-1],
	        color='black',linewidth=0.5)
	 
    ax[0].set_ylabel(ytitlee,fontsize=12)
	
# add second axis
    ax[1].set_position([0.1,0.12,0.85,0.20])
    if TypeTime == 'monthly':
	ax[1].set_xlim([TheMonths[0],TheMonths[MonCount-1]])
    else:
	ax[1].set_xlim([TheMonths[0],TheMonths[YrCount-1]])
    #ax.set_xticklabels([])    
    miny=(math.floor(10.*np.min(MSKallGB[3,:])))/10.  
    maxy=(math.ceil(10.*np.max(MSKallGB[3,:])))/10. 
    ax[1].set_ylim([miny,maxy])

    lowess = sm.nonparametric.lowess
    gots=np.where(allGB[3,:] > TheMDI)[0]
    
    if TypeTime == 'monthly':
        ax[1].plot(TheMonths[0:len(TheMonths)-1],MSKallGB[3,:],c='Black',linewidth=0.5)
        Smoothie=lowess(allGB[3,gots],gots,frac=(1./(len(TheMonths)-1))*(12*2))
	ax[1].plot(TheMonths[gots],Smoothie[:,1],c=cols[3],linewidth=1)
	
    else:
        ax[1].plot(TheMonths[0:len(TheMonths)-1],MSKallGB[3,:],c='Black',linewidth=1)
        Smoothie = lowess(allGB[3,gots],gots,frac=(1./(len(TheMonths)-1))*2)
        ax[1].plot(TheMonths[gots],Smoothie[:,1],c=cols[3],linewidth=2)

    if TypeTime == 'monthly':
        ax[1].hlines(0,TheMonths[0],TheMonths[MonCount-1],
	        color='black',linewidth=0.5)
    else:
	ax[1].hlines(0,TheMonths[0],TheMonths[YrCount-1],
	        color='black',linewidth=0.5)

    ax[1].annotate(ShortNameSource[3],xy=(0.02,0.8),xycoords='axes fraction',size=12,ha='left',color=cols[3])

    ytitlee='Difference ('+TheUnits+')'

    ax[1].set_ylabel(ytitlee,fontsize=12)
    ax[1].set_xlabel(xtitlee,fontsize=12)
     
    
# Figure Watermark and Labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    

    #plt.show()
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")

   # stop
     
    return #PlotTimeSeries
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# Loop through each gridbox

for gb in range(nGBs):
    # Can now establish the filename
    plottee=OUTFIL+'{:s}'.format(str(LongPoints[gb]))+'_'+'{:s}'.format(str(LatPoints[gb]))
    # Set up timeseries arrays
    if (timetype == 'monthly'):
        GBts=np.empty((nsources+1,nmons))
    else:
        GBts=np.empty((nsources+1,nyrs))
        
    GBts.fill(mdi)
    GBtrends=np.zeros((3,3)) # row for each source, columns for trend, 5thpct, 95thpct
    GBtrends.fill(mdi)

    GBcorrs=np.zeros((sum(range(nsources)))) # row for each source, columns for trend, 5thpct, 95thpct
    GBcorrs.fill(mdi)
    
# Loop through sources
    for ns in range(nsources+1):
        # set up tmparr for gridbox data
        tmpGB=np.empty(nmons)
	tmpGB.fill(mdi)
	
# Read in netCDF set of grids and extract gridbox
        tmpGB=ExtractGB(infilee[ns],LongPoints[gb],LatPoints[gb],tmpGB,styr,edyr,SourceNames,ns)
	bads=np.where(tmpGB < -999.)[0]
	tmpGB[bads]=mdi
	
# Reanomalies to climatology period
        tmpGB=CalcAnoms(tmpGB,stcl,edcl,mdi) # should modify tmpGB
	# may need to do something special with HadISDH raw to make it comparable with adj
	
# If timetype is annual then average to annual	
        if (timetype == 'annual'):
	    GBts[ns,:]=AverageAnnuals(tmpGB,mdi)
	else:
	    GBts[ns,:]=tmpGB
	        
# Fit a linear trend if desired
        if (fitlin) & (ns < nsources):
	    GBtrends[ns,:]=FitTrend(GBts[ns,:],timetype,mdi)

# Now make the fourth time series HadISDH adj-raw
    gots=np.where((GBts[0,:] > mdi) & (GBts[3,:] > mdi))[0]
    NEWts=np.empty(len(GBts[0,:]))
    NEWts.fill(mdi)
    NEWts[gots]=GBts[0,gots]-GBts[3,gots]
    # make mean of the diff series zero over the last two years of present data - no adj there!
    Meanie=np.mean(NEWts[gots[len(gots)-13:len(gots)]])
    NEWts[gots]=NEWts[gots]-Meanie
    GBts[3,:]=NEWts

# Now all sources are read in get corrs if desired
    if (addcor):
        GBcorrs=GetCorr(GBts[0:nsources,:],timetype,mdi,nsources)

	  
# Pass all sources and trends and otherinfo to plotter
    print('Plotting...',gb)
    PlotTimeSeries(plottee,GBts,LongPoints[gb],LatPoints[gb],GBtrends,GBcorrs,OtherInfo[gb],
                       unitees,nmons,nyrs,timetype,styr,edyr,mdi,nsources,Namey,SourceNames,ShortSourceNames)

#    if gb == 1:
#        stop()

print("And, we are done!")

