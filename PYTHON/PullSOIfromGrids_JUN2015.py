#!/usr/local/sci/bin/python

#***************************************
#11 Jun 2015 KMW - v1
# Takes any gridded field of mnnthly mean sea level pressure
# Can have long/lat/time manually or read in
# Pulls out gridbox closest to Darwin and Tahiti
# Calculates SOI for series
# Saves SOI time series to file

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PullSOIfromGrids_JUN2015.py
#
# REQUIRES
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

# Set up initial run choices - if these aren't set they will be taken from the NetCDF
StYr=1973
EdYr=2014
StMon=1
EdMon=12
NYrs=(EdYr+1)-StYr
NMons=NYrs*12

NLats=36
NLons=72
StLat=-90
StLon=-180
LatGB=180./NLats
LonGB=360./Nlons
Lats=np.arange(StLat+(LatGB/2.),StLat+180+(LatGB/2.),LatGB)
Lons=np.arange(StLon+(LonGB/2.),StLon+360+(LonGB/2.),LonGB)

# Files
INFIL='/data/local/hadkw/CMIP5/HadGEM2ES/ATMOSAMON/r2i1p1/psl/psl_Amon_HadGEM2-ES_historicalFULL_r2i1p1_1973-2014.nc'
OUTFIL='/data/local/hadkw/CMIP5/HadGEM2ES/STATISTICS/SOI_psl_Amon_HadGEM2-ES_historicalFULL_r2i1p1_1973-2014.nc'

# Lat and Lon of Darwin and Tahiti
Darwin=([130.8333,-12.4500])
Tahiti=([-149.4167,-17.6667])


#************************************************************************
# Subroutines
#************************************************************************
# ExtractDetails
def ExtractDetails(FileName,LongPoint,LatPoint,GBarr,TheLonCount,TheLatCount,TheLats,TheLons,TheStYr,TheEdYr,TheStMon,TheEdMon,TheYCount,TheMCount):
    ''' Read in netCDF '''
    ''' If any details are missing, find them, convert lats and lons as appropriate '''
    ''' Pull out data and convert to match converted lats and lons grid and pull out info '''
    ''' Find desired gridbox and pull out time series '''

    f=netcdf.netcdf_file(FileName,'r')

    var=f.variables['psl']
    pslvar=np.array(np.transpose(var.data[:]))

    if (TheStYr -- 0):
        var=f.variables['lat']
        latvar=np.array(np.transpose(var.data[:]))

        var=f.variables['lon']
        lonvar=np.array(np.transpose(var.data[:]))
    
        var=f.variables['time']
        timvar=np.array(np.transpose(var.data[:]))
        timunits=var.units
        Ypoint=int(timunits[11:15])
        Mpoint=int(timunits[16:18])
        Dpoint=int(timunits[19:21])
	
	# now use Datetime to find the actual start year/month and then get Year count and Month Count
	# REMEBER THAT MODELS HAVE 30 DAY MONTHS!!!! SO DATETIME WON'T WORK!!!
	
    # set up x axes
    if TheTimeType == 'monthly':
        TheMonths=[]
        yr=TheStYr
        mon=1
        for m in range(TheMCount):
            TheMonths.append(dt.date(yr,mon,1))
	    mon=mon+1
	    if mon == 13:
	        mon=1
	        yr=yr+1   
        TheMonths=np.array(TheMonths)	    
    else:
        TheMonths=[]
        yr=TheStYr
        mon=1
        for y in range(TheYCount):
            TheMonths.append(dt.date(yr,mon,1))
	    yr=yr+1   
        TheMonths=np.array(TheMonths)	    

 
    f.close()
    
    if (lonvar[0] > 100):	# convert to -180 to 180
        lonvar=lonvar
	pslvar=np.reshape()
	
    if (latvar[0] > 0): 	# convert to -90 to 90	
	latvar=latvar
	pslvar=np.reshape()
	
    # get the time details.	
    lnn=np.int(np.floor((LongPoint-(TheLons[0]))/abs(TheLons[1]-TheLons[0])))
    ltt=np.int(np.floor((LatPoint-(TheLats[0]))/abs(TheLats[1]-TheLats[0])))
    
    print(lnn,ltt)
    
    GBarr=tmpvar[lnn,ltt,:]

    return GBarr,TheLonCount,TheLatCount,TheLats,TheLons,TheStYr,TheEdYr,TheYCount,TheMCount # ExtractGB
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

