#!/usr/local/sci/bin/python

#***************************************
# 27 April 2015 KMW - v1
# For TEMPERATURE ONLY
# Plots a time series for any or all regions
# Option to overplot the same region from another source too 
# Option to add trends of each source
# Option to add correlation of each source
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotAllRegionsTimeSeries_APR2015.py
#
# REQUIRES
# RandomsRanges.py
# LinearTrends.py
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import numpy.ma as ma
import numpy as np
import sys, os
import scipy.stats
import struct
import os.path
import math
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
plotdiff=True

# What are the units?
unitees='$^{o}$C'  #'deg C'

# Number of sources to plot
nsources=6 	# this can be one to three

# Source Names for Title
Namey=['Globe (70$^{o}$S-70$^{o}$N)','N. Hemisphere (20$^{o}$N-70$^{o}$N)','Tropics (20$^{o}$S-20$^{o}$N)','S. Hemisphere (70$^{o}$S-20$^{o}$S)']
SourceNames=['HadISDH.2.0.1.2014p (adj)','HadISDH.2.0.1.2014p (raw)','CRUTEM4.3.0.0','GHCNM3','GISTEMP','BERKELEY EARTH']
ShortSourceNames=['HadISDH adj','HadISDH raw','CRUTEM','GHCNM','GISS','BERKELEY']

# Infiles
INFIL1='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_areaTS_19732014.nc'
INFIL2='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/HadISDH.landT.2.0.1.2014p_FLATgridRAW5by5_JAN2015_areaTS_19732014.nc'
INFIL3='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/CRUTEM4_T_areaTS_19732014.nc' #CRUTEM.4.3.0.0.anomalies_areaTS_18502014.nc
INFIL4='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/GHCNM_T_areaTS_19732014.nc' #GHCNM_18802014_areaTS_18802014.nc
INFIL5='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/GISS_T_areaTS_19732014.nc'
INFIL6='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/BERKELEY_T_areaTS_19732014.nc'
INFIL9='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/CRUTEM4_T_HadISDHMASKareaTS_19732014.nc' #CRUTEM.4.3.0.0.anomalies_areaTSMSK7605_19732014.nc
INFIL10='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/GHCNM_T_HadISDHMASKareaTS_19732014.nc' #GHCNM_18802014_areaTSMSK7605_19732014.nc
INFIL11='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/GISS_T_HadISDHMASKareaTS_19732014.nc'
INFIL12='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/BERKELEY_T_HadISDHMASKareaTS_19732014.nc'
infilee=list([INFIL1,INFIL2,INFIL3,INFIL4,INFIL5,INFIL6,INFIL1,INFIL2,INFIL9,INFIL10,INFIL11,INFIL12])

# What gridbox to pull out (use long (-180 to 180) and lat (-90 to 90) centres?
# If this is greater than one element then make multiple time series
nREGs=4
Regions=list(['glob_T_anoms','nhem_T_anoms','trop_T_anoms','shem_T_anoms'])
RegionsOTH=list(['glob_anoms','nhem_anoms','trop_anoms','shem_anoms'])
OtherInfo=list([])	# could be trends or correlations or ratios

# Outfiles
if (plotdiff):
    OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/AllRegionSourceTimeSeriesDIFFSHadISDH.landT.2.0.1.2014p' # add long box, lat box and .eps, .png
else:
    OUTFIL='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/AllRegionSourceTimeSeriesHadISDH.landT.2.0.1.2014p' # add long box, lat box and .eps, .png
    OUTFILC='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/AllCorrsMatrixHadISDH.landT.2.0.1.2014p' # correlation matrix for each region

# Set up time info
#styr=list([1973,1973,1973,1850,1880,1973,1973,1973,1973,1973,1973,1973,1973])	# actual start year, start year of source 1, start year of source 2 etc.
styr=list([1973,1973,1973,1973,1973,1973,1973,1973,1973,1973,1973,1973,1973])	# actual start year, start year of source 1, start year of source 2 etc.
edyr=list([2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014])	# actual end year, end year of source 1, end year of source 2 etc.
#styr=list([1973,1973,1973,1973,1973])	# actual start year, start year of source 1, start year of source 2 etc.
#edyr=list([2014,2014,2014,2014,2014])	# actual end year, end year of source 1, end year of source 2 etc.
nyrs=(edyr[0]-styr[0])+1
nmons=(nyrs)*12
climst=1981
climed=2010
stcl=climst-styr[0]
edcl=climed-styr[0]

# Set up variables
mdi=-1e30

#************************************************************************
# Subroutines
#************************************************************************
# ExtractREG
def ExtractREG(FileName,TheRegions,REGarr,YrStart,YrEnd,TheRCount,ThisSource):
    ''' Read in netCDF grid and pull out specific region time series '''

    st=(YrStart[0]-YrStart[ThisSource+1])*12
    # Check when source start is later than desired start
    if st < 0:
       st=0
    ed=((YrEnd[0]+1)-YrStart[ThisSource+1])*12
    
    print(st,ed)

    f=netcdf.netcdf_file(FileName,'r')
    for loo in range(TheRCount):
        var=f.variables[TheRegions[loo]]
        tmpvar=np.array(np.transpose(var.data[:]))
	# Can cope when source begins earlier or later than desired start
	REGarr[loo,len(REGarr[loo,:])-(ed-st):]=tmpvar[st:ed+1]    
    f.close()
    
    return REGarr # ExtractREG

#************************************************************************
# CALCANOMS
def CalcAnoms(REGarr,StClim,EdClim,TheMDI,TheRCount):
    ''' Calculate climatology for each month if 50% of months present '''
    ''' Get climate anomalies for all months '''
    for loo in range(TheRCount):
        tmparr=np.reshape(REGarr[loo,:],(len(REGarr[loo,:])/12,12))
        for mm in range(12):
            subarr=tmparr[:,mm]
	    climarr=tmparr[StClim:EdClim+1,mm]
	    climgots=np.where(climarr != TheMDI)[0]
	    gots=np.where(subarr != TheMDI)[0]
	    if len(climgots) >= 15:
	        subarr[gots]=subarr[gots]-np.mean(climarr[climgots])
	    else:
	        subarr[:]=TheMDI
	    tmparr[:,mm]=subarr
        REGarr[loo,:]=np.reshape(tmparr,np.size(tmparr))	
	
    return REGarr # CalcAnoms

#************************************************************************
# AVERAGEANNUALS
def AverageAnnuals(REGarr,TheMDI,TheRCount):
    ''' Calculate the annual mean anomalies where at least 75 % of months are present '''
    AnnArr=np.empty((TheRCount,len(REGarr[0,:])/12))
    AnnArr.fill(TheMDI)
    for loo in range(TheRCount):
        tmparr=np.reshape(REGarr[loo,:],(len(REGarr[loo,:])/12,12))
        for yy in range(len(AnnArr[loo,:])):
            subarr=tmparr[yy,:]
	    gots=np.where(subarr != TheMDI)[0]
	    if len(gots) >= 9:
	        AnnArr[loo,yy]=np.mean(subarr[gots])
	    
    return AnnArr # AverageAnnuals

#************************************************************************
# FITTREND
def FitTrend(REGarr,TypeTime,TheMDI,TheRCount):
    ''' Use MedianPairwise to fit a trend '''
    ''' return the 5th, median and 95th percentile slopes '''
    ''' If annual then multiply by 10, if monthly then multiply by 120 '''
    TrendStats=np.empty((TheRCount,3))
    for loo in range(TheRCount):
        TrendStats[loo,:]=MedianPairwise(REGarr[loo,:],TheMDI,TrendStats[loo,:])
        if TypeTime == 'monthly':
            TrendStats[loo,:]=np.array(TrendStats[loo,:])*120.
        else:
            TrendStats[loo,:]=np.array(TrendStats[loo,:])*10.
	
    return TrendStats # FitTrend

#************************************************************************
# GETCORR
def GetCorr(ALLarr,TypeTime,TheMDI,TheRCount,TheSCount):
    ''' Get simple correlation over months present '''
    CorrCount=sum(range(TheSCount)) # number of unique correlating pairs
    CorrStats=np.empty((CorrCount,TheRCount))
    point1=0	# source number
    point2=1 	# source number
    for cloo in range(CorrCount):
        print('Source points: ',point1,point2)
	for rloo in range(TheRCount):
	    tmp1=ALLarr[point1,rloo,:]
	    tmp2=ALLarr[point2,rloo,:]
	    gots=np.where((tmp1 > TheMDI) & (tmp2 > TheMDI))
            CorrStats[cloo,rloo]=np.corrcoef(tmp1[gots],tmp2[gots])[0,1]
        point2=point2+1
	if (point2 == TheSCount):	# loop through
	    point1=point1+1
	    point2=1+point1
	
    return CorrStats # GetCorr

#************************************************************************
# PlotTimeSeries
def PlotTimeSeries(TheFile,ALLarr,trendALL,corrALL,ALLarrMSK,trendALLMSK,corrALLMSK,ExtraInfo,TheUnits,
                   MonCount,YrCount,TypeTime,YrStart,YrEnd,TheMDI,SourceCount,
		   Titlee,ShortNameSource,TheRCount):
    ''' Plot a multi-panel and overlay the time series '''
    ''' Add a legend with colours and sources '''
    ''' Annotate with trends and otherinfo '''
    ''' This will be an 8 panel plot - Non-masked on the left, masked on the right '''
    ''' Save as png and eps '''
     
    # for clearer plotting
    MonCount=MonCount+1
    YrCount=YrCount+1    
    
    # Mask the time series
    MSKALLarr=ma.masked_where(ALLarr <= -999.,ALLarr)	# a fudge here as floating point precision issues in MDI
    MSKALLarrMSK=ma.masked_where(ALLarrMSK <= -999.,ALLarrMSK)	# a fudge here as floating point precision issues in MDI
    
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

    # set up number of panels and number of lines
    nplots=TheRCount
    print('PLOT NUMBERS: ',nplots)
    nlines=[]
    for n in range(nplots):
        nlines.append(SourceCount)
	
    Letteree=[]
    Letteree=LetterRange(0,nplots*2)
    
    # set up dimensions and plot - this is a 2 column nvar rows plot
    xpos=[]
    ypos=[]
    xfat=[]
    ytall=[]
    totalyspace=0.90	# start 0.08 end 0.98
    totalxspace=0.40	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.07)
	ypos.append(0.98-((n+1)*(totalyspace/nplots)))
	xfat.append(totalxspace)
	ytall.append(totalyspace/nplots)

            
    # set up dimensions and plot - this is a 3 column nvar rows plot

    if (SourceCount == 6):
        cols=['DimGrey','Firebrick','DodgerBlue','DarkOrange','MediumBlue','HotPink']			   
    else:
        cols=['Firebrick','DodgerBlue','DarkOrange','MediumBlue','HotPink']			   
    
    f,axarr=plt.subplots(TheRCount*2,figsize=(16,12),sharex=False)	#6,18

    for pp in range(nplots):
        print('Plot: ',pp)
	#axarr[pp].set_size(14)
	axarr[pp].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        if TypeTime == 'monthly':
	    axarr[pp].set_xlim([TheMonths[0],TheMonths[MonCount-1]])
        else:
	    axarr[pp].set_xlim([TheMonths[0],TheMonths[YrCount-1]])
	if pp < nplots-1:
	    axarr[pp].set_xticklabels([])    
        miny=(math.floor(10.*np.min(np.append(MSKALLarr,MSKALLarrMSK))))/10.  
        maxy=(math.ceil(10.*np.max(np.append(MSKALLarr,MSKALLarrMSK))))/10. 
        axarr[pp].set_ylim([miny,maxy])

        if TypeTime == 'monthly':
            for sn in range(SourceCount):
	        axarr[pp].plot(TheMonths[0:len(TheMonths)-1],MSKALLarr[sn,pp,:],c=cols[sn],linewidth=1)
        else:
            for sn in range(SourceCount):
	        axarr[pp].plot(TheMonths[0:len(TheMonths)-1],MSKALLarr[sn,pp,:],c=cols[sn],linewidth=2)
	
	axarr[pp].annotate(Letteree[pp]+') '+Titlee[pp],xy=(0.02,0.90),xycoords='axes fraction',size=16)

# get linear trend and annotate (does this work with masked arrays?)
        if trendALL[0,pp,0] != TheMDI:
	    scaletrend=([0.32,0.26,0.2,0.14,0.08,0.02])    #([0.88,0.82,0.76,0.70])
	    if (SourceCount == 5):
	        scaletrend=scaletrend=([0.26,0.2,0.14,0.08,0.02])
	    for sn in range(SourceCount):
	        linstr='{0:5.2f}'.format(trendALL[sn,pp,0])+' ('+'{0:5.2f}'.format(trendALL[sn,pp,1])+','+'{0:5.2f}'.format(trendALL[sn,pp,2])+') '+'{:s}'.format(TheUnits)+' dec$^{-1}$'
	        namstr=ShortNameSource[sn]+' ('+'{0:5.2f}'.format(np.std(MSKALLarr[sn,pp,:]))+')'
		#linstr="%19s %5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (NameSource[sn],trendGB[sn,0],trendGB[sn,1],trendGB[sn,2],TheUnits)
	        axarr[pp].annotate(namstr,xy=(0.5,scaletrend[sn]),xycoords='axes fraction',size=9,ha='left',color=cols[sn])
	        axarr[pp].annotate(linstr,xy=(0.72,scaletrend[sn]),xycoords='axes fraction',size=9,ha='left',color=cols[sn])

## get correlation and annotate (does this work with masked arrays?)
#        if corrALL[0,pp] != TheMDI:
#	    scaletrend=([0.02,0.08,0.14,0.2,0.26,0.32])
#	    point1=0
#	    point2=1
#	    for sn in range(len(corrALL[:,0])):
#	        linstr='r = '+'{0:6.3f}'.format(corrALL[sn,pp])
#	        #linstr="%19s %5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (NameSource[sn],trendGB[sn,0],trendGB[sn,1],trendGB[sn,2],TheUnits)
#	        axarr[pp].annotate(ShortNameSource[point1]+','+ShortNameSource[point2],xy=(0.6,scaletrend[sn]),xycoords='axes fraction',size=10,ha='left',color='Black')
#	        axarr[pp].annotate(linstr,xy=(0.85,scaletrend[sn]),xycoords='axes fraction',size=10,ha='left',color='Black')
#		point2=point2+1
#		if (point2 == SourceCount):	# loop through
#	            point1=point1+1
#	            point2=1+point1
	 
        # Plot a line at zero
	if TypeTime == 'monthly':
	    axarr[pp].hlines(0,TheMonths[0],TheMonths[MonCount-1],color='black',linewidth=0.5)
        else:
	    axarr[pp].hlines(0,TheMonths[0],TheMonths[YrCount-1],color='black',linewidth=0.5)
	
        axarr[pp].set_ylabel(ytitlee,fontsize=12)

    axarr[pp].set_xlabel(xtitlee,fontsize=12)


# NOW DO THE MASKED PLOTS
    # set up dimensions and plot - this is a 2 column nvar rows plot
    xpos=[]
    ypos=[]
    xfat=[]
    ytall=[]
    totalyspace=0.90	# start 0.08 end 0.98
    totalxspace=0.40	# start 0.12 end 0.98
    for n in range(nplots):
        xpos.append(0.57)
	ypos.append(0.98-((n+1)*(totalyspace/nplots)))
	xfat.append(totalxspace)
	ytall.append(totalyspace/nplots)

            
    # set up dimensions and plot - this is a 3 column nvar rows plot


#    f,axarr=plt.subplots(TheRCount*2,figsize=(16,12),sharex=False)	#6,18

    for pp in range(nplots):
        print('Plot: ',pp)
	#axarr[pp].set_size(14)
	axarr[pp+4].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        if TypeTime == 'monthly':
	    axarr[pp+4].set_xlim([TheMonths[0],TheMonths[MonCount-1]])
        else:
	    axarr[pp+4].set_xlim([TheMonths[0],TheMonths[YrCount-1]])
	if pp < nplots-1:
	    axarr[pp+4].set_xticklabels([])    
#        miny=(math.floor(10.*np.min(MSKALLarr)))/10.  
#        maxy=(math.ceil(10.*np.max(MSKALLarr)))/10. 
        axarr[pp+4].set_ylim([miny,maxy])

        if TypeTime == 'monthly':
            for sn in range(SourceCount):
	        axarr[pp+4].plot(TheMonths[0:len(TheMonths)-1],MSKALLarrMSK[sn,pp,:],c=cols[sn],linewidth=1)
        else:
            for sn in range(SourceCount):
	        axarr[pp+4].plot(TheMonths[0:len(TheMonths)-1],MSKALLarrMSK[sn,pp,:],c=cols[sn],linewidth=2)
	
	axarr[pp+4].annotate(Letteree[pp+4]+') '+Titlee[pp],xy=(0.02,0.9),xycoords='axes fraction',size=16)

# get linear trend and annotate (does this work with masked arrays?)
        if trendALLMSK[0,pp,0] != TheMDI:
	    scaletrend=([0.32,0.26,0.2,0.14,0.08,0.02])    #([0.88,0.82,0.76,0.70])
	    if (SourceCount == 5):
	        scaletrend=scaletrend=([0.26,0.2,0.14,0.08,0.02])
	    for sn in range(SourceCount):
	        if (sn > 1):
		    mush=' (masked)'
		else:
		    mush=''
		linstr='{0:5.2f}'.format(trendALLMSK[sn,pp,0])+' ('+'{0:5.2f}'.format(trendALLMSK[sn,pp,1])+','+'{0:5.2f}'.format(trendALLMSK[sn,pp,2])+') '+'{:s}'.format(TheUnits)+' dec$^{-1}$'
	        namstr=ShortNameSource[sn]+' ('+'{0:5.2f}'.format(np.std(MSKALLarrMSK[sn,pp,:]))+')'
		#linstr="%19s %5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (NameSource[sn],trendGB[sn,0],trendGB[sn,1],trendGB[sn,2],TheUnits)
	        axarr[pp+4].annotate(namstr,xy=(0.5,scaletrend[sn]),xycoords='axes fraction',size=9,ha='left',color=cols[sn])
	        axarr[pp+4].annotate(linstr,xy=(0.72,scaletrend[sn]),xycoords='axes fraction',size=9,ha='left',color=cols[sn])


## get correlation and annotate (does this work with masked arrays?)
#        if corrALL[0,pp] != TheMDI:
#	    scaletrend=([0.02,0.08,0.14,0.2,0.26,0.32])
#	    point1=0
#	    point2=1
#	    for sn in range(len(corrALL[:,0])):
#	        linstr='r = '+'{0:6.3f}'.format(corrALL[sn,pp])
#	        #linstr="%19s %5.2f (%5.2f to %5.2f) %s decade$^{-1}$ " % (NameSource[sn],trendGB[sn,0],trendGB[sn,1],trendGB[sn,2],TheUnits)
#	        axarr[pp+4].annotate(ShortNameSource[point1]+','+ShortNameSource[point2],xy=(0.6,scaletrend[sn]),xycoords='axes fraction',size=10,ha='left',color='Black')
#	        axarr[pp+4].annotate(linstr,xy=(0.85,scaletrend[sn]),xycoords='axes fraction',size=10,ha='left',color='Black')
#		point2=point2+1
#		if (point2 == SourceCount):	# loop through
#	            point1=point1+1
#	            point2=1+point1
	 
        # Plot a line at zero
	if TypeTime == 'monthly':
	    axarr[pp+4].hlines(0,TheMonths[0],TheMonths[MonCount-1],color='black',linewidth=0.5)
        else:
	    axarr[pp+4].hlines(0,TheMonths[0],TheMonths[YrCount-1],color='black',linewidth=0.5)
	
        axarr[pp+4].set_ylabel(ytitlee,fontsize=12)

    axarr[pp+4].set_xlabel(xtitlee,fontsize=12)

    
# Figure Watermark and Labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    

    #plt.show()
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")


     
    return #PlotTimeSeries

#************************************************************************
# PlotCorrMatrix
def PlotCorrMatrix(Filee,CorrAll,CorrAllMSK,TheRegions,TheNames):
    ''' This plots a correlation matrix for each region '''
    ''' Upper triangle is non-masked '''
    ''' Lower triangle is masked to HadISDH '''
    ''' Shades of YlOrRd from 0.8 to 1 '''
    
    # make matrix for each region
    RegCount=len(TheRegions)
    SourceCount=len(TheNames)
    GlobMat=np.empty((SourceCount,SourceCount))
    NHemMat=np.empty((SourceCount,SourceCount))
    TropMat=np.empty((SourceCount,SourceCount))
    SHemMat=np.empty((SourceCount,SourceCount))
    
    np.fill_diagonal(GlobMat,1.)
    np.fill_diagonal(NHemMat,1.)
    np.fill_diagonal(TropMat,1.)
    np.fill_diagonal(SHemMat,1.)

    print('FILLED DIAGONALS')
    
    point1=0
    point2=1
    for cc in range(len(CorrAll)):
        print(cc,point1,point2)
        GlobMat[point1,point2]=CorrAll[cc,0]
        NHemMat[point1,point2]=CorrAll[cc,1]
        TropMat[point1,point2]=CorrAll[cc,2]
        SHemMat[point1,point2]=CorrAll[cc,3]
        GlobMat[point2,point1]=CorrAllMSK[cc,0]
        NHemMat[point2,point1]=CorrAllMSK[cc,1]
        TropMat[point2,point1]=CorrAllMSK[cc,2]
        SHemMat[point2,point1]=CorrAllMSK[cc,3]
        point2=point2+1
	if (point2 == 6):
	    point1=point1+1
	    point2=point1+1

    print('GOT MATRICES')
	    
    # set up colours
    cmap=plt.get_cmap('YlOrRd')
        
    cmaplist=[cmap(i) for i in range(cmap.N)]
    print('CHOSEN COLOURS')
# remove the darkest and lightest (white and black) - and reverse
    for loo in range(30):
        cmaplist.remove(cmaplist[0])
    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    print('REFINED COLOURS')
    
    bounds=np.array([0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1,1.00001])

    strbounds=["%4.3g" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    print('SORTED COLOURS')
    
    # set up plot space and dimensions and plot - this is a 2 column nvar rows plot
#    xpos=list([0.12,0.62,0.12,0.62])
#    ypos=list([0.54,0.54,0.08,0.08])
#    xfat=list([0.36,0.36,0.36,0.36])
#    ytall=list([0.32,0.32,0.32,0.32])
    xpos=list([0.02,0.52,0.02,0.52])
    ypos=list([0.66,0.66,0.2,0.2])
    xfat=list([0.36,0.36,0.36,0.36])
    ytall=list([0.30,0.3,0.3,0.3])

    miny=0 
    maxy=6 
    xinds=range(7)
    yinds=range(6,-1,-1)
    ArrX,ArrY=np.meshgrid(xinds,yinds)

    Letteree=[]
    Letteree=LetterRange(0,len(TheRegions))

    f,axarr=plt.subplots(len(TheRegions),figsize=(12,13),sharex=False)	#6,18

    for pp in range(len(TheRegions)):
        print('Plot: ',pp)
	#axarr[pp].set_size(14)
	axarr[pp].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
        axarr[pp].set_xticklabels([])    
        axarr[pp].set_yticklabels([])    
        axarr[pp].set_ylim([miny,maxy])
        axarr[pp].set_xlim([miny,maxy])
        
	if (pp == 0):
	    ChosenMat=GlobMat
	if (pp == 1):
	    ChosenMat=NHemMat
	if (pp == 2):
	    ChosenMat=TropMat
	if (pp == 3):
	    ChosenMat=SHemMat
	
	print('MAKING THE GRIDS') 
	print(ChosenMat)   
        grids=axarr[pp].pcolor(ArrX,ArrY,ChosenMat,cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5) #,latlon='TRUE'

	for sn in range(len(TheNames)):
	    NamePos=sn*(1./6.)+(1./12.)
	    Topsies=-0.02
	    Sidesies=1.02
	    print(sn,NamePos,Topsies,Sidesies)
	    axarr[pp].annotate(TheNames[sn],xy=(NamePos,Topsies),xycoords='axes fraction',rotation=90,size=12,ha='center',va='top')
	    axarr[pp].annotate(TheNames[sn],xy=(Sidesies,1-NamePos),xycoords='axes fraction',rotation=0,size=12,ha='left',va='center')
	
	axarr[pp].annotate(Letteree[pp]+') '+TheRegions[pp],xy=(0.0,1.03),xycoords='axes fraction',size=16)
    
    # loop through plots (a to d) label sources
    
    # add colour bar along righthandside

    cbax=f.add_axes([0.1,0.05,0.8,0.03])
    cb=f.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=12) 
    cb.ax.set_xticklabels(strbounds)
    plt.figtext(0.5,0.01,'Correlation',size=16,ha='center')
    
    # save to file   
    #plt.show()
    plt.savefig(Filee+".eps")
    plt.savefig(Filee+".png")


    return #PlotCorrMatrix    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# Set up arrays
if (timetype == 'monthly'):
    REGts=np.empty((nsources,nREGs,nmons))
    REGtsMSK=np.empty((nsources,nREGs,nmons))
    DIFFts=np.empty((nsources-1,nREGs,nmons))
    DIFFtsMSK=np.empty((nsources-1,nREGs,nmons))
    moop='MON'
else:
    REGts=np.empty((nsources,nREGs,nyrs))
    REGtsMSK=np.empty((nsources,nREGs,nyrs))
    DIFFts=np.empty((nsources-1,nREGs,nyrs))
    DIFFtsMSK=np.empty((nsources-1,nREGs,nyrs))
    moop='ANN'

REGts.fill(mdi)
REGtsMSK.fill(mdi)

DIFFts.fill(mdi)
DIFFtsMSK.fill(mdi)

if not(plotdiff):
    ALLtrends=np.zeros((nsources,nREGs,3)) # row for each source, columns for trend, 5thpct, 95thpct
    ALLcorrs=np.zeros((sum(range(nsources)),nREGs)) # row for each source, columns for trend, 5thpct, 95thpct
    ALLtrendsMSK=np.zeros((nsources,nREGs,3)) # row for each source, columns for trend, 5thpct, 95thpct
    ALLcorrsMSK=np.zeros((sum(range(nsources)),nREGs)) # row for each source, columns for trend, 5thpct, 95thpct
else:
    ALLtrends=np.zeros((nsources-1,nREGs,3)) # row for each source, columns for trend, 5thpct, 95thpct
    ALLcorrs=np.zeros((sum(range(nsources-1)),nREGs)) # row for each source, columns for trend, 5thpct, 95thpct
    ALLtrendsMSK=np.zeros((nsources-1,nREGs,3)) # row for each source, columns for trend, 5thpct, 95thpct
    ALLcorrsMSK=np.zeros((sum(range(nsources-1)),nREGs)) # row for each source, columns for trend, 5thpct, 95thpct

ALLtrends.fill(mdi)
ALLcorrs.fill(mdi)
ALLtrendsMSK.fill(mdi)
ALLcorrsMSK.fill(mdi)

# Can now establish the filename
plottee=OUTFIL+'_'+moop
if not(plotdiff):
    plotteeC=OUTFILC+'_'+moop

# Loop through sources
for ns in range(nsources):
    # set up tmparr for gridbox data
    tmp=np.empty((nREGs,nmons))
    tmp.fill(mdi)
    tmpMSK=np.empty((nREGs,nmons))
    tmpMSK.fill(mdi)
	
# Read in netCDF set of grids and extract gridbox
    if (ns > 1):
        TheRegions=RegionsOTH
    else:
        TheRegions=Regions
    tmp=ExtractREG(infilee[ns],TheRegions,tmp,styr,edyr,nREGs,ns)
    bads=np.where(tmp < -999.)[0]
    tmp[bads]=mdi
    tmpMSK=ExtractREG(infilee[ns+6],TheRegions,tmpMSK,styr,edyr,nREGs,ns+6)
    bads=np.where(tmpMSK < -999.)[0]
    tmpMSK[bads]=mdi
        	
# Reanomalies to climatology period
    tmp=CalcAnoms(tmp,stcl,edcl,mdi,nREGs) # should modify tmpGB
    tmpMSK=CalcAnoms(tmpMSK,stcl,edcl,mdi,nREGs) # should modify tmpGB
	
# If timetype is annual then average to annual	
    if (timetype == 'annual'):
    	REGts[ns,:,:]=AverageAnnuals(tmp,mdi,nREGs)
        REGtsMSK[ns,:,:]=AverageAnnuals(tmpMSK,mdi,nREGs)
    else:
    	REGts[ns,:,:]=tmp
        REGtsMSK[ns,:,:]=tmpMSK
	        
    if not(plotdiff):
        # Fit a linear trend if desired
        if (fitlin):
            ALLtrends[ns,:,:]=FitTrend(REGts[ns,:,:],timetype,mdi,nREGs)
            ALLtrendsMSK[ns,:,:]=FitTrend(REGtsMSK[ns,:,:],timetype,mdi,nREGs)
	  
if not(plotdiff):    
    # Get correlations if desired (necesary for plotting correlation matrices)
    if (addcor):
        ALLcorrs=GetCorr(REGts,timetype,mdi,nREGs,nsources)
        ALLcorrsMSK=GetCorr(REGtsMSK,timetype,mdi,nREGs,nsources)
    # Pass all sources and trends and otherinfo to plotter
    print('Plotting...')
    PlotCorrMatrix(plotteeC,ALLcorrs,ALLcorrsMSK,Namey,ShortSourceNames)
    PlotTimeSeries(plottee,REGts,ALLtrends,ALLcorrs,REGtsMSK,ALLtrendsMSK,ALLcorrsMSK,OtherInfo,unitees,nmons,nyrs,timetype,styr,edyr,mdi,nsources,Namey,ShortSourceNames,nREGs)

else:
    # Now sort out difference arrays, trends and source names
    nsources=nsources-1
    NewShortSourceNames=ShortSourceNames[1:nsources+1]
    for ns in range(nsources):
        for nr in range(nREGs):
	    gots=np.where((REGts[0,nr,:] != mdi) & (REGts[ns+1,nr,:] != mdi))[0]
	    DIFFts[ns,nr,gots]=REGts[ns+1,nr,:]-REGts[0,nr,gots]
	    gots=np.where((REGtsMSK[0,nr,:] != mdi) & (REGtsMSK[ns+1,nr,:] != mdi))[0]
	    DIFFtsMSK[ns,nr,gots]=REGtsMSK[ns+1,nr,:]-REGtsMSK[0,nr,gots]
        # Fit a linear trend if desired
        if (fitlin):
            ALLtrends[ns,:,:]=FitTrend(DIFFts[ns,:,:],timetype,mdi,nREGs)
            ALLtrendsMSK[ns,:,:]=FitTrend(DIFFtsMSK[ns,:,:],timetype,mdi,nREGs)

    PlotTimeSeries(plottee,DIFFts,ALLtrends,ALLcorrs,DIFFtsMSK,ALLtrendsMSK,ALLcorrsMSK,OtherInfo,unitees,nmons,nyrs,timetype,styr,edyr,mdi,nsources,Namey,NewShortSourceNames,nREGs)
    
stop()

print("And, we are done!")

