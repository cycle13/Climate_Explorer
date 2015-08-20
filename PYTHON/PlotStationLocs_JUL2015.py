#!/usr/local/sci/bin/python

#***************************************
# 22 March 2014 KMW - v1
# Plots station locations for each variable including removed sub/super sats if desired  
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotStationLocs_JUL2015.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import scipy.stats
import struct
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from netCDF4 import Dataset


# RESTART VALUE
Restarter='------'				#'------'		#'681040'

Plottype='SATSSUBS' #'JUSTGOODS','SATSSUBS'

# Set up initial run choices
param='t'	# tw, q, e, rh, t, td, dpd
param2='T'	# Tw, q, e, RH, T, Td, DPD
homogtype='IDPHA'	# 'IDPHA','PHA','PHADPD'
nowmon='JUL'
nowyear='2015'
thenmon='JAN'
thenyear='2015'
version='2.0.1.2014p'
version2='2.0.1.2014p_SELECTmax1_most3'

# Set up directories and files
INDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'

#INLIST='Posthomog'+homogtype+param+'_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INLIST='HadISDH.land'+param2+'.2.0.1.2014p_SELECTmax1_most3_JUL2015.txt'
INSATS='Posthomog'+homogtype+param+'_satsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INSUBS='Posthomog'+homogtype+param+'_subzerosHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
OUTPLOT='StationLocations.land'+param2+'.'+version2+'_'+nowmon+nowyear

# Set up variables
ngoods=0	#set once file read in
nsats=0		#set once file read in
nsubs=0		#set once file read in

GoodWMOs=[]
GoodLats=[]
GoodLons=[]
SatWMOs=[]
SatLats=[]
SatLons=[]
SubWMOs=[]
SubLats=[]
SubLons=[]

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
# PlotNiceDotsMap
def PlotNiceDotsMap(TheFile,TheGLats,TheGLons,TheSTLats,TheSTLons,TheSBLats,TheSBLons,
                    TheGCount,TheSTCount,TheSBCount,TheLetter,TheNamee,ThePlot):
    ''' Plot all good lats and lons TheGLats, TheGLons '''
    ''' If there are sats (TheSTCount > 0) plot all sat lats TheSTLats, TheSTLons '''
    ''' If there are subzeros (TheSBCount > 0) plot all sat lats TheSBLats, TheSBLons '''
    ''' Label plot with totals '''
    ''' Save as png and eps '''
      
    # set up plot
 
    plt.clf()
    plt.figure(1,figsize=(8,3))
    plt.axes([0.05,0.05,0.9,0.9])
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
    
    # plot black dots for the goods
    x,y = m(TheGLons,TheGLats)
    m.scatter(x,y,s=7,marker='o',color="CornflowerBlue")
    
    # plot blue dots for the supersats
    if TheSTCount > 0:
        x,y = m(TheSTLons,TheSTLats)
        m.scatter(x,y,s=7,marker='o',color="Violet")

    # plot red dots for the subs
    if TheSBCount > 0:
        x,y = m(TheSBLons,TheSBLats)
        m.scatter(x,y,s=7,marker='o',color="Gold")
    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.9,TheLetter,size=16)
    plt.figtext(0.5,0.9,TheNamee,size=16,ha='center')
    
    print(ngoods,nsats,nsubs)
    totalgoods=TheGCount-(TheSTCount+TheSBCount)
    plt.annotate(str(totalgoods)+" good stations",xy=(0.52,0.265),xytext=None, xycoords='axes fraction',color="black")
    x,y = m(5,-40)
    m.scatter(x,y,s=15,marker='o',color='CornflowerBlue',)
    
    if (ThePlot == 'SATSSUBS'):
        plt.annotate(str(TheSTCount)+" supersaturated error stations",xy=(0.52,0.22),xytext=None, xycoords='axes fraction',color="black")
        x,y = m(5,-48)
        m.scatter(x,y,s=15,marker='o',color='Violet',)

        plt.annotate(str(TheSBCount)+" subzero error stations",xy=(0.52,0.175),xytext=None, xycoords='axes fraction',color="black")
        x,y = m(5,-56)
        m.scatter(x,y,s=15,marker='o',color='Gold',)

#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotNiceDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list
#MyTypes=("|S6","|S5","float","float","float","|S4","|S30","|S7")
MyTypes=("|S6","|S5","float","float","float","|S4","|S30")
#MyDelimiters=[6,5,9,10,7,4,30,7]
MyDelimiters=[6,5,9,10,7,4,30]
MyFile=INDIR+INLIST
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
GoodWMOs=np.array(RawData['f0'])
GoodLats=np.array(RawData['f2'])
GoodLons=np.array(RawData['f3'])
ngoods=len(GoodWMOs)

# If Plottype='SATSSUBS' then look at sats and subs too
nsats=0
nsubs=0

if (Plottype == 'SATSSUBS'):
    # get longs and lats of goods, supersats and subzeros if they exist
    if param != 't':
        print("Looking for sats...")
        MyTypes=("|S6","|S5","|S21")
        MyDelimiters=[6,5,21]
        MyFile=INDIR+INSATS
        RawData=ReadData(MyFile,MyTypes,MyDelimiters)
        SatWMOs=np.array(RawData['f0'])
        totsats=len(SatWMOs)
        SatLats=[]
        SatLons=[]
        for id in range(totsats):
            for findee in range(ngoods):
	        if SatWMOs[id] == GoodWMOs[findee]:
	            SatLats.append(GoodLats[findee])
		    SatLons.append(GoodLons[findee])
		    nsats=nsats+1
		    break

    if param == 'q' or param == 'e':
        print("Looking for subs...")
        MyTypes=("|S6","|S5","|S21")
        MyDelimiters=[6,5,21]
        MyFile=INDIR+INSUBS
        RawData=ReadData(MyFile,MyTypes,MyDelimiters)
        SubWMOs=np.array(RawData['f0'])
        totsubs=len(SubWMOs)
        SubLats=[]
        SubLons=[]
        for id in range(totsubs):
            for findee in range(ngoods):
	        if SubWMOs[id] == GoodWMOs[findee]:
	            SubLats.append(GoodLats[findee])
		    SubLons.append(GoodLons[findee])
		    nsubs=nsubs+1
		    break

# pass to plotter

if param == 't':
    Letty='a)'
    Namey='temperature'
if param == 'tw':
    Letty='a)'
    Namey='wet bulb temperature'
if param == 'td':
    Letty='b)'
    Namey='dew point temperature'
if param == 'q':
    Letty='c)'
    Namey='specific humidity'
if param == 'e':
    Letty='b)'
    Namey='vapour pressure'
if param == 'rh':
    Letty='d)'
    Namey='relative humidity'
if param == 'dpd':
    Letty='c)'
    Namey='dew point depression'
    
MyFile=OUTDIR+OUTPLOT
PlotNiceDotsMap(MyFile,GoodLats,GoodLons,
                        SatLats,SatLons,
			SubLats,SubLons,
			ngoods,nsats,nsubs,
			Letty,Namey,Plottype)
		
#    stop()

print("And, we are done!")

