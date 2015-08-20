#!/usr/local/sci/bin/python

#***************************************
# 22 March 2014 KMW - v1
# Plots station locations for each variable including removed sub/super sats  
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotTemporalCoverage_MAR2014.py
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
import os.path
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from netCDF4 import Dataset
from scipy.io import netcdf


# RESTART VALUE
Restarter='------'				#'------'		#'681040'

# Set up initial run choices
param='e'	# tw, q, e, rh, t, td, dpd
param2='e'	# Tw, q, e, RH, T, Td, DPD
homogtype='IDPHA'	# 'IDPHA','PHA','PHADPD'
nowmon='SEP'
nowyear='2014'
thenmon='JAN'
thenyear='2014'
version='2.0.0.2013p'
styr=1973
edyr=2013
nyrs=edyr-styr
nmons=(nyrs+1)*12

# Set up directories and files
if param=='q': 
    STATIONDIR='IDPHANETCDF/QDIR/'
elif param=='e': 
    STATIONDIR='IDPHANETCDF/EDIR/'
elif param=='rh': 
    STATIONDIR='IDPHANETCDF/RHDIR/'
elif param=='t': 
    STATIONDIR='IDPHANETCDF/TDIR/'
elif param=='tw': 
    STATIONDIR='IDPHANETCDF/TWDIR/'
elif param=='td': 
    STATIONDIR='IDPHANETCDF/TDDIR/'
elif param=='dpd': 
    STATIONDIR='PHANETCDF/DPDDIR/'

LISTDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/LISTS_DOCS/'
PLOTDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/'
DATADIR='/data/local/hadkw/HADCRUH2/UPDATE2013/MONTHLIES/HOMOG/'+STATIONDIR

INLIST='Posthomog'+homogtype+param+'_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INSATS='Posthomog'+homogtype+param+'_satsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INSUBS='Posthomog'+homogtype+param+'_subzerosHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
OUTPLOT='TemporalCoverageWMO.land'+param2+'.'+version+'_'+nowmon+nowyear
OUTLIST='TemporalCoverageWMO.land'+param2+'.'+version+'_'+nowmon+nowyear+'.txt'

DATAFILE='_homogJAN2014.nc'

# Set up variables
ngoods=0	#set once file read in
nsats=0		#set once file read in
nsubs=0		#set once file read in

GoodWMOs=[]
GoodWBANs=[]
SatWMOs=[]
SubWMOs=[]

GoodGots=[]	# will be an nStation(row) by nMonths(column) array of 0/1

TotalGots=np.zeros(nmons)

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
# SearchStation
def SearchStation(TheGots,TheWMOs,TheWBANs,TheFilePre,TheFilePost,TheParam,TheNCounts):
    ''' Loop through all stations. '''
    ''' Open station netCDF and fill row with 0s/1s for '''
    ''' month presences '''
    for ss in range(TheNCounts):
        MyNCFile=TheFilePre+TheWMOs[ss]+TheWBANs[ss]+TheFilePost
        f=netcdf.netcdf_file(MyNCFile,'r')
        if TheParam=='q': 
 	    var=f.variables['q_anoms']
        elif TheParam=='e':
            var=f.variables['e_anoms']
        elif TheParam=='rh':
	    var=f.variables['rh_anoms']
        elif TheParam=='t':
	    var=f.variables['t_anoms']
        elif TheParam=='tw':
	    var=f.variables['tw_anoms']
        elif TheParam=='td':
	    var=f.variables['td_anoms']
        elif TheParam=='dpd':
	    var=f.variables['dpd_anoms']	    
	vals=np.array(var.data)
        f.close()
        
#	print(vals[vals > 0.])
        TheGots[ss,vals > (-100.)]=1
#        print(TheWMOs[ss],vals[0:20],TheGots[ss,0:20])
#        stop()
   
    return TheGots # SearchStation

#************************************************************************
# WriteCoverage
def WriteCoverage(TheFile,TheWMOs,TheGots,TheNCount,TheMCount):
    ''' Write a text file with a row for each station '''
    ''' and a 0 or 1 in each column for whether that '''
    ''' monthy is present. '''

    for outt in range(TheNCount):
        if outt == 0:
	    goo=np.append(TheWMOs[outt],TheGots[outt,:])
	else:
	    goo=np.vstack((goo,np.append(TheWMOs[outt],TheGots[outt,:])))

    np.savetxt(TheFile,goo,fmt=list(np.append('%6s',['%3s']*TheMCount)),delimiter=' ')
    
    return # WriteCoverage

#************************************************************************
# PlotNiceCoverageGraph
def PlotNiceCoverageGraph(TheFile,TheWMOs,TheGots,TheTotalGots,
                          TheNCount,TheMCount,TheStYr,TheEdYr,
			  TheParam,TheLetter,TheNamey):
    ''' Plot lines for all months with data for all stations by WMO ID - left Yaxis'''
    ''' Colour by WMO ID block '''
    ''' Plot line for total month ob count - right Y axis '''
    ''' Save as png and eps '''
    
    # WMO Blocks
    # 0-199999 = Europe
    # 200000-379999 = Russia and Eastern Europe
    # 380000-439999 = Central Asia, Middle East, India/Pakistan
    # 440000-599850 = East Asia, China and Taiwan
    # 600000-689999 = Africa
    # 690000-749999 = USA, Canada
    # 760000-819999 = Central America
    # 820000-879999 = South America
    # 880000-899999 = Antarctica
    # 911000-919999 = Pacific Islands (inc. Hawaii)
    # 930000-949999 = Australasia (in. NZ)
    # 960000-988999 = Indonesia/Phillippines/Borneo
    
    collist=list(['DarkRed','Crimson','OrangeRed','DarkOrange','Gold',
             'DarkGreen','OliveDrab','MediumBlue','DeepSkyBlue',
	     'MidnightBlue','MediumSlateBlue','DarkViolet'])
    labblist=list([['Europe',0,199999],
              ['Russia/Eastern Europe',200000,379999],
	      ['Central and southern Asia/Middle East',380000,439999],
	      ['East Asia',440000,599999],
	      ['Africa',600000,689999],
	      ['North America',690000,749999],
	      ['Central America',760000,799999],
	      ['South America',800000,879999],
	      ['Antarctica',880000,899999],
	      ['Pacific Islands',910000,919999],
	      ['Australasia',930000,949999],
	      ['Indonesia/Philippines/Borneo',960000,999999]])
	  
    TotsArr=np.zeros((12,TheMCount),dtype=np.int) # 2d array for WMO month counts
    
    colpoints=[]
    Beginner=0	# this is a pointer to tell us when to make a Y axis label, store the location too
    BeginnerALT=10000	# this is a pointer to tell us when to make a Y axis label, store the location too
    FlipSwitch=100000
    YaxPointers=[]
    YaxPointers.append(0)
    YaxPointersALT=[]
    YaxPointersALT.append(0)
    
    for n in range(TheNCount):
        if TheWMOs[n].astype(int) <= labblist[0][2] and TheWMOs[n].astype(int) >= labblist[0][1]:
	    TotsArr[0,]=TotsArr[0,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[1][2] and TheWMOs[n].astype(int) >= labblist[1][1]: 
	    TotsArr[1,]=TotsArr[1,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[2][2] and TheWMOs[n].astype(int) >= labblist[2][1]: 
	    TotsArr[2,]=TotsArr[2,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[3][2] and TheWMOs[n].astype(int) >= labblist[3][1]: 
	    TotsArr[3,]=TotsArr[3,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[4][2] and TheWMOs[n].astype(int) >= labblist[4][1]: 
	    TotsArr[4,]=TotsArr[4,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[5][2] and TheWMOs[n].astype(int) >= labblist[5][1]: 
	    TotsArr[5,]=TotsArr[5,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[6][2] and TheWMOs[n].astype(int) >= labblist[6][1]: 
	    TotsArr[6,]=TotsArr[6,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[7][2] and TheWMOs[n].astype(int) >= labblist[7][1]: 
	    TotsArr[7,]=TotsArr[7,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[8][2] and TheWMOs[n].astype(int) >= labblist[8][1]: 
	    TotsArr[8,]=TotsArr[8,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[9][2] and TheWMOs[n].astype(int) >= labblist[9][1]: 
	    TotsArr[9,]=TotsArr[9,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[10][2] and TheWMOs[n].astype(int) >= labblist[10][1]: 
	    TotsArr[10,]=TotsArr[10,]+TheGots[n] 
        if TheWMOs[n].astype(int) <= labblist[11][2] and TheWMOs[n].astype(int) >= labblist[11][1]: 
	    TotsArr[11,]=TotsArr[11,]+TheGots[n] 

    print(TheNCount,len(colpoints))
    
    # set up axes
    TheYCount=TheMCount/12
    TheYears=np.reshape(range(TheStYr,TheStYr+TheYCount),TheYCount)
    #TheMonths=np.array(range(0,TheMCount))
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
#    print(TheMonths[0:10])
#    print(type(TheMonths))
    
    ytitlee1='No. of Stations in WMO Block'
    ytitlee2='Total No. of Stations'
    xtitlee='Years'
#    titlee=TheParam+' Observation Counts'
    titlee=TheNamey
            
    # set up plot
 
    fig = plt.figure(1,figsize=(8,6))
    plt.clf()	# needs to be called after figure!!! (create the figure, then clear the plot space)

    ax1=plt.subplot(1,1,1)

    ax1.set_xlim([TheMonths[0],TheMonths[TheMCount-1]])   

    ax1.set_xticks([TheMonths[2*12],TheMonths[7*12],TheMonths[12*12],
                TheMonths[17*12],TheMonths[22*12],TheMonths[27*12],
		TheMonths[32*12],TheMonths[37*12]])
    ax1.set_xticklabels([TheMonths[2*12].year,TheMonths[7*12].year,TheMonths[12*12].year,
                TheMonths[17*12].year,TheMonths[22*12].year,TheMonths[27*12].year,
		TheMonths[32*12].year,TheMonths[37*12].year],fontsize=17)

    ax1.tick_params(axis='both', which='major', labelsize=17)

    for wmo in range(12):
        ax1.plot(TheMonths,TotsArr[wmo,],color=collist[wmo],linewidth=4)
    
    ax1.set_title(titlee,size=18)
    
#    ax1.set_ylim([0,TheNCount])		#ax1.set_ylim([0,1000000])    
#    ax1.set_xlim([TheMonths[0],TheMonths[TheMCount-1]])   

    ax1.set_xlabel(xtitlee,size=18)
    ax1.set_ylabel(ytitlee1,size=18)
    
    ax2=ax1.twinx()
    ax2.set_ylim([0,4000])
    ax2.set_ylabel(ytitlee2,size=18)
    ax2.tick_params(axis='both', which='major', labelsize=17)
    ax2.plot(TheMonths,TheTotalGots,'black',linewidth=4)

    # add watermark and plot labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)

    plt.figtext(0.02,0.96,TheLetter,size=18)
    #plt.figtext(0.5,0.94,TheNamey,size=18,ha='center')
    
    plt.tight_layout()

    #plt.show()
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")

#    raw_input("stop")	# REALLY USEFUL TO INTERACT WITHIN SUBROUTINE ctrl C
    # plt.ion()
    # plt.show() can then zoom and save
     
    return #PlotNiceDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list
MyTypes=("|S6","|S5","float","float","float","|S4","|S30","|S7")
MyDelimiters=[6,5,9,10,7,4,30,7]
MyFile=LISTDIR+INLIST
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
GoodWMOs=np.array(RawData['f0'])
GoodWBANs=np.array(RawData['f1'])
ngoods=len(GoodWMOs)

print("Starting Goods: ",ngoods)

# get longs and lats of goods, supersats and subzeros if they exist
if param != 't':
    print("Looking for sats...")
    MyTypes=("|S6","|S5","|S21")
    MyDelimiters=[6,5,21]
    MyFile=LISTDIR+INSATS
    RawData=ReadData(MyFile,MyTypes,MyDelimiters)
    SatWMOs=np.array(RawData['f0'])
    nsats=len(SatWMOs)
    for id in range(nsats):
        for findee in range(ngoods):
	    if SatWMOs[id] == GoodWMOs[findee]:
	        ngoods=ngoods-1
		GoodWMOs=np.delete(GoodWMOs,findee,axis=0)
		GoodWBANs=np.delete(GoodWBANs,findee,axis=0)
		break

if param == 'q' or param == 'e':
    print("Looking for subs...")
    MyTypes=("|S6","|S5","|S21")
    MyDelimiters=[6,5,21]
    MyFile=LISTDIR+INSUBS
    RawData=ReadData(MyFile,MyTypes,MyDelimiters)
    SubWMOs=np.array(RawData['f0'])
    nsubs=len(SubWMOs)
    for id in range(nsubs):
        for findee in range(ngoods):
	    if SubWMOs[id] == GoodWMOs[findee]:
	        ngoods=ngoods-1
		GoodWMOs=np.delete(GoodWMOs,findee,axis=0)
		GoodWBANs=np.delete(GoodWBANs,findee,axis=0)
		break

print("Ending Goods: ",ngoods,nsats,nsubs)


# Loop through Good Stations and search for data presence
GoodGots=np.zeros((ngoods,nmons))
# if list file exists then read this in (quicker) else loop through data
MyListFile=LISTDIR+OUTLIST
if os.path.isfile(MyListFile):
    MyTypes=np.append("|S16",["int"]*nmons)
    MyDelimiters=np.append([6],[4]*nmons)
    RawData=ReadData(MyListFile,MyTypes,MyDelimiters)
    
    for n in range(ngoods):
        var=tuple(RawData[n])
        GoodGots[n,:]=var[1:nmons+1]
	
#    print(GoodGots[0,:])
else:
    MyFilePre=DATADIR
    MyFilePost=DATAFILE
    GoodGots=SearchStation(GoodGots,GoodWMOs,GoodWBANs,MyFilePre,MyFilePost,param,ngoods)
    
    GoodGots=GoodGots.astype(int)

    # Write this all out for future use
    WriteCoverage(MyListFile,GoodWMOs,GoodGots,ngoods,nmons)    

# get total obs per month
TotalGots[:]=np.sum(GoodGots,axis=0)	#sum all rows in each column

# pass to plotter
if param == 't':
    Letty='a)'
    Namey='temperature'
if param == 'tw':
    Letty='e)'
    Namey='wet bulb temperature'
if param == 'td':
    Letty='b)'
    Namey='dew point temperature'
if param == 'q':
    Letty='c)'
    Namey='specific humidity'
if param == 'e':
    Letty='f)'
    Namey='vapour pressure'
if param == 'rh':
    Letty='d)'
    Namey='relative humidity'
if param == 'dpd':
    Letty='g)'
    Namey='dew point depression'

MyFile=PLOTDIR+OUTPLOT
PlotNiceCoverageGraph(MyFile,GoodWMOs,GoodGots,TotalGots,
                      ngoods,nmons,styr,edyr,param2,
		      Letty,Namey)
		
#    stop()

print("And, we are done!")

