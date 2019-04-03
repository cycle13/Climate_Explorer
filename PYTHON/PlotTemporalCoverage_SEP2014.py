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
param='q'	# tw, q, e, rh, t, td, dpd
param2='q'	# Tw, q, e, RH, T, Td, DPD
homogtype='RAW'	# 'IDPHA','PHA','PHADPD', 'RAW'
nowmon='MAY'
nowyear='2018'
thenmon='JAN'
thenyear='2018'
version='4.0.0.2017f'
styr=1973
edyr=2017
nyrs=edyr-styr
nmons=(nyrs+1)*12

# Set up directories and files
if homogtype == 'RAW':
    STATIONDIR = 'NETCDF/'
else:
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

LISTDIR='/data/local/hadkw/HADCRUH2/UPDATE2017/LISTS_DOCS/'
PLOTDIR='/data/local/hadkw/HADCRUH2/UPDATE2017/IMAGES/BUILD/'
DATADIR='/data/local/hadkw/HADCRUH2/UPDATE2017/MONTHLIES/HOMOG/'+STATIONDIR

INLIST='Posthomog'+homogtype+param+'_anoms8110_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INSATS='Posthomog'+homogtype+param+'_anoms8110_satsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
INSUBS='Posthomog'+homogtype+param+'_anoms8110_subzerosHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
OUTPLOT='TemporalCoverage.land'+param2+'.'+version+'_'+nowmon+nowyear
OUTLIST='TemporalCoverage.land'+param2+'.'+version+'_'+nowmon+nowyear+'.txt'

DATAFILE='_anoms8110_homogJAN2018.nc'

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
	      ['Central America',750000,799999],
#	      ['Central America',760000,799999],
	      ['South America',800000,879999],
	      ['Antarctica',880000,899999],
	      ['Pacific Islands',900000,919999],
#	      ['Pacific Islands',910000,919999],
	      ['Australasia',920000,949999],
#	      ['Australasia',930000,949999],
	      ['Indonesia/Philippines/Borneo',950000,999999]])
#	      ['Indonesia/Philippines/Borneo',960000,999999]])
    
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
	    colpoints.append(collist[0])
        if TheWMOs[n].astype(int) <= labblist[1][2] and TheWMOs[n].astype(int) >= labblist[1][1]: 
	    colpoints.append(collist[1])
	    if Beginner == 0:
	        Beginner = 1	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[2][2] and TheWMOs[n].astype(int) >= labblist[2][1]: 
	    colpoints.append(collist[2])
	    if Beginner == 1:
	        Beginner = 2	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[3][2] and TheWMOs[n].astype(int) >= labblist[3][1]: 
	    colpoints.append(collist[3])
	    if Beginner == 2:
	        Beginner = 3	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[4][2] and TheWMOs[n].astype(int) >= labblist[4][1]: 
	    colpoints.append(collist[4])
	    if Beginner == 3:
	        Beginner = 4	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[5][2] and TheWMOs[n].astype(int) >= labblist[5][1]: 
	    colpoints.append(collist[5])
	    if Beginner == 4:
	        Beginner = 5	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[6][2] and TheWMOs[n].astype(int) >= labblist[6][1]: 
	    colpoints.append(collist[6])
	    if Beginner == 5:
	        Beginner = 6	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[7][2] and TheWMOs[n].astype(int) >= labblist[7][1]: 
	    colpoints.append(collist[7])
	    if Beginner == 6:
	        Beginner = 7	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[8][2] and TheWMOs[n].astype(int) >= labblist[8][1]: 
	    colpoints.append(collist[8])
	    if Beginner == 7:
	        Beginner = 8	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[9][2] and TheWMOs[n].astype(int) >= labblist[9][1]: 
	    colpoints.append(collist[9])
	    if Beginner == 8:
	        Beginner = 9	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[10][2] and TheWMOs[n].astype(int) >= labblist[10][1]: 
	    colpoints.append(collist[10])
	    if Beginner == 9:
	        Beginner = 10	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if TheWMOs[n].astype(int) <= labblist[11][2] and TheWMOs[n].astype(int) >= labblist[11][1]: 
	    colpoints.append(collist[11])
	    if Beginner == 10:
	        Beginner = 11	# now in next boundary so set up Y axis labels and pointers
		YaxPointers.append(n)
        if BeginnerALT != FlipSwitch and TheWMOs[n].astype(int) >= FlipSwitch:
            BeginnerALT=FlipSwitch
	    FlipSwitch=FlipSwitch+100000
            YaxPointersALT.append(n)

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
    
    ytitlee1='WMO ID'
    ytitlee2='No. of Stations'
    xtitlee='Years'
#    titlee=TheParam+' Observation Counts'
    titlee=TheNamey
            
    # set up plot
 
    fig = plt.figure(1,figsize=(8,11))
    plt.clf()	# needs to be called after figure!!! (create the figure, then clear the plot space)
#    ax1=fig.add_subplot(111,title=titlee)#
    ax1=fig.add_subplot(111)#

    for wmo in range(TheNCount):
#        print(TheGots[wmo,:])
#        print(TheMonths[TheGots[wmo,:] > 0])
#        print(TheMonths[TheGots[wmo] > 0])
#	print(np.repeat(TheWMOs[wmo].astype(int),len(TheMonths[TheGots[wmo] > 0])))

        plt.plot(TheMonths[TheGots[wmo] == 0],np.repeat(wmo,len(TheMonths[TheGots[wmo] == 0])),
	         color=colpoints[wmo],marker='.',markersize=0.01, ls='')	#,linewidth=2)	
#        plt.plot(TheMonths[TheGots[wmo] > 0],np.repeat(wmo,len(TheMonths[TheGots[wmo] > 0])),#
#	         color=colpoints[wmo],marker='.',markersize=0.01, ls='')	#,linewidth=2)	
    
    ax1.set_title(titlee,size=18)
    #ax1=plt.axes([0.15,0.1,0.85,0.9])
    ax1.set_ylim([0,TheNCount])		#ax1.set_ylim([0,1000000])    
    #plt.yticks(YaxPointers,['010000','200000','380000','440000','600000','690000','760000','820000','880000','910000','930000','960000'])
#    plt.yticks(np.append(YaxPointers[0:9],YaxPointers[11]),['010000','200000','380000','440000','600000','690000','760000','820000','880000','960000'])
    plt.yticks(YaxPointersALT,['010000','100000','200000','300000','400000','500000','600000','700000','800000','900000'])
    ax1.set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
    plt.xticks([TheMonths[2*12],TheMonths[7*12],TheMonths[12*12],
                TheMonths[17*12],TheMonths[22*12],TheMonths[27*12],
		TheMonths[32*12],TheMonths[37*12],TheMonths[42*12]])
    ax1.tick_params(axis='both', which='major', labelsize=17)
#    plt.xlim([TheStYr,TheStYr+TheYCount])
   

    ax1.set_xlabel(xtitlee,size=18)
    ax1.set_ylabel(ytitlee1,size=18)
    
    ax2=ax1.twinx()
    ax2.set_ylim([0,9000])
    ax2.set_ylabel(ytitlee2,size=18)
    ax2.tick_params(axis='both', which='major', labelsize=17)
    plt.plot(TheMonths,TheTotalGots,'black',linewidth=4)

    #for tl in ax1.get_yticklabels():
    #    tl.set_color('b')

    # add watermark and plot labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.02,0.97,TheLetter,size=18)
    #plt.figtext(0.5,0.94,TheNamey,size=18,ha='center')
    

    # annotate figures with WMO bands

    plt.tight_layout()

    #plt.show()
#    plt.savefig(TheFile+"INV.eps")
#    plt.savefig(TheFile+"INV.png")
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

MyFile=PLOTDIR+OUTPLOT
PlotNiceCoverageGraph(MyFile,GoodWMOs,GoodGots,TotalGots,
                      ngoods,nmons,styr,edyr,param2,
		      Letty,Namey)
		
#    stop()

print("And, we are done!")

