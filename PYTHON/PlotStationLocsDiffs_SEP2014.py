#!/usr/local/sci/bin/python

#***************************************
# 22 March 2014 KMW - v1
# Plots station locations common to all variables and then maps 
# for each variable including removed sub/super sats  
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotStationLocsDiffs_SEP2014.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import scipy.stats
import os.path
import struct
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from netCDF4 import Dataset


# RESTART VALUE
Restarter='------'				#'------'		#'681040'

# Set up initial run choices
param=['t','td','q','rh','tw','e','dpd']
param2=['T','Td','q','RH','Tw','e','DPD']
homogtype=['IDPHA','PHADPD','IDPHA','IDPHA','IDPHA','IDPHA','PHA']	# 'IDPHA','PHA','PHADPD'
nowmon='SEP'
nowyear='2014'
thenmon='JAN'
thenyear='2014'
version='2.0.0.2013p'

# Set up directories and files
INDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/LISTS_DOCS/'
OUTDIR='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/'

OUTPLOTALL='StationLocationsDiffs.landALL.'+version+'_'+nowmon+nowyear
ALLINS=['PosthomogIDPHAt_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt',
	'PosthomogPHADPDtd_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt',
	'PosthomogIDPHAq_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt',
	'PosthomogIDPHArh_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt',
	'PosthomogIDPHAtw_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt',
	'PosthomogIDPHAe_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt',
	'PosthomogPHAdpd_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt']

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

allgoods=np.zeros(7,dtype='i')
allsats=np.zeros(7,dtype='i')	
allsubs=np.zeros(7,dtype='i')
	
AllGoodDiffsWMOs=[[] for i in xrange(7)]	# list of lists for all variables
AllGoodDiffsLats=[[] for i in xrange(7)]
AllGoodDiffsLons=[[] for i in xrange(7)]
AllSatDiffsWMOs=[[] for i in xrange(7)]
AllSatDiffsLats=[[] for i in xrange(7)]
AllSatDiffsLons=[[] for i in xrange(7)]
AllSubDiffsWMOs=[[] for i in xrange(7)]
AllSubDiffsLats=[[] for i in xrange(7)]
AllSubDiffsLons=[[] for i in xrange(7)]

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
# PlotNiceDiffDotsAllMap
def PlotNiceDiffDotsAllMap(TheFile,TheGLats,TheGLons,TheSTLats,TheSTLons,TheSBLats,TheSBLons,
                    TheAllGLats,TheAllGLons,TheAllSTLats,TheAllSTLons,TheAllSBLats,TheAllSBLons,
                    TheGCount,TheSTCount,TheSBCount,TheAllGCount,TheAllSTCount,TheAllSBCount,TheLetter,TheNamee):
    ''' Plot all good lats and lons TheGLats, TheGLons '''
    ''' If there are sats (TheSTCount > 0) plot all sat lats TheSTLats, TheSTLons '''
    ''' If there are subzeros (TheSBCount > 0) plot all sat lats TheSBLats, TheSBLons '''
    ''' Label plot with totals '''
    ''' Save as png and eps '''
      
    # set up plot
    xpos=[0.05,0.01,0.51,0.01,0.51]
    ypos=[0.63,0.33,0.33,0.03,0.03]
    xfat=[0.9,0.48,0.48,0.48,0.48]
    ytall=[0.32,0.26,0.26,0.26,0.26]
    
    #f,axarr=plt.subplots(5,figsize=(12,14),sharex=False)	#6,18
    f=plt.figure(5,figsize=(12,13))	#6,18

    #axarr[0].set_position([xpos[0],ypos[0],xfat[0],ytall[0]])
    plt.axes([xpos[0],ypos[0],xfat[0],ytall[0]])
 
#    plt.clf()
#    plt.figure(1,figsize=(8,3))
#    plt.axes([0.05,0.05,0.9,0.95])
    
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
    m.scatter(x,y,s=1,marker='o',color="LightSkyBlue")
    
    plt.figtext(xpos[0]+0.15,ypos[0]+ytall[0],TheLetter[0],size=18)
    plt.figtext(0.5,ypos[0]+ytall[0]+0.01,TheNamee[0],size=18,ha='center')
    
    print(TheGCount,TheSTCount,TheSBCount)
#    totalgoods=TheGCount-(TheSTCount+TheSBCount)

    totalgoods=TheGCount
#    axarr[0].annotate(str(totalgoods)+" good stations",xy=(0.52,0.265),xytext=None, xycoords='axes fraction',color="black")
#    plt.annotate(str(totalgoods)+" good stations",xy=(0.52,0.265),xytext=None, xycoords='axes fraction',color="black")
    x,y = m(5,-40)
    m.scatter(x,y,s=15,marker='o',color='LightSkyBlue',)
    x,y = m(10,-40)
    plt.text(x,y,str(totalgoods)+" good",ha="left", va="center", size=14)    

    for plts in range(4):
#        axarr[plts+1].set_position([xpos[plts+1],ypos[plts+1],xfat[plts+1],ytall[plts+1]])
        plt.axes([xpos[plts+1],ypos[plts+1],xfat[plts+1],ytall[plts+1]])
 
    
        # plot map without continents and coastlines
        m = Basemap(projection='kav7',lon_0=0)
        # draw map boundary, transparent
        m.drawmapboundary()
        m.drawcoastlines()
        # draw paralells and medians, no labels
        m.drawparallels(np.arange(-90,90.,30.))
        m.drawmeridians(np.arange(-180,180.,60.))
    
        # plot black dots for the goods
        x,y = m(TheAllGLons[plts],TheAllGLats[plts])
        m.scatter(x,y,s=2,marker='o',color="CornflowerBlue")
    
        # plot blue dots for the supersats
        if TheAllSTCount[plts] > 0:
            x,y = m(TheAllSTLons[plts],TheAllSTLats[plts])
            m.scatter(x,y,s=5,marker='o',color="Violet")

        # plot red dots for the subs
        if TheAllSBCount[plts] > 0:
            x,y = m(TheAllSBLons[plts],TheAllSBLats[plts])
            m.scatter(x,y,s=5,marker='o',color="Gold")
    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
        plt.figtext(xpos[plts+1]+0.05,ypos[plts+1]+ytall[plts+1],TheLetter[plts+1],size=18)
        plt.figtext(xpos[plts+1]+0.24,ypos[plts+1]+ytall[plts+1]+0.01,TheNamee[plts+1],size=18,ha='center')
    
        print(TheAllGCount[plts],TheAllSTCount[plts],TheAllSBCount[plts])
        totalgoods=TheAllGCount[plts]-(TheAllSTCount[plts]+TheAllSBCount[plts])   
#        plt.annotate(str(totalgoods)+" good stations",xy=(0.52,0.265),xytext=None, xycoords='axes fraction',color="black")
        x,y = m(5,-40)
        m.scatter(x,y,s=15,marker='o',color='CornflowerBlue',)
        x,y = m(10,-40)
        plt.text(x,y,str(totalgoods)+" good",ha="left", va="center", size=12)    
    
        if TheAllSTCount[plts] > 0:
#            plt.annotate(str(TheAllSTCount[plts])+" supersaturated stations",xy=(0.52,0.22),xytext=None, xycoords='axes fraction',color="black")
            x,y = m(5,-48)
            m.scatter(x,y,s=15,marker='o',color='Violet',)
            x,y = m(10,-48)
            plt.text(x,y,str(TheAllSTCount[plts])+" supersaturated",ha="left", va="center", size=12)    

        if TheAllSBCount[plts] > 0:
#            plt.annotate(str(TheAllSBCount[plts])+" subzero stations",xy=(0.52,0.175),xytext=None, xycoords='axes fraction',color="black")
            x,y = m(5,-56)
            m.scatter(x,y,s=15,marker='o',color='Gold',)
            x,y = m(10,-56)
            plt.text(x,y,str(TheAllSBCount[plts])+" subzero",ha="left", va="center", size=12)    

#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotNiceDiffDotsAllMap
    
#************************************************************************
# PlotNiceDiffDotsMap
def PlotNiceDiffDotsMap(TheFile,TheGLats,TheGLons,TheSTLats,TheSTLons,TheSBLats,TheSBLons,
                    TheGCount,TheSTCount,TheSBCount,TheLetter,TheNamee):
    ''' Plot all good lats and lons TheGLats, TheGLons '''
    ''' If there are sats (TheSTCount > 0) plot all sat lats TheSTLats, TheSTLons '''
    ''' If there are subzeros (TheSBCount > 0) plot all sat lats TheSBLats, TheSBLons '''
    ''' Label plot with totals '''
    ''' Save as png and eps '''
      
    # set up plot
 
    plt.clf()
    plt.figure(1,figsize=(8,5))
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
    m.scatter(x,y,s=2,marker='o',color="CornflowerBlue")
    
    # plot blue dots for the supersats
    if TheSTCount > 0:
        x,y = m(TheSTLons,TheSTLats)
        m.scatter(x,y,s=5,marker='o',color="Violet")

    # plot red dots for the subs
    if TheSBCount > 0:
        x,y = m(TheSBLons,TheSBLats)
        m.scatter(x,y,s=5,marker='o',color="Gold")
    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.92,TheLetter,size=18)
    plt.figtext(0.5,0.95,TheNamee,size=18,ha='center')
    
    print(TheGCount,TheSTCount,TheSBCount)
    totalgoods=TheGCount-(TheSTCount+TheSBCount)
#    plt.annotate(str(totalgoods)+" good stations",xy=(0.52,0.265),xytext=None, xycoords='axes fraction',color="black")
    x,y = m(5,-40)
    m.scatter(x,y,s=15,marker='o',color='CornflowerBlue',)
    x,y = m(10,-40)
    plt.text(x,y,str(totalgoods)+" good",ha="left", va="center", size=14)    
    
    if TheSTCount > 0:
#        plt.annotate(str(TheSTCount)+" supersaturated error stations",xy=(0.52,0.22),xytext=None, xycoords='axes fraction',color="black")
        x,y = m(5,-48)
        m.scatter(x,y,s=15,marker='o',color='Violet',)
        x,y = m(10,-48)
        plt.text(x,y,str(TheSTCount)+" supersaturated",ha="left", va="center", size=14)    

    if TheSBCount > 0:
#        plt.annotate(str(TheSBCount)+" subzero error stations",xy=(0.52,0.175),xytext=None, xycoords='axes fraction',color="black")
        x,y = m(5,-56)
        m.scatter(x,y,s=15,marker='o',color='Gold',)
        x,y = m(10,-56)
        plt.text(x,y,str(TheSBCount)+" subzero",ha="left", va="center", size=14)    

#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotNiceDiffDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in all of the variable lists to fill arrays
for vari in range(7):
    INLIST='Posthomog'+homogtype[vari]+param[vari]+'_goodsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
    INSATS='Posthomog'+homogtype[vari]+param[vari]+'_satsHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
    INSUBS='Posthomog'+homogtype[vari]+param[vari]+'_subzerosHadISDH.'+version+'_'+thenmon+thenyear+'.txt'
 
    # read in station list
    MyTypes=("|S6","|S5","float","float","float","|S4","|S30","|S7")
    MyDelimiters=[6,5,9,10,7,4,30,7]
    MyFile=INDIR+INLIST
    RawData=ReadData(MyFile,MyTypes,MyDelimiters)
    AllGoodDiffsWMOs[vari]=np.array(RawData['f0'])
    AllGoodDiffsLats[vari]=np.array(RawData['f2'])
    AllGoodDiffsLons[vari]=np.array(RawData['f3'])
    allgoods[vari]=len(AllGoodDiffsWMOs[vari])

    # get longs and lats of goods, supersats and subzeros if they exist
    if param[vari] != 't':
        print("Looking for sats...")
        MyTypes=("|S6","|S5","|S21")
        MyDelimiters=[6,5,21]
        MyFile=INDIR+INSATS
        RawData=ReadData(MyFile,MyTypes,MyDelimiters)
        AllSatDiffsWMOs[vari]=np.array(RawData['f0'])
        allsats[vari]=len(AllSatDiffsWMOs[vari])
        AllSatDiffsLats[vari]=np.zeros(allsats[vari])
        AllSatDiffsLons[vari]=np.zeros(allsats[vari])
        for id in range(allsats[vari]):
            for findee in range(allgoods[vari]):
	        if AllSatDiffsWMOs[vari][id] == AllGoodDiffsWMOs[vari][findee]:
	            AllSatDiffsLats[vari][id]=AllGoodDiffsLats[vari][findee]
		    AllSatDiffsLons[vari][id]=AllGoodDiffsLons[vari][findee]
		    break
        # Now remove those Supersat stations from the Goods list
	getrids=[indi for indi,sts in enumerate(AllGoodDiffsWMOs[vari]) if np.array(np.where(np.array(AllSatDiffsWMOs[vari]) == sts)).size > 0]
#    stop()
        AllGoodDiffsWMOs[vari]=[sts for indi,sts in enumerate(AllGoodDiffsWMOs[vari]) if indi not in getrids]
        AllGoodDiffsLats[vari]=[sts for indi,sts in enumerate(AllGoodDiffsLats[vari]) if indi not in getrids]
        AllGoodDiffsLons[vari]=[sts for indi,sts in enumerate(AllGoodDiffsLons[vari]) if indi not in getrids]

    if param[vari] == 'q' or param[vari] == 'e':
        print("Looking for subs...")
        MyTypes=("|S6","|S5","|S21")
        MyDelimiters=[6,5,21]
        MyFile=INDIR+INSUBS
        RawData=ReadData(MyFile,MyTypes,MyDelimiters)
        AllSubDiffsWMOs[vari]=np.array(RawData['f0'])
        allsubs[vari]=len(AllSubDiffsWMOs[vari])
        AllSubDiffsLats[vari]=np.zeros(allsubs[vari])
        AllSubDiffsLons[vari]=np.zeros(allsubs[vari])
        for id in range(allsubs[vari]):
            for findee in range(allgoods[vari]):
	        if AllSubDiffsWMOs[vari][id] == AllGoodDiffsWMOs[vari][findee]:
	            AllSubDiffsLats[vari][id]=AllGoodDiffsLats[vari][findee]
		    AllSubDiffsLons[vari][id]=AllGoodDiffsLons[vari][findee]
		    break
        # Now remove those Subzero stations from the Goods list
	getrids=[indi for indi,sts in enumerate(AllGoodDiffsWMOs[vari]) if np.array(np.where(np.array(AllSubDiffsWMOs[vari]) == sts)).size > 0]
#    stop()
        AllGoodDiffsWMOs[vari]=[sts for indi,sts in enumerate(AllGoodDiffsWMOs[vari]) if indi not in getrids]
        AllGoodDiffsLats[vari]=[sts for indi,sts in enumerate(AllGoodDiffsLats[vari]) if indi not in getrids]
        AllGoodDiffsLons[vari]=[sts for indi,sts in enumerate(AllGoodDiffsLons[vari]) if indi not in getrids]

# Now find which stations are common to all and which are unique
# Start will smallest number of station (Tw) to get commons
for sts,statid in enumerate(AllGoodDiffsWMOs[4]):
    pointer=np.zeros(7)
    pointer[:]=-99
    pointer[4]=sts	# set the Tw station ID loc to non-missing
    for vv,vari in enumerate([0,1,2,3,5,6]):	# missing out Tw
        gots=np.where(np.array(AllGoodDiffsWMOs[vari]) == statid)
	if np.array(gots).size > 0:
	    pointer[vari]=gots[0] # set the variable station ID loc
	else:
	    break
    if pointer[6] >= 0:	# then we have continued the loop through all variable and have a common station
        # now add to commons list and remove from all others 	
        GoodWMOs.append(statid)
        GoodLats.append(AllGoodDiffsLats[4][sts])
        GoodLons.append(AllGoodDiffsLons[4][sts])

ngoods=len(GoodWMOs)

# Now map to all variables and remove
for vari in range(7):
    getrids=[indi for indi,sts in enumerate(AllGoodDiffsWMOs[vari]) if np.array(np.where(np.array(GoodWMOs) == sts)).size > 0]
#    stop()
    AllGoodDiffsWMOs[vari]=[sts for indi,sts in enumerate(AllGoodDiffsWMOs[vari]) if indi not in getrids]
    AllGoodDiffsLats[vari]=[sts for indi,sts in enumerate(AllGoodDiffsLats[vari]) if indi not in getrids]
    AllGoodDiffsLons[vari]=[sts for indi,sts in enumerate(AllGoodDiffsLons[vari]) if indi not in getrids]

# pass to plotter
Letty=['a)','b)','c)','d)','e)']
Namey=['all variables','temperature','dew point temperature','specific humdity','relative humidity']

MyFile=OUTDIR+OUTPLOTALL
PlotNiceDiffDotsAllMap(MyFile,GoodLats,GoodLons,SatLats,SatLons,SubLats,SubLons,
                       AllGoodDiffsLats,AllGoodDiffsLons,
                        AllSatDiffsLats,AllSatDiffsLons,
			AllSubDiffsLats,AllSubDiffsLons,
			ngoods,nsats,nsubs,allgoods,allsats,allsubs,
			Letty,Namey)

for vv,vari in enumerate([4,5,6]):

    if param[vari] == 'tw':
        Letty='a)'
        Namey='wet bulb temperature'
    if param[vari] == 'e':
        Letty='b)'
        Namey='vapour pressure'
    if param[vari] == 'dpd':
        Letty='c)'
        Namey='dew point depression'
    
    OUTPLOTEXT='StationLocationsDiffs.land'+param2[vari]+'.'+version+'_'+nowmon+nowyear
    MyFile=OUTDIR+OUTPLOTEXT
    PlotNiceDiffDotsMap(MyFile,AllGoodDiffsLats[vari],AllGoodDiffsLons[vari],
                        AllSatDiffsLats[vari],AllSatDiffsLons[vari],
			AllSubDiffsLats[vari],AllSubDiffsLons[vari],
			allgoods[vari],allsats[vari],allsubs[vari],
			Letty,Namey)
		
#    stop()

print("And, we are done!")

