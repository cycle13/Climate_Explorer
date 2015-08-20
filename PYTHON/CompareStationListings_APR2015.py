#!/usr/local/sci/bin/python

#***************************************
# 28th April 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# Reads in station file for CRUTEM4 (downloaded APR 2015) and for HadISDH.landT
# Tries to find same, unique, maybe stations
#       SAME IF AT LEAST ONE:
# 	Same lat and long (to a single decimal place - may have some 'dupicate' HadISDH stations here then) AND elevation within +/- 100m
#	Same Station name (and lat/lon within +/- 1 whole degrees
# 	MAYBE IF AT LEAST ONE:
#	Nearly Same Station name and lat/long within +/- 1 whole degrees
#	Same lat/lon but elevation diff by > 100m
# Plots station locs coloured by both, crutem only, hadisdh only, queried
# Outputs a list of these station groups

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 CompareStationListings_APR2015.py
#
# REQUIRES
# 
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
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
#from netCDF4 import Dataset
from scipy.io import netcdf
from scipy.stats.stats import pearsonr   
from RandomsRanges import LetterRange
from LinearTrends import MedianPairwise
import re	# regular expression stuff for character replacement within strings
from scipy.stats import itemfreq # a way of looking for unique values and their frequency of occurence

# Input Files
INFILEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/PosthomogIDPHAt_goodsHadISDH.2.0.1.2014p_JAN2015.txt'
#INFILEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CRUTEM4_StationList_APR2015.txt'
INFILEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CRUTEM4_StationList_APR2015.txt'
#INFILEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/ghcnm.tavg.v3.2.2.20150430.qca.inv'

# Output Files
OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CompareStationListings_CRUTEM4HadISDH_APR2015'
OUTLISTSAMEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsSAME_CRUTEM4HadISDH_APR2015.txt'
OUTLISTSAMEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsSAME_HadISDHCRUTEM4_APR2015.txt'
OUTLISTMAYBEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsMAYBE_CRUTEM4HadISDH_APR2015.txt'
OUTLISTMAYBEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsMAYBE_HadISDHCRUTEM4_APR2015.txt'
OUTLISTUNIQC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsUNIQ_CRUTEM4HadISDH_APR2015.txt'
OUTLISTUNIQH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsUNIQ_HadISDHCRUTEM4_APR2015.txt'
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CompareStationListings_GHCNM3CRUTEM4_APR2015'
#OUTLISTSAMEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsSAME_GHCNM3CRUTEM4_APR2015.txt'
#OUTLISTSAMEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsSAME_CRUTEM4GHCNM3_APR2015.txt'
#OUTLISTMAYBEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsMAYBE_GHCNM3CRUTEM4_APR2015.txt'
#OUTLISTMAYBEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsMAYBE_CRUTEM4GHCNM3_APR2015.txt'
#OUTLISTUNIQC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsUNIQ_GHCNM3CRUTEM4_APR2015.txt'
#OUTLISTUNIQH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsUNIQ_CRUTEM4GHCNM3_APR2015.txt'
#OUTPLOT='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/CompareStationListings_GHCNM3HadISDH_APR2015'
#OUTLISTSAMEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsSAME_GHCNM3HadISDH_APR2015.txt'
#OUTLISTSAMEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsSAME_HadISDHGHCNM3_APR2015.txt'
#OUTLISTMAYBEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsMAYBE_GHCNM3HadISDH_APR2015.txt'
#OUTLISTMAYBEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsMAYBE_HadISDHGHCNM3_APR2015.txt'
#OUTLISTUNIQC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsUNIQ_GHCNM3HadISDH_APR2015.txt'
#OUTLISTUNIQH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CompareStationListingsUNIQ_HadISDHGHCNM3_APR2015.txt'


# Variables
CRUTEM_count=0
HadISDH_count=0

CRUTEM_stations=list()	# list of arrays to store CRUTEM station info
HadISDH_stations=list()	# list of arrays to store HadISDH station info

CRUTEM_only=[]	# list of pointers to CRUTEM stations that are unique	
HadISDH_only=[]	# list of pointers to HadISDH stations that are unique	
Same_stations=[] # list of pointers to CRUTEM stations that are in HadISDH
MaybeHadISDH_stations=[] # list of pointers to HadISDH stations that are maybe also in CRUTEM4
MaybeCRUTEM_stations=[] # list of pointers to CRUTEM4 stations that are maybe also in HadISDH

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
# COMPUTE_JACCARD_INDEX
def compute_jaccard_index(set_1, set_2):
    ''' Each word must be converted to a set() '''
    ''' Word1 and Word2 are then compared for similarity on a scale of 0 to 1 '''
    ''' This should work on sets containing arrays too but not sure what the result would mean '''
    ''' NOTE TOKYO and KYOTO have a score of 1!!! '''
    ''' NOTE LNCN and LINCOLN have a score of 1!!! '''
    ''' http://love-python.blogspot.co.uk/2012/07/python-code-to-compute-jaccard-index.html '''
    set_1=set(set_1)
    set_2=set(set_2)
    n = len(set_1.intersection(set_2))

    return n / float(len(set_1) + len(set_2) - n) # COMPUTE_JACCARD_INDEX
#************************************************************************
# GetJaccardList
def GetJaccardList(list1,list2,jaccouts):	# fullouts,partouts
    ''' Loop through each word in list1 '''
    ''' compare with each word in list2 '''
    ''' find the maximum Jaccard score/scores and store '''
    ''' OLD: if 1.0 then add location of these relative to list2 to the fullouts list of lists'''
    ''' OLD: or set to -999 '''
    ''' OLD: if between 0.7 and 1.0 then add location of these relative to list2 to the partouts list of lists '''
    ''' OLD: or set to -999 '''
    for loo in range(len(list1)):
        tmpJACC=[compute_jaccard_index(list1[loo],list2[i]) for i in range(len(list2))]
        jaccouts.append(np.where(np.array(tmpJACC) == max(tmpJACC))[0])
#        if (max(tmpJACC) == 1.0): 
#            fullouts.append(np.where(np.array(tmpJACC) == max(tmpJACC))[0])	# this could be an array in some cases
#        else:
#            fullouts.append(-999)
#        if (max(tmpJACC) > 0.7) & (max(tmpJACC) < 1.0): 
#            partouts.append(np.where(np.array(tmpJACC) == max(tmpJACC))[0]) # this could be an array in some cases
#        else:
#            partouts.append(-999)

    return jaccouts # GETJACCARDLIST fullouts,partouts 
#************************************************************************
# GetMatchingStatsList
def GetMatchingStatsList(names1,lats1,lons1,elevs1,names2,lats2,lons2,elevs2,uniques,sames,maybes):
    ''' Loop through each location in lats1/lons1/elevs1 '''
    ''' Compare with location in lats2/lons2/elevs2 '''
    ''' Look also at full and part Jaccard Index for this pair '''
    ''' find the locations that match within +/- 0.1 degree and 100m SAME ''' 
    ''' find the locations that match within +/- 0.1 degree but not 100m MAYBE '''
    ''' find the identical station names and locations within +/- 1 degree SAME '''
    ''' find the nearly indentical station names and locations within 1 degree MAYBE '''
    for loo in range(len(names1)):
        foundone=0	# set to one if find a potential match
	for goo in range(len(names2)):
	    if (abs(lons1[loo] - lons2[goo]) <= 1.) & (abs(lats1[loo] - lats2[goo]) <= 1.):	# potential same or maybe
	        jaccscore=compute_jaccard_index(names1[loo],names2[goo])
		if (abs(lons1[loo] - lons2[goo]) <= 0.1) & (abs(lats1[loo] - lats2[goo]) <= 0.1): 	# potential same or maybe
                    if (abs(elevs1[loo] - elevs2[goo]) <= 100):
		        sames.append([loo,goo])	# this is a same, append pointer pairs 
		    else:
		        maybes.append([loo,goo])	# this is a same, append pointer pairs 
		    foundone=1
		elif (jaccscore == 1.):	
		    sames.append([loo,goo])	# this is a same, append pointer pairs 
		    foundone=1
		elif (jaccscore > 0.7):
		    maybes.append([loo,goo])	# this is a same, append pointer pairs 
		    foundone=1
        if foundone == 0:
	    uniques.append(loo)		# this is defo a unique station
		
    return uniques,sames,maybes # GETMATCHINGLOCSLIST
#************************************************************************
# PlotStationLocs
def PlotStationLocs(TheFile,C_only,H_only,SameCinH,SameHinC,C_maybe,H_maybe,C_list,H_list,TheLetter,TheNamee):
    ''' Plot three maps: unique, unique, sames, maybes '''
    ''' Label plot with totals '''
    ''' Save as png and eps '''

    # set up dimensions and plot - this is a 2 column nvar rows plot
    nplots=4
    xpos=[0.025,0.525,0.025,0.525]
    ypos=[0.525,0.525,0.025,0.025]
    xfat=[0.45,0.45,0.45,0.45]
    ytall=[0.45,0.45,0.45,0.45]
      
    f=plt.figure(4,figsize=(12,8))	#6,18

    plt.axes([xpos[0],ypos[0],xfat[0],ytall[0]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
        
    # plot blue dots for HadISDH
    x,y = m(H_list[0][H_only],H_list[1][H_only])	# long and lat
    m.scatter(x,y,s=5,marker='o',color="DodgerBlue")
    
    ax1=plt.axes([xpos[0],ypos[0]-0.02,xfat[0],ytall[0]*0.1],frameon=False) # map only
    ax1.set_ylim(0,1)
    ax1.set_xlim(0,1)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.annotate(str(len(H_only))+" HadISDH stations",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate(str(len(H_only))+" CRUTEM4 stations",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.05,0.6,markersize=10,marker='o',color='DodgerBlue')

    plt.figtext(xpos[0],ypos[0]+ytall[0],TheLetter[0],size=18)
    plt.figtext(xpos[0]+xfat[0]/2.,ypos[0]+ytall[0],TheNamee[0],size=18,ha='center')    

    plt.axes([xpos[1],ypos[1],xfat[1],ytall[1]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
    
    # plot red dots for CRUTEM
    x,y = m(C_list[0][C_only],C_list[1][C_only])	# long and lat
    m.scatter(x,y,s=2,marker='o',color="Firebrick")
    
    ax1=plt.axes([xpos[1],ypos[1]-0.02,xfat[1],ytall[1]*0.1],frameon=False) # map only
    ax1.set_ylim(0,1)
    ax1.set_xlim(0,1)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.annotate(str(len(C_only))+" CRUTEM4 stations",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate(str(len(C_only))+" GHCNM3 stations",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.05,0.6,markersize=10,marker='o',color='Firebrick')

    plt.figtext(xpos[1],ypos[1]+ytall[1],TheLetter[1],size=18)
    plt.figtext(xpos[1]+xfat[0]/2.,ypos[1]+ytall[1],TheNamee[1],size=18,ha='center')    


    plt.axes([xpos[2],ypos[2],xfat[2],ytall[2]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
    
    # plot red dots for sames
    pointees=[SameCinH[i][0] for i in range(len(SameCinH))]
    x,y = m(H_list[0][pointees],H_list[1][pointees])	# long and lat
    m.scatter(x,y,s=2,marker='o',color="Firebrick")
    
    # how many HadISDH stations with more than one match?
    uniqies=itemfreq(pointees)
    countdupsH=0
    for loo in range(len(uniqies)):
        if (uniqies[loo][1] > 1):
	    countdupsH=countdupsH+1
    pointees=[SameHinC[i][0] for i in range(len(SameHinC))]
    uniqies=itemfreq(pointees)
    countdupsC=0
    for loo in range(len(uniqies)):
        if (uniqies[loo][1] > 1):
	    countdupsC=countdupsC+1

    ax1=plt.axes([xpos[2],ypos[2]-0.02,xfat[2],ytall[2]*0.1],frameon=False) # map only
    ax1.set_ylim(0,1)
    ax1.set_xlim(0,1)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.annotate(str(len(SameCinH))+" Matching stations ("+str(countdupsH)+", "+str(countdupsC)+")",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate("HadISDH has "+str(countdupsH)+" multiple CRUTEM4 matches",xy=(0.07,0.51),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate("CRUTEM4 has "+str(countdupsC)+" multiple HadISDH matches",xy=(0.07,0.21),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.05,0.6,markersize=10,marker='o',color='Firebrick')

    plt.figtext(xpos[2],ypos[2]+ytall[2],TheLetter[2],size=18)
    plt.figtext(xpos[2]+xfat[0]/2.,ypos[2]+ytall[2],TheNamee[2],size=18,ha='center')    

    plt.axes([xpos[3],ypos[3],xfat[3],ytall[3]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
        
    # plot blue dots for maybes
    pointees=[C_maybe[i][0] for i in range(len(C_maybe))]
    x,y = m(H_list[0][pointees],H_list[1][pointees])	# long and lat
    m.scatter(x,y,s=5,marker='o',color="DodgerBlue")

    # how many HadISDH stations with more than one match?
    uniqies=itemfreq(pointees)
    countdupsH=0
    for loo in range(len(uniqies)):
        if (uniqies[loo][1] > 1):
	    countdupsH=countdupsH+1
    pointees=[H_maybe[i][0] for i in range(len(H_maybe))]
    uniqies=itemfreq(pointees)
    countdupsC=0
    for loo in range(len(uniqies)):
        if (uniqies[loo][1] > 1):
	    countdupsC=countdupsC+1

    ax1=plt.axes([xpos[3],ypos[3]-0.02,xfat[3],ytall[3]*0.1],frameon=False) # map only
    ax1.set_ylim(0,1)
    ax1.set_xlim(0,1)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.annotate(str(len(H_maybe))+" Potentially matching stations ("+str(countdupsH)+", "+str(countdupsC)+")",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate("HadISDH has "+str(countdupsH)+" multiple CRUTEM4 matches",xy=(0.07,0.51),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate("CRUTEM4 has "+str(countdupsC)+" multiple HadISDH matches",xy=(0.07,0.21),xytext=None, xycoords='axes fraction',color="Black",size=12)
    plt.plot(0.05,0.6,markersize=10,marker='o',color='DodgerBlue')

    plt.figtext(xpos[3],ypos[3]+ytall[3],TheLetter[3],size=18)
    plt.figtext(xpos[3]+xfat[0]/2.,ypos[3]+ytall[3],TheNamee[3],size=18,ha='center')    
    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotStationLocs
    
#************************************************************************
# WriteOut
def WriteOut(TheFile,TwoUp,Pointers,Datalist1,Datalist2=0):
    ''' Use Pointers to point to each station '''
    ''' Write station info out to file given '''
    ''' Input DataList components: '''
    ''' 4,3,1,0,2 = ids, names, lats, longs, elevs '''

    filee=open(TheFile,'a+')

    for outt in range(len(Pointers)):
        if TwoUp == 1:
#            filee.write('%11s %30s %7.3f %8.3f %7.1f \n' % (IDs[Pointers[outt]],Names[Pointers[outt]],
#	                Lats[Pointers[outt]],Longs[Pointers[outt]],Elevss[Pointers[outt]])
            filee.write(str('{:11s}'.format(Datalist1[4][Pointers[outt]])+' '+'{:30s}'.format(Datalist1[3][Pointers[outt]])+' '+'{:7.3f}'.format(Datalist1[1][Pointers[outt]])+' '+'{:8.3f}'.format(Datalist1[0][Pointers[outt]])+' '+'{:7.1f}'.format(Datalist1[2][Pointers[outt]])+'\n'))
        else:
            filee.write(str('{:11s}'.format(Datalist1[4][Pointers[outt][0]])+' '+'{:30s}'.format(Datalist1[3][Pointers[outt][0]])+' '+'{:7.3f}'.format(Datalist1[1][Pointers[outt][0]])+' '+'{:8.3f}'.format(Datalist1[0][Pointers[outt][0]])+' '+'{:7.1f}'.format(Datalist1[2][Pointers[outt][0]])+' '+'{:11s}'.format(Datalist2[4][Pointers[outt][1]])+' '+'{:30s}'.format(Datalist2[3][Pointers[outt][1]])+' '+'{:7.3f}'.format(Datalist2[1][Pointers[outt][1]])+' '+'{:8.3f}'.format(Datalist2[0][Pointers[outt][1]])+' '+'{:7.1f}'.format(Datalist2[2][Pointers[outt][1]])+'\n'))

    filee.close()

    return #WriteOut
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list for HadISDH
MyTypes=("|S11","float","float","float","|S4","|S30","|S7")
MyDelimiters=[11,9,10,7,4,30,7]
MyFile=INFILEH
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
HadISDH_stations.append(np.array(RawData['f2']))	# longitudes
HadISDH_stations.append(np.array(RawData['f1']))	# latitudes
HadISDH_stations.append(np.array(RawData['f3']))	# elevations
HadISDH_stations.append(np.array(RawData['f5']))	# names
HadISDH_stations.append(np.array(RawData['f0']))	# ids
HadISDH_count=len(HadISDH_stations[0])

#MyTypes=("|S6","|S3","|S20","|S15","float","float","int","int","int") # need to remove '--' characters from station name
#MyDelimiters=[6,3,20,15,7,7,6,5,5]
#MyFile=INFILEH
#RawData=ReadData(MyFile,MyTypes,MyDelimiters)
#HadISDH_stations.append(-(np.array(RawData['f5'])))	# longitudes
#HadISDH_stations.append(np.array(RawData['f4']))	# latitudes
#HadISDH_stations.append(np.array(RawData['f6']))	# elevations
#HadISDH_stations.append(np.array(RawData['f2']))	# names
#HadISDH_stations.append(np.array(RawData['f0']))	# ids
#HadISDH_count=len(HadISDH_stations[0])

MyTypes=("|S6","|S3","|S20","|S15","float","float","int","int","int") # need to remove '--' characters from station name
MyDelimiters=[6,3,20,15,7,7,6,5,5]
MyFile=INFILEC
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
CRUTEM_stations.append(-(np.array(RawData['f5'])))	# longitudes WHY ARE THESE IN REVERSE???
CRUTEM_stations.append(np.array(RawData['f4']))	# latitudes
CRUTEM_stations.append(np.array(RawData['f6']))	# elevations
CRUTEM_stations.append(np.array(RawData['f2']))	# names
CRUTEM_stations.append(np.array(RawData['f0']))	# ids
CRUTEM_count=len(CRUTEM_stations[0])

#MyTypes=("|S11","float","float","float","|S1","|S20","|S49") # need to remove '--' characters from station name
#MyDelimiters=[11,9,10,7,1,20,49]
#MyFile=INFILEC
#RawData=ReadData(MyFile,MyTypes,MyDelimiters)
#CRUTEM_stations.append(np.array(RawData['f2']))	# longitudes WHY ARE THESE IN REVERSE???
#CRUTEM_stations.append(np.array(RawData['f1']))	# latitudes
#CRUTEM_stations.append(np.array(RawData['f3']))	# elevations
#CRUTEM_stations.append(np.array(RawData['f5']))	# names
#CRUTEM_stations.append(np.array(RawData['f0']))	# ids
#CRUTEM_count=len(CRUTEM_stations[0])
print('Read in data...')

# remove '--' and WHITESPACE characters from each station name within CRUTEM and HadISDH
CRUTEM_stations[3]=[re.sub('[ -]','',CRUTEM_stations[3][i]) for i in range(CRUTEM_count)]
HadISDH_stations[3]=[re.sub('[ -]','',HadISDH_stations[3][i]) for i in range(HadISDH_count)]
print('Removed whitespace and dashes...')

# find station names that match exactly, or are similar - using JACCARD INDEX
##JACCARD_FULLCinH=list()        
##JACCARD_FULLHinC=list()        
##JACCARD_PARTCinH=list()        
##JACCARD_PARTHinC=list()        
##JACCARD_FULLHinC,JACCARD_PARTHinC=GetJaccardList(HadISDH_stations[3],CRUTEM_stations[3],JACCARD_FULLHinC,JACCARD_PARTHinC)
##JACCARD_FULLCinH,JACCARD_PARTCinH=GetJaccardList(CRUTEM_stations[3],HadISDH_stations[3],JACCARD_FULLCinH,JACCARD_PARTCinH)
#JACCARD_CinH=list()        
#JACCARD_HinC=list()        
#JACCARD_HinC=GetJaccardList(HadISDH_stations[3],CRUTEM_stations[3],JACCARD_HinC)
##JACCARD_CinH=GetJaccardList(CRUTEM_stations[3],HadISDH_stations[3],JACCARD_CinH)

# find the locations that match within +/- 0.1 degree and 100m
# find the locations that match within +/- 0.1 degree but not 100m
# find the identical station names and locations within +/- 1 degree
# find the nearly indentical station names and locations within 1 degree
CRUTEM_only=list()	# list of pointers to CRUTEM stations that are unique	
HadISDH_only=list()	# list of pointers to HadISDH stations that are unique	
Same_stationsH=list() # list of pointers to CRUTEM stations that are in HadISDH
Same_stationsC=list() # list of pointers to HadISDH stations that are in CRUTEM
MaybeHadISDH_stations=list() # list of pointers to HadISDH stations that are maybe also in CRUTEM4
MaybeCRUTEM_stations=list() # list of pointers to CRUTEM4 stations that are maybe also in HadISDH
HadISDH_only,Same_stationsH,MaybeCRUTEM_stations=GetMatchingStatsList(HadISDH_stations[3],HadISDH_stations[1],
                                                              HadISDH_stations[0],HadISDH_stations[2],
							      CRUTEM_stations[3],CRUTEM_stations[1],
                                                              CRUTEM_stations[0],CRUTEM_stations[2],
							      HadISDH_only,Same_stationsH,MaybeCRUTEM_stations)
CRUTEM_only,Same_stationsC,MaybeHadISDH_stations=GetMatchingStatsList(CRUTEM_stations[3],CRUTEM_stations[1],
	                                                      CRUTEM_stations[0],CRUTEM_stations[2],
							      HadISDH_stations[3],HadISDH_stations[1],
                                                              HadISDH_stations[0],HadISDH_stations[2],  
                                                              CRUTEM_only,Same_stationsC,MaybeHadISDH_stations)
print('Found matches...')

#stop()
# Write out these lists to file: 2 means output both sources, 1 means just 1
WriteOut(OUTLISTSAMEC,2,Same_stationsC,CRUTEM_stations,HadISDH_stations)
WriteOut(OUTLISTSAMEH,2,Same_stationsH,HadISDH_stations,CRUTEM_stations)
WriteOut(OUTLISTMAYBEC,2,MaybeHadISDH_stations,CRUTEM_stations,HadISDH_stations)
WriteOut(OUTLISTMAYBEH,2,MaybeCRUTEM_stations,HadISDH_stations,CRUTEM_stations)
WriteOut(OUTLISTUNIQC,1,CRUTEM_only,CRUTEM_stations)
WriteOut(OUTLISTUNIQH,1,HadISDH_only,HadISDH_stations)
print('Written out data...')
							      
# pass to plotter
Letty=['a)','b)','c)','d)']
Namey=['Unique HadISDH Stations','Unique CRUTEM4 Stations','Matching Stations','Potentially Matching Stations']
#Namey=['Unique CRUTEM4 Stations','Unique GHCNM3 Stations','Matching Stations','Potentially Matching Stations']
#Namey=['Unique HadISDH Stations','Unique GHCNM3 Stations','Matching Stations','Potentially Matching Stations']

PlotStationLocs(OUTPLOT,CRUTEM_only,HadISDH_only,Same_stationsH,Same_stationsC,
                        MaybeCRUTEM_stations,MaybeHadISDH_stations,CRUTEM_stations,HadISDH_stations,
			Letty,Namey)
		
#    stop()

print("And, we are done!")

