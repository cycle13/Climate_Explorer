#!/usr/local/sci/bin/python

#***************************************
# 3rd June 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# Given a list of gridbox centres (assuming 5by5 here)
# Given a list of stations (HadISDH, CRUTEM, GHCNM)
# Finds all stations within that gridbox
# Outputs list to file

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetStationsInGridBoxes_JUN2015.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import numpy as np
import sys, os

# Input Files
INFILEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/PosthomogIDPHAt_goodsHadISDH.2.0.1.2014p_JAN2015.txt'
INFILEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CRUTEM4_StationList_APR2015.txt'
INFILEG='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/ghcnm.tavg.v3.2.2.20150430.qca.inv'

# Output Files
#OUTLIST='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/StationsInGridboxes_HadISDH.landTCRUTEM4GHCNM3_JUN2015.txt'
OUTLIST='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/StationsInGridboxesCOR_HadISDH.landTCRUTEM4GHCNM3_JUN2015.txt'

# Variables
# gridbox lists
# What gridbox to pull out (use long (-180 to 180) and lat (-90 to 90) centres?
# If this is greater than one element then make multiple time series
# HadISDH vs CRUTEM
#nGBs=list([13,11,8])
nGBs=list([10,10,10])
TotalGBs=np.sum(nGBs)
WhichPair=np.concatenate((np.concatenate((np.repeat('H/C',nGBs[0]),np.repeat('H/G',nGBs[1]))),np.repeat('C/G',nGBs[2])))
#LongPoints=list([-57.500,157.500,167.500,-67.500,-67.500,-132.500,147.500,-62.500,-82.500,122.500,147.500,112.500,112.500,37.500,-172.500,177.500,167.500,-167.500,147.500,-67.500,-157.500,112.500,112.500,147.500,177.500,-167.500,-82.500,-132.500,-172.500,157.500,127.500,-57.500])
#LatPoints=list([7.500,7.500,-47.500,-52.500,-17.500,-22.500,-22.500,-42.500,-7.500,-27.500,-17.500,-67.500,-27.500,-67.500,57.500,-12.500,-47.500,62.500,-22.500,-52.500,17.500,-67.500,-27.500,-17.500,-12.500,62.500,-7.500,-22.500,57.500,7.500,-17.500,7.500])
LongPoints=list([107.5,102.5,-67.5,37.5,-72.5,-67.5,122.5,-62.5,-57.5,-177.5,17.5,-72.5,-67.5,-77.5,17.5,37.5,-37.5,12.5,-37.5,82.5,-77.5,17.5,137.5,37.5,177.5,-132.5,-62.5,157.5,-67.5,-77.5])
LatPoints=list([-7.5,-2.5,7.5,2.5,-17.5,-2.5,7.5,2.5,7.5,-17.5,-27.5,-17.5,7.5,2.5,2.5,-2.5,-7.5,-2.5,-12.5,22.5,2.5,-27.5,7.5,-2.5,-12.5,-22.5,-2.5,7.5,7.5,17.5])

CDRLabels=list(['HadISDH','CRUTEM','GHCNM'])

CRUTEM_count=0
GHCNM_count=0
HadISDH_count=0

GHCNM_stations=list()	# list of arrays to store GHCNM station info
CRUTEM_stations=list()	# list of arrays to store CRUTEM station info
HadISDH_stations=list()	# list of arrays to store HadISDH station info

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
# FindStations
def FindStations(StationList,Slat,Nlat,Wlon,Elon,FoundStationsList):
    ''' Loop through the given list of station info '''
    ''' Append to FoundStationList all info when a station is found within the box '''
    ''' Station lat/lon can be within or equal to gridbox boundaries '''
    # find lats and lons that are within the box
    # now find where these match (may be able to quadruple loop the where
    Gots=np.where((StationList[0] >= Wlon) & (StationList[0] <= Elon) & (StationList[1] >= Slat) & (StationList[1] <= Nlat))[0]
    
    if (len(Gots) > 0):
        for loo in range (5):
	    FoundStationsList.append(StationList[loo][Gots])
    
    return FoundStationsList # FindStations
#************************************************************************
# WriteOut
def WriteOut(TheFile,Lat,Lon,TheLabels,FoundH,FoundC,FoundG):
    ''' Append station info out to file given '''

    filee=open(TheFile,'a+')
    countH=0
    countC=0
    countG=0
    if (len(FoundH) > 0):
        countH=len(FoundH[0])
    if (len(FoundC) > 0):
        countC=len(FoundC[0])
    if (len(FoundG) > 0):
        countG=len(FoundG[0])
    
    filee.write(str('Long: '+'{:8.3f}'.format(Lon)+' Lat: '+'{:7.3f}'.format(Lat)+' HadISDH: '+'{:3d}'.format(countH)+' CRUTEM: '+'{:3d}'.format(countC)+' GHCNM: '+'{:3d}'.format(countG)+'\n'))
    
    if (countH > 0):
        for outt in range(len(FoundH[0])):
            filee.write(str('{:11s}'.format(FoundH[4][outt])+' '+'{:8.3f}'.format(FoundH[0][outt])+' '+'{:7.3f}'.format(FoundH[1][outt])+' '+'{:7.1f}'.format(FoundH[2][outt])+' '+'{:30s}'.format(FoundH[3][outt])+' '+'{:7s}'.format(TheLabels[0])+'\n'))

    if (countC > 0):
        for outt in range(len(FoundC[0])):
            filee.write(str('{:11s}'.format(FoundC[4][outt])+' '+'{:8.3f}'.format(FoundC[0][outt])+' '+'{:7.3f}'.format(FoundC[1][outt])+' '+'{:7.1f}'.format(FoundC[2][outt])+' '+'{:30s}'.format(FoundC[3][outt])+' '+'{:7s}'.format(TheLabels[1])+'\n'))

    if (countG > 0):
        for outt in range(len(FoundG[0])):
            filee.write(str('{:11s}'.format(FoundG[4][outt])+' '+'{:8.3f}'.format(FoundG[0][outt])+' '+'{:7.3f}'.format(FoundG[1][outt])+' '+'{:7.1f}'.format(FoundG[2][outt])+' '+'{:30s}'.format(FoundG[3][outt])+' '+'{:7s}'.format(TheLabels[2])+'\n'))

    filee.close()

    return #WriteOut
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list for HadISDH and append to list of arrays
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

# read in station list for CRUTEM4 and append to list of arrays
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

# read in station list for GHCNM3 and append to list of arrays
MyTypes=("|S11","float","float","float","|S1","|S20","|S49") # need to remove '--' characters from station name
MyDelimiters=[11,9,10,7,1,20,49]
MyFile=INFILEG
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
GHCNM_stations.append(np.array(RawData['f2']))	# longitudes WHY ARE THESE IN REVERSE???
GHCNM_stations.append(np.array(RawData['f1']))	# latitudes
GHCNM_stations.append(np.array(RawData['f3']))	# elevations
GHCNM_stations.append(np.array(RawData['f5']))	# names
GHCNM_stations.append(np.array(RawData['f0']))	# ids
GHCNM_count=len(GHCNM_stations[0])
print('Read in data...')

# Loop through each gridbox and then through each station list to find stations within
for gg in range(TotalGBs):
    latS=LatPoints[gg]-2.5
    latN=LatPoints[gg]+2.5
    lonW=LongPoints[gg]-2.5
    lonE=LongPoints[gg]+2.5
    
    Found_Hs=list() # list of arrays for ID, Lat, Lon, Elev, Name
    Found_Cs=list() # list of arrays for ID, Lat, Lon, Elev, Name
    Found_Gs=list() # list of arrays for ID, Lat, Lon, Elev, Name

    # Loop through HadISDH
    Found_Hs=FindStations(HadISDH_stations,latS,latN,lonW,lonE,Found_Hs)
    if (len(Found_Hs) > 0):
        count=len(Found_Hs[0])
    else:
        count=0
    print('HadISDH: ',count)
    
    # Loop through CRUTEM
    Found_Cs=FindStations(CRUTEM_stations,latS,latN,lonW,lonE,Found_Cs)    
    if (len(Found_Cs) > 0):
        count=len(Found_Cs[0])
    else:
        count=0
    print('CRUTEM: ',count)
    
    # Loop through GHCNM
    Found_Gs=FindStations(GHCNM_stations,latS,latN,lonW,lonE,Found_Gs)
    if (len(Found_Gs) > 0):
        count=len(Found_Gs[0])
    else:
        count=0
    print('GHCNM: ',count)

    # Write to file by appending
    print('Writing to file....',gg)
    WriteOut(OUTLIST,LatPoints[gg],LongPoints[gg],CDRLabels,Found_Hs,Found_Cs,Found_Gs)

#    stop()

print("And, we are done!")

