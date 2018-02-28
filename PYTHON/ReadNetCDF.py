#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 16 October 2015
# Last update: 16 October 2015
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This reads in a netCDF file and outputs numpy arrays or lists of numpy arrays:
#    GetGrid: of a time,lat,lon gridded dataset
#    GetField: of a lat,lon gridded field of e.g trends
#    GetTS: of a time, point/station dataset 
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Inbuilt:
# import numpy as np
# import scipy.stats
# import itertools
# from scipy.io import netcdf
## from netCDF4 import Dataset
# import netCDF4 as nc4  
# import pdb # pdb.set_trace() or c
#
# Kate's:
#
# -----------------------
# DATA
# -----------------------
# The code requires a netCDF file
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# GetGrid:
#    INPUTS:    
#    # A scalar string containing the filepath+filename
#    FileName = '/MyDir/MyFile.nc'
#    # A LIST of strings with the names of each variable to read in
#    ReadInfo = ['t_anoms'[ # ['t_anoms','t_abs']
#    OPTIONAL LatInfo and LonInfo lists to build LatList and LonList if read in isn't going to work
#    Default reads in lat and lon variables assuming LatInfo=['latitude'],LonInfo=['longitude']
#    Option 1: just provide the variable names to read in:
#    LatInfo = ['latitude'] # variable name, number of latitudes, start latitude
#    LonInfo = ['longitude'] # variable name, number of longitudes, start longitude
#    Option 2: provide number of boxes and start lat/lon:
#    LatInfo = [36,-87.5] # number of latitudes, start latitude
#    LonInfo = [72,-177.5] # number of longitudes, start longitude
#
#    python2.7 
#    from ReadNetCDF import GetGrid
#
#    # Premaking of arrays not strictly necessary:
#    TmpData = np.reshape(np.repeat(mdi,NMons*NLats*NLons),(NMons,NLats,NLons)) # entire field of data read in
#    LatList = [] # list of latitude gridbox centres
#    LonList = [] # list of longitude gridbox centres
#    # Multiple grid read in case at little more complicated, and still not necessary but nice?:
#    import itertools
#    TmpData=list(itertools.repeat(np.reshape(np.repeat(mdi,NMons*NLats*NLons),(NMons,NLats,NLons)),len(ReadInfo)))
#
#    TmpData,LatList,LonList = ReadNetCDF.GetGrid(FileName,ReadInfo,LatInfo,LonInfo)
#    # Multiple grid read in case requires unpacking of the list of np.arrays
#    TmpData1=TmpData[0] # etc.
#
# GetField:
#
# GetTS:
#
# -----------------------
# OUTPUT
# -----------------------
# GetGrid:
#    TmpData: a numpy array of the 3D gridded field as time, lat, lon or a list of numpy arrays to be unpacked
#    LatList: a numpy array of latitude gridbox centres
#    LonList: a numpy array of longitude gridbox centre
#
# GetField:
#
# GetTS:
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 2 (28th February 2018)
# ---------
#  
# Enhancements
#  
# Changes
# I've added a new GetGrid4 which does the same as GetGrid but works using netCDF4.
# This 'should' then mean that it can work with netCDF4 format and also automatically scale/offset data on read in
# Later I may make it such that you can pull out the variable entirely rather than just the data array so that
# you can access all attributes
#  
# Bug fixes
#
#
# Version 1 (16th October 2015)
# ---------
#  
# Enhancements
#  
# Changes
#  
# Bug fixes
#  
# -----------------------
# OTHER INFORMATION
# -----------------------
#    
#************************************************************************
# Functions
#************************************************************************
# GETGRID
def GetGrid(FileName,
            ReadInfo,
	    LatInfo = ['latitude'],
	    LonInfo = ['longitude']):
    ''' Open the NetCDF File
        Get the list of latitudes either from file or derive
        Get the list of longitudes either from file or derive
        Get the data (1 to 10 fields)
	INPUTS:
	FileName: string containing filepath/name
	ReadInfo: a LIST of string variable names (1+) for the grid to read in
	LatInfo: 
	    LatInfo=['latitude'[ # DEFAULT
	    LatInfo=['lat'[ # alternative LIST of string variable name
	    LatInfo=[36,-87.5] # list containing number of lats, a float for the start latitude
	LatInfo: 
	    LatInfo=['longitude'[ # DEFAULT
	    LatInfo=['lon'[ # alternative LIST of string variable name
	    LatInfo=[72,-177.5] # list containing number of lons, a float for the start longitude 
	OUTPUTS:
	LatList: a numpy array of latitude gridbox centres
	LonList: a numpy array of longitude gridbox centres
	TheData:
	    1 variable: a numpy array of time,lat,lon
	    2+ variables: a list of numpy arrays of time,lat,lon to be unpacked '''

    # Set up python imports
    import numpy as np
    import scipy.stats
    from scipy.io import netcdf
    import pdb # pdb.set_trace() or c
    
#    print(FileName,ReadInfo,LatInfo,LonInfo)

    # Open the netCDF file
    ncf=netcdf.netcdf_file(FileName,'r')

    # ncf.variables this lists the variable names
    
    # If LatInfo is only 1 element long then read in the variable, else calculate using given nlats, start_lat
    if (len(LatInfo) == 1):
        var=ncf.variables[LatInfo[0]]
        TheLatList=np.array(var.data)
    else:	
        if (LatInfo[1] < 0):
            TheLatList=np.arange(LatInfo[1], LatInfo[1]+180.,(180./LatInfo[0]))
        else:
            TheLatList=np.arange(LatInfo[1], LatInfo[1]-180.,-(180./LatInfo[0]))    

    # If LonInfo is only 1 element long then read in the variable, else calculate using given nlons, start_lon
    if (len(LonInfo) == 1):
        var=ncf.variables[LonInfo[0]]
        TheLonList=np.array(var.data)
    else:
        if (LonInfo[1] < 10):
            TheLonList=np.arange(LonInfo[1], LonInfo[1]+360.,(360./LonInfo[0]))
        else:
            TheLonList=np.arange(LonInfo[1], LonInfo[1]-360.,-(360./LonInfo[0]))    

    # If ReadInfo is only 1 element long then read into a numpy array, else make a list of arrays and then read in all
    if (len(ReadInfo) == 1):
        var=ncf.variables[ReadInfo[0]]
        TheData=np.array(var.data) # times, lats, lons
	#pdb.set_trace()
    else:
        # Initialise TheData as a list
	TheData=[]
	for loo in range(len(ReadInfo)):
            var=ncf.variables[ReadInfo[loo]]
	    #pdb.set_trace()
            TmpData=np.array(var.data) # times, lats, lons
	    TheData.append(TmpData)	

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # GetGrid

#************************************************************************
#************************************************************************
# GETGRID4
def GetGrid4(FileName,
            ReadInfo,
	    LatInfo = ['latitude'],
	    LonInfo = ['longitude']):
    ''' Open the NetCDF File
        Get the list of latitudes either from file or derive
        Get the list of longitudes either from file or derive
        Get the data (1 to 10 fields)
	INPUTS:
	FileName: string containing filepath/name
	ReadInfo: a LIST of string variable names (1+) for the grid to read in
	LatInfo: 
	    LatInfo=['latitude'[ # DEFAULT
	    LatInfo=['lat'[ # alternative LIST of string variable name
	    LatInfo=[36,-87.5] # list containing number of lats, a float for the start latitude
	LatInfo: 
	    LatInfo=['longitude'[ # DEFAULT
	    LatInfo=['lon'[ # alternative LIST of string variable name
	    LatInfo=[72,-177.5] # list containing number of lons, a float for the start longitude 
	OUTPUTS:
	LatList: a numpy array of latitude gridbox centres
	LonList: a numpy array of longitude gridbox centres
	TheData:
	    1 variable: a numpy array of time,lat,lon
	    2+ variables: a list of numpy arrays of time,lat,lon to be unpacked '''

    # Set up python imports
    import numpy as np
    import scipy.stats
    from scipy.io import netcdf
    #from netCDF4 import Dataset  
    import netCDF4 as nc4  
    import pdb # pdb.set_trace() or c
    
#    print(FileName,ReadInfo,LatInfo,LonInfo)

    # Open the netCDF file
    ncf=nc4.Dataset(FileName,'r')

    # ncf.variables this lists the variable names
    
    # If LatInfo is only 1 element long then read in the variable, else calculate using given nlats, start_lat
    if (len(LatInfo) == 1):
        var = ncf.variables[LatInfo[0]] # loads the data and attributes
        TheLatList = var[:]             # just pulls out the data as a numpy array
	#var.ncattrs() # prints the attributues
    else:	
        if (LatInfo[1] < 0):
            TheLatList=np.arange(LatInfo[1], LatInfo[1]+180.,(180./LatInfo[0]))
        else:
            TheLatList=np.arange(LatInfo[1], LatInfo[1]-180.,-(180./LatInfo[0]))    

    # If LonInfo is only 1 element long then read in the variable, else calculate using given nlons, start_lon
    if (len(LonInfo) == 1):
        var = ncf.variables[LonInfo[0]]
        TheLonList = var[:]
    else:
        if (LonInfo[1] < 10):
            TheLonList=np.arange(LonInfo[1], LonInfo[1]+360.,(360./LonInfo[0]))
        else:
            TheLonList=np.arange(LonInfo[1], LonInfo[1]-360.,-(360./LonInfo[0]))    

    # If ReadInfo is only 1 element long then read into a numpy array, else make a list of arrays and then read in all
    if (len(ReadInfo) == 1):
        var = ncf.variables[ReadInfo[0]]
        TheData = var[:] # times, lats, lons - THIS AUTOMATICALLY APPLIES SCALE AND OFFSET!!!
	#var.ncattrs() # prints the attributues
	#pdb.set_trace()
    else:
        # Initialise TheData as a list
	TheData = []
	for loo in range(len(ReadInfo)):
            var = ncf.variables[ReadInfo[loo]]
	    #var.ncattrs() # prints the attributues
	    #var.add_offset # prints the add_offset attribute
	    #pdb.set_trace()
            TmpData = var[:] # times, lats, lons - THIS AUTOMATICALLY APPLIES SCALE AND OFFSET!!!
	    TheData.append(TmpData)	

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # GetGrid4

#*********************************************************************************
