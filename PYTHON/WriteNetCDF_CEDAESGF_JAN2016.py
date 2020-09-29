#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 8 January 2016
# Last update: 6 April 2020
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This is a module called to write out data to a netCDF file ready for ingestion into the CEDA/ESGF 
# archive. It should be CF compliant to at least CF-1.0. Standard names are used where available and
# variable names are taken from PCMDI standard names (or designed to be PCMDI compatible). 
# 
# A table listing the header info is available for each type of file. For HadISDH they can be found 
# here: www.metoffice.gov.uk/hadobs/hadisdh/formattables.html
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Python:
# import numpy as np
# import sys, os
# import scipy.stats
# import struct
# import os.path
# import datetime as dt
# from datetime import datetime
# from matplotlib.dates import date2num,num2date
# from netCDF4 import Dataset
# from netCDF4 import stringtoarr # for putting strings in as netCDF variables
# from scipy.io import netcdf
# import pdb # pdb.set_trace() or c
#
# Kate's:
# from ReadNetCDF import GetGrid - written by Kate Willett, reads in any netCDF grid, can cope with multiple fields
# 
# -----------------------
# DATA
# -----------------------
# Any numpy arrays of data, latitude and longtitude info, mdi info etc:
#
# Filename: string containing filepath and filename e.g.,
#	/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/huss_HadISDH_HadOBS_19730101-20141231_v2-0-1-2014p.nc
# Dates: a four elemnt intger dictionary of {'StYr': integer e.g., 1973,
#                                            'StMon': integer between 1 and 12,
#                                            'EdYr': integer e.g., 2014,
# 					                         'EdMon': integer between 0 and 11}
# Latitudes: float array of latitude gridbox centres from North to South
# Longitudes: float array of longitude gridbox centres from West to East
# ClimPoints: two element list of start year and end year of climatology period e.g., (1976, 2005)
# DataObject: a list of data arrays to be written out
#	???
# DimObject: a list containing a list of dimension names and then dictionaries of attributes applying to each dimension
#       [['time','month','month_name','latitude','longitude','bound_pairs'], # all can be filled with other given info
#        [{'dimname':'time',
#          'standard_name': 'time',
#	   'long_name': 'time',
#	   'units': 'days since 1973-1-1 00:00:00',
#          'axis': 'T',
#          'calendar': 'gregorian',
#          'start_year','1973s',
#          etc...
# AttrObject: a list of dictionaries of attributes applying to each data array to be written out 
#	[[{'varname': 'huss',
#          'standard_name': 'specific_humidity',
#	   'long_name': 'near surface (~2m) specific humidity',
#          'cell_methods': 'time: mean(interval: 1 month comment: anomaly from climatology) area: mean where land (stations within gridbox)',
#	   'comment': 'gridbox mean monthly mean climate anomalies from stations'
#          'units': 'g/kg'
#          'scale_factor': '100'
# 
#       }],[{}],etc...]
# GlobAttrObject: a dictionary of global attributes
#	
#
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# Pass the filename and required variables to the code which will then create a netCDF file storing everything.
# These variables may have been set up using Convert_CEDAESGF_JAN2016.py, which in turns uses ReadNetCDF.py 
# to read in the data.
# 
# from WriteNetCDF_CEDAESGF_JAN2016 import WriteNCCF
# WriteNCCF(Filename,Dates,Latitudes,Longitudes,ClimPoints,DataObjectDimObject,AttrObject,GlobAttrObject)
# 
# -----------------------
# OUTPUT
# -----------------------
# A CEDA/ESGF compliant netCDF file with the provided filename
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 2 6th Mar 2020
# ---------
#  
# Enhancements
# Can now work with marine data - not completed yet
#  
# Changes
# Now python 3
# Now doesn't apply scale and offset and so saves as float32 with MDI=-1e30
#  
# Bug fixes
#  
# Version 1 8th Jan 2016
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
# This may be amended later to apply compression (netCDF4)
#
# Useful links on CF standards:
# Standard Names:
# http://cfconventions.org/Data/cf-standard-names/29/build/cf-standard-name-table.html
# CF Conventions:
# http://cfconventions.org/
# Cell Methods:
# http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch07s03.html
# PCMDI standard names:
# http://www-pcmdi.llnl.gov/projects/amip/DATASTDS/VARNAMES/main.html
#
#************************************************************************
#                                 START
#************************************************************************
# Set up python imports
import numpy as np
import sys, os
import scipy.stats
import struct
import os.path
import datetime as dt
from datetime import datetime
from matplotlib.dates import date2num,num2date
from netCDF4 import Dataset
from netCDF4 import stringtoarr # for putting strings in as netCDF variables
from scipy.io import netcdf
import pdb # pdb.set_trace() or c

from ReadNetCDF import GetGrid # written by Kate Willett, reads in any netCDF grid, can cope with multiple fields

# Set up hardwired variables
OLDMDI = -1e30
MDI = -999 # just check how this is going to work with storing everything as integers
MonthName = ['January   ',
             'February  ',
	     'March     ',
	     'April     ',
	     'May       ',
	     'June      ',
	     'July      ',
	     'August    ',
	     'September ',
	     'October   ',
	     'November  ',
	     'December  ']
MonthDays = [31,28,31,30,31,30,31,31,30,31,30,31]

#************************************************************************
# Subroutines
#************************************************************************
# MakeDaysSince
def MakeDaysSince(TheStYr,TheStMon,TheEdYr,TheEdMon):
    ''' Take counts of months since styr, stmn (assume 15th day of month) '''
    ''' Work out counts of days since styr,stmn, January - incl leap days '''
    ''' Also work out time boundaries 1st and last day of month '''
    ''' This can cope with incomplete years or individual months '''
    
    # set up arrays for month mid points and month bounds
    DaysArray=np.empty(((TheEdYr-TheStYr)+1)*((TheEdMon-TheStMon)+1))
    BoundsArray=np.empty((((TheEdYr-TheStYr)+1)*((TheEdMon-TheStMon)+1),2))
    
    # make a date object for each time point and subtract start date
    StartDate=datetime(TheStYr,TheStMon,1,0,0,0)	# January
    TheYear=TheStYr
    TheMonth=TheStMon
    for mm in range(len(DaysArray)):
        if (TheMonth < 12):
            DaysArray[mm]=(datetime(TheYear,TheMonth+1,1,0,0,0)-datetime(TheYear,TheMonth,1,0,0,0)).days/2. + (datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days
            BoundsArray[mm,0]=(datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days+1
            BoundsArray[mm,1]=(datetime(TheYear,TheMonth+1,1,0,0,0)-StartDate).days
        else:
            DaysArray[mm]=(datetime(TheYear+1,1,1,0,0,0)-datetime(TheYear,TheMonth,1,0,0,0)).days/2. + (datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days	
            BoundsArray[mm,0]=(datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days+1
            BoundsArray[mm,1]=(datetime(TheYear+1,1,1,0,0,0)-StartDate).days
        TheMonth=TheMonth+1
        if (TheMonth == 13):
            TheMonth=1
            TheYear=TheYear+1
            
    return DaysArray,BoundsArray

#************************************************************************
# WriteNCCF
def WriteNCCF(FileName,Dates,Latitudes,Longitudes,ClimPoints,DataObject,DimObject,AttrObject,GlobAttrObject):
    ''' Sort out the date/times to write out and time bounds '''
    ''' Sort out clim bounds '''
    ''' Sort out lat and long bounds '''
    ''' Convert variables using the obtained scale_factor and add_offset: stored_var=int((var-offset)/scale) '''
    ''' Write to file, set up given dimensions, looping through all potential variables and their attributes, and then the provided dictionary of global attributes '''
    
    # Sort out date/times to write out
    print(Dates)
    TimPoints,TimBounds = MakeDaysSince(Dates['StYr'],Dates['StMon'],Dates['EdYr'],Dates['EdMon'])
    nTims = len(TimPoints)
	
    # Sort out clim bounds - paired strings
    ClimBounds = np.empty((12,2),dtype='|S10')
    for mm in range(12):
        ClimBounds[mm,0] = str(ClimPoints[0])+'-'+str(mm+1)+'-'+str(1)
        ClimBounds[mm,1] = str(ClimPoints[1])+'-'+str(mm+1)+'-'+str(MonthDays[mm])
		
    # Sort out LatBounds and LonBounds
    LatBounds = np.empty((len(Latitudes),2),dtype='float')
    LonBounds = np.empty((len(Longitudes),2),dtype='float')
	
    LatBounds[:,0] = Latitudes - ((Latitudes[1]-Latitudes[0])/2.)
    LatBounds[:,1] = Latitudes + ((Latitudes[1]-Latitudes[0])/2.)

    LonBounds[:,0] = Longitudes - ((Longitudes[1]-Longitudes[0])/2.)
    LonBounds[:,1] = Longitudes + ((Longitudes[1]-Longitudes[0])/2.)	
	
    #pdb.set_trace()
    
#    # No need to convert float data using given scale_factor and add_offset to integers - done within writing program (packV = (V-offset)/scale
#    # Not sure what this does to float precision though...
#    # Change mdi into an integer -999 because these are stored as integers
#    for vv in range(len(DataObject)):
#        DataObject[vv][np.where(DataObject[vv] == OLDMDI)] = MDI

    # Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
    ncfw=Dataset(FileName,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first
    
    # Write out the global attributes
    if ('description' in GlobAttrObject):
        ncfw.description = GlobAttrObject['description']
	#print(GlobAttrObject['description'])
	
    if ('File_created' in GlobAttrObject):
        ncfw.File_created = GlobAttrObject['File_created']

    if ('Title' in GlobAttrObject):
        ncfw.Title = GlobAttrObject['Title']

    if ('Institution' in GlobAttrObject):
        ncfw.Institution = GlobAttrObject['Institution']

    if ('History' in GlobAttrObject):
        ncfw.History = GlobAttrObject['History']

    if ('Licence' in GlobAttrObject):
        ncfw.Licence = GlobAttrObject['Licence']

    if ('Project' in GlobAttrObject):
        ncfw.Project = GlobAttrObject['Project']

    if ('Processing_level' in GlobAttrObject):
        ncfw.Processing_level = GlobAttrObject['Processing_level']

    if ('Acknowledgement' in GlobAttrObject):
        ncfw.Acknowledgement = GlobAttrObject['Acknowledgement']

    if ('Source' in GlobAttrObject):
        ncfw.Source = GlobAttrObject['Source']

    if ('Comment' in GlobAttrObject):
        ncfw.Comment = GlobAttrObject['Comment']

    if ('References' in GlobAttrObject):
        ncfw.References = GlobAttrObject['References']

    if ('Creator_name' in GlobAttrObject):
        ncfw.Creator_name = GlobAttrObject['Creator_name']

    if ('Creator_email' in GlobAttrObject):
        ncfw.Creator_email = GlobAttrObject['Creator_email']

    if ('Version' in GlobAttrObject):
        ncfw.Version = GlobAttrObject['Version']

    if ('doi' in GlobAttrObject):
        ncfw.doi = GlobAttrObject['doi']

    if ('Conventions' in GlobAttrObject):
        ncfw.Conventions = GlobAttrObject['Conventions']

    if ('netcdf_type' in GlobAttrObject):
        ncfw.netcdf_type = GlobAttrObject['netcdf_type']
	
    # Loop through and set up the dimension names and quantities
    for vv in range(len(DimObject[0])):
        ncfw.createDimension(DimObject[0][vv],DimObject[1][vv])
	
    # Go through each dimension and set up the variable and attributes for that dimension if needed
    for vv in range(len(DimObject)-2): # ignore first two elements of the list but count all other dictionaries
        print(DimObject[vv+2]['var_name'])
	
	# NOt 100% sure this works in a loop with overwriting
	# initiate variable with name, type and dimensions
        MyVar = ncfw.createVariable(DimObject[vv+2]['var_name'],DimObject[vv+2]['var_type'],DimObject[vv+2]['var_dims'])
        
	# Apply any other attributes
        if ('standard_name' in DimObject[vv+2]):
            MyVar.standard_name = DimObject[vv+2]['standard_name']
	    
        if ('long_name' in DimObject[vv+2]):
            MyVar.long_name = DimObject[vv+2]['long_name']
	    
        if ('units' in DimObject[vv+2]):
            MyVar.units = DimObject[vv+2]['units']
		   	 
        if ('axis' in DimObject[vv+2]):
            MyVar.axis = DimObject[vv+2]['axis']

        if ('calendar' in DimObject[vv+2]):
            MyVar.calendar = DimObject[vv+2]['calendar']

        if ('start_year' in DimObject[vv+2]):
            MyVar.start_year = DimObject[vv+2]['start_year']

        if ('end_year' in DimObject[vv+2]):
            MyVar.end_year = DimObject[vv+2]['end_year']

        if ('start_month' in DimObject[vv+2]):
            MyVar.start_month = DimObject[vv+2]['start_month']

        if ('end_month' in DimObject[vv+2]):
            MyVar.end_month = DimObject[vv+2]['end_month']

        if ('bounds' in DimObject[vv+2]):
            MyVar.bounds = DimObject[vv+2]['bounds']

        if ('climatology' in DimObject[vv+2]):
            MyVar.climatology = DimObject[vv+2]['climatology']

        if ('point_spacing' in DimObject[vv+2]):
            MyVar.point_spacing = DimObject[vv+2]['point_spacing']
	
	# Provide the data to the variable
        if (DimObject[vv+2]['var_name'] == 'time'):
            MyVar[:] = TimPoints

        if (DimObject[vv+2]['var_name'] == 'bounds_time'):
            MyVar[:,:] = TimBounds

        if (DimObject[vv+2]['var_name'] == 'month'):
            for mm in range(12):
                MyVar[mm,:] = stringtoarr(MonthName[mm],10)

        if (DimObject[vv+2]['var_name'] == 'climbounds'):
            for mm in range(12):
                MyVar[mm,0,:] = stringtoarr(ClimBounds[mm,0],10)
                MyVar[mm,1,:] = stringtoarr(ClimBounds[mm,1],10)

        if (DimObject[vv+2]['var_name'] == 'latitude'):
            MyVar[:] = Latitudes

        if (DimObject[vv+2]['var_name'] == 'bounds_lat'):
            MyVar[:,:] = LatBounds

        if (DimObject[vv+2]['var_name'] == 'longitude'):
            MyVar[:] = Longitudes

        if (DimObject[vv+2]['var_name'] == 'bounds_lon'):
            MyVar[:,:] = LonBounds

    # Go through each variable and set up the variable attributes
    for vv in range(len(AttrObject)): # ignore first two elements of the list but count all other dictionaries

        print(AttrObject[vv]['var_name'])

        # NOt 100% sure this works in a loop with overwriting
	# initiate variable with name, type and dimensions
        MyVar = ncfw.createVariable(AttrObject[vv]['var_name'],AttrObject[vv]['var_type'],AttrObject[vv]['var_dims'],zlib=True,fill_value = AttrObject[vv]['_FillValue'])
        
	# Apply any other attributes
        if ('standard_name' in AttrObject[vv]):
            MyVar.standard_name = AttrObject[vv]['standard_name']
	    
        if ('long_name' in AttrObject[vv]):
            MyVar.long_name = AttrObject[vv]['long_name']
	    
# Too many issues with CF compliance
#        if ('cell_methods' in AttrObject[vv]):
#            MyVar.cell_methods = AttrObject[vv]['cell_methods']
	    
        if ('comment' in AttrObject[vv]):
            MyVar.comment = AttrObject[vv]['comment']
	    
        if ('units' in AttrObject[vv]):
            MyVar.units = AttrObject[vv]['units']
		   	 
        if ('axis' in AttrObject[vv]):
            MyVar.axis = AttrObject[vv]['axis']

#        if ('add_offset' in AttrObject[vv]):
#            MyVar.add_offset = AttrObject[vv]['add_offset']
#
#        if ('scale_factor' in AttrObject[vv]):
#            MyVar.scale_factor = AttrObject[vv]['scale_factor']

#        if ('valid_min' in AttrObject[vv]):
#            MyVar.valid_min = AttrObject[vv]['valid_min']#
#
#        if ('valid_max' in AttrObject[vv]):
#            MyVar.valid_max = AttrObject[vv]['valid_max']

#        if ('missing_value' in AttrObject[vv]):
#            MyVar.missing_value = AttrObject[vv]['missing_value']

#        if ('_FillValue' in AttrObject[vv]):
#	    MyVar._FillValue = AttrObject[vv]['_FillValue']

        if ('reference_period' in AttrObject[vv]):
            MyVar.reference_period = AttrObject[vv]['reference_period']

        if ('ancillary_variables' in AttrObject[vv]):
            MyVar.ancillary_variables = AttrObject[vv]['ancillary_variables']
	
	# Provide the data to the variable - depending on howmany dimensions there are
        if (len(AttrObject[vv]['var_dims']) == 1):
            MyVar[:] = DataObject[vv]
	    
        if (len(AttrObject[vv]['var_dims']) == 2):
            MyVar[:,:] = DataObject[vv]
	    
        if (len(AttrObject[vv]['var_dims']) == 3):
            MyVar[:,:,:] = DataObject[vv]
	    
	    
    ncfw.close()
   
    return # WriteNCCF

    
#************************************************************************
