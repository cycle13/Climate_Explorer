2231	 #!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 8 January 2016
# Last update: 8 January 2016
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This is program to read in data, read in appropriate dictionaries and output a new netCDF file 
# ready for ingestion into the CEDA/ESGF archive. It should be CF compliant to at least CF-1.0. Standard names are used where available and
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
# from scipy.io import netcdf
# import pdb # pdb.set_trace() or c
#
# Kate's:
# from ReadNetCDF import GetGrid - written by Kate Willett, reads in any netCDF grid, can cope with multiple fields
# from WriteNetCDF_CEDAESGF_JAN2016 import WriteNCCF - written by Kate Willett, outputs in CEDA/ESGF format
#
# -----------------------
# DATA
# -----------------------
# InPath: A filepath string linking to the netCDF file to read in
# InFile: A filename string of the netCDF file to read in
# OutPath: A filepath string linking to the netCDF file to output
# OutFile: A filename string of the netCDF file to write out
# InFilename: A string containing the filepath and filename of the input file
# OutFilename: A string containing the filepath and filename of the output file
# A netCDF file containing arrays of data which are then read into numpy arrays:
# 			data, latitude and longtitude info, mdi info etc:
# Dates: a four element intger dictionary of {'StYr': integer e.g., 1973,
#                                            'StMon': integer between 0 and 11,
#                                            'EdYr': integer e.g., 2014,
# 					                         'EdMon': integer between 0 and 11}
# Latitudes: float array of latitude gridbox centres from South to North
# Longitudes: float array of longitude gridbox centres from West to East
# ClimPoints: two element list of start year and end year of climatology period e.g., (1976, 2005)
# DataObject: a list of numpy data arrays to be written out
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
# GlobAttrObject: a list of two element lists of global attributes
#	
#
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# At line: ??? set 'MyChoice' to the variable name to set up the info and dictionaries
# you wish to read, then save and run the code:
# python2.7 Convert_CEDAESGF_JAN2016.py
# 
# The code reads all of the info and dictionaries for the set variable, reads in the data and 
# stores as a list of numpy data arrays. It then calls the write function to output the file:
#
# WriteNCCF(OutFilename,Dates,Latitudes,Longitudes,ClimPoints,DataObject,DimObject,AttrObject)
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
from scipy.io import netcdf
import pdb # pdb.set_trace() or c

from ReadNetCDF import GetGrid # written by Kate Willett, reads in any netCDF grid, can cope with multiple fields
from WriteNetCDF_CEDAESGF_JAN2016 import WriteNCCF

# Set up hardwired variables

#************************************************************************
MyChoice='q' # Choose your choice of dictionary: q, rh, e, t, td, tw, dpd

#-----------------------------------------------------------------------
# Set up for q
if (MyChoice == 'q'):
    # Files
    InPath='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
    InFile='HadISDH.landq.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf.nc'
    OutPath='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
    OutFile='huss_HadISDH_HadOBS_19730101-20141231_v2-0-1-2014p.nc'    
	
    # Dates
    StYr=1973
    StMon=1
    EdYr=2014
    EdMon=12
    ClimPoints=(1976,2005)
	
    # Space Dimensions
    LatInfo=[36,-87.5] # list of gridbox centres calc by ReadNetCDF prog
    LonInfo=[72,-177.5] # list of gridbox centres calc by ReadNetCDF prog
	
    # Data Object List - list of variable names to read in from the netCDF file
    DataObjectList=['q_anoms',
		    'q_abs',
		    'q_std',
		    'q_clims',
		    'q_combinederr',
		    'q_samplingerr',
		    'q_stationerr',
		    'q_obserr',
		    'q_climerr',
		    'q_adjerr',
		    'mean_n_stations',
	            'actual_n_stations',
		    'q_rbar',
		    'q_sbarSQ']
					
    # DimObject list
    DimObjectList=[['time','month','characters','latitude','longitude','bound_pairs'],
                   [((EdYr-StYr)*12)+EdMon,12,10,36,72,2],
	           dict([('var_type','f4'),
		   	 ('var_name','time'),
		   	 ('var_dims',('time',)),
		   	 ('standard_name','time'),
		   	 ('long_name','time'),
		   	 ('units','days since 1973-1-1 00:00:00'),
		   	 ('axis','T'),
		   	 ('calendar','gregorian'),
		   	 ('start_year',1973),
		   	 ('end_year',2014),
		   	 ('start_month',1),
		   	 ('end_month',12),
		   	 ('bounds','bounds_time')]),
		   dict([('var_type','i4'),
		   	 ('var_name','bounds_time'),
		   	 ('var_dims',('time','bound_pairs',)), 
		   	 ('standard_name','time'),
		   	 ('long_name','time period boundaries')]),
		   dict([('var_type','S1'),
		   	 ('var_name','month'),
		   	 ('var_dims',('month','characters',)), 
		   	 ('long_name','month of year')]),
		   dict([('var_type','S1'),
		   	 ('var_name','climbounds'),
		   	 ('var_dims',('month','bound_pairs','characters')), 
		   	 ('long_name','climatology period boundaries')]),
	           dict([('var_type','f4'),
		   	 ('var_name','lat'),
		   	 ('var_dims',('latitude',)), 
		   	 ('standard_name','latitude'),
		   	 ('long_name','gridbox centre latitude'),
		   	 ('units','degrees_north'),
		   	 ('point_spacing','even'),
		   	 ('bounds','bounds_lat')]),
		   dict([('var_type','f4'),
		   	 ('var_name','bounds_lat'),
		   	 ('var_dims',('latitude','bound_pairs',)), 
		   	 ('standard_name','latitude'),
		   	 ('long_name','latitude gridbox boundaries')]),
	           dict([('var_type','f4'),
		   	 ('var_name','lon'),
		   	 ('var_dims',('longitude',)), 
		   	 ('standard_name','longitude'),
		   	 ('long_name','gridbox centre longitude'),
		   	 ('units','degrees_east'),
		   	 ('point_spacing','even'),
		   	 ('bounds','bounds_lon')]),
		   dict([('var_type','f4'),
		   	 ('var_name','bounds_lon'),
		   	 ('var_dims',('longitude','bound_pairs',)), 
		   	 ('standard_name','longitude'),
		   	 ('long_name','longitude gridbox boundaries')])]
	
    # AttrObject list
    AttrObjectList=[dict([('var_type','i4'),
			  ('var_name','hussa'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~2m) specific humidity anomaly'),
			  ('cell_methods','time: mean (interval: 1 month comment: anomaly from climatology) area: mean where land (stations within gridbox)'),
			  ('comment','gridbox mean monthly mean climate anomaly from stations'),
			  ('units','g/kg'),
                          ('add_offset',-100.0), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999),
			  ('reference_period','1976, 2005'),
			  ('ancillary_variables','stdunc sampunc stnunc measunc adjunc climunc')]),
	            dict([('var_type','i4'),
			  ('var_name','huss'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('standard_name','specific_humidity'),
			  ('long_name','near surface (~2m) specific humidity'),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where land (stations within gridbox)'),
			  ('comment','gridbox mean monthly mean from stations'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999),
			  ('ancillary_variables','stdunc sampunc stnunc measunc adjunc climunc')]),
	            dict([('var_type','i4'),
			  ('var_name','std'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~2m) specific humidity standard deviation'),
			  ('cell_methods','time: mean (interval: 1 month) area: variance where land (stations within gridbox)'),
			  ('comment','gridbox standard deviation of monthly mean climate anomaly from stations'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','clm'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','near surface (~2m) specific humidity climatology'),
			  ('cell_methods','time: mean (interval: 1 month comment: over 30 year climatology period) area: mean where land (stations within gridbox)'),
			  ('comment','gridbox mean of monthly mean from stations'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999),
			  ('reference_period','1976, 2005')]),
	            dict([('var_type','i4'),
			  ('var_name','stdunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated combined 1 sigma uncertainty for gridbox'),
			  ('comment','gridbox mean monthly station uncertainty and gridbox sampling uncertainty combined in quadrature assumed uncorrelated'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','sampunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 1 sigma sampling uncertainty for gridbox'),
			  ('cell_methods','area: mean where land (stations within gridbox)'),
			  ('comment','gridbox sampling uncertainty (Jones et al 1997) based on spatio-temporal station presence and intersite correlation assumed uncorrelated'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','stnunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 1 sigma station uncertainty for gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where land (stations within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly measurement, adjustment and climatology uncertainty combined in quadrature for each station and then in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','measunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 1 sigma measurement uncertainty for gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where land (stations within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly measurement uncertainty for each station combined in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','climunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 1 sigma climatology uncertainty for gridbox'),
			  ('cell_methods','area: mean where land (stations within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly climatology uncertainty for each station combined in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','adjunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 1 sigma adjustment uncertainty for gridbox'),
			  ('cell_methods','area: mean where land (stations within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly adjustment (applied and missed) uncertainty for each station combined in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999),
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','meanstncount'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean number of stations within gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: sum where land (stations within gridbox)'),
			  ('units','1'),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','stncount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','actual number of stations within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where land (stations within gridbox)'),
			  ('units','1'),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','rbar'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','intersite correlation (rbar)'),
			  ('comment','intersite correlation for each gridbox following Jones et al 1997 (rbar)'),
			  ('units','1'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-999)]),
	            dict([('var_type','i4'),
			  ('var_name','sbar2'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean gridbox variance (sbar2)'),
			  ('comment','mean variance over all stations in gridbox following Jones et al 1997 (rbar)'),
			  ('units','g/kg'),
			  ('add_offset',-100.), # storedval=int((var-offset)/scale)
			  ('scale_factor',0.01),
			  ('valid_min',0),
			  ('valid_max',30000), # max integer
			  ('missing_value',-999), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-999)])]  

    # Set up Global Attribute List of Lists
    GlobAttrObjectList=dict([['File_created',datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')], # Is there a call for time stamping?
			     ['Description','HadISDH monthly mean land surface specific humidity climate monitoring product \
from 1973 onwards. Stations have been quality controlled, homogenised and averaged over 5deg by 5deg gridboxes (no smoothing \
or interpolation). Gridbox 1 sigma uncertainty estimates are provided that represent spatio-temporal sampling within the \
gridbox and combined station level uncertainties for measurement, climatology and homogenisation (applied and missed adjustments).'],
			     ['Title','HadISDH monthly mean land surface specific humidity climate monitoring product'], 
			     ['Institution','Met Office Hadley Centre (UK), \
National Centres for Environmental Information (USA), \
National Physical Laboratory (UK), \
Climatic Research Unit (UK), \
Maynooth University (Rep. of Ireland)'],
			     ['History','See Willett et al., (2014) REFERENCE for more information. \
See www.metoffice.gov.uk/hadobs/hadisdh/ for more information and related data and figures. \
See also http://catalogue.ceda.ac.uk/uuid/251474c7b09449d8b9e7aeaf1461858f for more information and related data.\
Follow @metofficeHadOBS to keep up to date with Met Office Hadley Centre HadOBS dataset developements. \
See hadisdh.blogspot.co.uk for HadISDH updates, bug fixes and explorations.'], 
			     ['Licence','HadISDH is distributed under the Non-Commercial Government Licence: \
http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/non-commercial-government-licence.htm. \
The data are freely available for any non-comercial use with attribution to the data providers. Please cite \
Willett et al.,(2014) with a link to the REFERENCE provided in the REFERENCE attribute.'],
			     ['Project', 'HadOBS: Met Office Hadley Centre Climate Monitoring Data-product www.metoffice.gov.uk/hadobs'],
			     ['Processing_level','Hourly station data selected for length and continuity, quality controlled, averaged \
to monthly means, adjusted to remove inhomogeneity and then averaged over 5deg by 5deg gridboxes.'],
			     ['Acknowledgment','Kate Willett, Robert Dunn and David Parker were supported by the Joint DECC/Defra \
Met Office Hadley Centre Climate Programme (GA01101). Phil D. Jones has been supported by the Office of Science (BER), US Department \
of Energy, Grant no. DE-SC0005689. Stephanie Bell and Michael de Podesta were supported by the UK National Measurement System Programme \
for Engineering and Flow Metrology, and by the MeteoMet Project of the European Metrology Research Programme. We also thank the internal \
reviewers at NCDC, NPL and the Met Office for their thorough reviews and advice. We thank Adrian Simmons and Andreas Becker for their \
thoughtful and constructive reviews.'],
			     ['Source','HadISD.1.0.3.2014f (Dunn et al. 2012) from the National Centres for Environmental Information \
International Surface Database (ISD): www.ncdc.noaa.gov/isd'],
			     ['Comment',''],
			     ['References','Willett, K. M., Dunn, R. J. H., Thorne, P. W., Bell, S., de Podesta, M., Parker, D. E., \
Jones, P. D. and Williams, Jr., C. N.: HadISDH land surface multi-variable humidity and temperature record for climate monitoring, \
Clim. Past, 10, 1983-2006, doi:10.5194/cp-10-1983-2014, 2014'],
			     ['Creator_name','Kate Willett'],
			     ['Creator_email','kate.willett@metoffice.gov.uk'],
			     ['Version','v2.0.1.2014p'],
			     ['doi',''], # This needs to be filled in
			     ['Conventions','CF-1.6'],
			     ['netCDF_type','NETCDF4_CLASSIC']]) # May now actually be higher than this - can we say NetCDF4 compliant?
						  
						  
#------------------------------------------------------------------------	
# Set up for rh
#if (MyChoice == 'rh'):	

#------------------------------------------------------------------------	
# Set up for e
#if (MyChoice == 'e'):	

#------------------------------------------------------------------------	
# Set up for td
#if (MyChoice == 'td'):	

#------------------------------------------------------------------------	
# Set up for tw
#if (MyChoice == 'tw'):	

#------------------------------------------------------------------------	
# Set up for t
#if (MyChoice == 't'):	

#------------------------------------------------------------------------	
# Set up for dpd
#if (MyChoice == 'dpd'):	


#************************************************************************
# Main Program
#************************************************************************
# Get filenames
InFileName = InPath+InFile
OutFileName = OutPath+OutFile

# Set up dates dictionary
Dates = dict([('StYr',StYr),('StMon',StMon),('EdYr',EdYr),('EdMon',EdMon)])

# Read in data objects to make a list of numpy arrays DataObject and also supply Latitudes and Longitudes
DataObject,Latitudes,Longitudes = GetGrid(InFileName,DataObjectList,LatInfo,LonInfo)

# Call writing module
WriteNCCF(OutFileName,Dates,Latitudes,Longitudes,ClimPoints,DataObject,DimObjectList,AttrObjectList,GlobAttrObjectList)

# Finish
print("And we are done!")
