#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 8 January 2016
# Last update: 6 March 2020
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
#>module load scitools/default-current
# python3 Convert_CEDAESGF_JAN2016.py
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
# Version 2 6th Mar 2020
# ---------
#  
# Enhancements
# Can now work with marine data too - not completed yet
#  
# Changes
# Now Python 3
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
from scipy.io import netcdf
import pdb # pdb.set_trace() or c

from ReadNetCDF import GetGrid4 # written by Kate Willett, reads in any netCDF grid, can cope with multiple fields
from WriteNetCDF_CEDAESGF_JAN2016 import WriteNCCF

# Set up hardwired variables

#************************************************************************
MyChoice='q' # Choose your choice of dictionary: q, rh, e, t, td, tw, dpd

Domain = 'land'
#Domain = 'marine'
#Domain = 'blend'

# Generic Things
# VERSION
if (Domain == 'land'):
    
    version = 'v4-2-0-2019f'
    verstring = '4.2.0.2019f'
    DomainHGT = '2'
    DomainOBS = 'stations'

elif (Domain == 'marine'):

    version = 'v1-0-0-2019f'
    verstring = '1.0.0.2019f'
#    version = 'v1-0-0-2018f'
#    verstring = '1.0.0.2018f'
    DomainHGT = '10'
    DomainOBS = 'ships'

elif (Domain == 'blend'):

    version = 'v1-0-0-2019f'
    verstring = '1.0.0.2019f'
    DomainHGT = '2/10'
    DomainOBS = 'observations'

#nowmon = 'FEB'
#nowyear = '2019'
nowmon = 'JAN'
nowyear = '2020'

# Dates
StYr=1973
StMon=1
EdYr=2019
EdMon=12
ClimPoints=(1981,2010)
#ClimPoints=(1976,2005)
climbo = str(ClimPoints[0])[2:4]+str(ClimPoints[1])[2:4]

# Files
InPath='/data/users/hadkw/WORKING_HADISDH/UPDATE'+str(EdYr)+'/STATISTICS/GRIDS/'
OutPath='/data/users/hadkw/WORKING_HADISDH/UPDATE'+str(EdYr)+'/STATISTICS/GRIDS/'

# Space Dimensions
LatInfo=[36,-87.5] # list of gridbox centres calc by ReadNetCDF prog
LonInfo=[72,-177.5] # list of gridbox centres calc by ReadNetCDF prog

# Set up Global Attributes for  List of Lists
if (Domain == 'land'):
    Description = 'climate monitoring product \
               from 1973 onwards. Stations have been quality controlled, homogenised and averaged over 5deg by 5deg gridboxes (no smoothing \
               or interpolation). Gridbox 2 sigma uncertainty estimates are provided that represent spatio-temporal sampling within the \
               gridbox and combined station level uncertainties for measurement, climatology and homogenisation (applied and missed adjustments).'
    Institution = 'Met Office Hadley Centre (UK), \
               National Centres for Environmental Information (USA), \
               National Physical Laboratory (UK), \
               Climatic Research Unit (UK), \
               Maynooth University (Rep. of Ireland)'
    History = 'See Willett et al., (2014) REFERENCE for more information. \
           See www.metoffice.gov.uk/hadobs/hadisdh/ for more information and related data and figures. \
           Follow @metofficeHadOBS to keep up to date with Met Office Hadley Centre HadOBS dataset developements. \
           See hadisdh.blogspot.co.uk for HadISDH updates, bug fixes and explorations.'
# Non-commercial license
#    Licence = 'HadISDH is distributed under the Non-Commercial Government Licence: \
#           http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/non-commercial-government-licence.htm. \
#           The data are freely available for any non-commercial use with attribution to the data providers. Please cite \
#           Willett et al.,(2014) and Smith et al., (2011) with a link to the REFERENCES provided in the REFERENCE attribute.'
# Open Government License
    Licence = 'HadISDH is distributed under the Open Government Licence: \
               http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/. \
               The data are freely available for use with attribution to the data providers. Please cite \
               references provided in the REFERENCE attribute.'
    Project = 'HadOBS: Met Office Hadley Centre Climate Monitoring Data-product www.metoffice.gov.uk/hadobs'
    Processing_level = 'Hourly station data selected for length and continuity, quality controlled, averaged \
                    to monthly means, adjusted to remove inhomogeneity and then averaged over 5deg by 5deg gridboxes.'
    Acknowledgement = 'Kate Willett, Robert Dunn and David Parker were supported by the Joint UK BEIS/Defra \
                  Met Office Hadley Centre Climate Programme (GA01101). Phil D. Jones has been supported by the Office of Science (BER), US Department \
                  of Energy, Grant no. DE-SC0005689. Stephanie Bell and Michael de Podesta were supported by the UK National Measurement System Programme \
                  for Engineering and Flow Metrology, and by the MeteoMet Project of the European Metrology Research Programme. We also thank the internal \
                  reviewers at NCDC, NPL and the Met Office for their thorough reviews and advice. We thank Adrian Simmons and Andreas Becker for their \
                  thoughtful and constructive reviews.'
# THIS BIT NEEDS EDITING EVERY YEAR **************************
    Source = 'HadISD.3.1.0.2019f (Dunn et al. 2016) from the National Centres for Environmental Information \
          International Surface Database (ISD) (Smith et al., 2011) from www.ncdc.noaa.gov/isd'
    Comment = ''
    References = 'Willett, K. M., Dunn, R. J. H., Thorne, P. W., Bell, S., de Podesta, M., Parker, D. E., \
              Jones, P. D. and Williams, Jr., C. N.: HadISDH land surface multi-variable humidity and temperature record for climate monitoring, \
              Clim. Past, 10, 1983-2006, doi:10.5194/cp-10-1983-2014, 2014, \
	      Dunn, R. J. H., et al. 2016: Expanding HadISD: quality-controlled, sub-daily \
	      station data from 1931, Geoscientific Instrumentation, Methods and Data Systems, 5, 473-491. \
	      Smith, A., N. Lott, and R. Vose, 2011: The Integrated Surface Database: Recent \
	      Developments and Partnerships. Bulletin of the American Meteorological Society, \
	      92, 704-708, doi:10.1175/2011BAMS3015.1'
    Creator_name = 'Kate Willett'
    Creator_email = 'kate.willett@metoffice.gov.uk'
    Conventions = 'CF-1.6'
    netCDF_type = 'NETCDF4_CLASSIC'

elif (Domain == 'marine'):
    Description = 'climate monitoring product \
               from 1973 onwards. Ship observations have been quality controlled, bias adjusted (height to 10m, unaspirated screens) and averaged over 5deg by 5deg gridboxes (no smoothing \
               or interpolation). Gridbox 2 sigma uncertainty estimates are provided that represent spatio-temporal sampling within the \
               gridbox and combined observation level uncertainties for measurement, climatology, whole number reporting and bias adjustment.'
    Institution = 'Met Office Hadley Centre (UK), \
               National Oceanography Centre Southampton (UK)'
    History = 'See Willett et al., (2020) REFERENCE for more information. \
           See www.metoffice.gov.uk/hadobs/hadisdh/ for more information and related data and figures. \
           Follow @metofficeHadOBS to keep up to date with Met Office Hadley Centre HadOBS dataset developements. \
           See hadisdh.blogspot.co.uk for HadISDH updates, bug fixes and explorations.'
## Non-commercial license
#    Licence = 'HadISDH is distributed under the Non-Commercial Government Licence: \
#           http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/non-commercial-government-licence.htm. \
#           The data are freely available for any non-comercial use with attribution to the data providers. Please cite \
#           Willett et al.,(2020) and Freeman et al., (2017) with a link to the REFERENCES provided in the REFERENCE attribute.'
# Open Government License
    Licence = 'HadISDH is distributed under the Open Government Licence: \
               http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/. \
               The data are freely available for use with attribution to the data providers. Please cite \
               references provided in the REFERENCE attribute.'
    Project = 'HadOBS: Met Office Hadley Centre Climate Monitoring Data-product www.metoffice.gov.uk/hadobs'
    Processing_level = 'Hourly ship observation data quality controlled, bias adjusted, averaged \
                    to monthly means, then averaged over 5deg by 5deg gridboxes.'
    Acknowledgement = 'Kate Willett, Robert Dunn and John Kennedy were supported by the Joint UK BEIS/Defra \
                  Met Office Hadley Centre Climate Programme (GA01101).'
# THIS BIT NEEDS EDITING EVERY YEAR **************************
    Source = 'ICOADS.3.0.0 (1973-2014) and ICOADS.3.0.1 (2015 onwards) (Freeman et al. 2017) from icoads.noaa.gov.'
    Comment = ''
    References = 'Willett, K. M.,, R. J. H. Dunn, J. J. Kennedy, and D. I. Berry, in review.: Development of the HadISDH marine \
              humidity climate monitoring dataset. Earth System Science Datasets. \
	      Freeman, E., S.D. Woodruff, S.J. Worley, S.J. Lubker, E.C. Kent, W.E. Angel, D.I . Berry, P. Brohan, R. Eastman, L. Gates, W. \
	      Gloeden, Z. Ji, J. Lawrimore, N.A. Rayner, G. Rosenhagen, and S.R. Smith, 2017: ICOADS Release 3.0: A major update to the \
	      historical marine climate record. Int. J. Climatol. (CLIMAR-IV Special Issue), 37, 2211-2237 (doi:10.1002/joc.4775).'
    Creator_name = 'Kate Willett'
    Creator_email = 'kate.willett@metoffice.gov.uk'
    Conventions = 'CF-1.6'
    netCDF_type = 'NETCDF4_CLASSIC'

elif (Domain == 'blend'):
    Description = 'climate monitoring product \
               from 1973 onwards. Ship and station observations have been quality controlled, buddy checked (marine), homogenised (land), bias adjusted (marine - height to 10m, unaspirated screens) \
	       and averaged over 5deg by 5deg gridboxes (no smoothing \
               or interpolation). Gridbox 2 sigma uncertainty estimates are provided that represent spatio-temporal sampling within the \
               gridbox and combined observation level uncertainties.'
    Institution = 'Met Office Hadley Centre (UK), \
               National Centres for Environmental Information (USA), \
               National Physical Laboratory (UK), \
               Climatic Research Unit (UK), \
               Maynooth University (Rep. of Ireland) \
               National Oceanography Centre Southampton (UK)'
    History = 'See Willett et al., (2014, 2020) REFERENCE for more information. \
           See www.metoffice.gov.uk/hadobs/hadisdh/ for more information and related data and figures. \
           Follow @metofficeHadOBS to keep up to date with Met Office Hadley Centre HadOBS dataset developements. \
           See hadisdh.blogspot.co.uk for HadISDH updates, bug fixes and explorations.'
## Non-commercial license
#    Licence = 'HadISDH is distributed under the Non-Commercial Government Licence: \
#           http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/non-commercial-government-licence.htm. \
#           The data are freely available for any non-comercial use with attribution to the data providers. Please cite \
#           Willett et al.,(2020) and Freeman et al., (2017) with a link to the REFERENCES provided in the REFERENCE attribute.'
# Open Government License
    Licence = 'HadISDH is distributed under the Open Government Licence: \
               http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/. \
               The data are freely available for use with attribution to the data providers. Please cite \
               references provided in the REFERENCE attribute.'
    Project = 'HadOBS: Met Office Hadley Centre Climate Monitoring Data-product www.metoffice.gov.uk/hadobs'
    Processing_level = 'Hourly ship and station observation data quality controlled, buddy checked (marine), \
                       homogenised (land), bias adjusted (marine), averaged \
                    to monthly means, and averaged over 5deg by 5deg gridboxes.'
    Acknowledgement = 'Kate Willett, Robert Dunn, John Kennedy and David Parker were supported by the Joint UK BEIS/Defra \
                  Met Office Hadley Centre Climate Programme (GA01101).'
# THIS BIT NEEDS EDITING EVERY YEAR **************************
    Source = 'ICOADS.3.1.0 (1973-2014) and ICOADS.3.0.1 (2015 onwards) (Freeman et al. 2017) from icoads.noaa.gov, \
          HadISD.3.1.2018f (Dunn et al., 2016), International Surface Database (ISD) (Smith et al, 2011) from www.ncdc.noaa.gov/isd'
    Comment = ''
    References = 'Willett, K. M.,, R. J. H. Dunn, J. J. Kennedy, and D. I. Berry, in review.: Development of the HadISDH marine \
              humidity climate monitoring dataset. Earth System Science Datasets. \
	      Freeman, E., S.D. Woodruff, S.J. Worley, S.J. Lubker, E.C. Kent, W.E. Angel, D.I . Berry, P. Brohan, R. Eastman, L. Gates, W. \
	      Gloeden, Z. Ji, J. Lawrimore, N.A. Rayner, G. Rosenhagen, and S.R. Smith, 2017: ICOADS Release 3.0: A major update to the \
	      historical marine climate record. Int. J. Climatol. (CLIMAR-IV Special Issue), 37, 2211-2237 (doi:10.1002/joc.4775). \
	      Willett, K. M., Dunn, R. J. H., Thorne, P. W., Bell, S., de Podesta, M., Parker, D. E., \
              Jones, P. D. and Williams, Jr., C. N.: HadISDH land surface multi-variable humidity and temperature record for climate monitoring, \
              Clim. Past, 10, 1983-2006, doi:10.5194/cp-10-1983-2014, 2014, \
	      Dunn, R. J. H., et al. 2016: Expanding HadISD: quality-controlled, sub-daily \
	      station data from 1931, Geoscientific Instrumentation, Methods and Data Systems, 5, 473-491. \
	      Smith, A., N. Lott, and R. Vose, 2011: The Integrated Surface Database: Recent \
	      Developments and Partnerships. Bulletin of the American Meteorological Society, \
	      92, 704-708, doi:10.1175/2011BAMS3015.1'
    Creator_name = 'Kate Willett'
    Creator_email = 'kate.willett@metoffice.gov.uk'
    Conventions = 'CF-1.6'
    netCDF_type = 'NETCDF4_CLASSIC'

# Variable IDS and Names DIct
NameList = dict([('q',['q','g/kg','specific_humidity','specific humidity','hussa','huss']),
                ('rh',['RH','%rh','relative_humidity','relative humidity','hursa','hurs']),
		('e',['e','hPa','water_vapor_partial_pressure_in_air','vapor pressure','vpsa','vps']),
		('t',['T','deg C','air_temperature','air temperature','tasa','tas']),
		('tw',['Tw','deg C','wet_bulb_temperature','wetbulb temperature','twsa','tws']),
		('td',['Td','deg C','dew_point_temperature','dewpoint temperature','tdsa','tds']),
		('dpd',['DPD','deg C','dew_point_depression','dewpoint depression','dpdsa','dpds'])])

# Data Object List - list of variable names to read in from the netCDF file
if (Domain == 'land'):
    DataObjectList=[MyChoice+'_anoms',
		MyChoice+'_abs',
		MyChoice+'_std',
		MyChoice+'_clims',
		MyChoice+'_combinederr', # (2 sigma)
		MyChoice+'_samplingerr', # (2 sigma)
		MyChoice+'_stationerr', # (2 sigma)
		MyChoice+'_obserr', # (2 sigma)
		MyChoice+'_climerr',
		MyChoice+'_adjerr',
		'mean_n_stations',
	        'actual_n_stations',
		MyChoice+'_rbar',
		MyChoice+'_sbarSQ']

elif (Domain == 'marine'):						
    DataObjectList=[MyChoice+'_anoms',
		MyChoice+'_abs',
		#MyChoice+'_std',
		MyChoice+'_clims',
		MyChoice+'_clim_std', # monthly climatological standard deviation,
		MyChoice+'_abs_uHGT', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uHGT', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uSCN', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uSCN', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uC', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uC', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uR', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uR', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uM', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uM', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uOBS', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uOBS', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uSAMP', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uSAMP', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_uSAMP_n_grids', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_usbarSQ', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_usbarSQ', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_usbarSQ_n_grids', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_urbar', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_urbar', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_abs_uFULL', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_anoms_uFULL', # monthyl mean <var> adjustment height uncertainty (2 sigma)
		MyChoice+'_n_grids', # number of 1by1 daily grids within gridbox
	        MyChoice+'_n_obs', # number of observations within gridbox
		MyChoice+'_clims_n_grids', # number of 1by1 daily grids within gridbox climatology
	        MyChoice+'_clims_n_obs', # number of observations within gridbox climatology
		MyChoice+'_clim_std_n_grids', # number of 1by1 daily grids within gridbox climatological standard deviation
	        MyChoice+'_clim_std_n_obs'] # number of observations within gridbox climatological standard deviation

elif (Domain == 'blend'):
    DataObjectList=[MyChoice+'_anoms',
		MyChoice+'_abs',
		MyChoice+'_clims',
		MyChoice+'_land_std',
		MyChoice+'_marine_clims_std',
		MyChoice+'_anoms_obserr', # (2 sigma)
		MyChoice+'_abs_obserr', # (2 sigma)
		MyChoice+'_anoms_samplingerr', # (2 sigma)
		MyChoice+'_abs_samplingerr', # (2 sigma)
		'marine_actual_pseudo_stations',
		'marine_mean_pseudo_stations',
		MyChoice+'_anoms_combinederr', # (2 sigma)
		MyChoice+'_abs_combinederr', # (2 sigma)
		'land_mean_n_stations',
	        'land_actual_n_stations',
		'marine_n_grids',
		'marine_n_obs',
		'marine_clim_n_grids',
		'marine_clim_n_obs',
		'marine_climstd_n_grids',
		'marine_climstd_n_obs']
    				    
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
    		     ('start_year',StYr),
    		     ('end_year',EdYr),
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
    		     ('var_name','latitude'),
    		     ('var_dims',('latitude',)), 
    		     ('standard_name','latitude'),
    		     ('long_name','gridbox centre latitude'),
    		     ('units','degrees_north'),
    		     ('axis','Y'),
    		     ('point_spacing','even'),
    		     ('bounds','bounds_lat')]),
    	       dict([('var_type','f4'),
    		     ('var_name','bounds_lat'),
    		     ('var_dims',('latitude','bound_pairs',)), 
    		     ('standard_name','latitude'),
    		     ('long_name','latitude gridbox boundaries')]),
    	       dict([('var_type','f4'),
    		     ('var_name','longitude'),
    		     ('var_dims',('longitude',)), 
    		     ('standard_name','longitude'),
    		     ('long_name','gridbox centre longitude'),
    		     ('units','degrees_east'),
    		     ('axis','X'),
    		     ('point_spacing','even'),
    		     ('bounds','bounds_lon')]),
    	       dict([('var_type','f4'),
    		     ('var_name','bounds_lon'),
    		     ('var_dims',('longitude','bound_pairs',)), 
    		     ('standard_name','longitude'),
    		     ('long_name','longitude gridbox boundaries')])]

    # AttrObject list
    # These were all i4 and I was using the offset and scale factor but the output netcdf didn't work with ncview in that the missing
    # data identifier wasn't recognised. I'm hoping its ok to use zlib = True and save as floats rather than faff with scale and offset
    # Its MUCH MUCH smaller at least
    # I've also switched the MDI to -1e30
    # I've also not written out the valid_min and valid_max as these are fairly meaningless
    # I suppose they should be the physical max and min of the variable
if (Domain == 'land'):

    AttrObjectList=[dict([('var_type','float32'),
			  ('var_name',NameList[MyChoice][4]),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' anomaly'),
			  ('cell_methods','time: mean (interval: 1 month comment: anomaly from climatology) area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean monthly mean climate anomaly from '+DomainOBS),
			  ('units', NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1])),		#'1976, 2005'),
			  ('ancillary_variables','stdunc sampunc stnunc measunc adjunc climunc')]),
	            dict([('var_type','float32'),
			  ('var_name',NameList[MyChoice][5]),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('standard_name',NameList[MyChoice][2]),
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean monthly mean from '+DomainOBS),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('ancillary_variables','stdunc sampunc stnunc measunc adjunc climunc')]),
	            dict([('var_type','float32'),
			  ('var_name','std'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' standard deviation'),
			  ('cell_methods','time: mean (interval: 1 month) area: variance where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox standard deviation of monthly mean climate anomaly from '+DomainOBS),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','clm'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' climatology'),
			  ('cell_methods','time: mean (interval: 1 month comment: over 30 year climatology period) area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean of 30 yr climatological monthly mean from '+DomainOBS),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1]))]),		#'1976, 2005')]),
	            dict([('var_type','float32'),
			  ('var_name','stdunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated combined 2 sigma uncertainty for gridbox'),
			  ('comment','gridbox mean monthly station uncertainty and gridbox sampling uncertainty combined in quadrature assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','sampunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma sampling uncertainty for gridbox'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox sampling uncertainty (Jones et al 1997) based on spatio-temporal station presence and intersite correlation assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','stnunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma observation uncertainty for gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where '+Domain+' ('+DomainOBS+' within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly measurement, adjustment and climatology uncertainty combined in quadrature for each '+DomainOBS+' and then in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','measunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma measurement uncertainty for gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where '+Domain+' ('+DomainOBS+' within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly measurement uncertainty for each observation combined in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','climunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma climatology uncertainty for gridbox'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly climatology uncertainty for each observation combined in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','adjunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma adjustment uncertainty for gridbox'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox combined in quadtrature)'),
			  ('comment','gridbox mean monthly adjustment (applied and missed) uncertainty for each observation combined in quadrature over the gridbox assumed to be uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','i4'),
			  ('var_name','meanstncount'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean number of '+DomainOBS+' within gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: sum where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','stncount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','actual number of '+DomainOBS+' within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','float32'),
			  ('var_name','rbar'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','intersite correlation (rbar)'),
			  ('comment','intersite correlation for each gridbox following Jones et al 1997 (rbar)'),
			  ('units','1'),
			  ('missing_value',-1e30), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','sbar2'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean gridbox variance (sbar2)'),
			  ('comment','mean variance over all observations in gridbox following Jones et al 1997 (sbar2)'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1e30)])]  

elif (Domain == 'marine'):
    #marine attributes
    AttrObjectList=[dict([('var_type','float32'),
			  ('var_name',NameList[MyChoice][4]),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' anomaly'),
			  ('cell_methods','time: mean (interval: 1 month comment: anomaly from climatology) area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean monthly mean climate anomaly from '+DomainOBS),
			  ('units', NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1])),		#'1976, 2005'),
			  ('ancillary_variables','abs_stdunc anoms_stdunc abs_sampunc anoms_stdunc abs_obsunc anoms_obsunc')]),
	            dict([('var_type','float32'),
			  ('var_name',NameList[MyChoice][5]),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('standard_name',NameList[MyChoice][2]),
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]),
			  ('cell_methods','time: mean (interval: 1 month) area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean monthly mean from '+DomainOBS),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('ancillary_variables','abs_stdunc anoms_stdunc abs_sampunc anoms_stdunc abs_obsunc anoms_obsunc')]),
#	            dict([('var_type','float32'),
#			  ('var_name','std'),
#			  ('var_dims',('time','latitude','longitude',)), 
#			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' standard deviation'),
#			  ('cell_methods','time: mean (interval: 1 month) area: variance where '+Domain+' ('+DomainOBS+' within gridbox)'),
#			  ('comment','gridbox standard deviation of monthly mean climate anomaly from '+DomainOBS),
#			  ('units',NameList[MyChoice][1]),
#			  ('missing_value',-1e30),
#			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','clm'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' climatology'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox) time: mean (interval: 1 month comment: over 30 year climatology period)'),
			  ('comment','30 year monthly mean of gridbox mean'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1]))]),		#'1976, 2005')]),
	            dict([('var_type','float32'),
			  ('var_name','clmstd'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' climatological standard deviations'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox) time: standard deviation of monthly means (interval: 1 month comment: over 30 year climatology period) '),
			  ('comment','30 year standard deviation of gridbox monthly mean'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1]))]),		#'1976, 2005')]),
	            dict([('var_type','float32'),
			  ('var_name','abs_hgtadjunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','correlated 2 sigma uncertainty for ship height bias adjustments for actual values'),
			  ('comment','gridbox mean monthly ship height bias adjustment uncertainty combined in quadrature assuming correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_hgtadjunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','correlated 2 sigma uncertainty for ship height bias adjustments for anomaly values'),
			  ('comment','gridbox mean monthly ship height bias adjustment uncertainty combined in quadrature assuming correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_instadjunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','correlated 2 sigma uncertainty for instrument bias adjustments for actual values'),
			  ('comment','gridbox mean monthly instrument bias adjustment uncertainty combined in quadrature assuming correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_instadjunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','correlated 2 sigma uncertainty for instrument bias adjustments for anomaly values'),
			  ('comment','gridbox mean monthly instrument bias adjustment uncertainty combined in quadrature assuming correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_clmunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','correlated 2 sigma uncertainty for climatology for actual values'),
			  ('comment','gridbox mean monthly climatology uncertainty combined in quadrature assuming correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_clmunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','correlated 2 sigma uncertainty for climatology for anomaly values'),
			  ('comment','gridbox mean monthly climatology uncertainty combined in quadrature assuming correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_wholeunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma uncertainty for whole number reporting for actual values'),
			  ('comment','gridbox mean monthly whole number uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_wholeunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma uncertainty for whole number reporting for anomaly values'),
			  ('comment','gridbox mean monthly whole number reporting uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_measunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma uncertainty for measurement for actual values'),
			  ('comment','gridbox mean monthly measurement uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_measunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma uncertainty for measurement for anomaly values'),
			  ('comment','gridbox mean monthly measurement uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_obsunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma combined observation uncertainty for actual values'),
			  ('comment','gridbox mean monthly combined observation uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_obsunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma combined observations uncertainty for anomaly values'),
			  ('comment','gridbox mean monthly combined observations uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_sampunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma sampling uncertainty for gridbox actual values'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox sampling uncertainty (Jones et al 1997) based on spatio-temporal station presence and intersite correlation assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_sampunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma sampling uncertainty for gridbox anomaly values'),
			  ('cell_methods','area: mean where '+Domain+' ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox sampling uncertainty (Jones et al 1997) based on spatio-temporal station presence and intersite correlation assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','i4'),
			  ('var_name','pseudostncount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','number of pseudo stations within gridbox'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_sbarsq'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','gridbox mean pseudo-station variance (sbarSQ for sampling uncertainty) for gridbox actual values'),
			  ('comment','mean variance over all observations in gridbox following Jones et al 1997 (sbarSQ)'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_sbarsq'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','gridbox mean pseudo-station variance (sbarSQ for sampling uncertainty) for gridbox anomaly values'),
			  ('comment','mean variance over all observations in gridbox following Jones et al 1997 (sbarSQ)'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','i4'),
			  ('var_name','meanpseudostncount'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean number of pseudo stations within gridbox'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_rbar'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','intersite correlation (rbar) for actual values'),
			  ('comment','intersite correlation for each gridbox following Jones et al 1997 (rbar)'),
			  ('units','1'),
			  ('missing_value',-1e30), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_rbar'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','intersite correlation (rbar) for anomalies'),
			  ('comment','intersite correlation for each gridbox following Jones et al 1997 (rbar)'),
			  ('units','1'),
			  ('missing_value',-1e30), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_stdunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated combined 2 sigma uncertainty for gridbox actual values'),
			  ('comment','gridbox mean monthly observation uncertainty and gridbox sampling uncertainty combined in quadrature assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_stdunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated combined 2 sigma uncertainty for gridbox anomalu values'),
			  ('comment','gridbox mean monthly observation uncertainty and gridbox sampling uncertainty combined in quadrature assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','i4'),
			  ('var_name','gridcount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','number of 1by1 daily grids within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where 1by1 daily grids within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','obscount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','number of observations within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where observations within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','clmgridcount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of 1by1 daily grids within gridbox climatology'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where 1by1 daily grids within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','clmobscount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of observations within gridbox climatology'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where observations within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','clmstdgridcount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of 1by1 daily grids within gridbox climatological standard deviations'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where 1by1 daily grids within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','clmstdobscount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of observations within gridbox climatological standard deviations'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where observations within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)])]  

elif (Domain == 'blend'):
    #marine attributes    				    
    AttrObjectList=[dict([('var_type','float32'),
			  ('var_name',NameList[MyChoice][4]),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' anomaly'),
			  ('cell_methods','time: mean (interval: 1 month comment: anomaly from climatology) area: mean ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean monthly mean climate anomaly from '+DomainOBS),
			  ('units', NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1])),		#'1976, 2005'),
			  ('ancillary_variables','abs_stdunc anoms_stdunc abs_sampunc anoms_stdunc abs_obsunc anoms_obsunc')]),
	            dict([('var_type','float32'),
			  ('var_name',NameList[MyChoice][5]),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('standard_name',NameList[MyChoice][2]),
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]),
			  ('cell_methods','time: mean (interval: 1 month) area: mean ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox mean monthly mean from '+DomainOBS),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('ancillary_variables','abs_stdunc anoms_stdunc abs_sampunc anoms_stdunc abs_obsunc anoms_obsunc')]),
	            dict([('var_type','float32'),
			  ('var_name','clm'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' climatology'),
			  ('cell_methods','area: mean ('+DomainOBS+' within gridbox) time: mean (interval: 1 month comment: over 30 year climatology period)'),
			  ('comment','30 year monthly mean of gridbox mean'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1]))]),		#'1976, 2005')]),
	            dict([('var_type','float32'),
			  ('var_name','land_std'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','near surface (~2m) '+NameList[MyChoice][3]+' standard deviation'),
			  ('cell_methods','time: mean (interval: 1 month) area: variance where land ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox standard deviation of monthly mean climate anomaly from '+DomainOBS),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','marine_clmstd'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','near surface (~'+DomainHGT+'m) '+NameList[MyChoice][3]+' climatological standard deviations'),
			  ('cell_methods','area: mean where ocean ('+DomainOBS+' within gridbox) time: standard deviation of monthly means (interval: 1 month comment: over 30 year climatology period) '),
			  ('comment','30 year standard deviation of gridbox monthly mean'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30),
			  ('reference_period',str(ClimPoints[0])+', '+str(ClimPoints[1]))]),		#'1976, 2005')]),
	            dict([('var_type','float32'),
			  ('var_name','abs_obsunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma combined observation uncertainty for actual values'),
			  ('comment','gridbox mean monthly combined observation uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_obsunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma combined observations uncertainty for anomaly values'),
			  ('comment','gridbox mean monthly combined observations uncertainty combined in quadrature assuming no correlation'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_sampunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma sampling uncertainty for gridbox actual values'),
			  ('cell_methods','area: mean ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox sampling uncertainty (Jones et al 1997) based on spatio-temporal station presence and intersite correlation assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_sampunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated 2 sigma sampling uncertainty for gridbox anomaly values'),
			  ('cell_methods','area: mean ('+DomainOBS+' within gridbox)'),
			  ('comment','gridbox sampling uncertainty (Jones et al 1997) based on spatio-temporal station presence and intersite correlation assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_pseudostncount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','number of pseudo stations within gridbox'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_meanpseudostncount'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean number of pseudo stations within gridbox'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','float32'),
			  ('var_name','abs_stdunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated combined 2 sigma uncertainty for gridbox actual values'),
			  ('comment','gridbox mean monthly observation uncertainty and gridbox sampling uncertainty combined in quadrature assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','float32'),
			  ('var_name','anoms_stdunc'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','uncorrelated combined 2 sigma uncertainty for gridbox anomalu values'),
			  ('comment','gridbox mean monthly observation uncertainty and gridbox sampling uncertainty combined in quadrature assumed uncorrelated'),
			  ('units',NameList[MyChoice][1]),
			  ('missing_value',-1e30),
			  ('_FillValue',-1e30)]),
	            dict([('var_type','i4'),
			  ('var_name','land_meanstncount'),
			  ('var_dims',('latitude','longitude',)), 
			  ('long_name','mean number of '+DomainOBS+' within gridbox'),
			  ('cell_methods','time: mean (interval: 1 month) area: sum ('+DomainOBS+' within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','land_stncount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','actual number of '+DomainOBS+' within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum ('+DomainOBS+' within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_gridcount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','number of 1by1 daily grids within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where 1by1 daily grids within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_obscount'),
			  ('var_dims',('time','latitude','longitude',)), 
			  ('long_name','number of observations within gridbox'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where observations within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_clmgridcount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of 1by1 daily grids within gridbox climatology'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where 1by1 daily grids within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_clmobscount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of observations within gridbox climatology'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where observations within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_clmstdgridcount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of 1by1 daily grids within gridbox climatological standard deviations'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where 1by1 daily grids within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)]),
	            dict([('var_type','i4'),
			  ('var_name','marine_clmstdobscount'),
			  ('var_dims',('month','latitude','longitude',)), 
			  ('long_name','number of observations within gridbox climatological standard deviations'),
			  ('cell_methods','time: sum (interval: 1 month) area: sum where observations within gridbox)'),
			  ('units','1'),
			  ('missing_value',-1), # DOUBLE CHECK WHETHER THIS IS PERMISSABLE FOR INT?
			  ('_FillValue',-1)])]  

# Set up Global Attribute List of Lists
GlobAttrObjectList=dict([['File_created',datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')], # Is there a call for time stamping?
			     ['Description','HadISDH monthly mean '+Domain+' surface '+NameList[MyChoice][3]+' '+Description],
			     ['Title','HadISDH monthly mean '+Domain+' surface '+NameList[MyChoice][3]+' climate monitoring product'], 
			     ['Institution',Institution],
			     ['History','See http://catalogue.ceda.ac.uk/uuid/251474c7b09449d8b9e7aeaf1461858f for more information and related data.'+History], 
			     ['Licence',Licence],
			     ['Project', Project],
			     ['Processing_level',Processing_level],
			     ['Acknowledgement',Acknowledgement],
			     ['Source',Source],
			     ['Comment',''],
			     ['References',References],
			     ['Creator_name',Creator_name],
			     ['Creator_email',Creator_email],
			     ['Version',version],
			     ['doi',''], # This needs to be filled in
			     ['Conventions',Conventions],
			     ['netCDF_type',netCDF_type]]) 


#-----------------------------------------------------------------------
# Set up for q
if (MyChoice == 'q'):
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.landq.'+verstring+'_FLATgridIDPHA5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marineq.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blendq.'+verstring+'_FLATgridIDPHABClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='huss-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    
    								  						  
#------------------------------------------------------------------------	
# Set up for rh
if (MyChoice == 'rh'):	
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.landRH.'+verstring+'_FLATgridIDPHA5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marineRH.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blendRH.'+verstring+'_FLATgridIDPHABClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='hurs-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    
		
#------------------------------------------------------------------------	
# Set up for e
if (MyChoice == 'e'):	
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.lande.'+verstring+'_FLATgridIDPHA5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marinee.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blende.'+verstring+'_FLATgridIDPHABClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='vps-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    
        	
#------------------------------------------------------------------------	
# Set up for td
if (MyChoice == 'td'):	
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.landTd.'+verstring+'_FLATgridPHADPD5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marineTd.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blendTd.'+verstring+'_FLATgridPHADPDBClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='tds-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    
	
#------------------------------------------------------------------------	
# Set up for tw
if (MyChoice == 'tw'):	
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.landTw.'+verstring+'_FLATgridIDPHA5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marineTw.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blendTw.'+verstring+'_FLATgridIDPHABClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='tws-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    

#------------------------------------------------------------------------	
# Set up for t
if (MyChoice == 't'):	
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.landT.'+verstring+'_FLATgridIDPHA5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marineT.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blendT.'+verstring+'_FLATgridIDPHABClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='tas-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    
						
#------------------------------------------------------------------------	
# Set up for dpd
if (MyChoice == 'dpd'):	
    # Files
    if (Domain == 'land'):
        InFile='HadISDH.landDPD.'+verstring+'_FLATgridPHA5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'marine'):
        InFile='HadISDH.marineDPD.'+verstring+'_BClocalSHIP5by5both_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    elif (Domain == 'blend'):
        InFile='HadISDH.blendDPD.'+verstring+'_FLATgridPHABClocalSHIPboth5by5_anoms'+climbo+'_'+nowmon+nowyear+'_cf.nc'
    OutFile='dpds-'+Domain+'_HadISDH_HadOBS_'+version+'_'+str(StYr)+'0101-'+str(EdYr)+'1231.nc'    
						
#************************************************************************
# Main Program
#************************************************************************
# Get filenames
InFileName = InPath+InFile
OutFileName = OutPath+OutFile

# Set up dates dictionary
Dates = dict([('StYr',StYr),('StMon',StMon),('EdYr',EdYr),('EdMon',EdMon)])

# Read in data objects to make a list of numpy arrays DataObject and also supply Latitudes and Longitudes
#DataObject,Latitudes,Longitudes = GetGrid(InFileName,DataObjectList,LatInfo,LonInfo)
DataObject,Latitudes,Longitudes = GetGrid4(InFileName,DataObjectList,LatInfo,LonInfo)

# Check data mdis - land counts are currently -1e30 and need to be -1
if (Domain == 'land'):
    DataObject[10][np.where(DataObject[10] <= 0)] = -1
    DataObject[11][np.where(DataObject[11] <= 0)] = -1
#pdb.set_trace()
#ADD A CATCH HERE TO CHANGE MDI

# Call writing module
WriteNCCF(OutFileName,Dates,Latitudes,Longitudes,ClimPoints,DataObject,DimObjectList,AttrObjectList,GlobAttrObjectList)

# Finish
print("And we are done!")
