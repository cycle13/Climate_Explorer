#!/usr/bin/python

# Code to read in years of hourly high resolution ERA5 t, td and land_sea_mask data
# Converts to q, RH, e, tw and DPD
# Aggregates to daily average
# Regrids to 1by1 degree and shifts to -179.5 ot 179.5 and 180 lats from 89.5 to -89.5 (was 181 lats!!!)
# Outputs as netCDF

# Later code will read in and convert to pentads, monthlies, anomalies etc and combine to make complete record or
# append the latest year

# This can also spin through all years or cope with annual updates
#*******************************************
# START
#*******************************************
import os
import datetime as dt
import calendar
import numpy as np
import sys
import time
import pdb
import iris
import iris.coord_categorisation
from iris.coords import DimCoord
from iris.cube import Cube
import cf_units

import CalcHums as ch

#import utils

#sys.path.append('/data/users/rdunn/reanalyses/code/era5/cdsapi-0.1.4')
#sys.path.append('/data/users/hadkw/WORKING_HADISDH/UPDATE2019/PROGS/PYTHON/cdsapi-0.1.4')
sys.path.append('/home/h04/hadkw/HadISDH_Code/HADISDH_BUILD/cdsapi-0.1.4')
import cdsapi

# Set up directory
#DataLoc = '/data/users/hadkw/WORKING_HADISDH/UPDATE2019/OTHERDATA/ERA5/'
DataLoc = '/scratch/hadkw/UPDATE2020/OTHERDATA/ERA5/'
print(DataLoc)

"""
Butchered from 

http://fcm1.metoffice.com/projects/utils/browser/CM_ML/trunk/NAO_Precip_Regr/get_era5_uwind.py
and /data/users/rdunn/reanalysis/era5/get_era.py

"""

#****************************************
def retrieve(year, variable, month, ndays):
    '''
    Use ECMWF API to get the data

    4.5GB per month --> 55GB per year, 50mins per month of processing
    '''

    if variable == "2m_temperature":
        varlist = ["2m_temperature"]
#        varlist = ["2m_temperature", "land_sea_mask"] # if you want to download both at once
    elif variable == "2m_dewpoint_temperature":
        varlist = ["2m_dewpoint_temperature"]
#        varlist = ["2m_dewpoint_temperature", "land_sea_mask"]
    elif variable == "surface_pressure":
        varlist = ["surface_pressure"]
    elif variable == "land_sea_mask":
        varlist = ["land_sea_mask"]
    else:
        print("please provide correct variable to download")
        return
    
    days = ["{:2d}".format(d+1) for d in range(ndays)]

    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':varlist,
            'year':"{}".format(year),
            'month':"{:02d}".format(month),
            'day':days,
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00',
            ]
        },
        os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable))
        )

#    time.sleep(5) # to allow any writing process to finish up.
#    # make a "success" file
#    with open(os.path.join(DataLoc, "{}{:02d}_{}_success.txt".format(year, month, variable)), "w") as outfile:
#        
#        outfile.write("Success {}".format(dt.datetime.now()))

    return # retreive

#****************************************
def check_files(year, variable, month, ndays):
    ''' This reads in the t, td and p files and checks for full download '''
    ''' If the last hour field has identical values for every lat and lon box then it has failed '''
    ''' A failed file is removed and program aborts '''
    ''' Program will need to be restarted '''
    
    action = 'contrinue' # output if file is ok
    
    test_cube = iris.load(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))
    # convert from list to cube
    test_cube = test_cube[0]
    test_data = test_cube.data[-1,:,:]
    
    # Are all values the same - np.unique() has a length of 1 if so
    if (len(np.unique(test_cube.data[-1,:,:])) == 1):
        
        # remove failed files
        os.remove(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))

        action = 'retrieve' # either retry download or exit code depending on stage in process
        ## exit the program
        #sys.exit('Incomplete file download')
	
    return action

#****************************************
def convert(year, month, ndays, remove=False):
    """
    Now need to:
    - convert to q, RH , e, tw, DPD
    - aggregate to daily averages
    - regrid to 1by1 gridboxes
    """

    MDI = -999.

# Set up null_cube with desired gridding format to use as a template 
# Does this have to have the same time dimensions?   
#    ndays = np.int(p_cube.data[:,0,0] / 24)
    
    time = DimCoord(np.arange(ndays*24),
                        standard_name = 'time',
                        units = 'hours')
    latitude = DimCoord(np.linspace(89.5, -89.5, 180),
#    latitude = DimCoord(np.linspace(90, -90, 181),
                        standard_name = 'latitude',
			long_name = 'gridbox centre latitude',
                        units = 'degrees_north')
    longitude = DimCoord(np.linspace(-179.5, 179.5, 360),
#    longitude = DimCoord(np.linspace(0, 359, 360),
                         standard_name='longitude',
			 long_name = 'gridbox centre longitude',
                        units = 'degrees_east')
    null_cube = Cube(np.zeros((ndays*24, 180, 360), np.float32),
                 dim_coords_and_dims=[(time, 0),
		                     (latitude, 1),
                                     (longitude, 2)])
    print('Check null_cube for new grid')
#    pdb.set_trace()

## START OF LSM************************************************
#    # read in land_sea_mask
#    variable = "land_sea_mask"
#    lsm_cube = iris.load(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))
#    #pdb.set_trace()
#    # convert from list to cube
#    lsm_cube = lsm_cube[0]#
#
## REgrid to 1by1 degree - cash the source, template, gridding type for later use - faster
#    regridder = iris.analysis.Linear().regridder(lsm_cube, null_cube)
#    lsm_cube_1by1 = regridder(lsm_cube)
#    print('Check lsm_cube_1by1 for new grid')
##    pdb.set_trace()#
#
#    # remove old cube
#    lsm_cube = 0
#
#    lsm_cube_1by1 = lsm_cube_1by1[0,:,:] 
##    lsm_cube_1by1_field = lsm_cube_1by1.extract(iris.Constraint(time=0)) 
#    lsm_cube_1by1.units = "1"
#    print(lsm_cube_1by1)
#    print('Check lsm_cube_1by1 for 2m_temperature')
#    #pdb.set_trace()
#
## output
#    iris.save(lsm_cube_1by1, os.path.join(DataLoc, "{}{:02d}_{}.nc".format(year, month, variable)), zlib=True)
#    print('Check lsm_cube_1by1 output')
#    pdb.set_trace()
## END OF LSM************************************************************

    # read in t, td and sp (may be VERY LARGE
    variable = "2m_temperature"
    t_cube = iris.load(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))
    #pdb.set_trace()
    # convert from list to cube
    t_cube = t_cube[0]

# REgrid to 1by1 degree - cash the source, template, gridding type for later use - faster
    regridder = iris.analysis.Linear().regridder(t_cube, null_cube)
    t_cube_1by1 = regridder(t_cube)
    print('Check t_cube_1by1 for new grid')
#    pdb.set_trace()

    # remove old cube
    t_cube = 0

    t_cube_1by1.data -= 273.15 # convert to C
    t_cube_1by1.units = "degreesC"
    print('Check t_cube_1by1 for 2m_temperature')
    #pdb.set_trace()

    variable = "2m_dewpoint_temperature"
    td_cube = iris.load(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))
    # convert from list to cube
    td_cube = td_cube[0]

# REgrid to 1by1 degree - cash the source, template, gridding type for later use - faster
    td_cube_1by1 = regridder(td_cube)
    print('Check td_cube_1by1 for new grid')
#    pdb.set_trace()

    # remove old cube
    td_cube = 0

    td_cube_1by1.data -= 273.15 # convert to C
    td_cube_1by1.units = "degreesC"
    print('Check td_cube_1by1 for 2m_dewpoint_temperature')
#    pdb.set_trace()

    variable = "surface_pressure"
    p_cube = iris.load(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))
    # convert from list to cube
    p_cube = p_cube[0]

# REgrid to 1by1 degree - cash the source, template, gridding type for later use - faster
    p_cube_1by1 = regridder(p_cube)
    print('Check p_cube_1by1 for new grid')
#    pdb.set_trace()

    # remove old cube
    p_cube = 0

    p_cube_1by1.data /= 100. # convert to C
    p_cube_1by1.units = "hPa"
    print('Check p_cube_1by1 for surface_pressure')
#    pdb.set_trace()

#    # if it contains 2 cubes where we have downloaded mask and wish to mask to land or sea....
#    if len(p_cubelist) == 2:
#        # extract both cubes
#        pcube1 = p_cubelist[0]
#        pcube2 = p_cubelist[1]#
#
#        masked1, = np.where(pcube1.data.mask[:, 0, 0] == True)
#        masked2, = np.where(pcube2.data.mask[:, 0, 0] == True)
#
#        # use locations of masks to overwrite
#        tp_cube = pcube1[:]
#        tp_cube.data[masked1] = pcube2.data[masked1]
#        tp_cube.var_name = "tp"
#
#    # else it's just a single cube, so easier to deal with
#    elif len(p_cubelist) == 1:#
#
#        tp_cube = p_cubelist[0]
#        tp_cube.var_name = "tp"


# No masking internally within this code...
# Process q
# Copy the t_cube and then change some of the fields?
    variable = 'specific_humidity'
    q_cube = t_cube_1by1.copy()
    q_cube.fill_value = MDI # not sure whether we're doing -999 yet if saving as integer
    q_cube.units = cf_units.Unit("g kg-2")
    q_cube.var_name = "q2m"
    q_cube.long_name = "2 metre specific humidity"

# Populate the q data
    q_cube.data = ch.sh(td_cube_1by1.data,t_cube_1by1.data,p_cube_1by1.data,roundit=False)
    print('Check q_cube for new data')
#    pdb.set_trace()
    
    ## mask all regions which are 100% ocean
    #cube.data[lsm.data == 0] = utils.MDI
    #cube.data = np.ma.masked_where(lsm.data == 0, cube.data)
    #cube.data.fill_value = utils.MDI

# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(q_cube, "time", name="day_of_month")
    
    q_cube_day = q_cube.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    q_cube = 0
    q_cube_day.remove_coord("day_of_month")
    q_cube_day.units = cf_units.Unit("g kg-2")
    print('Check q_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(q_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    q_cube_day=0
    print('Check q_cube_1by1 output')
#    pdb.set_trace()
    
# Process RH
# Copy the t_cube and then change some of the fields?
    variable = 'relative_humidity'
    rh_cube = t_cube_1by1.copy()
    rh_cube.fill_value = MDI # not sure whether we're doing -999 yet if saving as integer
    rh_cube.units = cf_units.Unit("%")
    rh_cube.var_name = "rh2m"
    rh_cube.long_name = "2 metre relative humidity"

# Populate the q data
    rh_cube.data = ch.rh(td_cube_1by1.data,t_cube_1by1.data,p_cube_1by1.data,roundit=False)
    print('Check rh_cube for new data')
#    pdb.set_trace()
    
    ## mask all regions which are 100% ocean
    #cube.data[lsm.data == 0] = utils.MDI
    #cube.data = np.ma.masked_where(lsm.data == 0, cube.data)
    #cube.data.fill_value = utils.MDI

# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(rh_cube, "time", name="day_of_month")
    
    rh_cube_day = rh_cube.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    rh_cube = 0
    rh_cube_day.remove_coord("day_of_month")
    rh_cube_day.units = cf_units.Unit("%")
    print('Check rh_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(rh_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    rh_cube_day=0
    print('Check rh_cube_1by1 output')
    #pdb.set_trace()

# Process e
# Copy the t_cube and then change some of the fields?
    variable = 'vapour_pressure'
    e_cube = t_cube_1by1.copy()
    e_cube.fill_value = MDI # not sure whether we're doing -999 yet if saving as integer
    e_cube.units = cf_units.Unit("hPa")
    e_cube.var_name = "e2m"
    e_cube.long_name = "2 metre vapour pressure"

# Populate the q data
    e_cube.data = ch.vap(td_cube_1by1.data,t_cube_1by1.data,p_cube_1by1.data,roundit=False)
    print('Check e_cube for new data')
#    pdb.set_trace()
    
    ## mask all regions which are 100% ocean
    #cube.data[lsm.data == 0] = utils.MDI
    #cube.data = np.ma.masked_where(lsm.data == 0, cube.data)
    #cube.data.fill_value = utils.MDI

# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(e_cube, "time", name="day_of_month")
    
    e_cube_day = e_cube.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    e_cube = 0
    e_cube_day.remove_coord("day_of_month")
    e_cube_day.units = cf_units.Unit("hPa")
    print('Check e_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(e_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    e_cube_day=0
    print('Check e_cube_1by1 output')
#    pdb.set_trace()

# Process tw
# Copy the t_cube and then change some of the fields?
    variable = 'wetbulb_temperature'
    tw_cube = t_cube_1by1.copy()
    tw_cube.fill_value = MDI # not sure whether we're doing -999 yet if saving as integer
    tw_cube.units = cf_units.Unit("degrees C")
    tw_cube.var_name = "tw2m"
    tw_cube.long_name = "2 metre wetbulb temperature"

# Populate the q data
    tw_cube.data = ch.wb(td_cube_1by1.data,t_cube_1by1.data,p_cube_1by1.data,roundit=False)
    print('Check tw_cube for new data')
#    pdb.set_trace()
    
    ## mask all regions which are 100% ocean
    #cube.data[lsm.data == 0] = utils.MDI
    #cube.data = np.ma.masked_where(lsm.data == 0, cube.data)
    #cube.data.fill_value = utils.MDI

# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(tw_cube, "time", name="day_of_month")
    
    tw_cube_day = tw_cube.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    tw_cube = 0
    tw_cube_day.remove_coord("day_of_month")
    tw_cube_day.units = cf_units.Unit("degrees C")
    print('Check tw_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(tw_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    tw_cube_day=0
    print('Check tw_cube_1by1 output')
#    pdb.set_trace()

# Process dpd
# Copy the t_cube and then change some of the fields?
    variable = 'dewpoint_depression'
    dpd_cube = t_cube_1by1.copy()
    dpd_cube.fill_value = MDI # not sure whether we're doing -999 yet if saving as integer
    dpd_cube.units = cf_units.Unit("degrees C")
    dpd_cube.var_name = "dpd2m"
    dpd_cube.long_name = "2 metre dewpoint depression"

# Populate the q data
    dpd_cube.data = ch.dpd(td_cube_1by1.data,t_cube_1by1.data,roundit=False)
    print('Check dpd_cube for new data')
#    pdb.set_trace()
    
    ## mask all regions which are 100% ocean
    #cube.data[lsm.data == 0] = utils.MDI
    #cube.data = np.ma.masked_where(lsm.data == 0, cube.data)
    #cube.data.fill_value = utils.MDI

# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(dpd_cube, "time", name="day_of_month")
    
    dpd_cube_day = dpd_cube.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    dpd_cube = 0
    dpd_cube_day.remove_coord("day_of_month")
    dpd_cube_day.units = cf_units.Unit("degrees C")
    print('Check dpd_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(dpd_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    dpd_cube_day=0
    print('Check dpd_cube_1by1 output')
#    pdb.set_trace()

# Process Td
    variable = '2m_dewpoint_temperature'
# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(td_cube_1by1, "time", name="day_of_month")
    
    td_cube_day = td_cube_1by1.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    td_cube_1by1 = 0
    td_cube_day.remove_coord("day_of_month")
    td_cube_day.units = cf_units.Unit("degrees C")
    td_cube_day.var_name = "td2m"
    print('Check td_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(td_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    td_cube_day=0
    print('Check td_cube_1by1 output')
#    pdb.set_trace()

# Process T
    variable = '2m_temperature'
# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(t_cube_1by1, "time", name="day_of_month")
    
    t_cube_day = t_cube_1by1.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    t_cube_1by1 = 0
    t_cube_day.remove_coord("day_of_month")
    t_cube_day.units = cf_units.Unit("degrees C")
    t_cube_day.var_name = "t2m"
    print('Check t_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(t_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    t_cube_day=0
    print('Check t_cube_1by1 output')
#    pdb.set_trace()

# Process P
    variable = 'surface_pressure'
# Aggregate to daily
    # add a "day" indicator to allow aggregation
    iris.coord_categorisation.add_day_of_month(p_cube_1by1, "time", name="day_of_month")
    
    p_cube_day = p_cube_1by1.aggregated_by(["day_of_month"], iris.analysis.MEAN)
    p_cube_1by1 = 0
    p_cube_day.remove_coord("day_of_month")
    p_cube_day.units = cf_units.Unit("hPa")
    p_cube_day.var_name = "p2m"
    print('Check p_cube for daily averages')
#    pdb.set_trace()

# output
    iris.save(p_cube_day, os.path.join(DataLoc, "{}{:02d}_daily_{}.nc".format(year, month, variable)), zlib=True)
    p_cube_day=0
    print('Check p_cube_1by1 output')
#    pdb.set_trace()

#    # append precipitation cube to temperature one
#    cubelist += [tp_cube]

    # remove input files
    if remove:
        for variable in ["2m_temperature", "2m_dewpoint_temperature", "surface_pressure"]:
#        for variable in ["2m_temperature", "2m_dewpoint_temperature", "surface_pressure", "land_sea_mask"]:
            os.remove(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))

    return # combine
    
#****************************************
if __name__ == "__main__":

    import argparse

#    # Set up directory
#    DataLoc = '/data/users/hadkw/WORKING_HADISDH/UPDATE2019/OTHERDATA/ERA5/'
#    print(DataLoc)

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', dest='start', action='store', default=1979, type=int,
                        help='Start year [1979]')
    parser.add_argument('--end', dest='end', action='store', default=2019, type=int,
                        help='End year [2019]')
    parser.add_argument('--remove', dest='remove', action='store_true', default=False,
                        help='Remove hourly and monthly files, default = False')
 
    args = parser.parse_args()         

    print(args)

    for year in np.arange(args.start, args.end+1):
        print(year)

        if not os.path.exists(os.path.join(DataLoc, "{}12_daily_surface_pressure.nc".format(year))):
                    
            for month in np.arange(1, 13):

                print("{} - {}".format(year, month))
                # get number of days
                ndays = calendar.monthrange(year, month)[1]
                
                if not os.path.exists(os.path.join(DataLoc, "{}{:02d}_daily_surface_pressure.nc".format(year, month))):
                    
                    for variable in ["2m_temperature", "2m_dewpoint_temperature", "surface_pressure"]:
                                    
			# If there is a file present then check it first
                        action = 'retrieve' # instruction as to whether to retrieve (exit code if second try doesn't work), not bother (carry on)

                        if os.path.exists(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable))):

                            print("{} - {} - {} already downloaded so testing".format(year, month, variable))    
                            action = check_files(year, variable, month, ndays)
                            print('Code to', action)    
                        
                        if (action == 'retrieve'):
			
                            retrieve(year, variable, month, ndays) 
                            action = 'continue'
                            action = check_files(year, variable, month, ndays)
			    
                        if (action == 'retrieve'): # then the download has failed for some reason so stop the code and restart    
                        
                            sys.exit('Failed download from CDS')   
			
                    convert(year, month, ndays, remove = args.remove)
            
                else:
                                    
                    print("{} - {} already downloaded".format(year, month))    

        else:
                                    
            print("{} already downloaded".format(year))    

#*******************************************
# END
#*******************************************
