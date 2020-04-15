#!/usr/bin/python

# Code to read in ERA5 land_sea_mask data and save as 1by1 single field
# Outputs as netCDF


# This can also spin through all years or cope with annual updates
#*******************************************
# START
#*******************************************
import os
import datetime as dt
import calendar
import iris
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
sys.path.append('/data/users/hadkw/WORKING_HADISDH/UPDATE2019/PROGS/PYTHON/cdsapi-0.1.4')
import cdsapi

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

    variable == "land_sea_mask"
    varlist = ["land_sea_mask"]

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

    time.sleep(5) # to allow any writing process to finish up.
    # make a "success" file
    with open(os.path.join(DataLoc, "{}{:02d}_{}_success.txt".format(year, month, variable)), "w") as outfile:
        
        outfile.write("Success {}".format(dt.datetime.now()))

    return # retreive

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
    latitude = DimCoord(np.linspace(89.5, -89.5, 180), # DIFF FROM OTHER ERA5 DOWNLOAD
#    latitude = DimCoord(np.linspace(90, -90, 181),
                        standard_name = 'latitude',
			long_name = 'gridbox centre latitude',
                        units = 'degrees_north')
    longitude = DimCoord(np.linspace(-179.5, 179.5, 360), # DIFF FROM OTHER ERA5 DOWNLOAD
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

# START OF LSM************************************************
    # read in land_sea_mask
    variable = "land_sea_mask"
    lsm_cube = iris.load(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))
    #pdb.set_trace()
    # convert from list to cube
    lsm_cube = lsm_cube[0]#

# REgrid to 1by1 degree - cash the source, template, gridding type for later use - faster
    regridder = iris.analysis.Linear().regridder(lsm_cube, null_cube)
    lsm_cube_1by1 = regridder(lsm_cube)
    print('Check lsm_cube_1by1 for new grid')
#    pdb.set_trace()#

    # remove old cube
    lsm_cube = 0

    lsm_cube_1by1 = lsm_cube_1by1[0,:,:] 
#    lsm_cube_1by1_field = lsm_cube_1by1.extract(iris.Constraint(time=0)) 
    lsm_cube_1by1.units = "1"
    print(lsm_cube_1by1)
    print('Check lsm_cube_1by1 for land_sea_mask')
    #pdb.set_trace()

# output
    iris.save(lsm_cube_1by1, os.path.join(DataLoc, "{}{:02d}_hourly_{}_ERA5.nc".format(year, month, variable)), zlib=True)
    print('Check lsm_cube_1by1 output')
#    pdb.set_trace()
# END OF LSM************************************************************

    # remove input files
    if remove:
        for variable in ["2m_temperature", "2m_dewpoint_temperature", "surface_pressure"]:
#        for variable in ["2m_temperature", "2m_dewpoint_temperature", "surface_pressure", "land_sea_mask"]:
            os.remove(os.path.join(DataLoc, "{}{:02d}_hourly_{}.nc".format(year, month, variable)))

    return # combine
    
#****************************************
if __name__ == "__main__":

    import argparse

    # Set up directory
    DataLoc = '/data/users/hadkw/WORKING_HADISDH/UPDATE2019/OTHERDATA/'
    print(DataLoc)

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

        for month in np.arange(1, 2): # I think all months are the same

            print("{} - {}".format(year, month))

            # get number of days
            ndays = calendar.monthrange(year, month)[1]
               
            for variable in ["land_sea_mask"]:
                                
                retrieve(year, variable, month, ndays)
                convert(year, month, ndays, remove = args.remove)

#*******************************************
# END
#*******************************************
