#!/usr/bin/python
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


# Set up null_cube with desired gridding format to use as a template 
# DAILY MODEL DATA - CHANGE DEPENDING ON WHETHER YOU ARE WORKING WITH HOURLY, DAILY, MONTHLY
ndays = 360 # if 1 year of data
nmons = 12 # if 1 year of data
nyrs = 1 # if 1 year of data
time = DimCoord(np.arange(ndays),
                        standard_name = 'time',
                        units = 'days')
# Set up to regrid to HadISDH like			
latitude = DimCoord(np.linspace(-87.5, 87.5, 36),
#    latitude = DimCoord(np.linspace(90, -90, 181),
                        standard_name = 'latitude',
			long_name = 'gridbox centre latitude',
                        units = 'degrees_north')
longitude = DimCoord(np.linspace(-177.5, 177.5, 72),
#    longitude = DimCoord(np.linspace(0, 359, 360),
                         standard_name='longitude',
			 long_name = 'gridbox centre longitude',
                        units = 'degrees_east')
null_cube = Cube(np.zeros((ndays, 36, 72), np.float32),
                 dim_coords_and_dims=[(time, 0),
		                     (latitude, 1),
                                     (longitude, 2)])
#    print('Check null_cube for new grid')
#    pdb.set_trace()

# read in t, td and sp (may be VERY LARGE
variable = "hursa" # hursa relative humidity
# if the model file contains only one variable then this should work
old_cube = iris.load(filename)
# convert from list to cube
old_cube = old_cube[0]
# if the model file contains more than one then you will need to tell if which one to read in - the following uses the long_name but I think you can use other things
#old_cube = iris.load(filename, file_long_name)[0]    # I think this already converts from a list to a cube
#pdb.set_trace()

# REgrid to 5by5 degree - cash the source, template, gridding type for later use - faster
regridder = iris.analysis.Linear().regridder(old_cube, null_cube)
new_cube = regridder(old_cube)
print('Check new_cube for new grid then type c and <return> to carry on')
pdb.set_trace()

# output
iris.save(new_cube, outfilename, zlib=True)
