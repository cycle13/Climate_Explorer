#!/usr/bin/env python
import calendar
# import ecmwfapi
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

def retrieve_era_interim():
    yearStart = 1979                        # adjust to your requirements - as of 2017-07 only 2010-01-01 onwards is available
    yearEnd = 2017                          # adjust to your requirements
    months = [1,2,3,4,5,6,7,8,9,10,11,12]   # adjust to your requirements
    
    years = range(yearStart, yearEnd+1)
    print 'Years: ',years
    decades = list(set([divmod(i, 10)[0] for i in years]))
    decades = [x * 10 for x in decades]
    decades.sort()
    print 'Decades:', decades
    
    # loop through decades and create a month list
    for d in decades:
        requestDates=''
        for y in years:
            if ((divmod(y,10)[0])*10) == d:
                for m in months:
                    requestDates = requestDates+str(y)+(str(m)).zfill(2)+'01/'
        requestDates = requestDates[:-1]
        print 'Requesting dates: ', requestDates
        target = 'era_interim_moda_%d'% (d)    # specifies the output file name
        print 'Output file: ', target
        era_interim_request(requestDates, d, target)

# the actual data request
def era_interim_request(requestDates, decade, target):
    '''
        Change the keywords below to adapt to your needs.
        The easiest way to do this is:
        1. go to http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/
        2. click through to the parameters you want, for any date
        3. click 'View the MARS request'
        '''
    server.retrieve({
                    "class": "ei",              # do not change
                    "dataset": "interim",       # do not change
                    'expver': '1',              # do not change
                    'stream': 'moda',           # monthly means of daily means
                    'type': 'an',               # analysis (versus forecast, fc)
                    'levtype': 'sfc',           # surface data (versus pressure level, pl, and model level, ml)
                    'param': '165.128/166.128',         # here: sea surface temperature (param 34) and mean sea level pressure (param 151)
                    'grid': '5/5',        # horizontal resolution of output in degrees lat/lon
                    'format': 'netcdf',         # get output in netcdf; only works with regular grids; for GRIB remove this line
                    'date': requestDates,       # dates, set automatically from above
                    'decade': decade,           # decade set automatically from above
                    'target': target            # output file name, set automatically from above
                    })

if __name__ == '__main__':
    retrieve_era_interim()
