#!/usr/local/sci/bin/python
# PYTHON3
#
# test command line arguements
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
#
# # set up all editables as you wish within the script 
# > module load scitools/experimental-current
# > python PlotAverageMap_cartopy_JUL2018_timeplotaverage.py
#
# or with some/all command line arguments
# > module load scitools/experimental-current
# >python MakeObsCountList_APR2016.py --year1 '1973' --year2 '1973' --month1 '01' --month2 '01' --typee 'ERAclimNBC'
#
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# Set up python imports
import numpy as np
import sys, os, getopt
import struct
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)


def main(argv):
    ' Some info about what this code does and any things we might want to specify in the command line when running the program '
    #************************************************************************
    # INPUT PARAMETERS AS STRINGS!!!!
    year1 = '-999' # input integer years either as integers or strings
    year2 = '-999' # input integer years either as integers or strings
    month1 = '01' # input integer months 1-12 either as integers or strings with or without left 0 fill
    month2 = '12' # input integer months 1-12 either as integers or strings with or without left 0 fill
    switch = 'all'

    try:
        opts, args = getopt.getopt(argv, "hi:",
	                           ["year1=","year2=","month1=","month2=","switch="])
    except getopt.GetoptError:
        print('Usage (as strings) arg_test.py --year1 <1973> --year2 <1973> --month1 <01> --month2 <12> --switch <all>')
        sys.exit(2)
    pdb.set_trace()
    for opt, arg in opts:
        if opt == "--year1":
            try:
                year1 = int(arg)
            except:
                sys.exit("Failed: year1 not an integer")
        elif opt == "--year2":
            try:
                year2 = int(arg)
            except:
                sys.exit("Failed: year2 not an integer")
        elif opt == "--month1":
            try:
                month1 = int(arg)
            except:
                sys.exit("Failed: month1 not an integer")
            else:
                if (month1 < 1) | (month1 > 12):
                    sys.exit('month1 too high or too low, should be 1-12')	       
        elif opt == "--month2":
            try:
                month2 = int(arg)
            except:
                sys.exit("Failed: month2 not an integer")
            else:
                if (month2 < 1) | (month2 > 12):
                    sys.exit('month2 too high or too low, should be 1-12')	       
        elif opt == "--switch":
            try:
                switch = str(arg)
            except:
                sys.exit("Failed: switch not a string")
            else:
                if (switch  not in ['all','ship','buoys']):
                    sys.exit('switch invalid - should be one of [all, ship, buoys]')	       

    assert year1 != '-999' and year2 != '-999', "Year not specified."

    print(year1, year2, month1, month2, switch)
#    pdb.set_trace()
    return

#************************************************************************
# MAIN PROGRAM
#************************************************************************

if __name__ == '__main__':
    main(sys.argv[1:])
