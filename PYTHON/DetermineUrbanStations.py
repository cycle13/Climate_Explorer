#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 27 April 2021 
# (adapted from David Parker's station_9pixels_landuse_ISTI.py)
# Last update: 27 April 2021
# Location: /home/h04/hadkw/HadISDH_Code/CLIMEXPLORER/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Code reads in a list of station lats and longs and matches with pixels from a
# land use satellite dataset to determine likely land use in a variety of ways:
# - a) Urban = 1, Not urban = 0
# - b) average of pixel + 8 surrounding pixels urban = 1.0, not urban = 0.0
# - c) as b) but weight 0.6 [centre], 0.15/(1+cos(lat)) [E&W], 
#   0.15*cos(lat)/(1+cos(lat)) [N&S], 0.025 [corners]
# - d) Reduced Urban Signal by surrounding water, desert or sparse vegation = 
#   as c) * function of neighbouring 8 pixels being water, barren or sparse veg
#   if all surrounding pixels are water, desert or sparse veg then urban signal 
#   is reduced to 0.
# -----------------------
# DATA
# -----------------------
# A station file, code needs to be set up to read it
# - HadISD 
#   /data/users/hadkw/WORKING_HADISDH/UPDATE<YYYY>/LISTS_DOCS/
#    HadISD.<versiondots>_candidate_stations_details.txt
#
# A land use pixel dataset, code needs to be set up to read it
# - ESACCI_LandCover
#  /project/LandCoverCCI/V1/ESACCI-LC-Map-nc/
#   ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-CurrentVersion-inlandwater.nc
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# > module load scitools/default-current
# > python DetermineUrbanStations.py
# -----------------------
# OUTPUT
# -----------------------
# A statopm list:
# - HadISD:
#   /data/users/hadkw/WORKING_HADISDH/UPDATE<YYYY>/LISTS_DOCS/
#    HadISD.<versiondots>_candidate_stations_details_landuse.txt
#   Station_ID Station_Name Latitude Longitude Elevation Land_use Urban_a 
#   Urban_b Urban_c Urban_d Pixel_NW Pixel_N Pixel_NE Pixel_W Pixel_E Pixel_SW 
#   Pixel_S Pixel_SE
#   f"{x:>12} {x:<31} {x:7.3f} {x:8.3f} {x:7.1f} {x:4} {x:1} {x:8.2f} {x:8.2f} 
#   {x:8.2f} {x:4} {x:4} {x:4} {x:4} {x:4} {x:4} {x:4} {x:4}"
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 (27th April 2021)
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
# LAND USE FLAG MEANINGS. VALUES ARE IN BYTE FORMAT, RANGE 0-220
# Land use code /project/LandCoverCCI/user-tool/current/resources/LCCS_class_definitions.csv: \n"
#  0 no_data\n"
#  10 cropland_rainfed\n"
#  11 cropland_rainfed_herbaceous_cover\n"
#  12 cropland_rainfed_tree_or_shrub_cover\n"
#  20 cropland_irrigated\n"
#  30 mosaic_cropland\n"
#  40 mosaic_natural_vegetation\n"
#  50 tree_broadleaved_evergreen_closed_to_open\n"
#  60 tree_broadleaved_deciduous_closed_to_open\n"
#  61 tree_broadleaved_deciduous_closed\n"
#  62 tree_broadleaved_deciduous_open\n"
#  70 tree_needleleaved_evergreen_closed_to_open\n"
#  71 tree_needleleaved_evergreen_closed\n"
#  72 tree_needleleaved_evergreen_open\n"
#  80 tree_needleleaved_deciduous_closed_to_open\n"
#  81 tree_needleleaved_deciduous_closed\n"
#  82 tree_needleleaved_deciduous_open\n"
#  90 tree_mixed\n"
#  100 mosaic_tree_and_shrub\n"
#  110 mosaic_herbaceous\n"
#  120 shrubland\n"
#  121 shrubland_evergreen\n"
#  122 shrubland_deciduous\n"
#  130 grassland\n"
#  140 lichens_and_mosses\n"
#  150 sparse_vegetation\n" #*
#  152 sparse_shrub\n" #*
#  153 sparse_herbacious\n" #*
#  160 tree_cover_flooded_fresh_or_brakish_water\n"
#  170 tree_cover_flooded_saline_water\n"
#  180 shrub_or_herbaceous_cover_flooded\n"
#  190 urban\n"
#  200 bare_areas\n" #*
#  201 bare_areas_consolidated\n" #*
#  202 bare_areas_unconsolidated\n" #*
#  210 water\n" #*
#  211 water(other? - not listed in definitions file but present in netCDF)\n" #*
#  220 snow_and_ice\n"
#************************************************************************
# Set up python imports
import iris
import numpy as np
import pdb
import netCDF4 as nc4  

# Set up global variables
InStation = "/data/users/hadkw/WORKING_HADISDH/UPDATE2020/LISTS_DOCS/HadISD.3.1.2.202101p_candidate_stations_details.txt"  
InLandUse = "/project/LandCoverCCI/V1/ESACCI-LC-Map-nc/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-CurrentVersion-inlandwater.nc"

OutStation = "/data/users/hadkw/WORKING_HADISDH/UPDATE2020/LISTS_DOCS/HadISD.3.1.2.202101p_candidate_stations_details_landuse.txt"  

#************************************************************************
# SUBROUTINES
#************************************************************************

#************************************************************************
# MAIN
#************************************************************************
cubeg = iris.load(InLandUse) # This loads the global file 
# cube[0] is not what we want 
cubeg1 = cubeg[1] 

#ncf = nc4.Dataset(InLandUse,'r')
#LCLatValues = ncf.variables['lat'][:]
#LCLongValues = ncf.variables['lon'][:]
#pdb.set_trace()
#var = ncf.variables['lccs_class'][:]  

LCLatValues = cubeg1.coord('latitude').points
nlats = len(cubeg1.coord('latitude').points) # nlats is 64800
LCLongValues = cubeg1.coord('longitude').points 
# Latitude decreases (moves south) but longitude increases (moves east) with 
# increasing array element. The increments are approximately 0.002777 degrees 
# of latitude or longitude but are not absolutely fixed.
nlongs = len(cubeg1.coord('longitude').points) # nlongs is 129600
# These are the maximum and minimum latitude increments in the land use dataset 
dy_max = 0.00278473 
dy_min = 0.0027771
# These are the maximum and minimum longitude increments in the land use dataset 
dx_max = 0.00279236 
dx_min = 0.0027771

station_out = open(OutStation, 'a')
station_out.write(
    f"ESACCI-LC-Map-nc/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-CurrentVersion-inlandwater.nc"
    f"Urban_a = 1.0 if urban, 0.0 if not"
    f"Urban_b = mean of 9 pixels around target (1.0 if urban, 0.0 if not)"
    f"Urban_c = cosine weighted mean of 9 pixels around target (1.0 if urban, 0.0 if not)"
    f"Urban_d = as for c but reduced by cosine weighted function of presence of desert, water or sparse vegetation in surrounding 8 pixels"
    f"Land use code LCCS_class_definitions.csv: \n"
    f"  0 no_data\n"
    f"  10 cropland_rainfed\n"
    f"  11 cropland_rainfed_herbaceous_cover\n"
    f"  12 cropland_rainfed_tree_or_shrub_cover\n"
    f"  20 cropland_irrigated\n"
    f"  30 mosaic_cropland\n"
    f"  40 mosaic_natural_vegetation\n"
    f"  50 tree_broadleaved_evergreen_closed_to_open\n"
    f"  60 tree_broadleaved_deciduous_closed_to_open\n"
    f"  61 tree_broadleaved_deciduous_closed\n"
    f"  62 tree_broadleaved_deciduous_open\n"
    f"  70 tree_needleleaved_evergreen_closed_to_open\n"
    f"  71 tree_needleleaved_evergreen_closed\n"
    f"  72 tree_needleleaved_evergreen_open\n"
    f"  80 tree_needleleaved_deciduous_closed_to_open\n"
    f"  81 tree_needleleaved_deciduous_closed\n"
    f"  82 tree_needleleaved_deciduous_open\n"
    f"  90 tree_mixed\n"
    f"  100 mosaic_tree_and_shrub\n"
    f"  110 mosaic_herbaceous\n"
    f"  120 shrubland\n"
    f"  121 shrubland_evergreen\n"
    f"  122 shrubland_deciduous\n"
    f"  130 grassland\n"
    f"  140 lichens_and_mosses\n"
    f"  150 sparse_vegetation\n" #*
    f"  152 sparse_shrub\n" #*
    f"  153 sparse_herbacious\n" #*
    f"  160 tree_cover_flooded_fresh_or_brakish_water\n"
    f"  170 tree_cover_flooded_saline_water\n"
    f"  180 shrub_or_herbaceous_cover_flooded\n"
    f"  190 urban\n"
    f"  200 bare_areas\n" #*
    f"  201 bare_areas_consolidated\n" #*
    f"  202 bare_areas_unconsolidated\n" #*
    f"  210 water\n" #*
    f"  211 water(other? - not listed in definitions file but present in netCDF)\n" #*
    f"  220 snow_and_ice\n"
)
station_out.write(
    f"Station_ID Station_Name Latitude Longitude Elevation Land_use Urban_a "
    f"Urban_b Urban_c Urban_d Pixel_NW Pixel_N Pixel_NE Pixel_W Pixel_E "
    f"Pixel_SW Pixel_S Pixel_SE\n"
)
#   f"{x:>12} {x:<31} {x:7.3f} {x:8.3f} {x:7.1f} {x:4} {x:1} {x:8.2f} {x:8.2f} 
#   {x:8.2f} {x:4} {x:4} {x:4} {x:4} {x:4} {x:4} {x:4} {x:4}"
station_out.close()

with open(InStation, 'r+b') as station_in:

    for line in station_in:

        StationID = line.decode('utf8')[0:12]   
        StationName = line.decode('utf8')[13:43]   
        Latitude = float(line.decode('utf8')[44:51])   
        Longitude = float(line.decode('utf8')[52:60])   
        Elevation = float(line.decode('utf8')[61:68])   

        RadianLatitude = np.radians(Latitude) 

        # Treat the South pole separately later in this script
        if Latitude > -89.5: 

            # Maximum possible latitude array-element corresponding to the 
            # station given the range of latitude increments            
            imx = int(1.0 + (90.0 - Latitude) / dy_min)
            # Minimum possible latitude array-element corresponding to the 
            # station given the range of latitude increments
            imn = int(-1.0 + (90.0 - Latitude) / dy_max) 
            # Maximum possible longitude array-element corresponding to the 
            # station given the range of longitude increments            
            jmx = int(1.0 + (Longitude + 180) / dx_min) 
            # Minimum possible longitude array-element corresponding to the 
            # station given the range of longitude increments            
            jmn = int(-1.0 + (Longitude + 180) / dx_max) 
      
            print(StationID, Latitude, Longitude, imn, imx, jmn, jmx) 

            # Now iterate through the land use tiles, finding the latitude and 
            # longitude ranges of each one until a tile matches the station 
            # location. 
            # To reduce this loop, run from imn = (90-Latitude) / dy_max to 
            # imx = (90-Latitude) / dy_min 
            # and from jmn = (Longitude+180) / dx_max to 
            # jmx = (Longitude+180) / dx_min, each one to the nearest integer 
            # mins(maxs) rounded down(up).

            # Get only the data for the range of array-elements imn through 
            # imx, jmn through jmx
            cubegd = cubeg1[imn:imx+1, jmn:jmx+1].data   
            # Set a 'not found' value for the station land use. 
            # This value does not occur in the input land-use dataset; if it 
            # is output, the station has not been found!            
            stlanduse = -256 
            # Array of land use in station pixel and in surrounding pixels
            pxlu9 = np.zeros((3,3), dtype = int)  
            # Array of urban indicators (1.0 urban [-66], 0.0 not urban) in 
            # station pixel and in surrounding pixels            
            pxur9 = np.zeros((3,3), dtype = float)  
            # Array of Water/Desert indicators (1.0 water or desert [150, 152, 
            # 153, 200, 201, 202, 210, 211?], 0.0 not water or desert) in 
            #station pixel and in surrounding pixels            
            pxwd9 = np.zeros((3,3), dtype = float)  
      
            for i in range(imn, imx+1):
                for j in range(jmn, jmx+1): 
                    latmax = 90.0 # North Pole
                    latmin = -89.0 # most southerly latitude considered here
                    longmax = 180.0 # Dateline at eastern edge
                    longmin = -180.0 # Dateline at western edge
          
                    if i != 0: 
                        latmax = LCLatValues[i] + 0.5 * \
                                 (LCLatValues[i-1]  - LCLatValues[i])
          
                    if i != nlats-1: 
                        latmin = LCLatValues[i] - 0.5 * \
                                (LCLatValues[i]  - LCLatValues[i+1])  

                    if j != nlongs-1: 
                        longmax = LCLongValues[j] + 0.5 * \
                                  (LCLongValues[j+1]  - LCLongValues[j])

                    if j != 0: 
                        longmin = LCLongValues[j] - 0.5 * \
                                  (LCLongValues[j]  - LCLongValues[j-1]) 

                    if (Latitude <= latmax) and \
                       (Latitude >= latmin) and \
                       (Longitude <= longmax) and \
                       (Longitude >= longmin): 
            
                        for ii in range(3):
              
                            for jj in range(3):
                
                                pxlu9[ii,jj] = cubegd[i-imn+ii-1, j-jmn+jj-1]
                       
                                # pxur9 gives the urban status of each of the 
                                #9 pixels centred on the station
                                if pxlu9[ii,jj] == 190: 
                                    pxur9[ii,jj] = 1.0  

                                # pxwd9 gives the water/desert status of each of the 
                                #9 pixels centred on the station
                                if (pxlu9[ii,jj] == 150) or \
                                   (pxlu9[ii,jj] == 151) or \
                                   (pxlu9[ii,jj] == 152) or \
                                   (pxlu9[ii,jj] == 200) or \
                                   (pxlu9[ii,jj] == 201) or \
                                   (pxlu9[ii,jj] == 202) or \
                                   (pxlu9[ii,jj] == 210) or \
                                   (pxlu9[ii,jj] == 211): # water or desert or sparse veg pixels

                                    pxwd9[ii,jj] = 1.0  
#                        pdb.set_trace()
                        break

        # Fix for South pole: -36 denotes snow and ice.
        elif Latitude <= -89.5:          
            
            pxlu9[:,:] = 220
            pxur9[:,:] = 0.0

        # This is the land use in the station pixel
        print("Station Land Use Pixel: ",pxlu9[1,1])    

        # Calculate alternative measures of urbanness using surrounding pixels
        coslat = np.cos(RadianLatitude)  
        # e.g. at 60N coslat=0.5: NorthSouthWeight=0.5, EastWestWeight=0.1 
        # are the weights of the north & south, east & west side boxes
        NorthSouthWeight = 0.15 * coslat / (1.0 + coslat)
        EastWestWeight = 0.15 / (1.0 + coslat) 
        # Average urbanization in the 9-pixel block centred on the station
        StationUrbanAverage = np.mean(pxur9)
        # Cosine weighted average urbanization (centre 0.6; sides 0.3; corners 
        # 0.1) in the 9-pixel block centred on the station. 
        # East West weighted more than North South sides by factor 1/coslat    
        StationWeightedUrbanAverage = (0.025 * pxur9[0,0] + 
                                       NorthSouthWeight * pxur9[0,1] + 
                                       0.025 * pxur9[0,2] + 
                                       EastWestWeight * pxur9[1,0] + 
                                       0.6 * pxur9[1,1] + 
                                       EastWestWeight * pxur9[1,2] + 
                                       0.025 * pxur9[2,0] + 
                                       NorthSouthWeight * pxur9[2,1] + 
                                       0.025 * pxur9[2,2]) 
        # Indicator for summing the prevalence of water or desert in 
        # surrounding pixels.
        WaterDesert = 0.0 
        # Count neighbouring pixels with either water or desert, with 
        # weighting the same as above, thus having a maximum total weight 
        # of 0.4. Do not include target pixel!
        WaterDesertWeight = (0.025 * pxwd9[0,0] + 
                             NorthSouthWeight * pxwd9[0,1] + 
                             0.025 * pxwd9[0,2] + 
                             EastWestWeight * pxwd9[1,0] + 
                             EastWestWeight * pxwd9[1,2] + 
                             0.025 * pxwd9[2,0] + 
                             NorthSouthWeight * pxwd9[2,1] + 
                             0.025 * pxwd9[2,2]) 
        # A proportion (total weight of water plus desert plus sparse veg 
        # pixels * 2.5) is subtracted from the urban influence coming to 
        # zero if city is fully surrounded by water or desert
        WaterDesertReducedUrbanAverage = 0.0
        if StationWeightedUrbanAverage > 0.0: 
            WaterDesertReducedUrbanAverage = StationWeightedUrbanAverage * \
                                             (1 - 2.5 * WaterDesertWeight)

        station_out = open(OutStation, 'a')
        station_out.write(
            f"{StationID:>12} {StationName:<31} {Latitude:7.3f} "
            f"{Longitude:8.3f} {Elevation:7.1f} {pxlu9[1,1]:4} "
            f"{pxur9[1,1]:1} {StationUrbanAverage:8.2f} "
            f"{StationWeightedUrbanAverage:8.2f} "
            f"{WaterDesertReducedUrbanAverage:8.2f} {pxlu9[0,0]:4} "
            f"{pxlu9[0,1]:4} {pxlu9[0,2]:4} {pxlu9[1,0]:4} "
            f"{pxlu9[1,2]:4} {pxlu9[2,0]:4} {pxlu9[2,1]:4} "
            f"{pxlu9[2,2]:4}\n"
        )
        station_out.close()
