#!/usr/local/sci/bin/python
#***************************************
#
#   Plot uncertainty maps (with coverage completeness)
#
#   23 July 2013 RJHD
#***************************************
#************************************************************************
#                    SVN Info
# $Rev:: 36                                            $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2015-03-27 10:03:26 +0000 (Fri, 27 Mar 2015) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import numpy as np
import datetime as dt
import netCDF4 as ncdf
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import matplotlib as mpl
import matplotlib.patches

import math
import os
from mpl_toolkits.basemap import Basemap
from plt_ts_utils import *


# functions
#***************************************
def MedianPairwiseSlopes(xdata,ydata, mdi):
    '''
    Calculate the median of the pairwise slopes - assumes no missing values

    :param array xdata: x array
    :param array ydata: y array
    :param float mdi: missing data indicator
    :returns: float of slope
    '''
    # run through all pairs, and obtain slope
    slopes=[]
    for i in range(len(xdata)):
        for j in range(i+1,len(xdata)):
            if ydata[i] >  mdi and ydata[j] > mdi:
                slopes+=[(ydata[j]-ydata[i])/(xdata[j]-xdata[i])]

    mpw=np.median(np.array(slopes))

    # copied from median_pairwise.pro methodology (Mark McCarthy)
    slopes.sort()

    good_data=np.where(ydata != mdi)[0]

    n=len(ydata[good_data])

    dof=n*(n-1)/2
    w=math.sqrt(n*(n-1)*((2.*n)+5.)/18.)

    rank_upper=((dof+1.96*w)/2.)+1
    rank_lower=((dof-1.96*w)/2.)+1

    if rank_upper >= len(slopes): rank_upper=len(slopes)-1
    if rank_upper < 0: rank_upper=0
    if rank_lower < 0: rank_lower=0

    upper=slopes[int(rank_upper)]
    lower=slopes[int(rank_lower)]
    

    return mpw,lower,upper      # MedianPairwiseSlopes

#*********************************************
def mad(data, median = False):    
    ''' Calculate the MAD of the data '''
    
    if median:
        mad = np.ma.mean(np.ma.abs(data - np.ma.median(data)))
        
    else:        
        mad = np.ma.mean(np.ma.abs(data - np.ma.mean(data)))

    return mad # mean_absolute_deviation

#***************************************
def DifferencingCalculation(data, verbose=False):
    '''
    Calculate trend for index for each gridbox
    return the range, significance and the detrended data

    :param array data: masked array
    :param bool verbose: extra print to screen
    '''

    nyears,nlat,nlon=data.shape

    # convert YYYYMMDD to YYYY
    years=np.arange(nyears)

    # set up arrays
    difference=np.zeros([nlat,nlon])
    difference.fill(FLTMDI)
    significance=np.zeros([nlat,nlon])


    # run through each grid box
    for t,lat in enumerate(range(nlat)):
        if verbose: print "%i/%i" % (t, nlat)
        for lon in range(nlon):
            # get the timeseries for this gridbox
            gridbox=data[:,lat,lon]
            
            # if there is data
            if np.max(gridbox)!=np.min(gridbox):

                # need to do completeness check with spread through decades

                # first 20 and last 20 average

                first = np.ma.median(gridbox[:20])
                last = np.ma.median(gridbox[-20:])

                difference[lat,lon] = last - first
                
                first_std = mad(gridbox[:20], median = True)
                last_std = mad(gridbox[-20:], median = True)

                if np.abs(difference[lat,lon]) > (first_std + last_std):
                    # standard deviations mean differences don't overlap
                    significance[lat,lon]=1
                else:
                    significance[lat,lon]=0
                
                if verbose:
                    plt.clf()
                    plt.plot(years[np.where(gridbox.mask != True)],gridbox.compressed(),'bo')
                    plt.axhline(y=first, xmax = 20./len(gridbox), color='r')
                    plt.axhline(y=first+first_std, xmax = 20./len(gridbox), color='r', ls=':')
                    plt.axhline(y=first-first_std, xmax = 20./len(gridbox), color='r', ls=':')
                    plt.axhline(y=last, xmin = (len(gridbox)-20.)/len(gridbox), color='r')
                    plt.axhline(y=last+last_std, xmin = (len(gridbox)-20.)/len(gridbox), color='r', ls=':')
                    plt.axhline(y=last-last_std, xmin = (len(gridbox)-20.)/len(gridbox), color='r', ls=':')
                    plt.xlim(0, len(gridbox))
                    plt.show()
                         

    return difference, significance # Differencing Calculation

#***************************************
def TrendingCalculation(data,verbose=False):
    '''
    Calculate trend for index for each gridbox
    return the range, significance and the detrended data

    :param array data: masked array
    :param bool verbose: extra print to screen
    '''

    nyears,nlat,nlon=data.shape

    # convert YYYYMMDD to YYYY
    years=np.arange(nyears)

    # set up arrays
    trend=np.zeros([nlat,nlon])
    trend.fill(FLTMDI)
    sigma=np.zeros([nlat,nlon])
    sigma.fill(FLTMDI)
    significance=np.zeros([nlat,nlon])

    detrended_data=np.zeros(data.shape)
    detrended_data.fill(FLTMDI)

    # run through each grid box
    for t,lat in enumerate(range(nlat)):
        if verbose: print "%i/%i" % (t, nlat)
        for lon in range(nlon):
            # get the timeseries for this gridbox
            gridbox=data[:,lat,lon]
            
            # if there is data
            if np.max(gridbox)!=np.min(gridbox):

                # need to do completeness check with spread through decades

                # median pairwise slopes
                mpw,mpw_l,mpw_u=MedianPairwiseSlopes(years,gridbox.filled(FLTMDI),FLTMDI)


                # get trend, sigma and significance
                trend[lat,lon]=mpw

                sigma[lat,lon]=(mpw_u-mpw_l)/2./2.  # 95% ~ 2 sigma full width

                if ((mpw > 0) and (mpw_l > 0)) or ((mpw < 0 and mpw_u < 0)):
                    # "significant" as both slope is clearly different from zero
                    significance[lat,lon]=1
                else:
                    significance[lat,lon]=0

                # detrend

                # find point on line
                ypivot=np.mean(gridbox.compressed())
                xpivot=np.mean(years[np.where(gridbox.mask != True)])

                # get linear trend values at each year
                lineary=[mpw*(yr-xpivot)+ypivot for yr in years[np.where(gridbox.mask != True)]]


                detrended=gridbox.compressed()-lineary
                
                if False:
                    plt.clf()
                    plt.plot(years[np.where(gridbox.mask != True)],gridbox.compressed(),'bo')
                    plt.plot(years[np.where(gridbox.mask != True)],lineary, 'r-')
                    plt.show()
                         

                detrended_data[np.where(gridbox.mask != True),lat,lon]=detrended

    return trend, sigma, significance, detrended_data # Trending Calculation


#***************************************
def FixVersionCompleteness(true_version_completeness,compl_text):
    '''
    scale the true version completeness to the plotted one

    :param array true_version_completeness: array of completeness
    :param list compl_text: list of str mapping to use
    '''
    plot_version_completeness=np.zeros(true_version_completeness.shape)
    plot_version_completeness.fill(-99) # use as missing

    locs=np.where(true_version_completeness!=0)
    plot_version_completeness[locs]=-88 # use as only 1 model 

    if len(compl_text[0])==len(compl_text[1])==1:
        for cr,compl_value in enumerate(compl_text):
            # extract the value
            locs=np.where(true_version_completeness==int(compl_value))

            plot_version_completeness[locs]=cr

    else:

        for cr,compl_range in enumerate(compl_text):
            # extract the range.
            compl_range=compl_range.split("-")
            if len(compl_range)==2:
                lower=int(compl_range[0])
                upper=int(compl_range[1])
            
                for value in range(lower,upper+1):
                    locs=np.where(true_version_completeness==value)
                    plot_version_completeness[locs]=cr


            else:
                locs=np.where(true_version_completeness==int(compl_range[0]))
                plot_version_completeness[locs]=cr
                
            
    return plot_version_completeness # FixVersionCompleteness

#***************************************
def PlotMap(lons,lats,index,version_completeness,mean_correlations,title,compl_text,
            vmin=0.1,
            vmax=1.001,
            log=False,
            extend='neither',
            nsteps=10,
            reverse=False,
            significance=[]):
    '''
    Plot the maps with significance shown by grey highlighting around the boxes

    :param array lons: longitude array
    :param array lats: latitude array
    :param str index: name of index
    :param array version_completeness: how many of the methods had a value at this grid box
    :param array mean_correlations: average correlation of the boxes
    :param str title: title for plot
    :param str compl_text: text for completeness axis of colourbar
    :param float vmin: min value for colourbar
    :param float vmax: max value for colourbar
    :param boolean log: log scale for colourbar
    :param str extend: extend colourbar option
    :param int nsteps: number of bits in colourbar
    :param boolean reverse: to reverse the colourscale
    :param array significance: array for significant boxes
    '''
    fig=plt.figure(1,figsize=(10,6))
    ax=fig.add_axes([0.05,0.15,0.90,0.90])

    m = Basemap(projection='robin',lon_0=0,resolution='l',ax=ax)
    m.drawcoastlines(color='black', linewidth=.1)

    plt.title(index)

    Lons,Lats=np.meshgrid(np.append(lons,180),lats)
    pltlons,pltlats=m(Lons, Lats)
    pltlons,pltlats=np.array(pltlons),np.array(pltlats) 
 
    # roll through specified colour maps
    #  no point doing case of only 1 version fills boxes - hence the "+2" below and n-1 maps.
    # allow reversal of colormap 
    if reverse:
        cmaps=[plt.cm.Greens_r,plt.cm.Blues_r,plt.cm.Purples_r,plt.cm.Reds_r]
    else:
        cmaps=[plt.cm.Greens,plt.cm.Blues,plt.cm.Purples,plt.cm.Reds]

    # purples are add in - remove if only need three bins
    if len(compl_text)==3:
        cmaps=[cmaps[0],cmaps[1],cmaps[3]]

    for cm,cmap in enumerate(cmaps): 
      
        # discretise colormap - http://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
        cmaplist=[cmap(i) for i in range(cmap.N)]

        if cm!=plt.cm.Reds and reverse:
            cmaplist=cmaplist[56:] #  reduce darkness of the green/blue/purple
        elif cm!=plt.cm.Reds:
            cmaplist=cmaplist[:200] #  reduce darkness of the green/blue/purple
        cmap=cmap.from_list('this_cmap',cmaplist, cmap.N)  
     
        # allow for logarithmic colorbar
        if log:
            bounds=np.logspace(np.log10(vmin),np.log10(vmax),nsteps)
            strbounds=["%4.2g" % i for i in bounds]
            norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N) # LogNorm(vmin,vmax)
        else:
            bounds=np.linspace(vmin,vmax,nsteps)
            strbounds=["%4.2g" % i for i in bounds]
            norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)

        # plot the uncertainty map
        # note - masking - so mask all bar value of interest
        plot_mean_correlations=np.ma.masked_where(version_completeness != cm,mean_correlations)
        grids=ax.pcolor(pltlons,pltlats, plot_mean_correlations[:-1,:], cmap=cmap,norm=norm)

        if significance != []:
            plot_mean_correlations=np.ma.masked_where(significance == 0,plot_mean_correlations)
            grids=ax.pcolor(pltlons,pltlats, plot_mean_correlations[:-1,:], cmap=cmap,norm=norm, edgecolor='0.2',linewidth=0.5)

        cbax=fig.add_axes([0.15,0.07+(cm*0.03),0.70,0.03])
        cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds, extend=extend)

        if log:
            strbounds[0]=r'$\leq$'+strbounds[0]
            strbounds[-1]=r'$\geq$'+strbounds[-1]
        cb.ax.set_xticklabels(strbounds)

        for t in cb.ax.get_xticklabels():
            if cm==0:
                t.set_fontsize(7)
            else:
                t.set_fontsize(0)

        if cm==len(cmaps)-1:
            plt.title(title)

        cbax.text(1.02,0.5,compl_text[cm],va='center',size=8)
        cbax.text(-0.02,0.5,compl_text[cm],va='center',ha='right',size=8)

    plt.figtext(0.95,0.12,"realisations",size=9,ha='center')
    plt.figtext(0.95,0.10,"with",size=9,ha='center')
    plt.figtext(0.95,0.08,"coverage",size=9,ha='center')

    
    # finish off the map
    m.drawmeridians(np.arange(-180, 180, 60), color='#bbbbbb',labels=[1,0,0,1])
    m.drawparallels(np.arange(-90, 90, 30), color='#bbbbbb',labels=[1,0,0,1])
    m.drawmapboundary(zorder=10, fill_color='none')

    return # PlotMap


#***************************************
# fixed variables

DATANAME='JackKnife_'

COMPLETENESS=['2-34','35-67','68-100']

THRESHOLDS=['25']#,'50','75']
ITERATIONS=100

DATALOCATION='/data/local/rdunn/HadEX2/hadex2/gridded_2013/'
TSLOCATION='/data/local/rdunn/HadEX2/hadex2/timeseries_2013/'

PLOTLOCATION='/data/local/rdunn/HadEX2/hadex2/images_2013/'
MASKFILE='/data/local/rdunn/HadEX2/hadex2/h2_mask.msk'
INTMDI=-99.9
FLTMDI=-1.e30

INDICES=['CDD','CSDI','CWD','DTR','ETR','FD','GSL','ID','PRCPTOT','R10mm','R20mm','R95p','R95pTOT','R99p','R99pTOT','Rx1day','Rx5day','SDII','SU','TN90p','TN10p','TNx','TNn','TX90p','TX10p','TXx','TXn','WSDI','TR']
#INDICES = ['TX90p','TX10p','TN90p','TN10p','SDII','ETR']
STARTYEAR='1901'
ENDYEAR='2010'
STARTTIME=1951 # only last 60 years to increase coverage
ENDTIME=2010

#INDICES=['TR']

for thd,threshold in enumerate(THRESHOLDS):


    for idx,index in enumerate(INDICES):
        print index

        # need to sort magic numbers for lat/lon
        # set up blank arrays
        all_data=np.ma.zeros([ITERATIONS,ENDTIME-STARTTIME+1,73,96])
        all_trends=np.ma.zeros([ITERATIONS,73,96])
        all_completeness_masks=np.ma.zeros([ITERATIONS,73,96])
        all_differences=np.ma.zeros([ITERATIONS,73,96])

        # read in each file
        for iteration in range(ITERATIONS+1):

            if iteration==0:
                infilename='HadEX2_'+index+'_1901-2010_h2_mask_m4.nc'
            else:
                infilename=DATANAME+index+'_1901-2010_h2_mask_m4_'+threshold+'_'+str(iteration)+'.nc'

            ANNUALFILE=DATALOCATION+infilename

            month='Ann'

            # read the data
            try:
                ncfile=ncdf.Dataset(ANNUALFILE,'r')

                timevar=ncfile.variables['time'] # array of YYYMMDD
                latvar=ncfile.variables['lat'] # array of lats
                lonvar=ncfile.variables['lon'] # array of lons

                annualvar=ncfile.variables[month] # array of arrays

                times=timevar[:]
                lats=latvar[:]
                lons=lonvar[:]

                # restrict to times of interest
                years=np.array([int(str(y)[0:4]) for y in times])
                start_loc=np.where((years >= STARTTIME) & (years <= ENDTIME))[0]
                times=times[start_loc[0]:start_loc[-1]+1]

                annual_data=annualvar[:]
                annual_data=annual_data[start_loc[0]:start_loc[-1]+1,:,:]
            except IOError,RuntimeError:
                print "File not found: ",filename
                raise IOError,RuntimeError


            # have to test for nans here, and mask them out
            #nan_loc=np.where(np.isnan(annual_data.data) == True)
            #annual_data.mask[nan_loc]=True

            # adjust longitudes and latitudes by half

            lons=lons-(lons[1]-lons[0])/2.
            lats=lats-(lats[1]-lats[0])/2.


            # roll longitudes to be from 0->360 to -180->180

            lons=np.roll(lons, -(len(lons)/2))
            west=np.where((lons >180) & (lons <= 360))
            lons[west]=lons[west]-360
            lons[0]=-lons[0] # fix the 180deg

            # now roll data
            annual_data=np.roll(annual_data, -(len(lons)/2), axis=2)
            # this changes the missing value to -1e20, so fix mdi value
            mdi=annual_data.fill_value

            # find which boxes have x% of years with data
            yearly_completeness_mask=CompletenessCheck(annual_data, times, mdi, ENDTIME, STARTTIME, threshold=0.66)


            # find trends (MPW) & return detrended data
            annual_trend, annual_sigma, annual_signif, annual_detrended=TrendingCalculation(annual_data)

            # do difference between 1991-2010 and 1951-1970
            annual_difference, annual_diff_signif=DifferencingCalculation(annual_data)

            if iteration==0:
                significance=annual_signif
                diff_signif = annual_diff_signif
            else:
                all_trends[iteration-1,:,:]=annual_trend
                all_data[iteration-1,:,:,:]=annual_detrended # annual_data
                all_completeness_masks[iteration-1,:,:]=yearly_completeness_mask
                all_differences[iteration-1,:,:]=annual_difference
            


            print ANNUALFILE #, np.max(annual_data) #  status monitor

        # gone through all files
        all_data=all_data.filled(FLTMDI)
        all_completeness_masks=np.array(all_completeness_masks)

        # completeness across versions
        true_version_completeness=np.zeros(yearly_completeness_mask.shape)
        true_version_completeness[:,:]=ITERATIONS
        true_version_completeness=true_version_completeness-np.sum(all_completeness_masks,axis=0)
        # contains number of times a grid box is filled across all versions with the

        version_completeness=FixVersionCompleteness(true_version_completeness, COMPLETENESS)

        #*************************
        # trend ranges
        #*************************

        all_trends=np.ma.masked_where(all_trends == FLTMDI, all_trends)
        trend_range=np.abs(np.std(all_trends,axis=0)/np.mean(all_trends,axis=0))
        trend_range=np.ma.masked_where(true_version_completeness == 0,trend_range)
        trend_range=np.ma.masked_where(version_completeness == -99,trend_range)

        difference_range =np.abs(np.std(all_differences,axis=0)/np.mean(all_differences,axis=0))
        difference_range=np.ma.masked_where(true_version_completeness == 0,difference_range)
        difference_range=np.ma.masked_where(version_completeness == -99,difference_range)

        #*************************
        # start plot - relative size of trend
        #*************************

        plt.clf()
        # remove zero's from trend_range - as logarithmic scale can't plot them
        vmin=0.01
        vmax = 10.
        PlotMap(lons,lats,index,version_completeness,trend_range,r'$\sigma$(trends)/$\mu$(trends)',COMPLETENESS,vmin=vmin,vmax=vmax,reverse=True,log=True, significance=significance)

        watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
        #plt.figtext(0.01,0.01,watermarkstring,size=6)

        # save figure
        plt.savefig(PLOTLOCATION+"jack_"+index+"_normalised-trends.png", dpi=150)
        plt.savefig(PLOTLOCATION+"jack_"+index+"_normalised-trends.eps")

	plt.close()

        #*************************
        # start plot - relative size of differences between early and late
        #*************************

        plt.clf()
        # remove zero's from trend_range - as logarithmic scale can't plot them
        vmin=0.01
        vmax=10.
        small=np.where(difference_range < vmin)
        trend_range[small]=vmin/10.
        PlotMap(lons,lats,index,version_completeness,difference_range,r'$\sigma$(differences)/$\mu$(differences)',COMPLETENESS,vmin=vmin,vmax=vmax,reverse=True,log=True, significance=diff_signif)

        watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
        #plt.figtext(0.01,0.01,watermarkstring,size=6)
        plt.figtext(0.01,0.95,'(c)')

        # save figure
        plt.savefig(PLOTLOCATION+"jack_"+index+"_normalised-differences.png", dpi=150)
        plt.savefig(PLOTLOCATION+"jack_"+index+"_normalised-differences.eps")
	plt.close()
 #************************************************************************
#                                 END
#************************************************************************
