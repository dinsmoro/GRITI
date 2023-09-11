#GOAL: Plot only  keogram
#RD on 3/01/21
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from Code.subfun_figFitter import figFitter
import astropy.coordinates as coord
from astropy.time import Time, conf
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.subfun_sunAlsoRises_location import sunAlsoRises_location
    

def GRITI_keo_plot_wSun(keoData,timeUnique,plotLimValu,colorMap,plotLatRange,plotLongRange,latMillstone,longMillstone,
        dates,avg_anyAngle,avg_anyAngle_Width,avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name,
        avg_anyAngle_dataType,avg_anyAngle_plotLabel,FONT_titleFM,FONT_axisTickFM,FONT_axisLabelFM):
    
    #-----Unpack-----
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum']; #unpack
    dateRange_dayNum_full = dates['date range full dayNum']; #unpack
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(plotLimValu) == 1 ):
        plotLimValu = np.array( (-plotLimValu,plotLimValu) ); #make it a vector
    #END IF

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ):
        fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        divider = make_axes_locatable(ax); #prep to add an axis
        cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        ax = [ax]; #wrap it into a list to make things work
    else:
        import matplotlib.gridspec as gridspec
        fig = plt.figure();
        gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
        gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=20, ncols=1, subplot_spec = gridr[0], hspace = 5.0)
        gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
        # gridr.update(hspace=0.05); # set the spacing between axes.
        fig.add_subplot(gridr1[:]); #RTI plots are 2 tall
        # gridr.update(hspace=0.80); # set the spacing between axes.
        fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
        ax = fig.axes; #get a list of the axes
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        divider = make_axes_locatable(ax[0]); #prep to add an axis
        cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        divider2 = make_axes_locatable(ax[1]); #prep to add an axis
        cax2 = divider2.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        cax2.set_visible(False); #make it sneaky
    #END IF
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    for i in range(0,len(ax)):
        ax[i].set_aspect('auto'); #remove aspect ratios
    #END IF
    
    #[(timeUnique - dateRange_dayNum_zeroHr[1])*24 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (np.append(timeUnique,timeUnique[-1]+np.median(np.diff(timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax[0].pcolormesh(pltHelprX, pltHelprY,  keoData.T ,vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(FONT_axisTickFM); #yee
        #END FOR tick
        #cbar.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #they changed how the code works, this doesn't work anymore
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #now it's this
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    else:
        im = ax[0].pcolormesh(pltHelprX, pltHelprY,  keoData.T ,cmap=colorMap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(FONT_axisTickFM); #yee
        #END FOR tick
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #END IF
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    string_Title = avg_anyAngle_dataType+' Avg\'d on '+str(np.round(avg_anyAngle,2))+' deg Angle w/ '+ \
        str(np.round(avg_anyAngle_Width,2))+' arcdeg Width (Sun Location is Orange Line)'; #create mecha title
    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[0].set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' [hr]',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    # xAxisTicks = np.arange( (np.round((np.min(timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),2)) , \
    #         (np.round((np.max(timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),2)) + 4 , \
    #         4); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    # ax[0].set_xlim( (xAxisTicks[0],xAxisTicks[-1]) ); #set x axis ticks
    # ax[0].set_xticks(xAxisTicks); #set x axis ticks
    xAxisTicks, xAxisLims = GRITI_plotHelper_axisizerTime((timeUnique-dateRange_dayNum_zeroHr[1]*86400)/3600,ax=ax[0]); #automagic time ticks here
    
    avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    else:
        if(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Latitude'): #if Y axis is latitude, use latitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/13; #just goes for it if it's a super tiny range
        elif(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude'): #if Y axis is longitude, use longitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/13; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))+avg_anyAngle_Range_Chunks_Long_Plot_autoTick,avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax[0].set_ylim( (np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))) ); #set x axis ticks
    ax[0].set_yticks(yAxisTicks); #set x axis ticks
    
    
    #Now drawing line of interest
    # if( (settings_map['coord type'] == 'geo') & ('keo coord type' in settings_dataSpecific) ): #pull a switcheroo if needed
    #     if( settings_dataSpecific['keo coord type'] == 'mag' ):
    #         import copy
    #         import datetime
    #         from Code.subfun_convertToMag import convert_to_mag
            
    #         if( 'pierceAlt' in settings_dataSpecific ):
    #             keo_alt = settings_dataSpecific['pierceAlt'];
    #         else:
    #             keo_alt = 120.; #default, great for auroral zone stuff (like field aligned currents)
    #         #END IF
    #         timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
    #         kk = np.where(np.int64(data_timeUnique[timeIndex]/86400) == dates['date range full dayNum'][:,1])[0].item(); #get where the year is gonna be
    #         time4mag_hr = np.int32(np.mod(data_timeUnique[timeIndex],86400)//3600); #get hours
    #         time4mag_min = np.int32(np.mod(data_timeUnique[timeIndex],86400)//60-time4mag_hr*60); #get the minutes
    #         time4mag_sec = np.int32(np.mod(data_timeUnique[timeIndex],86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
    #         time4mag = datetime.datetime(dates['date range full'][kk,0],dates['date range full'][kk,1],dates['date range full'][kk,2], \
    #                                      hour = time4mag_hr, minute = time4mag_min, second = time4mag_sec); #date time object for aacgmv2    
            
    #         # avgPt_coords = copy.deepcopy(avgPt_coords); #copy big time
    #         [latMillstone, longMillstone] = convert_to_mag(latMillstone, longMillstone, keo_alt, time4mag); #convert
    #     #END IF
    # #END IF
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
            ax[0].plot( np.linspace(np.min((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=1); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
            ax[0].plot( np.linspace(np.min((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=1); #plots a point with a black line
        #END IF
    #END IF
    
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ):
        #Now draw sun position
        #yes I just abs wing making these solutions every time hoping to find the codegrail each time, I don't tho
        #!NO YEAR SUPPORT YET!
        #time to make some strings of time
        timeString = [[] for i in range(0,timeUnique.size)]; #preallocate list to hold the strings
        timeTime = [[] for i in range(0,timeUnique.size)]; #preallocate list to hold the strings
        #---Build time things needed---
        timeMonDay = subfun_dayNum_to_date(np.array([np.ones(timeUnique.size,dtype=np.int64)*dateRange_dayNum_full[0,0],np.int64(timeUnique/86400)]).T); #get the month and day
        timeTempHrFloat = (timeUnique/86400 - np.int64(timeUnique/86400))*24; #hr, get the decimal bit
        timeTempHr = np.int32(timeTempHrFloat); #hr, get the hours
        timeTempMinFloat = (timeTempHrFloat - timeTempHr)*60; #min, get the decimal bit
        timeTempMin = np.int32(timeTempMinFloat); #min, get the minutes
        timeTempSec = np.int32((timeTempMinFloat - timeTempMin)*60); #sec, get the sec (we won't go sub-sec)
        #---fix float32 time issues (we know we only operate on whole seconds)---
        timeTempMin[timeTempSec == 59] += 1; #increment
        timeTempSec[timeTempSec == 59] = 0; #reset to 0
        timeTempHr[timeTempMin == 60] += 1; #increment
        timeTempMin[timeTempMin == 60] = 0; #reset to 0
        timeTempHr[timeTempHr == 24] = 0; #reset to 0
        for i in range(0,timeUnique.size):
            timeString[i] = str(dateRange_dayNum_full[0,0])+'-'+str(timeMonDay[i,1]).zfill(2)+'-'+str(timeMonDay[i,2]).zfill(2)+'T'+ \
                str(timeTempHr[i]).zfill(2)+':'+str(timeTempMin[i]).zfill(2)+':'+str(timeTempSec[i]).zfill(2)+'.000'; #make the strings
        #END FOR i
        for i in range(0,timeUnique.size):
            #uses ISOT format 2000-01-01T00:00:00.000
            with conf.set_temp('use_fast_parser', 'force'):
                try:
                    timeTime[i] = Time(timeString[i], format='isot', scale='utc'); #convert
                except ValueError as errorz:
                    print(errorz); #that's a whoopsie (should never happen)
                #END TRY
            #END WITH
        #END FOR i
        sunPlot = np.empty(timeUnique.size,dtype=np.float32); #preallocate
        for i in range(0,timeUnique.size):
            sun = coord.EarthLocation.from_geocentric(coord.get_sun(timeTime[i]).transform_to(coord.ITRS).x,coord.get_sun(timeTime[i]).transform_to(coord.ITRS).y,coord.get_sun(timeTime[i]).transform_to(coord.ITRS).z); #trainwreck of a call to get the sun in lat/long
            if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
                sunPlot[i] = sun.lon.value; #get that data
            else:
                sunPlot[i] = sun.lat.value; #get that data
            #END IF
        #END FOR i
        # ax[0].plot( (timeUnique - dateRange_dayNum_zeroHr[1])*24, sunPlot,c='xkcd:orange',linewidth=2.5); #plots a point with an orange line
        indexes2plot = np.hstack( (np.arange(0,timeUnique.size,10),(timeUnique.size-1,)) ); #plot select stars
        ax[0].plot( ((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600)[indexes2plot], sunPlot[indexes2plot], c='xkcd:orange', marker='*', markersize=20, linestyle='None'); #plots a point with an orange line
        # ax[0].plot( np.linspace(np.min((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
        #         np.tile(-72.68,2) , #Y latitude OR longitude arcdeg
        #         c='xkcd:black',linewidth=1,linestyle='--'); #plots a point with a black line
    else:
        sunSubSolar_loc = sunAlsoRises_location(dateRange_dayNum_full,timeIndexes=timeUnique); #calc sun locations
        indexes2plot = np.hstack( (np.arange(0,timeUnique.size,10),(timeUnique.size-1,)) ); #plot select stars
        ax[1].plot( ((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600)[indexes2plot], sunSubSolar_loc['long'][indexes2plot], c='xkcd:orange', marker='*', markersize=10, linestyle='None'); #plots a point with an orange line
        ax[1].set_yticklabels([]); #remove y labels for day nite plot
        ax[1].set_yticks([]); #remove y tick marks for day nite plot
        # ax[1].set_ylim( (0,1) ); #set y lims to 0 and 1
        # ax[1].set_ylabel('Subsolar Long. [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
        ax[1].set_xticks(xAxisTicks); #set x axis ticks
        ax[1].set_xlim( xAxisLims ); #set x axis limits
        ax[1].set_xticklabels([]); #remove y labels for day nite plot
        ax[1].tick_params(axis='x',direction='in'); #ticks in
        # ax[1].spines['left'].set_visible(False); #turn off box lines
        # ax[1].spines['right'].set_visible(False); #turn off box lines
        # ax[1].spines['top'].set_visible(False); #turn off box lines
        # ax[1].plot( np.linspace(np.min((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
        #         np.tile(-72.68,2) , #Y latitude OR longitude arcdeg
        #         c='xkcd:black',linewidth=1,linestyle='--'); #plots a point with a black line
    #END IF    
    
    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.050, right = 0.945, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
       
#END DEF no return, plot is the return!




