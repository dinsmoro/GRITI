#GOAL: Plot only  keogram
#RD on 3/02/21
#
#INPUT: buncha stuff
#OUTPUT: rolled keo array

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from subfun_figFitter import figFitter
import astropy.coordinates as coord
from astropy.time import Time, conf
from subfun_dayNum_to_date import subfun_dayNum_to_date

def GRITI_keo_plot_sunCentered(keoData,timeUnique,plotLimValu,colorMap,plotLatRange,plotLongRange,latMillstone,longMillstone,
        dates,avg_anyAngle,avg_anyAngle_Width,avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name,
        avg_anyAngle_dataType,avg_anyAngle_plotLabel,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):
    
    #-----Unpack-----
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum']; #unpack
    dateRange_dayNum_full = dates['date range full dayNum']; #unpack
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(plotLimValu) == 1 ):
        plotLimValu = np.array( (-plotLimValu,plotLimValu) ); #make it a vector
    #END IF
    
    #-----Get the sun times-----
    #!NO YEAR SUPPORT YET!
    #time to make some strings of time
    timeString = [[] for i in range(0,timeUnique.size)]; #preallocate list to hold the strings
    timeTime = [[] for i in range(0,timeUnique.size)]; #preallocate list to hold the strings
    #---Build time things needed---
    timeMonDay = subfun_dayNum_to_date(np.array([np.ones(timeUnique.size,dtype=np.int64)*dateRange_dayNum_full[0,0],np.int64(timeUnique)]).T); #get the month and day
    timeTempHrFloat = (timeUnique - np.int64(timeUnique))*24; #hr, get the decimal bit
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

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    
    #[(timeUnique - dateRange_dayNum_zeroHr[1])*24 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (timeUnique - dateRange_dayNum_zeroHr[1])*24, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    #---now let's roll this so the sun is at the center, always---
    goalIndex = np.where( np.abs(avg_anyAngle_Range_Chunks_Long_Plot-0) == np.min(np.abs(avg_anyAngle_Range_Chunks_Long_Plot-0)) )[0].item() - 1; #set the goal index (closet to 0)
    midDegs = np.diff(avg_anyAngle_Range_Chunks_Long_Plot)/2 + avg_anyAngle_Range_Chunks_Long_Plot[:-1]; #degc, middle degrees
    for i in range(0,timeUnique.size):
        indxr = np.where( np.abs(sunPlot[i] - midDegs) == np.min(np.abs(sunPlot[i] - midDegs)) )[0].item(); #get closest index
        indxrDelta = goalIndex - indxr; #get the delta between the index and the goal index
        keoData[i,:] = np.roll(keoData[i,:], indxrDelta); #roll it to the goal index
    #END FOR i
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax.pcolormesh(pltHelprX, pltHelprY,  keoData.T ,vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        #cbar.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #they changed how the code works, this doesn't work anymore
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #now it's this
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    else:
        im = ax.pcolormesh(pltHelprX, pltHelprY,  keoData.T ,cmap=colorMap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #END IF
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    string_Title = avg_anyAngle_dataType+' Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
        str(np.round(avg_anyAngle_Width,2))+' arcdeg (Keo Pixels Under Sun Loc Aligned to Orange Line)'; #create mecha title
    ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax.set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' [hr]',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax.set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(timeUnique)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(timeUnique)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(timeUnique)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(timeUnique)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax.set_xticks(xAxisTicks); #set x axis ticks
    
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
    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)),avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax.set_yticks(yAxisTicks); #set x axis ticks
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= midDegs[goalIndex]) & (np.max(plotLongRange) >= midDegs[goalIndex]) ): #only plot if it's in the long range specified
            ax.plot( np.linspace(np.min((timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),np.max((timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),10,endpoint=True) , #X time hr
                    np.tile(midDegs[goalIndex],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:orange',linewidth=1); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= midDegs[goalIndex]) & (np.max(plotLatRange) >= midDegs[goalIndex]) ): #only plot if it's in the lat range specified
            ax.plot( np.linspace(np.min((timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),np.max((timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),10,endpoint=True) , #X time hr
                    np.tile(midDegs[goalIndex],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:orange',linewidth=1); #plots a point with a black line
        #END IF
    #END IF
    
    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.050, right = 0.945, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
    
    return keoData
#END DEF




