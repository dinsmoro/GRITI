#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from subfun_figFitter import figFitter

def GRITI_keo_plot(keo, data_timeUnique, data_timeRef, dates, \
        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
        FLG_fancyPlot = 0):
    
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
            
    #--- Unpack ---
    PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
    FONT_titleFM = settings_plot['font title FM'];
    FONT_axisTick = settings_plot['font axis tick'];
    FONT_axisLabelFM = settings_plot['font axis label FM'];
    journal_dpi = settings_plot['journal dpi'];
    plotLongRange = settings_map['long range'];
    plotLatRange = settings_map['lat range'];
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum'];
    dateRange_zeroHr = dates['date range zero hr'];
    dateRange_zeroHr_monthName = dates['date range zero hr month name'];
    dateRange_zeroHr_dayPostfix = dates['date range zero hr day post fix'];
    plotLimValu = settings_dataSpecific['keo plot lim'];
    avgPt_coords = settings_map['site coords'];
    plotLatRange = settings_map['lat range']; #unpack
    plotLongRange = settings_map['long range']; #unpack
    plotLimValu = settings_dataSpecific['keo plot lim']; #unpack
    keo_angle = settings_dataSpecific['keo angle'];
    keo_width = settings_dataSpecific['keo width'];
    # keo_N = settings_dataSpecific['keo N'];
    keo_name = settings_dataSpecific['keo name'];
    keo_range_plotLatLong_name = settings_dataSpecific['keo plot latlong name'];
    keo_range_plotLatLong_chunks = settings_dataSpecific['keo plot latlong chunks'];
    keo_colormap = settings_dataSpecific['keo colormap'];
    keo_plotLabel = keo_name+settings_dataSpecific['keo units'];
            
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(plotLimValu) == 1 ):
        plotLimValu = np.array( (-plotLimValu,plotLimValu) ); #make it a vector
    #END IF

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    if( FLG_fancyPlot == 0 ):
        fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    #END IF
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    
    #[(data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , keo_range_plotLatLong_chunks] ,
    pltHelprX, pltHelprY = np.meshgrid( (data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                keo_range_plotLatLong_chunks);
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax.pcolormesh(pltHelprX, pltHelprY,  keo.T ,vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),cmap=keo_colormap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        if( FLG_fancyPlot == 0 ):
            cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
            cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
        else:
            cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),11)); #create useful tick marks
            cax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
        #END IF
        cbar.set_label(keo_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        #cbar.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #they changed how the code works, this doesn't work anymore
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #now it's this
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    else:
        im = ax.pcolormesh(pltHelprX, pltHelprY,  keo.T ,cmap=keo_colormap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(keo_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #END IF
    
    if( FLG_fancyPlot == 0 ): #only title non-fancy plot
        #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
        #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
        #        ' arcdeg, Line Shows '+keo_range_plotLatLong_name+' of Millstone Hill Zenith Beam'; #create mecha title
        string_Title = keo_name+' Averaged on Angle of '+str(np.round(keo_angle,2)).rstrip('0').rstrip('.')+' deg and Width of '+ \
            str(np.round(keo_width,2)).rstrip('0').rstrip('.')+' arcdeg'; #create mecha title
        ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    #END IF
    ax.set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax.set_ylabel(keo_range_plotLatLong_name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    #makes plotting look better
    xAxisLims = np.array(((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600)); #get those xAxis limits, ensures everything is aligned
    
    if( (np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[0] = np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    if( (np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[1] = np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    # ax.set_xlim(xAxisLims); #set the xlims now
    yAxisLims = np.array(ax.get_ylim()); #get the current y axis limits
    if( (np.round(np.min(keo_range_plotLatLong_chunks)) - np.min(keo_range_plotLatLong_chunks)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[0] = np.round(np.min(keo_range_plotLatLong_chunks)); #force the rounded value
    #END IF
    if( (np.round(np.max(keo_range_plotLatLong_chunks)) - np.max(keo_range_plotLatLong_chunks)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[1] = np.round(np.max(keo_range_plotLatLong_chunks)); #force the rounded value
    #END IF
    ax.set_ylim(yAxisLims); #set the ylims now
    
    if( FLG_fancyPlot == 0 ): #different aspect ratios require different spacing assumptions
        keo_range_plotLatLong_chunks_autoTick = (np.ceil(np.max(keo_range_plotLatLong_chunks)) - np.floor(np.min(keo_range_plotLatLong_chunks)))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( keo_range_plotLatLong_chunks_autoTick > 25 ):
            keo_range_plotLatLong_chunks_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 10 ):
            keo_range_plotLatLong_chunks_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 5 ):
            keo_range_plotLatLong_chunks_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 2 ):
            keo_range_plotLatLong_chunks_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 1 ):
            keo_range_plotLatLong_chunks_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            keo_range_plotLatLong_chunks_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
        else:
            if(keo_range_plotLatLong_name == 'Latitude'): #if Y axis is latitude, use latitude
                keo_range_plotLatLong_chunks_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/13; #just goes for it if it's a super tiny range
            elif(keo_range_plotLatLong_name == 'Longitude'): #if Y axis is longitude, use longitude
                keo_range_plotLatLong_chunks_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/13; #just goes for it if it's a super tiny range
            #END IF
        #END IF
    else:
        keo_range_plotLatLong_chunks_autoTick = (np.ceil(np.max(keo_range_plotLatLong_chunks)) - np.floor(np.min(keo_range_plotLatLong_chunks)))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( keo_range_plotLatLong_chunks_autoTick > 20 ):
            keo_range_plotLatLong_chunks_autoTick = 30; #sets the tick setting to 30 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 10 ):
            keo_range_plotLatLong_chunks_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 5 ):
            keo_range_plotLatLong_chunks_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 2 ):
            keo_range_plotLatLong_chunks_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick > 1 ):
            keo_range_plotLatLong_chunks_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
        elif( keo_range_plotLatLong_chunks_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            keo_range_plotLatLong_chunks_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
        else:
            if(keo_range_plotLatLong_name == 'Latitude'): #if Y axis is latitude, use latitude
                keo_range_plotLatLong_chunks_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/17; #just goes for it if it's a super tiny range
            elif(keo_range_plotLatLong_name == 'Longitude'): #if Y axis is longitude, use longitude
                keo_range_plotLatLong_chunks_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/17; #just goes for it if it's a super tiny range
            #END IF
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(keo_range_plotLatLong_chunks)),np.ceil(np.max(keo_range_plotLatLong_chunks)),keo_range_plotLatLong_chunks_autoTick ),2); #creates y ticks automagically
    ax.set_yticks(yAxisTicks); #set x axis ticks
    ax.set_ylim(yAxisLims); #set the ylims now [gotta do it after ticks in case they're a little too far]
    
    time_autoTick = (np.ceil((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.floor((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( time_autoTick > 24 ):
        time_autoTick = 24; #sets the tick setting to 15 arcdegrees per tick
    elif( time_autoTick > 12 ):
        time_autoTick = 12; #sets the tick setting to 15 arcdegrees per tick
    elif( time_autoTick > 6 ):
        time_autoTick = 6; #sets the tick setting to 10 arcdegrees per tick
    elif( time_autoTick > 4 ):
        time_autoTick = 4; #sets the tick setting to 5 arcdegrees per tick
    elif( time_autoTick > 2 ):
        time_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( time_autoTick >= 1 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        time_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    else:
        time_autoTick = ((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600 - (np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600)/13; #just goes for it if it's a super tiny range
    #END IF
    xAxisTicks = np.arange( (np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) , \
            (np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) + time_autoTick , \
            time_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax.set_xticks(xAxisTicks); #set x axis ticks
    ax.set_xlim(xAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    
    #Now drawing line of interest
    if( keo_range_plotLatLong_name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= avgPt_coords[0,1]) & (np.max(plotLongRange) >= avgPt_coords[0,1]) ): #only plot if it's in the long range specified
            ax.plot( np.linspace(np.min((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
                    np.ones(2)*avgPt_coords[0,1] , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmoller); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= avgPt_coords[0,0]) & (np.max(plotLatRange) >= avgPt_coords[0,0]) ): #only plot if it's in the lat range specified
            ax.plot( np.linspace(np.min((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
                    np.ones(2)*avgPt_coords[0,0] , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmoller); #plots a point with a black line
        #END IF
    #END IF
    
    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.050, right = 0.945, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    if( FLG_fancyPlot == 0 ):
        plt.show(); #req to make plot show up
    else:
        fig.savefig(settings_paths['fancyPlots']+'\\'+settings_dataSpecific['keo data type']+'_keo.png'); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
       
#no return, plot is the return!




