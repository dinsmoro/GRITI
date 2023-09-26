#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
from os.path import join as path_join
from Code.subfun_sunAlsoRises import sunAlsoRises
from Code.subfun_sunAlsoRises_timeZone import sunAlsoRises_timeZone
from Code.subfun_textNice import textNice
from Code.subfun_figFitter import figFitter
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.GRITI_plotHelper_axisizerLatLong import GRITI_plotHelper_axisizerLatLong

def GRITI_keo_plot(keo, data_timeUnique, data_timeRef, dates, \
        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
        timeCutout = None, FLG_drawSun=False, FLG_fancyPlot = 0, highlighter = False, settings_config =  None):
    
    if( FLG_fancyPlot >= 1 ):
        if( np.any(timeCutout != None) ):
            print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo_cutout IN fancyPlot FOLDER'); #report since you won't see anything
        elif( FLG_drawSun == True ):
            print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo_wSun IN fancyPlot FOLDER'); #report since you won't see anything
        else:
            print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo IN fancyPlot FOLDER'); #report since you won't see anything
        #END IF
    #END IF
    letteringPositionX = -0.105; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.88; #set the X position of the lettering (e.g., a. b. c. ...)
            
    #--- Unpack ---
    PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
    FONT_titleFM = settings_plot['font title FM'];
    FONT_axisTickFM = settings_plot['font axis tick FM'];
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
    coordType = settings_map['coord type']; #unpack
    plotLimValu = settings_dataSpecific['keo plot lim']; #unpack
    keo_angle = settings_dataSpecific['keo angle'];
    keo_width = settings_dataSpecific['keo width'];
    # keo_N = settings_dataSpecific['keo N'];
    keo_name = settings_dataSpecific['keo labels'];
    keo_plotLatLong_name = settings_dataSpecific['keo plot latlong name'];
    keo_plotLatLong_chunks = settings_dataSpecific['keo plot latlong chunks'];
    keo_colormap = settings_dataSpecific['keo colormap'];
    keo_plotLabel = keo_name+settings_dataSpecific['keo units'];
    
    if( 'keo coord type' in settings_dataSpecific ):
        coordType = settings_dataSpecific['keo coord type']; #use this if provided
    #END IF
    
    if( 'degree label' in settings_map ):
        latlong_unitName = settings_map['degree label']; #use supplied label
    else:
        latlong_unitName = '°';
    #END IF
    if( latlong_unitName != '' ):
        latlong_unitName_bracketed = settings_plot['unit L']+latlong_unitName+settings_plot['unit R']; #bracket it
    else:
        latlong_unitName_bracketed = ''; #nada
    #END IF
    if( 'indicate direction' in settings_map ):
        if( settings_map['indicate direction'][0] == True ):
            if(keo_plotLatLong_name == 'Latitude'):
                keo_plotLatLong_dirAdder = settings_map['indicate direction'][1]['lat']; #get the lat
            else:
                keo_plotLatLong_dirAdder = settings_map['indicate direction'][1]['long']; #get the long
            #END IF
        else:
            keo_plotLatLong_dirAdder = ''; #nada
        #END IF
    else:
        keo_plotLatLong_dirAdder = ''; #nada
    #END IF
    
    
    if( settings_plot['save file type'].lower() != '.pdf' ):
        plot_rasterize = False; #don't need to bother
    else:
        plot_rasterize = True; #it needs it
    #END IF
            
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(plotLimValu) == 1 ):
        plotLimValu = np.array( (-plotLimValu,plotLimValu) ); #make it a vector
    #END IF
    
    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    if( 'day nite shading' in settings_dataSpecific ):
        if( settings_dataSpecific['day nite shading'] == -1 ):
            figNumRows = 2; #need 2 rows for the day nite line thing
            import timezonefinder
            from urllib.request import urlopen
            from Code.subfun_strstr import strstr
            import pytz
            from datetime import datetime
        else:
            figNumRows = 1;
        #END IF
    else:
        figNumRows = 1;
    #END IF
    if( FLG_fancyPlot == 0 ):
        if( figNumRows == 1 ):
            fig, ax = plt.subplots(nrows=figNumRows, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
            ax = [ax]; #wrap into list for ez later
        else:
            import matplotlib.gridspec as gridspec
            fig = plt.figure();
            gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
            gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[0], hspace = 6.0)
            gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
            # gridr.update(hspace=0.05); # set the spacing between axes.
            fig.add_subplot(gridr1[0]); #keogram is 20 tall
            # gridr.update(hspace=0.80); # set the spacing between axes.
            fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
            ax = fig.axes; #get a list of the axes
            ax[1].set_aspect('auto');
        #END IF
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        if( figNumRows == 1 ):
            fig, ax = plt.subplots(nrows=figNumRows, ncols=1, figsize=(14,8.5), dpi=settings_plot['journal dpi']); #use instead of fig because it inits an axis too (I think I dunno)
            ax = [ax]; #wrap into list for ez later
        else:
            import matplotlib.gridspec as gridspec
            fig = plt.figure(figsize=(14,10.5),dpi=settings_plot['journal dpi']);
            # gridr = gridspec.GridSpec(nrows=7, ncols=1, figure=fig);
            # fig.subplots_adjust(hspace=0.75)
            gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
            gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[0], hspace = 5.0)
            gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
            # gridr.update(hspace=0.05); # set the spacing between axes.
            fig.add_subplot(gridr1[0]); #RTI plots are 2 tall
            # gridr.update(hspace=0.80); # set the spacing between axes.
            fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
            ax = fig.axes; #get a list of the axes
            ax[1].set_aspect('auto');
        #END IF
    #END IF
    
    divider = make_axes_locatable(ax[0]); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    
    #[(data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , keo_plotLatLong_chunks] ,
    pltHelprX, pltHelprY = np.meshgrid( (np.append(data_timeUnique,data_timeUnique[-1]+np.median(np.diff(data_timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                keo_plotLatLong_chunks);
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax[0].pcolormesh(pltHelprX, pltHelprY,  keo.T ,vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),cmap=keo_colormap, linewidth=0, rasterized=plot_rasterize); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        # if( FLG_fancyPlot == 0 ):
        #     cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        #     cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
        # else:
        #     cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),11)); #create useful tick marks
        #     cax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
        # #END IF
        cax_ticks = np.arange(np.min(plotLimValu),np.max(plotLimValu)+0.1,0.1); #busted out to deal with 0 being -0
        cax_ticks[np.isclose(cax_ticks, 0)] = 0; #enforce true 0
        cax.yaxis.set_ticks(cax_ticks); #create useful tick marks
        cax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
        cbar.set_label(keo_plotLabel); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(FONT_axisTickFM); #yee
        #END FOR tick
        #cbar.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #they changed how the code works, this doesn't work anymore
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu)); #now it's this
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    else:
        im = ax[0].pcolormesh(pltHelprX, pltHelprY,  keo.T ,cmap=keo_colormap, linewidth=0, rasterized=plot_rasterize); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(keo_plotLabel); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(FONT_axisTickFM); #yee
        #END FOR tick
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #END IF
    
    if( np.any(highlighter != False) ):        
        pltr_time = np.copy(data_timeUnique); #copy time
        k = np.logical_not(highlighter); #remove times not highlighted
        pltr_time = np.delete(pltr_time,k,axis=0); #delete em

        k = np.isclose(np.diff(pltr_time), np.median(np.diff(data_timeUnique))); #find contiguous time bits
        pltr_time = (pltr_time - dateRange_dayNum_zeroHr[1]*86400)/3600; #for plottin, adjust to hrs
        kj = np.where(~k)[0]+1; #add 1 for index b/c diff
        if(kj[0] != 0):
            kj = np.insert(kj, 0, 0); #insert 0
        #END IF
        if(kj[-1] != (pltr_time.size-1) ):
            kj = np.append(kj, pltr_time.size-1); #insert end size
        #END IF
        for j in range(0,kj.size-1):
            ax[0].axvspan(pltr_time[kj[j]], pltr_time[kj[j+1]-1], ymin=0, ymax=1, alpha=0.45, edgecolor='none', facecolor='xkcd:brick red');
        #END FOR j
    #END IF
    
    if( FLG_fancyPlot == 0 ): #only title non-fancy plot
        #    string_Title = 'TEC Averaged on Angle of '+textNice(np.round(avg_anyAngle,2))+' deg and Width of '+ \
        #        textNice(np.round(avg_anyAngle_Width,2))+' '+latlong_unitName, Avg Step # = '+textNice(avg_anyAngle_N)+ \
        #        ' '+latlong_unitName, Line Shows '+keo_plotLatLong_name+' of Millstone Hill Zenith Beam'; #create mecha title
        string_Title = keo_name+' Averaged on Angle of '+textNice(np.round(keo_angle,2))+' ° and Width of '+ \
            textNice(np.round(keo_width,2))+' arc°'; #create mecha title - for debug so hardcode unit
        ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    #END IF
    if( settings_dataSpecific['use local time'] == False ):
        ax[0].set_xlabel('Time in UT '+settings_plot['unit L']+'hr'+settings_plot['unit R']+' | 0 Hr on '+dateRange_zeroHr_monthName+' '+textNice(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' (Day '+textNice(dateRange_dayNum_zeroHr[1])+'), '+textNice(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    else:
        if( 'use local time override' in settings_dataSpecific ):
            if( settings_dataSpecific['use local time override'][0] == True ):
                avgPt_coords_hur = [[settings_dataSpecific['use local time override'][1][0], settings_dataSpecific['use local time override'][1][1]]]; #override!
                plotLatRange_hur = [avgPt_coords_hur[0][0]-8,avgPt_coords_hur[0][0]+8]; #make something up
                plotLongRange_hur = [avgPt_coords_hur[0][1]-8,avgPt_coords_hur[0][1]+8]; #make something up
            else:
                avgPt_coords_hur = avgPt_coords; #alias
                plotLatRange_hur = plotLatRange; #alias
                plotLongRange_hur = plotLongRange; #alias
            #END IF
        else:
            avgPt_coords_hur = avgPt_coords; #alias
            plotLatRange_hur = plotLatRange; #alias
            plotLongRange_hur = plotLongRange; #alias
        #END IF
        localTime_timeZone, localTime_UTCoffset = sunAlsoRises_timeZone(avgPt_coords_hur, plotLatRange_hur, plotLongRange_hur, dateRange_zeroHr);
        # ax[0].set_xlabel('Local Time [hr] | '+localTime_timeZone+' = UT '+localTime_UTCoffset+' hr | 0 hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' (Day '+str(dateRange_dayNum_zeroHr[1])+'), '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        ax[0].set_xlabel('Local Time '+settings_plot['unit L']+'hr'+settings_plot['unit R']+' | 0 hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label, salih specific request
    #END IF
    if( coordType == 'geo' ):
        ax[0].set_ylabel(keo_plotLatLong_name+' '+latlong_unitName_bracketed+keo_plotLatLong_dirAdder,fontproperties=FONT_axisLabelFM); #set the y axis label
    elif( coordType == 'mag' ):
        ax[0].set_ylabel(keo_plotLatLong_name+' (Geomag) '+latlong_unitName_bracketed+keo_plotLatLong_dirAdder,fontproperties=FONT_axisLabelFM); #set the y axis label
    #END IF
    
    if( 'day nite shading' in settings_dataSpecific ):
        if( settings_dataSpecific['day nite shading'] < 0 ):
            if( 'day nite shading lettering' in settings_dataSpecific ):
                if( settings_dataSpecific['day nite shading lettering'] == True ):
                    letteringB = ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[0].transAxes); #print the text saying the day or nite
                #END IF
            else:
                letteringB = ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[0].transAxes); #print the text saying the day or nite
            #END IF
        #END IF
    #END IF
    
    # #makes plotting look better
    # xAxisLims = np.array(((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600)); #get those xAxis limits, ensures everything is aligned
    # if( (np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 5*10**-2 ): #fix x axis stuff (time)
    #     xAxisLims[0] = np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    # #END IF
    # if( (np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 5*10**-2 ): #fix x axis stuff (time)
    #     xAxisLims[1] = np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    # #END IF
    # # ax.set_xlim(xAxisLims); #set the xlims now
    
    # yAxisLims = np.array(ax.get_ylim()); #get the current y axis limits
    # if( (np.round(np.min(keo_plotLatLong_chunks)) - np.min(keo_plotLatLong_chunks)) < 5*10**-2 ): #fix y axis stuff (lat or longitude)
    #     yAxisLims[0] = np.round(np.min(keo_plotLatLong_chunks)); #force the rounded value
    # #END IF
    # if( (np.round(np.max(keo_plotLatLong_chunks)) - np.max(keo_plotLatLong_chunks)) < 5*10**-2 ): #fix y axis stuff (lat or longitude)
    #     yAxisLims[1] = np.round(np.max(keo_plotLatLong_chunks)); #force the rounded value
    # #END IF
    # ax.set_ylim(yAxisLims); #set the ylims now
    
    # if( FLG_fancyPlot == 0 ): #different aspect ratios require different spacing assumptions
    #     keo_plotLatLong_chunks_autoTick = (np.ceil(np.max(keo_plotLatLong_chunks)) - np.floor(np.min(keo_plotLatLong_chunks)))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    #     if( keo_plotLatLong_chunks_autoTick > 25 ):
    #         keo_plotLatLong_chunks_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 10 ):
    #         keo_plotLatLong_chunks_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 5 ):
    #         keo_plotLatLong_chunks_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 2 ):
    #         keo_plotLatLong_chunks_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 1 ):
    #         keo_plotLatLong_chunks_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    #         keo_plotLatLong_chunks_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    #     else:
    #         if(keo_plotLatLong_name == 'Latitude'): #if Y axis is latitude, use latitude
    #             keo_plotLatLong_chunks_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/13; #just goes for it if it's a super tiny range
    #         elif(keo_plotLatLong_name == 'Longitude'): #if Y axis is longitude, use longitude
    #             keo_plotLatLong_chunks_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/13; #just goes for it if it's a super tiny range
    #         #END IF
    #     #END IF
    # else:
    #     keo_plotLatLong_chunks_autoTick = (np.ceil(np.max(keo_plotLatLong_chunks)) - np.floor(np.min(keo_plotLatLong_chunks)))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    #     if( keo_plotLatLong_chunks_autoTick > 20 ):
    #         keo_plotLatLong_chunks_autoTick = 30; #sets the tick setting to 30 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 10 ):
    #         keo_plotLatLong_chunks_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 5 ):
    #         keo_plotLatLong_chunks_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 2 ):
    #         keo_plotLatLong_chunks_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick > 1 ):
    #         keo_plotLatLong_chunks_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    #     elif( keo_plotLatLong_chunks_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    #         keo_plotLatLong_chunks_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    #     else:
    #         if(keo_plotLatLong_name == 'Latitude'): #if Y axis is latitude, use latitude
    #             keo_plotLatLong_chunks_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/17; #just goes for it if it's a super tiny range
    #         elif(keo_plotLatLong_name == 'Longitude'): #if Y axis is longitude, use longitude
    #             keo_plotLatLong_chunks_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/17; #just goes for it if it's a super tiny range
    #         #END IF
    #     #END IF
    # #END IF
    # yAxisTicks = np.round(np.arange( np.floor(np.min(keo_plotLatLong_chunks)),np.ceil(np.max(keo_plotLatLong_chunks)),keo_plotLatLong_chunks_autoTick ),2); #creates y ticks automagically
    # ax.set_yticks(yAxisTicks); #set x axis ticks
    # ax.set_ylim(yAxisLims); #set the ylims now [gotta do it after ticks in case they're a little too far]
    
    if( FLG_fancyPlot == 0 ): #different aspect ratios require different spacing assumptions
        GRITI_plotHelper_axisizerLatLong(keo_plotLatLong_chunks,ax=ax[0],axDir='y',tickNumGoal=17);
    else:
        GRITI_plotHelper_axisizerLatLong(keo_plotLatLong_chunks,ax=ax[0],axDir='y',tickNumGoal=13);
    #END IF
    
    # time_autoTick = (np.ceil((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.floor((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    # if( time_autoTick > 24 ):
    #     time_autoTick = 24*(time_autoTick//24); #sets the tick setting to 15 arcdegrees per tick
    # elif( time_autoTick > 12 ):
    #     time_autoTick = 12; #sets the tick setting to 15 arcdegrees per tick
    # elif( time_autoTick > 6 ):
    #     time_autoTick = 6; #sets the tick setting to 10 arcdegrees per tick
    # elif( time_autoTick > 4 ):
    #     time_autoTick = 4; #sets the tick setting to 5 arcdegrees per tick
    # elif( time_autoTick > 2 ):
    #     time_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    # elif( time_autoTick >= 1 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    #     time_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    # else:
    #     time_autoTick = ((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600 - (np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600)/13; #just goes for it if it's a super tiny range
    # #END IF
    # xAxisTicks = np.arange( (np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) , \
    #         (np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(data_timeRef)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) + time_autoTick , \
    #         time_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    # ax.set_xticks(xAxisTicks); #set x axis ticks
    # ax.set_xlim(xAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    if( np.any(timeCutout == None) ):
        timezAxisTicks, timezAxisLims = GRITI_plotHelper_axisizerTime((data_timeRef-dateRange_dayNum_zeroHr[1]*86400)/3600,ax=ax[0]); #automagic time ticks here
    else:
        timezAxisTicks, timezAxisLims = GRITI_plotHelper_axisizerTime((data_timeRef-dateRange_dayNum_zeroHr[1]*86400)/3600,ax=ax[0],FLG_manualLims=timeCutout); #automagic time ticks here
        ax[0].format_coord = plot_cursor_relabeler(ax[0]); #hopefully make it work
    #END IF
           
    
    #Now drawing line of interest
    if( (settings_map['coord type'] == 'geo') & ('keo coord type' in settings_dataSpecific) ): #pull a switcheroo if needed
        if( settings_dataSpecific['keo coord type'] == 'mag' ):
            import copy
            import datetime
            from Code.subfun_convertToMag import convert_to_mag
            
            if( 'pierceAlt' in settings_dataSpecific ):
                keo_alt = settings_dataSpecific['pierceAlt'];
            else:
                keo_alt = 120.; #default, great for auroral zone stuff (like field aligned currents)
            #END IF
            timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
            kk = np.where(np.int64(data_timeUnique[timeIndex]/86400) == dates['date range full dayNum'][:,1])[0].item(); #get where the year is gonna be
            time4mag_hr = np.int32(np.mod(data_timeUnique[timeIndex],86400)//3600); #get hours
            time4mag_min = np.int32(np.mod(data_timeUnique[timeIndex],86400)//60-time4mag_hr*60); #get the minutes
            time4mag_sec = np.int32(np.mod(data_timeUnique[timeIndex],86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
            time4mag = datetime.datetime(dates['date range full'][kk,0],dates['date range full'][kk,1],dates['date range full'][kk,2], \
                                         hour = time4mag_hr, minute = time4mag_min, second = time4mag_sec); #date time object for aacgmv2    
            
            avgPt_coords = copy.deepcopy(avgPt_coords); #copy big time
            [avgPt_coords[0,0], avgPt_coords[0,1]] = convert_to_mag(avgPt_coords[0,0], avgPt_coords[0,1], keo_alt, time4mag); #convert
        #END IF
    #END IF
    if( keo_plotLatLong_name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= avgPt_coords[0,1]) & (np.max(plotLongRange) >= avgPt_coords[0,1]) ): #only plot if it's in the long range specified
            ax[0].plot( np.linspace(np.min((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
                    np.ones(2)*avgPt_coords[0,1] , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmoller); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= avgPt_coords[0,0]) & (np.max(plotLatRange) >= avgPt_coords[0,0]) ): #only plot if it's in the lat range specified
            ax[0].plot( np.linspace(np.min((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((data_timeRef - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
                    np.ones(2)*avgPt_coords[0,0] , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmoller); #plots a point with a black line
        #END IF
    #END IF
    
    
    #----- Plot Shading to represent 'night' -----
    if( 'day nite shading' in settings_dataSpecific ):
        if( settings_dataSpecific['day nite shading'] == 1 ):
            if( (np.min(plotLatRange) <= avgPt_coords[0,0]) & (np.max(plotLatRange) >= avgPt_coords[0,0]) ): #only plot if it's in the lat range specified
                if( (avgPt_coords[0,0] > 45) & (np.min(plotLatRange) < 45) ):
                    latToUse = 45; #cap at 50 deg lat to keep it from getting too zesty at the poles
                elif( (avgPt_coords[0,0] < -45) & (np.max(plotLatRange) < -45) ):
                    latToUse = -45; #cap at 50 deg lat to keep it from getting too zesty at the poles
                else:
                    latToUse = avgPt_coords[0,0];
                #END IF
            else:
                latToUse = np.mean(plotLatRange); #just take avg of plt area I guess
            #END IF
            if( (np.min(plotLongRange) <= avgPt_coords[0,1]) & (np.max(plotLongRange) >= avgPt_coords[0,1]) ): #only plot if it's in the long range specified
                longToUse = avgPt_coords[0,1];
            else:
                longToUse = np.mean(plotLongRange); #just take avg of plt area I guess
            #END IF
            
            (niteTimes_sunRise, niteTimes_sunSet, niteTimes_dateRange_fullPad) = sunAlsoRises(dates['date range full'],latToUse,longToUse); #call sunrise/set function

            niteTime_dateRange_dayNum_fullPad = subfun_date_to_dayNum(niteTimes_dateRange_fullPad); #convert to dayNum
            niteTimes_sunRise = (niteTimes_sunRise + niteTime_dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hrs, center around zero hr and convert ot hrs
            niteTimes_sunSet = (niteTimes_sunSet + niteTime_dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hrs, center around zero hr and convert ot hrs
            niteTimes_sunRise = niteTimes_sunRise[1:]; #remove 1st
            niteTimes_sunSet = niteTimes_sunSet[:-1]; #remove last
            #FIFTH STEP: PLOT THIS STUFF
            for j in range(0,niteTimes_sunSet.size):
                if(keo_plotLatLong_name == 'Latitude'): #if Y axis is latitude, use latitude
                    recta = Rectangle((niteTimes_sunSet[j], np.min(plotLatRange)), niteTimes_sunRise[j]-niteTimes_sunSet[j], np.diff(plotLatRange).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.20); #make that patch
                    ax[0].add_patch(recta); #add on that patch
                else: #otherwise longitude
                    recta = Rectangle((niteTimes_sunSet[j], np.min(plotLongRange)), niteTimes_sunRise[j]-niteTimes_sunSet[j], np.diff(plotLongRange).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.20); #make that patch
                    ax[0].add_patch(recta); #add on that patch
                #END IF
            #END FOR j
        elif( settings_dataSpecific['day nite shading'] == 2 ):
            for j in range(0,dates['date range zero hr hours'].size):
                doubleKeo_niteTimes_per = dates['date range zero hr hours'][j] + np.asarray(settings_dataSpecific['day nite shading times']); #get the nite time ranges
                if(keo_plotLatLong_name == 'Latitude'): #if Y axis is latitude, use latitude
                    recta = Rectangle((doubleKeo_niteTimes_per[0], np.min(plotLatRange)), np.diff(doubleKeo_niteTimes_per).item(), np.diff(plotLatRange).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.20); #make that patch
                    ax[0].add_patch(recta); #add on that patch
                else: #otherwise longitude
                    recta = Rectangle((doubleKeo_niteTimes_per[0], np.min(plotLongRange)), np.diff(doubleKeo_niteTimes_per).item(), np.diff(plotLongRange).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.20); #make that patch
                    ax[0].add_patch(recta); #add on that patch
                #END IF
            #END FOR j
        elif( settings_dataSpecific['day nite shading'] == -1 ):
            if( (np.min(plotLatRange) <= settings_map['site coords'][0,0]) & (np.max(plotLatRange) >= settings_map['site coords'][0,0]) ): #only plot if it's in the lat range specified
                if( (settings_map['site coords'][0,0] > 45) & (np.min(plotLatRange) < 45) ):
                    latToUse = 45; #cap at 50 deg lat to keep it from getting too zesty at the poles
                elif( (settings_map['site coords'][0,0] < -45) & (np.max(plotLatRange) < -45) ):
                    latToUse = -45; #cap at 50 deg lat to keep it from getting too zesty at the poles
                else:
                    latToUse = settings_map['site coords'][0,0];
                #END IF
            else:
                latToUse = np.mean(plotLatRange); #just take avg of plt area I guess
            #END IF
            if( (np.min(plotLongRange) <= settings_map['site coords'][0,1]) & (np.max(plotLongRange) >= settings_map['site coords'][0,1]) ): #only plot if it's in the long range specified
                longToUse = settings_map['site coords'][0,1];
            else:
                longToUse = np.mean(plotLongRange); #just take avg of plt area I guess
            #END IF
            
            #----- DAYNITE STUFF IF NEEDED -----
            #FIRST: GET SUNRISE/SUNSET TIMES
            (dayNite_sunrise, dayNite_sunset, daynites_dateRange_fullPad) = sunAlsoRises(dates['date range full'],latToUse,longToUse); #call sunrise/set function
        
            daynite_dateRange_dayNum_fullPad = subfun_date_to_dayNum(daynites_dateRange_fullPad); #convert to dayNum
            dayNite_sunrise = (dayNite_sunrise + daynite_dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hrs, center around zero hr and convert ot hrs
            dayNite_sunset = (dayNite_sunset + daynite_dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hrs, center around zero hr and convert ot hrs
            dayNite_sunrise = dayNite_sunrise[1:]; #remove 1st
            dayNite_sunset = dayNite_sunset[:-1]; #remove last
            
            #SECOND: GET TIME ZONE FOR LOCAL TIME
            tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
            dayNite_timeZoneID = tf.certain_timezone_at(lat=settings_map['site coords'][0][0], lng=settings_map['site coords'][0][1]); #use it to find the time zone
            if dayNite_timeZoneID is None:
                #use geonames site as a backup
                url = 'http://api.geonames.org/timezone?lat='+str(settings_map['site coords'][0][0])+'&lng='+str(settings_map['site coords'][0][1])+'&username='+settings_config['login GeoNames TimeZone']['user']; #create link for lat/long
                webpage = urlopen(url).read(); #get the raw HTML and read it
                try:
                    charset = webpage.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                    if( charset is None ):
                        charset = 'utf-8'; #assume utf-8
                    #END IF
                except:
                    charset = 'utf-8'; #assume utf-8
                #END TRY
                webpage = webpage.decode(charset); #"decode" the HTML content so it's legible
                # index_start = strstr(webpage,'<dstOffset>')[0]; #get where dstOffset is
                # index_end = strstr(webpage,'</dstOffset>')[0]; #get where dstOffset is
                # timeZone_offset_str = webpage[index_start+11:index_end]; #get dst offset
                # timeZone_offset = np.float64(timeZone_offset_str); #and convert to number
                # index_start = strstr(webpage,'<gmtOffset>')[0]; #get where UT offset is
                # index_end = strstr(webpage,'</gmtOffset>')[0]; #get where UT offset is
                # dayNite_UToffset_str = webpage[index_start+11:index_end]; #get UT offset
                # dayNite_UToffset = np.float64(dayNite_UToffset_str); #and convert to number
                index_start = strstr(webpage,'<timezoneId>')[0]; #get where time zone ID is
                index_end = strstr(webpage,'</timezoneId>')[0]; #get where time zone ID is
                dayNite_timeZoneID = webpage[index_start+12:index_end]; #get time zone ID
            #END IF
            
            #THIRD: TIME ZONE STUFF
            timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
            timeZoneObj_UTC = pytz.timezone('UTC'); #make a timezone object
            
            timeZone_datetime = timeZoneObj_UTC.localize(datetime.strptime(str(dates['date range zero hr'][0])+'-'+str(dates['date range zero hr'][1])+'-'+str(dates['date range zero hr'][2])+'T00:00:00.000', \
                '%Y-%m-%dT%H:%M:%S.%f')).astimezone(timeZoneObj); #create datetime object, set it to the UTC time zone (which it is, datetime just doesn't know), then convert it to local time zone
            timeZone_offset_str = timeZone_datetime.isoformat()[strstr(timeZone_datetime.isoformat(),'-')[-1]:]; #get the DST time offset
            if( strstr(timeZone_offset_str,':').size > 0 ):
                if( timeZone_offset_str[strstr(timeZone_offset_str,':')[0]+1:] == '00' ):
                    timeZone_offset_str = timeZone_offset_str[:strstr(timeZone_offset_str,':')[0]]; #remove the :00 if it's just that
                    if( timeZone_offset_str[1] == '0' ):
                        timeZone_offset_str = timeZone_offset_str.replace('0',''); #remove the 0 that was extraneous
                    #END IF
                    timeZone_offset = np.int64(timeZone_offset_str); #get the number version
                #END IF
                else:
                    timeZone_offset = np.float64(timeZone_offset_str[:strstr(timeZone_offset_str,':')[0]]) + \
                        np.int64(timeZone_offset_str[strstr(timeZone_offset_str,':')[0]+1:])/60; #convert to an hour decimal
                #END IF
            #END IF
            timeZone_name = timeZone_datetime.tzname(); #get the time zone name (like 'EST' or 'EDT' depending on standard or daylight savings time)
            timeZone_DSTnUTCOffset = timeZone_datetime.dst().total_seconds()/3600; #get the time offset
            if( np.mod(timeZone_DSTnUTCOffset,1) == 0 ):
                timeZone_DSTnUTCOffset = np.int64(timeZone_DSTnUTCOffset); #convert to integer
            #END IF
                
            #FOURTH STEP: PREP FOR PLOTTING BY ALIGNING TIMES, MAKING PLOT VARIABLES
            dayNite_sunrise += timeZone_offset; #make local time
            dayNite_sunset += timeZone_offset; #make local time
        
            xTime = np.sort( np.concatenate( (dayNite_sunrise,dayNite_sunset) ) ); #hr, xvar to plot against
            yDayNite = np.ones(xTime.shape); #prep if day or night, set all to day
            for i in range(0,dayNite_sunset.size):
                yDayNite[xTime == dayNite_sunset[i]] = 0; #set sunset times to sunset
            #END FOR i
            xTime = xTime.repeat(2); #interleave repeated values
            yDayNite = np.roll(yDayNite.repeat(2),1) #interleave repeated values and circular shift by 1
            yDayNite[0] = yDayNite[1]; #set that to match (for plotting niceness)
            yDayNite[-1] = yDayNite[-2]; #set that to match (for plotting niceness)
            
            # xTime[ xTime < np.min(timezAxisLims)+timeZone_offset ] = np.min(timezAxisLims)+timeZone_offset; #limit the xTime to the plotted times
            # xTime[ xTime > np.max(timezAxisLims)+timeZone_offset ] = np.max(timezAxisLims)+timeZone_offset;
            # xTime[0] = xTime[1]; #set that to match (for plotting niceness)
            # xTime[-1] = xTime[-2]; #set that to match (for plotting niceness)
            yDayNite = np.delete(yDayNite, (xTime > np.max(timezAxisLims)+timeZone_offset) | (xTime < np.min(timezAxisLims)+timeZone_offset) ); #remove out of bounds stuff
            xTime = np.delete(xTime, (xTime > np.max(timezAxisLims)+timeZone_offset) | (xTime < np.min(timezAxisLims)+timeZone_offset) ); #remove out of bounds stuff
            #fill in edges so they just end
            yDayNite = np.append(np.insert(yDayNite,0,yDayNite[0]),yDayNite[-1]);
            xTime = np.append(np.insert(xTime,0,np.min(timezAxisLims)+timeZone_offset),np.max(timezAxisLims)+timeZone_offset);
            
            
            #FIFTH STEP: ACTUALLY PLOTTING
            divider2 = make_axes_locatable(ax[1]); #prep to add an axis
            cax2 = divider2.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
            cax2.set_visible(False); #mkae it invisible so it matches the other plots in width
            ax[1].plot(xTime,yDayNite,color='xkcd:black',linewidth=settings_plot['line width']['thicc'], antialiased=True);
            if( timeZone_DSTnUTCOffset != 0 ):
                strang = timeZone_offset_str+' (Daylight Savings) '+timeZone_name+' Time Zone';
            else:
                 strang = timeZone_offset_str+' '+timeZone_name+' Time Zone';
            #END IF
            # if( yDayNite[-1] == yDayNite[-2] ):
            #     #hacky offset so day or night isn't printed off the plot (or skipped)
            #     offset = 3;
            # else:
            #     offset = 1;
            # #END IF
            xTime_diff_median = np.median(np.diff(xTime)[::2]); #average day/nite,improve alg for higher alts by including yDayNite but keyboard tiny so not now
            for i in range(0,xTime.size,2):
                if( xTime_diff_median*.5 < (xTime[i+1]-xTime[i]) ):
                    if( yDayNite[i] == 1 ):
                        ax[1].text( (xTime[i+1]-xTime[i])/2+xTime[i], \
                            0.25,'Day', color='k', ha='center', va='baseline', fontproperties=settings_plot['font title FM']); #print the text saying the day or nite
                    else:
                        ax[1].text( (xTime[i+1]-xTime[i])/2+xTime[i], \
                           0.25, 'Night', color='k', ha='center', va='baseline', fontproperties=settings_plot['font title FM']); #print the text saying the day or nite
                    #END IF
                #END IF
            #END FOR i
            ax[1].set_xticks(timezAxisTicks+timeZone_offset); #set x axis ticks
            if( settings_dataSpecific['use local time'] == False ):
                ax[1].set_xlabel('Local Time '+settings_plot['unit L']+'hr'+settings_plot['unit R']+' | '+strang,fontproperties=settings_plot['font axis label FM']);
            else:
                ax[1].set_xticklabels([]); #remove x labels for day nite plost
            #END IF
            # ax[1].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
            ax[1].set_yticklabels([]); #remove y labels for day nite plost
            ax[1].set_yticks([]); #remove y tick marks for day nite plot
            ax[1].set_xlim( (timezAxisLims[0]+timeZone_offset,timezAxisLims[1]+timeZone_offset) ); #set x axis limits
            ax[1].set_ylim( (0,1) ); #set y lims to 0 and 1
            ax[1].spines['left'].set_visible(False); #turn off box lines
            ax[1].spines['right'].set_visible(False); #turn off box lines
            ax[1].spines['top'].set_visible(False); #turn off box lines
            # ax[1].grid(b=True, which='major', axis='x', color='xkcd:light grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
            # ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[1].transAxes); #print the text saying the day or nite
        
            if( settings_dataSpecific['day nite shading'] < 0 ):
                if( 'day nite shading lettering' in settings_dataSpecific ):
                    if( settings_dataSpecific['day nite shading lettering'] == True ):
                        letteringB = ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[1].transAxes); #print the text saying the day or nite
                    #END IF
                else:
                    letteringB = ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[1].transAxes); #print the text saying the day or nite
                #END IF
            #END IF
        #END IF
    #END IF
    
    if( FLG_drawSun ):
        from Code.subfun_sunAlsoRises_location import sunAlsoRises_location
        # from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
        # import astropy.coordinates as coord
        # from astropy.time import Time, conf
        # #Now draw sun position
        # #yes I just abs wing making these solutions every time hoping to find the codegrail each time, I don't tho
        # #!NO YEAR SUPPORT YET!
        # #time to make some strings of time
        # timeString = [[] for i in range(0,data_timeUnique.size)]; #preallocate list to hold the strings
        # timeTime = [[] for i in range(0,data_timeUnique.size)]; #preallocate list to hold the strings
        # #---Build time things needed---
        # timeMonDay = subfun_dayNum_to_date(np.array([np.ones(data_timeUnique.size,dtype=np.int64)*dates['date range full dayNum'][0,0],np.int64(data_timeUnique/86400)]).T); #get the month and day
        # timeTempHrFloat = (data_timeUnique/86400 - np.int64(data_timeUnique/86400))*24; #hr, get the decimal bit
        # timeTempHr = np.int32(timeTempHrFloat); #hr, get the hours
        # timeTempMinFloat = (timeTempHrFloat - timeTempHr)*60; #min, get the decimal bit
        # timeTempMin = np.int32(timeTempMinFloat); #min, get the minutes
        # timeTempSec = np.int32((timeTempMinFloat - timeTempMin)*60); #sec, get the sec (we won't go sub-sec)
        # #---fix float32 time issues (we know we only operate on whole seconds)---
        # timeTempMin[timeTempSec == 59] += 1; #increment
        # timeTempSec[timeTempSec == 59] = 0; #reset to 0
        # timeTempHr[timeTempMin == 60] += 1; #increment
        # timeTempMin[timeTempMin == 60] = 0; #reset to 0
        # timeTempHr[timeTempHr == 24] = 0; #reset to 0
        # for i in range(0,data_timeUnique.size):
        #     timeString[i] = str(dates['date range full dayNum'][0,0])+'-'+str(timeMonDay[i,1]).zfill(2)+'-'+str(timeMonDay[i,2]).zfill(2)+'T'+ \
        #         str(timeTempHr[i]).zfill(2)+':'+str(timeTempMin[i]).zfill(2)+':'+str(timeTempSec[i]).zfill(2)+'.000'; #make the strings
        # #END FOR i
        # for i in range(0,data_timeUnique.size):
        #     #uses ISOT format 2000-01-01T00:00:00.000
        #     with conf.set_temp('use_fast_parser', 'force'):
        #         try:
        #             timeTime[i] = Time(timeString[i], format='isot', scale='utc'); #convert
        #         except ValueError as errorz:
        #             print(errorz); #that's a whoopsie (should never happen)
        #         #END TRY
        #     #END WITH
        # #END FOR i
        # sunPlot = np.empty(data_timeUnique.size,dtype=np.float32); #preallocate
        # sunPlot_other = np.empty(data_timeUnique.size,dtype=np.float32); #preallocate
        # for i in range(0,data_timeUnique.size):
        #     sun = coord.EarthLocation.from_geocentric(coord.get_sun(timeTime[i]).transform_to(coord.ITRS).x,coord.get_sun(timeTime[i]).transform_to(coord.ITRS).y,coord.get_sun(timeTime[i]).transform_to(coord.ITRS).z); #trainwreck of a call to get the sun in lat/long
        #     if( keo_plotLatLong_name == 'Longitude' ): #if true, longitude
        #         sunPlot[i] = sun.lon.value; #get that data
        #         sunPlot_other[i] = sun.lat.value; #get that data
        #     else:
        #         sunPlot[i] = sun.lat.value; #get that data
        #         sunPlot_other[i] = sun.lon.value; #get that data
        #     #END IF
        # #END FOR i
        sunSubSolar_loc = sunAlsoRises_location(dates['date range full padded dayNum'],timeIndexes=data_timeUnique); #calc sun locations
        if( keo_plotLatLong_name == 'Longitude' ): #if true, longitude
            sunPlot = sunSubSolar_loc['long'];
            sunPlot_other = sunSubSolar_loc['lat'];
        else:
            sunPlot = sunSubSolar_loc['lat'];
            sunPlot_other = sunSubSolar_loc['long'];
        #END IF
        if( (settings_dataSpecific['keo coord type'] == 'mag') ):
            import datetime
            from Code.subfun_convertToMag import convert_to_mag
            
            if( 'pierceAlt' in settings_dataSpecific ):
                keo_alt = settings_dataSpecific['pierceAlt'];
            else:
                keo_alt = 120.; #default, great for auroral zone stuff (like field aligned currents)
            #END IF
            time4mag = [None for i in range(data_timeUnique.size)]; #make same size as data
            for i in range(data_timeUnique.size):
                # timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
                kk = np.where(np.int64(data_timeUnique[i]/86400) == dates['date range full padded dayNum'][:,1])[0].item(); #get where the year is gonna be
                time4mag_hr = np.int32(np.mod(data_timeUnique[i],86400)//3600); #get hours
                time4mag_min = np.int32(np.mod(data_timeUnique[i],86400)//60-time4mag_hr*60); #get the minutes
                time4mag_sec = np.int32(np.mod(data_timeUnique[i],86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
                time4mag[i] = datetime.datetime(dates['date range full padded'][kk,0],dates['date range full padded'][kk,1],dates['date range full padded'][kk,2], \
                                             hour = time4mag_hr, minute = time4mag_min, second = time4mag_sec); #date time object for aacgmv2    
            #END FOR i
            keo_alt = np.ones(data_timeUnique.size)*keo_alt; #make same size as data
            if( keo_plotLatLong_name == 'Longitude' ): #if true, longitude
                [sunPlot_other, sunPlot] = convert_to_mag(sunPlot_other, sunPlot, keo_alt, time4mag); #convert
            else:
                [sunPlot, sunPlot_other] = convert_to_mag(sunPlot, sunPlot_other, keo_alt, time4mag); #convert to mag
            #END IF
        #END IF
        # ax.plot( (data_timeUnique - dateRange_dayNum_zeroHr[1])*24, sunPlot,c='xkcd:orange',linewidth=2.5); #plots a point with an orange line
        if( keo_plotLatLong_name == 'Longitude' ): #if true, longitude
            indexes2plot = np.hstack( (np.arange(0,data_timeUnique.size,35),(data_timeUnique.size-1,)) ); #plot select stars
            ax[0].plot( ((data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600)[indexes2plot], sunPlot[indexes2plot], c='xkcd:orange', marker='*', markersize=20, linestyle='None'); #plots a point with an orange line
        else:
            indexes2plot = np.hstack( (np.arange(0,data_timeUnique.size,35),(data_timeUnique.size-1,)) ); #plot select stars
            ax[0].plot( ((data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600)[indexes2plot], sunPlot[indexes2plot], c='xkcd:orange', marker='*', markersize=20, linestyle='None'); #plots a point with an orange line
        #END IF
        # ax[0].plot( np.linspace(np.min((data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((data_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),2,endpoint=True) , #X time hr
        #         np.tile(-72.68,2) , #Y latitude OR longitude arcdeg
        #         c='xkcd:black',linewidth=1,linestyle='--'); #plots a point with a black line
    #END IF
    
    if( settings_dataSpecific['use local time'] == True ):
        labelz = ax[0].get_xticklabels(); #get them labels
        #--- edit ticks ---
        if( localTime_UTCoffset != 0 ):
            for jk in range (0,len(labelz)):
                middleLad = timezAxisTicks[jk] + np.int64(localTime_UTCoffset);
                if( np.isclose(middleLad, np.int64(np.round(middleLad))) ):
                    middleLad = np.int64(np.round(middleLad)); #it can be simply an int64
                #END IF
                labelz[jk].set_text(textNice(middleLad));
            #END FOR jk
        #END IF
        #--- applky ticks---
        ax[0].set_xticklabels(labelz); #apply new labels
    #END IF
    
    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.050, right = 0.945, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    if( FLG_fancyPlot != 0 ):
        if( np.any(timeCutout != None) ):
            figFileName = path_join(settings_paths['fancyPlots'],settings_dataSpecific['keo data type']+'_keo_cutout'+settings_plot['save file type']);
        elif( FLG_drawSun == True ):
            figFileName = path_join(settings_paths['fancyPlots'],settings_dataSpecific['keo data type']+'_keo_wSun'+settings_plot['save file type']); #save the figure
        else:
            figFileName = path_join(settings_paths['fancyPlots'],settings_dataSpecific['keo data type']+'_keo'+settings_plot['save file type']); #save the figure
        #END IF
        # if( settings_plot['save file type'].lower() != '.pdf' ):
        fig.savefig(figFileName); #save the figure
        # else:
        #     from matplotlib.backends.backend_pdf import PdfPages
        #     pp = PdfPages(figFileName);
        #     pp.savefig(fig);
        #     pp.close();
        #END IF
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF

def plot_cursor_relabeler(ax): #inspired by https://stackoverflow.com/a/21585524/2403531
    #ax_under needs to be a list, even if just [ax0]
    #ax_names should be in the order of [top, most under, less under, least under] (e.g. 1 longer than ax_under)
    def format_coord(x, y):
        # x, y are data coordinates, they appear mystically
        # returnString = '';
        # display_coord = ax.transData.transform((x,y)); # convert to display coords
        # ax_coord = ax.transData.inverted().transform(display_coord); #fire up the inverter & convert back to data coords with respect to ax
        # coords = []; #prep list
        # for j in range(0,len(ax_under)):
        #     ax_coord = ax_under[j].transData.inverted().transform(display_coord); #fire up the inverter & convert back to data coords with respect to ax
        #     print(str(ax_coord))
        #     print(str(len(ax_under)))
        #     print(str(ax_names))
        #     # coords.append(ax_coord); #create list of coords
        #     returnString = returnString + ax_names[j]+' ({:.3f}, {:.3f}) | '.format(ax_coord[0], ax_coord[1]); #tack it on
        # #END FOR j
        # returnString = ' | '+ax_names[-1]+': ({:.3f}, {:.3f})'.format(x, y); #end of the string        
        # return (returnString)
        return ('time={:.3f}, lat/long={:.3f}'.format(x, y))
    #END DEF
    return format_coord
#END DEF


