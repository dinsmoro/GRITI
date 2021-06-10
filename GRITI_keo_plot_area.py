"""
GOAL: Plot area that average alg with average over
RD on 6/11/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tick
# import cartopy as cartopy #cartopy replaces basemap b/c it's updated
# import os
from subfun_figFitter import figFitter

def GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
        FLG_fancyPlot = 0):
    #----Start plotting-----
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo_area IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    
    #--- Unpack ---
    PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    # PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    # PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
    journal_dpi = settings_plot['journal dpi'];
    plotLongRange = settings_map['long range'];
    plotLatRange = settings_map['lat range'];
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum'];
    keo_angle = settings_dataSpecific['keo angle'];
    keo_width = settings_dataSpecific['keo width'];
    keo_N = settings_dataSpecific['keo N'];
    plotLimValu = settings_dataSpecific['keo plot lim'];
    avgPt_coords = settings_map['site coords'];
    
    #THIS IS TO VIEW DATA AVG RANGE BEING TAKEN
    if( FLG_fancyPlot == 0 ):
        fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        
        map_lat_autoTick = settings_map['lat autotick']; #get the auto ticks
        map_long_autoTick = settings_map['lat autotick']; #get the auto ticks
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
        
        map_lat_autoTick = settings_map['lat autotick fancy']; #get the auto ticks
        map_long_autoTick = settings_map['lat autotick fancy']; #get the auto ticks
    #END IF
    
    ax = plt.axes(projection=settings_map['projection']); #redefine the axis to be a geographical axis [AFTER cax else there's ERRORS oof]
    ax.set_aspect('auto');

    cax = ax.inset_axes((1.02, 0, 0.02, 1)); #make a color bar axis
    # divider = make_axes_locatable(ax); #prep to add an axis
    # cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #---ADD GRID LINES, SET PLOTTING AREA---
    # gl = ax.gridlines(linewidth=PLOT_lineWidthSmoller, color='xkcd:black', alpha=0.5, linestyle='--', draw_labels=True); #draw some well-described gridlines
    # gl.xlabels_top = False; #turn off all, let ticks be handled by set_xticks
    # gl.xlabels_bottom = False; #turn off all, let ticks be handled by set_xticks
    # gl.ylabels_right = False; #turn off all, let ticks be handled by set_yticks
    # gl.ylabels_left = False; #turn off all, let ticks be handled by set_yticks
    # gl.xlocator = tick.FixedLocator(np.arange(np.min(plotLongRange),np.max(plotLongRange)+map_long_autoTick,map_long_autoTick)); #this works ok, but be consistent use set_xticks
    # gl.ylocator = tick.FixedLocator(np.arange(np.min(plotLatRange),np.max(plotLatRange)+map_lat_autoTick,map_lat_autoTick)); #this doesn't plot -90 and 90 labels, but is req to get the gridlines right
    ax.set_xticks(np.arange(np.min(plotLongRange),np.max(plotLongRange)+map_long_autoTick,map_long_autoTick),crs=settings_map['projection']); #gotta plot ticks with this to get -90 and 90
    ax.set_yticks(np.arange(np.min(plotLatRange),np.max(plotLatRange)+map_lat_autoTick,map_lat_autoTick),crs=settings_map['projection']); #gotta plot ticks with this to get -90 and 90
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    ax.set_extent(plotLongRange + plotLatRange); #set the plot extent, set at end - x and y ticks can extend plot area so this'll reign it in
    
    #---DRAW SOME COASTLINES, MAYBE COLOR IN SOME STUFF---
    fig.canvas.draw(); #key for all instances
    if( plt.isinteractive() == True ): #only needs fig.canvas.draw() for interactive, flush and show and pause aren't req
        fig.canvas.flush_events(); #only needed for interactive plotting
        plt.show(); #req to make plot show up
        plt.pause(0.01); #this is oddly highly required for plotting in an IDE interactively (seems 0.01 is lowest I can go)
    #END IF
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get info on the size of the plot area to know what geographic scale to use
    mapper_resolution = np.max( [np.abs(plotLatRange[0]-plotLatRange[1])/bbox.height , np.abs(plotLongRange[0]-plotLongRange[1])/bbox.width] ); #degc, max extent covered in the plot
    if( mapper_resolution > 20 ): #arbitrary numbers
        mapper_resolution = '110m'; #the resolution to use for plotting geographical features
    elif( mapper_resolution > 0.5 ): #arbitrary numbers
        mapper_resolution = '50m'; #the resolution to use for plotting geographical features
    else:
        #otherwise if the deg/in for the plot is super small use the highest detail possible
        mapper_resolution = '10m'; #the resolution to use for plotting geographical features
    #END IF
    # if( settings_map['world color'] == True ):
    #     ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', mapper_resolution, edgecolor='face', facecolor=settings_map['land color'], alpha=0.75,zorder=75)); #idk what these calls really mean
    #     ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'ocean', mapper_resolution, edgecolor='face', facecolor=settings_map['water color'], alpha=0.75,zorder=75)); #idk what these calls really mean
    # #END IF
    ax.coastlines(resolution=mapper_resolution, color='xkcd:black',zorder=75); #draw the coastlines
    #reinforce axis limits after this drawing stuff
    ax.set_xlim(plotLongRange); #set x limits
    ax.set_ylim(plotLatRange); #set y limits
    
    #---ACTUALLY PLOT REAL STUFF HERE---
    # (_, counts) = np.unique(data_time, return_counts=True); #get the counts of the data
    # timeIndex = np.where( counts >= data_time.size/data_timeUnique.size )[0][0]; #get an index where there's a lot of data
    timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
    timeHr = np.int16((data_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600); #calc the hour
    timeMin = np.int16(np.abs((data_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)/60); #calc the minute
    timeSec = np.int16(np.round(np.abs((data_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)/60 - timeMin)); #calc the second
    if( timeSec == 60 ):
        timeSec = 0; #roll bak to  0
        timeMin += 1; #increment by 1
    #END IF
    if( timeMin == 60 ):
        timeMin = 0; #roll back to 0
        timeHr += 1; #increment by 1
    #END IF
    
    k = np.where( data_time == data_timeUnique[timeIndex])[0]; #gets during a time period
    
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['keo scatter size'],c=data_data[k],cmap=settings_dataSpecific['keo colormap'], vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),zorder=100);
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(settings_dataSpecific['keo name']+settings_dataSpecific['keo units']); #tabel the colorbar
        cbar.ax.tick_params(labelsize=settings_plot['font axis tick']);
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu));
        # cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),11)); #create useful tick marks
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    else:
        im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['keo scatter size'],c=data_data[k],cmap=settings_dataSpecific['keo colormap'],zorder=100);
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(settings_dataSpecific['keo name']+settings_dataSpecific['keo units']); #tabel the colorbar
        cbar.ax.tick_params(labelsize=settings_plot['font axis tick']);
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #END IF
    if( FLG_fancyPlot == 0 ):
        string_Title = settings_dataSpecific['keo name']+' Keo w/ Angle='+str(np.round(keo_angle,2))+' deg, Width='+ \
            str(np.round(keo_width,2))+' arcdeg, Step #='+str(keo_N)+ \
            ', Step Width='+ \
            str(np.round(np.sqrt( (temp_Longs_up[0,0] - temp_Longs_up[0,1])**2 + (temp_Lats_up[0,0] - temp_Lats_up[0,1])**2 ),2))+ \
            ' arcdeg | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' at '+str(timeHr)+':'+str(timeMin)+':'+str(timeSec); #create mecha title
        ax.set_title(string_Title,fontproperties=settings_plot['font title FM']); #set the title
    #END IF
    
    #draw a box where the data will be gathered
    temp_mapCoords = ( np.hstack( [np.linspace(keo_range[0,1],keo_range[1,1],200) , \
        np.linspace(keo_range[1,1],keo_range[3,1],200) , \
        np.linspace(keo_range[3,1],keo_range[2,1],200) , \
        np.linspace(keo_range[2,1],keo_range[0,1],200)] ) , \
        np.hstack( [np.linspace(keo_range[0,0],keo_range[1,0],200) , \
        np.linspace(keo_range[1,0],keo_range[3,0],200) , \
        np.linspace(keo_range[3,0],keo_range[2,0],200) , \
        np.linspace(keo_range[2,0],keo_range[0,0],200)] ) ); #convert to the geographic map coords
    ax.plot( temp_mapCoords[0],  #X longitude arcdeg
        temp_mapCoords[1],  #Y latitude arcdeg
        c='xkcd:fuchsia',linewidth=PLOT_lineWidthThicc, zorder=90);
    
    if( (np.arange(10,keo_N,10).size > 10) & (np.arange(10,keo_N,10).size <= 30) ):
        #a good zone - not too few, not too many
        temp_arangeStart = 10;
        temp_arangeSpacing = 10;    
    else:
        #otherwise need different start and spacing
        if( keo_N > 10 ):
            temp_arangeStart = np.int64(np.round(keo_N/10)); #scale it
            temp_arangeSpacing = np.int64(np.round(keo_N/10)); #scale it
        else:
            temp_arangeStart = 1; #min it at 1
            temp_arangeSpacing = 1; #min it at 1
        #END IF
    #END IF
            
    for i in np.arange(temp_arangeStart,keo_N,temp_arangeSpacing):
        temp_mapCoords = ( np.linspace( temp_Long_List[i,0],temp_Long_List[i,3],200 ) , \
            np.linspace( temp_Lat_List[i,0],temp_Lat_List[i,3],200 ) ); #convert to the geographic map coords
        
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c='xkcd:fuchsia',linewidth=PLOT_lineWidthSmol, zorder=90);
    #END FOR i
    
    #plot a * where the ISR is
    if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
        temp_mapCoords = (avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=settings_map['site marker type'], color=settings_map['site marker color'], markersize=settings_map['site marker size'], zorder=150);
    #END IF
    
    figFitter(fig); #fit the fig fast
    if( FLG_fancyPlot == 0 ):
        plt.show(); #req to make plot show up
    else:
        fig.savefig(settings_paths['fancyPlots']+'\\'+settings_dataSpecific['keo data type']+'_keo_area.png'); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF