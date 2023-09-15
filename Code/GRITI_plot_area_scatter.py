"""
GOAL: Plot area that average alg with average over
RD on 6/11/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tickr
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
import cartopy.feature as cfeature #cartopy replaces basemap b/c it's updated
import os
from Code.subfun_figFitter import figFitter

def GRITI_plot_area_scatter(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, \
        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
        dates = None, boxArea = False, FLG_useAllTimes=False, FLG_fancyPlot = 0):
    #----Start plotting-----
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: '+settings_dataSpecific['name'].replace(' ','-')+'_area_scatter IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    
    #--- Unpack ---
    PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    # PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    # PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
    journal_dpi = settings_plot['journal dpi'];
    plotLongRange = settings_map['long range'];
    plotLatRange = settings_map['lat range'];
    plotLimValu = settings_dataSpecific['plot lim'];
    avgPt_coords = settings_map['site coords'];
    
    #allow for no dates to be declared
    if( dates is None ):
        dateRange_dayNum_zeroHr = np.array((2077, np.int16(np.round(np.median(data_timeRef/86400))))); #recreate the "zero hr" day 
    else:
        dateRange_dayNum_zeroHr = dates['date range zero hr dayNum'];
    #END IF
    
    
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
    ax.axis('off'); #cartopy makes a new axis
    ax = plt.axes(projection=settings_map['projection']); #redefine the axis to be a geographical axis [AFTER cax else there's ERRORS oof]

    cax = ax.inset_axes((1.12, 0, 0.04, 1)); #make a color bar axis
    # divider = make_axes_locatable(ax); #prep to add an axis
    # cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #---ADD GRID LINES, SET PLOTTING AREA---
    if( settings_map['projection name'] == 'npstere' ):
        #Cartopy has a long way to go to match basemap, wow
        import matplotlib.path as mpath
        ax.set_aspect('equal');
        map_lat_autoTick = 15;
        map_long_autoTick = 30;
        plotLongTicks = np.arange(np.min(plotLongRange)+map_long_autoTick,np.max(plotLongRange)+map_long_autoTick,map_long_autoTick); #get the long ticks, will reuse
        plotLatTicks = np.arange(np.min(plotLatRange),np.max(plotLatRange),map_lat_autoTick); #get the lat ticks, will reuse
        # no labels, they're abs borked with the calls to limit the plot range
        gl = ax.gridlines(xlocs=plotLongTicks,ylocs=plotLatTicks,
            linewidth=PLOT_lineWidthSmoller, color='xkcd:black', alpha=0.5, linestyle='--', draw_labels=False, zorder=90); #draw some well-described gridlines
        
        #--- following based on https://stackoverflow.com/a/61986546/2403531 ---
        #this makes it be the lat range limits but makes it square
        ax.set_extent(plotLongRange + np.flip(plotLatRange).tolist(),ccrs.PlateCarree()); #set the extent, only platecarree makes anything remotely work with this call
        
        #this makes it not square again b/c platecarree makes it square
        theta = np.linspace(0, 2*np.pi, 360);
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T;
        center, radius = [0.5, 0.5], 0.5;
        circle = mpath.Path(verts * radius + center);
        ax.set_boundary(circle, transform=ax.transAxes);  #without this; get rect bound
        
        
        # for label alignment
        # alignmentVert = 'center' # also bottom, top
        # alignmentHoriz = 'center' # right, left
        
        #-longitude labels-
        for i in range(0,plotLongTicks.size):
            plotLatTemp = np.min(plotLatRange); #may need minor adjustments to keep verything nice
            if( plotLongTicks[i] == 180 ):
                alignmentHoriz = 'center';
                alignmentVert = 'bottom';
            elif( plotLongTicks[i] > 0 ):
                alignmentHoriz = 'left';
                alignmentVert = 'center';
                if( plotLongTicks[i] < 45 ):
                    plotLatTemp += -0.85; #give it a bump
                elif( (plotLongTicks[i] >= 45) & (plotLongTicks[i] < 75) ):
                    plotLatTemp += -0.5; #give it a bump
                elif( (plotLongTicks[i] >= 105) & (plotLongTicks[i] < 135) ):
                    plotLatTemp += -0.5; #give it a bump
                elif( (plotLongTicks[i] >= 135) & (plotLongTicks[i] < 165) ):
                    plotLatTemp += -0.85; #give it a bump
                #END IF
            elif( plotLongTicks[i] < 0 ):
                alignmentHoriz = 'right';
                alignmentVert = 'center';
                if( plotLongTicks[i] > -60 ):
                    plotLatTemp += -1; #give it a bump
                #END IF
            elif(  plotLongTicks[i] == 0 ):
                alignmentHoriz = 'center';
                alignmentVert = 'top';
                plotLatTemp += -0.15; #smol adjust
            #END IF
            long_proj, lat_proj = ax.projection.transform_point(plotLongTicks[i], plotLatTemp, ccrs.Geodetic()); #puts the longitude labels at the minimum latitude
            ax.text(long_proj, lat_proj, str(plotLongTicks[i])+u'\u00B0', \
                va=alignmentVert, ha=alignmentHoriz, color='xkcd:black',zorder=100);
        #END FOR i
        
        #-latitude labels-
        for i in range(0,plotLatTicks.size):
            long_proj, lat_proj = ax.projection.transform_point(135, plotLatTicks[i], ccrs.Geodetic()); #puts the latitude marks on the -45 longitude line
            ax.text(long_proj, lat_proj, str(plotLatTicks[i])+u'\u00B0', \
                va='center', ha='center', color='xkcd:black',zorder=100);
        #END FOR i
    else:
        ax.set_aspect('auto');
        ax.set_xticks(np.arange(np.min(plotLongRange),np.max(plotLongRange)+map_long_autoTick,map_long_autoTick),crs=settings_map['projection']); #gotta plot ticks with this to get -90 and 90
        ax.set_yticks(np.arange(np.min(plotLatRange),np.max(plotLatRange)+map_lat_autoTick,map_lat_autoTick),crs=settings_map['projection']); #gotta plot ticks with this to get -90 and 90
        # gl.xformatter = LONGITUDE_FORMATTER
        # gl.yformatter = LATITUDE_FORMATTER
        # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
        # gl.xlocator = tick.FixedLocator(np.arange(np.min(plotLongRange),np.max(plotLongRange)+map_long_autoTick,map_long_autoTick)); #this works ok, but be consistent use set_xticks
        # gl.ylocator = tick.FixedLocator(np.arange(np.min(plotLatRange),np.max(plotLatRange)+map_lat_autoTick,map_lat_autoTick)); #this doesn't plot -90 and 90 labels, but is req to get the gridlines right
        ax.set_extent(plotLongRange + plotLatRange); #set the plot extent, set at end - x and y ticks can extend plot area so this'll reign it in
    #END IF
    
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
    ax.add_feature(cfeature.OCEAN, facecolor=settings_map['water color'],alpha=0.5,zorder=0)
    ax.add_feature(cfeature.LAND, facecolor=settings_map['land color'],alpha=0.75,zorder=0)
    #reinforce axis limits after this drawing stuff
    if( settings_map['projection name'] != 'mill' ):
        ax.set_xlim(plotLongRange); #set x limits
        ax.set_ylim(plotLatRange); #set y limits
    #END IF
    
    #---ACTUALLY PLOT REAL STUFF HERE---
    if( FLG_useAllTimes == False ):
        (_, counts) = np.unique(data_time, return_counts=True); #get the counts of the data
        # timeIndex = np.where( counts >= data_time.size/data_timeUnique.size )[0][0]; #get an index where there's a lot of data
        if( np.all(counts == counts[0]) ):
            data_mean = np.mean(data_data.reshape(data_data.size//counts[0],counts[0]),axis=1); #get the mean at each time stamp
            timeIndex = np.where( np.min(np.abs(data_mean - np.mean(data_mean))) == np.abs(data_mean - np.mean(data_mean)) )[0][0]; #get an index where there's a lot of data
        else:
            timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
        #END IF
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
        
        k = np.where( (data_time == data_timeUnique[timeIndex]) & (data_data >= 1) )[0]; #gets during a time period
    else:
        if( data_time.size > 20000 ):
            k = np.zeros(data_time.size, dtype=np.bool_); #use just some
            k_limiter = np.random.randint(0, high=data_time.size, size=20000, dtype=np.int64);
            k[k_limiter] = 1; #use just some
        else:
            k = np.ones(data_time.size, dtype=np.bool_); #use em all
        #END IF
        
        #for compatability later
        timeIndex = np.where( np.abs(np.int16(np.round(np.median(data_timeRef/86400))) - data_timeRef/86400) == np.min(np.abs(np.int16(np.round(np.median(data_timeRef/86400))) - data_timeRef/86400)) )[0][0]; #get the index closest
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
    #END IF
    
    if( settings_map['projection name'] != 'npstere' ):
        titleOffset = 1; #offset for the title
    else:
        #polar needs a title offset
        titleOffset = 1.05; #offset for the title
    #END IF
    
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['scatter size'],c=data_data[k],cmap=settings_dataSpecific['colormap'], vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),zorder=10,transform=ccrs.PlateCarree());
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_major_formatter(tickr.FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(settings_dataSpecific['name']+settings_dataSpecific['units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
        #END FOR tick
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu));
        # cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),11)); #create useful tick marks
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    else:
        im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['scatter size'],c=data_data[k],cmap=settings_dataSpecific['colormap'],zorder=10,transform=ccrs.PlateCarree());
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(settings_dataSpecific['name']+settings_dataSpecific['units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
        #END FOR tick
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #END IF
    temp_yTicks = np.insert(cbar.ax.get_yticks(),0,1); #get em, tack on a 1
    temp_yLim = np.asarray(cbar.mappable.get_clim()); #get em, convert to array to make a copy prob
    temp_yLim[0] = 1; #force the bottom number to be 1
    cbar.mappable.set_clim(vmin=temp_yLim[0], vmax=temp_yLim[-1]);
    cax.yaxis.set_ticks(temp_yTicks); #create useful tick marks
    
    if( FLG_fancyPlot == 0 ):
        string_Title = settings_dataSpecific['name']+' Scatter Plot | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' at '+str(timeHr)+':'+str(timeMin).zfill(2)+':'+str(timeSec).zfill(2); #create mecha title
        ax.set_title(string_Title,fontproperties=settings_plot['font title FM'],y=titleOffset); #set the title
    #END IF
    
    if( np.all(boxArea != False) ):
        #draw a box where the data will be gathered
        temp_mapCoords = ( np.hstack( [np.linspace(boxArea[0,1],boxArea[1,1],200) , \
            np.linspace(boxArea[1,1],boxArea[3,1],200) , \
            np.linspace(boxArea[3,1],boxArea[2,1],200) , \
            np.linspace(boxArea[2,1],boxArea[0,1],200)] ) , \
            np.hstack( [np.linspace(boxArea[0,0],boxArea[1,0],200) , \
            np.linspace(boxArea[1,0],boxArea[3,0],200) , \
            np.linspace(boxArea[3,0],boxArea[2,0],200) , \
            np.linspace(boxArea[2,0],boxArea[0,0],200)] ) ); #convert to the geographic map coords
        ax.plot( temp_mapCoords[0],  #X longitude arcdeg
            temp_mapCoords[1],  #Y latitude arcdeg
            c='xkcd:fuchsia', linewidth=PLOT_lineWidthThicc, transform=ccrs.PlateCarree(), zorder=500); #
    #END IF
                    
    # #plot a * where the ISR is
    # if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
    #     temp_mapCoords = (avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
    #     ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=settings_map['site marker type'], color=settings_map['site marker color'], markersize=settings_map['site marker size'], zorder=150);
    # #END IF
    
    figFitter(fig, tryLim = 2); #fit the fig fast
    if( FLG_fancyPlot != 0 ):
        fig.savefig(os.path.join(settings_paths['fancyPlots'],settings_dataSpecific['name'].replace(' ','-')+'_area_scatter.png')); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF