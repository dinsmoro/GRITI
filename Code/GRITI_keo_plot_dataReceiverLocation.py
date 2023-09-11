"""
GOAL: Plot area that average alg with average over
RD on 6/11/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

#!!! THIS IS CURRENTLY UNFEASIBLE WITH CURRENT DATA PRODUCT TYPE !!

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tickr
# import cartopy as cartopy #cartopy replaces basemap b/c it's updated
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
# import os
import datetime
from subfun_figFitter import figFitter
from GRITI_plotHelper_area_init import GRITI_plotHelper_area_init

def GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
        FLG_fancyPlot = 0):
    #----Start plotting-----
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo_dataReceiverLocation IN fancyPlot FOLDER'); #report since you won't see anything
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
    coordType = settings_map['coord type']; #unpack
    if( 'lat long words' in list(settings_dataSpecific.keys()) ): #this is a late addon, added in a way that isn't important if it isn't there
        plot_latLongWords = settings_dataSpecific['lat long words'];
    else:
        plot_latLongWords = False;
    #END IF
      
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
    
    #--- Make a cartopy map ---
    kk = np.where(np.int64(data_timeUnique[timeIndex]/86400) == dates['date range full dayNum'][:,1])[0].item(); #get where the year is gonna be
    time4mag_hr = np.int32(np.mod(data_timeUnique[timeIndex],86400)//3600); #get hours
    time4mag_min = np.int32(np.mod(data_timeUnique[timeIndex],86400)//60-time4mag_hr*60); #get the minutes
    time4mag_sec = np.int32(np.mod(data_timeUnique[timeIndex],86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
    time4mag = datetime.datetime(dates['date range full'][kk,0],dates['date range full'][kk,1],dates['date range full'][kk,2], \
                                 hour = time4mag_hr, minute = time4mag_min, second = time4mag_sec); #date time object for aacgmv2    
    fig, ax, cax = GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_map, settings_plot, FLG_fancyPlot, time4mag=time4mag, alt4mag=120.);
    
    # k = np.where( data_time == data_timeUnique[timeIndex])[0]; #gets during a time period
    
    # k = np.where( np.unique(np.vstack((data_lat,data_long)), return_index=False, return_inverse=False, return_counts=False, axis=0) );
    
    if( np.any(np.isinf(plotLimValu)) == False ):
        im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['keo scatter size'],c=data_data[k],cmap=settings_dataSpecific['keo colormap'], vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),zorder=100,transform=ccrs.PlateCarree());
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_major_formatter(tickr.FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(settings_dataSpecific['keo labels']+settings_dataSpecific['keo units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
        #END FOR tick
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu));
        # cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),11)); #create useful tick marks
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    else:
        im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['keo scatter size'],c=data_data[k],cmap=settings_dataSpecific['keo colormap'],zorder=100,transform=ccrs.PlateCarree());
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(settings_dataSpecific['keo labels']+settings_dataSpecific['keo units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
        #END FOR tick
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #END IF
    if( FLG_fancyPlot == 0 ):
        string_Title = settings_dataSpecific['keo labels']+' Keo w/ Angle='+str(np.round(keo_angle,2))+' deg, Width='+ \
            str(np.round(keo_width,2))+' arcdeg, Step #='+str(keo_N)+ \
            ', Step Width='+ \
            str(np.round(np.sqrt( (temp_Longs_up[0,0] - temp_Longs_up[0,1])**2 + (temp_Lats_up[0,0] - temp_Lats_up[0,1])**2 ),2))+ \
            ' arcdeg | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' at '+str(timeHr)+':'+str(timeMin)+':'+str(timeSec); #create mecha title
        ax.set_title(string_Title,fontproperties=settings_plot['font title FM']); #set the title
    #END IF
    if( plot_latLongWords == True ):
        if( coordType == 'geo' ):
            ax.set_xlabel('Longitude [arcdeg]');
            ax.set_ylabel('Latitude [arcdeg]');
        elif( coordType == 'mag' ):
            ax.set_xlabel('Longitude (Geomag) [arcdeg]');
            ax.set_ylabel('Latitude (Geomag) [arcdeg]');
        #END IF
    #END IF
    
    #draw a box where the data will be gathered
    lineNum = 200; # number of pts to split line into, if using non-mercator more may be needed to make it look less bad on curvy stuff
    temp_mapCoords = ( np.hstack( [np.linspace(keo_range[0,1],keo_range[1,1],lineNum) , \
        np.linspace(keo_range[1,1],keo_range[3,1],lineNum) , \
        np.linspace(keo_range[3,1],keo_range[2,1],lineNum) , \
        np.linspace(keo_range[2,1],keo_range[0,1],lineNum)] ) , \
        np.hstack( [np.linspace(keo_range[0,0],keo_range[1,0],lineNum) , \
        np.linspace(keo_range[1,0],keo_range[3,0],lineNum) , \
        np.linspace(keo_range[3,0],keo_range[2,0],lineNum) , \
        np.linspace(keo_range[2,0],keo_range[0,0],lineNum)] ) ); #convert to the geographic map coords
    ax.plot( temp_mapCoords[0],  #X longitude arcdeg
        temp_mapCoords[1],  #Y latitude arcdeg
        c='xkcd:fuchsia',linewidth=PLOT_lineWidthThicc, zorder=90, transform=ccrs.PlateCarree());
    
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
        temp_mapCoords = ( np.linspace( temp_Long_List[i,0],temp_Long_List[i,3],lineNum ) , \
            np.linspace( temp_Lat_List[i,0],temp_Lat_List[i,3],lineNum ) ); #convert to the geographic map coords
        
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c='xkcd:fuchsia',linewidth=PLOT_lineWidthSmol, zorder=90, transform=ccrs.PlateCarree());
    #END FOR i
    
    #plot a * where the ISR is
    if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
        temp_mapCoords = (avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=settings_map['site marker type'], color=settings_map['site marker color'], markersize=settings_map['site marker size'], zorder=150, transform=ccrs.PlateCarree());
    #END IF
    
    figFitter(fig); #fit the fig fast
    if( FLG_fancyPlot == 0 ):
        plt.show(); #req to make plot show up
    else:
        fig.savefig(settings_paths['fancyPlots']+'\\'+settings_dataSpecific['keo data type']+'_keo_dataReceiverLocation'+settings_plot['save file type']); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF