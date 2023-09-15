"""
GOAL: Plot data density that average alg with average over
RD on 6/11/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.ticker as tick
# import cartopy as cartopy #cartopy replaces basemap b/c it's updated
# import os
from Code.GRITI_plotHelper_area_init import GRITI_plotHelper_area_init
from Code.subfun_figFitter import figFitter

def GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
        FLG_fancyPlot = 0):
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: '+settings_dataSpecific['keo data type']+'_keo_area_dataDensity IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    
    #----Start plotting-----
    #Unpack line widths
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
    # dateRange_dayNum_zeroHr = dates['date range zero hr dayNum'];
    # keo_angle = settings_dataSpecific['keo angle'];
    # keo_width = settings_dataSpecific['keo width'];
    keo_N = settings_dataSpecific['keo N'];
    # plotLimValu = settings_dataSpecific['keo plot lim'];
    avgPt_coords = settings_map['site coords'];
    
    #--- Make a cartopy map ---
    fig, ax, cax = GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_map, settings_plot, FLG_fancyPlot);
    
    #--- Get info on how big the plot is ---
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get info on the size of the plot area to know what geographic scale to use
    mapper_resolution = np.max( [np.abs(plotLatRange[0]-plotLatRange[1])/bbox.height , np.abs(plotLongRange[0]-plotLongRange[1])/bbox.width] ); #degc, max extent covered in the plot
    if( mapper_resolution > 20 ): #arbitrary numbers
        # mapper_resolution = '110m'; #the resolution to use for plotting geographical features
        binner_scaler = 0.5; #scaler
    elif( mapper_resolution > 0.5 ): #arbitrary numbers
        # mapper_resolution = '50m'; #the resolution to use for plotting geographical features
        binner_scaler = 1; #scaler
    else:
        #otherwise if the deg/in for the plot is super small use the highest detail possible
        # mapper_resolution = '10m'; #the resolution to use for plotting geographical features
        binner_scaler = 5; #scaler
    #END IF
    
    #---ACTUALLY PLOT REAL STUFF HERE---    
    # im = ax.scatter(data_long,data_lat,s=20,facecolors='none', edgecolors='xkcd:blue green');
    keoDensity, keoLong, keoLat = np.histogram2d(data_long, data_lat, bins=(np.int64(np.diff(plotLongRange).item()*binner_scaler), np.int64(np.diff(plotLatRange).item()*binner_scaler)) );
    keoDensity[keoDensity == 0] = np.nan; #set 0's to NaN so they don't get colored in (more obvious no data)
    pltHelprX, pltHelprY = np.meshgrid( keoLong, keoLat);
    from matplotlib.axes import Axes
    from cartopy.mpl.geoaxes import GeoAxes
    GeoAxes._pcolormesh_patched = Axes.pcolormesh; #cartopy keeps getting rekt by reimplementing matplotlib, why is this the new way
    im = ax.pcolormesh(pltHelprX, pltHelprY,  keoDensity.T ,cmap=settings_dataSpecific['keo colormap']); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar.set_label('Data Density [# Pts. Per Pixel]'); #tabel the colorbar
    for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
        tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
    #END FOR tick
    cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    cbarLims = cbar.mappable.get_clim(); #get the current colorbar lims
    cbar.mappable.set_clim(vmin=1, vmax=np.max(cbarLims)); #force 1 for lower end
    cbarTicks = cax.get_yticks(); #get the current ticks
    cax.yaxis.set_ticks(np.hstack( ((1,),cbarTicks) )); #create useful tick marks with 1 at the end
    
    if( FLG_fancyPlot == 0 ):
        string_Title = settings_dataSpecific['keo data type']+' Data Density from '+ \
            str(dates['date range'][0,1])+'/'+str(dates['date range'][0,2])+'/'+str(dates['date range'][0,0])+' to '+ \
            str(dates['date range'][-1,1])+'/'+str(dates['date range'][-1,2])+'/'+str(dates['date range'][-1,0])+' (M/D/YR)'; #create mecha title
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
        fig.savefig(os.path.join(settings_paths['fancyPlots'],settings_dataSpecific['keo data type']+'_keo_area_dataDensity.png')); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF