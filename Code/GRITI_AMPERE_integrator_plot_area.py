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
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
# import os
import datetime
import copy
import aacgmv2
from Code.subfun_figFitter import figFitter
from Code.GRITI_plotHelper_area_init import GRITI_plotHelper_area_init
from Code.GRITI_plotHelper_axisizerLatLong import GRITI_plotHelper_axisizerLatLong
from Code.GRITI_movieMaker_subfun_dataGridder import GRITI_movieMaker_subfun_dataGridder

def GRITI_AMPERE_integrator_plot_area(AMPERE_data, dates, settings_AMPERE, \
        settings_map, settings_plot, settings_paths, \
        AMPERE_integrateArea_time = False, plot_terrainDraw = False, 
        plot_latLongWords = True, FLG_gridded = True, 
        FLG_sunLoc = False, FLG_sunSpin = False, settings_config_tz = None,
        FLG_fancyPlot = 0):
    #----Start plotting-----
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: '+settings_AMPERE['data type']+'_integrate_area IN fancyPlot FOLDER'); #report since you won't see anything
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
    avgPt_coords = settings_map['site coords'];
    coordType = settings_map['coord type']; #unpack
    if( 'altitude' in settings_AMPERE ):
        AMPERE_alt = settings_AMPERE['altitude'];
    else:
        AMPERE_alt = 120.; #default, great for auroral zone stuff (like field aligned currents)
    #END IF
    
    if( settings_plot['save file type'].lower() != '.pdf' ):
        plot_rasterize = False; #don't need to bother
    else:
        plot_rasterize = True; #it needs it
    #END IF
    
    if( FLG_gridded == True ):
        FLG_wordColorDub = False; #prevents world from being colored in, can't see it with pcolor
    else:
        FLG_wordColorDub = True; #helps allow world to be colored in
    #END IF
        
    #---Make box to show integration boundary---
    if( settings_AMPERE['integrate method'] == 0 ):
        #regular
        integrateLatRange = plotLatRange;
        integrateLongRange = plotLongRange;
    elif( settings_AMPERE['integrate method'] == 1 ):
        #max is 90 now within plot range
        integrateLatRange = copy.deepcopy(plotLatRange);
        integrateLatRange[np.where(integrateLatRange == np.max(integrateLatRange))[0].item()] = 90;
        integrateLongRange = plotLongRange;
    elif( settings_AMPERE['integrate method'] == 2 ):
        #max is user-defined now within plot range
        integrateLatRange = copy.deepcopy(plotLatRange);
        integrateLatRange[np.where(integrateLatRange == np.max(integrateLatRange))[0].item()] = settings_AMPERE['integrate method lat val'];
        integrateLongRange = plotLongRange;
    elif( settings_AMPERE['integrate method'] == 3 ):
        #360 long, plot lat range
        integrateLatRange = plotLatRange;
        integrateLongRange = [-180,180];
    elif( settings_AMPERE['integrate method'] == 4 ):
        #user-defined min to 90 max within plot range
        integrateLatRange = copy.deepcopy(plotLatRange);
        integrateLatRange[np.where(integrateLatRange == np.max(integrateLatRange))[0].item()] = 90;
        integrateLatRange[np.where(integrateLatRange == np.min(integrateLatRange))[0].item()] = settings_AMPERE['integrate method lat val'];
        integrateLongRange = plotLongRange;
    elif( settings_AMPERE['integrate method'] == 5 ):
        #user-defined min to 90 max for 360 long
        integrateLatRange = copy.deepcopy(plotLatRange);
        integrateLatRange[np.where(integrateLatRange == np.max(integrateLatRange))[0].item()] = 90;
        integrateLatRange[np.where(integrateLatRange == np.min(integrateLatRange))[0].item()] = settings_AMPERE['integrate method lat val'];
        integrateLongRange = [-180,180];
    elif( settings_AMPERE['integrate method'] == 6 ):
        #user-defined max for 360 long
        integrateLatRange = copy.deepcopy(plotLatRange);
        integrateLatRange[np.where(integrateLatRange == np.max(integrateLatRange))[0].item()] = settings_AMPERE['integrate method lat val'];
        integrateLatRange[np.where(integrateLatRange == np.min(integrateLatRange))[0].item()] = 0;
        integrateLongRange = [-180,180];
    elif( settings_AMPERE['integrate method'] == 7 ):
        #radius party
        import sys
        print('add this')
        sys.crash();
    #END IF  
    integrate_range = np.concatenate( [ [np.max(integrateLatRange), np.min(integrateLatRange), np.max(integrateLatRange), np.min(integrateLatRange)], [np.max(integrateLongRange) , np.max(integrateLongRange), np.min(integrateLongRange), np.min(integrateLongRange)]]).reshape(2,4).T; #arcdeg, record pts for use
    
    #not super sure this actually fixes it, rather just makes it go through truncation idk been a while
    integrate_range[np.where(integrate_range[:,1] > 180),1]  = 180.; #fix over 180 so it flips over to the other side
    integrate_range[np.where(integrate_range[:,1] < -180),1]  = -180.; #fix less than -180 so it flips over to the other side
    integrate_range[np.where(integrate_range[:,0] > 90),0]  = 90.; #fix over 90 so it flips over to the other side
    integrate_range[np.where(integrate_range[:,0] < -90),0]  = -90.; #fix less than -90 so it flips over to the other side
      
    #---ACTUALLY PLOT REAL STUFF HERE---
    #count-based (useless for AMPERE)
    # (_, counts) = np.unique(AMPERE_data['time'], return_counts=True); #get the counts of the data
    # timeIndex = np.where( counts >= AMPERE_data['time'].size/AMPERE_data['time unique'].size )[0][0]; #get an index where there's a lot of data
    #integration based (useful for AMPERE)
    # timeIndex = np.where( AMPERE_data['integrated'] == AMPERE_data['integrated'].max() )[0][0]; #get an index where there's a lot
    timeIndex = np.where( np.abs(-15.4*3600 - (AMPERE_data['time unique']-dateRange_dayNum_zeroHr[1]*86400)) == np.min(np.abs(-15.4*3600 - (AMPERE_data['time unique']-dateRange_dayNum_zeroHr[1]*86400))) )[0][0]; #get an index where there's a lot of data
    #time index override (useful for AMPERE)
    # if( AMPERE_integrateArea_time != False ):
    #     timeRef = AMPERE_integrateArea_time; #use supplied
    # else:
    #     timeRef = AMPERE_data['time unique'][0]; #just choose 1st one yolo
    # #END IF
    # timeIndex = np.where( np.abs(timeRef - (AMPERE_data['time unique']-dateRange_dayNum_zeroHr[1]*86400)) == np.min(np.abs(timeRef - (AMPERE_data['time unique']-dateRange_dayNum_zeroHr[1]*86400))) )[0][0]; #get an index where there's a lot of data
    timeHr = np.int16((AMPERE_data['time unique'][timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600); #calc the hour
    timeMin = np.int16(np.round((timeHr - (AMPERE_data['time unique'][timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600)*60)); #calc the minute
    timeSec = np.int16(np.round(((timeHr - (AMPERE_data['time unique'][timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600)*60 - timeMin)*60)); #calc the second
    if( timeSec == 60 ):
        timeSec = 0; #roll bak to  0
        timeMin += 1; #increment by 1
    #END IF
    if( timeMin == 60 ):
        timeMin = 0; #roll back to 0
        timeHr += 1; #increment by 1
    #END IF
    
    #--- Make a cartopy map ---
    kk = np.where(np.int64(AMPERE_data['time unique'][timeIndex]/86400) == dates['date range full dayNum'][:,1])[0].item(); #get where the year is gonna be
    time4mag_hr = np.int32(np.mod(AMPERE_data['time unique'][timeIndex],86400)//3600); #get hours
    time4mag_min = np.int32(np.mod(AMPERE_data['time unique'][timeIndex],86400)//60-time4mag_hr*60); #get the minutes
    time4mag_sec = np.int32(np.mod(AMPERE_data['time unique'][timeIndex],86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
    time4mag = datetime.datetime(dates['date range full'][kk,0],dates['date range full'][kk,1],dates['date range full'][kk,2], \
                                 hour = time4mag_hr, minute = time4mag_min, second = time4mag_sec); #date time object for aacgmv2
    
    settings_map = copy.deepcopy(settings_map); #deep copy
    settings_map['projection'] = ccrs.NorthPolarStereo(); #always use this
    settings_map['projection name'] = 'npstere'; #no south support yet...
    # settings_map['coord type'] = settings_AMPERE['integrate method coord type'];
    plotLatRange = [40, 90]; #redefine
    plotLongRange = [-180, 180]; #redefine
    fig, ax, cax = GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_map, settings_plot, FLG_fancyPlot, time4mag=time4mag, alt4mag=AMPERE_alt, FLG_wordColorDub=FLG_wordColorDub);
    
    if( plot_terrainDraw == True ):
        prevAspect = ax.get_aspect();
        ax.stock_img();
        ax.set_aspect(prevAspect); #reset aspect ratio
        import cartopy.feature as cfeature
        ax.add_feature(cfeature.BORDERS, edgecolor='xkcd:black', zorder=45);
        # ax.coastlines(resolution=mapper_resolution, color='xkcd:black',zorder=75); #draw the coastlines
        #none of the stamen stuff worked at all
        # from cartopy.io.img_tiles import Stamen
        # import copy
        # stamen_terrain = Stamen('terrain'); #get the thing fired up
        # settings_map = copy.deepcopy(settings_map); #copy map settings
        # settings_map['projection'] = stamen_terrain.crs; #it has to be this apparently for it to show up
        # fig, ax, cax = GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_map, settings_plot, FLG_fancyPlot, time4mag=time4mag, alt4mag=AMPERE_alt);
        # ax.add_image(stamen_terrain, 3, zorder=0); #draw it on
    #END IF
    
    
    
    k = AMPERE_data['time'] == AMPERE_data['time unique'][timeIndex]; #gets during a time period
    
    if( np.all(np.isinf(settings_AMPERE['plot lim'])) == False ):
        if( np.any(np.isinf(settings_AMPERE['plot lim'])) ): #deal with infs
            settings_AMPERE = copy.deepcopy(settings_AMPERE); #deep copy
            if( np.isinf(settings_AMPERE['plot lim'][0]) ):
                if( np.sign(settings_AMPERE['plot lim'][0]) > 0 ):
                    settings_AMPERE['plot lim'][0] = np.max(AMPERE_data[settings_AMPERE['data type']][k]);
                else:
                    settings_AMPERE['plot lim'][0] = np.min(AMPERE_data[settings_AMPERE['data type']][k]);
                #END IF
            #END IF
            if( np.isinf(settings_AMPERE['plot lim'][1]) ):
                if( np.sign(settings_AMPERE['plot lim'][1]) > 0 ):
                    settings_AMPERE['plot lim'][1] = np.max(AMPERE_data[settings_AMPERE['data type']][k]);
                else:
                    settings_AMPERE['plot lim'][1] = np.min(AMPERE_data[settings_AMPERE['data type']][k]);
                #END IF
            #END IF
        #END IF
        if( np.all(np.asarray(settings_AMPERE['plot lim']) >= 0) ): #indicates positive-only so enforce lower limit to prevent plotting (does not activate for symmetric-around-0 data)
            k = np.where(k & (AMPERE_data[settings_AMPERE['data type']] >= np.min(settings_AMPERE['plot lim'])))[0]; #enforce lower limit, get as indexes to improve speed
        else:
            k = np.where(k)[0]; #just get indexes b/c faster
        #END IF
        if( FLG_gridded == True ):
            AMPERE_lat_delta = np.median(np.diff(np.unique(AMPERE_data['lat'])));
            if( np.isclose(AMPERE_lat_delta,np.int64(AMPERE_lat_delta)) ):
                AMPERE_lat_delta = np.int64(AMPERE_lat_delta); #convert to integer if it's an integer
            elif( AMPERE_lat_delta < 1e-4 ):
                AMPERE_lat_delta = 1; #override
            #END IF
            AMPERE_long_delta = np.median(np.diff(np.unique(AMPERE_data['long'])));
            if( np.isclose(AMPERE_long_delta,np.int64(AMPERE_long_delta)) ):
                AMPERE_long_delta = np.int64(AMPERE_long_delta); #convert to integer if it's an integer
            elif( AMPERE_long_delta < 1e-4 ):
                AMPERE_long_delta = 15; #override
            #END IF
            gif_Grid_Lat_AMPERE = np.arange(np.min(plotLatRange),np.max(plotLatRange)+AMPERE_lat_delta,AMPERE_lat_delta); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
            gif_Grid_Long_AMPERE = np.arange(np.min(plotLongRange),np.max(plotLongRange)+AMPERE_long_delta,AMPERE_long_delta); #degc, create long points specilized for AMPERE
                
            gif_Grid_AMPERE = GRITI_movieMaker_subfun_dataGridder(AMPERE_data['lat'][k], AMPERE_data['long'][k], AMPERE_data[settings_AMPERE['data type']][k], gif_Grid_Lat_AMPERE,gif_Grid_Long_AMPERE,gif_Grid_Lat_AMPERE.size-1,gif_Grid_Long_AMPERE.size-1,AMPERE_lat_delta,AMPERE_long_delta,2,101,8).T; #101 disables the data rejection stuff b/c AMPERE doesn't need it
            
            AMPERE_latLongMapped = np.meshgrid( gif_Grid_Long_AMPERE, gif_Grid_Lat_AMPERE); #helps the pcolor work
            im = ax.pcolormesh(AMPERE_latLongMapped[0], AMPERE_latLongMapped[1], gif_Grid_AMPERE, vmin=np.min(settings_AMPERE['plot lim']), vmax=np.max(settings_AMPERE['plot lim']), cmap=settings_AMPERE['colormap'], zorder=100, transform=ccrs.PlateCarree(), rasterized=plot_rasterize); # pseudocolor plot "stretched" to the grid
        else:
            im = ax.scatter(AMPERE_data['long'][k], AMPERE_data['lat'][k], s=settings_AMPERE['scatter size'], c=AMPERE_data[settings_AMPERE['data type']][k], cmap=settings_AMPERE['colormap'], vmin=np.min(settings_AMPERE['plot lim']), vmax=np.max(settings_AMPERE['plot lim']), zorder=100,transform=ccrs.PlateCarree());
        #END IF
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        # cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f')); #force a rounded format
        if( np.all(np.mod(cbar.get_ticks(),1) == 0) ):
            cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.0f')); #force a rounded format
        else:
            cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
        #END IF
        cbar.set_label(settings_AMPERE['labels'][settings_AMPERE['data type']]+settings_AMPERE['units'][settings_AMPERE['data type']]); #tabel the colorbar
        cbar.ax.tick_params(labelsize=settings_plot['font axis tick']);
        cbar.mappable.set_clim(vmin=np.min(settings_AMPERE['plot lim']), vmax=np.max(settings_AMPERE['plot lim']));
        # cax.yaxis.set_ticks(np.linspace(np.min(settings_AMPERE['plot lim']),np.max(settings_AMPERE['plot lim']),5)); #create useful tick marks
        # cax.yaxis.set_ticks(np.linspace(np.min(settings_AMPERE['plot lim']),np.max(settings_AMPERE['plot lim']),11)); #create useful tick marks
        
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
        
        if( np.all(np.asarray(settings_AMPERE['plot lim']) >= 0) ):
            #ensure min tick is set
            caxTicks = cbar.get_ticks();
            caxTicks[0] = np.min(settings_AMPERE['plot lim']); #enforce
            cbar.set_ticks(caxTicks);
        #END IF
    else:
        if( FLG_gridded == True ):            
            AMPERE_lat_delta = np.median(np.diff(np.unique(AMPERE_data['lat'])));
            if( np.isclose(AMPERE_lat_delta,np.int64(AMPERE_lat_delta)) ):
                AMPERE_lat_delta = np.int64(AMPERE_lat_delta); #convert to integer if it's an integer
            elif( AMPERE_lat_delta < 1e-4 ):
                AMPERE_lat_delta = 1; #override
            #END IF
            AMPERE_long_delta = np.median(np.diff(np.unique(AMPERE_data['long'])));
            if( np.isclose(AMPERE_long_delta,np.int64(AMPERE_long_delta)) ):
                AMPERE_long_delta = np.int64(AMPERE_long_delta); #convert to integer if it's an integer
            elif( AMPERE_long_delta < 1e-4 ):
                AMPERE_long_delta = 15; #override
            #END IF
            gif_Grid_Lat_AMPERE = np.arange(np.min(plotLatRange),np.max(plotLatRange)+AMPERE_lat_delta,AMPERE_lat_delta); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
            gif_Grid_Long_AMPERE = np.arange(np.min(plotLongRange),np.max(plotLongRange)+AMPERE_long_delta,AMPERE_long_delta); #degc, create long points specilized for AMPERE
               
            gif_Grid_AMPERE = GRITI_movieMaker_subfun_dataGridder(AMPERE_data['lat'][k],AMPERE_data['long'][k],AMPERE_data[settings_AMPERE['data type']][k],gif_Grid_Lat_AMPERE,gif_Grid_Long_AMPERE,gif_Grid_Lat_AMPERE.size-1,gif_Grid_Long_AMPERE.size-1,AMPERE_lat_delta,AMPERE_long_delta,2,101,8).T; #101 disables the data rejection stuff b/c AMPERE doesn't need it
            
            AMPERE_latLongMapped = np.meshgrid( gif_Grid_Long_AMPERE, gif_Grid_Lat_AMPERE); #helps the pcolor work
            im = ax.pcolormesh(AMPERE_latLongMapped[0], AMPERE_latLongMapped[1], gif_Grid_AMPERE, cmap=settings_AMPERE['colormap'], zorder=100, transform=ccrs.PlateCarree(), rasterized=plot_rasterize); # pseudocolor plot "stretched" to the grid
        else:
            im = ax.scatter(AMPERE_data['long'][k],AMPERE_data['lat'][k],s=settings_AMPERE['scatter size'],c=AMPERE_data[settings_AMPERE['data type']][k],cmap=settings_AMPERE['colormap'],zorder=100,transform=ccrs.PlateCarree());
        #END IF
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(settings_AMPERE['labels'][settings_AMPERE['data type']]+settings_AMPERE['units'][settings_AMPERE['data type']]); #tabel the colorbar
        cbar.ax.tick_params(labelsize=settings_plot['font axis tick']);
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #END IF
    if( FLG_fancyPlot == 0 ):
        string_Title = 'Integration Method '+str(settings_AMPERE['integrate method'])+' | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' @ '+str(timeHr)+':'+str(timeMin).zfill(2)+':'+str(timeSec).zfill(2); #create mecha title
        if( 'stere' not in settings_map['projection name'] ):
            ax.set_title(string_Title,fontproperties=settings_plot['font title FM']); #set the title
        else:
            ax.set_title(string_Title,fontproperties=settings_plot['font title FM'],y=1.075); #set the title
        #END IF
    #END IF
    if( 'stere' not in settings_map['projection name'] ):
        if( plot_latLongWords == True ):
            if( coordType == 'geo' ):
                ax.set_xlabel('Longitude [arcdeg]');
                ax.set_ylabel('Latitude [arcdeg]');
            elif( coordType == 'mag' ):
                ax.set_xlabel('Longitude (Geomag) [arcdeg]');
                ax.set_ylabel('Latitude (Geomag) [arcdeg]');
            #END IF
        #END IF
    #END IF
    
    #draw a box where the data will be gathered
    lineNum = 200; # number of pts to split line into, if using non-mercator more may be needed to make it look less bad on curvy stuff
    if( ('stere' in settings_map['projection name']) & ((integrate_range[1,1]-integrate_range[3,1]) == 360) ):
        temp_mapCoords = [ 
            np.linspace(integrate_range[1,1],integrate_range[3,1],lineNum) , \
            np.linspace(integrate_range[2,0],integrate_range[0,0],lineNum) ]; #convert to the geographic map coords
    else:
        temp_mapCoords = [ np.hstack( [np.linspace(integrate_range[0,1],integrate_range[1,1],lineNum) , \
            np.linspace(integrate_range[1,1],integrate_range[3,1],lineNum) , \
            np.linspace(integrate_range[3,1],integrate_range[2,1],lineNum) , \
            np.linspace(integrate_range[2,1],integrate_range[0,1],lineNum)] ) , \
            np.hstack( [np.linspace(integrate_range[0,0],integrate_range[1,0],lineNum) , \
            np.linspace(integrate_range[1,0],integrate_range[3,0],lineNum) , \
            np.linspace(integrate_range[3,0],integrate_range[2,0],lineNum) , \
            np.linspace(integrate_range[2,0],integrate_range[0,0],lineNum)] ) ]; #convert to the geographic map coords
    #END IF
        
    if( settings_AMPERE['integrate method coord type'] != settings_map['coord type'] ):
        if( settings_AMPERE['integrate method coord type'] != 'mag' ):
            aacgmv2_method = 'G2A'; #geo to mag
        else:
            aacgmv2_method = 'A2G'; #mag to geo
        #END IF
        for jj in range(0,temp_mapCoords[0].size):
            [temp_mapCoords[1][jj], temp_mapCoords[0][jj], _] = aacgmv2.convert_latlon(temp_mapCoords[1][jj], temp_mapCoords[0][jj], AMPERE_alt, time4mag, method_code=aacgmv2_method); #converts from geographic to geomagnetic (AACGMv2) or vice-versa
        #END FOR jj
        if( ((integrate_range[1,1]-integrate_range[3,1]) == 360) ):
            kj = np.where(np.abs(np.diff(temp_mapCoords[0])) > (np.mean(np.abs(np.diff(temp_mapCoords[0]))) + np.std(np.abs(np.diff(temp_mapCoords[0])))*3))[0]+1; #get outliers (+1 for diff)
            temp_mapCoords[0] = np.insert(temp_mapCoords[0], kj, np.nan); #nan the gap so it doesn't get weird
            temp_mapCoords[1] = np.insert(temp_mapCoords[1], kj, np.nan); #nan the gap so it doesn't get weird
        #END IF
    #END IF
    if( isinstance(settings_AMPERE['colormap'],str) ):
        areaBounder_color = 'xkcd:fuchsia'; #stands out
        loc_color = settings_map['site marker color']; #default
    else:
        if( np.any(np.all(np.abs(settings_AMPERE['colormap'].colors - np.array([237,17,217])/255) < 0.15, axis=1)) ):
            areaBounder_color = 'xkcd:vermillion'; #if colormap uses fuschia-like color use redish instead b/c looks different to colorblind
            loc_color = 'xkcd:dark teal'; #not more purple
        else:
            areaBounder_color = 'xkcd:fuchsia'; #stands out
            loc_color = settings_map['site marker color']; #default
       #END IF
    #END IF
    ax.plot( temp_mapCoords[0],  #X longitude arcdeg
        temp_mapCoords[1],  #Y latitude arcdeg
        c=areaBounder_color,linewidth=PLOT_lineWidthThicc, zorder=180, transform=ccrs.PlateCarree());
    # if( ('stere' in settings_map['projection name']) & ((integrate_range[1,1]-integrate_range[3,1]) == 360) ):
    #     ax.fill( temp_mapCoords[0],  #X longitude arcdeg
    #         temp_mapCoords[1],  #Y latitude arcdeg
    #         c=areaBounder_color, facecolor=areaBounder_color, edgecolor=areaBounder_color, alpha=0.5, zorder=180, transform=ccrs.PlateCarree());
    # #END IF
    
    #plot a * where the ISR is
    if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
        temp_mapCoords = (avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=settings_map['site marker type'], color=loc_color, markersize=settings_map['site marker size'], zorder=150, transform=ccrs.PlateCarree());
        #END IF
    #END IF
    
    FLG_sunLoc = True
    if( FLG_sunLoc == True ):
        from Code.subfun_sunAlsoRises_location import sunAlsoRises_location
        sunSubSolar_loc = sunAlsoRises_location(np.expand_dims(np.array( (dateRange_dayNum_zeroHr[0], np.int16(AMPERE_data['time unique'][timeIndex]/86400)) ).astype(np.int16),axis=0),timeIndexes=AMPERE_data['time unique'][timeIndex],timeZone='UTC');
        sunSubSolar_loc['lat'] = sunSubSolar_loc['lat'][0]; #undo
        sunSubSolar_loc['long'] = sunSubSolar_loc['long'][0];
        
        if( settings_map['coord type'] == 'mag' ):
            [sunSubSolar_loc['lat'], sunSubSolar_loc['long'], _] = aacgmv2.convert_latlon(sunSubSolar_loc['lat'], sunSubSolar_loc['long'], AMPERE_alt, time4mag, method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
        #END IF

        if( FLG_sunSpin == 0):
            x = (1.05*0.5*np.cos((sunSubSolar_loc['long']-90)*np.pi/180))+0.5; #geoMap coordinate 
            y = (1.05*0.5*np.sin((sunSubSolar_loc['long']-90)*np.pi/180))+0.5; #geoMap coordinate 
            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=FONT_axisTickFM);
            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
            hSun = ax.add_artist(circleSun); #plot the sun
        else:
            #this does not work in cartopy!! so don't
            hSun = ax.set_theta_offset((sunSubSolar_loc['long']-90)*np.pi/180); #turn the whole plot so top is where the sun is
            #I'm not sure if I can get this working easily - not a lot of optons.
        #END IF
    #END IF
    
    if( 'stere' not in settings_map['projection name'] ):
        GRITI_plotHelper_axisizerLatLong(plotLatRange,ax=ax,axDir='y',tickNumGoal=13,tickReducer=0,FLG_extendLims=False); #auto tick the thing
        GRITI_plotHelper_axisizerLatLong(plotLongRange,ax=ax,axDir='x',tickReducer=3.5,FLG_extendLims=False);
    #END IF
    
    figFitter(fig); #fit the fig fast
    if( FLG_fancyPlot != 0 ):
        fig.savefig(os.path.join(settings_paths['fancyPlots'],settings_AMPERE['data type']+'_integrate_area'+settings_plot['save file type'])); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF