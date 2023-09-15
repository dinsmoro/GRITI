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
# import cartopy as cartopy #cartopy replaces basemap b/c it's updated
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
import os
import copy
import datetime
from Code.GRITI_movieMaker_subfun_dataGridder import GRITI_movieMaker_subfun_dataGridder
from Code.subfun_figFitter import figFitter
from Code.GRITI_plotHelper_area_init import GRITI_plotHelper_area_init
from Code.GRITI_plotHelper_axisizerLatLong import GRITI_plotHelper_axisizerLatLong

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
    coordType = settings_map['coord type']; #unpack
    if( 'lat long words' in list(settings_dataSpecific.keys()) ): #this is a late addon, added in a way that isn't important if it isn't there
        plot_latLongWords = settings_dataSpecific['lat long words'];
        plot_latLongUnitsOnLabels = True;
    else:
        plot_latLongWords = False;
        plot_latLongUnitsOnLabels = True;
    #END IF
    if( 'pierceAlt' in settings_dataSpecific ):
        keo_alt = settings_dataSpecific['pierceAlt'];
    else:
        keo_alt = 120.; #default, great for auroral zone stuff (like field aligned currents)
    #END IF
    if( 'keo gridder' in settings_dataSpecific ):
        FLG_gridder = settings_dataSpecific['keo gridder']; #use this
    else:
        FLG_gridder = False;
    #END IF
    if( 'keo coord type' not in settings_dataSpecific ):
        settings_dataSpecific = copy.deepcopy(settings_dataSpecific); #deep copy to prevent editing previous
        settings_dataSpecific['keo coord type'] = coordType; #set this for later
    #END IF
    if( 'terrain draw' in settings_dataSpecific ):
        FLG_wordColorDub = settings_dataSpecific['terrain draw']; #set this for later
    else:
        FLG_wordColorDub = False;
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
            latLong_lat_dirAdder = settings_map['indicate direction'][1]['lat']; #get the lat
            latLong_long_dirAdder = settings_map['indicate direction'][1]['long']; #get the long
        else:
            latLong_lat_dirAdder = ''; #nada
            latLong_long_dirAdder = ''; #nada
        #END IF
    else:
        latLong_lat_dirAdder = ''; #nada
        latLong_long_dirAdder = ''; #nada
    #END IF
    
    if( settings_plot['save file type'].lower() != '.pdf' ):
        plot_rasterize = False; #don't need to bother
    else:
        plot_rasterize = True; #it needs it
    #END IF
      
    #---ACTUALLY PLOT REAL STUFF HERE---
    (_, counts) = np.unique(data_time, return_counts=True); #get the counts of the data
    if( np.all(counts == counts[0]) ): #this means that its a set grid of data so counts not useful
        timeIndex = np.nansum(data_data.reshape(-1,counts[0]), axis=1); #reuse illegal style, since it's a consistent grid regrid it so that sum can happen
        timeIndex = np.where( timeIndex == np.max(timeIndex) )[0][0]; #get an index where there's a lot of intensity
    else:
        # timeIndex = np.where( counts >= data_time.size/data_timeUnique.size )[0][0]; #get an index where there's a lot of data
        timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
    #END IF
    # timeIndex = np.where( np.abs(-15.4*3600 - (data_timeUnique-dateRange_dayNum_zeroHr[1]*86400)) == np.min(np.abs(-15.4*3600 - (data_timeUnique-dateRange_dayNum_zeroHr[1]*86400))) )[0][0]; #get an index where there's a lot of data
    timeIndex = np.where( np.abs(104640 - (data_timeUnique-dateRange_dayNum_zeroHr[1]*86400)) == np.min(np.abs(104640 - (data_timeUnique-dateRange_dayNum_zeroHr[1]*86400))) )[0][0]; #get an index where there's a lot of data [29:04, didn't divide well]
    # timeIndex = np.where( np.abs(data_timeRef[0] - data_timeUnique) == np.min(np.abs(data_timeRef[0] - data_timeUnique)) )[0][0]; #get an index where there's a lot of data
    timeHr = np.int16((data_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600); #calc the hour
    timeMin = np.int16(np.round(((data_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)*60)); #calc the minute
    timeSec = np.int16(np.round((((data_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)*60 - timeMin)*60)); #calc the second
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
    if( 'keo polar mode' in settings_dataSpecific ):
        if( settings_dataSpecific['keo polar mode'] == 1 ):
            settings_mapCOPY = copy.deepcopy(settings_map); #deep copy to prevent editing previous
            settings_mapCOPY['projection'] = ccrs.NorthPolarStereo(); #always use this
            settings_mapCOPY['projection name'] = 'npstere'; #no south support yet...
            if( 'keo coord type' in settings_dataSpecific ):
                if( settings_dataSpecific['keo coord type'] == 'mag' ):
                    settings_mapCOPY['coord type'] = 'mag'; #set that
                #END IF
            #END IF
            plotLatRange = [40, 90]; #redefine
            plotLongRange = [-180, 180]; #redefine
        else:
            settings_mapCOPY = settings_map; #smooth out stuff with an alias
        #END IF
    else:
        settings_mapCOPY = settings_map; #smooth out stuff with an alias
    #END IF
    fig, ax, cax = GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_mapCOPY, settings_plot, FLG_fancyPlot, time4mag=time4mag, alt4mag=keo_alt, FLG_wordColorDub=FLG_wordColorDub);
    
    if( 'border draw' in settings_dataSpecific ):
        if( (settings_dataSpecific['border draw'] == True) & (settings_dataSpecific['keo coord type'] == 'geo') ):
            # ax.stock_img();
            # ax.set_aspect('auto'); #reset aspect ratio
            import cartopy.feature as cfeature
            ax.add_feature(cfeature.BORDERS, edgecolor='xkcd:black', zorder=75);
            # ax.coastlines(resolution=mapper_resolution, color='xkcd:black',zorder=75); #draw the coastlines
            #none of the stamen stuff worked at all
            # from cartopy.io.img_tiles import Stamen
            # import copy
            # stamen_terrain = Stamen('terrain'); #get the thing fired up
            # settings_mapCopy = copy.deepcopy(settings_map); #copy map settings
            # settings_mapCopy['projection'] = stamen_terrain.crs; #it has to be this apparently for it to show up
            # fig, ax, cax = GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_mapCopy, settings_plot, FLG_fancyPlot, time4mag=time4mag, alt4mag=keo_alt);
            # ax.add_image(stamen_terrain, 3, zorder=0); #draw it on
        #END IF
    #END IF
    
    k = data_time == data_timeUnique[timeIndex]; #gets during a time period
    
    if( np.all(np.isinf(plotLimValu)) == False ):
        if( np.any(np.isinf(plotLimValu)) ): #deal with infs
            plotLimValu = np.copy(plotLimValu); #copy
            if( np.isinf(plotLimValu[0]) ):
                if( np.sign(plotLimValu[0]) > 0 ):
                    plotLimValu[0] = np.max(data_data[k]);
                else:
                    plotLimValu[0] = np.min(data_data[k]);
                #END IF
            #END IF
            if( np.isinf(plotLimValu[1]) ):
                if( np.sign(plotLimValu[1]) > 0 ):
                    plotLimValu[1] = np.max(data_data[k]);
                else:
                    plotLimValu[1] = np.min(data_data[k]);
                #END IF
            #END IF
        #END IF
        if( np.all(np.asarray(plotLimValu) >= 0) ): #indicates positive-only so enforce lower limit to prevent plotting (does not activate for symmetric-around-0 data)
            k = np.where(k & (data_data >= np.min(plotLimValu)))[0]; #enforce lower limit, get as indexes to improve speed
        else:
            k = np.where(k)[0]; #just get indexes b/c faster
        #END IF

        if( FLG_gridder == True ):
            data_lat_delta = np.median(np.diff(np.unique(data_lat)));
            if( np.isclose(data_lat_delta,np.int64(data_lat_delta)) ):
                data_lat_delta = np.int64(data_lat_delta); #convert to integer if it's an integer
            elif( data_lat_delta < 1e-4 ):
                data_lat_delta = 1; #override
            #END IF
            data_long_delta = np.median(np.diff(np.unique(data_long)));
            if( np.isclose(data_long_delta,np.int64(data_long_delta)) ):
                data_long_delta = np.int64(data_long_delta); #convert to integer if it's an integer
            elif( data_long_delta < 1e-4 ):
                data_long_delta = 15; #override
            #END IF
            gif_Grid_Lat = np.arange(np.min(plotLatRange),np.max(plotLatRange)+data_lat_delta,data_lat_delta); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
            gif_Grid_Long = np.arange(np.min(plotLongRange),np.max(plotLongRange)+data_long_delta,data_long_delta); #degc, create long points specilized for AMPERE
            pltHelprX, pltHelprY = np.meshgrid( gif_Grid_Long, gif_Grid_Lat); #helps the pcolor work
            gif_Grid = GRITI_movieMaker_subfun_dataGridder(data_lat[k],data_long[k],data_data[k],gif_Grid_Lat,gif_Grid_Long,gif_Grid_Lat.size-1,gif_Grid_Long.size-1,data_lat_delta,data_long_delta,2,101,8).T; #101 disables the data rejection stuff b/c AMPERE doesn't need it
            im = ax.pcolormesh(pltHelprX, pltHelprY, gif_Grid, vmin=np.min(plotLimValu), vmax=np.max(plotLimValu), cmap=settings_dataSpecific['keo colormap'], zorder=100, transform=ccrs.PlateCarree(), rasterized=plot_rasterize); # pseudocolor plot "stretched" to the grid
        else:
            im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['keo scatter size'],c=data_data[k],cmap=settings_dataSpecific['keo colormap'], vmin=np.min(plotLimValu), vmax=np.max(plotLimValu),zorder=100,transform=ccrs.PlateCarree());
        #END IF
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        # cax.yaxis.set_major_formatter(tickr.FormatStrFormatter('%.2f')); #force a rounded format
        # if( np.all(np.mod(cbar.get_ticks(),1) == 0) ):
        #     cax.yaxis.set_major_formatter(tickr.FormatStrFormatter('%.0f')); #force a rounded format
        # else:
        #     cax.yaxis.set_major_formatter(tickr.FormatStrFormatter('%.1f')); #force a rounded format
        # #END IF
        cbar.set_label(settings_dataSpecific['keo labels']+settings_dataSpecific['keo units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
        #END FOR tick
        cbar.mappable.set_clim(vmin=np.min(plotLimValu), vmax=np.max(plotLimValu));
        # cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),5)); #create useful tick marks
        # cax.yaxis.set_ticks(np.linspace(np.min(plotLimValu),np.max(plotLimValu),11)); #create useful tick marks
        cax_ticks = np.arange(np.min(plotLimValu),np.max(plotLimValu)+0.1,0.1); #busted out to deal with 0 being -0
        cax_ticks[np.isclose(cax_ticks, 0)] = 0; #enforce true 0
        cax.yaxis.set_ticks(cax_ticks); #create useful tick marks
        cax.yaxis.set_major_formatter(tickr.FormatStrFormatter('%.1f')); #force a rounded format
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
        
        if( np.all(np.asarray(plotLimValu) >= 0) ):
            #ensure min tick is set
            caxTicks = cbar.get_ticks();
            caxTicks[0] = np.min(plotLimValu); #enforce
            cbar.set_ticks(caxTicks);
        #END IF
    else:
        if( FLG_gridder == True ):
            data_lat_delta = np.median(np.diff(np.unique(data_lat)));
            if( np.isclose(data_lat_delta,np.int64(data_lat_delta)) ):
                data_lat_delta = np.int64(data_lat_delta); #convert to integer if it's an integer
            elif( data_lat_delta < 1e-4 ):
                data_lat_delta = 1; #override
            #END IF
            data_long_delta = np.median(np.diff(np.unique(data_long)));
            if( np.isclose(data_long_delta,np.int64(data_long_delta)) ):
                data_long_delta = np.int64(data_long_delta); #convert to integer if it's an integer
            elif( data_long_delta < 1e-4 ):
                data_long_delta = 15; #override
            #END IF
            gif_Grid_Lat = np.arange(np.min(plotLatRange),np.max(plotLatRange)+data_lat_delta,data_lat_delta); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
            gif_Grid_Long = np.arange(np.min(plotLongRange),np.max(plotLongRange)+data_long_delta,data_long_delta); #degc, create long points specilized for AMPERE
            pltHelprX, pltHelprY = np.meshgrid( gif_Grid_Long, gif_Grid_Lat); #helps the pcolor work
            gif_Grid = GRITI_movieMaker_subfun_dataGridder(data_lat[k],data_long[k],data_data[k],gif_Grid_Lat,gif_Grid_Long,gif_Grid_Lat.size-1,gif_Grid_Long.size-1,data_lat_delta,data_long_delta,2,101,8).T; #101 disables the data rejection stuff b/c AMPERE doesn't need it
            im = ax.pcolormesh(pltHelprX, pltHelprY, gif_Grid, cmap=settings_dataSpecific['keo colormap'],zorder=100, transform=ccrs.PlateCarree(), rasterized=plot_rasterize); # pseudocolor plot "stretched" to the grid
        else:
            im = ax.scatter(data_long[k],data_lat[k],s=settings_dataSpecific['keo scatter size'],c=data_data[k],cmap=settings_dataSpecific['keo colormap'],zorder=100,transform=ccrs.PlateCarree());
        #END IF
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(settings_dataSpecific['keo labels']+settings_dataSpecific['keo units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings_plot['font axis tick FM']); #yee
        #END FOR tick
        cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #END IF
    if( FLG_fancyPlot == 0 ):
        string_Title = 'Keo Angle='+str(np.round(keo_angle,2))+'°, Width='+ \
            str(np.round(keo_width,2))+'arc°, S#='+str(keo_N)+ \
            ', SWidth='+ \
            str(np.round(np.sqrt( (temp_Longs_up[0,0] - temp_Longs_up[0,1])**2 + (temp_Lats_up[0,0] - temp_Lats_up[0,1])**2 ),2))+ \
            'arc° | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' @ '+str(timeHr)+':'+str(timeMin).zfill(2)+':'+str(timeSec).zfill(2); #create mecha title
        if( 'stere' not in settings_mapCOPY['projection name'] ):
            ax.set_title(string_Title,fontproperties=settings_plot['font title FM']); #set the title
        else:
            ax.set_title(string_Title,fontproperties=settings_plot['font title FM'],y=1.075); #set the title
        #END IF
    #END IF
    if( 'stere' not in settings_mapCOPY['projection name'] ):
        if( plot_latLongWords == True ):
            if( coordType == 'geo' ):
                if( (latlong_unitName_bracketed == '') & (latLong_long_dirAdder != '') ):
                   #if no units and just the directions, drop the latitude/longitude
                   ax.set_xlabel(latLong_long_dirAdder.replace('\n',''),fontproperties=settings_plot['font axis label FM']);
                   ax.set_ylabel(latLong_lat_dirAdder.replace('\n',''),fontproperties=settings_plot['font axis label FM']);
                else:
                    ax.set_xlabel('Longitude '+latlong_unitName_bracketed+latLong_long_dirAdder,fontproperties=settings_plot['font axis label FM']);
                    ax.set_ylabel('Latitude '+latlong_unitName_bracketed+latLong_lat_dirAdder,fontproperties=settings_plot['font axis label FM']);
                #END IF
            elif( coordType == 'mag' ):
                ax.set_xlabel('Longitude (Geomag) '+latlong_unitName_bracketed+latLong_long_dirAdder,fontproperties=settings_plot['font axis label FM']);
                ax.set_ylabel('Latitude (Geomag) '+latlong_unitName_bracketed+latLong_lat_dirAdder,fontproperties=settings_plot['font axis label FM']);
            #END IF
        #END IF
    #END IF
    
    #draw a box where the data will be gathered
    if( isinstance(settings_dataSpecific['keo colormap'],str) ):
        keogramLine_color = 'xkcd:fuchsia'; #stands out
        keogramPt_color = settings_map['site marker color']; #default
    else:
        if( np.any(np.all(np.abs(settings_dataSpecific['keo colormap'].colors - np.array([237,17,217])/255) < 0.15, axis=1)) ):
            keogramLine_color = 'xkcd:vermillion'; #if colormap uses fuschia-like color use red instead b/c looks different to colorblind
            keogramPt_color = 'xkcd:dark cyan'; #not more purple
        else:
            keogramLine_color = 'xkcd:fuchsia'; #stands out
            keogramPt_color = settings_map['site marker color']; #default
        #END IF
    #END IF
    lineNum = 200; # number of pts to split line into, if using non-mercator more may be needed to make it look less bad on curvy stuff
    if( ('stere' in settings_mapCOPY['projection name']) & (np.abs(360-np.abs((keo_range[1,1]-keo_range[3,1]))) < 0.1) ):
        temp_mapCoords = [ 
            np.linspace(keo_range[1,1],keo_range[3,1],lineNum) , \
            np.linspace(keo_range[2,0],keo_range[0,0],lineNum) ]; #convert to the geographic map coords
        kj = np.where(np.abs(np.diff(temp_mapCoords[0])) > 355)[0]+1; #get outliers (+1 for diff)
        temp_mapCoords[0] = np.insert(temp_mapCoords[0], kj, np.nan); #nan the gap so it doesn't get weird
        temp_mapCoords[1] = np.insert(temp_mapCoords[1], kj, np.nan); #nan the gap so it doesn't get weird
    else:
        temp_mapCoords = [ np.hstack( [np.linspace(keo_range[0,1],keo_range[1,1],lineNum) , \
            np.linspace(keo_range[1,1],keo_range[3,1],lineNum) , \
            np.linspace(keo_range[3,1],keo_range[2,1],lineNum) , \
            np.linspace(keo_range[2,1],keo_range[0,1],lineNum)] ) , \
            np.hstack( [np.linspace(keo_range[0,0],keo_range[1,0],lineNum) , \
            np.linspace(keo_range[1,0],keo_range[3,0],lineNum) , \
            np.linspace(keo_range[3,0],keo_range[2,0],lineNum) , \
            np.linspace(keo_range[2,0],keo_range[0,0],lineNum)] ) ]; #convert to the geographic map coords
        if( np.abs(360-np.abs((keo_range[1,1]-keo_range[3,1]))) < 0.1 ):
            kj = np.where(np.abs(np.diff(temp_mapCoords[0])) > 355)[0]+1; #get outliers (+1 for diff)
            temp_mapCoords[0] = np.insert(temp_mapCoords[0], kj, np.nan); #nan the gap so it doesn't get weird
            temp_mapCoords[1] = np.insert(temp_mapCoords[1], kj, np.nan); #nan the gap so it doesn't get weird
        #END IF
    #END IF
    ax.plot( temp_mapCoords[0],  #X longitude arcdeg
        temp_mapCoords[1],  #Y latitude arcdeg
        c=keogramLine_color,linewidth=PLOT_lineWidthThicc, zorder=190, transform=ccrs.PlateCarree()); #fuchsia
    
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
            c=keogramLine_color,linewidth=PLOT_lineWidthSmol, zorder=190, transform=ccrs.PlateCarree()); #fuchsia
    #END FOR i
    if( ('stere' in settings_mapCOPY['projection name']) & (np.abs(360-np.abs((keo_range[1,1]-keo_range[3,1]))) < 0.1) ):
        #360 fix again, draws line on -180/180 spot since for non-sphere they don't touch so it's NBD, but on sphere they do touch so it looks weird w/o a line
        temp_mapCoords = ( np.linspace( 180,180,lineNum ) , \
            np.linspace( temp_Lat_List[0,0],temp_Lat_List[0,3],lineNum ) ); #convert to the geographic map coords
        
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c=keogramLine_color,linewidth=PLOT_lineWidthSmol, zorder=190, transform=ccrs.PlateCarree()); #fuchsia
    #END IF
    
    #plot a * where the ISR is
    if( (settings_map['coord type'] == 'geo') & ('keo coord type' in settings_dataSpecific) ): #pull a switcheroo if needed
        if( settings_dataSpecific['keo coord type'] == 'mag' ):
            from Code.subfun_convertToMag import convert_to_mag
            avgPt_coords = copy.deepcopy(avgPt_coords); #copy big time
            [avgPt_coords[0,0], avgPt_coords[0,1]] = convert_to_mag(avgPt_coords[0,0], avgPt_coords[0,1], keo_alt, time4mag); #convert
        #END IF
    #END IF
    if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
        temp_mapCoords = (avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=settings_map['site marker type'], color=keogramPt_color, markersize=settings_map['site marker size'], zorder=150, transform=ccrs.PlateCarree());
    #END IF
    
    if( 'stere' in settings_mapCOPY['projection name'] ):
        from Code.subfun_sunAlsoRises_location import sunAlsoRises_location
        import aacgmv2
        sunSubSolar_loc = sunAlsoRises_location(np.expand_dims(np.array( (dateRange_dayNum_zeroHr[0], np.int16(data_timeUnique[timeIndex]/86400)) ).astype(np.int16),axis=0),timeIndexes=data_timeUnique[timeIndex],timeZone='UTC');
        sunSubSolar_loc['lat'] = sunSubSolar_loc['lat'][0]; #undo for 1 value
        sunSubSolar_loc['long'] = sunSubSolar_loc['long'][0];
        
        if( settings_mapCOPY['coord type'] == 'mag' ):
            [sunSubSolar_loc['lat'], sunSubSolar_loc['long'], _] = aacgmv2.convert_latlon(sunSubSolar_loc['lat'], sunSubSolar_loc['long'], keo_alt, time4mag, method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
        #END IF

        FLG_sunSpin = 0; #disabled for now, cartopy can't spin
        if( FLG_sunSpin == 0):
            x = (1.05*0.5*np.cos((sunSubSolar_loc['long']-90)*np.pi/180))+0.5; #geoMap coordinate 
            y = (1.05*0.5*np.sin((sunSubSolar_loc['long']-90)*np.pi/180))+0.5; #geoMap coordinate 
            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=FONT_axisTickFM);
            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
            ax.add_artist(circleSun); #plot the sun
        else:
            #this does not work in cartopy!! so don't
            ax.set_theta_offset((sunSubSolar_loc['long']-90)*np.pi/180); #turn the whole plot so top is where the sun is
            #I'm not sure if I can get this working easily - not a lot of optons.
        #END IF
    #END IF
    
    if( 'stere' not in settings_mapCOPY['projection name'] ):
        if( FLG_fancyPlot == 0 ):
            tickNumGoal_x = 28; #for x axis ticks
            tickNumGoal_y = 17; #for the y axis ticks
        else:
            tickNumGoal_x = 22; #for x axis ticks
            tickNumGoal_y = 13; #for the y axis ticks
        #END IF
        GRITI_plotHelper_axisizerLatLong(plotLatRange,ax=ax,axDir='y',tickNumGoal=tickNumGoal_y,tickReducer=0,FLG_extendLims=False,FLG_labelUnits=plot_latLongUnitsOnLabels); #auto tick the thing
        GRITI_plotHelper_axisizerLatLong(plotLongRange,ax=ax,axDir='x',tickNumGoal=tickNumGoal_x,tickReducer=3.5,FLG_extendLims=False,FLG_labelUnits=plot_latLongUnitsOnLabels);
    #END IF
    
    figFitter(fig); #fit the fig fast
    if( FLG_fancyPlot == 0 ):
        # plt.show(); #req to make plot show up
        pass; #it's not
    else:
        fig.savefig(os.path.join(settings_paths['fancyPlots'],settings_dataSpecific['keo data type']+'_keo_area'+settings_plot['save file type'])); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF
