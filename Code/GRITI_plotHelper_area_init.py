"""
GOAL: Initiate a Cartopy plot
RD on 6/11/21

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tick
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
import cartopy.feature as cfeature
from Code.subfun_textNice import textNice
from Code.subfun_pltPause import pause as pltPause

def GRITI_plotHelper_area_init(plotLatRange, plotLongRange, settings_map, settings_plot, FLG_fancyPlot, 
                               FLG_wordColorDub=True, numSubplots=False, time4mag=None, alt4mag=120., 
                               figAlreadyInfo=None, FLG_double_cax=False, FLG_latLabel=True, FLG_longLabel=True, FLG_forcePolarLatLabels=False,
                               figsizer_override=None, figPPI_override=None):
    
    if( figAlreadyInfo == None ):
        #Initialize a map plot
        
        #this only matters on ioff fancy plots
        if( np.all(figsizer_override == None) ):
            figsizer = (14,8.5);
        else:
            figsizer = figsizer_override;
        #END IF
        if( np.all(figPPI_override == None) ):
            figPPI = settings_plot['journal dpi'];
        else:
            figPPI = figPPI_override;
        #END IF
        
        if( numSubplots != False ):
            if( np.isscalar(numSubplots) ):
                numRows = numSubplots;
                numCols = 1;
            else:
                numRows = numSubplots[0];
                numCols = numSubplots[1];
            #END IF
        #END IF        
        if( FLG_fancyPlot == 0 ):
            if( numSubplots != False ):
                fig, ax = plt.subplots(nrows=numRows,ncols=numCols, subplot_kw={'projection': settings_map['projection']}); #use instead of fig because it inits an axis too (I think I dunno)
                ax = ax.ravel(); #flatten
            else:
                fig, ax = plt.subplots(subplot_kw={'projection': settings_map['projection']}); #use instead of fig because it inits an axis too (I think I dunno)
                ax = [ax]; #list it
            #END IF
            figManager = fig.canvas.manager; #req to maximize
            figManager.window.showMaximized(); #force maximized
            
            map_lat_autoTick = settings_map['lat autotick']; #get the auto ticks
            map_long_autoTick = settings_map['lat autotick']; #get the auto ticks
        else:
            plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
            if( numSubplots != False ):
                fig, ax = plt.subplots(nrows=numRows,ncols=numCols,figsize=figsizer,dpi=figPPI, subplot_kw={'projection': settings_map['projection']}); #use instead of fig because it inits an axis too (I think I dunno)
                ax = ax.ravel(); #flatten
            else:
                fig, ax = plt.subplots(figsize=figsizer,dpi=figPPI_override, subplot_kw={'projection': settings_map['projection']}); #use instead of fig because it inits an axis too (I think I dunno)
                ax = [ax]; #list it
            #END IF
            map_lat_autoTick = settings_map['lat autotick fancy']; #get the auto ticks
            map_long_autoTick = settings_map['lat autotick fancy']; #get the auto ticks
        #END IF
    else:
        fig = figAlreadyInfo['fig'];
        ax = figAlreadyInfo['ax'];
        if( np.size(ax) != 1 ):
            if( isinstance(ax,list) | isinstance(ax,tuple) ):
                numSubplots = len(ax); #list of axes, get the num
                numRows = numSubplots;
                numCols = 1;
            else:
                numSubplots = ax.shape; #list of axes
                numRows = numSubplots[0];
                numCols = numSubplots[1];
                ax = ax.ravel(); #flatten
            #END IF
        else:
            ax = [ax];
        #END IF
        if( figAlreadyInfo['fancy plot'] == 0 ):
            map_lat_autoTick = settings_map['lat autotick']; #get the auto ticks
            map_long_autoTick = settings_map['lat autotick']; #get the auto ticks
        else:
            map_lat_autoTick = settings_map['lat autotick fancy']; #get the auto ticks
            map_long_autoTick = settings_map['lat autotick fancy']; #get the auto ticks
        #END IF
    #END IF
    
    for i in range(0,len(ax)):
        # ax[i].axis('off'); #cartopy makes a new axis
        # ax[i] = plt.axes(projection=settings_map['projection']); #redefine the axis to be a geographical axis [AFTER cax else there's ERRORS oof]
        ax[i].set_aspect('auto');
            
        #---ADD GRID LINES, SET PLOTTING AREA---
        if( 'stere' not in settings_map['projection name'] ):
            if( figAlreadyInfo == None ):
                if( FLG_double_cax == False ):
                    cax = ax[-1].inset_axes((1.02, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax.yaxis.label.set_text('holder');
                    cax.yaxis.set_label_position('right')
                    cax.yaxis.set_ticks_position('right'); #simulate real deal
                else:
                    cax = [None, None]; #prep a list
                    cax_i = 0;
                    cax[cax_i] = ax[-1].inset_axes((1.02, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax[cax_i].yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax[cax_i].yaxis.label.set_text('holder');
                    cax[cax_i].yaxis.set_label_position('right')
                    cax[cax_i].yaxis.set_ticks_position('right'); #simulate real deal
                    cax_i += 1;
                    cax[cax_i] = ax[-1].inset_axes((-0.02, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax[cax_i].yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax[cax_i].yaxis.label.set_text('holder');
                    cax[cax_i].yaxis.set_label_position('left')
                    cax[cax_i].yaxis.set_ticks_position('left'); #simulate real deal
                #END IF
            #$END IF
            # divider = make_axes_locatable(ax); #prep to add an axis
            # cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
            
            # gl = ax.gridlines(linewidth=PLOT_lineWidthSmoller, color='xkcd:black', alpha=0.5, linestyle='--', draw_labels=True); #draw some well-described gridlines
            # gl.xlabels_top = False; #turn off all, let ticks be handled by set_xticks
            # gl.xlabels_bottom = False; #turn off all, let ticks be handled by set_xticks
            # gl.ylabels_right = False; #turn off all, let ticks be handled by set_yticks
            # gl.ylabels_left = False; #turn off all, let ticks be handled by set_yticks
            # gl.xlocator = tick.FixedLocator(np.arange(np.min(plotLongRange),np.max(plotLongRange)+map_long_autoTick,map_long_autoTick)); #this works ok, but be consistent use set_xticks
            # gl.ylocator = tick.FixedLocator(np.arange(np.min(plotLatRange),np.max(plotLatRange)+map_lat_autoTick,map_lat_autoTick)); #this doesn't plot -90 and 90 labels, but is req to get the gridlines right
            ax[i].set_xticks(np.arange(np.floor(np.min(plotLongRange)),np.ceil(np.max(plotLongRange))+map_long_autoTick,map_long_autoTick),crs=settings_map['projection']); #gotta plot ticks with this to get -90 and 90
            ax[i].set_yticks(np.arange(np.floor(np.min(plotLatRange)),np.ceil(np.max(plotLatRange))+map_lat_autoTick,map_lat_autoTick),crs=settings_map['projection']); #gotta plot ticks with this to get -90 and 90
            # gl.xformatter = LONGITUDE_FORMATTER
            # gl.yformatter = LATITUDE_FORMATTER
            # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
            # ax.set_extent(plotLongRange + plotLatRange); #set the plot extent, set at end - x and y ticks can extend plot area so this'll reign it in
            #--- following based on https://stackoverflow.com/a/61986546/2403531 ---
            #this makes it be the lat range limits but makes it square
        else:
            if( (figAlreadyInfo == None) & FLG_latLabel & FLG_longLabel ):
                if( FLG_double_cax == False ):
                    cax = ax[-1].inset_axes((1.135, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax.yaxis.label.set_text('holder');
                    cax.yaxis.set_label_position('right')
                    cax.yaxis.set_ticks_position('right'); #simulate real deal
                else:
                    cax = [None, None]; #prep a list
                    cax_i = 0;
                    cax[cax_i] = ax[-1].inset_axes((1.145, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax[cax_i].yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax[cax_i].yaxis.label.set_text('holder');
                    cax[cax_i].yaxis.set_label_position('right')
                    cax[cax_i].yaxis.set_ticks_position('right'); #simulate real deal
                    cax_i += 1;
                    cax[cax_i] = ax[-1].inset_axes((-0.19, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax[cax_i].yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax[cax_i].yaxis.label.set_text('holder');
                    cax[cax_i].yaxis.set_label_position('left')
                    cax[cax_i].yaxis.set_ticks_position('left'); #simulate real deal
                #END IF
                FLG_cartopyLabels = True; #use em
            elif( (figAlreadyInfo == None) & (( not FLG_latLabel ) | ( not FLG_longLabel )) ):
                if( FLG_double_cax == False ):
                    cax = ax[-1].inset_axes((1.05, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax.yaxis.label.set_text('holder');
                    cax.yaxis.set_label_position('right')
                    cax.yaxis.set_ticks_position('right'); #simulate real deal
                else:
                    cax = [None, None]; #prep a list
                    cax_i = 0;
                    cax[cax_i] = ax[-1].inset_axes((1.05, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax[cax_i].yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax[cax_i].yaxis.label.set_text('holder');
                    cax[cax_i].yaxis.set_label_position('right')
                    cax[cax_i].yaxis.set_ticks_position('right'); #simulate real deal
                    cax_i += 1;
                    cax[cax_i] = ax[-1].inset_axes((-0.08, 0, 0.02, 1)); #make a color bar axis on last axis
                    cax[cax_i].yaxis.label.set_font_properties(settings_plot['font axis label FM']);
                    cax[cax_i].yaxis.label.set_text('holder');
                    cax[cax_i].yaxis.set_label_position('left')
                    cax[cax_i].yaxis.set_ticks_position('left'); #simulate real deal
                #END IF
                FLG_cartopyLabels = False; #make my own I guess
            elif( (figAlreadyInfo != None) & FLG_latLabel & FLG_longLabel ):
                FLG_cartopyLabels = True; #just go for it
            elif( (figAlreadyInfo != None) & (( not FLG_latLabel ) | ( not FLG_longLabel )) ):
                FLG_cartopyLabels = False; #nop
            #END IF
            # divider = make_axes_locatable(ax); #prep to add an axis
            # cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
            
            from Code.GRITI_plotHelper_area_cartopyStereLabeler import GRITI_plotHelper_area_cartopyStereLabeler
            ax[i].set_extent(plotLongRange + np.flip(plotLatRange).tolist(),ccrs.PlateCarree()); #set the extent, only platecarree makes anything remotely work with this call
            ax[i].set_aspect('equal')
            #from https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html
            import matplotlib.path as mpath
            # Compute a circle in axes coordinates, which we can use as a boundary
            # for the map. We can pan/zoom as much as we like - the boundary will be
            # permanently circular.
            theta = np.linspace(0, 2*np.pi, 100);
            center, radius = [0.5, 0.5], 0.5; #does not support south hemi 1-np.min(plotLatRange)/90
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T;
            circle = mpath.Path(verts * radius + center);
            ax[i].set_boundary(circle, transform=ax[i].transAxes);
            gl = ax[i].gridlines(draw_labels=FLG_cartopyLabels, zorder=375, linewidth=1.0, color='black', alpha=.75, linestyle='--', xlabel_style={'fontproperties':settings_plot['font axis tick FM']}, ylabel_style={'fontproperties':settings_plot['font axis tick FM']}); #draw gridlines was zorder=375
            ticks_long = np.flip(np.arange(180,-180,-30)); #long ticks
            ticks_lat = np.arange(45,90,15); #lat ticks
            gl.xlocator = tick.FixedLocator(ticks_long);
            gl.ylocator = tick.FixedLocator(ticks_lat);
            if( FLG_forcePolarLatLabels == True ):
                # cartopy draw doesn't seem to be able to draw latitude labels yet, so still do it myself
                GRITI_plotHelper_area_cartopyStereLabeler(ax[i],ticks_lat,ticks_long,settings_plot,np.min(plotLatRange), FLG_latLabel=True, FLG_longLabel=False);  
            #END IF
            if( ( not FLG_latLabel ) | ( not FLG_longLabel ) ):
                #only need to do it myself if both false since cartopy added in labels finally to stere
                GRITI_plotHelper_area_cartopyStereLabeler(ax[i],ticks_lat,ticks_long,settings_plot,np.min(plotLatRange), FLG_latLabel=FLG_latLabel, FLG_longLabel=FLG_longLabel);            
            #END IF
        #END IF
    #END FOR i

    
    #---DRAW SOME COASTLINES, MAYBE COLOR IN SOME STUFF---
    fig.canvas.draw(); #key for all instances
    if( plt.isinteractive() == True ): #only needs fig.canvas.draw() for interactive, flush and show and pause aren't req
        fig.canvas.flush_events(); #only needed for interactive plotting
        # plt.show(); #req to make plot show up
        pltPause(0.01); #this is oddly highly required for plotting in an IDE interactively (seems 0.01 is lowest I can go)
    #END IF
    
    for i in range(0,len(ax)):
        bbox = ax[i].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get info on the size of the plot area to know what geographic scale to use
        mapper_resolution = np.max( [np.abs(plotLatRange[0]-plotLatRange[1])/bbox.height , np.abs(plotLongRange[0]-plotLongRange[1])/bbox.width] ); #degc, max extent covered in the plot
        if( mapper_resolution > 25 ): #arbitrary numbers
            mapper_resolution = '110m'; #the resolution to use for plotting geographical features
        elif( mapper_resolution > 0.5 ): #arbitrary numbers
            mapper_resolution = '50m'; #the resolution to use for plotting geographical features
        else:
            #otherwise if the deg/in for the plot is super small use the highest detail possible
            mapper_resolution = '10m'; #the resolution to use for plotting geographical features
        #END IF
        if( (settings_map['world color'] == True) & (FLG_wordColorDub == True) & (settings_map['coord type'] == 'geo') ):
            ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', mapper_resolution, edgecolor='face', facecolor=settings_map['land color'], alpha=0.75,zorder=50)); #idk what these calls really mean
            ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', mapper_resolution, edgecolor='face', facecolor=settings_map['water color'], alpha=0.75,zorder=49)); #idk what these calls really mean
        #END IF
        if( settings_map['coord type'] == 'geo' ):
            ax[i].coastlines(resolution=mapper_resolution, color='xkcd:black',zorder=275); #draw the coastlines
        elif( settings_map['coord type'] == 'mag' ):
            cc = cfeature.NaturalEarthFeature('physical', 'coastline', mapper_resolution, color='xkcd:black',zorder=75); #get a feature to mess with
            geom_mag = []; #prep a list
            for geom in cc.geometries():
                geom_mag.append(convert_to_mag(geom, time4mag, alt4mag));
                # geom_mag.append(geom);
                # coords = list(geom.coords);
            #END FOR geom
            cc_mag = cfeature.ShapelyFeature(geom_mag, settings_map['projection'], color='xkcd:black',zorder=75);
            for geom in cc_mag.geometries():
                try:
                    ax[i].plot(*geom.coords.xy, color='xkcd:black', linewidth=1.0, zorder=275, transform=ccrs.PlateCarree());
                except NotImplementedError:
                    pass; #avoids "NotImplementedError: Multi-part geometries do not provide a coordinate sequence"
                    #I don't really know or care to figure out what's going on in the geom stuff (it is obtuse) so when something errors just skip it I guess - seems only a few things do this
            #END FOR geom
            # ax.add_feature(cc_mag);
        #END IF
        #reinforce axis limits after this drawing stuff
        #reinforce axis limits after this drawing stuff
        if( settings_map['projection name'] == 'mill' ):
            ax[i].set_xlim(plotLongRange); #set x limits
            ax[i].set_ylim(plotLatRange); #set y limits
        #END IF
    #END IF
    
    if( numSubplots == False ):
        ax = ax[0]; #pull out of list
    else:
        if( numCols != 1 ):
            ax = ax.reshape(numRows,numCols); #reform to original shape for later use (ezier this way)
        #END IF
    #END IF
    
    if( figAlreadyInfo == None ): #only need to return if it didn't exist before
        return fig, ax, cax
    else:
        return ax
    #END IF
#END DEF

def convert_to_mag(geom, time4mag, alt4mag): #based on fabulously thorough code at https://gis.stackexchange.com/a/291293
    import aacgmv2 #install with: pip install aacgmv2 [need buildtools what a pain]
    # from math import isnan as isnan
    if geom.is_empty:
        return geom
    #END IF
    if geom.has_z:
        def convert_to_mag_doer(coords, time4mag, _):
            for long, lat, alt4mag in coords:
                [lat_mag, long_mag, alt_mag] = aacgmv2.convert_latlon(lat, long, alt4mag, time4mag, method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
                yield (long_mag, lat_mag, alt_mag)
            #END FOR long, lat, alt_mag
        #END DEF
    else:
        def convert_to_mag_doer(coords, time4mag, alt4mag):
            for long, lat in coords:
                [lat_mag, long_mag, _] = aacgmv2.convert_latlon_arr(lat, long, alt4mag, time4mag, method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
                # tryCntr = 0;
                # while( (isnan(long_mag.item()) | isnan(lat_mag.item())) & (tryCntr < 5) ):
                #     #recalc at higher altitude
                #     [long_mag, lat_mag, _] = aacgmv2.convert_latlon_arr(lat, long, alt4mag+200, time4mag, method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
                #     tryCntr += 1 #increment
                # #END IF
                yield (long_mag.item(), lat_mag.item())
            #END FOR long, lat
        #END DEF
    #END IF
    # Process coordinates from each supported geometry type
    if geom.geom_type in ('Point', 'LineString', 'LinearRing'): #updated from geom.type to geom.geom_type due to depreciation warning (noted just in case older versions lack geom_type)
        converted = list(convert_to_mag_doer(geom.coords, time4mag, alt4mag)); #list of tuples
        converted_numpy = np.asarray(converted); #to numpy for easy peasy
        kj = np.where(np.abs(np.diff(converted_numpy[:,0])) > 355)[0]+1; #get outliers (+1 for diff)
        converted_numpy = np.insert(converted_numpy, kj, np.nan, axis=0); #nan the gap so it doesn't get weird
        converted = converted_numpy.tolist(); #all list now, but should be OK
        return type(geom)(converted)
    elif geom.geom_type == 'Polygon':
        ring = geom.exterior
        shell = type(ring)(list(convert_to_mag_doer(ring.coords, time4mag, alt4mag)));
        holes = list(geom.interiors);
        for pos, ring in enumerate(holes):
            holes[pos] = type(ring)(list(convert_to_mag_doer(ring.coords, time4mag, alt4mag)));
        #END FOR pos, ring
        return type(geom)(shell, holes)
    elif geom.geom_type.startswith('Multi') or geom.type == 'GeometryCollection':
        # Recursive call
        return type(geom)([convert_to_mag(part, time4mag, alt4mag) for part in geom.geoms])
    else:
        raise ValueError('Type %r not recognized' % geom.type)
    #END IF
#END DEF