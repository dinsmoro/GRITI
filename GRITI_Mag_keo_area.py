"""
GOAL: Plot area that average alg with average over
RD on 6/11/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
import cartopy as cartopy #cartopy replaces basemap b/c it's updated
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER #use these to make it like say N or S next to numbers on the plot ticks
from subfun_figFitter import figFitter

def GRITI_Mag_keo_area(data, dates, settings, keo_width, keo_range, temp_lats_up, temp_longs_up, temp_lat_list, temp_long_list):
    
    def strike(text): #from https://stackoverflow.com/questions/25244454/python-create-strikethrough-strikeout-overstrike-string-type @pdw
        result = ''
        for c in text:
            result = result + c + '\u0336'
        return result
    #END DEF
    
    #==============UNPACK SETTINGS==============
    keo_angle = settings['Mag']['keo angle orig'];
    # keo_width = settings['Mag']['keo width orig'];
    keo_N = settings['Mag']['keo N'];
    # keo_45vsLatLong = settings['Mag']['keo 45 lat/long'];
    map_plotLatRange = settings['Mag']['lat range'];
    map_plotLongRange = settings['Mag']['long range'];
    
    if( settings['Mag']['keo set stations'] != 1 ):
        siteNames = data['Mag']['site names']; #get the site names
    else:
        siteNames = settings['Mag']['keo set stations names']; #set the site names to whatever the user had
    #END IF
    removeList = []; #prep
    for j in range(0,len(siteNames)):
        if( ~((data['Mag'][siteNames[j]]['lat'] >= np.min(settings['Mag']['lat range'])) & \
                (data['Mag'][siteNames[j]]['lat'] <= np.max(settings['Mag']['lat range'])) & \
                (data['Mag'][siteNames[j]]['long'] >= np.min(settings['Mag']['long range'])) & \
                (data['Mag'][siteNames[j]]['long'] <= np.max(settings['Mag']['long range']))) ):
            removeList.append(siteNames[j]); #prep for removal b/c not within the range we want
        #END IF
    #END FOR j
    #solving lists takes a lot of slowww lists ohw ell
    for j in range(0,len(removeList)):
        siteNames.remove(removeList[j]); #remove the stuff we don't need
    #END FOR j
    siteNamesIndex = np.where(np.in1d(data['Mag']['site names'],siteNames))[0]; #get indexes [keeps site coloring consistent even if all are not in use]
    
    #plot help with autotick calculating
    map_long_autoTick = (np.max(map_plotLongRange) - np.min(map_plotLongRange))/25; #tries to split the longitude range into 25 parts (based off of 360/15+1)
    if( map_long_autoTick > 14 ):
        map_long_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    elif( map_long_autoTick > 10 ):
        map_long_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( map_long_autoTick > 5 ):
        map_long_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( map_long_autoTick > 2 ):
        map_long_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( map_long_autoTick > 1 ):
        map_long_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( map_long_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        map_long_autoTick = 1; #sets the tick setting to 1 arcdegree per tick
    else:
        map_long_autoTick = (np.max(map_plotLongRange) - np.min(map_plotLongRange))/15; #just goes for it if it's a super tiny range
    #END IF
    map_lat_autoTick = (np.max(map_plotLatRange) - np.min(map_plotLatRange))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( map_lat_autoTick > 10 ):
        map_lat_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( map_lat_autoTick > 5 ):
        map_lat_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( map_lat_autoTick > 2 ):
        map_lat_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( map_lat_autoTick > 1 ):
        map_lat_autoTick = 2; #sets the tick setting to 2 arcdegrees per tick
    elif( map_lat_autoTick > 0.75 ): #0.75 because 10/13 = 0.76something and it sounded good for enough 1 arcdeg ticks
        map_lat_autoTick = 1; #sets the tick setting to 1 arcdegree per tick
    else:
        map_lat_autoTick = (np.max(map_plotLatRange) - np.min(map_plotLatRange))/13; #just goes for it if it's a super tiny range
    #END IF

    #---START TO PLOT HERE---
    fig = plt.figure(); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    ax = fig.add_subplot(1,1,1, projection=settings['map']['projection']);
    # divider = make_axes_locatable(ax[0]); #prep to add an axis
    # dividerKeo = make_axes_locatable(ax[1]); #prep to add an axis
    ax.set_aspect('auto');
    
    #---ADD GRID LINES, SET PLOTTING AREA---
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True); #draw some well-described gridlines
    gl.xlabels_top = False; #turn off all, let ticks be handled by set_xticks
    gl.xlabels_bottom = False; #turn off all, let ticks be handled by set_xticks
    gl.ylabels_right = False; #turn off all, let ticks be handled by set_yticks
    gl.ylabels_left = False; #turn off all, let ticks be handled by set_yticks
    gl.xlocator = mticker.FixedLocator(np.arange(np.min(map_plotLongRange),np.max(map_plotLongRange)+map_long_autoTick,map_long_autoTick)); #this works ok, but be consistent use set_xticks
    gl.ylocator = mticker.FixedLocator(np.arange(np.min(map_plotLatRange),np.max(map_plotLatRange)+map_lat_autoTick,map_lat_autoTick)); #this doesn't plot -90 and 90 labels, but is req to get the gridlines right
    ax.set_xticks(np.arange(np.min(map_plotLongRange),np.max(map_plotLongRange)+map_long_autoTick,map_long_autoTick),crs=settings['map']['projection']); #gotta plot ticks with this to get -90 and 90
    ax.set_yticks(np.arange(np.min(map_plotLatRange),np.max(map_plotLatRange)+map_lat_autoTick,map_lat_autoTick),crs=settings['map']['projection']); #gotta plot ticks with this to get -90 and 90
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    ax.set_extent(map_plotLongRange + map_plotLatRange); #set the plot extent, set at end - x and y ticks can extend plot area so this'll reign it in
    
    #---DRAW SOME COASTLINES, MAYBE COLOR IN SOME STUFF---
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get info on the size of the plot area to know what geographic scale to use
    mapper_resolution = np.max( [np.abs(map_plotLatRange[0]-map_plotLatRange[1])/bbox.height , np.abs(map_plotLongRange[0]-map_plotLongRange[1])/bbox.width] ); #degc, max extent covered in the plot
    if( mapper_resolution > 20 ): #arbitrary numbers
        mapper_resolution = '110m'; #the resolution to use for plotting geographical features
    elif( mapper_resolution > 0.5 ): #arbitrary numbers
        mapper_resolution = '50m'; #the resolution to use for plotting geographical features
    else:
        #otherwise if the deg/in for the plot is super small use the highest detail possible
        mapper_resolution = '10m'; #the resolution to use for plotting geographical features
    #END IF
    if( settings['map']['world color'] == True ):
        ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', mapper_resolution, edgecolor='face', facecolor=settings['map']['land color'], alpha=0.75)); #idk what these calls really mean
        ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'ocean', mapper_resolution, edgecolor='face', facecolor=settings['map']['water color'], alpha=0.75)); #idk what these calls really mean
    #END IF
    ax.coastlines(resolution=mapper_resolution, color='xkcd:black'); #draw the coastlines
    #reinforce axis limits after this drawing stuff
    ax.set_xlim(map_plotLongRange); #set x limits
    ax.set_ylim(map_plotLatRange); #set y limits
    
    
    #---ACTUALLY PLOT REAL STUFF HERE---
    sN = []; #prep empty list of text handles
    mN = []; #prep empty list of text handles
    for j in range(0,len(siteNames)):
        mN.append( ax.plot( data['Mag'][siteNames[j]]['long'], data['Mag'][siteNames[j]]['lat'], color=settings['plot']['color'][siteNamesIndex[j]], marker=settings['map']['marker type'], markersize=settings['map']['marker size'], linewidth=0, zorder=50, transform=settings['map']['projection'] ) ); #plot the sites
        sN.append( ax.text( data['Mag'][siteNames[j]]['long'], data['Mag'][siteNames[j]]['lat'], siteNames[j]+'\n  '+str(siteNamesIndex[j])+'', zorder=55, transform=settings['map']['projection']) ); #write name of site
    #END FOR j
    if( settings['Mag']['keo set stations'] == 1 ):
        sitesRedacted = np.setxor1d(data['Mag']['site names'],settings['Mag']['keo set stations names']);
        siteNamesIndex = np.where(np.in1d(data['Mag']['site names'],sitesRedacted))[0]; #get indexes [keeps site coloring consistent even if all are not in use]
        sNR = []; #prep empty list of text handles
        mNR = []; #prep empty list of text handles
        for j in range(0,len(sitesRedacted)):
            mNR.append( ax.plot( data['Mag'][sitesRedacted[j]]['long'], data['Mag'][sitesRedacted[j]]['lat'], color=settings['plot']['color'][siteNamesIndex[j]], marker=settings['map']['marker type'], markersize=settings['map']['marker size'], linewidth=0, zorder=50, transform=settings['map']['projection'] ) ); #plot the sites
            sNR.append( ax.text( data['Mag'][sitesRedacted[j]]['long'], data['Mag'][sitesRedacted[j]]['lat'], strike(sitesRedacted[j])+'\n  '+str(siteNamesIndex[j])+'', zorder=55, transform=settings['map']['projection'], c='xkcd:grey') ); #write name of site
        #END FOR j
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
        c='xkcd:fuchsia',linewidth=settings['plot']['line width']['thicc'], 
        zorder=85, transform=settings['map']['projection']);
    
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
        temp_mapCoords = ( np.linspace( temp_long_list[i,0],temp_long_list[i,3],200 ) , \
            np.linspace( temp_lat_list[i,0],temp_lat_list[i,3],200 ) ); #convert to the geographic map coords
        
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c='xkcd:fuchsia',linewidth=settings['plot']['line width']['smol'],
            zorder=85, transform=settings['map']['projection']);
    #END FOR i
    
    string_Title = 'NR Canada Magnetic Observatories'+' Keo w/ Angle='+str(np.round(keo_angle,2))+'°, Keo Width='+ \
        str(np.round(keo_width,2))+' arcdeg, Step #='+str(keo_N)+ \
        ', Step Width='+ \
        str(np.round(np.sqrt( (temp_longs_up[0,0] - temp_longs_up[0,1])**2 + (temp_lats_up[0,0] - temp_lats_up[0,1])**2 ),2))+ \
        ' arcdeg'; #create mecha title
    ax.set_title(string_Title,fontproperties=settings['plot']['font title FM']); #set the title
    
    #---FINALIZE PLOTTING---
    figFitter(fig); #fit that fig fast
    # fig.subplots_adjust(left = 0.045, right = 0.980, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
