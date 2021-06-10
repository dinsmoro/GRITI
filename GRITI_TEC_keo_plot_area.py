"""
GOAL: Plot area that average alg with average over
RD on 6/11/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tick
import os
from subfun_figFitter import figFitter

def GRITI_TEC_keo_plot_area(TEC_time,TEC_dTEC,\
        TEC_lat,TEC_long,TEC_timeUnique, TEC_plotLimValu, time_Ref, \
        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, avgPt_coords, avg_anyAngle_Range, \
        avg_anyAngle, avg_anyAngle_Width, avg_anyAngle_N, avg_anyAngle_dataType, avg_anyAngle_plotLabel, \
        plotLatRange,plotLongRange,plotLatRange_autoTick,plotLongRange_autoTick, \
        dateRange_dayNum_zeroHr, gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size, \
        colorMap, avg_anyAngle_polarMode,geoMap_projectionStyle,BasemapFixDir, \
        FONT_titleFM, FONT_axisLabelFM, FONT_axisTickFM, FONT_axisTick):
    #----Start plotting-----
    os.environ["PROJ_LIB"] = BasemapFixDir; #hack beacuse people are awful coders
    from mpl_toolkits.basemap import Basemap #import here because need that fix
    
    #THIS IS TO VIEW DATA AVG RANGE BEING TAKEN
    fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    if( (np.abs(np.diff(plotLatRange)).item() > 30) | (np.abs(np.diff(plotLongRange)).item() > 30) ):
        areaThreshold = 10000; #bigger area thresh prevents drawing lil lakes
    else:
        #if close, should draw small stuff (else things like PR don't get drawn)
        areaThreshold = 1000; #smaller area thresh allows for PR to get drawn
    #END IF
    
    if( avg_anyAngle_polarMode == 0):
        #mill for square Mercator style
        #robin for oval shape
        geoMap = Basemap(projection=geoMap_projectionStyle, lat_0=np.mean(plotLatRange), lon_0=np.mean(plotLongRange), #projection type, and I think lat_0/lon_0 are the centers?
            resolution = 'i', area_thresh = areaThreshold, ax=ax, #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
            llcrnrlon=np.float32(plotLongRange[0]), llcrnrlat=np.float32(plotLatRange[0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
            urcrnrlon=np.float32(plotLongRange[1]), urcrnrlat=np.float32(plotLatRange[1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
                    
        #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
        geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(plotLongRange)),np.ceil(np.max(plotLongRange))+1,plotLongRange_autoTick),2), 
            labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='w' ); #adds the labels but keeps the lines invisible
        geoMap.drawparallels(np.round(np.arange(np.floor(np.min(plotLatRange)),np.ceil(np.max(plotLatRange))+1,plotLatRange_autoTick),2), 
            labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='w' ); #adds the labels but keeps the lines invisible
        #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
        #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax.set_aspect('auto');
        
        # fig.subplots_adjust(left = 0.038, right = 0.945, top = 0.96, bottom = 0.05); #sets padding to small numbers for minimal white space
    else: 
        if( np.mean(plotLatRange) >= 0 ): #north pole projection
            geoMap = Basemap(projection=geoMap_projectionStyle,boundinglat=np.float32(plotLatRange[0]), lat_0=90, lon_0=np.mean(plotLongRange), #projection type, and I think lat_0/lon_0 are the centers?
                resolution = 'i', area_thresh = areaThreshold, ax=ax, round=True); #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
            #the other calls used above don't seem to do anything - boundinglat seems to be the limiter here
            parallelTicks = np.round(np.arange(90,np.min(plotLatRange),-15),0);
        else: #south pole projection
            geoMap = Basemap(projection=geoMap_projectionStyle,boundinglat=np.float32(plotLatRange[0]), lat_0=-90, lon_0=np.mean(plotLongRange), #projection type, and I think lat_0/lon_0 are the centers?
                resolution = 'i', area_thresh = areaThreshold, ax=ax, round=True); #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
            #the other calls used above don't seem to do anything - boundinglat seems to be the limiter here
            parallelTicks = np.round(np.arange(-90,np.max(plotLatRange),15),0);
        #END IF
                        
        #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
        geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(plotLongRange)),np.ceil(np.max(plotLongRange))+1,30),0), 
            labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
        geoMap.drawparallels(parallelTicks, 
            labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
        #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
        #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
        for j in np.arange(0,360,30): #longitude labels
            x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
            y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
            if( j > 180 ):
                angle = j-360; #deg, flip to negative
            else:
                angle = j; #deg, angle is OK
            #END IF
            if( angle == 180):
                y = y - 0.01; #small nudge, cause this one is too close
            #END IF
            if( angle == -60):
                x = x - 0.003; #small nudge, cause this one is too close
            #END IF
            ax.text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=FONT_axisTickFM)
        #END FOR j
        latPts = np.roll(parallelTicks,-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
        latPts[-1] = np.min(plotLatRange); #degc, last is the final latitude value
        latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
        for j in range(0,latPts.size): #latitude labels
            x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
            y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
            ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=FONT_axisTickFM)
            #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=FONT_axisTickFM)
        #END FOR j

        #Remove the aspect ratio from the basemap so it fills the screen better
        ax.set_aspect('equal');
        geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
        geoCircle.set_clip_on(False); #prevent weird things where the circle is clipped off
        
        # fig.subplots_adjust(left = 0.035, right = 0.95, top = 0.96, bottom = 0.03); #sets padding to small numbers for minimal white space
    #END IF
         
    geoMap.drawcoastlines();
    #map.drawcountries()
    #map.fillcontinents(color='coral')
    #map.drawmapboundary()
    
    #xticks( np.round(np.arange(np.floor(np.min(plotLongRange)),np.ceil(np.max(plotLongRange))+1,plotLongRange_autoTick),2) ); #creates x ticks automagically
    #yticks( np.round(np.arange(np.floor(np.min(plotLatRange)),np.ceil(np.max(plotLatRange))+1,plotLatRange_autoTick),2) ); #creates y ticks automagically
    
    # (_, counts) = np.unique(TEC_float[:,TEC_time], return_counts=True); #get the counts of the data
    # timeIndex = np.where( counts >= TEC_float[:,TEC_time].size/TEC_timeUnique.size )[0][0]; #get an index where there's a lot of data
    timeIndex = np.where( np.abs(time_Ref[0] - TEC_timeUnique) == np.min(np.abs(time_Ref[0] - TEC_timeUnique)) )[0][0]; #get an index where there's a lot of data
    timeHr = np.int16((TEC_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600); #calc the hour
    timeMin = np.int16(np.abs((TEC_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)/60); #calc the minute
    timeSec = np.int16(np.round(np.abs((TEC_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)/60 - timeMin)); #calc the second
    if( timeSec == 60 ):
        timeSec = 0; #roll bak to  0
        timeMin += 1; #increment by 1
    #END IF
    if( timeMin == 60 ):
        timeMin = 0; #roll back to 0
        timeHr += 1; #increment by 1
    #END IF
    
    k = np.where( TEC_time == TEC_timeUnique[timeIndex])[0]; #gets during a time period
    TEC_latLongMapped = geoMap(TEC_long[k],TEC_lat[k]); #convert the lat/long arcdeg to the current map coordinates
    
    if( np.any(np.isinf(TEC_plotLimValu)) == False ):
        im = ax.scatter(TEC_latLongMapped[0],TEC_latLongMapped[1],s=20,c=TEC_dTEC[k],cmap=colorMap, vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    else:
        im = ax.scatter(TEC_latLongMapped[0],TEC_latLongMapped[1],s=20,c=TEC_dTEC[k],cmap=colorMap);
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #END IF
    string_Title = avg_anyAngle_dataType+' Avging w/ Angle='+str(np.round(avg_anyAngle,2))+' deg, Avg Width='+ \
        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Step #='+str(avg_anyAngle_N)+ \
        ', Step Width='+ \
        str(np.round(np.sqrt( (temp_Longs_up[0,0] - temp_Longs_up[0,1])**2 + (temp_Lats_up[0,0] - temp_Lats_up[0,1])**2 ),2))+ \
        ' arcdeg | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' at '+str(timeHr)+':'+str(timeMin)+':'+str(timeSec); #create mecha title
    ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    #draw a box where the data will be gathered
    temp_mapCoords = geoMap( np.hstack( [np.linspace(avg_anyAngle_Range[0,1],avg_anyAngle_Range[1,1],200) , \
        np.linspace(avg_anyAngle_Range[1,1],avg_anyAngle_Range[3,1],200) , \
        np.linspace(avg_anyAngle_Range[3,1],avg_anyAngle_Range[2,1],200) , \
        np.linspace(avg_anyAngle_Range[2,1],avg_anyAngle_Range[0,1],200)] ) , \
        np.hstack( [np.linspace(avg_anyAngle_Range[0,0],avg_anyAngle_Range[1,0],200) , \
        np.linspace(avg_anyAngle_Range[1,0],avg_anyAngle_Range[3,0],200) , \
        np.linspace(avg_anyAngle_Range[3,0],avg_anyAngle_Range[2,0],200) , \
        np.linspace(avg_anyAngle_Range[2,0],avg_anyAngle_Range[0,0],200)] ) ); #convert to the geographic map coords
    ax.plot( temp_mapCoords[0],  #X longitude arcdeg
        temp_mapCoords[1],  #Y latitude arcdeg
        c='xkcd:fuchsia',linewidth=4);
    
    if( (np.arange(10,avg_anyAngle_N,10).size > 10) & (np.arange(10,avg_anyAngle_N,10).size <= 30) ):
        #a good zone - not too few, not too many
        temp_arangeStart = 10;
        temp_arangeSpacing = 10;    
    else:
        #otherwise need different start and spacing
        if( avg_anyAngle_N > 10 ):
            temp_arangeStart = np.int64(np.round(avg_anyAngle_N/10)); #scale it
            temp_arangeSpacing = np.int64(np.round(avg_anyAngle_N/10)); #scale it
        else:
            temp_arangeStart = 1; #min it at 1
            temp_arangeSpacing = 1; #min it at 1
        #END IF
    #END IF
            
    for i in np.arange(temp_arangeStart,avg_anyAngle_N,temp_arangeSpacing):
        temp_mapCoords = geoMap( np.linspace( temp_Long_List[i,0],temp_Long_List[i,3],200 ) , \
            np.linspace( temp_Lat_List[i,0],temp_Lat_List[i,3],200 ) ); #convert to the geographic map coords
        
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c='xkcd:fuchsia',linewidth=1);
    #END FOR i
    
    #plot a * where the ISR is
    if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
        temp_mapCoords = geoMap(avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=gif_Millstone_Marker, color=gif_Millstone_Marker_Color, markersize=gif_Millstone_Marker_Size, zorder=50);
    #END IF
    
    figFitter(fig); #fit the fig fast
    plt.show(); #req to make plot show up