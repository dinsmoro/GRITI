"""
GOAL: Plot area that average alg with average over
RD on 2/16/21

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tick
import os
import matplotlib.patches as pat
from subfun_figFitter import figFitter

def GRITI_TEC_keo_fancyPlot_area_explainer(TEC_time,TEC_dTEC,\
        TEC_lat,TEC_long,TEC_timeUnique, TEC_plotLimValu, time_Ref, \
        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, avgPt_coords, avg_anyAngle_Range, \
        avg_anyAngle, avg_anyAngle_Width, avg_anyAngle_N, avg_anyAngle_dataType, avg_anyAngle_plotLabel, \
        plotLatRange,plotLongRange,plotLatRange_autoTick,plotLongRange_autoTick_Crunched, \
        dateRange_dayNum_zeroHr, gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size, \
        colorMap, avg_anyAngle_polarMode,geoMap_projectionStyle,BasemapFixDir, \
        FONT_titleFM, FONT_axisLabelFM, FONT_axisTickFM, FONT_axisTick, \
        PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi):
    print('MAKING FANCY PLOT: TEC_keo_plot_area_explainer IN fancyPlot FOLDER'); #report since you won't see anything
    
    #-----Start plotting-----
    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth['double plus']; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth['regular plus']; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth['regular']; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth['smol']; #get the line widths
    
    FLG_zoomer = 1; #user gotta set this, zooms in around the * or doesn't
    if( FLG_zoomer == 0 ):
        gif_Millstone_Marker_Size = gif_Millstone_Marker_Size;
        TEC_scatterSize = 20;
    else:
        gif_Millstone_Marker_Size = 30;
        TEC_scatterSize = 80;
    #END IF
    
    if( FLG_zoomer == 1 ):
        #----- Zoom in on plot area -----
        if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
            #--- Do a mini ray trace to find where to zoom in on ---
            isIn = False; #set to start
            i = 0; #prep cntr
            while( (isIn == False) & (i < avg_anyAngle_N) ):
                #--- Step 1, Calc Ref Vectors ---
                V12 = np.array( (temp_Long_List[i,1]-temp_Long_List[i,0],temp_Lat_List[i,1]-temp_Lat_List[i,0]) ); #get ref vect
                M12 = np.sqrt(V12[0]**2 + V12[1]**2); #get the mag
                V14 = np.array( (temp_Long_List[i,3]-temp_Long_List[i,0],temp_Lat_List[i,3]-temp_Lat_List[i,0]) ); #get ref vect
                M14 = np.sqrt(V14[0]**2 + V14[1]**2); #get the mag
                V32 = np.array( (temp_Long_List[i,1]-temp_Long_List[i,2],temp_Lat_List[i,1]-temp_Lat_List[i,2]) ); #get ref vect
                M32 = np.sqrt(V32[0]**2 + V32[1]**2); #get the mag
                V34 = np.array( (temp_Long_List[i,3]-temp_Long_List[i,2],temp_Lat_List[i,3]-temp_Lat_List[i,2]) ); #get ref vect
                M34 = np.sqrt(V34[0]**2 + V34[1]**2); #get the mag
                #--- Step 2, Calc Ref "Angles", don't take the acos --- (these should always turn out to be 90 in this instance)
                cos1214 = np.dot(V12,V14)/(M12*M14); #get the angle between the two vectors, don't take the acos for speed
                cos3234 = np.dot(V32,V34)/(M32*M34); #get the angle between the two vectors, don't take the acos for speed
                #--- Step 3, Calc Arb. Point Vectors ---
                V1p = np.array( (avgPt_coords[0,1]-temp_Long_List[i,0],avgPt_coords[0,0]-temp_Lat_List[i,0]) ); #get ref vect
                M1p =np.sqrt(V1p[0]**2 + V1p[1]**2); #get the mag
                V3p = np.array( (avgPt_coords[0,1]-temp_Long_List[i,2],avgPt_coords[0,0]-temp_Lat_List[i,2]) ); #get ref vect
                M3p = np.sqrt(V3p[0]**2 + V3p[1]**2); #get the mag
                #--- Step 4, Calc Arb. Point "Angles", don't take the acos ---
                cos121p = np.dot(V12,V1p)/(M12*M1p); #get the angle between the two vectors, don't take the acos for speed
                cos141p = np.dot(V14,V1p)/(M14*M1p); #get the angle between the two vectors, don't take the acos for speed
                cos323p = np.dot(V32,V3p)/(M32*M3p); #get the angle between the two vectors, don't take the acos for speed
                cos343p = np.dot(V34,V3p)/(M34*M3p); #get the angle between the two vectors, don't take the acos for speed
                #--- Step 5, Check if Arb. Point is inside the rectangle ---
                isIn = (cos121p >= cos1214) & (cos141p >= cos1214) & \
                    (cos323p >= cos3234) & (cos343p >= cos3234)
                
                i += 1; #increment cntr
            #END WHILE
            isIn = i-1; #get the place where it was in
            plotLatRange = [np.min(temp_Lat_List[isIn,:]), np.max(temp_Lat_List[isIn,:])]; #arcdeg, overwrite the plot limits
            plotLongRange = [np.min(temp_Long_List[isIn,:]), np.max(temp_Long_List[isIn,:])]; #arcdeg, overwrite the plot limits
        else:
            print('ERROR: Plot area of '+str(plotLatRange)+' lat '+str(plotLongRange)+' long does not include the defined ISR point at '+str(avgPt_coords[0,:])+ 'lat, long. Plot needs that. Quitting via crash.');
            import sys #just to fake crash
            sys.crash();
        #END IF
    #END IF
    
    plotLatRange_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( plotLatRange_autoTick > 10 ):
        plotLatRange_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( plotLatRange_autoTick > 5 ):
        plotLatRange_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( plotLatRange_autoTick > 2 ):
        plotLatRange_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( plotLatRange_autoTick > 1 ):
        plotLatRange_autoTick = 2; #sets the tick setting to 2 arcdegrees per tick
    elif( plotLatRange_autoTick > 0.75 ): #0.75 because 10/13 = 0.76something and it sounded good for enough 1 arcdeg ticks
        plotLatRange_autoTick = 1; #sets the tick setting to 1 arcdegree per tick
    else:
        plotLatRange_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/14; #just goes for it if it's a super tiny range
    #END IF
    plotLongRange_autoTick_Crunched = (np.max(plotLongRange) - np.min(plotLongRange))/13; #tries to split the longitude range into 25 parts (based off of 360/15+1)
    if( plotLongRange_autoTick_Crunched > 25 ):
        plotLongRange_autoTick_Crunched = 30; #sets the tick setting to 15 arcdegrees per tick
    elif( plotLongRange_autoTick_Crunched > 10 ):
        plotLongRange_autoTick_Crunched = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( plotLongRange_autoTick_Crunched > 5 ):
        plotLongRange_autoTick_Crunched = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( plotLongRange_autoTick_Crunched > 2 ):
        plotLongRange_autoTick_Crunched = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( plotLongRange_autoTick_Crunched > 1 ):
        plotLongRange_autoTick_Crunched = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( plotLongRange_autoTick_Crunched >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        plotLongRange_autoTick_Crunched = 1; #sets the tick setting to 1 arcdegree per tick
    else:
        plotLongRange_autoTick_Crunched = (np.max(plotLongRange) - np.min(plotLongRange))/7; #just goes for it if it's a super tiny range
    #END IF
    
    
    os.environ["PROJ_LIB"] = BasemapFixDir; #hack beacuse people are awful coders
    from mpl_toolkits.basemap import Basemap #import here because need that fix
    
    #THIS IS TO VIEW DATA AVG RANGE BEING TAKEN
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
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
        geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(plotLongRange)),np.ceil(np.max(plotLongRange))+plotLongRange_autoTick_Crunched,plotLongRange_autoTick_Crunched),2), 
            labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='w' ); #adds the labels but keeps the lines invisible
        geoMap.drawparallels(np.round(np.arange(np.floor(np.min(plotLatRange)),np.ceil(np.max(plotLatRange))+plotLatRange_autoTick,plotLatRange_autoTick),2), 
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
    
    #Fix basemap ticks
    ticks_lat = np.round(np.arange(np.floor(np.min(plotLatRange)),np.ceil(np.max(plotLatRange))+plotLatRange_autoTick,plotLatRange_autoTick),2);
    ticks_lat = np.delete(ticks_lat,np.where((ticks_lat > np.max(plotLatRange))|(ticks_lat < np.min(plotLatRange))) ); #remove out-of-boudns
    ticks_long = np.round(np.arange(np.floor(np.min(plotLongRange)),np.ceil(np.max(plotLongRange))+plotLongRange_autoTick_Crunched,plotLongRange_autoTick_Crunched),2);
    ticks_long = np.delete(ticks_long,np.where((ticks_long > np.max(plotLongRange))|(ticks_long < np.min(plotLongRange))) ); #remove out-of-boudns
    
    ticks_long, _ = geoMap(ticks_long, np.zeros(ticks_long.shape)); #convert the lat/long arcdeg to the current map coordinates
    _, ticks_lat = geoMap(np.zeros(ticks_lat.shape), ticks_lat); #convert the lat/long arcdeg to the current map coordinates
    
    ax.set_xticks(ticks_long); #set the tick mark locations
    ax.set_yticks(ticks_lat); #set the tick mark locations
    ax.tick_params(axis='both',which='major'); #turn on the ticks
    ax.xaxis.set_ticklabels([]); #remove the labels (Basemap already has it labeled)
    ax.yaxis.set_ticklabels([]); #remove the labels (Basemap already has it labeled)
    
    # (_, counts) = np.unique(TEC_float[:,TEC_time], return_counts=True); #get the counts of the data
    # timeIndex = np.where( counts >= TEC_float[:,TEC_time].size/TEC_timeUnique.size )[0][0]; #get an index where there's a lot of data
    # timeHr = np.int16((TEC_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600); #calc the hour
    # timeMin = np.int16(np.abs((TEC_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)*60); #calc the minute
    # timeSec = np.int16(np.round((np.abs((TEC_timeUnique[timeIndex]-dateRange_dayNum_zeroHr[1]*86400)/3600 - timeHr)*60 - timeMin)*60)); #calc the second
    timeIndex = np.where( np.abs(time_Ref[0] - TEC_timeUnique) == np.min(np.abs(time_Ref[0] - TEC_timeUnique)) )[0][0]; #get an index where there's a lot of data
    
    k = np.where( TEC_time == TEC_timeUnique[timeIndex])[0]; #gets during a time period
    TEC_latLongMapped = geoMap(TEC_long[k],TEC_lat[k]); #convert the lat/long arcdeg to the current map coordinates
    
    im = ax.scatter(TEC_latLongMapped[0],TEC_latLongMapped[1],s=TEC_scatterSize,c=TEC_dTEC[k],cmap=colorMap, vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.f')); #force a rounded format
    cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar.ax.tick_params(labelsize=FONT_axisTick) 
    cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu))
    cax.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax.yaxis.label.set_font_properties(FONT_axisLabelFM)
    # string_Title = avg_anyAngle_dataType+' Avging w/ Angle='+str(np.round(avg_anyAngle,2))+' deg, Avg Width='+ \
    #     str(np.round(avg_anyAngle_Width,2))+' arcdeg, Step #='+str(avg_anyAngle_N)+ \
    #     ', Step Width='+ \
    #     str(np.round(np.sqrt( (temp_Longs_up[0,0] - temp_Longs_up[0,1])**2 + (temp_Lats_up[0,0] - temp_Lats_up[0,1])**2 ),2))+ \
    #     ' arcdeg | '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' at '+str(timeHr)+':'+str(timeMin)+':'+str(timeSec); #create mecha title
    # ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
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
        c='xkcd:fuchsia',linewidth=PLOT_lineWidthThicc);
    
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
            c='xkcd:fuchsia',linewidth=PLOT_lineWidthSmol);
    #END FOR i
    
    #---for methodsx paper to show where the points are--
    for i in range(0,avg_anyAngle_Range.shape[0]):
        temp_mapCoords = geoMap(avg_anyAngle_Range[i,1],avg_anyAngle_Range[i,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker='o', color='xkcd:black', markersize=10, zorder=50,clip_on=False);
    #END FOR i
    for i in np.arange(temp_arangeStart,avg_anyAngle_N,temp_arangeSpacing):
        temp_mapCoords = geoMap( np.array( (temp_Long_List[i,0],temp_Long_List[i,3]) ) , \
            np.array( (temp_Lat_List[i,0],temp_Lat_List[i,3]) ) ); #convert to the geographic map coords
        ax.plot(temp_mapCoords[0][0],temp_mapCoords[1][0],marker='X', color='xkcd:black', markersize=13, zorder=50,clip_on=False);
        ax.plot(temp_mapCoords[0][1],temp_mapCoords[1][1],marker='X', color='xkcd:black', markersize=13, zorder=50,clip_on=False);
    #END FOR i
    if( FLG_zoomer == 1 ):
        plot_xDist = np.diff(plotLongRange).item(); #arcdeg, get the long delta
        plot_yDist = np.diff(plotLatRange).item(); #arcdeg, get the lat delta
        temp_mapCoords = geoMap(temp_Long_List[isIn,0]+plot_xDist*.0195,temp_Lat_List[isIn,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], '$P_1$', 
            color='xkcd:black',horizontalalignment='left',verticalalignment='center',fontproperties=FONT_axisLabelFM , clip_on=False);
        temp_mapCoords = geoMap(temp_Long_List[isIn,1],temp_Lat_List[isIn,1]-plot_yDist*.0175); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], r'$P_2$', 
            color='xkcd:black',horizontalalignment='center',verticalalignment='top',fontproperties=FONT_axisLabelFM , clip_on=False);
        temp_mapCoords = geoMap(temp_Long_List[isIn,2]-plot_xDist*.0195,temp_Lat_List[isIn,2]); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], '$P_3$', 
            color='xkcd:black',horizontalalignment='right',verticalalignment='center',fontproperties=FONT_axisLabelFM , clip_on=False);
        temp_mapCoords = geoMap(temp_Long_List[isIn,3],temp_Lat_List[isIn,3]+plot_yDist*.0175); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], r'$P_4$', 
            color='xkcd:black',horizontalalignment='center',verticalalignment='bottom',fontproperties=FONT_axisLabelFM , clip_on=False);
        
    
        temp_mapCoords = geoMap(temp_Long_List[isIn,0]+plot_xDist*.10*np.cos(np.linspace(-avg_anyAngle*np.pi/180,avg_anyAngle*np.pi/180,num=200,endpoint=True)),
            temp_Lat_List[isIn,0]+plot_xDist*.10*np.sin(np.linspace(-avg_anyAngle*np.pi/180,avg_anyAngle*np.pi/180,num=200,endpoint=True))); #convert the lat/long arcdeg to the current map coordinates
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c='xkcd:black',linewidth=PLOT_lineWidthSmol,linestyle=':',clip_on=False,zorder=0);
        temp_mapCoords = geoMap(temp_Long_List[isIn,0]+plot_xDist*.1025,temp_Lat_List[isIn,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], r'$\theta_{12,14}$', 
            color='xkcd:black',horizontalalignment='left',verticalalignment='center',fontproperties=FONT_axisLabelFM , clip_on=False);
        
        temp_mapCoords = geoMap(temp_Long_List[isIn,2]+plot_xDist*.10*np.cos(np.linspace(-avg_anyAngle*np.pi/180-np.pi,avg_anyAngle*np.pi/180-np.pi,num=200,endpoint=True)),
            temp_Lat_List[isIn,2]+plot_xDist*.10*np.sin(np.linspace(-avg_anyAngle*np.pi/180-np.pi,avg_anyAngle*np.pi/180-np.pi,num=200,endpoint=True))); #convert the lat/long arcdeg to the current map coordinates
        ax.plot( temp_mapCoords[0] , #X longitude arcdeg
            temp_mapCoords[1] , #Y latitude arcdeg
            c='xkcd:black',linewidth=PLOT_lineWidthSmol,linestyle=':',clip_on=False,zorder=0);
        temp_mapCoords = geoMap(temp_Long_List[isIn,2]-plot_xDist*.1025,temp_Lat_List[isIn,2]); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], r'$\theta_{32,34}$', 
            color='xkcd:black',horizontalalignment='right',verticalalignment='center',fontproperties=FONT_axisLabelFM , clip_on=False);
        
        temp_mapCoords = geoMap((temp_Long_List[isIn,0]+temp_Long_List[isIn,1])/2-plot_xDist*.007,(temp_Lat_List[isIn,0]+temp_Lat_List[isIn,1])/2); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], '$\overrightarrow{V}_{12}$', 
            color='xkcd:black',horizontalalignment='right',verticalalignment='bottom',fontproperties=FONT_axisLabelFM , clip_on=False);
        
        temp_mapCoords = geoMap((temp_Long_List[isIn,0]+temp_Long_List[isIn,3])/2,(temp_Lat_List[isIn,0]+temp_Lat_List[isIn,3])/2); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], '$\overrightarrow{V}_{14}$', 
            color='xkcd:black',horizontalalignment='right',verticalalignment='top',fontproperties=FONT_axisLabelFM , clip_on=False);
        
        temp_mapCoords = geoMap((temp_Long_List[isIn,2]+temp_Long_List[isIn,1])/2,(temp_Lat_List[isIn,2]+temp_Lat_List[isIn,1])/2); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], '$\overrightarrow{V}_{32}$', 
            color='xkcd:black',horizontalalignment='left',verticalalignment='bottom',fontproperties=FONT_axisLabelFM , clip_on=False);
        
        temp_mapCoords = geoMap((temp_Long_List[isIn,2]+temp_Long_List[isIn,3])/2,(temp_Lat_List[isIn,2]+temp_Lat_List[isIn,3])/2); #convert the lat/long arcdeg to the current map coordinates
        ax.text( temp_mapCoords[0],temp_mapCoords[1], '$\overrightarrow{V}_{34}$', 
            color='xkcd:black',horizontalalignment='left',verticalalignment='top',fontproperties=FONT_axisLabelFM , clip_on=False);
    #END IF
    
    #plot a * where the ISR is
    if( (avgPt_coords[0,0] <= np.max(plotLatRange)) & (avgPt_coords[0,0] >= np.min(plotLatRange)) &(avgPt_coords[0,1] <= np.max(plotLongRange)) & (avgPt_coords[0,1] >= np.min(plotLongRange)) ): #only plot if the * is within the range of plotting
        temp_mapCoords = geoMap(avgPt_coords[0,1],avgPt_coords[0,0]); #convert the lat/long arcdeg to the current map coordinates
        ax.plot(temp_mapCoords[0],temp_mapCoords[1],marker=gif_Millstone_Marker, color=gif_Millstone_Marker_Color, markersize=gif_Millstone_Marker_Size, zorder=50);
    #END IF
    
    figFitter(fig); #fit the fig fast
    # if( np.all(np.array(plotLongRange) < 0) == True ):
    #     fig.subplots_adjust(left = 0.095, right = 0.9178, top = 0.980, bottom = 0.052); #sets padding to small numbers for minimal white space
    # else:
    #     fig.subplots_adjust(left = 0.095, right = 0.9178, top = 0.980, bottom = 0.042); #sets padding to small numbers for minimal white space
    # #END IF
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig(folder[3]+'\\TEC_avgKeo_area_explainer.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff