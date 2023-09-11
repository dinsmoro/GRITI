"""
GOAL: Plot only ISR POPL HP RTI
RD on 6/4/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tick
import timezonefinder
from datetime import datetime, timedelta
import pytz
# from astroplan import Observer
# import astropy.units as astroUnits
# from astropy.time import Time
from urllib.request import urlopen
from Code.subfun_strstr import strstr
from Code.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
from Code.subfun_figFitter import figFitter
from Code.subfun_figLetteringFitter import subfun_figLetteringFitter
from Code.subfun_sunAlsoRises import sunAlsoRises
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.GRITI_plotHelper_axisizerLatLong import GRITI_plotHelper_axisizerLatLong

def GRITI_TEC_keo_fancyPlot_TEC_wDayNite(vTECChunked_anyAngleAvg,TEC_timeUnique,TEC_plotLimValu, \
        colorMap,plotLatRange,plotLongRange,latMillstone,longMillstone,dateRange_dayNum_zeroHr,time_Ref,latLong_ref, \
        avg_anyAngle,avg_anyAngle_Width,avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name, \
        avg_anyAngle_dataType,avg_anyAngle_plotLabel,dateRange_zeroHr,dateRange_zeroHr_monthName,\
        dateRange_zeroHr_dayPostfix, dateRange_dayNum_full, dateRange_full,\
        FONT_grandioseFM, FONT_titleFM,FONT_axisTick,FONT_axisLabelFM, PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi, 
        settings_plot, settings_dataSpecific):
    
    print('MAKING FANCY PLOT: TEC_avgAnyAngle_plot_wDayNite IN fancyPlot FOLDER'); #report since you won't see anything
    
    letteringPositionX = -0.065; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.88; #set the X position of the lettering (e.g., a. b. c. ...)
    # letteringPositionX = -0.089; #more for bigger font, special use
    
    xAxisLims = ((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #get those xAxis limits, ensures everything is aligned
    
    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth['double plus']; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth['regular plus']; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth['regular']; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth['smol']; #get the line widths
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(TEC_plotLimValu) == 1 ):
        TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    #END IF
    
    #prep the plot
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig = plt.figure(figsize=(14,10.5),dpi=journal_dpi);
    # gridr = gridspec.GridSpec(nrows=7, ncols=1, figure=fig);
    # fig.subplots_adjust(hspace=0.75)
    gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
    gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=20, ncols=1, subplot_spec = gridr[0], hspace = 5.0)
    gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
    # gridr.update(hspace=0.05); # set the spacing between axes.
    fig.add_subplot(gridr1[:]); #RTI plots are 2 tall
    # gridr.update(hspace=0.80); # set the spacing between axes.
    fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
    ax = fig.axes; #get a list of the axes
    
    #Remove the aspect ratio so it fills the screen better
    ax[0].set_aspect('auto'); #set to auto for all axes
    ax[1].set_aspect('auto'); #set to auto for all axes
    
    #!!!-----Plot ISR POPL HP results as a RTI-----!!!
    #Plot just the ISR POPL HP results
    #ZENITH STUFF
    divider0 = make_axes_locatable(ax[0]); #prep to add an axis
    cax0 = divider0.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    pltHelprX, pltHelprY = np.meshgrid( (np.append(TEC_timeUnique,TEC_timeUnique[-1]+np.median(np.diff(TEC_timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im0 = ax[0].pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar0 = fig.colorbar(im0, cax=cax0, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax0.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
    cbar0.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar0.ax.tick_params(labelsize=FONT_axisTick);
    cbar0.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cax0.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax0.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    # string_Title = 'Zenith POPL with High-pass Filter '+str(filter_cutoffPeriod)+' hr Cutoff'; #create mecha title
    # ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    # ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    if( settings_dataSpecific['day nite only local'] == False ):
        ax[0].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    #END IF
    ax[0].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    # #makes plotting look better
    # xAxisLims = np.array(xAxisLims); #get the current x axis limits
    # if( (np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
    #     xAxisLims[0] = np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    # #END IF
    # if( (np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
    #     xAxisLims[1] = np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    # #END IF
    # yAxisLims = np.array(ax[0].get_ylim()); #get the current y axis limits
    # if( (np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)) - np.min(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
    #     yAxisLims[0] = np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    # #END IF
    # if( (np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.max(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
    #     yAxisLims[1] = np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    # #END IF
    # ax[0].set_ylim(yAxisLims); #set the ylims no
    
    # time_autoTick = (np.ceil((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.floor((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    # if( time_autoTick > 24 ):
    #     time_autoTick = 24; #sets the tick setting to 15 arcdegrees per tick
    # elif( time_autoTick > 12 ):
    #     time_autoTick = 12; #sets the tick setting to 15 arcdegrees per tick
    # elif( time_autoTick > 6 ):
    #     time_autoTick = 6; #sets the tick setting to 10 arcdegrees per tick
    # elif( time_autoTick > 4 ):
    #     time_autoTick = 4; #sets the tick setting to 5 arcdegrees per tick
    # elif( time_autoTick > 2 ):
    #     time_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    # elif( time_autoTick >= 1 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    #     time_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    # else:
    #     time_autoTick = ((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600 - (np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600)/13; #just goes for it if it's a super tiny range
    # #END IF
    # xAxisTicks = np.arange( (np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) , \
    #         (np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) + time_autoTick , \
    #         time_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    # ax[0].set_xticks(xAxisTicks); #set x axis ticks
    # # ax[0].set_xticklabels([]); #if statement to remove x axis labels except for the last line
    # ax[0].set_xlim( xAxisLims ); #set x axis limits
    
    # avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    # if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 20 ):
    #     avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    # elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 ):
    #     avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    # elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 ):
    #     avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    # elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 ):
    #     avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    # elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 ):
    #     avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    # elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    #     avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    # else:
    #     if(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Latitude'): #if Y axis is latitude, use latitude
    #         avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/17; #just goes for it if it's a super tiny range
    #     elif(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude'): #if Y axis is longitude, use longitude
    #         avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/17; #just goes for it if it's a super tiny range
    #     #END IF
    # #END IF
    # yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))+avg_anyAngle_Range_Chunks_Long_Plot_autoTick,avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    # ax[0].set_yticks(yAxisTicks); #set x axis ticks
    # if( (np.abs(np.max(yAxisTicks)-np.max(yAxisLims)) <= avg_anyAngle_Range_Chunks_Long_Plot_autoTick/2) & (np.max(yAxisTicks) > np.max(yAxisLims)) ):
    #     yAxisLims[1] = np.max(yAxisTicks); #set the yAxisLims to be the tick if its worth making the tick happen
    # #END IF
    # if( (np.abs(np.min(yAxisTicks)-np.min(yAxisLims)) <= avg_anyAngle_Range_Chunks_Long_Plot_autoTick/2) & (np.min(yAxisTicks) < np.min(yAxisLims)) ):
    #     yAxisLims[0] = np.min(yAxisTicks); #set the yAxisLims to be the tick if its worth making the tick happen
    # #END IF
    # ax[0].set_ylim(yAxisLims); #set the ylims no
    
    xAxisTicks, xAxisLims = GRITI_plotHelper_axisizerTime((time_Ref-dateRange_dayNum_zeroHr[1]*86400)/3600,ax=ax[0]); #automagic time ticks here
    # if( FLG_fancyPlot == 0 ): #different aspect ratios require different spacing assumptions
    #     GRITI_plotHelper_axisizerLatLong(keo_plotLatLong_chunks,ax=ax,axDir='y',tickNumGoal=17);
    # else:
    GRITI_plotHelper_axisizerLatLong(avg_anyAngle_Range_Chunks_Long_Plot,ax=ax[0],axDir='y',tickNumGoal=13);
    #END IF
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
            ax[0].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
            ax[0].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    #END IF   
    if( np.any(np.abs(np.array(plotLongRange)) >= 100) & (avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude') ): #if true, longitude
        letteringA = ax[0].text( letteringPositionX-0.01, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text saying the day or nite
    else:
        letteringA = ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text saying the day or nite
    #END IF
    #~~~~Special arrows for only world-wide longitude plot~~~~
    if( (plotLongRange[0] == -180) and (plotLongRange[1] == 180) and (plotLatRange[0] == -90) and (plotLatRange[1] == 90) and \
            (avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude')  ):
        xlimsDiff = np.diff(ax[0].get_xlim()).item(); #get the x lims
        ylimsDiff = np.diff(ax[0].get_ylim()).item(); #get the y lims
        arrowXoffset = xlimsDiff*0.045; #goal is to get a 45 degree angle arrow despite the disparate scales
        arrowYoffset = ylimsDiff*0.045; #goal is to get a 45 degree angle arrow despite the disparate scales.185
        
        #-----1st arrow set-----
        arrowX = -1.0; #set the X position
        arrowY = -60; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        arrowX = -1.0; #set the X position
        arrowY = -150; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY-arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
            
        #-----2nd arrow set----- 
        arrowX = 2.75; #set the X position
        arrowY = 155; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        arrowX = 2.75; #set the X position
        arrowY = 70; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY-arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        
        # arrowX = 5.5; #set the X position
        # arrowY = -40; #set the Y position
        # ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
        #     arrowprops=dict(width=PLOT_lineWidthDoublePlus,
        #     shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        
        #-----3rd arrow set----- 
        arrowX = 10.6; #set the X position
        arrowY = 45; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        arrowX = 10.6; #set the X position
        arrowY = -20; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY-arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        
        #-----4th arrow set----- 
        arrowX = 13.5; #set the X position
        arrowY = -69; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        arrowX = 13.5; #set the X position
        arrowY = -135; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY-arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        
        #-----5th arrow set----- 
        arrowX = 33; #set the X position
        arrowY = 38; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        arrowX = 33; #set the X position
        arrowY = -20; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY-arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        #-----5th arrow set----- 
        arrowX = 40; #set the X position
        arrowY = -60; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY+arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
        arrowX = 40; #set the X position
        arrowY = -130; #set the Y position
        ax[0].annotate(s='', xy=(arrowX,arrowY), xytext=(arrowX+arrowXoffset,arrowY-arrowYoffset), \
            arrowprops=dict(width=PLOT_lineWidthDoublePlus,
            shrinkA=0,shrinkB=0, color='xkcd:bright purple') ); #use annotate because direct arrow call was doing real dump stuff
    #END IF
    
    
    #!!!DAY NITE PLOT STUFF!!!
    # Xaxisvar_min = (time_Ref[0] - dateRange_dayNum_zeroHr[1]*86400)/3600; #UT hr, min time to compare to
    # Xaxisvar_max = (time_Ref[-1] - dateRange_dayNum_zeroHr[1]*86400)/3600; #UT hr, max time to compare to
#         #below is alt. full time
#         Xaxisvar_min = (min(timeUnique) - dateRange_zeroHr(2))*24; #UT hr, min time to compare to
#         Xaxisvar_max = (max(timeUnique) - dateRange_zeroHr(2))*24; #UT hr, min time to compare to

    #FIRST BATTLE: REGION WHERE LAT/LONG IS   
    if( (np.min(plotLatRange) <= latLong_ref[0][0]) & (np.max(plotLatRange) >= latLong_ref[0][0]) & \
       (np.min(plotLongRange) <= latLong_ref[0][1]) & (np.max(plotLongRange) >= latLong_ref[0][1]) ):
        # Only use latLong_ref if its within the plot area, otherwise ditch it b/c it's set wronk
        latLong_use = (latLong_ref[0][0],latLong_ref[0][1]);
    else:
        print('WARNING: latLong_ref[0] '+str(np.round(latLong_ref[0][0],2)).rstrip('0').rstrip('.')+' lat | '+str(np.round(latLong_ref[0][1],2)).rstrip('0').rstrip('.')+\
              ' long is not within lat range of '+str(np.min(plotLatRange))+' to '+str(np.max(plotLatRange))+\
              ' | long range of '+str(np.min(plotLongRange))+' to '+str(np.max(plotLongRange))+'. Ignoring latLong_ref[0] and using mean of plotLat/LongRange.'+\
              '('+str(np.round(np.mean(plotLatRange),2)).rstrip('0').rstrip('.')+' lat | '+str(np.round(np.mean(plotLongRange),2)).rstrip('0').rstrip('.')+' long)'); # Report a warning on a wrong setting
        latLong_use = (np.mean(plotLatRange),np.mean(plotLongRange)); # Set to mean of plotLat/LongRange
    #END IF    
    tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
    dayNite_timeZoneID = tf.certain_timezone_at(lat=latLong_use[0], lng=latLong_use[1]); #use it to find the time zone
    if dayNite_timeZoneID is None:
        #use geonames site as a backup
        url = 'http://api.geonames.org/timezone?lat='+str(latLong_use[0])+'&lng='+str(latLong_use[1])+'&username=razzluhdzuul'; #create link for lat/long
        webpage = urlopen(url).read(); #get the raw HTML and read it
        try:
            charset = webpage.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
            if( charset is None ):
                charset = 'utf-8'; #assume utf-8
            #END IF
        except:
            charset = 'utf-8'; #assume utf-8
        #END TRY
        webpage = webpage.decode(charset); #"decode" the HTML content so it's legible
        # index_start = strstr(webpage,'<dstOffset>')[0]; #get where dstOffset is
        # index_end = strstr(webpage,'</dstOffset>')[0]; #get where dstOffset is
        # dayNite_DSToffset_str = webpage[index_start+11:index_end]; #get dst offset
        # dayNite_DSToffset = np.float64(dayNite_DSToffset_str); #and convert to number
        # index_start = strstr(webpage,'<gmtOffset>')[0]; #get where UT offset is
        # index_end = strstr(webpage,'</gmtOffset>')[0]; #get where UT offset is
        # dayNite_UToffset_str = webpage[index_start+11:index_end]; #get UT offset
        # dayNite_UToffset = np.float64(dayNite_UToffset_str); #and convert to number
        index_start = strstr(webpage,'<timezoneId>')[0]; #get where time zone ID is
        index_end = strstr(webpage,'</timezoneId>')[0]; #get where time zone ID is
        dayNite_timeZoneID = webpage[index_start+12:index_end]; #get time zone ID
    #END IF
    
    # dateRange_dates = sFUN_dayNumber_to_Date_MULTIPLE(dateRange); #convert day# range to yr/mon/day range with all inbetween filled in
    # dateRange_extended = sFUN_dayNumber_to_Date_MULTIPLE(dateRange,1); #convert day# range to yr/day# flushed out
    # dayNite_DSTactive = isdst(datetime(dateRange_dates[1,1],dateRange_dates[1,2],dateRange_dates[:,3],'TimeZone',dayNite_timeZoneID)); #get DST active or not for each day
    
    #SECOND STEP: CALC LOCAL SUNRISE/SUNSET TIMES
    #will adjust for DST later - calcs done in UT/GMT
    # based on calc steps in https://www.mathworks.com/examples/matlab/community/21093-estimating-sunrise-and-sunset
    # long_corrected = 4*(latLong_ref[1] - 15*dayNite_UToffset); #calc corrected longitude, for sunrise/sunset time

    # dayNite_B = 360*(dateRange_dayNum_full[:,1] - 81)/365; #some sort of angle based on days and stuff
    # dayNite_EoT_corrected = 9.87*np.sin(2*dayNite_B*np.pi/180) - 7.53*np.cos(dayNite_B*np.pi/180) - 1.5*np.sin(dayNite_B*np.pi/180); #eq for Time Correction
    # dayNite_solar_corrected = long_corrected + dayNite_EoT_corrected; #min, solar time correction - for noon

    # dayNite_solar_declination = np.arcsin(np.sin(23.45*np.pi/180)*np.sin(360*(dateRange_dayNum_full[:,1] - 81)/365)*np.pi/180); #deg, solar declination

    # dayNite_sunrise = 12 - np.arccos(-np.tan(latLong_ref[0]*np.pi/180)*np.tan(dayNite_solar_declination*np.pi/180))/15 - dayNite_solar_corrected/60; #hr, sunrise time
    # dayNite_sunset = 12 + np.arccos(-np.tan(latLong_ref[0]*np.pi/180)*np.tan(dayNite_solar_declination*np.pi/180))/15 - dayNite_solar_corrected/60; #hr, sunrise time

    dayNite_sunrise, dayNite_sunset, dateRange_fullPad = sunAlsoRises(dateRange_full,latLong_use[0],latLong_use[1]); # Get the sunrise and sunset times
    dateRange_dayNum_fullPad = subfun_date_to_dayNum(dateRange_fullPad); #convert
    timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
    # timeZoneObj_UTC = pytz.timezone('UTC'); #make a timezone object
    #Daylight savings local time fix
    
    #Pad dateRange_full to account for the timezone offset (a day/night change might happen the next local day, so pad for it)
    # tzString = datetime.strptime(str(dateRange_full[0,:]), '[%Y\t%m\t%d]').astimezone(timeZoneObj).strftime('%z').rstrip('0'); #get the timezone string
    # if( tzString[0] == '+' ):
    #     tempTime = datetime.strptime(str(dateRange_full[-1,:]), '[%Y\t%m\t%d]') + timedelta(days=1);
    #     dateRange_fullPad = np.vstack(( dateRange_full, np.array( ( np.int16(tempTime.strftime('%Y')), np.int16(tempTime.strftime('%m')), np.int16(tempTime.strftime('%d'))) ) ));
    #     dateRange_dayNum_fullPad = np.vstack(( dateRange_dayNum_full, np.array( ( np.int16(tempTime.strftime('%Y')), np.int16(tempTime.strftime('%j'))) ) ));
    # else:
    #     tempTime = datetime.strptime(str(dateRange_full[0,:]), '[%Y\t%m\t%d]') - timedelta(days=1);
    #     dateRange_fullPad = np.vstack(( np.array( ( np.int16(tempTime.strftime('%Y')), np.int16(tempTime.strftime('%m')), np.int16(tempTime.strftime('%d'))) ), dateRange_full ));
    #     dateRange_dayNum_fullPad = np.vstack(( np.array( ( np.int16(tempTime.strftime('%Y')), np.int16(tempTime.strftime('%j'))) ), dateRange_dayNum_full ));
    # #END IF
    
    # for i in range(0,dateRange_fullPad.shape[0]):
    #     dateRange_timeObj = Time(str(dateRange_fullPad[i,0])+'-'+str(dateRange_fullPad[i,1])+'-'+str(dateRange_fullPad[i,2])+'T12:00:00', format='isot', scale='utc'); #make an astropy time object
        
    #     #Run through the days
    #     dayNite_sunriseTemp = ISR.sun_rise_time(dateRange_timeObj, which='nearest', horizon=0*astroUnits.deg, n_grid_points=150).to_value('isot'); #get the sunrise time
    #     dayNite_sunriseObj = timeZoneObj_UTC.localize(datetime.strptime(dayNite_sunriseTemp, '%Y-%m-%dT%H:%M:%S.%f')).astimezone(timeZoneObj); #create datetime object, set it to the UTC time zone (which it is, datetime just doesn't know), then convert it to local time zone
    #     #convert to decimal hour
    #     dayNite_sunriseTemp = dayNite_sunriseObj.isoformat(); #convert to a string that's easily readable
    #     index_start = strstr(dayNite_sunriseTemp,'T')[0]; #get where the hour starts
    #     index_colon = strstr(dayNite_sunriseTemp,':'); #get where the colons are
    #     index_dash = strstr(dayNite_sunriseTemp,'-')[-1]; #get where the last -, end of time
    #     if( index_colon[1] > index_dash ):
    #         #then it is +UTC for local instead of -UTC - so search for last +
    #         index_dash = strstr(dayNite_sunriseTemp,'+')[-1]; #get where the last -, end of time
    #     #END IF
    #     dayNite_sunrise[i] = dateRange_dayNum_fullPad[i,1] + np.int64(dayNite_sunriseTemp[index_start+1:index_colon[0]])/24 + \
    #         np.int64(dayNite_sunriseTemp[index_colon[0]+1:index_colon[1]])/1440 + np.float64(dayNite_sunriseTemp[index_colon[1]+1:index_dash])/86400; #make into day units
        
    #     #Run through the days
    #     dayNite_sunsetTemp = ISR.sun_set_time(dateRange_timeObj, which='nearest', horizon=0*astroUnits.deg, n_grid_points=150).to_value('isot'); #get the sunrise time
    #     dayNite_sunsetObj = timeZoneObj_UTC.localize(datetime.strptime(dayNite_sunsetTemp, '%Y-%m-%dT%H:%M:%S.%f')).astimezone(timeZoneObj); #create datetime object, set it to the UTC time zone (which it is, datetime just doesn't know), then convert it to local time zone
    #     #convert to decimal hour
    #     dayNite_sunsetTemp = dayNite_sunsetObj.isoformat(); #convert to a string that's easily readable
    #     index_start = strstr(dayNite_sunsetTemp,'T')[0]; #get where the hour starts
    #     index_colon = strstr(dayNite_sunsetTemp,':'); #get where the colons are
    #     index_dash = strstr(dayNite_sunsetTemp,'-')[-1]; #get where the last -, end of time
    #     if( index_colon[1] > index_dash ):
    #         #then it is +UTC for local instead of -UTC - so search for last +
    #         index_dash = strstr(dayNite_sunsetTemp,'+')[-1]; #get where the last -, end of time
    #     #END IF
    #     dayNite_sunset[i] = dateRange_dayNum_fullPad[i,1] + np.int64(dayNite_sunsetTemp[index_start+1:index_colon[0]])/24 + \
    #         np.int64(dayNite_sunsetTemp[index_colon[0]+1:index_colon[1]])/1440 + np.float64(dayNite_sunsetTemp[index_colon[1]+1:index_dash])/86400; #make into day units
    # #END FOR i
    timeZoneObj_zeroHr = timeZoneObj.localize(datetime.strptime(str(dateRange_zeroHr), '[%Y\t%m\t%d]')); #time zone info at zero hr
    # dayNite_DSToffset_str = timeZoneObj_zeroHr.strftime('%z'); #get the DST time offset
    # if( strstr(dayNite_DSToffset_str,':').size > 0 ):
    #     if( dayNite_DSToffset_str[strstr(dayNite_DSToffset_str,':')[0]+1:] == '00' ):
    #         dayNite_DSToffset_str = dayNite_DSToffset_str[:strstr(dayNite_DSToffset_str,':')[0]]; #remove the :00 if it's just that
    #         if( dayNite_DSToffset_str[1] == '0' ):
    dayNite_DSToffset_str = timeZoneObj_zeroHr.strftime('%z').replace('0',''); #remove the 0 that was extraneous
            #END IF
    dayNite_DSToffset = np.int64(dayNite_DSToffset_str); #get the number version
        #END IF
        # else:
        #     dayNite_DSToffset = np.float64(dayNite_DSToffset_str[:strstr(dayNite_DSToffset_str,':')[0]]) + \
        #         np.int64(dayNite_DSToffset_str[strstr(dayNite_DSToffset_str,':')[0]+1:])/60; #convert to an hour decimal
        # #END IF
    #END IF
    dayNite_timeZoneName = timeZoneObj_zeroHr.tzname(); #get the time zone name (like 'EST' or 'EDT' depending on standard or daylight savings time)
    dateNite_DSTnUTCOffset = timeZoneObj_zeroHr.dst().total_seconds()/3600; #get the time offset
    if( np.mod(dateNite_DSTnUTCOffset,1) == 0 ):
        dateNite_DSTnUTCOffset = np.int64(dateNite_DSTnUTCOffset); #convert to integer
    #END IF

    #THIRD STEP: PREP FOR PLOTTING BY ALIGNING TIMES, MAKING PLOT VARIABLES
    Xaxisvar_min = np.min(xAxisLims) + dayNite_DSToffset; #hr local, UT time of -12 conv. to local
    Xaxisvar_max = np.max(xAxisLims) + dayNite_DSToffset; #hr local, UT time of -12 conv. to local

    dayNite_sunrise = (dayNite_sunrise + dayNite_DSToffset/24 + dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hr, convert so 0 hr is in the middle (and convert from days to hours)
    dayNite_sunset = (dayNite_sunset + dayNite_DSToffset/24 + dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hr, convert so 0 hr is in the middle (and convert from days to hours)

    xTime = np.sort( np.concatenate( (dayNite_sunrise,dayNite_sunset) ) ); #hr, xvar to plot against
    kmin = xTime < Xaxisvar_min; # Calc these for dealing with dayNite continuation
    kmax = xTime > Xaxisvar_max;
    xTime[ kmin ] = Xaxisvar_min; #limit the xTime to the plotted times
    xTime[ kmax ] = Xaxisvar_max;
    yDayNite = np.ones(xTime.shape); #prep if day or night, set all to day
    for i in range(0,dayNite_sunset.size):
        yDayNite[xTime == dayNite_sunset[i]] = 0; #set sunset times to sunset
    #END FOR i
    yDayNite[kmin] = np.abs(yDayNite[np.where(kmin)[0][-1]+1]-1); #set time below the minimum plotted value to be the other thing
    yDayNite[kmax] = yDayNite[np.where(kmax)[0][0]-1]; #set time above the maximum plotted value to be the same (so there's no weird day shift on the edge of the plot due to stuff)
    xTime = xTime.repeat(2); #interleave repeated values
    yDayNite = np.roll(yDayNite.repeat(2),1) #interleave repeated values and circular shift by 1
    yDayNite[0] = yDayNite[1]; #set that to match (for plotting niceness)
    yDayNite[-1] = yDayNite[-2]; #set that to match (for plotting niceness)
    # if( (xTime[0] > Xaxisvar_min) ):
    #     #if this is so, gotta append the Xaxisvar_min
    #     if( np.any(np.isin(dayNite_sunrise,xTime[0])) ):
    #         #if true, then edge time is a sunrise - so the appended time would be night until that sunrise time
    #         yDayNite[0] = 0; #set to sunset
    #         yDayNite = np.insert(yDayNite,0,np.zeros(2)); #add a new value that's also 0 to make plotting work for the new xTime value
    #     else:
    #         #if false is sunset - so appended time would be day until that sunset time
    #         yDayNite[0] = 1; #set to sunrise
    #         yDayNite = np.insert(yDayNite,0,np.ones(2)); #add a new value that's also 0 to make plotting work for the new xTime value
    #     #END IF
    #     xTime = np.insert(xTime,0,np.tile(Xaxisvar_min,2)); #deal with time edges
    # #END IF
    # if( (xTime[-1] < Xaxisvar_max) ):
    #     #if this is so, gotta append the Xaxisvar_min
    #     if( np.any(np.isin(dayNite_sunset,xTime[-1])) ):
    #         #if true, then edge time is a sunset - so it is sunset until the end of the plotting time
    #         yDayNite[-1] = 0; #set to sunset
    #         yDayNite = np.append(yDayNite,np.zeros(2)); #add a new value that's also 0 to make plotting work for the new xTime value
    #     else:
    #         #if false is sunrise - so it is sunrise until the end of the plotting time
    #         yDayNite[-1] = 1; #set to sunrise
    #         yDayNite = np.append(yDayNite,np.ones(2)); #add a new value that's also 0 to make plotting work for the new xTime value
    #     #END IF
    #     xTime = np.append(xTime,np.tile(Xaxisvar_max,2)); #deal with time edges
    # #END IF
        
    

    # xTime[0] = xTime[1]; #set that to match (for plotting niceness)
    # xTime[-1] = xTime[-2]; #set that to match (for plotting niceness)
    
    #FOURTH STEP: ACTUALLY PLOTTING
    divider1 = make_axes_locatable(ax[1]); #prep to add an axis
    cax1 = divider1.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    cax1.set_visible(False); #mkae it invisible so it matches the other plots in width
    ax[1].plot(xTime,yDayNite,color='xkcd:black',linewidth=PLOT_lineWidthThicc, antialiased=True);
    if( dateNite_DSTnUTCOffset != 0 ):
        strang = dayNite_DSToffset_str+' (Daylight Savings) '+dayNite_timeZoneName+' Time Zone';
    else:
         strang = dayNite_DSToffset_str+' '+dayNite_timeZoneName+' Time Zone';
    #END IF
    # if( yDayNite[-1] == yDayNite[-2] ):
    #     #hacky offset so day or night isn't printed off the plot (or skipped)
    #     offset = 3;
    # else:
    #     offset = 1;
    # #END IF
    for i in range(1,xTime.size-2,2):
        if(xTime[i+1]-xTime[i] > 5 ):
            if( yDayNite[i] == 1 ):
                ax[1].text( (xTime[i+1]-xTime[i])/2+xTime[i], \
                    0.25,'Day', color='k', horizontalalignment='center', fontproperties=FONT_titleFM); #print the text saying the day or nite
            else:
                ax[1].text( (xTime[i+1]-xTime[i])/2+xTime[i], \
                   0.25, 'Night', color='k', horizontalalignment='center', fontproperties=FONT_titleFM); #print the text saying the day or nite
            #END IF
        #END IF
    #END FOR i
    if( settings_dataSpecific['day nite only local'] == False ):
        ax[1].set_xlabel('Local Time [hr] | '+strang,fontproperties=FONT_axisLabelFM);
    else:
        ax[0].set_xlabel('Local Time [hr] | '+dayNite_timeZoneName+' = UT '+dayNite_DSToffset_str+' hr | 0 hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' (Day '+str(dateRange_dayNum_zeroHr[1])+') '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        ax[1].set_xticklabels([]); #remove y labels for day nite plot
        plt.draw(); #req so that labels are current
        labelz = ax[0].get_xticklabels(); #get them labels
        tickz = xAxisTicks+dayNite_DSToffset; #prepare new tick labels (not in str form)
        for i in range(0,len(labelz)):
            labelz[i].set_text(str(tickz[i]).rstrip('0').rstrip('.')); #set the new tick labels
        #END FOR i
        ax[0].set_xticklabels(labelz); #set the new labelz (may not be necessary, seems labelz is already linked but w/e)
    #END IF
    ax[1].set_xticks(xAxisTicks+dayNite_DSToffset); #set x axis ticks
    # ax[1].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1]*86400)/3600) ); #set x axis limits
    ax[1].set_yticklabels([]); #remove y labels for day nite plot
    ax[1].set_yticks([]); #remove y tick marks for day nite plot
    ax[1].set_xlim( (xAxisLims[0]+dayNite_DSToffset,xAxisLims[1]+dayNite_DSToffset) ); #set x axis limits
    ax[1].set_ylim( (0,1) ); #set y lims to 0 and 1
    ax[1].spines['left'].set_visible(False); #turn off box lines
    ax[1].spines['right'].set_visible(False); #turn off box lines
    ax[1].spines['top'].set_visible(False); #turn off box lines
    # ax[2].grid(b=True, which='major', axis='x', color='xkcd:light grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
    letteringB = ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=FONT_grandioseFM, transform=ax[1].transAxes); #print the text saying the day or nite
    
    figFitter(fig); #fit that fig fast
    # if( np.any(np.abs(np.array(plotLongRange)) >= 100) & (avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude') ): #if true, longitude
    #     if( np.any(avg_anyAngle_Range_Chunks_Long_Plot < 0) ):
    #         fig.subplots_adjust(left = 0.098, right = 0.915, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    #     else:
    #         fig.subplots_adjust(left = 0.092, right = 0.915, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    #     #END IF
    # elif( np.any(np.array(plotLatRange) < 0) & (avg_anyAngle_Range_Chunks_Long_Plot_Name != 'Longitude') ):
    #     fig.subplots_adjust(left = 0.080, right = 0.915, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
        
    # elif( np.any(avg_anyAngle_Range_Chunks_Long_Plot < 0) & (avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude') ):
    #     fig.subplots_adjust(left = 0.080, right = 0.915, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    # else:
    #     fig.subplots_adjust(left = 0.065, right = 0.915, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    # #END IF
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    
    #-----Attempt to detect a./b. overlap with negative signs-----
    subfun_figLetteringFitter(fig, ax); #call a function that makes the a./b./.. lettering fit automagically
    
    fig.savefig(folder[3]+'\TEC_avgKeo_wDayNite'+settings_plot['save file type']); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff