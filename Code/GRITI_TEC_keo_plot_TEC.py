#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.subfun_textNice import textNice
from Code.subfun_figFitter import figFitter

def GRITI_TEC_keo_plot_TEC(vTECChunked_anyAngleAvg,TEC_timeUnique,TEC_plotLimValu,colorMap,plotLatRange,plotLongRange,latMillstone,longMillstone,dateRange_dayNum_zeroHr,avg_anyAngle,avg_anyAngle_Width,avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name,avg_anyAngle_dataType,avg_anyAngle_plotLabel,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):

    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(TEC_plotLimValu) == 1 ):
        TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    #END IF

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (np.append(TEC_timeUnique,TEC_timeUnique[-1]+np.median(np.diff(TEC_timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    if( np.any(np.isinf(TEC_plotLimValu)) == False ):
        im = ax.pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),5)); #create useful tick marks
        cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        #cbar.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu)); #they changed how the code works, this doesn't work anymore
        cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu)); #now it's this
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    else:
        im = ax.pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,cmap=colorMap); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #END IF
    #    string_Title = 'TEC Averaged on Angle of '+textNice(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #        textNice(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+textNice(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    string_Title = avg_anyAngle_dataType+' Averaged on Angle of '+textNice(np.round(avg_anyAngle,2))+' deg and Width of '+ \
        textNice(np.round(avg_anyAngle_Width,2))+' arcdeg'; #create mecha title
    ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax.set_xlabel('Time in UT - 0 Hr on Day '+textNice(dateRange_dayNum_zeroHr[1])+', '+textNice(dateRange_dayNum_zeroHr[0])+' [hr]',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax.set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    # xAxis_min = np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600);
    # xAxis_max = np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600);
    
    # xAxis_delta = (xAxis_max - xAxis_min)/19; #19 is hard-coded maximum number of tick marks (and labels) on the x axis that fit
    # if( xAxis_delta <= 12 ): #slice it up based on time
    #     if( xAxis_delta < 1 ):
    #         xAxis_delta = 1;            
    #     elif( xAxis_delta < 2 ):
    #         xAxis_delta = 2;  
    #     elif( xAxis_delta < 4 ):
    #         xAxis_delta = 4;  
    #     elif( xAxis_delta < 6 ):
    #         xAxis_delta = 6;  
    #     elif( xAxis_delta <= 12 ):
    #         xAxis_delta = 12;
    #     #END IF
    # else:
    #     xAxis_delta = 24*(xAxis_delta//24); #just do a multiple of 24
    # #END IF
    
    # xAxisTicks = np.arange( (np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),xAxis_delta)) , \
    #         (np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),xAxis_delta)) + xAxis_delta , \
    #         xAxis_delta); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    # ax.set_xticks(xAxisTicks); #set x axis ticks
    # ax.set_xlim((xAxis_min,xAxis_max)); #set x axis limits
    GRITI_plotHelper_axisizerTime((TEC_timeUnique-dateRange_dayNum_zeroHr[1]*86400)/3600,ax=ax); #automagic time ticks here
    
    avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    else:
        if(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Latitude'): #if Y axis is latitude, use latitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/13; #just goes for it if it's a super tiny range
        elif(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude'): #if Y axis is longitude, use longitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/13; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)),avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax.set_yticks(yAxisTicks); #set y axis ticks
    ax.set_ylim((np.min(yAxisTicks),np.max(yAxisTicks))); #set y axis limits
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
            ax.plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=1); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
            ax.plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=1); #plots a point with a black line
        #END IF
    #END IF
    
    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.050, right = 0.945, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
       
#no return, plot is the return!




