#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tick
from subfun_figFitter import figFitter

def GRITI_TEC_keo_fancyPlot_TEC(vTECChunked_anyAngleAvg,TEC_timeUnique,TEC_plotLimValu, \
        colorMap,plotLatRange,plotLongRange,latMillstone,longMillstone,dateRange_dayNum_zeroHr,time_Ref,latLong_ref, \
        avg_anyAngle,avg_anyAngle_Width,avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name, \
        avg_anyAngle_dataType,avg_anyAngle_plotLabel,dateRange_zeroHr,dateRange_zeroHr_monthName,dateRange_zeroHr_dayPostfix,\
        FONT_titleFM,FONT_axisTick,FONT_axisLabelFM, PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi):
    print('MAKING FANCY PLOT: TEC_avgAnyAngle_plot IN fancyPlot FOLDER'); #report since you won't see anything
    
    #----Start plotting-----
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
    
    xAxisLims = np.array(((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600)); #get those xAxis limits, ensures everything is aligned

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im = ax.pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
    cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar.ax.tick_params(labelsize=FONT_axisTick);
    cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cax.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    # string_Title = avg_anyAngle_dataType+' Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #     str(np.round(avg_anyAngle_Width,2))+' arcdeg'; #create mecha title
    # ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax.set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax.set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    #makes plotting look better
    if( (np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[0] = np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    if( (np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[1] = np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    # ax.set_xlim(xAxisLims); #set the xlims now
    yAxisLims = np.array(ax.get_ylim()); #get the current y axis limits
    if( (np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)) - np.min(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[0] = np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    if( (np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.max(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[1] = np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    ax.set_ylim(yAxisLims); #set the ylims now
    
    time_autoTick = (np.ceil((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.floor((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( time_autoTick > 24 ):
        time_autoTick = 24; #sets the tick setting to 15 arcdegrees per tick
    elif( time_autoTick > 12 ):
        time_autoTick = 12; #sets the tick setting to 15 arcdegrees per tick
    elif( time_autoTick > 6 ):
        time_autoTick = 6; #sets the tick setting to 10 arcdegrees per tick
    elif( time_autoTick > 4 ):
        time_autoTick = 4; #sets the tick setting to 5 arcdegrees per tick
    elif( time_autoTick > 2 ):
        time_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( time_autoTick >= 1 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        time_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    else:
        time_autoTick = ((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600 - (np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600)/13; #just goes for it if it's a super tiny range
    #END IF
    xAxisTicks = np.arange( (np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) , \
            (np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) + time_autoTick , \
            time_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax.set_xticks(xAxisTicks); #set x axis ticks
    ax.set_xlim(xAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    
    avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 20 ):
        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 30 arcdegrees per tick
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
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/17; #just goes for it if it's a super tiny range
        elif(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude'): #if Y axis is longitude, use longitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/17; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))+avg_anyAngle_Range_Chunks_Long_Plot_autoTick,avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax.set_yticks(yAxisTicks); #set x axis ticks
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
            ax.plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
            ax.plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    #END IF
    
    figFitter(fig); #fit the fig fast
    # if( np.any(np.abs(np.array(plotLongRange)) >= 100) & (avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude') ): #if true, longitude
    #     fig.subplots_adjust(left = 0.092, right = 0.9178, top = 0.980, bottom = 0.085); #sets padding to small numbers for minimal white space
    # elif( np.any(np.array(plotLatRange) < 0) & (avg_anyAngle_Range_Chunks_Long_Plot_Name != 'Longitude') ):
    #     fig.subplots_adjust(left = 0.080, right = 0.9178, top = 0.980, bottom = 0.085); #sets padding to small numbers for minimal white space
    # else:
    #     fig.subplots_adjust(left = 0.060, right = 0.9178, top = 0.980, bottom = 0.085); #sets padding to small numbers for minimal white space
    # #END IF
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig(folder[3]+'\TEC_avgKeo.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff
       
#no return, plot is the return!




