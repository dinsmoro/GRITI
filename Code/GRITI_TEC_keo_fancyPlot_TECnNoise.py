#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tick
from Code.GRITI_TEC_keo_noise import GRITI_TEC_keo_noise
from Code.GRITI_TEC_randomSynth import GRITI_TEC_randomSynth
from Code.subfun_figFitter import figFitter

def GRITI_TEC_keo_fancyPlot_TECnNoise(vTECChunked_anyAngleAvg,TEC_timeUnique,TEC_plotLimValu, \
        TEC_float,locFloat_lat, locFloat_long, locFloat_time, time_Ref, \
        noise_background_mean, noise_background_stdev, Re, avg_anyAngle_N, avg_anyAngle_45vsLatLong, \
        wave_latRange, wave_longRange, wave_N, wave_angle, wave_phase, wave_waveLength, wave_period, wave_amp, \
        colorMap,plotLatRange,plotLongRange,latMillstone,longMillstone,dateRange_dayNum_zeroHr, \
        avg_anyAngle,avg_anyAngle_Width,avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name, \
        avg_anyAngle_dataType,avg_anyAngle_plotLabel,dateRange_zeroHr,dateRange_zeroHr_monthName,dateRange_zeroHr_dayPostfix,\
        Zenith_time,  \
        FONT_grandioseFM, FONT_titleFM,FONT_axisTick,FONT_axisLabelFM, PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi):
    print('MAKING FANCY PLOT: TEC_avgAnyAngle_plot_wNoise IN fancyPlot FOLDER'); #report since you won't see anything
    
    letteringPositionX = -0.085; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.92; #set the X position of the lettering (e.g., a. b. c. ...)
    
    def cbarFormatter(x, y):
        cbarVal = '{:1.0e}'.format(x).replace('+',''); #get it in basic scientific format
        if( cbarVal[cbarVal.find('e')+1] == '0' ):
            cbarVal = cbarVal[:cbarVal.find('e')+1]+cbarVal[cbarVal.find('e')+2:]; #remove the 0 in '2e05' or something like that
        #END IF
        if( cbarVal[0] == '0' ):
            cbarVal = '0'; #just set to 0
        #END IF
        return cbarVal
    #END DEF
    
    #shared x-axis limits are important here
    xAxisLims = ((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(Zenith_time)-dateRange_dayNum_zeroHr[1]*86400)/3600); #get those xAxis limits, ensures everything is aligned
    
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

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig, ax = plt.subplots(nrows=2, ncols=1,figsize=(14,10.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    fig.subplots_adjust(hspace=0.080);
    divider0 = make_axes_locatable(ax[0]); #prep to add an axis
    cax0 = divider0.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    divider1 = make_axes_locatable(ax[1]); #prep to add an axis
    cax1 = divider1.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , avg_anyAngle_Range_Chunks_Long_Plot] ,
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
    
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    # string_Title = avg_anyAngle_dataType+' Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
    #     str(np.round(avg_anyAngle_Width,2))+' arcdeg'; #create mecha title
    # ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    # ax.set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[0].set_xlim(xAxisLims); #set the common x-axis limits
    #makes plotting look better
    # xAxisLims = np.array(ax[0].get_xlim()); #get the current x axis limits
    if( (np.round(xAxisLims[0]) - xAxisLims[0]) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[0] = np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    if( (np.round(xAxisLims[1]) - xAxisLims[1]) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[1] = np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    ax[0].set_xlim(xAxisLims); #set the xlims now
    yAxisLims = np.array(ax[0].get_ylim()); #get the current y axis limits
    if( (np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)) - np.min(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[0] = np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    if( (np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.max(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[1] = np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    ax[0].set_ylim(yAxisLims); #set the ylims now
    
    time_autoTick = (np.ceil((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.floor((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
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
        time_autoTick = ((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600 - (np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600)/17; #just goes for it if it's a super tiny range
    #END IF
    xAxisTicks = np.arange( (np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) , \
            (np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) + time_autoTick , \
            time_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    ax[0].set_xlim(xAxisLims); #set the xlims now again
    
    avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)))/17; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 20 ):
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
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/17; #just goes for it if it's a super tiny range
        elif(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude'): #if Y axis is longitude, use longitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/17; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))+avg_anyAngle_Range_Chunks_Long_Plot_autoTick,avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax[0].set_yticks(yAxisTicks); #set x axis ticks
    ax[0].set_ylim(yAxisLims); #set the ylims now
    
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
    ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text saying the day or nite
    
    #!!!-----Plot TEC noise verison of above plot-----!!!
    #Prep the TEC noise data
    TEC_noise = GRITI_TEC_randomSynth(TEC_float.shape[0],TEC_float[:,locFloat_lat],TEC_float[:,locFloat_long],TEC_float[:,locFloat_time], \
        noise_background_mean,noise_background_stdev,Re,dateRange_zeroHr, \
        plotLatRange,plotLongRange,0,0, \
        wave_latRange,wave_longRange,wave_N,wave_angle,wave_phase,wave_waveLength,wave_period,wave_amp, \
        FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,TEC_plotLimValu,2,FLG_plotStuff=0); #replace the delta-vTEC data with random data OR random data with synth waves embedded
    #Average the data
    avg_anyAngle = avg_anyAngle*180/np.pi; #convert to deg from rad again
    vTECChunked_anyAngleAvg, avg_anyAngle, avg_anyAngle_Width,  avg_anyAngle_Range_Chunks_Long_Plot, avg_anyAngle_Range_Chunks_Long_Plot_Name = \
        GRITI_TEC_keo_noise(plotLatRange,plotLongRange,TEC_noise,TEC_float,TEC_timeUnique,TEC_plotLimValu,
            colorMap,0,locFloat_time,locFloat_lat,locFloat_long, time_Ref, avg_anyAngle,avg_anyAngle_N,avg_anyAngle_Width,
            avg_anyAngle_45vsLatLong,0,0,dateRange_dayNum_zeroHr,0,
            0,0, 0, 0, 0,
            FONT_titleFM,FONT_axisTick,0,0,0,
            avg_anyAngle_dataType,avg_anyAngle_plotLabel,0,PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi,
            avg_anyAngle_polarMode=0,FLG_disablePlot=1); #average the noise data
        
        
    #Do other stuff
    pltHelprX, pltHelprY = np.meshgrid( (TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im1 = ax[1].pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax1.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
    cbar1.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar1.ax.tick_params(labelsize=FONT_axisTick);
    cbar1.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cax1.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax1.yaxis.label.set_font_properties(FONT_axisLabelFM);
    
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle*180/np.pi,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    # string_Title = avg_anyAngle_dataType+' Averaged on Angle of '+str(np.round(avg_anyAngle*180/np.pi,2))+' deg and Width of '+ \
    #     str(np.round(avg_anyAngle_Width,2))+' arcdeg'; #create mecha title
    # ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    # ax.set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    #makes plotting look better
    ax[1].set_xlim(xAxisLims); #set the xlims now
    ax[1].set_ylim(yAxisLims); #set the ylims now
    
    xAxisTicks = np.arange( (np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) , \
            (np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1]*86400)/3600),time_autoTick)) + time_autoTick , \
            time_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks
    ax[1].set_xlim(xAxisLims); #set the xlims now again

    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))+avg_anyAngle_Range_Chunks_Long_Plot_autoTick,avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax[1].set_yticks(yAxisTicks); #set x axis ticks
    ax[1].set_ylim(yAxisLims); #set the ylims now
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
            ax[1].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
            ax[1].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    #END IF
    ax[1].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label

    ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=FONT_grandioseFM, transform=ax[1].transAxes); #print the text saying the day or nite
    
    figFitter(fig); #fit the fig fast
    # if( np.any(np.abs(np.array(plotLongRange)) >= 100) & (avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude') ): #if true, longitude
    #     fig.subplots_adjust(left = 0.092, right = 0.911, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    # elif( np.any(np.array(plotLatRange) < 0) & (avg_anyAngle_Range_Chunks_Long_Plot_Name != 'Longitude') ):
    #     fig.subplots_adjust(left = 0.075, right = 0.911, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    # else:
    #     fig.subplots_adjust(left = 0.075, right = 0.911, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    # #END IF
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig(folder[3]+'TEC_avgKeo&Noise.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff
       
#no return, plot is the return!




