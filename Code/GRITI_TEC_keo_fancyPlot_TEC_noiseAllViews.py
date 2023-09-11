"""
GOAL: Plot only ISR POPL HP RTI
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Code.subfun_monthNum_to_word import subfun_monthNum_to_word
from Code.GRITI_TEC_randomSynth import GRITI_TEC_randomSynth
from Code.GRITI_TEC_keo_noise import GRITI_TEC_keo_noise
from Code.GRITI_TEC_avgPt import GRITI_TEC_avgPt
from Code.GRITI_TEC_avgPt_timeMatch import GRITI_TEC_avgPt_timeMatch
from Code.subfun_figFitter import figFitter

def GRITI_TEC_keo_fancyPlot_TEC_noiseAllViews(time_Ref, Re,
        geoMap_projectionStyle, BasemapFixDir, colorMap,
        plotLatRange, plotLongRange,  plotLatRange_zoomZone, plotLongRange_zoomZone,
        plotLatRange_autoTick,plotLongRange_autoTick,plotLongRange_autoTick_Crunched, 
        TEC_timeUnique, TEC_float, locFloat_time, locFloat_lat, locFloat_long, locFloat_dTEC,
        TEC_plotLimValu, noise_background_mean, noise_background_stdev, avgPt_TECnoise_iterations,
        avgPt_pointRadius, avg_anyAngle_45vsLatLong, avg_anyAngle_plotLabel,
        gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size,
        dataReject, dataRejectOrig, dataRejectLimit, dataRejectLimitOrig, dataRejectMax,
        Zenith_time, Zenith_height, Zenith_POPL_hp, MISA_time, MISA_height, MISA_POPL_hp,
        pointAltitude, avgPt_ISRavgAlt, filter_cutoffPeriod,
        avgPt_coords,time_cutout_range,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,
        dateRange_zeroHr, dateRange_zeroHr_monthName, dateRange_zeroHr_dayPostfix, \
        FONT_grandioseFM, FONT_titleFM,FONT_axisTickFM,FONT_axisTick,FONT_axisLabelFM,
        PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi):
    
    print('MAKING FANCY PLOT: TEC_avgAnyAngle_plot_noiseAllViews IN fancyPlot FOLDER !!THIS WILL TAKE A WHILE!!'); #report since you won't see anything

    letteringPositionX = -0.1025; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.92; #set the X position of the lettering (e.g., a. b. c. ...)
    
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
    
    #-----Plot ISR POPL HP results as a RTI-----
    #Plot just the ISR POPL HP results
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig = plt.figure(figsize=(14,18.3),dpi=journal_dpi);
    gridr = gridspec.GridSpec(nrows=7, ncols=1, figure=fig);
    gridr.update(hspace=0.22); # set the spacing between axes. 
    fig.add_subplot(gridr[0]); #noise time series
    fig.add_subplot(gridr[1:3]); #keogram east coast
    fig.add_subplot(gridr[3:5]); #keogram world 90
    fig.add_subplot(gridr[5:7]); #keogram world 0
    ax = fig.axes; #get a list of the axes
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    ax[2].set_aspect('auto');
    ax[3].set_aspect('auto');
    
    
    print('noiseAllViews fancyPlot #1 - Time Series'); #report
    #----Prep time series data-----
    #Now for the ISR stuff    
    Zenith_height_atISRavgAlt = np.where( np.min(np.abs( Zenith_height - pointAltitude )) == np.abs( Zenith_height - pointAltitude ) )[0][0]; #get the index of where Zenith_height is closest to pointAltitude (in km)
    Zenith_height_atISRavgAltIndexes = np.array( ( np.where(np.min(np.abs( (Zenith_height[Zenith_height_atISRavgAlt]-avgPt_ISRavgAlt) - Zenith_height)) == np.abs(Zenith_height[Zenith_height_atISRavgAlt]-avgPt_ISRavgAlt - Zenith_height) )[0][0] , \
                                    np.where(np.min(np.abs( (Zenith_height[Zenith_height_atISRavgAlt]+avgPt_ISRavgAlt) - Zenith_height)) == np.abs(Zenith_height[Zenith_height_atISRavgAlt]+avgPt_ISRavgAlt - Zenith_height) )[0][0] ) ); #get the upper and lower indexes of height to average (+/- avgPt_ISRavgAlt in km)
    Zenith_POPL_hp_altAvgd = np.mean( Zenith_POPL_hp[Zenith_height_atISRavgAltIndexes[0]:Zenith_height_atISRavgAltIndexes[1]+1,:] , axis=0 ); #average +/- avgPt_ISRavgAlt in km around the center altitude pointAltitude in km
    
    MISA_height_atISRavgAlt = np.where( np.min(np.abs( MISA_height - pointAltitude )) == np.abs( MISA_height - pointAltitude ) )[0][0]; #get the index of where MISA_height is closest to pointAltitude (in km)
    MISA_height_atISRavgAltIndexes = np.array( ( np.where(np.min(np.abs( (MISA_height[MISA_height_atISRavgAlt]-avgPt_ISRavgAlt) - MISA_height)) == np.abs(MISA_height[MISA_height_atISRavgAlt]-avgPt_ISRavgAlt - MISA_height) )[0][0] , 
                                    np.where(np.min(np.abs( (MISA_height[MISA_height_atISRavgAlt]+avgPt_ISRavgAlt) - MISA_height)) == np.abs(MISA_height[MISA_height_atISRavgAlt]+avgPt_ISRavgAlt - MISA_height) )[0][0] ) ); #get the upper and lower indexes of height to average (+/- avgPt_ISRavgAlt in km)
    MISA_POPL_hp_altAvgd = np.mean( MISA_POPL_hp[MISA_height_atISRavgAltIndexes[0]:MISA_height_atISRavgAltIndexes[1]+1,:] , axis=0 ); #average +/- avgPt_ISRavgAlt in km around the center altitude pointAltitude in km
    
    #cut out the times
    time_cutout_indexesZ = np.array( ( np.where(np.min(np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) )) == np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) ) )[0][0] , \
        np.where(np.min(np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) )) == np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) ) )[0][0] ) ); #get the indexes for that time cutout range
    Zenith_time_cutOut = Zenith_time[time_cutout_indexesZ[0]:time_cutout_indexesZ[1]+1];
    Zenith_POPL_hp_altAvgd_cutOut = Zenith_POPL_hp_altAvgd[time_cutout_indexesZ[0]:time_cutout_indexesZ[1]+1];
    
    time_cutout_indexes = np.array( ( np.where(np.min(np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) )) == np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) ) )[0][0] , 
        np.where(np.min(np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) )) == np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) ) )[0][0] ) ); #get the indexes for that time cutout range
    MISA_time_cutOut = MISA_time[time_cutout_indexes[0]:time_cutout_indexes[1]+1];
    MISA_POPL_hp_altAvgd_cutOut = MISA_POPL_hp_altAvgd[time_cutout_indexes[0]:time_cutout_indexes[1]+1];
    
    avgPt_TECnoise_timeMatch_HP_cutOut_avg = np.zeros( (avgPt_TECnoise_iterations, time_cutout_indexesZ[1]-time_cutout_indexesZ[0]+1) ); #preallocate
    pointRadiusAngular = (avgPt_pointRadius/Re)*180/np.pi; #get the angular radius to get a small subset of points to deal with
    k = ((avgPt_coords[0,0]-pointRadiusAngular <= TEC_float[:,locFloat_lat]) & (avgPt_coords[0,0]+pointRadiusAngular >= TEC_float[:,locFloat_lat])) & \
        ((avgPt_coords[0,1]-pointRadiusAngular <= TEC_float[:,locFloat_long]) & (avgPt_coords[0,1]+pointRadiusAngular >= TEC_float[:,locFloat_long])); #get only east coast to lower calcs needed
    for i in range(0,avgPt_TECnoise_iterations):
        TEC_noise = GRITI_TEC_randomSynth(TEC_float[k,locFloat_lat].shape[0],TEC_float[k,locFloat_lat],TEC_float[k,locFloat_long],TEC_float[k,locFloat_time], \
            noise_background_mean,noise_background_stdev,Re,dateRange_zeroHr, \
            plotLatRange,plotLongRange,0,0, \
            0,0,0,0,0,0,0,0, \
            FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,TEC_plotLimValu,1,FLG_plotStuff=0); #replace the delta-vTEC data with random data 
        
        avgPt_TECnoise, _, avgPt_TECnoise_time, _, _, _  = \
            GRITI_TEC_avgPt(TEC_timeUnique,TEC_float[k,locFloat_lat],TEC_float[k,locFloat_long],TEC_float[k,locFloat_time],TEC_noise, \
            avgPt_coords[0,:],avgPt_pointRadius,Re,dateRange_dayNum_zeroHr, \
            dataReject,dataRejectOrig,dataRejectLimit,dataRejectLimitOrig,dataRejectMax,FLG_report=0); #average points in a radius
            
        _, avgPt_TECnoise_timeMatch_HP, avgPt_vTEC_timeMatch_time = GRITI_TEC_avgPt_timeMatch(avgPt_TECnoise,avgPt_TECnoise_time,Zenith_time,dateRange_dayNum_zeroHr,filter_cutoffPeriod=filter_cutoffPeriod);
    
        avgPt_TECnoise_timeMatch_HP_cutOut = avgPt_TECnoise_timeMatch_HP[time_cutout_indexesZ[0]:time_cutout_indexesZ[1]+1];
        avgPt_vTEC_timeMatch_time_cutOut = avgPt_vTEC_timeMatch_time[time_cutout_indexesZ[0]:time_cutout_indexesZ[1]+1];
        
        pwr_TECZN = np.sqrt(1/avgPt_TECnoise_timeMatch_HP_cutOut.size*np.sum(avgPt_TECnoise_timeMatch_HP_cutOut**2)); #estimate power of signal
        
        avgPt_TECnoise_timeMatch_HP_cutOut_avg[i,:] = 1/pwr_TECZN*avgPt_TECnoise_timeMatch_HP_cutOut; #record
    #END IF
    avgPt_TECnoise_timeMatch_HP_cutOut = np.mean(avgPt_TECnoise_timeMatch_HP_cutOut_avg,axis=0); #get output
    # avgPt_TECnoise_timeMatch_HP_cutOut = avgPt_TECnoise_timeMatch_HP_cutOut_avg[2,:]; #get output

    #----PLOT TIME SERIES DATA----    
    pwr_TECZN = np.sqrt(1/avgPt_TECnoise_timeMatch_HP_cutOut.size*np.sum(avgPt_TECnoise_timeMatch_HP_cutOut**2)); #estimate power of signal
    # pwr_TECZN = 1; #since it was made of averaged power signals, don't re-boost it
    Zenith_pwr = np.sqrt(1/Zenith_POPL_hp_altAvgd_cutOut.size*np.sum(Zenith_POPL_hp_altAvgd_cutOut**2)); #estimate power of signal
    MISA_pwr = np.sqrt(1/MISA_POPL_hp_altAvgd_cutOut.size*np.sum(MISA_POPL_hp_altAvgd_cutOut**2)); #estimate power of signal
    
    #Real quick side move to calc correlation coefficients
    # R_ZvTECZN = np.corrcoef(1/Zenith_pwr*Zenith_POPL_hp_altAvgd_cutOut,1/pwr_TECZN*avgPt_TECnoise_timeMatch_HP_cutOut)[0,1];
    # R_MvTECZN = np.corrcoef(1/MISA_pwr*MISA_POPL_hp_altAvgd_cutOut,1/pwr_TECZN*avgPt_TECnoise_timeMatch_HP_cutOut)[0,1];
    # print('ZvTECZN: '+str(np.round(R_ZvTECZN,2))+' | R_MvTECZN: '+str(np.round(R_MvTECZN,2))); #report
    
    line1, = ax[0].plot( (avgPt_vTEC_timeMatch_time_cutOut - dateRange_dayNum_zeroHr[1]*86400)/3600 , 1/pwr_TECZN*avgPt_TECnoise_timeMatch_HP_cutOut , color='xkcd:fire engine red' ,linewidth=PLOT_lineWidthRegularPlus, linestyle='-', marker='o'); # plot some time series
    line2, = ax[0].plot( (Zenith_time_cutOut - dateRange_dayNum_zeroHr[1]*86400)/3600 , 1/Zenith_pwr*Zenith_POPL_hp_altAvgd_cutOut , color='xkcd:electric blue' ,linewidth=PLOT_lineWidthRegularPlus, linestyle='--'); # plot some time series
    line3, = ax[0].plot( (MISA_time_cutOut - dateRange_dayNum_zeroHr[1]*86400)/3600 , 1/MISA_pwr*MISA_POPL_hp_altAvgd_cutOut , color='xkcd:vivid green' ,linewidth=PLOT_lineWidthRegularPlus, linestyle ='-'); # plot some time series
    
    ax[0].legend((line2, line3, line1), ('Zenith N$_{e^-}$', 'MISA N$_{e^-}$', 'delta-vTEC [Noise]') ,loc='lower right', ncol=3, prop=FONT_axisLabelFM);
    # string_Title = 'delta-vTEC ISR Interval Matched AVG''d High-passed AND Zenith/MISA POPL HP AVG''d +/-'+str(avgPt_ISRavgAlt)+' km around '+str(pointAltitude)+' km at '+str(np.round(avgPt_coords[0],2))+' deg lat/'+str(np.round(avgPt_coords[1],2))+'deg long'; #create mecha title
    # ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
#    ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('Amplitude',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1]*86400)/3600),2)) , \
            (np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1]*86400)/3600) - np.mod(np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1]*86400)/3600),2)) + 2 , \
            1); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    ax[0].set_xlim( ((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1]*86400)/3600) ); #set x axis limits
    ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text saying the day or nite
    ax[0].spines['right'].set_visible(False); #turn off box lines
    ax[0].spines['top'].set_visible(False); #turn off box lines
    
        
    xAxisLims = ((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #get those xAxis limits, ensures everything is aligned
    #makes plotting look better
    if( (np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[0] = np.round((np.min(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    if( (np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) - (np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600) < 10**-2 ): #fix x axis stuff (time)
        xAxisLims[1] = np.round((np.max(time_Ref)-dateRange_dayNum_zeroHr[1]*86400)/3600); #force the rounded value
    #END IF
    
    print('noiseAllViews fancyPlot #2 - East Coast Keo'); #report
    #-----Prep East Coast TEC Keogram-----
    k = ((plotLatRange_zoomZone[0] <= TEC_float[:,locFloat_lat]) & (plotLatRange_zoomZone[1] >= TEC_float[:,locFloat_lat])) & \
        ((plotLongRange_zoomZone[0] <= TEC_float[:,locFloat_long]) & (plotLongRange_zoomZone[1] >= TEC_float[:,locFloat_long])); #get only east coast to lower calcs needed
    TEC_noise = GRITI_TEC_randomSynth(TEC_float[k,locFloat_lat].shape[0],TEC_float[k,locFloat_lat],TEC_float[k,locFloat_long],TEC_float[k,locFloat_time], \
        noise_background_mean,noise_background_stdev,Re,dateRange_zeroHr, \
        plotLatRange,plotLongRange,0,0, \
        0,0,0,0,0,0,0,0, \
        FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,TEC_plotLimValu,1,FLG_plotStuff=0); #replace the delta-vTEC data with random data 
    (vTECChunked_anyAngleAvg, avg_anyAngle,avg_anyAngle_Width, \
    avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name) = \
        GRITI_TEC_keo_noise(plotLatRange_zoomZone,plotLongRange_zoomZone,TEC_noise,TEC_float[k,:],TEC_timeUnique,\
            TEC_plotLimValu,'jet',locFloat_dTEC,locFloat_time,locFloat_lat,locFloat_long,90, \
            200,360,avg_anyAngle_45vsLatLong,avgPt_coords,geoMap_projectionStyle,\
            dateRange_dayNum_zeroHr,plotLatRange_autoTick,plotLongRange_autoTick,plotLongRange_autoTick_Crunched, gif_Millstone_Marker, gif_Millstone_Marker_Color, \
            gif_Millstone_Marker_Size,FONT_titleFM,FONT_axisTick,FONT_axisTickFM,FONT_axisLabelFM,BasemapFixDir,\
            'delta-vTEC','delta-vTEC [TECU]',0,PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi, \
            avg_anyAngle_polarMode=0,FLG_disablePlot=1);
    
    #-----Plot TEC East Coast results as a Keogram-----
    #Plot just the TEC
    divider1 = make_axes_locatable(ax[1]); #prep to add an axis
    cax1 = divider1.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (np.append(TEC_timeUnique,TEC_timeUnique[-1]+np.median(np.diff(TEC_timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im1 = ax[1].pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax1.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
    cbar1.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar1.ax.tick_params(labelsize=FONT_axisTick);
    cbar1.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cax1.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax1.yaxis.label.set_font_properties(FONT_axisLabelFM);

    # ax[1].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    #makes plotting look better
    yAxisLims = np.array(ax[1].get_ylim()); #get the current y axis limits
    if( (np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)) - np.min(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[0] = np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    if( (np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.max(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[1] = np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    ax[1].set_ylim(yAxisLims); #set the ylims now
    
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
    ax[1].set_xticks(xAxisTicks); #set x axis ticks
    ax[1].set_xlim(xAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    
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
    ax[1].set_yticks(yAxisTicks); #set x axis ticks
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= avgPt_coords[0,1]) & (np.max(plotLongRange) >= avgPt_coords[0,1]) ): #only plot if it's in the long range specified
            ax[1].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(avgPt_coords[0,1],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= avgPt_coords[0,0]) & (np.max(plotLatRange) >= avgPt_coords[0,0]) ): #only plot if it's in the lat range specified
            ax[1].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(avgPt_coords[0,0],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    #END IF
    ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=FONT_grandioseFM, transform=ax[1].transAxes); #print the text saying the day or nite


    #get TEC noise that covers the entire USA now
    TEC_noise = GRITI_TEC_randomSynth(TEC_float.shape[0],TEC_float[:,locFloat_lat],TEC_float[:,locFloat_long],TEC_float[:,locFloat_time], \
        noise_background_mean,noise_background_stdev,Re,dateRange_zeroHr, \
        plotLatRange,plotLongRange,0,0, \
        0,0,0,0,0,0,0,0, \
        FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,TEC_plotLimValu,1,FLG_plotStuff=0); #replace the delta-vTEC data with random data 

    print('noiseAllViews fancyPlot #3 - USA 0 Deg Keo'); #report
    #-----Prep USA 0 degree average TEC Keogram-----
    (vTECChunked_anyAngleAvg, avg_anyAngle,avg_anyAngle_Width, \
    avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name) = \
        GRITI_TEC_keo_noise(plotLatRange,plotLongRange,TEC_noise,TEC_float,TEC_timeUnique,\
            TEC_plotLimValu,'jet',locFloat_dTEC,locFloat_time,locFloat_lat,locFloat_long,0, \
            200,360,avg_anyAngle_45vsLatLong,avgPt_coords,geoMap_projectionStyle,\
            dateRange_dayNum_zeroHr,plotLatRange_autoTick,plotLongRange_autoTick,plotLongRange_autoTick_Crunched, gif_Millstone_Marker, gif_Millstone_Marker_Color, \
            gif_Millstone_Marker_Size,FONT_titleFM,FONT_axisTick,FONT_axisTickFM,FONT_axisLabelFM,BasemapFixDir,\
            'delta-vTEC','delta-vTEC [TECU]',0,PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi, \
            avg_anyAngle_polarMode=0,FLG_disablePlot=1);
    
    #-----Plot USA 0 degree average results as a Keogram-----
    #Plot just the TEC
    divider2 = make_axes_locatable(ax[2]); #prep to add an axis
    cax2 = divider2.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (np.append(TEC_timeUnique,TEC_timeUnique[-1]+np.median(np.diff(TEC_timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im2 = ax[2].pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar2 = fig.colorbar(im2, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax2.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
    cbar2.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar2.ax.tick_params(labelsize=FONT_axisTick);
    cbar2.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cax2.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax2.yaxis.label.set_font_properties(FONT_axisLabelFM);

    # ax[2].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    #makes plotting look better
    yAxisLims = np.array(ax[2].get_ylim()); #get the current y axis limits
    if( (np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)) - np.min(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[0] = np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    if( (np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.max(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[1] = np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    ax[2].set_ylim(yAxisLims); #set the ylims now
    
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
    ax[2].set_xticks(xAxisTicks); #set x axis ticks
    ax[2].set_xlim(xAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    
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
    ax[2].set_yticks(yAxisTicks); #set x axis ticks
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= avgPt_coords[0,1]) & (np.max(plotLongRange) >= avgPt_coords[0,1]) ): #only plot if it's in the long range specified
            ax[2].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(avgPt_coords[0,1],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= avgPt_coords[0,0]) & (np.max(plotLatRange) >= avgPt_coords[0,0]) ): #only plot if it's in the lat range specified
            ax[2].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(avgPt_coords[0,0],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    #END IF
    ax[2].text( letteringPositionX, letteringPositionY, 'c.', color='r', fontproperties=FONT_grandioseFM, transform=ax[2].transAxes); #print the text saying the day or nite

    print('noiseAllViews fancyPlot #4 - USA 90 Deg Keo'); #report
    #-----Prep USA 90 degree average TEC Keogram-----
    (vTECChunked_anyAngleAvg, avg_anyAngle,avg_anyAngle_Width, \
    avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name) = \
        GRITI_TEC_keo_noise(plotLatRange,plotLongRange,TEC_noise,TEC_float,TEC_timeUnique,\
            TEC_plotLimValu,'jet',locFloat_dTEC,locFloat_time,locFloat_lat,locFloat_long,90, \
            200,360,avg_anyAngle_45vsLatLong,avgPt_coords,geoMap_projectionStyle,\
            dateRange_dayNum_zeroHr,plotLatRange_autoTick,plotLongRange_autoTick,plotLongRange_autoTick_Crunched, gif_Millstone_Marker, gif_Millstone_Marker_Color, \
            gif_Millstone_Marker_Size,FONT_titleFM,FONT_axisTick,FONT_axisTickFM,FONT_axisLabelFM,BasemapFixDir,\
            'delta-vTEC','delta-vTEC [TECU]',0,PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi, \
            avg_anyAngle_polarMode=0,FLG_disablePlot=1);
    
    #-----Plot USA 0 degree average results as a Keogram-----
    #Plot just the TEC
    divider3 = make_axes_locatable(ax[3]); #prep to add an axis
    cax3 = divider3.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (np.append(TEC_timeUnique,TEC_timeUnique[-1]+np.median(np.diff(TEC_timeUnique))) - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im3 = ax[3].pcolormesh(pltHelprX, pltHelprY,  vTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar3 = fig.colorbar(im3, cax=cax3, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax3.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f')); #force a rounded format
    cbar3.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar3.ax.tick_params(labelsize=FONT_axisTick);
    cbar3.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
    cax3.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),11)); #create useful tick marks
    cax3.yaxis.label.set_font_properties(FONT_axisLabelFM);

    ax[3].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[3].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    #makes plotting look better
    yAxisLims = np.array(ax[3].get_ylim()); #get the current y axis limits
    if( (np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)) - np.min(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[0] = np.round(np.min(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    if( (np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.max(avg_anyAngle_Range_Chunks_Long_Plot)) < 10**-2 ): #fix y axis stuff (lat or longitude)
        yAxisLims[1] = np.round(np.max(avg_anyAngle_Range_Chunks_Long_Plot)); #force the rounded value
    #END IF
    ax[3].set_ylim(yAxisLims); #set the ylims now
    
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
    ax[3].set_xticks(xAxisTicks); #set x axis ticks
    ax[3].set_xlim(xAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    
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
    ax[3].set_yticks(yAxisTicks); #set x axis ticks
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= avgPt_coords[0,1]) & (np.max(plotLongRange) >= avgPt_coords[0,1]) ): #only plot if it's in the long range specified
            ax[3].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(avgPt_coords[0,1],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= avgPt_coords[0,0]) & (np.max(plotLatRange) >= avgPt_coords[0,0]) ): #only plot if it's in the lat range specified
            ax[3].plot( np.linspace(np.min((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),np.max((TEC_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(avgPt_coords[0,0],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=PLOT_lineWidthSmol); #plots a point with a black line
        #END IF
    #END IF
    ax[3].text( letteringPositionX, letteringPositionY, 'd.', color='r', fontproperties=FONT_grandioseFM, transform=ax[3].transAxes); #print the text saying the day or nite

    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.092, right = 0.912, top = 0.988, bottom = 0.045); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig(folder[3]+'TEC_timeSeries&avgKeo_allNoise.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff