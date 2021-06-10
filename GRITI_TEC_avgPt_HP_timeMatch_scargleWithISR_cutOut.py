"""
GOAL: Plot only Lomb-Scargle of ISR SNR HP at 300 km as well as ISR SNR and ISR SNR HP for comparsion and looking
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from subfun_lombscargle import subfun_lombscargle

def GRITI_TEC_avgPt_HP_timeMatch_scargleWithISR_cutOut(avgPt_vTEC_timeMatch_HP_cutOut,avgPt_vTEC_timeMatch_time_cutOut,Zenith_time_cutOut,Zenith_SNR_hp_altAvgd_cutOut,MISA_time_cutOut,MISA_SNR_hp_altAvgd_cutOut,pointAltitude,avgPt_ISRavgAlt,time_cutout_range,filter_cutoffPeriod,avgPt_pointRadius,plot_Period_Lim,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,PLOT_lineWidth):

    avgPt_vTEC_scargPeriod, avgPt_vTEC_scargPower, avgPt_vTEC_scarggf = subfun_lombscargle((avgPt_vTEC_timeMatch_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , avgPt_vTEC_timeMatch_HP_cutOut); #scargle that data
    avgPt_vTEC_scargPeriod = avgPt_vTEC_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
    
    Zenith_scargPeriod, Zenith_scargPower, Zenith_scarggf = subfun_lombscargle((Zenith_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , Zenith_SNR_hp_altAvgd_cutOut); #scargle that data
    Zenith_scargPeriod = Zenith_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
    
    MISA_scargPeriod, MISA_scargPower, MISA_scarggf = subfun_lombscargle((MISA_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , MISA_SNR_hp_altAvgd_cutOut); #scargle that data
    MISA_scargPeriod = MISA_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
    
    #-----Plot ISR SNR HP results as a RTI-----
    #Plot just the ISR SNR HP results
    fig, ax = plt.subplots(nrows=3, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    ax[2].set_aspect('auto');
    
    #~~~~~~~~~~ delta-vTEC STUFF STUFF ~~~~~~~~~~
    #-----delta-vTEC lomb-scargle of delta-vTEC avg'd 50 km radius and matched to ISR time interval and high-passed-----
    ax[0].plot( avgPt_vTEC_scargPeriod, avgPt_vTEC_scargPower ,color='xkcd:fire engine red' ,linewidth=PLOT_lineWidth, linestyle='-'); #plot
    ax[0].plot( avgPt_vTEC_scargPeriod, np.tile(avgPt_vTEC_scarggf,np.size(avgPt_vTEC_scargPeriod)) , color="xkcd:grey" ); #plot
    
#    ax[0].set_xlabel("Periods (min)",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('delta-vTEC Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    string_Title = 'Scargle of delta-vTEC Avg''d in a '+str(avgPt_pointRadius)+' km radius & HP''d w/ '+str(filter_cutoffPeriod)+' hr cutoff & Time Avg''d to match ISR Interval - Time Cutout '+str(np.min(time_cutout_range))+' to '+str(np.max(time_cutout_range))+' hrs'; #create mecha title
    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    ax[0].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    if( np.all(np.isnan(avgPt_vTEC_scargPower)) == 0 ): #avoids an error
        ax[0].set_ylim( (0, np.max(avgPt_vTEC_scargPower)+0.1*np.max(avgPt_vTEC_scargPower) ) ); #set x axis limits
    #END IF
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    
    #~~~~~~~~~~ ZENITH STUFF ~~~~~~~~~~
    #-----ZENITH lomb-scargle of SNR high-passed and averaged +/-25 km at 300km-----
    ax[1].plot( Zenith_scargPeriod, Zenith_scargPower ,color='xkcd:electric blue' ,linewidth=PLOT_lineWidth, linestyle='-'); #plot
    ax[1].plot( Zenith_scargPeriod, np.tile(Zenith_scarggf,np.size(Zenith_scargPeriod)) , color="xkcd:grey" ); #plot
    
#    ax[1].set_xlabel("Periods (min)",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel('Zenith Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    string_Title = 'Scargle of Zenith SNR HP''d w/ '+str(filter_cutoffPeriod)+' hr cutoff & Avg''d +/- '+str(avgPt_ISRavgAlt)+' around '+str(pointAltitude)+' km altitude - Time Cutout '+str(np.min(time_cutout_range))+' to '+str(np.max(time_cutout_range))+' hrs'; #create mecha title
    ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    ax[1].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    if( np.all(np.isnan(Zenith_scargPower)) == 0 ): #avoids an error
        ax[1].set_ylim( (0, np.max(Zenith_scargPower)+0.1*np.max(Zenith_scargPower) ) ); #set x axis limits
    #END IF
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks

    
    #~~~~~~~~~~ MISA STUFF ~~~~~~~~~~
    #-----MISA lomb-scargle of SNR high-passed and averaged +/- 25km at 300km-----
    ax[2].plot( MISA_scargPeriod, MISA_scargPower, color='xkcd:vivid green' ,linewidth=PLOT_lineWidth, linestyle ='-' ); #plot
    ax[2].plot( MISA_scargPeriod, np.tile(MISA_scarggf,np.size(MISA_scargPeriod)) , color="xkcd:grey" ); #plot
    
    ax[2].set_xlabel("Periods (min)",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2].set_ylabel('MISA Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    string_Title = 'Scargle of Zenith SNR HP''d w/ '+str(filter_cutoffPeriod)+' hr cutoff & Avg''d +/- '+str(avgPt_ISRavgAlt)+' around '+str(pointAltitude)+' km altitude - Time Cutout '+str(np.min(time_cutout_range))+' to '+str(np.max(time_cutout_range))+' hrs'; #create mecha title
    ax[2].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    ax[2].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    if( np.all(np.isnan(MISA_scargPower)) == 0 ): #avoids an error
        ax[2].set_ylim( (0, np.max(MISA_scargPower)+0.1*np.max(MISA_scargPower) ) ); #set x axis limits
    #END IF
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[2].set_xticks(xAxisTicks); #set x axis ticks
    
    #final plot adjusting stuff
    fig.subplots_adjust(left = 0.065, right = 0.975, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up