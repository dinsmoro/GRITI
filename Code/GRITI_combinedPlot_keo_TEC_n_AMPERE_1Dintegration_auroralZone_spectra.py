#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from Code.GRITI_TEC_avgPt_timeMatch import GRITI_TEC_avgPt_timeMatch
from Code.subfun_lombscargle import subfun_lombscargle
from Code.subfun_highpass import subfun_highpass
from scipy import signal
import warnings

def GRITI_combinedPlot_keo_TEC_n_AMPERE_1Dintegration_auroralZone_spectra(vTECChunked_anyAngleAvg,TEC_timeUnique,TEC_plotLimValu,colorMap, \
        AMPERE_data,AMPERE_timeUnique,locAMPERE_time,locAMPERE_lat,locAMPERE_long,AMPERE_plot_index,AMPERE_plot_indexes,AMPERE_plot_labels, \
        plotLatRange,plotLongRange,latMillstone,longMillstone,dateRange_dayNum_zeroHr,avg_anyAngle,avg_anyAngle_Width, \
        avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name,avg_anyAngle_dataType,avg_anyAngle_plotLabel, \
        time_Ref,time_Reference,dateRange,dateRange_zeroHr,dateRange_zeroHr_monthName,dateRange_zeroHr_dayPostfix,time_cutout_range_delay, \
        PLOT_lineWidth,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,settings,FLG_avg_anyAngle_Scargle_FFT):    
    #Unpack settings
    settings_spectra = settings['spectra']; #get the spectra settings
    
    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth[0]; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth[1]; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth[2]; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth[3]; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth[4]; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth[5]; #get the line widths
    #Unpack everything else
    plot_Period_Lim = settings_spectra['period limit max']; #get the period limit
    plot_period_min = settings_spectra['period limit min']; #get the period limit
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(TEC_plotLimValu) == 1 ):
        TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    #END IF
    
    #Prep to spectrum
    #-----Get TEC data line-----
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        k = np.where( np.min(np.abs(longMillstone-avg_anyAngle_Range_Chunks_Long_Plot)) == np.abs(longMillstone-avg_anyAngle_Range_Chunks_Long_Plot))[0].item(); #get
    else: #ottherwise latitude
        k = np.where( np.min(np.abs(latMillstone-avg_anyAngle_Range_Chunks_Long_Plot)) == np.abs(latMillstone-avg_anyAngle_Range_Chunks_Long_Plot))[0].item(); #get
    #END IF
    vTECChunked_anyAngleAvg_lined = vTECChunked_anyAngleAvg[:,k]; #get just a line of data at the desired lat or long
    
    #-----Plot AMPERE results as a 1D line-----
    AMPERE_timeUnique_hr = (AMPERE_timeUnique - dateRange_dayNum_zeroHr[1])*24; #hr, convert to hr with 0 hr at specified day
    
    if( np.mod(np.round(np.min(AMPERE_timeUnique_hr)),2) == 0 ):
        AMPERE_time_hr_axis_min = np.round(np.min(AMPERE_timeUnique_hr)); #is even, good to go
    else:
        AMPERE_time_hr_axis_min = np.round(np.min(AMPERE_timeUnique_hr))+1; #is odd, make even
    #END IF
    if( np.mod(np.round(np.max(AMPERE_timeUnique_hr)),2) == 0 ):
        AMPERE_time_hr_axis_max = np.round(np.max(AMPERE_timeUnique_hr)); #is even, good to go
    else:
        AMPERE_time_hr_axis_max = np.round(np.max(AMPERE_timeUnique_hr))-1; #is odd, make even
    #END IF
    
    if( (np.min(plotLatRange) >= 0) & (np.max(plotLatRange) >= 0) ):
        #northern hemisphere
        kInRange = AMPERE_data[:,locAMPERE_lat] >= 0; #get data in the range
    else:
        #southern hemisphere
        kInRange = AMPERE_data[:,locAMPERE_lat] <= 0; #get data in the range
    #END IF
    
    AMPERE_jouleHeating_integrate = np.zeros( AMPERE_timeUnique_hr.size , dtype=np.float64); #prep integrated joule heating
    for i in range(AMPERE_timeUnique_hr.size):
        k = AMPERE_timeUnique[i] == AMPERE_data[:,locAMPERE_time]; #get the right time
        AMPERE_jouleHeating_integrate[i] = np.sum((AMPERE_data[k&kInRange,AMPERE_plot_index])); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
    #END FOR i
    AMPERE_jouleHeating_integrate = np.log10(AMPERE_jouleHeating_integrate); #log it
    
    #make sure FFT can happen if it is on
    if( (FLG_avg_anyAngle_Scargle_FFT == 1) & ((np.isnan(AMPERE_jouleHeating_integrate).sum() > 0) | (np.isnan(vTECChunked_anyAngleAvg_lined).sum() > 0)) ):
        #if there are data gaps, data needs to be scargled
        # FLG_avg_anyAngle_Scargle_FFT = 0;
        #also NaNs need to be yeeted
        if( np.isnan(AMPERE_jouleHeating_integrate).sum() > 0 ):
            k = np.logical_not(np.isnan(AMPERE_jouleHeating_integrate));
            AMPERE_jouleHeating_integrate = AMPERE_jouleHeating_integrate[k]; #remove NaNs
            AMPERE_timeUnique_hr = AMPERE_timeUnique_hr[k]; #remove NaNs
        elif( np.isnan(vTECChunked_anyAngleAvg_lined).sum() > 0 ):
            k = np.logical_not(np.isnan(vTECChunked_anyAngleAvg_lined));
            vTECChunked_anyAngleAvg_lined = vTECChunked_anyAngleAvg_lined[k]; #remove NaNs
            TEC_timeUnique = TEC_timeUnique[k]; #remove NaNs
        #END IF
    #END IF
    #Force TEC onto AMPERE data cadence
    if( np.isclose(np.median(np.diff(TEC_timeUnique)),np.median(np.diff(AMPERE_timeUnique))) == False ):
        #Match the data in the 1st input (and its time in the 2nd input) to the time scale given in the 3rd input time and return that data and that data's highpassed form
        _, vTECChunked_anyAngleAvg_lined_timeMatch_HP, TEC_timeUnique_timeMatch = GRITI_TEC_avgPt_timeMatch(vTECChunked_anyAngleAvg_lined,TEC_timeUnique,AMPERE_timeUnique-time_cutout_range_delay/24,dateRange_dayNum_zeroHr,filter_cutoffPeriod=settings_spectra['filter cutoff period']);    
    else:
        #-----Highpass the data to keep the power within the period range we want-----
        vTECChunked_anyAngleAvg_lined_timeMatch_HP = subfun_highpass((TEC_timeUnique - dateRange_dayNum_zeroHr[1])*24, vTECChunked_anyAngleAvg_lined, filter_cutoffPeriod=settings_spectra['filter cutoff period'], filter_order=settings_spectra['filter order'], windowType=settings_spectra['window type'], axisToUse=1);
    #END IF
    #-----Highpass the data to keep the power within the period range we want-----
    AMPERE_jouleHeating_integrate_HP = subfun_highpass(AMPERE_timeUnique_hr, AMPERE_jouleHeating_integrate, filter_cutoffPeriod=settings_spectra['filter cutoff period'], filter_order=settings_spectra['filter order'], windowType=settings_spectra['window type'], axisToUse=1)

    #-----now scargle or FFT it-----
    if( FLG_avg_anyAngle_Scargle_FFT == 0 ):
        spectraName = 'Scargle'; #set the name for plotting
        #Scargle TEC
        vTEC_anyAngleAvg_scargPeriod, vTEC_anyAngleAvg_scargPower, vTEC_anyAngleAvg_scarggf = subfun_lombscargle((TEC_timeUnique_timeMatch - dateRange_dayNum_zeroHr[1])*24 , vTECChunked_anyAngleAvg_lined_timeMatch_HP); #scargle that data
        vTEC_anyAngleAvg_scargPeriod = vTEC_anyAngleAvg_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        #Scargle AMPERE
        AMPERE_integrated_scargPeriod, AMPERE_integrated_scargPower, AMPERE_integrated_scarggf = subfun_lombscargle(AMPERE_timeUnique_hr , AMPERE_jouleHeating_integrate); #scargle that data
        AMPERE_integrated_scargPeriod = AMPERE_integrated_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
    else:
        #otherwise FFT
        spectraName = 'FFT'; #set the name for plotting
        window = np.hamming(settings_spectra['windowLength']);
        #TEC TIME
        pwr_TEC = np.sqrt(1/vTECChunked_anyAngleAvg_lined_timeMatch_HP.size*np.sum(vTECChunked_anyAngleAvg_lined_timeMatch_HP**2)); #estimate power of signal
        Fs = 1/(np.median(np.diff(TEC_timeUnique_timeMatch))*24*60); #min, time delta in freq form
        [freqs_TEC,Cxx_TEC] = signal.welch(1/pwr_TEC*vTECChunked_anyAngleAvg_lined_timeMatch_HP ,window=window,noverlap=settings_spectra['noverlap'],nfft=settings_spectra['nfft']['6min'],fs=Fs);
        #AMPERE TIME
        pwr_AMPERE = np.sqrt(1/AMPERE_jouleHeating_integrate_HP.size*np.sum(AMPERE_jouleHeating_integrate_HP**2)); #estimate power of signal
        Fs = 1/(np.median(np.diff(AMPERE_timeUnique_hr))*60); #min, time delta in freq form
        [freqs_AMPERE,Cxx_AMPERE] = signal.welch(1/pwr_AMPERE*AMPERE_jouleHeating_integrate_HP ,window=window,noverlap=settings_spectra['noverlap'],nfft=settings_spectra['nfft']['6min'],fs=Fs);
    #END IF
    
    #-----Plot the spectra-----
    warnings.filterwarnings("ignore", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
    fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    
    #~~~~~~~~~~ delta-vTEC STUFF STUFF ~~~~~~~~~~
    if( FLG_avg_anyAngle_Scargle_FFT == 0 ):
        pTEC, = ax.plot( vTEC_anyAngleAvg_scargPeriod, vTEC_anyAngleAvg_scargPower ,color='xkcd:fire engine red' ,linewidth=PLOT_lineWidthRegular, linestyle='-'); #plot
        ax.plot( vTEC_anyAngleAvg_scargPeriod, np.tile(vTEC_anyAngleAvg_scarggf,np.size(vTEC_anyAngleAvg_scargPeriod)) , color="xkcd:grey" ); #plot
        # TEC_minPeriod = np.min(vTEC_anyAngleAvg_scargPeriod); #get the min period for the TEC
    else:
        pTEC, = ax.plot(1/freqs_TEC,(Cxx_TEC),color='xkcd:fire engine red',linewidth=PLOT_lineWidthRegular, linestyle='-');
        # TEC_minPeriod = np.min(1/freqs_TEC); #get the min period for the TEC
    #END IF

    # ax.set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    # if( np.all(np.isnan(vTEC_anyAngleAvg_scargPower)) == 0 ): #avoids an error
    #     ax.set_ylim( (0, np.max(vTEC_anyAngleAvg_scargPower)+0.1*np.max(vTEC_anyAngleAvg_scargPower) ) ); #set x axis limits
    # #END IF
    # xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    # ax.set_xticks(xAxisTicks); #set x axis ticks
    
    #~~~~~~~~~~ AMPERE STUFF STUFF ~~~~~~~~~~
    AMPERE_plot_label = AMPERE_plot_labels[np.where(AMPERE_plot_indexes == AMPERE_plot_index)[0].item()]; #get the label
    AMPERE_plot_label_noUnits = AMPERE_plot_label[0:AMPERE_plot_label.find('[')-1]; #remove the (units)
    
    if( FLG_avg_anyAngle_Scargle_FFT == 0 ):
        pAMPERE, = ax.plot( AMPERE_integrated_scargPeriod, AMPERE_integrated_scargPower ,color='xkcd:cerulean' ,linewidth=PLOT_lineWidthRegular, linestyle='--'); #plot
        ax.plot( AMPERE_integrated_scargPeriod, np.tile(AMPERE_integrated_scarggf,np.size(AMPERE_integrated_scargPeriod)) , color="xkcd:grey", linestyle='--' ); #plot
        # AMPERE_minPeriod = np.min(AMPERE_integrated_scargPeriod); #get the min period for the TEC
    else:
        pAMPERE, = ax.plot(1/freqs_AMPERE,(Cxx_AMPERE),color='xkcd:cerulean',linewidth=PLOT_lineWidthRegular, linestyle='--');
        # AMPERE_minPeriod = np.min(1/freqs_AMPERE); #get the min period for the TEC
    #END IF
    # plot_period_min = np.min( (TEC_minPeriod,AMPERE_minPeriod) ); #get the min period
    
    ax.set_xlabel("Periods (min)",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax.set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    string_Title = spectraName+' of delta-vTEC Avg''d at '+str(np.round(latMillstone,2))+' arcdeg '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' and Log of '+AMPERE_plot_label; #create mecha title
    ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    # if( np.all(np.isnan(AMPERE_integrated_scargPower)) == 0 ): #avoids an error
    #     ax.set_ylim( (0, np.max(AMPERE_integrated_scargPower)+0.1*np.max(AMPERE_integrated_scargPower) ) ); #set x axis limits
    # #END IF
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax.set_xticks(xAxisTicks); #set x axis ticks
    ax.set_xlim( (plot_period_min, plot_Period_Lim) ); #set x axis limits
    
    ax.legend([pTEC,pAMPERE],
        ['TEC',
         AMPERE_plot_label_noUnits,],
        loc='upper left');
    
    fig.subplots_adjust(left = 0.068, right = 0.935, top = 0.96, bottom = 0.075 , hspace = 0.205); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
    warnings.filterwarnings("default", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
       
#no return, plot is the return!




