"""
GOAL: Plot only lomb-scargle of an OMNI data type
RD on 5/11/19

INPUT: buncha OMNI stuff
option: OMNI_plot_scargle_highpassOption, chooses if 0 no high-pass on data, 1 both high-pass and no high-pass on data, or 2 high-pass on data only (default 1)
option: plot_Period_Lim, (minutes) period to limit plot by (default 120)
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from subfun_lombscargle import subfun_lombscargle
from subfun_date_to_dayNum import subfun_date_to_dayNum

def GRITI_OMNI_plot_scargle(OMNI_data,OMNI_timeUnique,OMNI_dict,OMNI_dictPlot,OMNI_plot_scargle_name,filter_cutoffPeriod,dateRange,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisLabelFM,settings,plot_Period_Lim=120,OMNI_plot_scargle_highpassOption=1):
    #unpack settings
    settings_spectra = settings['spectra'];    
    
    #----OMNI can be longer than the dateRange due to time shift plotting, don't want that here----
    dateRange_dayNum = subfun_date_to_dayNum(dateRange); #get the full date range in daynum format
    k = (OMNI_timeUnique < dateRange_dayNum[0,1]) | (OMNI_timeUnique >= dateRange_dayNum[-1,1]+1); #get stuff outside the date range
    #!! NO YEAR SUPPORT !!
    OMNI_timeUnique = np.delete(OMNI_timeUnique,k); #delete em    
    OMNI_data = np.delete(OMNI_data,k,axis=0); #delete em
    
    if( OMNI_plot_scargle_highpassOption == 0 ): #only original data, no high-passed data
    
        #-----Plot an OMNI data index's spectrum-----
        OMNI_timeUnique_hr = (OMNI_timeUnique - dateRange_dayNum_zeroHr[1])*24; #hr, convert to hr with 0 hr at specified day
        
        OMNI_data_scargPeriod, OMNI_data_scargPower, OMNI_data_scarggf = subfun_lombscargle(OMNI_timeUnique_hr , OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]]); #scargle that data
        OMNI_data_scargPeriod = OMNI_data_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        
        #OMNI plot prep
        OMNI_plot_scargle_label = OMNI_dictPlot[OMNI_dict[OMNI_plot_scargle_name]]; #get the label
        OMNI_plot_scargle_labelNoUnits = OMNI_plot_scargle_label[0:OMNI_plot_scargle_label.find('[')-1]; #remove the (units)
        
        #Start the OMNI scargle plot
        fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        
        ax.plot( OMNI_data_scargPeriod, OMNI_data_scargPower ); #plot
        ax.plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
        ax.set_xlabel("Periods (min)"+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
            '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax.set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax.set_xlim( (0, plot_Period_Lim) ); #set x axis limits
        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
            ax.set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
        else:
            ax.set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
        #END IF 
        
        string_Title = OMNI_plot_scargle_labelNoUnits+' - Lomb-Scargle Periodogram'; #create mecha title
        ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
        
        fig.subplots_adjust(left = 0.050, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
        plt.show(); #req to make plot show up    
        
    elif( OMNI_plot_scargle_highpassOption == 1 ): #combines original data scargle and high-passed data scargle
        from subfun_highpass import subfun_highpass
        from subfun_lowpass import subfun_lowpass
        
        #-----Plot an OMNI data index's spectrum-----
        OMNI_timeUnique_hr = (OMNI_timeUnique - dateRange_dayNum_zeroHr[1])*24; #hr, convert to hr with 0 hr at specified day
        
        #Start the OMNI scargle plot
        fig, ax = plt.subplots(nrows=2,ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        
        #~~~~~HIGH-PASSED ONE~~~~~
        OMNI_data_hp = subfun_highpass(OMNI_timeUnique_hr,OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]],filter_cutoffPeriod=filter_cutoffPeriod); #high-pass that data
        OMNI_data_hp = subfun_lowpass(OMNI_timeUnique_hr,OMNI_data_hp,filter_cutoffPeriod=25/60); #low-pass that data
        
        OMNI_data_scargPeriod, OMNI_data_scargPower, OMNI_data_scarggf = subfun_lombscargle(OMNI_timeUnique_hr , OMNI_data_hp); #scargle that data
        OMNI_data_scargPeriod = OMNI_data_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        
        #OMNI plot prep
        OMNI_plot_scargle_label = OMNI_dictPlot[OMNI_dict[OMNI_plot_scargle_name]]; #get the label
        OMNI_plot_scargle_labelNoUnits = OMNI_plot_scargle_label[0:OMNI_plot_scargle_label.find('[')-1]; #remove the (units)

        ax[0].plot( OMNI_data_scargPeriod, OMNI_data_scargPower ); #plot
        ax[0].plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
#        ax[0].set_xlabel("Periods (min)"+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
#            '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax[0].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax[0].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
            ax[0].set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
        else:
            ax[0].set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
        #END IF         
        string_Title = OMNI_plot_scargle_labelNoUnits+' w/ High-pass Cutoff Period of '+str(filter_cutoffPeriod)+' hrs - Lomb-Scargle Periodogram'; #create mecha title
        ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
        
        #~~~~~NO FILTER ONE~~~~~
        OMNI_data_scargPeriod, OMNI_data_scargPower, OMNI_data_scarggf = subfun_lombscargle(OMNI_timeUnique_hr , OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]]); #scargle that data
        OMNI_data_scargPeriod = OMNI_data_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        
        ax[1].plot( OMNI_data_scargPeriod, OMNI_data_scargPower ); #plot
        ax[1].plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
        ax[1].set_xlabel("Periods (min)"+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
            '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax[1].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax[1].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
            ax[1].set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
        else:
            ax[1].set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
        #END IF         
        string_Title = OMNI_plot_scargle_labelNoUnits+' - Lomb-Scargle Periodogram'; #create mecha title
        ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
        
        fig.subplots_adjust(left = 0.050, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
        plt.show(); #req to make plot show up  
        
    elif( OMNI_plot_scargle_highpassOption == 2 ): #only high-passed data
        from subfun_highpass import subfun_highpass
        
        #-----Plot an OMNI data index's spectrum-----
        OMNI_timeUnique_hr = (OMNI_timeUnique - dateRange_dayNum_zeroHr[1])*24; #hr, convert to hr with 0 hr at specified day
        
        OMNI_data_hp = subfun_highpass(OMNI_timeUnique_hr,OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]],filter_cutoffPeriod=filter_cutoffPeriod); #high-pass that data
        
        OMNI_data_scargPeriod, OMNI_data_scargPower, OMNI_data_scarggf = subfun_lombscargle(OMNI_timeUnique_hr , OMNI_data_hp); #scargle that data
        OMNI_data_scargPeriod = OMNI_data_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        
        #OMNI plot prep
        OMNI_plot_scargle_label = OMNI_dictPlot[OMNI_dict[OMNI_plot_scargle_name]]; #get the label
        OMNI_plot_scargle_labelNoUnits = OMNI_plot_scargle_label[0:OMNI_plot_scargle_label.find('[')-1]; #remove the (units)
        
        #Start the OMNI scargle plot
        fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        
        ax.plot( OMNI_data_scargPeriod, OMNI_data_scargPower ); #plot
        ax.plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
        ax.set_xlabel("Periods (min)"+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
            '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax.set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax.set_xlim( (0, plot_Period_Lim) ); #set x axis limits
        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
            ax.set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
        else:
            ax.set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
        #END IF         
        string_Title = OMNI_plot_scargle_labelNoUnits+' w/ High-pass Cutoff Period of '+str(filter_cutoffPeriod)+' hrs - Lomb-Scargle Periodogram'; #create mecha title
        ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
        
        fig.subplots_adjust(left = 0.050, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
        plt.show(); #req to make plot show up    

    elif( OMNI_plot_scargle_highpassOption == 3 ): #combines original data scargle and delta data scargle
        from scipy.signal import savgol_filter
        
        #-----Plot an OMNI data index's spectrum-----
        OMNI_timeUnique_hr = (OMNI_timeUnique - dateRange_dayNum_zeroHr[1])*24; #hr, convert to hr with 0 hr at specified day
        
        #Start the OMNI scargle plot
        fig, ax = plt.subplots(nrows=2,ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        
        #~~~~~DELTA-DATA ONE~~~~~
        #prep the Sav-Gol filter for debiasing
        windowLen_savGol = np.int64(np.round(settings_spectra['savgol filter period']/(np.median(np.diff(OMNI_timeUnique))*86400))); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
        #from conversations with AC ^
        if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
            windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
        #END IF
        polyYvals = savgol_filter(OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]],windowLen_savGol, settings_spectra['savgol filter order'] ); #filter it up
        OMNI_data_delta = OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]] - polyYvals;
        
        OMNI_data_scargPeriod, OMNI_data_scargPower, OMNI_data_scarggf = subfun_lombscargle(OMNI_timeUnique_hr , OMNI_data_delta); #scargle that data
        OMNI_data_scargPeriod = OMNI_data_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        
        #OMNI plot prep
        OMNI_plot_scargle_label = OMNI_dictPlot[OMNI_dict[OMNI_plot_scargle_name]]; #get the label
        OMNI_plot_scargle_labelNoUnits = OMNI_plot_scargle_label[0:OMNI_plot_scargle_label.find('[')-1]; #remove the (units)

        ax[0].plot( OMNI_data_scargPeriod, OMNI_data_scargPower ); #plot
        ax[0].plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
#        ax[0].set_xlabel("Periods (min)"+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
#            '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax[0].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax[0].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
            ax[0].set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
        else:
            ax[0].set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
        #END IF         
        string_Title = OMNI_plot_scargle_labelNoUnits+' w/ Sav-Gol Delta Order '+str(settings_spectra['savgol filter order'])+' and Period of '+str(settings_spectra['savgol filter period'])+' min - Lomb-Scargle Periodogram'; #create mecha title
        ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
        
        #~~~~~NO FILTER ONE~~~~~
        OMNI_data_scargPeriod, OMNI_data_scargPower, OMNI_data_scarggf = subfun_lombscargle(OMNI_timeUnique_hr , OMNI_data[:,OMNI_dict[OMNI_plot_scargle_name]]); #scargle that data
        OMNI_data_scargPeriod = OMNI_data_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
        
        ax[1].plot( OMNI_data_scargPeriod, OMNI_data_scargPower ); #plot
        ax[1].plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
        ax[1].set_xlabel("Periods (min)"+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
            '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax[1].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax[1].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
            ax[1].set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
        else:
            ax[1].set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
        #END IF         
        string_Title = OMNI_plot_scargle_labelNoUnits+' - Lomb-Scargle Periodogram'; #create mecha title
        ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
        
        fig.subplots_adjust(left = 0.050, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
        plt.show(); #req to make plot show up     
    #END IF