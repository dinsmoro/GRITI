"""
GOAL: Plot only Lomb-Scargle of ISR POPL HP at 300 km as well as ISR POPL and ISR POPL HP for comparsion and looking
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
from scipy import signal
import matplotlib.pyplot as plt
from subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_ISR_Haystack_plot_POPL_FFTSet(Zenith_time,Zenith_POPL,Zenith_POPL_hp,Zenith_filtHeight,MISA_time,MISA_POPL,MISA_POPL_hp,MISA_filtHeight,pointAltitude,plot_Period_Lim,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):

    window = np.hamming(np.int64(Zenith_time.size/1.5));
    nooverlap = np.int64(Zenith_time.size/3);
    nfft = 2**12 #2^12=4096
    Fs_Z = 1/(np.median(np.diff(Zenith_time))*24*60); #min, zenith time delta in freq form
    Fs_M = 1/(np.median(np.diff(MISA_time))*24*60); #min, MISA time delta in freq form
    
    pwr_Z = np.sqrt(1/Zenith_POPL_hp[Zenith_filtHeight,:].size*np.sum(Zenith_POPL_hp[Zenith_filtHeight,:]**2)); #estimate power of signal
    pwr_M = np.sqrt(1/MISA_POPL_hp[MISA_filtHeight,:].size*np.sum(MISA_POPL_hp[MISA_filtHeight,:]**2)); #estimate power of signal

    [freqs_Z,Cxx_Z] = signal.welch(1/pwr_Z*Zenith_POPL_hp[Zenith_filtHeight,:] ,window=window,noverlap=nooverlap,nfft=nfft,fs=Fs_Z);
    [freqs_M,Cxx_M] = signal.welch(1/pwr_M*MISA_POPL_hp[MISA_filtHeight,:] ,window=window,noverlap=nooverlap,nfft=nfft,fs=Fs_M);    
    
    #-----Plot ISR POPL HP results as a RTI-----
    #Plot just the ISR POPL HP results
    fig, ax = plt.subplots(nrows=3, ncols=2); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0,0].set_aspect('auto');
    ax[0,1].set_aspect('auto');
    ax[1,0].set_aspect('auto');
    ax[1,1].set_aspect('auto');
    ax[2,0].set_aspect('auto');
    ax[2,1].set_aspect('auto');
    
    #~~~~~~~~~~ ZENITH STUFF ~~~~~~~~~~
    #-----ZENITH POPL at 300km-----
    ax[0,0].plot( (Zenith_time-dateRange_dayNum_zeroHr[1])*24, Zenith_POPL[Zenith_filtHeight,:] ); #plot
    
    string_Title = 'Zenith Beam - '+subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]); #create mecha title
    ax[0,0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[0,0].set_ylabel('POPL at '+str(pointAltitude)+ ' km',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0,0].set_xticks(xAxisTicks); #set x axis ticks
    ax[0,0].set_xlim( ((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    #-----ZENITH POPL high-passed at 300km-----
    ax[1,0].plot( (Zenith_time-dateRange_dayNum_zeroHr[1])*24, Zenith_POPL_hp[Zenith_filtHeight,:] ); #plot
    
    ax[1,0].set_ylabel('High-pass POPL at '+str(pointAltitude)+ ' km',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1,0].set_xticks(xAxisTicks); #set x axis ticks
    ax[1,0].set_xlim( ((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    #-----ZENITH spectra of POPL high-passed at 300km-----
    ax[2,0].plot(1/freqs_Z,Cxx_Z); # ,color='xkcd:deep red',linewidth=1.5, linestyle='--'
#    ax[2,0].plot( Zenith_scargPeriod, np.tile(Zenith_scarggf,np.size(Zenith_scargPeriod)) , color="xkcd:grey" ); #plot
    
    ax[2,0].set_xlabel("Periods [mi]",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2,0].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[2,0].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[2,0].set_xticks(xAxisTicks); #set x axis ticks
    
    #~~~~~~~~~~ MISA STUFF ~~~~~~~~~~
    #-----MISA POPL at 300km-----
    ax[0,1].plot( (MISA_time-dateRange_dayNum_zeroHr[1])*24, MISA_POPL[MISA_filtHeight,:] ); #plot
    
    string_Title = 'MISA Beam - Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' (hr)'; #create mecha title
    ax[0,1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[0,1].set_ylabel('POPL at '+str(pointAltitude)+ ' km',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0,1].set_xticks(xAxisTicks); #set x axis ticks
    ax[0,1].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    #-----MISA POPL high-passed at 300km-----
    ax[1,1].plot( (MISA_time-dateRange_dayNum_zeroHr[1])*24, MISA_POPL_hp[MISA_filtHeight,:] ); #plot
    
    ax[1,1].set_ylabel('High-pass POPL at '+str(pointAltitude)+ ' km',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1,1].set_xticks(xAxisTicks); #set x axis ticks
    ax[1,1].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    #-----MISA lomb-scargle of POPL high-passed at 300km-----
    ax[2,1].plot(1/freqs_M,Cxx_M); #,color='xkcd:goldenrod',linewidth=1.5, linestyle='-.'
#    ax[2,1].plot( MISA_scargPeriod, np.tile(MISA_scarggf,np.size(MISA_scargPeriod)) , color="xkcd:grey" ); #plot
    
    ax[2,1].set_xlabel("Periods [min]",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2,1].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[2,1].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[2,1].set_xticks(xAxisTicks); #set x axis ticks
    
    #final plot adjusting stuff
    fig.subplots_adjust(left = 0.065, right = 0.975, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up