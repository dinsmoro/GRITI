"""
GOAL: Plot only Lomb-Scargle of ISR POPL HP at 300 km as well as ISR POPL and ISR POPL HP for comparsion and looking
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from subfun_lombscargle import subfun_lombscargle
from subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_ISR_Haystack_plot_POPL_scargleSet(Zenith_time,Zenith_POPL,Zenith_POPL_hp,Zenith_filtHeight,MISA_time,MISA_POPL,MISA_POPL_hp,MISA_filtHeight,pointAltitude,plot_Period_Lim,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):

    Zenith_scargPeriod, Zenith_scargPower, Zenith_scarggf = subfun_lombscargle((Zenith_time - dateRange_dayNum_zeroHr[1])*24 , Zenith_POPL_hp[Zenith_filtHeight,:]); #scargle that data
    Zenith_scargPeriod = Zenith_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
    
    MISA_scargPeriod, MISA_scargPower, MISA_scarggf = subfun_lombscargle((MISA_time - dateRange_dayNum_zeroHr[1])*24 , MISA_POPL_hp[MISA_filtHeight,:]); #scargle that data
    MISA_scargPeriod = MISA_scargPeriod*60; #min, adjust the period out from hrs to minutes (since hrs goes in)
    
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
    
    #-----ZENITH lomb-scargle of POPL high-passed at 300km-----
    ax[2,0].plot( Zenith_scargPeriod, Zenith_scargPower ); #plot
    ax[2,0].plot( Zenith_scargPeriod, np.tile(Zenith_scarggf,np.size(Zenith_scargPeriod)) , color="xkcd:grey" ); #plot
    
    ax[2,0].set_xlabel("Periods [min]",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2,0].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[2,0].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[2,0].set_xticks(xAxisTicks); #set x axis ticks
    if( np.all(np.isnan(Zenith_scargPower)) == 0 ): #avoids an error
        ax[2,0].set_ylim( (0, np.max(Zenith_scargPower)+0.1*np.max(Zenith_scargPower) ) ); #set x axis limits
    #END IF
    
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
    ax[2,1].plot( MISA_scargPeriod, MISA_scargPower ); #plot
    ax[2,1].plot( MISA_scargPeriod, np.tile(MISA_scarggf,np.size(MISA_scargPeriod)) , color="xkcd:grey" ); #plot
    
    ax[2,1].set_xlabel("Periods [min]",fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2,1].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[2,1].set_xlim( (0, plot_Period_Lim) ); #set x axis limits
    xAxisTicks = np.arange( 0, plot_Period_Lim+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[2,1].set_xticks(xAxisTicks); #set x axis ticks
    if( np.all(np.isnan(MISA_scargPower)) == 0 ): #avoids an error
        ax[2,1].set_ylim( (0, np.max(MISA_scargPower)+0.1*np.max(MISA_scargPower) ) ); #set x axis limits
    #END IF
    
    #final plot adjusting stuff
    fig.subplots_adjust(left = 0.065, right = 0.975, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up