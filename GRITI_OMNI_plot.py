"""
GOAL: Plot only OMNI
RD on 5/11/19

INPUT: buncha OMNI stuff
OUTPUT: no vars, just a plot is made

opt: 0 (default) - plots regular
    1 - plots with vertical lines dividing the days
    2 - plots with vertical lines dividing the days, and the dates written just above the bottom axis
    3 - does the same as above, but uses day number format instead of dates
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from subfun_date_to_dayNum import subfun_date_to_dayNum
from subfun_figFitter import figFitter

def GRITI_OMNI_plot(OMNI_data,OMNI_timeUnique,OMNI_dict,OMNI_dictPlot,OMNI_plotSet_name, \
        time_Ref,time_Reference,dateRange_full,dateRange_zeroHr,dateRange_dayNum_zeroHr,dateRange_zeroHr_monthName, \
        dateRange_zeroHr_dayPostfix,FONT_titleFM,FONT_axisLabelFM,opt=0,FLG_fancyPlot=0):

    #----OMNI can be longer than the dateRange due to time shift plotting, don't want that here----
    dateRange_dayNum_full = subfun_date_to_dayNum(dateRange_full); #get the full date range in daynum format
    k = (OMNI_timeUnique/86400 < dateRange_dayNum_full[0,1]) | (OMNI_timeUnique/86400 >= dateRange_dayNum_full[-1,1]+1); #get stuff outside the date range
    #!! NO YEAR SUPPORT !!
    OMNI_timeUnique = np.delete(OMNI_timeUnique,k); #delete em    
    OMNI_data = np.delete(OMNI_data,k,axis=0); #delete em

    #-----Plot OMNI results as several time series-----
    OMNI_timeUnique_hr = (OMNI_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
    
    if( np.mod(np.round(np.min(OMNI_timeUnique_hr)),2) == 0 ):
        OMNI_time_hr_axis_min = np.round(np.min(OMNI_timeUnique_hr)); #is even, good to go
    else:
        OMNI_time_hr_axis_min = np.round(np.min(OMNI_timeUnique_hr))+1; #is odd, make even
    #END IF
    if( np.mod(np.round(np.max(OMNI_timeUnique_hr)),2) == 0 ):
        OMNI_time_hr_axis_max = np.round(np.max(OMNI_timeUnique_hr)); #is even, good to go
    else:
        OMNI_time_hr_axis_max = np.round(np.max(OMNI_timeUnique_hr))-1; #is odd, make even
    #END IF
    
    xAxisTicks = np.arange(OMNI_time_hr_axis_min,OMNI_time_hr_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    
    if( opt >= 1 ):
        dateRange_dayNum_full = subfun_date_to_dayNum(dateRange_full); #call function to get the date range into dayNumber form (easy to work with)
    #END IF
    
    #Start the OMNI plot
    fig, ax = plt.subplots(nrows=len(OMNI_plotSet_name), ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #plot each one systematically
    for i in range(0,len(OMNI_plotSet_name)):
        ax[i].set_aspect('auto'); #Remove the aspect ratio from the basemap so it fills the screen better
        ax[i].plot( OMNI_timeUnique_hr, OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]] , linewidth=0.8 ); #plot
        
        if( (np.abs((np.min(OMNI_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.min(time_Ref))/3600 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'OMNI')) ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
            ax[i].plot( np.repeat( (np.min(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , 10) , np.linspace(np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),num=10), linewidth=1.75, color='r'); #plot red lines showing ISR data time
        if( (np.abs((np.max(OMNI_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.max(time_Ref))/3600 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'OMNI')) ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
            ax[i].plot( np.repeat( (np.max(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , 10) , np.linspace(np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),num=10), linewidth=1.75, color='r'); #plot red lines showing ISR data time
        #END IF
        
        ax[i].set_xticks(xAxisTicks); #set x axis ticks
        if( i != (len(OMNI_plotSet_name)-1)):
            ax[i].set_xticklabels([]); #if statement to remove x axis labels except for the last line
            ax[i].tick_params(axis="x",direction="in");
        #END IF
        ax[i].set_xlim( OMNI_time_hr_axis_min , OMNI_time_hr_axis_max ); #set y axis limits
        
        ax[i].set_ylabel(OMNI_dictPlot[OMNI_dict[OMNI_plotSet_name[i]]],fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax[i].set_ylim( np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]) , np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]) ); #set y axis limits
        
        ax[i].grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on
        
        if( opt >= 1 ):
            for j in range(0,dateRange_full.shape[0]):
                #run through each day, print the line
                if( j != 0 ):
                    ax[i].plot( np.repeat( (dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24, 10) , np.linspace(np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),num=10), linewidth=1.00, color='k', linestyle='--'); #plot black dotted lines showing the date separation
                #END IF
                if( (i == 0) & (opt >= 2) ): #only print on the top plot
                    if( opt == 2 ):
                        ax[i].text((dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-(np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]))*.15, str(dateRange_full[j,1])+'/'+str(dateRange_full[j,2])+'/'+str(dateRange_full[j,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
                    elif( opt == 3 ):
                        ax[i].text((dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-(np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]))*.15, str(dateRange_dayNum_full[j,1])+', '+str(dateRange_dayNum_full[j,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
                    #END IF
                #END IF
            #END FOR j
        #END IF
    
    #END FOR i
    
    #do things for the whole plot to finish up now
    string_Title = 'OMNI-Sourced Data (1 min resolution) for '+str(dateRange_full[0,1])+'/'+str(dateRange_full[0,2])+ \
        '/'+str(dateRange_full[0,0])+' to '+str(dateRange_full[-1,1])+ \
        '/'+str(dateRange_full[-1,2])+'/'+str(dateRange_full[-1,0])+ \
        ' (M/D/Y)'; #create mecha title
    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    ax[i].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    
    figFitter(fig); #fit that fit fast
    # fig.subplots_adjust(left = 0.050, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up