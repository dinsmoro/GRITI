"""
GOAL: Plot only Kp
RD on 5/11/19

INPUT: buncha Kp stuff
OUTPUT: no vars, just a plot is made

opt: 0 (default) - plots regular
    1 - plots with vertical lines dividing the days
    2 - plots with vertical lines dividing the days, and the dates written just above the bottom axis
    3 - does the same as above, but uses day number format instead of dates
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_figFitter import figFitter

def GRITI_Kp_plot(Kp_data,Kp_time,time_Ref,time_Reference, \
        dateRange_full,dateRange_zeroHr,dateRange_dayNum_zeroHr,dateRange_zeroHr_monthName,dateRange_zeroHr_dayPostfix, \
        FONT_titleFM,FONT_axisLabelFM,opt=0,FLG_fancyPlot=0):
    
    #----Kp can be longer than the dateRange due to time shift plotting, don't want that here----
    dateRange_dayNum_full = subfun_date_to_dayNum(dateRange_full); #get the full date range in daynum format
    k = np.unique(np.round(Kp_time/86400))[0:-1];
    k = (k < dateRange_dayNum_full[0,1]) | (k > dateRange_dayNum_full[-1,1]); #get stuff outside the date range
    Kp_data = np.delete(Kp_data,k,axis=0); #delete em
    k = (Kp_time/86400 <= dateRange_dayNum_full[0,1]) | (Kp_time/86400 > dateRange_dayNum_full[-1,1]+1); #get stuff outside the date range
    #!! NO YEAR SUPPORT !!
    Kp_time = np.delete(Kp_time,k); #delete em
    
    #-----Plot Kp results as a time series-----
    Kp_data_plot = np.repeat(Kp_data.flatten(),2); #replicate teh Kp for plotting
    Kp_time_plot = (np.repeat(Kp_time,2) - dateRange_dayNum_zeroHr[1]*86400)/3600; #hrs at 0 hr, replicate the hours for plotting purposes
    Kp_time_plot = np.hstack( (np.min(Kp_time_plot)-3,Kp_time_plot[0:len(Kp_time_plot)-1]) ); #readjust so the hour range matches what is real (e.g. Kp lasts 3 hr, so 0-3 hr is same Kp value)
    
    #Plot just the Kp results
    fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    ax.plot( Kp_time_plot, Kp_data_plot,linewidth=3 ); #plot
    
    string_Title = 'Kp Index for '+str(dateRange_full[0,1])+'/'+str(dateRange_full[0,2])+ \
        '/'+str(dateRange_full[0,0])+' to '+str(dateRange_full[-1,1])+ \
        '/'+str(dateRange_full[-1,2])+'/'+str(dateRange_full[-1,0])+ \
        ' (M/D/Y)'; #create mecha title
    ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax.set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax.set_ylabel('Kp Index, 3 Hr Resolution',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    if( np.mod(np.min(Kp_time_plot),2) == 0 ):
        Kp_time_axis_min = np.min(Kp_time_plot); #is even, good to go
    else:
        Kp_time_axis_min = np.min(Kp_time_plot)+1; #is odd, make even
    #END IF
    
    if( np.mod(np.max(Kp_time_plot),2) == 0 ):
        Kp_time_axis_max = np.max(Kp_time_plot); #is even, good to go
    else:
        Kp_time_axis_max = np.max(Kp_time_plot)-1; #is odd, make even
    #END IF
    xAxisTicks = np.arange(Kp_time_axis_min,Kp_time_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    ax.set_xticks(xAxisTicks); #set x axis ticks
    ax.set_xlim( Kp_time_axis_min , Kp_time_axis_max ); #set y axis limits
    
    yAxisTicks = np.arange(0,np.max(Kp_data)+1/3,1/3); #creates y ticks automagically [Kp goes in 1/3 steps]
    ax.set_yticks(yAxisTicks); #set x axis ticks
    ax.set_ylim( 0 , np.max(Kp_data)+1/3 ); #set y axis limits
    
    ytickLabelsNew = []; #prep list
    for i in range(0,yAxisTicks.size):
        labelInt = np.int64(yAxisTicks[i]); #get the integer portion
        labelFrac = yAxisTicks[i] - labelInt; #get fractional portion  $\frac{a}{b}$
        labelText = '';
        if( (labelInt != 0) | ((labelInt == 0) & np.isclose(labelFrac,0)) ):
            labelText += str(labelInt)+' '; #add in the integer
        #END IF
        if( ~np.isclose(labelFrac,0) ):
            labelText += '$\\frac{'+str(np.int64(np.round(labelFrac*3)))+'}{3}$'; #add in the fraction
        #END IF
        labelText = labelText.rstrip(' '); #remove extra spaces if there are any
        # ytickLabels[i].set_text(labelText); #set the text
        ytickLabelsNew.append(labelText); #append
    #END FOR i
    ax.set_yticklabels(ytickLabelsNew); #set the labels
    
    ax.grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on
    
    if( (np.abs((np.min(Kp_time_plot)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.min(time_Ref))/3600 >= 0.25) & (time_Reference != 'Kp') ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
        ax.plot( np.repeat( (np.min(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , np.arange(0,(np.max(Kp_data)+1),0.5).size) , np.arange(0,(np.max(Kp_data)+1),0.5), linewidth=1.75, color='r'); #plot red lines showing ISR data time
    if( (np.abs((np.max(Kp_time_plot)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.max(time_Ref))/3600 >= 0.25) & (time_Reference != 'Kp') ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
        ax.plot( np.repeat( (np.max(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , np.arange(0,(np.max(Kp_data)+1),0.5).size) , np.arange(0,(np.max(Kp_data)+1),0.5), linewidth=1.75, color='r'); #plot red lines showing ISR data time
    #END IF
    
    if( opt >= 1 ):
        dateRange_dayNum_full = subfun_date_to_dayNum(dateRange_full); #call function to get the date range into dayNumber form (easy to work with)
        for i in range(0,dateRange_full.shape[0]):
            #run through each day, print the line
            if( i != 0 ):
                ax.plot( np.repeat( (dateRange_dayNum_full[i,1] - dateRange_dayNum_zeroHr[1])*24, np.arange(0,(np.max(Kp_data)+1),0.5).size) , np.arange(0,(np.max(Kp_data)+1),0.5), linewidth=1.00, color='k', linestyle='--'); #plot black dotted lines showing the date separation
            #END IF
            if( opt == 2 ):
                ax.text((dateRange_dayNum_full[i,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(Kp_data)+0.20, str(dateRange_full[i,1])+'/'+str(dateRange_full[i,2])+'/'+str(dateRange_full[i,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
            elif( opt == 3 ):
                ax.text((dateRange_dayNum_full[i,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(Kp_data)+0.20, str(dateRange_dayNum_full[i,1])+', '+str(dateRange_dayNum_full[i,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
            #END IF
        #END FOR i
    #END IF
    
    figFitter(fig); #fit that fig fast
    # fig.subplots_adjust(left = 0.045, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up