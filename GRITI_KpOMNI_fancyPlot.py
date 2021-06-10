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
import matplotlib.gridspec as gridspec
from subfun_date_to_dayNum import subfun_date_to_dayNum
from subfun_figFitter import figFitter


def GRITI_KpOMNI_fancyPlot(Kp_data,Kp_time, OMNI_data,OMNI_timeUnique, \
        OMNI_dict, OMNI_dictPlot, OMNI_plotSet_name, time_Ref,time_Reference, \
        dateRange_full, dateRange_zeroHr, dateRange_dayNum_zeroHr, dateRange_zeroHr_monthName, \
        dateRange_zeroHr_dayPostfix, FONT_grandioseFM, FONT_titleFM, FONT_axisLabelFM, \
        PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi,opt=2):
    print('MAKING FANCY PLOT: KpOMNI_fancyPlot IN fancyPlot FOLDER'); #report since you won't see anything
    
    letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.90; #set the X position of the lettering (e.g., a. b. c. ...)

    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth['double plus']; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth['regular plus']; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth['regular']; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth['smol']; #get the line widths
    
    #----Kp can be longer than the dateRange due to time shift plotting, don't want that here----
    dateRange_dayNum_full = subfun_date_to_dayNum(dateRange_full); #get the full date range in daynum format
    k = np.unique(np.round(Kp_time/86400))[0:-1];
    k = (k < dateRange_dayNum_full[0,1]) | (k > dateRange_dayNum_full[-1,1]); #get stuff outside the date range
    Kp_data = np.delete(Kp_data,k,axis=0); #delete em
    k = (Kp_time/86400 <= dateRange_dayNum_full[0,1]) | (Kp_time/86400 > dateRange_dayNum_full[-1,1]+1); #get stuff outside the date range
    #!! NO YEAR SUPPORT !!
    Kp_time = np.delete(Kp_time,k); #delete em
    #----OMNI can be longer than the dateRange due to time shift plotting, don't want that here----
    k = (OMNI_timeUnique/86400 < dateRange_dayNum_full[0,1]) | (OMNI_timeUnique/86400 >= dateRange_dayNum_full[-1,1]+1); #get stuff outside the date range
    #!! NO YEAR SUPPORT !!
    OMNI_timeUnique = np.delete(OMNI_timeUnique,k); #delete em    
    OMNI_data = np.delete(OMNI_data,k,axis=0); #delete em
    
    # if( opt >= 1 ):
    #     dateRange_dayNum_full = subfun_date_to_dayNum(dateRange_full); #call function to get the date range into dayNumber form (easy to work with)
    # #END IF
    
    #-----PREP Kp STUFF-----
    Kp_data_plot = np.repeat(Kp_data.flatten(),2); #replicate teh Kp for plotting
    Kp_time_plot = (np.repeat(Kp_time,2) - dateRange_dayNum_zeroHr[1]*86400)/3600; #hrs at 0 hr, replicate the hours for plotting purposes
    Kp_time_plot = np.hstack( (np.min(Kp_time_plot)-3,Kp_time_plot[0:len(Kp_time_plot)-1]) ); #readjust so the hour range matches what is real (e.g. Kp lasts 3 hr, so 0-3 hr is same Kp value)
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
    xAxisTicksKp = np.arange(Kp_time_axis_min,Kp_time_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    
    #-----PREP OMNI STUFF-----
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
    xAxisTicksMax = 25; #maximum number of ticks that can fit on the bottom
    xAxisTicksSteps = np.array( (1,2,4,6,12,24) ); #allowed steps to choose from
    xAxisTicksNum = (OMNI_time_hr_axis_max-OMNI_time_hr_axis_min)/xAxisTicksSteps+1; #calc how many tick marks there will be
    xAxisTicksStep = np.where(xAxisTicksNum <= xAxisTicksMax)[0][0]; #get the step that is good enough
    xAxisTicksOMNI = np.arange(OMNI_time_hr_axis_min,OMNI_time_hr_axis_max+xAxisTicksSteps[xAxisTicksStep],xAxisTicksSteps[xAxisTicksStep]); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    
    #-----PLOT THE Kp AND OMNI STUFF-----
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig = plt.figure(figsize=(14,14.5),dpi=journal_dpi);
    gridr = gridspec.GridSpec(nrows=len(OMNI_plotSet_name)+1, ncols=1, figure=fig);
    gridr.update(hspace=0.05); # set the spacing between axes. 
    fig.add_subplot(gridr[0:1]);
    for i in range(1,len(OMNI_plotSet_name)+1):
        fig.add_subplot(gridr[i]); #dynamically add axes
    #END FOR i
    ax = fig.axes; #get a list of the axes
    
    #Start the Kp portion
    #Remove the aspect ratio so it fills the screen better
    ax[0].set_aspect('auto'); #set to auto for all axes
    ax[0].plot( Kp_time_plot, Kp_data_plot,linewidth=PLOT_lineWidthDoublePlus, antialiased=True); #plot
    
    # string_Title = 'Kp Index for '+str(dateRange_full[0,1])+'/'+str(dateRange_full[0,2])+ \
    #     '/'+str(dateRange_full[0,0])+' to '+str(dateRange_full[-1,1])+ \
    #     '/'+str(dateRange_full[-1,2])+'/'+str(dateRange_full[-1,0])+ \
    #     ' (M/D/Y)'; #create mecha title
    # ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    # ax.set_xlabel('Time in UT (hr) - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('Kp Index',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[0].set_xticks(xAxisTicksKp); #set x axis ticks
    ax[0].set_xlim( Kp_time_axis_min , Kp_time_axis_max ); #set y axis limits
    ax[0].set_xticklabels([]); #if statement to remove x axis labels except for the last line
    ax[0].tick_params(axis="x",direction="in");
    ax[0].grid(b=True, which='major', axis='both', color='xkcd:light grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
    ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text labelling the lettering a. b. c. ect.
    
#    yAxisTicks = ; #creates y ticks automagically
#    ax.set_yticks(yAxisTicks); #set x axis ticks
    ax[0].set_ylim( 0 , np.max(Kp_data)+0.5 ); #set y axis limits
    
    if( (np.abs((np.min(Kp_time_plot)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.min(time_Ref))/3600 >= 0.25) & (time_Reference != 'Kp') ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
        ax[0].plot( np.repeat( (np.min(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , np.arange(0,(np.max(Kp_data)+1),0.5).size) , np.arange(0,(np.max(Kp_data)+1),0.5), linewidth=PLOT_lineWidthPlus, color='r', antialiased=True); #plot red lines showing ISR data time
    if( (np.abs((np.max(Kp_time_plot)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.max(time_Ref))/3600 >= 0.25) & (time_Reference != 'Kp') ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
        ax[0].plot( np.repeat( (np.max(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , np.arange(0,(np.max(Kp_data)+1),0.5).size) , np.arange(0,(np.max(Kp_data)+1),0.5), linewidth=PLOT_lineWidthPlus, color='r', antialiased=True); #plot red lines showing ISR data time
    #END IF
    if( opt >= 1 ):
        for i in range(0,dateRange_full.shape[0]):
            #run through each day, print the line
            if( i != 0 ):
                ax[0].plot( np.repeat( (dateRange_dayNum_full[i,1] - dateRange_dayNum_zeroHr[1])*24, np.arange(0,(np.max(Kp_data)+1),0.5).size) , np.arange(0,(np.max(Kp_data)+1),0.5), linewidth=PLOT_lineWidthRegularPlus, color='k', linestyle='--', antialiased=True); #plot black dotted lines showing the date separation
            #END IF
            if( opt == 2 ):
                ax[0].text((dateRange_dayNum_full[i,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(Kp_data)+0.05, str(dateRange_full[i,1])+'/'+str(dateRange_full[i,2])+'/'+str(dateRange_full[i,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
            elif( opt == 3 ):
                ax[0].text((dateRange_dayNum_full[i,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(Kp_data)+0.05, str(dateRange_dayNum_full[i,1])+', '+str(dateRange_dayNum_full[i,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
            #END IF
        #END FOR i
    #END IF
    
    #Start the OMNI portion
    #plot each one systematically
    for i in range(0,len(OMNI_plotSet_name)):
        ax[i+1].set_aspect('auto'); #Remove the aspect ratio from the basemap so it fills the screen better
        ax[i+1].plot( OMNI_timeUnique_hr, OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]] , linewidth=PLOT_lineWidthSmol, antialiased=True); #plot
        
        if( (np.abs((np.min(OMNI_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.min(time_Ref))*24 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'OMNI')) ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
            ax[i+1].plot( np.repeat( (np.min(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , 10) , np.linspace(np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),num=10), linewidth=PLOT_lineWidthPlus, color='r', antialiased=True); #plot red lines showing ISR data time
        if( (np.abs((np.max(OMNI_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.max(time_Ref))*24 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'OMNI')) ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
            ax[i+1].plot( np.repeat( (np.max(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , 10) , np.linspace(np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),num=10), linewidth=PLOT_lineWidthPlus, color='r', antialiased=True); #plot red lines showing ISR data time
        #END IF
        
        ax[i+1].set_xticks(xAxisTicksOMNI); #set x axis ticks
        if( i != (len(OMNI_plotSet_name)-1)):
            ax[i+1].set_xticklabels([]); #if statement to remove x axis labels except for the last line
            ax[i+1].tick_params(axis="x",direction="in");
        #END IF
        ax[i+1].set_xlim( OMNI_time_hr_axis_min , OMNI_time_hr_axis_max ); #set y axis limits
        
        ax[i+1].set_ylabel(OMNI_dictPlot[OMNI_dict[OMNI_plotSet_name[i]]],fontproperties=FONT_axisLabelFM); #set the y axis label
        
        ax[i+1].set_ylim( np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]) , np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]) ); #set y axis limits
        
        ax[i+1].grid(b=True, which='major', axis='both', color='xkcd:light grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
        
        if( 98+i < 123): #avoids weird stuff in a weird case
            if( (i == 0) | (i == 1) ):
                #adjusted y positioning for b. and c.
                ax[i+1].text( letteringPositionX, letteringPositionY+.05, chr(98+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i+1].transAxes); #print the text labelling the lettering a. b. c. ect.
                #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
            else:
                ax[i+1].text( letteringPositionX, letteringPositionY, chr(98+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i+1].transAxes); #print the text labelling the lettering a. b. c. ect.
                #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
            #END IF
        #END IF
        
        if( opt >= 1 ):
            for j in range(0,dateRange_full.shape[0]):
                #run through each day, print the line
                if( j != 0 ):
                    ax[i+1].plot( np.repeat( (dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24, 10) , np.linspace(np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]),num=10), linewidth=PLOT_lineWidthRegularPlus, color='k', linestyle='--', antialiased=True); #plot black dotted lines showing the date separation
                #END IF
                # if( (i == 0) & (opt >= 2) ): #only print on the top plot
                #     if( opt == 2 ):
                #         ax[i+1].text((dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-(np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]))*.15, str(dateRange_full[j,1])+'/'+str(dateRange_full[j,2])+'/'+str(dateRange_full[j,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
                #     elif( opt == 3 ):
                #         ax[i+1].text((dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-(np.max(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]])-np.min(OMNI_data[:,OMNI_dict[OMNI_plotSet_name[i]]]))*.15, str(dateRange_dayNum_full[j,1])+', '+str(dateRange_dayNum_full[j,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
                #     #END IF
                # #END IF
            #END FOR j
        #END IF
    
    #END FOR i
    
    #do things for the whole plot to finish up now
    # string_Title = 'OMNI-Sourced Data (1 min resolution) for '+str(dateRange_full[0,1])+'/'+str(dateRange_full[0,2])+ \
    #     '/'+str(dateRange_full[0,0])+' to '+str(dateRange_full[-1,1])+ \
    #     '/'+str(dateRange_full[-1,2])+'/'+str(dateRange_full[-1,0])+ \
    #     ' (M/D/Y)'; #create mecha title
    # ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    ax[i+1].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label

    figFitter(fig); #fit that fig fast
    # fig.subplots_adjust(left = 0.10, right = 0.982, top = 0.990, bottom = 0.050); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig('fancyPlots\Kp&OMNI.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff