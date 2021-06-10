#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

def GRITI_combinedPlot_keo_TEC_n_AMPERE_1Dintegration(sTECChunked_anyAngleAvg,TEC_timeUnique,TEC_plotLimValu,colorMap, \
        AMPERE_data,AMPERE_timeUnique,locAMPERE_time,locAMPERE_lat,locAMPERE_long,AMPERE_plot_index,AMPERE_plot_indexes,AMPERE_plot_labels, \
        FLG_AMPERE_upTo90,plotLatRange,plotLongRange,latMillstone,longMillstone,dateRange_dayNum_zeroHr,avg_anyAngle,avg_anyAngle_Width, \
        avg_anyAngle_Range_Chunks_Long_Plot,avg_anyAngle_Range_Chunks_Long_Plot_Name,avg_anyAngle_dataType,avg_anyAngle_plotLabel, \
        time_Ref,time_Reference,dateRange,dateRange_zeroHr,dateRange_zeroHr_monthName,dateRange_zeroHr_dayPostfix,time_cutout_range_delay, \
        PLOT_lineWidth,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):

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
    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    divider = make_axes_locatable(ax[0]); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    divider2 = make_axes_locatable(ax[1]); #prep to add an axis
    cax2 = divider2.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    cax2.set_visible(False); #mkae it invisible so it matches the other plots in width
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    
    #[(TEC_timeUnique - dateRange_dayNum_zeroHr[1])*24 , avg_anyAngle_Range_Chunks_Long_Plot] ,
    pltHelprX, pltHelprY = np.meshgrid( (TEC_timeUnique - dateRange_dayNum_zeroHr[1])*24, \
                avg_anyAngle_Range_Chunks_Long_Plot);
    im = ax[0].pcolormesh(pltHelprX, pltHelprY,  sTECChunked_anyAngleAvg.T ,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu),cmap=colorMap); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cax.yaxis.set_ticks(np.linspace(np.min(TEC_plotLimValu),np.max(TEC_plotLimValu),5)); #create useful tick marks
    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
    cbar.set_label(avg_anyAngle_plotLabel); #tabel the colorbar
    cbar.ax.tick_params(labelsize=FONT_axisTick);
    #cbar.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu)); #they changed how the code works, this doesn't work anymore
    cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu)); #now it's this
    cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle*180/np.pi,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+avg_anyAngle_Range_Chunks_Long_Plot_Name+' of Millstone Hill Zenith Beam'; #create mecha title
    string_Title = avg_anyAngle_dataType+' Averaged on Angle of '+str(np.round(avg_anyAngle,2))+' deg and Width of '+ \
        str(np.round(avg_anyAngle_Width,2))+' arcdeg'; #create mecha title
    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    # ax[0].set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' [hr]',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel(avg_anyAngle_Range_Chunks_Long_Plot_Name+' [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(TEC_timeUnique)-dateRange_dayNum_zeroHr[1])*24),4)) , \
            (np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(TEC_timeUnique)-dateRange_dayNum_zeroHr[1])*24),4)) + 4 , \
            4); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    
    avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot)) - np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 ):
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
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLatRange) - np.min(plotLatRange))/13; #just goes for it if it's a super tiny range
        elif(avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude'): #if Y axis is longitude, use longitude
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (np.max(plotLongRange) - np.min(plotLongRange))/13; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(avg_anyAngle_Range_Chunks_Long_Plot)),np.ceil(np.max(avg_anyAngle_Range_Chunks_Long_Plot))+avg_anyAngle_Range_Chunks_Long_Plot_autoTick,avg_anyAngle_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
    ax[0].set_yticks(yAxisTicks); #set x axis ticks
    
    #Now drawing line of interest
    if( avg_anyAngle_Range_Chunks_Long_Plot_Name == 'Longitude' ): #if true, longitude
        if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
            ax[0].plot( np.linspace(np.min((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),np.max((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),10,endpoint=True) , #X time hr
                    np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=1); #plots a point with a black line
        #END IF
    else: #else latitude
        if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
            ax[0].plot( np.linspace(np.min((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),np.max((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),10,endpoint=True) , #X time hr
                    np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=1); #plots a point with a black line
        #END IF
    #END IF
    ax[0].set_xlim([np.min(xAxisTicks),np.max(xAxisTicks)] ); #force axis to reach ends (if it plots 47.99995 it won't put a 48)
    ax[0].set_ylim([np.min(yAxisTicks),np.max(yAxisTicks)] ); #force axis to reach ends (if it plots 47.99995 it won't put a 48)
    
    
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
    
    #prep to plot
    if( FLG_AMPERE_upTo90 == 1 ):
        plotLatRange = [np.min(plotLatRange),90]; #force up to 90 for averaging
    #END IF

    
    kInRange = (AMPERE_data[:,locAMPERE_lat] >= np.min(plotLatRange)) & (AMPERE_data[:,locAMPERE_lat] <= np.max(plotLatRange)) & \
        (AMPERE_data[:,locAMPERE_long] >= np.min(plotLongRange)) & (AMPERE_data[:,locAMPERE_long] <= np.max(plotLongRange)); #get data in the range
    
    AMPERE_jouleHeating_integrate = np.zeros( AMPERE_timeUnique_hr.size , dtype=np.float64); #prep integrated joule heating
    for i in range(AMPERE_timeUnique_hr.size):
        k = AMPERE_timeUnique[i] == AMPERE_data[:,locAMPERE_time]; #get the right time
        AMPERE_jouleHeating_integrate[i] = np.sum((AMPERE_data[k&kInRange,AMPERE_plot_index])); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
    #END FOR i
    
    AMPERE_plot_label = AMPERE_plot_labels[np.where(AMPERE_plot_indexes == AMPERE_plot_index)[0][0]]; #get the label
    AMPERE_plot_label_noUnits = AMPERE_plot_label[0:AMPERE_plot_label.find('[')-1]; #remove the (units)
            
    ax[1].plot( AMPERE_timeUnique_hr + time_cutout_range_delay, AMPERE_jouleHeating_integrate , linewidth=PLOT_lineWidthRegularPlus ); #plot
    
    if( (np.abs((np.min(AMPERE_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1]) - np.min(time_Ref))*24 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'AMPERE')) ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
        ax[1].plot( np.repeat( (np.min(time_Ref) - dateRange_dayNum_zeroHr[1])*24 , 10) , np.linspace(np.min(AMPERE_jouleHeating_integrate),np.max(AMPERE_jouleHeating_integrate),num=10), linewidth=1.75, color='r'); #plot red lines showing ISR data time
    if( (np.abs((np.max(AMPERE_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1]) - np.max(time_Ref))*24 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'AMPERE')) ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
        ax[1].plot( np.repeat( (np.max(time_Ref) - dateRange_dayNum_zeroHr[1])*24 , 10) , np.linspace(np.min(AMPERE_jouleHeating_integrate),np.max(AMPERE_jouleHeating_integrate),num=10), linewidth=1.75, color='r'); #plot red lines showing ISR data time
    #END IF
    
    xAxisTicks = np.arange(AMPERE_time_hr_axis_min,AMPERE_time_hr_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks

    ax[1].set_xlim( AMPERE_time_hr_axis_min , AMPERE_time_hr_axis_max ); #set y axis limits
    
    ax[1].set_ylabel(AMPERE_plot_label,fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax[1].set_ylim( np.min(AMPERE_jouleHeating_integrate) , np.max(AMPERE_jouleHeating_integrate) ); #set y axis limits
    
    ax[1].grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on        
    
    string_Title = 'Integrated AMPERE '+AMPERE_plot_label_noUnits+' in the Nothern Hemisphere ['+str(time_cutout_range_delay)+' hr delay]'; #create mecha title
    ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    
    ax[1].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    
    fig.subplots_adjust(left = 0.068, right = 0.935, top = 0.96, bottom = 0.075 , hspace = 0.205); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
       
#no return, plot is the return!




