"""
GOAL: Plot only SuperMAG
RD on 5/11/19

INPUT: buncha SuperMAG stuff
OUTPUT: no vars, just a plot is made

opt: 0 (default) - plots regular
    1 - plots with vertical lines dividing the days
    2 - plots with vertical lines dividing the days, and the dates written just above the bottom axis
    3 - does the same as above, but uses day number format instead of dates
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from Code.subfun_textNice import textNice
from Code.subfun_figFitter import figFitter

def GRITI_SuperMAG_plot(SuperMAG_data, settings_SuperMAG, dates, settings_plot, \
        opt=0, time_Ref=None, time_Reference=None, settings_paths=None):
    
    #-----Unpack-----
    dateRange_full = dates['date range full']; #get the full date range in date format
    dateRange_dayNum_full = dates['date range full dayNum']; #get the full date range in daynum format
    dateRange_zeroHr = dates['date range zero hr'];
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum'];
    dateRange_zeroHr_monthName = dates['date range zero hr month name'];
    dateRange_zeroHr_dayPostfix = dates['date range zero hr day post fix'];
    SuperMAG_timeUnique = SuperMAG_data['time unique']; #get the time unique
    SuperMAG_plotSet = settings_SuperMAG['plot set']; #get the plot set to plot
    FONT_titleFM = settings_plot['font title FM'];
    FONT_axisLabelFM = settings_plot['font axis label FM'];
    PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
    PLOT_lineWidthSmollerest = settings_plot['line width']['smollerest']; #get the line widths

    #----SuperMAG can be longer than the dateRange due to time shift plotting, don't want that here----
    kr = (SuperMAG_timeUnique < dateRange_dayNum_full[0,1]*86400) | (SuperMAG_timeUnique >= (dateRange_dayNum_full[-1,1]+1)*86400); #get stuff outside the date range
    kr_not = (SuperMAG_timeUnique >= dateRange_dayNum_full[0,1]*86400) & (SuperMAG_timeUnique < (dateRange_dayNum_full[-1,1]+1)*86400); #get stuff outside the date range
    #!! NO YEAR SUPPORT !!
    SuperMAG_timeUnique = np.delete(SuperMAG_timeUnique,kr); #delete em    

    #-----Plot SuperMAG results as several time series-----
    SuperMAG_timeUnique_hr = (SuperMAG_timeUnique - dateRange_dayNum_zeroHr[1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
    
    if( np.mod(np.round(np.min(SuperMAG_timeUnique_hr)),2) == 0 ):
        SuperMAG_time_hr_axis_min = np.round(np.min(SuperMAG_timeUnique_hr)); #is even, good to go
    else:
        SuperMAG_time_hr_axis_min = np.round(np.min(SuperMAG_timeUnique_hr))+1; #is odd, make even
    #END IF
    if( np.mod(np.round(np.max(SuperMAG_timeUnique_hr)),2) == 0 ):
        SuperMAG_time_hr_axis_max = np.round(np.max(SuperMAG_timeUnique_hr)); #is even, good to go
    else:
        SuperMAG_time_hr_axis_max = np.round(np.max(SuperMAG_timeUnique_hr))-1; #is odd, make even
    #END IF
    
    xAxisTicks = np.arange(SuperMAG_time_hr_axis_min,SuperMAG_time_hr_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
        
    if( 'fancy plot' in settings_plot ):
        FLG_fancyPlot = settings_plot['fancy plot']; #read it
        if( (settings_paths == None) & (FLG_fancyPlot >= 1) ):
            print('WARNING in GRITI_SuperMAG_plot: FLG_fancyPlot >= 1 but no settings_paths provided. Disabling FancyPlot(TM) since a path is needed to save the file to a location.');
            FLG_fancyPlot = 0; #set it
        #END IF
    else:
        FLG_fancyPlot = 0; #set it
    #END IF
    
    #Start the SuperMAG plot    
    #--- Prep Plot ---
    if( FLG_fancyPlot == 0 ):
        fig, ax = plt.subplots(nrows=len(SuperMAG_plotSet), ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        journal_dpi = settings_plot['journal dpi'];
        print('MAKING FANCY PLOT: SuperMAG_products IN fancyPlot FOLDER'); #report since you won't see anything
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(nrows=len(SuperMAG_plotSet), ncols=1,figsize=(14,14.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    #END IF
    
    #plot each one systematically
    for i in range(0,len(SuperMAG_plotSet)):
        SuperMAG_dataToPlot = SuperMAG_data[SuperMAG_plotSet[i]][~kr]; #get the data to plot
        ax[i].set_aspect('auto'); #Remove the aspect ratio from the basemap so it fills the screen better
        ax[i].plot( SuperMAG_timeUnique_hr, SuperMAG_dataToPlot, linewidth=PLOT_lineWidthSmollerest ); #plot
        
        if( (time_Reference != None) & np.all(time_Ref != None) ):
            if( (np.abs((np.min(SuperMAG_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.min(time_Ref))/3600 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'OMNI') & (time_Reference != 'SuperMAG')) ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
                ax[i].plot( np.repeat( (np.min(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , 2) , (np.min(SuperMAG_dataToPlot),np.max(SuperMAG_dataToPlot)), linewidth=PLOT_lineWidthRegular, color='r'); #plot red lines showing ISR data time
            if( (np.abs((np.max(SuperMAG_timeUnique_hr)/24 + dateRange_dayNum_zeroHr[1])*86400 - np.max(time_Ref))/3600 >= 0.25) & ((time_Reference != 'Kp') & (time_Reference != 'OMNI') & (time_Reference != 'SuperMAG')) ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
                ax[i].plot( np.repeat( (np.max(time_Ref) - dateRange_dayNum_zeroHr[1]*86400)/3600 , 2) , (np.min(SuperMAG_dataToPlot),np.max(SuperMAG_dataToPlot)), linewidth=PLOT_lineWidthRegular, color='r'); #plot red lines showing ISR data time
            #END IF
        #END IF
        
        ax[i].set_xticks(xAxisTicks); #set x axis ticks
        if( i != (len(SuperMAG_plotSet)-1)):
            ax[i].set_xticklabels([]); #if statement to remove x axis labels except for the last line
            ax[i].tick_params(axis="x",direction="in");
        #END IF
        ax[i].set_xlim( SuperMAG_time_hr_axis_min , SuperMAG_time_hr_axis_max ); #set y axis limits
        
        if( len(SuperMAG_plotSet) < 6 ):
            ax[i].set_ylabel(settings_SuperMAG['labels'][SuperMAG_plotSet[i]]+settings_SuperMAG['units'][SuperMAG_plotSet[i]],fontproperties=FONT_axisLabelFM); #set the y axis label
        else:
            ax[i].set_ylabel(settings_SuperMAG['labels'][SuperMAG_plotSet[i]],fontproperties=FONT_axisLabelFM); #set the y axis label [drop units for dense plot]
        #END IF
        
        ax[i].set_ylim( np.nanmin(SuperMAG_dataToPlot) , np.nanmax(SuperMAG_dataToPlot) ); #set y axis limits
        
        ax[i].grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on
        
        if( opt >= 1 ):
            for j in range(0,dateRange_full.shape[0]):
                #run through each day, print the line
                if( j != 0 ):
                    ax[i].plot( np.repeat( (dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24, 2) , (np.min(SuperMAG_dataToPlot),np.max(SuperMAG_dataToPlot)), linewidth=PLOT_lineWidthSmoller, color='k', linestyle='--'); #plot black dotted lines showing the date separation
                #END IF
                if( (i == 0) & (opt >= 2) ): #only print on the top plot
                    if( opt == 2 ):
                        ax[i].text((dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(SuperMAG_dataToPlot)-(np.max(SuperMAG_dataToPlot)-np.min(SuperMAG_dataToPlot))*.15, str(dateRange_full[j,1])+'/'+str(dateRange_full[j,2])+'/'+str(dateRange_full[j,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
                    elif( opt == 3 ):
                        ax[i].text((dateRange_dayNum_full[j,1] - dateRange_dayNum_zeroHr[1])*24+10, np.max(SuperMAG_dataToPlot)-(np.max(SuperMAG_dataToPlot)-np.min(SuperMAG_dataToPlot))*.15, str(dateRange_dayNum_full[j,1])+', '+str(dateRange_dayNum_full[j,0]), fontproperties=FONT_axisLabelFM); #print the text saying the days
                    #END IF
                #END IF
            #END FOR j
        #END IF
    
    #END FOR i
    
    #do things for the whole plot to finish up now
    if( FLG_fancyPlot == 0 ):
        string_Title = 'SuperMAG-Sourced Data ('+textNice(SuperMAG_data['data rate']/60)+' min resolution) for '+str(dateRange_full[0,1])+'/'+str(dateRange_full[0,2])+ \
            '/'+str(dateRange_full[0,0])+' to '+str(dateRange_full[-1,1])+ \
            '/'+str(dateRange_full[-1,2])+'/'+str(dateRange_full[-1,0])+ \
            ' (M/D/Y)'; #create mecha title
        ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    #END IF
    
    ax[i].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    
    #====== Standard Figure Fit & Saving =====
    figFitter(fig); #fit that fit fast
    if( FLG_fancyPlot == 0 ):
        plt.show(); #req to make plot show up
    else:
        tryLim = 3; #number of times to try
        tryNum = 0; #try counter
        FLG_win = False; #prep flag
        while( (tryNum < tryLim) & (FLG_win == False) ):
            try:
                fig.savefig(settings_paths['fancyPlots']+'\\'+'SuperMAG_products'+settings_plot['save file type']); #save the figure
                plt.close(); #close figure b/c it lurks apparently
                plt.ion(); #re-enable it for later stuff
                FLG_win = True;
            except PermissionError:
                print('WARNING in GRITI_SuperMAG_plot: Try #'+str(tryNum+1)+'/'+str(tryLim)+', error saving a FancyPlot(TM): PermissionError for the file '+settings_paths['fancyPlots']+'\\'+'SuperMAG_products'+settings_plot['save file type']+' Waiting 3 seconds.');
                try:
                    time.sleep(3);
                except:
                    import time
                    time.sleep(3);
                #END TRY
                tryNum += 1; #increment failure
            #END TRY
        #END WHILE
        if( tryNum == tryLim ):
            print('\nWARNING in GRITI_SuperMAG_plot: Error saving a FancyPlot(TM): PermissionError for the file '+settings_paths['fancyPlots']+'\\'+'SuperMAG_products'+settings_plot['save file type']+' Disabling FancyPlot(TM), plot will likely be real messed.');
            plt.ion(); #re-enable it for later stuff
            plt.show(); #req to make plot show up
        #END IF
    #END IF