import numpy as np
import matplotlib.pyplot as plt
from Code.GRITI_AMPERE_integrator import GRITI_AMPERE_integrator
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.subfun_figFitter import figFitter
    
def GRITI_AMPERE_integrator_plot(AMPERE_data, timeRef, dates, settings_AMPERE, \
        settings_map, settings_plot, settings_paths, FLG_fancyPlot = 0, highlighter = False):
    
    if( FLG_fancyPlot >= 1 ):
        print('MAKING FANCY PLOT: '+settings_AMPERE['data type']+'_integrate IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    
    #Unpack everything else
    plotLatRange = settings_map['lat range'];
    plotLongRange = settings_map['long range'];
    FONT_titleFM = settings_plot['font title FM'];
    FONT_axisLabelFM = settings_plot['font axis label FM'];
    journal_dpi = settings_plot['journal dpi'];
    
    #Unpack line widths
    PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    
    #----- Integrate AMPERE Data -----    
    # AMPERE_integrated = GRITI_AMPERE_integrator(AMPERE_data, dates, settings_AMPERE, plotLatRange, plotLongRange, settings_AMPERE['integrate method'], settings_AMPERE['integrate method lat val'], AMPERE_integrateMethod_log=settings_AMPERE['integrate method log']); #integrate with the integrator function
    AMPERE_integrated = AMPERE_data['integrated']

    # #--- Time match to 6 minutes if needed ---
    # if( np.isclose(AMPERE_data['data rate'],360.) == False ):
    #     sixMin_timeUnique = np.arange(dates['date range zero hr hour bounds'][0]*3600,dates['date range zero hr hour bounds'][1]*3600,360); #sec, arange time stamps in 6 minute steps
    #     AMPERE_integrated, AMPERE_timeUnique_matched = subfun_timeMatch(AMPERE_integrated, (AMPERE_data['time unique']-dates['date range zero hr dayNum'][1]*86400), sixMin_timeUnique, timeMatch_delta=AMPERE_data['data rate'], FLG_useSum=1); #time match alg to align to 6 minute cadence, add because it's a count (?)
    #     AMPERE_timeUnique_hr = (AMPERE_timeUnique_matched)/3600; #hr, convert to hr with 0 hr at specified day
    # else:
    AMPERE_timeUnique_hr = (AMPERE_data['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
    # #END IF
    
    #--- Prepare axis limits ---
    # if( np.mod(np.round(np.min(AMPERE_timeUnique_hr)),2) == 0 ):
    #     AMPERE_time_hr_axis_min = np.round(np.min(AMPERE_timeUnique_hr)); #is even, good to go
    # else:
    #     AMPERE_time_hr_axis_min = np.round(np.min(AMPERE_timeUnique_hr))+1; #is odd, make even
    # #END IF
    # if( np.mod(np.round(np.max(AMPERE_timeUnique_hr)),2) == 0 ):
    #     AMPERE_time_hr_axis_max = np.round(np.max(AMPERE_timeUnique_hr)); #is even, good to go
    # else:
    #     AMPERE_time_hr_axis_max = np.round(np.max(AMPERE_timeUnique_hr))-1; #is odd, make even
    # #END IF
    
    #-----BEGIN THE PLOTTING!------
    AMPERE_plot_label = settings_AMPERE['labels'][settings_AMPERE['data type']]+settings_AMPERE['units'][settings_AMPERE['data type']]; #get the label
    AMPERE_plot_label_noUnits = settings_AMPERE['labels'][settings_AMPERE['data type']]; #remove the (units)
    
    #Start the AMPERE plot
    if( FLG_fancyPlot == 0 ):
        fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    #END IF
    
    ax.set_aspect('auto'); #Remove the aspect ratio from the basemap so it fills the screen better
    ax.plot( AMPERE_timeUnique_hr, AMPERE_integrated , linewidth=PLOT_lineWidthRegular ); #plot
    
    if( np.any(highlighter != False) ):        
        pltr_time = np.copy(AMPERE_data['time unique']); #copy time
        k = np.logical_not(highlighter); #remove times not highlighted
        pltr_time = np.delete(pltr_time,k,axis=0); #delete em

        k = np.isclose(np.diff(pltr_time), np.median(np.diff(AMPERE_data['time unique']))); #find contiguous time bits
        pltr_time = (pltr_time - dates['date range zero hr dayNum'][1]*86400)/3600; #for plottin, adjust to hrs
        kj = np.where(~k)[0]+1; #add 1 for index b/c diff
        if(kj[0] != 0):
            kj = np.insert(kj, 0, 0); #insert 0
        #END IF
        if(kj[-1] != (pltr_time.size-1) ):
            kj = np.append(kj, pltr_time.size-1); #insert end size
        #END IF
        for j in range(0,kj.size-1):
            ax.axvspan(pltr_time[kj[j]], pltr_time[kj[j+1]-1], ymin=0, ymax=1, alpha=0.45, edgecolor='none', facecolor='xkcd:brick red');
        #END FOR j
    #END IF
    
    if( (np.abs((np.min(AMPERE_timeUnique_hr)*3600 + dates['date range zero hr dayNum'][1]*86400) - np.min(timeRef))/3600 >= 0.25) & ((settings_plot['time ref'] != 'Kp') & (settings_plot['time ref'] != 'OMNI') & (settings_plot['time ref'] != 'SuperMAG') & (settings_plot['time ref'] != 'AMPERE')) ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
        ax.plot( np.repeat( (np.min(timeRef) - dates['date range zero hr dayNum'][1]*86400)/3600 , 10) , np.linspace(np.min(AMPERE_integrated),np.max(AMPERE_integrated),num=10), linewidth=1.75, color='r'); #plot red lines showing ISR data time
    if( (np.abs((np.max(AMPERE_timeUnique_hr)*3600 + dates['date range zero hr dayNum'][1]*86400) - np.max(timeRef))/3600 >= 0.25) & ((settings_plot['time ref'] != 'Kp') & (settings_plot['time ref'] != 'OMNI') & (settings_plot['time ref'] != 'SuperMAG') & (settings_plot['time ref'] != 'AMPERE')) ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
        ax.plot( np.repeat( (np.max(timeRef) - dates['date range zero hr dayNum'][1]*86400)/3600 , 10) , np.linspace(np.min(AMPERE_integrated),np.max(AMPERE_integrated),num=10), linewidth=1.75, color='r'); #plot red lines showing ISR data time
    #END IF
    
    # xAxisTicks = np.arange(AMPERE_time_hr_axis_min,AMPERE_time_hr_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    # ax.set_xticks(xAxisTicks); #set x axis ticks
    
    # ax.set_xlim( AMPERE_time_hr_axis_min , AMPERE_time_hr_axis_max ); #set y axis limits
    
    GRITI_plotHelper_axisizerTime(AMPERE_timeUnique_hr,ax=ax); #do all the x axis stuff
    
    ax.set_ylabel('Integrated '+AMPERE_plot_label,fontproperties=FONT_axisLabelFM); #set the y axis label
    
    ax.set_ylim( np.min(AMPERE_integrated) , np.max(AMPERE_integrated) ); #set y axis limits
    
    ax.grid(visible=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on        
    
    if( FLG_fancyPlot == 0 ):
        string_title = 'Integrated AMPERE '+AMPERE_plot_label_noUnits; #prep the title
        if( settings_AMPERE['integrate method'] == 0 ):
            string_title = string_title + ' within keo area'; #add to mecha title
        elif( settings_AMPERE['integrate method'] == 1 ):
            string_title = string_title + ' within keo long & up to pole'; #add to mecha title
        elif( settings_AMPERE['integrate method'] == 2 ):
            string_title = string_title + ' within keo long & up to '+str(settings_AMPERE['integrate method lat val'])+' degc lat'; #add to mecha title
        elif( settings_AMPERE['integrate method'] == 3 ):
            if( (np.min(plotLatRange) < 0) & (np.max(plotLatRange) >= 0) ):
                string_title = string_title + ' both Hemispheres'; #add to mecha title
            else:
                if( (np.min(plotLatRange) >= 0) & (np.max(plotLatRange) >= 0) ):
                    #northern hemisphere
                    string_title = string_title + ' Northern Hemisphere'; #add to mecha title
                else:
                    #southern hemisphere
                    string_title = string_title + ' Southern Hemisphere'; #add to mecha title
                #END IF
            #END IF
        #END IF
        string_title += ' '+str(dates['date range'][0,1])+'/'+str(dates['date range'][0,2])+ \
        '/'+str(dates['date range'][0,0])+' to '+str(dates['date range'][-1,1])+ \
        '/'+str(dates['date range'][-1,2])+'/'+str(dates['date range'][-1,0])+ \
        ' (M/D/Y)'; #create mecha title
        ax.set_title(string_title,fontproperties=FONT_titleFM); #set the title
    #END IF
    
    ax.set_xlabel('Time in UT [hr] - 0 Hr on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr day post fix']+' | Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    
    figFitter(fig); #fit the fig fast
    if( FLG_fancyPlot != 0 ):
        fig.savefig(settings_paths['fancyPlots']+'\\'+settings_AMPERE['data type']+'_integrate'+settings_plot['save file type']); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF