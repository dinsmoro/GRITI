#GOAL: Plot only  keogram
#RD on 3/02/21
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
import warnings
from subfun_filter import subfun_filter
from subfun_spectra import subfun_spectra
from subfun_figFitter import figFitter
    
def GRITI_spectral_analysisPlot(dataz, timez, timezUnit, dataRate, dataRateUnit, plotTimeUnit, \
        filtMethod, spectraMethod, dates, settings_spectra, settings_plot, settings_paths, datazType , titleOverride=None, reduceWindow=0, \
        FLG_fancyPlot = 0):
    
    if( FLG_fancyPlot == 1 ):
        print('MAKING FANCY PLOT: spectra_'+datazType+' IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    
    #----- Unpack -----
    dateRange = dates['date range'];
    FONT_axisLabelFM = settings_plot['font axis label FM'];
    FONT_titleFM = settings_plot['font title FM'];
    journal_dpi = settings_plot['journal dpi'];
    
    #--- Prepare for plotting multiple things on a plot ---
    if( isinstance(dataz,list) == True ):
        if( isinstance(dataz[0],list) == True ):
            datazOverallLen = len(dataz); #get the list-of-a-list size
            datazLen = len(dataz[0]); #get the data length
        else:
            datazOverallLen = 1; #set the list-of-a-list size
            datazLen = len(dataz); #get the data length
            #need to wrap everything as a list-of-a-list so that the code doesn't fail
            dataz = [dataz]; #wrap it
            timez = [timez]; #wrap it
            timezUnit = [timezUnit]; #wrap it
            dataRate = [dataRate]; #wrap it
            dataRateUnit =  [dataRateUnit]; #wrap it
            filtMethod = [filtMethod]; #wrap it
            spectraMethod = [spectraMethod]; #wrap it
            datazType = [datazType]; #wrap it
        #END IF
    else:
        #need to wrap everything as a list-of-a-list so that the code doesn't fail
        datazOverallLen = 1; #set the list-of-a-list size
        datazLen = 1; #only 1
        dataz = [[dataz]]; #wrap it
        timez = [[timez]]; #wrap it
        timezUnit = [[timezUnit]]; #wrap it
        dataRate = [[dataRate]]; #wrap it
        dataRateUnit =  [[dataRateUnit]]; #wrap it
        filtMethod = [[filtMethod]]; #wrap it
        spectraMethod = [[spectraMethod]]; #wrap it
        datazType = [[datazType]]; #wrap it
    #END IF
    
    if( plotTimeUnit == 'sec' ):
        plotTimeAdj = 1; #converter to plot time unit
    elif( plotTimeUnit == 'min' ):
        plotTimeAdj = 60; #converter to plot time unit
    elif( plotTimeUnit == 'hr' ):
        plotTimeAdj = 3600; #converter to plot time unit
    elif( plotTimeUnit == 'day' ):
        plotTimeAdj = 86400; #converter to plot time unit
    #END IF
    
    #PLOT IT UP
    warnings.filterwarnings("ignore", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
    
    #Start the spectra plot
    if( FLG_fancyPlot == 0 ):
        fig, ax = plt.subplots(nrows=datazOverallLen, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(nrows=datazOverallLen, ncols=1, figsize=(14,8.5), dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
    #END IF
    if( datazOverallLen == 1 ):
        ax = [ax]; #wrap it because lists 
    #END IF
    
    cntr = 0; #cntr for line colors
    for j in range(0,datazOverallLen):
        pZ = []; #prep
        lZ = []; #prep
        for i in range(0,datazLen):
            #----- Time Prep -----
            if( timezUnit[j][i] == 'min' ):
                timez[j][i] = timez[j][i]*60; #convert to sec
            elif( timezUnit[j][i] == 'hr' ):
                timez[j][i] = timez[j][i]*3600; #convert to sec
            elif( timezUnit == 'day' ):
                timez[j][i] = timez[j][i]*86400; #convert to sec
            #END IF
            if( dataRate[j][i] == None ):
                dataRate[j][i] = np.median(np.diff(timez)); #estimate the data rate
            #END IF
            if( dataRateUnit[j][i] == 'min' ):
                dataRate[j][i] = np.int64(dataRate[j][i]*60); #convert to sec
            elif( dataRateUnit[j][i] == 'hr' ):
                dataRate[j][i] = np.int64(dataRate[j][i]*3600); #convert to sec
            elif( dataRateUnit[j][i] == 'day' ):
                dataRate[j][i] = np.int64(dataRate[j][i]*86400); #convert to sec
            #END IF
            
            dataz_filt = subfun_filter(dataz[j][i], timez[j][i], filtMethod[j][i], settings_spectra, dataRate = dataRate[j][i], reduceWindow = reduceWindow); #filter
            
            if( spectraMethod[j][i].lower() == 'fft' ):
                Cxx, period = subfun_spectra( dataz_filt, timez[j][i], spectraMethod[j][i], settings_spectra, dataRate = dataRate[j][i], reduceWindow = reduceWindow); #get spectra
            elif( spectraMethod[j][i].lower().replace('-','').replace(' ','') == 'lombscargle'):
                Cxx, period, _ = subfun_spectra( dataz_filt, timez[j][i], spectraMethod[j][i], settings_spectra, dataRate = dataRate[j][i], reduceWindow = reduceWindow); #get spectra
            #END IF
            
            # Fs = 1/dataRate; #1/sec, time delta in freq form
            
            pT, = ax[j].plot( period/plotTimeAdj, Cxx, color=settings_plot['color'][cntr], linewidth=settings_plot['line width']['plus'], linestyle=settings_plot['line style'][cntr] ); #plot
        #        ax[j].plot( OMNI_data_scargPeriod, np.tile(OMNI_data_scarggf,np.size(OMNI_data_scargPeriod)) , color="xkcd:grey" ); #plot
        
            pZ.append(pT); #add onto the plot list
            lZ.append(datazType[j][i]); #add onto the legend list
            cntr += 1; #increment
        #END FOR i
    
        ax[j].set_ylabel('Normalized Power',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        xAxisTicks = np.int64(np.arange( 0, settings_spectra['period limit max']+10*60, 10*60)/plotTimeAdj); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
        ax[j].set_xticks(xAxisTicks); #set x axis ticks
        ax[j].set_xlim( (settings_spectra['period limit min']/plotTimeAdj, settings_spectra['period limit max']/plotTimeAdj) ); #set x axis limits
        ax[j].legend(pZ, lZ, loc='upper right');
    #        if( OMNI_data_scarggf < np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ):
    #            ax[j].set_ylim( (0, np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim])+0.1*np.max(OMNI_data_scargPower[OMNI_data_scargPeriod<=plot_Period_Lim]) ) ); #set x axis limits
    #        else:
    #            ax[j].set_ylim( (0, OMNI_data_scarggf+0.1*OMNI_data_scarggf ) ); #set x axis limits
    #        #END IF 
        if( FLG_fancyPlot == 0 ):
            if( titleOverride != None ):
                if( titleOverride == '&' ): #if it's just & then join with &
                    string_title = ' & '.join(list(set(datazType)))+' - '+' & '.join(list(set(filtMethod)))+' '+' & '.join(list(set(spectraMethod)))+' Power Spectra'; #create mecha title
                else:
                    string_title = titleOverride; #use this for the title
                #END IF
                ax[j].set_title(string_title,fontproperties=FONT_titleFM); #set the title
            #END IF
        #END IF
    #END FOR j
    ax[j].set_xlabel('Periods ['+plotTimeUnit+']'+' for Date Range '+str(dateRange[0,1])+'/'+str(dateRange[0,2])+ \
        '/'+str(dateRange[0,0])+' to '+str(dateRange[-1,1])+ '/'+str(dateRange[-1,2])+'/'+str(dateRange[-1,0])+ ' (M/D/Y)',fontproperties=FONT_axisLabelFM); #set the x axis label
    
    figFitter(fig); #fit the fig fast
    # fig.subplots_adjust(left = 0.050, right = 0.985, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    if( FLG_fancyPlot == 0 ):
        plt.show(); #req to make plot show up
    else:
        fig.savefig(settings_paths['fancyPlots']+'\\'+'spectra_6min_'+datazType+'.png'); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF  
    
    warnings.filterwarnings("default", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
#END DEF