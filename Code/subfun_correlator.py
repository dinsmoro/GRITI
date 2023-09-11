"""
Finds the optimal correlation between two signals sig1 & sig2
sig1 & sig2 must be on the same cadence.
- modes:
    corr : calcs correlation coeff for entirety of sig1, sig2 - time is not needed
    shift : calcs optimal correlation coeff by shifting sig2 around in time. sig1 & sig2 must be aligned in time and the same length.
    range : calcs optimal correlation coeff for a specific time range for sig1 by shifting sig2 around in time. sig1 & sig2 must be at the same cadence and sig2 needs to have enough data to be able to shift around the time range.
        timeRange - the time range in np.array((t1,t2)) or (t1,t2) or [t1,t2] format where t1 is less than t2
    interval : 
        timeInterval - time length required to be used as an interval for correlation
        intervalType - 'pos', 'neg', 'both' -> checks for intervals with positive/negative/pos+neg stretches that are at least timeInterval long
        intervalStepScaler - % of data to set as the min length to fit a linear fit to, default is 1% of total data size
        intervalFitScaler - R-Squared ejector for fitting linear fit to data, if new rSq is X times more than previous max it quits (not sure how used this is, never profiled to see what that alg does other than it works ( ͡° ͜ʖ ͡°))
        intervalFitEndEncourager - Fit coefficient that encourages the ends of the linear fits to connect, less makes them connect less and more makes them connect more to the detriment of fitting the data (though how detrimental idk)
        intervalFitMinRSq - 0 to 1, Minimum rSq value to use for interval time checks, idea is if its below X rSq it can't be trusted to be used as a positive or negative incidence
    interval manual :
        timeInterval - the time ranges to be used as intervals
- dataRate can be supplied or it will be estimated from time
- timeLimit is how long to check for correlation (e.g., 3 hours in the time-units provided)
- FLG_interpGaps will interpolate time gaps in sig1 and sig2 automatically if set to True
- reportDivisor will divide the optimal time offset value by the value FLG_reportDivisor to make it the units you want it to be, use [60,' min'] to convert sec to minutes and report the units of minutes
- reportRounder rounds the calc'd corr coeffs to the decimal precision given (default 3 decimal places)
- FLG_shiftDir ['pos','neg','both'] tells to investigate only positive, negative, or both positive and negative time shifts
- FLG_plot plots the correlation coefficients vs time in addition to reporting the corr coeff/time that maximizes the correlation coeff [does not do anything for mode='corr' b/c it's just a number]
- settings_plot needs to be provided if FLG_plot is on and it holds general plot settings that are unpacked at the beginning of the plot alg
- plotName is the name for the plot data (like "Joule Heating & PC(N)")
- FLG_fancyPlot plots the correletion coeff plot as a fancyPlot(TM)
"""
#this was written in a feverdream (incl. subfun_peakFinder) driven by michelle branch's enitre discography on repeat [it's like 4 albums]
#(best song is a shootout between fault line and I'd rather be in love, I'd know)
import numpy as np
from Code.subfun_filter import subfun_filter
from scipy.ndimage.filters import uniform_filter1d
from Code.subfun_textNice import textNice
from Code.subfun_strstr import strstr
import matplotlib.pyplot as plt
from Code.subfun_figFitter import figFitter
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
# from mpl_toolkits.axes_grid1 import make_axes_locatable

def subfun_correlator(sig1, sig2, mode='corr', sig1_filt=None, sig2_filt=None, sig1_noise=None, sig2_noise=None, time1=None, time2=None, dataRate=None, timeLimit=None, timeRange=None, timeInterval=None, \
                      intervalType='pos', intervalStepScaler=1, intervalFitScaler=10, intervalFitEndEncourager=1000, intervalFitMinRSq=0.5, \
                      FLG_interpGaps=False, reportDivisor=[1,'time-units'], reportRounder=3, FLG_shiftDir='both', FLG_plot=False, settings_plot=None, \
                      settings_paths=None, settings_spectra=None, plotName='', FLG_enableText=True, FLG_fancyPlot=0):
    #--- Make sure needed info is provided ---
    if( (mode != 'corr') & np.all(time2 == None)):
        error_text = 'Error in subfun_correlator: Mode is "'+mode+'" but no time2 input was provided. Returning an error.';
        if( FLG_enableText ):
            print(error_text); #report error
        #END IF
        return error_text #very helpful
    #END IF
    if( ((mode == 'range') | (mode == 'interval') | (mode == 'interval manual')) & np.all(time1 == None) ):
        error_text = 'Error in subfun_correlator: Mode is "'+mode+'" but no time1 input was provided. Returning an error.';
        if( FLG_enableText ):
            print(error_text); #report error
        #END IF
        return error_text #very helpful
    #END IF
    if(  (mode == 'corr') & (sig1.size != sig2.size) ):
        error_text = 'Error in subfun_correlator: sig1 and sig2 are not the same sizes ('+textNice(sig1.size)+' and '+textNice(sig2.size)+', respectively) and mode "'+mode+'" requires that. Returning an error.';
        if( FLG_enableText ):
            print(error_text); #report error
        #END IF
        return error_text #very helpful
    #END IF
    if( (mode == 'range') & np.all(timeRange == None) ):
        error_text = 'Error in subfun_correlator: Mode is "'+mode+'" but no timeRange was provided. Returning an error.';
        if( FLG_enableText ):
            print(error_text); #report error
        #END IF
        return error_text #very helpful
    #END IF
    if( (mode == 'interval') & np.all(timeInterval == None) ):
        error_text = 'Error in subfun_correlator: Mode is "'+mode+'" but no timeInterval was provided. Returning an error.';
        if( FLG_enableText ):
            print(error_text); #report error
        #END IF
        return error_text #very helpful
    #END IF
    if( (FLG_interpGaps == True) & (np.all(time1 == None) | np.all(time2 == None)) ):
        if( np.all(time1 == None) & np.all(time2 == None) ):
            error_text = 'Error in subfun_correlator: FLG_interpGaps is [True] but no time1 or time2 were provided. Returning an error.';
            if( FLG_enableText ):
                print(error_text); #report error
            #END IF
        elif( np.all(time1 == None) & np.all(time2 != None) ):
            error_text = 'Error in subfun_correlator: FLG_interpGaps is [True] but no time1 was provided. Returning an error.';
            if( FLG_enableText ):
                print(error_text); #report error
            #END IF
        else:
            error_text = 'Error in subfun_correlator: FLG_interpGaps is [True] but no time2 was provided. Returning an error.';
            if( FLG_enableText ):
                print(error_text); #report error
            #END IF
        #END IF
        return error_text #very helpful
    #END IF
    
    #--- Say it's go time ---
    if( FLG_enableText ):
        if( plotName == '' ):
            print('-- Correlation Calcing ---');
        else:
            print('-- On: '+plotName.replace('$\mathregular{','').replace('}$','')+' ---');
        #END IF
    #END IF

    #--- Make sure folder is made if needed ---
    if( FLG_fancyPlot >= 1 ):
        import os
        path_corrFancyPlot = os.path.join(settings_paths['fancyPlots'],'Correlations');
        if( not os.path.isdir(path_corrFancyPlot) ):
            os.makedirs(path_corrFancyPlot); #make it if not here
        #END IF
        journal_dpi = settings_plot['journal dpi'];
        if( plotName != '' ):
            plotName_fileName = '_'+plotName.replace(' ','-').replace('$\mathregular{_','').replace('}$','').replace('/','-');
        else:
            plotName_fileName = plotName.replace(' ','-').replace('$\mathregular{_','').replace('}$','').replace('/','-');
        #END IF
        
    #END IF

    #--- Prepare return variable ---
    corrRet = {}; #prep
    
    #--- Get time info together if needed ---
    if( (mode != 'corr') | (FLG_interpGaps == True) ): #non-corr modes need time info to do their thing
        #--- Estimate dataRate if not avail ---
        if( dataRate == None ):
            dataRate = np.median(np.diff(time2)); #$TIME, delta of time between readings
            if( np.isclose(np.mod(dataRate,1),0.0) ):
                dataRate = np.int64(dataRate); #if it's an integer, make it one to improve accuracy/alignment w/o floating point errors
            #END IF
        #END IF
        
        #--- Estimate time limit if not avail ---
        if( timeLimit == None ):
            timeLimit = dataRate*30; #$TIME, time limit set to 30 times the data rate to explore enough space but not too much (gives 3 hrs of time to check out with a 6 min cadence, seems reasonable)
        #END IF
        
        #--- Create time shift array ---
        if( FLG_shiftDir == 'pos' ):
            timeShift = np.arange(0, timeLimit+dataRate, dataRate); #$TIME, create a variable that covers the time shift
        elif( FLG_shiftDir == 'neg' ):
            timeShift = -np.flip(np.arange(0, timeLimit+dataRate, dataRate)); #make a negative array by flipping and negating the positive array
        elif( FLG_shiftDir == 'both' ):
            timeShift = np.concatenate((np.flip(np.arange(-dataRate, -timeLimit-dataRate, -dataRate)),np.arange(0, timeLimit+dataRate, dataRate))); #negative option is built like this to always ensure there's a 0
        else:
            error_text = 'Error in subfun_correlator: FLG_shiftDir is ['+str(FLG_shiftDir)+'] but whatever that is, it\'s not "pos", "neg", or "both". Returning an error.';
            if( FLG_enableText ):
                print(error_text); #report error
            #END IF
            return error_text #very helpful
        #END IF
    #END IF
    
    #--- Interpolate Gaps if Needed ---
    if( FLG_interpGaps == True ):
        if( np.any(np.diff(time1) > dataRate) ):
            if( FLG_enableText ):
                print('Warning in subfun_correlator: There are gaps in time1 that will be interpolated over.');
            #END IF
            sig1 = subfun_filter( sig1, 'interp', dataTime = time1, dataRate = dataRate);
            time1 = np.arange(time1[0],time1[-1]+dataRate,dataRate); #make full time range
        #END IF
        if( np.any(np.diff(time2) > dataRate) ):
            if( FLG_enableText ):
                print('Warning in subfun_correlator: There are gaps in time2 that will be interpolated over.');
            #END IF
            sig2 = subfun_filter( sig2, 'interp', dataTime = time2, dataRate = dataRate);
            time2 = np.arange(time2[0],time2[-1]+dataRate,dataRate); #make full time range
        #END IF
    #END IF
    
    #--- Filter if Needed ---
    if( sig1_filt != None ):
        sig1 = subfun_filter( np.copy(sig1), sig1_filt, dataTime = time1, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 1, FLG_reportNaNs = False); #filter (or not)
    #END IF
    if( sig2_filt != None ):
        sig2 = subfun_filter( np.copy(sig2), sig2_filt, dataTime = time2, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 1, FLG_reportNaNs = False); #filter (or not)
    #END IF
    
    #--- Calculate correlation coefficients for the various modes ---
    if( mode == 'corr' ):
        #CALC PWRS
        pwr_sig1 = np.sqrt(1/sig1.size*np.sum(sig1**2)); #estimate power of signal
        pwr_sig2 = np.sqrt(1/sig2.size*np.sum(sig2**2)); #estimate power of signal
        corr = np.corrcoef(1/pwr_sig1*sig1,1/pwr_sig2*sig2)[0,1];
        corr_textResults = 'Corr coeff: '+str(textNice(np.round(corr,reportRounder)));
        if( FLG_enableText ):
            print(corr_textResults);#print resupts
        #END IF
        corrRet['corr'] = corr; #record
        corrRet['data rate'] = dataRate; #record
        corrRet['plot name'] = plotName; #record
        corrRet['text results'] = corr_textResults; #record
    
    elif( mode == 'shift' ):
        #--- ROLLING TIME shift TO FIND OPTIMAL CORR COEFF ---
        corr = np.zeros(timeShift.size); #preallocate
        for i in range(0,timeShift.size):
            #--- Adjust by the time offset ---
            time_adj = time2 + timeShift[i]; #adjust timeRange by a timeShift for sig2
            
            #sig2 time adjustments
            time_adj_indexes = np.array( ( np.where(np.min(np.abs( time_adj - time2[0] )) == np.abs( time_adj - time2[0] ) )[0][0] , \
                np.where(np.min(np.abs( time_adj - time2[-1] )) == np.abs( time_adj - time2[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
            sig2_adj = sig2[time_adj_indexes[0]:time_adj_indexes[1]+1];
            # if( sig2_filt != None):
            #     if( strstr(sig2_filt.lower().replace('-','').replace(' ',''),'0mean').size > 0 ):
            #         sig2_adj = subfun_filter( np.copy(sig2_adj), '0 mean'); #filter (or not)
            #     #END IF
            # #END IF
            
            #sig1 time adjustments to match size of sig2
            if( timeShift[i] != 0 ):
                sig1_adj = sig1[:sig2_adj.size]; #cut off the end of sig1 to match the size of sig2 (with pos/neg time shift the sig1[0] wants to match sig2_adj[0] and to keep size same so cut off end, future me - I am 95% confident in this)
                # if( sig1_filt != None):
                #     if( strstr(sig1_filt.lower().replace('-','').replace(' ',''),'0mean').size > 0 ):
                #         sig1_adj = subfun_filter( np.copy(sig1_adj), '0 mean'); #filter (or not)
                #     #END IF
                # #END IF
            else:
                sig1_adj = sig1; #timeShift is 0 and its regular
            #END IF
            
            #CALC PWRS
            pwr_sig1_adj = np.sqrt(1/sig1_adj.size*np.sum(sig1_adj**2)); #estimate power of signal
            pwr_sig2_adj = np.sqrt(1/sig2_adj.size*np.sum(sig2_adj**2)); #estimate power of signal
            
            #Real quick side move to calc correlation coefficients
            corr[i] = np.corrcoef(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj*sig2_adj)[0,1];
        #END FOR i
        k = np.where( np.abs(corr) == np.nanmax(np.abs(corr)))[0]; #get maximum values
        corr_textResults = 'Corr coeff max: '+textNice(np.round(corr[k],reportRounder))+' at time shift: '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1];
        if( FLG_enableText ):
            print(corr_textResults);#print resupts
        #END IF
        if( k.size > 0 ):
            corrRet['corr'] = corr[k]; #record
            corrRet['time shift'] = timeShift[k]; #record
        else:
            corrRet['corr'] = np.array((np.nan,));
            corrRet['time shift'] = np.array((np.nan,));
        #END IF
        corrRet['corr vect'] = corr; #record
        corrRet['time shift vect'] = timeShift; #record
        corrRet['data rate'] = dataRate; #record
        corrRet['plot name'] = plotName; #record
        corrRet['text results'] = corr_textResults; #record
        
        if( (FLG_plot == True) & (k.size > 0) ):
            #--- Declare & Unpack ---
            if( (FLG_fancyPlot >= 1) & FLG_enableText ):
                print('MAKING FANCY PLOT: corr_'+mode+plotName_fileName+'_highlite'+' IN fancyPlot FOLDER'); #report since you won't see anything
            #END IF
            # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
            PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
            PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
            # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
            # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
            # PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
            # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
            FONT_titleFM = settings_plot['font title FM'];
            # FONT_axisTick = settings_plot['font axis tick'];
            FONT_axisLabelFM = settings_plot['font axis label FM'];
            
            if( FLG_fancyPlot == 0 ): #make a loop that can do fancy and non-fancy in one (correlation calcs are expensive when done en masse)
                FLG_fancyPlot_runner = np.array( (0,) );
            elif( FLG_fancyPlot == 1 ):
                FLG_fancyPlot_runner = np.array( (0,1) );
            elif( FLG_fancyPlot == 2 ):
                FLG_fancyPlot_runner = np.array( (1,) );
            #END IF
            
            #--- Prep Plot ---
            for figr in range(0, FLG_fancyPlot_runner.size):
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
                    figManager = fig.canvas.manager; #req to maximize
                    figManager.window.showMaximized(); #force maximized
                else:
                    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                    fig, ax = plt.subplots(nrows=2, ncols=1,figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
                #END IF
                fig.subplots_adjust(hspace=0);
                # divider = make_axes_locatable(ax); #prep to add an axis
                # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
                
                #Remove the aspect ratio from the basemap so it fills the screen better
                ax[0].set_aspect('auto');
                ax[1].set_aspect('auto');
                
                #--- Actual Plotting ---
                time_adj = time2 + timeShift[k]; #adjust timeRange by a timeShift for sig2
                time_adj_indexes = np.array( ( np.where(np.min(np.abs( time_adj - time2[0] )) == np.abs( time_adj - time2[0] ) )[0][0] , \
                    np.where(np.min(np.abs( time_adj - time2[-1] )) == np.abs( time_adj - time2[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
                ax[0].plot(time1/reportDivisor[0], sig1, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                ax[0].plot(time1[:sig2[time_adj_indexes[0]:time_adj_indexes[1]+1].size]/reportDivisor[0], sig1[:sig2[time_adj_indexes[0]:time_adj_indexes[1]+1].size], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                
                ax[1].plot(time2/reportDivisor[0], sig2, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                ax[1].plot((time2[time_adj_indexes[0]:time_adj_indexes[1]+1])/reportDivisor[0], sig2[time_adj_indexes[0]:time_adj_indexes[1]+1], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                #END FOR i
                
                #--- Axis and Titles and Stuff ---
                if( FLG_fancyPlot_runner[figr] == 0 ): #only title non-fancy plot
                    string_Title = mode.capitalize()+' Method, 2nd Plot has Time Shift of '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1];
                    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
                #END IF
                ax[1].set_xlabel('Time ['+reportDivisor[1]+']',fontproperties=FONT_axisLabelFM); #set the x axis label
                ax[0].set_ylabel(plotName[:plotName.find(' & ')],fontproperties=FONT_axisLabelFM); #set the y axis label
                ax[1].set_ylabel(plotName[plotName.find(' & ')+3:],fontproperties=FONT_axisLabelFM); #set the y axis label
                #END IF
                
                #nice axis ticks
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],FLG_removeLabels=True,FLG_tickDirIn=True);
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True);
                else:
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=True,FLG_tickDirIn=True);
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=False,FLG_tickDirIn=True);
                #END IF
                
                figFitter(fig); #fit the fig fast
                if( FLG_fancyPlot_runner[figr] != 0 ):
                    fig.savefig(os.path.join(path_corrFancyPlot,'corr_'+mode+plotName_fileName+'_highlite'+'.png')); #save the figure
                    plt.close(); #close figure b/c it lurks apparently
                    plt.ion(); #re-enable it for later stuff
                #END IF
            #END FOR figr
        elif( (FLG_plot == True) & (k.size == 0) ):
            print('WARNING in subfun_correlator: No plots produced because no max found (likely due to one of the data inputs being all NaNs or times not aligning at all.');
        #END IF
    
    elif( mode == 'range' ):
        #--- ROLLING TIME shift TO FIND OPTIMAL CORR COEFF ---
        corr = np.ones(timeShift.size)*np.nan; #preallocate
        time_adj_indexes = np.array( ( np.where(np.min(np.abs( time1 - timeRange[0] )) == np.abs( time1 - timeRange[0] ) )[0][0] , \
                np.where(np.min(np.abs( time1 - timeRange[-1] )) == np.abs( time1 - timeRange[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
        sig1_adj = sig1[time_adj_indexes[0]:time_adj_indexes[1]+1];
        #enforce a local 0 mean
        # if( sig1_filt != None):
        #     if( strstr(sig1_filt.lower().replace('-','').replace(' ',''),'0mean').size > 0 ):
        #         sig1_adj = subfun_filter( np.copy(sig1_adj), '0 mean'); #filter (or not)
        #     #END IF
        # #END IF
        sig1_adj = subfun_filter( np.copy(sig1_adj), '0 mean'); #filter (or not)
        pwr_sig1_adj = np.sqrt(1/sig1_adj.size*np.sum(sig1_adj**2)); #estimate power of signal
        sig1_adj_normd = 1/pwr_sig1_adj*sig1_adj; #calc once b/c static
                
        # #--- fix for timeRange edges extending past
        # if( np.any(timeRange < time1[0]) ):
        #     if( FLG_enableText ):
        #         print('WARNING in subfun_correlator: time1[0] is '+str(time1[0])+' which is more than timeRange[0] '+str(timeRange[timeRange < time1[0]])+' so timeRange[0] will be adjusted to match time1[0]\'s limit so that the analysis can work properly (vector size mismatch would NaN all analysis otherwise).');#print resupts
        #     #END IF
        #     timeRange = timeRange.copy(); #copy protect
        #     timeRange[timeRange < time1[0]] = time1[0]; #adjust the edge of the time range if it extends before time1 (it'll cause time2 to fail otherwise b/c of a vector size mismatch)
        # #END IF
        # if( np.any(timeRange > time1[-1]) ):
        #     if( FLG_enableText ):
        #         print('WARNING in subfun_correlator: time1[-1] is '+str(time1[-1])+' which is less than timeRange[-1] '+str(timeRange[timeRange > time1[-1]])+' so timeRange[-1] will be adjusted to match time1[-1]\'s limit so that the analysis can work properly (vector size mismatch would NaN all analysis otherwise).');#print resupts
        #     #END IF
        #     timeRange = timeRange.copy(); #copy protect
        #     timeRange[timeRange > time1[-1]] = time1[-1]; #adjust the edge of the time range if it extends before time1 (it'll cause time2 to fail otherwise b/c of a vector size mismatch)
        # #END IF
        if( np.any(timeRange < time1[0]) | np.any(timeRange > time1[-1]) ):
            #instead of fix, better to be consistent and just nan if it's not gonna work right
            k = np.empty( shape=(0, ) ); #make a nothing array for later
            if( ('min' in reportDivisor[1]) & (reportDivisor[0] == 60) ):
                corr_textResults = 'Corr coeff is NaN for time range: '+textNice(np.min(timeRange)/3600)+' to '+textNice(np.max(timeRange)/3600)+' hr because one of those was below time1[0] ('+str(time1[0])+') or above time1[-1] ('+str(time1[-1])+').';
            else:
                corr_textResults = 'Corr coeff is NaN for time range: '+textNice(np.min(timeRange)/reportDivisor[0])+' to '+textNice(np.max(timeRange)/reportDivisor[0])+' '+reportDivisor[1]+' because one of those was below time1[0] ('+str(time1[0])+') or above time1[-1] ('+str(time1[-1])+').';
            #END IF
        else:
            for i in range(0,timeShift.size):
                #--- Adjust by the time offset ---
                time_adj = time2 + timeShift[i]; #adjust timeRange by a timeShift for sig2
                
                #sig2 time adjustments
                time_adj_indexes = np.array( ( np.where(np.min(np.abs( time_adj - timeRange[0] )) == np.abs( time_adj - timeRange[0] ) )[0][0] , \
                    np.where(np.min(np.abs( time_adj - timeRange[-1] )) == np.abs( time_adj - timeRange[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
                sig2_adj = sig2[time_adj_indexes[0]:time_adj_indexes[1]+1];
                # if( sig2_filt != None):
                #     if( strstr(sig2_filt.lower().replace('-','').replace(' ',''),'0mean').size > 0 ):
                #         sig2_adj = subfun_filter( np.copy(sig2_adj), '0 mean'); #filter (or not)
                #     #END IF
                # #END IF
                sig2_adj = subfun_filter( np.copy(sig2_adj), '0 mean'); #filter (or not)
                
                if( sig2_adj.size == sig1_adj.size ): 
                    #CALC PWRS
                    pwr_sig2_adj = np.sqrt(1/sig2_adj.size*np.sum(sig2_adj**2)); #estimate power of signal
                    
                    #Real quick side move to calc correlation coefficients
                    corr[i] = np.corrcoef(sig1_adj_normd,1/pwr_sig2_adj*sig2_adj)[0,1];
                #END IF
            #END FOR i
            k = np.where( np.abs(corr) == np.nanmax(np.abs(corr)))[0]; #get maximum values
            if( ('min' in reportDivisor[1]) & (reportDivisor[0] == 60) ):
                #time range more useful in hrs if minutes requested for everything else
                corr_textResults = 'Corr coeff max: '+textNice(np.round(corr[k],reportRounder))+' at time shift: '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1]+'\nfor time range: '+textNice(np.min(timeRange)/3600)+' to '+textNice(np.max(timeRange)/3600)+' hr';
            else:
                corr_textResults = 'Corr coeff max: '+textNice(np.round(corr[k],reportRounder))+' at time shift: '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1]+'\nfor time range: '+textNice(np.min(timeRange)/3600)+' to '+textNice(np.max(timeRange)/3600)+' hr';
            #END IF
        #END IF
        if( FLG_enableText ):
            print(corr_textResults);#print resupts
        #END IF
        if( k.size > 0 ):
            corrRet['corr'] = corr[k]; #record
            corrRet['time shift'] = timeShift[k]; #record
        else:
            corrRet['corr'] = np.array((np.nan,));
            corrRet['time shift'] = np.array((np.nan,));
        #END IF
        corrRet['time range'] = timeRange; #record
        corrRet['data rate'] = dataRate; #record
        corrRet['corr vect'] = corr; #record
        corrRet['time shift vect'] = timeShift; #record
        corrRet['plot name'] = plotName; #record
        corrRet['text results'] = corr_textResults; #record
        
        if( (FLG_plot == True) & (k.size > 0) ):
            #--- Declare & Unpack ---
            if( FLG_fancyPlot >= 1 ):
                if( ('min' in reportDivisor[1]) & (reportDivisor[0] == 60) ):
                    #time range more useful in hrs if minutes requested for everything else
                    plotName_fileName = plotName_fileName+'_'+textNice(np.min(timeRange)/3600)+'to'+textNice(np.max(timeRange)/3600)+'hr';
                else:
                    plotName_fileName = plotName_fileName+'_'+textNice(np.min(timeRange)/reportDivisor[0])+'to'+textNice(np.max(timeRange)/reportDivisor[0])+reportDivisor[1];
                #END IF
                if( FLG_enableText ):
                    print('MAKING FANCY PLOT: corr_'+mode+plotName_fileName+'_highlite'+' IN fancyPlot FOLDER'); #report since you won't see anything
                #END IF
            #END IF
            # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
            PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
            PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
            # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
            # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
            # PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
            # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
            FONT_titleFM = settings_plot['font title FM'];
            # FONT_axisTick = settings_plot['font axis tick'];
            FONT_axisLabelFM = settings_plot['font axis label FM'];
            
            if( FLG_fancyPlot == 0 ): #make a loop that can do fancy and non-fancy in one (correlation calcs are expensive when done en masse)
                FLG_fancyPlot_runner = np.array( (0,) );
            elif( FLG_fancyPlot == 1 ):
                FLG_fancyPlot_runner = np.array( (0,1) );
            elif( FLG_fancyPlot == 2 ):
                FLG_fancyPlot_runner = np.array( (1,) );
            #END IF
            
            #--- Prep Plot ---
            for figr in range(0, FLG_fancyPlot_runner.size):
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
                    figManager = fig.canvas.manager; #req to maximize
                    figManager.window.showMaximized(); #force maximized
                else:
                    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                    fig, ax = plt.subplots(nrows=2, ncols=1,figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
                #END IF
                fig.subplots_adjust(hspace=0);
                # divider = make_axes_locatable(ax); #prep to add an axis
                # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
                
                #Remove the aspect ratio from the basemap so it fills the screen better
                ax[0].set_aspect('auto');
                ax[1].set_aspect('auto');
                
                #--- Actual Plotting ---
                time_adj_indexes = np.array( ( np.where(np.min(np.abs( time1 - timeRange[0] )) == np.abs( time1 - timeRange[0] ) )[0][0] , \
                    np.where(np.min(np.abs( time1 - timeRange[-1] )) == np.abs( time1 - timeRange[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
                sig1_adj = sig1[time_adj_indexes[0]:time_adj_indexes[1]+1];
                ax[0].plot(time1/reportDivisor[0], sig1, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                ax[0].plot(time1[time_adj_indexes[0]:time_adj_indexes[1]+1]/reportDivisor[0], sig1_adj, linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                
                time_adj = time2 + timeShift[k]; #adjust timeRange by a timeShift for sig2
                time_adj_indexes = np.array( ( np.where(np.min(np.abs( time_adj - timeRange[0] )) == np.abs( time_adj - timeRange[0] ) )[0][0] , \
                    np.where(np.min(np.abs( time_adj - timeRange[-1] )) == np.abs( time_adj - timeRange[-1] ) )[0][0] ) )
                ax[1].plot(time2/reportDivisor[0], sig2, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                ax[1].plot((time2[time_adj_indexes[0]:time_adj_indexes[1]+1])/reportDivisor[0], sig2[time_adj_indexes[0]:time_adj_indexes[1]+1], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                
                #--- Axis and Titles and Stuff ---
                if( FLG_fancyPlot_runner[figr] == 0 ): #only title non-fancy plot
                    string_Title = mode.capitalize()+' Method, Time Range of '+textNice(np.min(timeRange)/reportDivisor[0])+' to '+textNice(np.max(timeRange)/reportDivisor[0])+' '+reportDivisor[1]+'\n2nd Plot has Time Shift of '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1];
                    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
                #END IF
                ax[1].set_xlabel('Time ['+reportDivisor[1]+']',fontproperties=FONT_axisLabelFM); #set the x axis label
                ax[0].set_ylabel(plotName[:plotName.find(' & ')],fontproperties=FONT_axisLabelFM); #set the y axis label
                ax[1].set_ylabel(plotName[plotName.find(' & ')+3:],fontproperties=FONT_axisLabelFM); #set the y axis label
                #END IF
                
                if( time1.size != time2.size ):
                    #deal with some time issues (time2 can be longer than time1 for edges of time1)
                    if( timeShift[k] > 0 ):
                        #if timeShift positive it will use back-in-time time2 data so check the front edge of time1
                        if( (timeRange[0] - timeShift[k]) < time1[0] ):
                            ax[0].set_xlim( ((timeRange[0] - timeShift[k])/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                            ax[1].set_xlim( ((timeRange[0] - timeShift[k])/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                        else:
                            ax[0].set_xlim( (time1[0]/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                            ax[1].set_xlim( (time1[0]/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                        #END IF
                    elif( timeShift[k] < 0 ):
                        #if timeShift negative it will use forward-in-time time2 data so check the back edge of time1
                        if( (timeRange[-1] - timeShift[k]) > time1[-1] ):
                            ax[0].set_xlim( (time1[0]/reportDivisor[0], (timeRange[-1] - timeShift[k])/reportDivisor[0]) );
                            ax[1].set_xlim( (time1[0]/reportDivisor[0], (timeRange[-1] - timeShift[k])/reportDivisor[0]) );
                        else:
                            ax[0].set_xlim( (time1[0]/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                            ax[1].set_xlim( (time1[0]/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                        #END IF
                    else:
                        ax[0].set_xlim( (time1[0]/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                        ax[1].set_xlim( (time1[0]/reportDivisor[0], time1[-1]/reportDivisor[0]) );
                    #END IF
                #END IF
                
                #nice axis ticks
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],FLG_removeLabels=True,FLG_tickDirIn=True,FLG_manualLims=ax[0].get_xlim());
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True,FLG_manualLims=ax[1].get_xlim());
                else:
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=True,FLG_tickDirIn=True,FLG_manualLims=ax[0].get_xlim());
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=False,FLG_tickDirIn=True,FLG_manualLims=ax[1].get_xlim());
                #END IF
                
                figFitter(fig); #fit the fig fast
                if( FLG_fancyPlot_runner[figr] != 0 ):
                    fig.savefig(os.path.join(path_corrFancyPlot,'corr_'+mode+plotName_fileName+'_highlite'+'.png')); #save the figure
                    plt.close(); #close figure b/c it lurks apparently
                    plt.ion(); #re-enable it for later stuff
                #END IF
            #END FOR figr
        elif( (FLG_plot == True) & (k.size == 0) ):
            print('WARNING in subfun_correlator: No plots produced because no max found (likely due to one of the data inputs being all NaNs or times not aligning at all.');
        #END IF
        
    elif( mode == 'interval' ):
        #--- CALC INTERVALS on SIG1 ---
        #want when sig1 is increasing or decreasing or both depending on intervalType
        sig1_size = sig1.size; #calc it once
        sig1_smoothed = sig1; #no smooth
        # sig1_smoothed = subfun_filter(sig1, time1, 'Sav-Gol Smooth', {'savgol filter period':dataRate*intervalSmoother,'savgol filter order':1}, dataRate = dataRate); #smooth with savgol filter
        # intervalSmoother=20;
        # intervalNoiseRemover_fluctMax=12;
        # interval_diff = np.diff(sig1_smoothed)
        # interval_diff = np.insert(interval_diff,0,interval_diff[0]); #calc the diff of sig1
        # interval_diff[np.isclose(interval_diff,0)] = 0; #force close to 0 to be 0
        # # interval_diff_sign = np.sign(interval_diff); #get the sign of the diff
        # interval_mvAvg = uniform_filter1d(sig1_smoothed, size=intervalSmoother//2, mode='nearest');
        # interval_diff_normd = interval_diff*np.abs(interval_mvAvg); #normalize the diff
        # interval_diff_normd[np.isclose(interval_diff_normd,0)] = 0; #force close to 0 to be 0
        # interval_diff_normd_absmvAvg = uniform_filter1d(np.abs(interval_diff_normd), size=intervalSmoother*4, mode='constant');
        # interval_mvAvg_wide = uniform_filter1d(sig1_smoothed, size=intervalSmoother*4, mode='nearest');
        
        # interval_diff_normd_sign = np.sign(interval_diff_normd); #find signs of diff
        # interval_diff_normd_signChange = np.insert(np.diff(interval_diff_normd_sign),0,0); #find sign changes
        # possible_max_min = np.where(interval_diff_normd_signChange != 0)[0]; #get some stuff
        # possible_max_min_padded = np.append(np.insert(possible_max_min,0,0),sig1_size); #pad with 0 for start and size for end
        # possible_max_min_classify = np.zeros( possible_max_min.size, dtype=np.int32); #prep - 1 is max, -1 is min, 0 is nothing
        # for i in range(1, possible_max_min_padded.size-1):
        #     if( interval_diff_normd_sign[possible_max_min_padded[i]] == 0 ): #plateau detection
        #         #--- engage plateau scan ---
        #         j_prev = -1; #prep
        #         FLG_edge = False; #prep
        #         while( (interval_diff_normd_sign[possible_max_min_padded[i]+j_prev] == interval_diff_normd_sign[possible_max_min_padded[i]]) & (FLG_edge == False) ):
        #             if( i+j_prev == 0 ):
        #                 FLG_edge = True; #hit edge, note it
        #                 j_prev += 1; #offset comming increment
        #             #END IF
        #             j_prev += -1; #increment
        #         #END WHILE
        #         sign_prev = interval_diff_normd_sign[possible_max_min_padded[i]+j_prev]; #record it
                
        #         j_post = 1; #prep
        #         FLG_edge = False; #prep
        #         while( (interval_diff_normd_sign[possible_max_min_padded[i]+j_post] == interval_diff_normd_sign[possible_max_min_padded[i]]) & (FLG_edge == False) ):
        #             if( i+j_post == (sig1_size-1) ):
        #                 FLG_edge = True; #hit edge, note it
        #                 j_post += -1; #offset comming increment
        #             #END IF
        #             j_post += 1; #increment
        #         #END WHILE
        #         sign_post = interval_diff_normd_sign[possible_max_min_padded[i]+j_post]; #record it
                
        #         if( sign_prev > interval_diff_normd_sign[possible_max_min_padded[i]] > sign_post ):
        #             possible_max_min_classify[i-1] = 1; #classify as a max
        #             possible_max_min[i-1] += (possible_max_min_padded[i]+j_post-1-(possible_max_min_padded[i]+j_prev+1))//2; #if there's lots of 0's put the new center at the middle of the 0's
        #         elif( sign_prev < interval_diff_normd_sign[possible_max_min_padded[i]] < sign_post ):
        #             possible_max_min_classify[i-1] = -1; #classify as a min
        #             possible_max_min[i-1] += (possible_max_min_padded[i]+j_post-1-(possible_max_min_padded[i]+j_prev+1))//2; #if there's lots of 0's put the new center at the middle of the 0's
        #         #END IF
        #     elif( interval_diff_normd_sign[possible_max_min_padded[i]] == -1 ): #positive detection
        #         if( (interval_diff_normd_sign[possible_max_min_padded[i]-1] == 1) & (-1 == interval_diff_normd_sign[possible_max_min_padded[i]]) ):
        #             possible_max_min_classify[i-1] = 1; #classify as a max
        #         #END IF
        #     else: #negative detection
        #         if( (interval_diff_normd_sign[possible_max_min_padded[i]-1] == -1) & (1 == interval_diff_normd_sign[possible_max_min_padded[i]]) ):
        #             possible_max_min_classify[i-1] = -1; #classify as a min
        #         #END IF
        #     #END IF
        # #END FOR i
        # possible_max_min_classify_max = possible_max_min_classify == 1;
        # possible_max_min_classify_min = possible_max_min_classify == -1;        
        
        # plt.figure()
        # plt.plot(np.arange(0,sig1.size,1),sig1_smoothed)
        # # plt.plot(np.arange(0,sig1.size,1),sig1)
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_max]-1],sig1_smoothed[possible_max_min[possible_max_min_classify_max]-1],c='xkcd:red',s=50)
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_min]-1],sig1_smoothed[possible_max_min[possible_max_min_classify_min]-1],c='xkcd:green',s=50)
        
        #--- Remove max/mins that are small fluctuations ---
        # from scipy import interpolate
        from scipy.optimize import curve_fit, minimize
        from scipy.stats import linregress

        
        #declare function and constraints
        def funcFirst(x, b, m): #this one's for curve_fit b/c b/c
            return m*x + b
        #END DEF
        # def func(x,params):
        #     return params[1]*x + params[0]
        # #END DEF       
        # # def constraint_strt(params):
        # #     return y[0] - func(x[0],params)
        # # #END DEF
        # def constraint_strtStop(params, x, y):
        #     return y - func(x,params)
        # #END DEF
        # def goal_lstSqrs(params, x, y):
        #     y_try = func(x,params);
        #     return np.sum((y_try - y)**2) #ls sq error goal
        # #END DEF
        # constraints = [{'type':'eq', 'fun': constraint_strt},
        #                {'type':'eq', 'fun': constraint_stop}];
        # constraints_noStrt = [{'type':'eq', 'fun': constraint_stop, 'args':[x,y]}]; #use for 1st b/c we have nothing to connect to
        
        x = np.arange(0,sig1_size,1);
        stop_min = np.int64(np.round(sig1_smoothed.size*intervalStepScaler/100)); #declare min step size
        if( stop_min < 2 ):
            stop_min = 2; #prep
        #END IF
        step_size = np.int64(np.round(sig1_smoothed.size*intervalStepScaler/1000)); #declare min step size (intervalStepScaler%/10)
        step_minNum = np.int64(np.round(sig1_smoothed.size/step_size*intervalStepScaler/100)); #min number of steps to take
        strt = 0; #set
        stop = stop_min; #set
        fit_params = []; #prep
        fit = []; #prep
        rSq = []; #prep
        x_vals = []; #prep
        y_vals = []; #prep
        time_vals = []; #prep (for later)
        firstGuess = curve_fit(funcFirst, x, sig1_smoothed)[0];
        cntr_wide = 0; #prep
        while( (strt < sig1_size) & ((stop-strt) > 1) ):
            #--- while loop to find best rSq ---
            # if( strt != 0 ):
            #     constraints_toUse = [{'type':'eq', 'fun': constraint_stop, 'args':[x_vals[cntr_wide-1][-1],fit[cntr_wide-1][-1]]}]; #use for 1st b/c we have nothing to connect to
            # else:
            #     constraints_toUse = []; #this is special for 1st b/c nothing to align to 
            # #END IF
            fit_params_try = [];
            fit_try = [];
            rSq_try = []; #prep list
            cntr = 0; #prep
            FLG_eject = False;
            stop_try = [stop];
            while( ((cntr < step_minNum) | (FLG_eject == False)) & (stop_try[cntr] <= sig1_size) ):  
                fit_params_try.append(curve_fit(funcFirst, x[strt:stop_try[cntr]], sig1_smoothed[strt:stop_try[cntr]])[0]);
                #END TRY
                # if( cntr == 0 ):
                #     firstGuess = curve_fit(funcFirst, x, sig1_smoothed)[0]; #get unconstrained 1st guess
                # #END IF
                # fit_params_try.append(minimize(goal_lstSqrs, x0=firstGuess, args=(x[strt:stop_try[cntr]],sig1_smoothed[strt:stop_try[cntr]]), constraints=constraints_toUse)['x']);                     
                fit_try.append(fit_params_try[cntr][1]*x[strt:stop_try[cntr]] + fit_params_try[cntr][0]);
                # fit_try.append(func(x[strt:stop_try[cntr]],fit_params_try[cntr]));
                
                _, _, r_val, _, _ = linregress(sig1_smoothed[strt:stop_try[cntr]], fit_try[cntr]);
                rSq_try.append(r_val**2);
                if( rSq_try[cntr-1]*intervalFitScaler < np.max(rSq_try) ):
                    FLG_eject = True; #eject
                #END IF
                cntr += 1; #increment
                stop_try.append(stop_try[cntr-1]+step_size); #increment
            #END WHILE
            k = np.where( rSq_try == np.max(rSq_try))[0][0]; #get the maximal
            fit_params.append(fit_params_try[k]); #record
            fit.append(fit_try[k]); #record     
            rSq.append(rSq_try[k]); #record
            x_vals.append(x[strt:stop_try[k]]);
            y_vals.append(sig1_smoothed[strt:stop_try[k]]);
            time_vals.append(time1[strt:stop_try[k]]);
            strt = stop_try[k]; #move up the start
            stop = strt + stop_min; #add on the minimum stop range
            
            if( stop > sig1_size ):
                stop = sig1_size; #cap it
            #END IF
            cntr_wide += 1; #Increment
        #END WHILE
        x_vals_strtStop = np.zeros( (len(fit), 2), dtype=np.int64 );
        for i in range(0,len(fit)):
            x_vals_strtStop[i,0] = x_vals[i][0];
            x_vals_strtStop[i,1] = x_vals[i][-1];
        #END FOR i
        
        # plt.figure()
        # plt.plot(x,np.abs(interval_diff_normd))
        # plt.scatter(x[possible_max_min[possible_max_min_classify_max]-1],np.abs(interval_diff_normd[possible_max_min[possible_max_min_classify_max]-1]),c='xkcd:red',s=50)
        # plt.scatter(x[possible_max_min[possible_max_min_classify_min]-1],np.abs(interval_diff_normd[possible_max_min[possible_max_min_classify_min]-1]),c='xkcd:green',s=50)
        
        # plt.figure()
        # plt.plot(np.arange(0,sig1.size,1),sig1_smoothed)
        # for i in range(0,len(fit)):
        #     plt.plot(x_vals[i],fit[i])
        # #END FOR i
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_max]-1],sig1_smoothed[possible_max_min[possible_max_min_classify_max]-1],c='xkcd:red',s=50)
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_min]-1],sig1_smoothed[possible_max_min[possible_max_min_classify_min]-1],c='xkcd:green',s=50)
        
        def func(x_vals,fit_params):
            y_try = [None for i in range(0,len(fit_params))]; #preallocate
            for i in range(0,len(fit_params)):
                y_try[i] = fit_params[i][1]*x_vals[i] + fit_params[i][0]; #calc the values as we gooo
            #END FOR i
            return y_try
        #END DEF  
        def funcSingle(x,params):
            return params[1]*x + params[0]
        #END DEF 
        def goal_lstSqrs_list(fit_params_flat, x_vals, y_vals, sizer_per, intervalFitEndEncourager):
            fit_params = conv_flat2list(fit_params_flat, sizer_per); #make into a list
            y_try = func(x_vals, fit_params); #get the deal
            errLSSQ = 0; #prep
            for i in range(0,len(fit_params)):
                errLSSQ += np.sum((y_try[i] - y_vals[i])**2); #ls sq error goal
                #--- bonus, make ends meet ---
                if( (i != 0) & (i != len(fit_params)-1) ):
                    errLSSQ += ((funcSingle(x_vals[i-1][-1], fit_params[i-1]) - funcSingle(x_vals[i][0], fit_params[i]))**2 + \
                                (funcSingle(x_vals[i+1][0], fit_params[i+1]) - funcSingle(x_vals[i][-1], fit_params[i]))**2)*intervalFitEndEncourager; #encourage ends to meet
                elif( i == 0 ):
                    errLSSQ += ((funcSingle(x_vals[i+1][0], fit_params[i+1]) - funcSingle(x_vals[i][-1], fit_params[i]))**2)*intervalFitEndEncourager; #encourage ends to meet
                else:
                    errLSSQ += ((funcSingle(x_vals[i-1][-1], fit_params[i-1]) - funcSingle(x_vals[i][0], fit_params[i]))**2)*intervalFitEndEncourager; #encourage ends to meet
                #END IF
            #END FOR i
            return errLSSQ
        #END DEF
        #--- these are b/c we can't minimize in lists ---
        def conv_list2flat(fit_params):
            #convert to a flat array so basic functions can work with it
            sizer = 0; #preallocate
            sizer_per = np.empty(len(fit_params),dtype=np.int64); #preallocate
            for i in range(0,len(fit_params)):
                sizer += fit_params[i].size; #get size of fit_params
                sizer_per[i] = fit_params[i].size; #record individual sizes
            #END FOR i
            fit_params_flat = np.empty(sizer); #preallocate
            cntr = 0; #prep
            for i in range(0,len(fit_params)):
                fit_params_flat[cntr:cntr+sizer_per[i]] = fit_params[i]; #pull out the stuff
                cntr += sizer_per[i]; #increment
            #END FOR i
            return fit_params_flat, sizer_per
        #END DEF
        def conv_flat2list(fit_params_flat, sizer_per):
            fit_params = [None for i in range(0,sizer_per.size)]; #create unholy combo
            cntr = 0; #prep
            for i in range(0,sizer_per.size):
                fit_params[i] = fit_params_flat[cntr:cntr+sizer_per[i]]; #put it in
                cntr += sizer_per[i]; #increment
            #END FOR i
            return fit_params
        #END DEF
        #constraints don't work and make nothing happen
        # def constraint_strtStop(fit_params_flat, x_vals, sizer_per, idx1, idx1_idx, idx2, idx2_idx):
        #     fit_params = conv_flat2list(fit_params_flat, sizer_per); #convert to list
        #     return funcSingle(x_vals[idx1][idx1_idx],fit_params[idx1]) - funcSingle(x_vals[idx2][idx2_idx],fit_params[idx2])
        # #END DEF
        # #--- build constraints ---
        # constraints_toUse = []; #prep list
        # for i in range(0,len(fit_params)):
        #     if( (i != 0) & (i != len(fit_params)-1) ):
        #         constraints_toUse.append({'type':'eq', 'fun': constraint_strtStop, 'args': (x_vals, sizer_per, i, 0, i-1, -1)}); #append a start limit
        #         constraints_toUse.append({'type':'eq', 'fun': constraint_strtStop, 'args': (x_vals, sizer_per, i, -1, i+1, 0)}); #append an end limit
        #     elif( i == 0 ):
        #         constraints_toUse.append({'type':'eq', 'fun': constraint_strtStop, 'args': (x_vals, sizer_per, i, -1, i+1, 0)}); #append an end limit
        #     else:
        #         constraints_toUse.append({'type':'eq', 'fun': constraint_strtStop, 'args': (x_vals, sizer_per, i, 0, i-1, -1)}); #append a start limit
        #     #END IF
        # #END FOR i
        
        
        fit_params_flat, sizer_per = conv_list2flat(fit_params);
        fit_params_mecha_flat = minimize(goal_lstSqrs_list, x0=fit_params_flat, args=(x_vals, y_vals, sizer_per, intervalFitEndEncourager) )['x'];
        # fit_params_mecha_flat = minimize(goal_lstSqrs_list, x0=fit_params_flat, args=(x_vals, y_vals, sizer_per), constraints=constraints_toUse )['x'];
        fit_params_mecha = conv_flat2list(fit_params_mecha_flat, sizer_per); #make into a list
        fit_mecha = func(x_vals,fit_params_mecha); #calc fit vals
        
        # plt.figure()
        # plt.plot(np.arange(0,sig1.size,1),sig1_smoothed)
        # for i in range(0,len(fit_mecha)):
        #     plt.plot(x_vals[i],fit_mecha[i])
        # #END FOR i
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_max]-1],sig1_smoothed[possible_max_min[possible_max_min_classify_max]-1],c='xkcd:red',s=50)
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_min]-1],sig1_smoothed[possible_max_min[possible_max_min_classify_min]-1],c='xkcd:green',s=50)
        
        rSq_mecha = np.empty( len(fit_mecha) );
        fit_slope_mecha = np.empty( len(fit_mecha) );
        for i in range(0,len(fit_mecha)):
            _, _, r_val, _, _ = linregress(y_vals[i], fit_mecha[i]);
            rSq_mecha[i] = r_val**2;
            fit_slope_mecha[i] = fit_params_mecha[i][1]; #pull out the slope only
        #END FOR i

        # plt.figure()
        # plt.plot(np.arange(0,sig1.size,1),sig1)
        # for i in range(0,len(fit_mecha)):
        #     plt.plot(x_vals[i],fit_mecha[i])
        # #END FOR i
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_max]-1],sig1[possible_max_min[possible_max_min_classify_max]-1],c='xkcd:red',s=50)
        # plt.scatter(np.arange(0,sig1.size,1)[possible_max_min[possible_max_min_classify_min]-1],sig1[possible_max_min[possible_max_min_classify_min]-1],c='xkcd:green',s=50)

        
        # intervalNoiseRemover_fluctMax = intervalNoiseRemover_fluctMax/100; #convert from percent to decimal
        # max_min = possible_max_min[possible_max_min_classify != 0]-1; #get the actual max_mins
        # max_min_classify = possible_max_min_classify[possible_max_min_classify != 0]; #get the classified actual max_mins as well
        # max_min_padded = np.append(np.insert(max_min,0,0),sig1_size-1); #pad with 0 for start and size for end
        # pts_inbetween = np.diff(possible_max_min_padded);
        # sig1_smoothed_extra = uniform_filter1d(sig1_smoothed, size=np.int64(np.round(intervalSmoother*intervalSmootherScaler)), mode='nearest');
        # interval_diff_normd_smoothed = uniform_filter1d(interval_diff_normd, size=np.int64(np.round(intervalSmoother*4)), mode='nearest');
        # # max_min_keep_padded = np.ones( max_min.size+2, dtype=np.bool_); #preallocate keep array
        # FLG_remove = False; #prep
        # slopeLogic = False; #reset logic
        # cntr = 1; #prep hwile cntr
        # while( cntr < max_min_padded.size-1 ):
        #     if( max_min_padded[cntr] == 911 ):
        #         sys.crash()
        #     #END IF
        #     #--- mean test ---
        #     val_prev = sig1_smoothed[max_min_padded[cntr-1]];
        #     val_curr = sig1_smoothed[max_min_padded[cntr]];
        #     val_post = sig1_smoothed[max_min_padded[cntr+1]];
        #     mean_absMean = np.abs(val_prev+val_curr+val_post)/3; #get total absolute mean
        #     # slope_curr = (val_curr-val_prev)/(max_min_padded[cntr]-max_min_padded[cntr-1]);
        #     slope_curr = (sig1_smoothed[max_min_padded[cntr+1]]-sig1_smoothed[max_min_padded[cntr-1]])/(max_min_padded[cntr+1]-max_min_padded[cntr-1]);
        #     k = np.where( (x_vals_strtStop[:,0] < max_min_padded[cntr]) & (x_vals_strtStop[:,1] >= max_min_padded[cntr]) )[0].item();
        #     if( max_min_classify[cntr-1] == 1 ):
        #         #it's a max
        #         slope_max = np.max(sig1_smoothed[x_vals_strtStop[k,0]:x_vals_strtStop[k,1]+1]); #get the max in the range
        #         if( slope_max == val_curr ):
        #             slopeLogic = True; #keep based being the maximum in the range
        #         #END IF
        #     else:
        #         #it's a min
        #         slope_min = np.min(sig1_smoothed[x_vals_strtStop[k,0]:x_vals_strtStop[k,1]+1]); #get the max in the range
        #         if( slope_min == val_curr ):
        #             slopeLogic = True; #keep based being the maximum in the range
        #         #END IF
        #     #END IF
        #     slope_currFitted = fit_params[np.where( (x_vals_strtStop[:,0] < max_min_padded[cntr]) & (x_vals_strtStop[:,1] >= max_min_padded[cntr]) )[0].item()][0]; #get the fitted slope from above
        #     if( ~((np.abs(val_curr-val_prev)/mean_absMean >= intervalNoiseRemover_fluctMax) & \
        #           (np.abs(val_curr-val_post)/mean_absMean >= intervalNoiseRemover_fluctMax)) & \
        #         (slopeLogic == True) ):
        #         cntr += 1; #increment temporarially
        #         if( (cntr+1) < max_min_padded.size-1 ): #don't go over
        #             val_prev = sig1_smoothed[max_min_padded[cntr-1]];
        #             val_curr = sig1_smoothed[max_min_padded[cntr]];
        #             val_post = sig1_smoothed[max_min_padded[cntr+1]];
        #             mean_absMean = np.abs(val_prev+val_curr+val_post)/3; #get total mean
        #             # slope_currFwd = (val_curr-val_prev)/(max_min_padded[cntr]-max_min_padded[cntr-1]);
        #             slope_currFwd = (sig1_smoothed[max_min_padded[cntr+1]]-sig1_smoothed[max_min_padded[cntr-1]])/(max_min_padded[cntr+1]-max_min_padded[cntr-1]);
        #             slope_currFwdFitted = fit_params[np.where( (x_vals_strtStop[:,0] < max_min_padded[cntr]) & (x_vals_strtStop[:,1] >= max_min_padded[cntr]) )[0].item()][0]; #get the fitted slope from above
        #             # min_max_curr = max_min_classify[cntr-1]; #get min or max
        #             if( ~((np.abs(val_curr-val_prev)/mean_absMean >= intervalNoiseRemover_fluctMax) & \
        #                   (np.abs(val_curr-val_post)/mean_absMean >= intervalNoiseRemover_fluctMax)) | \
        #                 (slope_currFwd<slope_currFwdFitted) ):
        #                 #small fluctuations are removed
        #                 #slope is used to detect ~temporary~ plateaus
        #                 FLG_remove = True;
        #             #END IF
        #         #END IF
        #         cntr += -1; #remove the testing increment
        #     #END IF
        #     # #--- even farther mean test ---
        #     # if( FLG_doubleCheck == True ):
        #     #     cntr_fwd = cntr + 2; #check next max or min
        #     #     if( cntr_fwd < max_min_padded.size-1 ): #don't go over
        #     #         val_prev = sig1_smoothed[max_min_padded[cntr-1]];
        #     #         val_curr = sig1_smoothed[max_min_padded[cntr]];
        #     #         slope_prev2curr = (val_curr-val_prev)/(max_min_padded[cntr]-max_min_padded[cntr-1]);
        #     #         val_prev = sig1_smoothed[max_min_padded[cntr_fwd-1]];
        #     #         val_curr = sig1_smoothed[max_min_padded[cntr_fwd]];
        #     #         slope_prev2curr_fwd = (val_curr-val_prev)/(max_min_padded[cntr_fwd]-max_min_padded[cntr_fwd-1]);
        #     #         if( ~(((max_min_classify[cntr-1] == 1) & (slope_prev2curr > slope_prev2curr_fwd)) | \
        #     #              ((max_min_classify[cntr-1] == -1) & (slope_prev2curr < slope_prev2curr_fwd))) ):
        #     #             #middle-increases that aren't caught before are removed
        #     #             FLG_remove = True;
        #     #         #END IF
        #     #     #END IF
        #     #     FLG_doubleCheck = False; #reset flag
        #     # #END IF
            
            
        #     # #--- slope test ---
        #     # slope_prev = interval_diff_normd[max_min_padded[cntr-1]];
        #     # slope_prev_mean = np.mean(interval_diff_normd[max_min_padded[cntr-1]:max_min_padded[cntr]])
        #     # slope_curr = interval_diff_normd[max_min_padded[cntr]];
        #     # slope_curr_smoothed = interval_diff_normd_smoothed[max_min_padded[cntr]];
        #     # slope_post = interval_diff_normd[max_min_padded[cntr+1]];
        #     # slope_post_mean = np.mean(interval_diff_normd[max_min_padded[cntr]+1:max_min_padded[cntr+1]+1])
        #     # if( np.abs(interval_diff_normd[max_min_padded[cntr]])/np.abs(interval_diff_normd_smoothed[max_min_padded[cntr]]) >= intervalNoiseRemover_fluctMax*4 ):
        #     #     cntr += -1; #keep on truckin' by reverting the increment
        #     #     FLG_remove = False; #turn off flag in case it was turned on previously
        #     # else:
        #     #     #instantaneous slope is very different from running-mean slope which means its likely spurrious
        #     #     FLG_remove = True;
        #     # #END IF
            
        #     #--- remove alg ---
        #     if( FLG_remove == True ):
        #         max_min_padded = np.delete(max_min_padded,(cntr,cntr+1)); #remove
        #         max_min_classify = np.delete(max_min_classify,(cntr-1,cntr)); #remove
        #         # max_min_keep[cntr-2] = False; #set 2nd pt to reject b/c min/max comes in pairs
        #         # max_min_keep[cntr-1] = False; #set 1st pt to reject
        #         cntr += -1; #revert to do the same cntr # again b/c we just removed 2 pts from the set
        #         FLG_remove = False; #reset flag
        #     #END IF
        #     cntr += 1; #increment to next
        #     slopeLogic = False; #reset logic
        # #END WHILE
        
        # max_min_final = max_min_padded[1:-1]; #keep just the stuff that's good
        # max_min_final_classify = max_min_classify; #keep just the stuff that's good
        
        
        # max_min_final_classify_max = max_min_final_classify == 1;
        # max_min_final_classify_min = max_min_final_classify == -1;    
        
        # plt.figure()
        # plt.plot(np.arange(0,sig1.size,1),sig1_smoothed)
        # # plt.plot(np.arange(0,sig1.size,1),sig1)
        # plt.scatter(np.arange(0,sig1.size,1)[max_min_final[max_min_final_classify_max]],sig1_smoothed[max_min_final[max_min_final_classify_max]],c='xkcd:red',s=50)
        # plt.scatter(np.arange(0,sig1.size,1)[max_min_final[max_min_final_classify_min]],sig1_smoothed[max_min_final[max_min_final_classify_min]],c='xkcd:green',s=50)
        
        #--- find time interval req'd ---
        #identify positives
        if( intervalType == 'pos' ):
            kj = np.where((rSq_mecha > intervalFitMinRSq) & (fit_slope_mecha > 0))[0]; #get positive slopes with a good enough confidence
            kk = np.zeros( kj.size, dtype=np.bool_ );
            for i in range(0,kj.size):
                if( timeInterval <= dataRate*x_vals[kj[i]].size ):
                    kk[i] = True; #keep it
                #END FOR i
            #END FOR i
            kj = kj[kk]; #keep only the right stuff
        elif( intervalType == 'neg' ):
            kj = np.where((rSq_mecha > intervalFitMinRSq) & (fit_slope_mecha < 0))[0]; #get negative slopes with a good enough confidence
            kk = np.zeros( kj.size, dtype=np.bool_ );
            for i in range(0,kj.size):
                if( timeInterval <= dataRate*x_vals[kj[i]].size ):
                    kk[i] = True; #keep it
                #END FOR i
            #END FOR i
            kj = kj[kk]; #keep only the right stuff
        else: #otherwise both
            kj = np.where(rSq_mecha > intervalFitMinRSq)[0];
            kk = np.zeros( kj.size, dtype=np.bool_ );
            for i in range(0,kj.size):
                if( timeInterval <= dataRate*x_vals[kj[i]].size ):
                    kk[i] = True; #keep it
                #END FOR i
            #END FOR i
            kj = kj[kk]; #keep only the right stuff
        #END IF
        
        sizer = 0; #prep
        for i in range(0,kj.size):
            sizer += x_vals[kj[i]].size; #get the total size
        #END FOR i
        sig1_only = np.empty( sizer ); #preallocate
        # time1_only = np.empty( sizer ); #preallocate
        cntr = 0; #prep cntr
        for i in range(0,kj.size):
            sig1_only[cntr:cntr+x_vals[kj[i]].size] = y_vals[kj[i]]; #get the sig1
            # time1_only[cntr:cntr+x_vals[kj[i]].size] = time_vals[kj[i]]; #get the time1
            cntr += x_vals[kj[i]].size; #increment
        #END FOR i
        
        #--- ROLLING TIME shift TO FIND OPTIMAL CORR COEFF ---
        corr = np.zeros(timeShift.size); #preallocate
        #CALC PWRS
        pwr_sig1_only = np.sqrt(1/sig1_only.size*np.sum(sig1_only**2)); #estimate power of signal
        sig1_only_normd = 1/pwr_sig1_only*sig1_only; #norm it once b/c doesn't change
        for i in range(0,timeShift.size):
            #--- Adjust by the time offset ---
            sig2_adj = np.zeros( sizer ); #preallocate
            cntr = 0; #prep cntr
            for j in range(0,kj.size):
                # k_strt = np.where( np.min(np.abs(time2 - (time_vals[kj[j]][0] + timeShift[i]))) == np.abs(time2 - (time_vals[kj[j]][0] + timeShift[i])) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                # k_end = np.where( np.min(np.abs(time2 - (time_vals[kj[j]][-1] + timeShift[i]))) == np.abs(time2 - (time_vals[kj[j]][-1] + timeShift[i])) )[0][0];
                k_strt = np.where( np.min(np.abs( (time2+timeShift[i]) - time_vals[kj[j]][0])) == np.abs( (time2+timeShift[i]) - time_vals[kj[j]][0]) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                k_end = np.where( np.min(np.abs( (time2+timeShift[i]) - time_vals[kj[j]][-1])) == np.abs( (time2+timeShift[i]) - time_vals[kj[j]][-1]) )[0][0];
                if( (k_strt.size != 0) & (k_end.size != 0) & (k_strt != k_end) & (x_vals[kj[j]].size == (k_end-k_strt+1)) ):
                    sig2_adj[cntr:cntr+x_vals[kj[j]].size] = sig2[k_strt.item():k_end.item()+1]; #get the sig2 adjusted
                else:
                    sig2_adj[cntr:cntr+x_vals[kj[j]].size] = np.nan; #set to NaN
                #END IF
                cntr += x_vals[kj[j]].size; #increment
            #END FOR j
            jk = ~np.isnan(sig2_adj);
            if( jk.sum() > 0 ):
                sig2_adj = sig2_adj[jk]; #remove NaNs from shifts that were out of range
                
                if( jk.sum() == sig1_only_normd.size ):
                    sig1_only_normd_temp = sig1_only_normd; #keep it real
                else:
                    sig1_only_normd_temp = sig1_only_normd[jk]; #also remove this data to match size and time
                #END IF
                
                #CALC PWRS
                pwr_sig2_adj = np.sqrt(1/sig2_adj.size*np.sum(sig2_adj**2)); #estimate power of signal
                
                #Real quick side move to calc correlation coefficients
                corr[i] = np.corrcoef(sig1_only_normd_temp,1/pwr_sig2_adj*sig2_adj)[0,1];
            else:
                corr[i] = np.nan; #if timeshift is totally out of range and there's no data at all - return nan
            #END IF
        #END FOR i
        k = np.where( np.abs(corr) == np.nanmax(np.abs(corr)))[0]; #get maximum values
        if( intervalType == 'pos' ):
            string_intervalType = 'positive';
        elif( intervalType == 'neg' ):
            string_intervalType = 'negative';
        else:
            string_intervalType = 'pos+neg';
        #END IF
        corr_textResults = 'Corr coeff max: '+textNice(np.round(corr[k],reportRounder))+' at time shift: '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1]+' for any time interval that is '+string_intervalType+' for longer than '+textNice(timeInterval/reportDivisor[0])+' '+reportDivisor[1];
        if( FLG_enableText ):
            print(corr_textResults);#print resupts
        #END IF
        if( k.size > 0 ):
            corrRet['corr'] = corr[k]; #record
            corrRet['time shift'] = timeShift[k]; #record
        else:
            corrRet['corr'] = np.array((np.nan,));
            corrRet['time shift'] = np.array((np.nan,));
        #END IF
        corrRet['time interval'] = timeInterval; #record
        corrRet['time interval type'] = intervalType; #record
        corrRet['data rate'] = dataRate; #record
        corrRet['corr vect'] = corr; #record
        corrRet['time shift vect'] = timeShift; #record
        corrRet['plot name'] = plotName; #record
        corrRet['text results'] = corr_textResults; #record
        
        if( (FLG_plot == True) & (k.size > 0) ):
            #--- Declare & Unpack ---
            if( (FLG_fancyPlot >= 1) & FLG_enableText ):
                print('MAKING FANCY PLOT: corr_'+mode+plotName_fileName+'_highlite'+' IN fancyPlot FOLDER'); #report since you won't see anything
            #END IF
            # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
            PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
            PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
            # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
            # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
            # PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
            # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
            FONT_titleFM = settings_plot['font title FM'];
            # FONT_axisTick = settings_plot['font axis tick'];
            FONT_axisLabelFM = settings_plot['font axis label FM'];
            
            if( FLG_fancyPlot == 0 ): #make a loop that can do fancy and non-fancy in one (correlation calcs are expensive when done en masse)
                FLG_fancyPlot_runner = np.array( (0,) );
            elif( FLG_fancyPlot == 1 ):
                FLG_fancyPlot_runner = np.array( (0,1) );
            elif( FLG_fancyPlot == 2 ):
                FLG_fancyPlot_runner = np.array( (1,) );
            #END IF
            
            #--- Prep Plot ---
            for figr in range(0, FLG_fancyPlot_runner.size):
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
                    figManager = fig.canvas.manager; #req to maximize
                    figManager.window.showMaximized(); #force maximized
                else:
                    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                    fig, ax = plt.subplots(nrows=2, ncols=1,figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
                #END IF
                fig.subplots_adjust(hspace=0);
                # divider = make_axes_locatable(ax); #prep to add an axis
                # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
                
                #Remove the aspect ratio from the basemap so it fills the screen better
                ax[0].set_aspect('auto');
                ax[1].set_aspect('auto');
                
                #--- Actual Plotting ---
                ax[0].plot(time1/reportDivisor[0], sig1, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                for i in range(0,kj.size):
                    ax[0].plot(time_vals[kj[i]]/reportDivisor[0], y_vals[kj[i]], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                #END FOR i
                for i in range(0,len(fit_mecha)):
                    ax[0].plot(time_vals[i]/reportDivisor[0],fit_mecha[i],linewidth=PLOT_lineWidthPlus,linestyle='--',zorder=3);
                #END FOR i
                
                ax[1].plot(time2/reportDivisor[0], sig2, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                for i in range(0,kj.size):
                    # k_strt = np.where( np.min(np.abs(time2 - (time_vals[kj[i]][0] + timeShift[k]))) == np.abs(time2 - (time_vals[kj[i]][0] + timeShift[k])) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                    # k_end = np.where( np.min(np.abs(time2 - (time_vals[kj[i]][-1] + timeShift[k]))) == np.abs(time2 - (time_vals[kj[i]][-1] + timeShift[k])) )[0][0];
                    k_strt = np.where( np.min(np.abs( (time2 + timeShift[k]) - time_vals[kj[i]][0])) == np.abs( (time2 + timeShift[k]) - time_vals[kj[i]][0]) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                    k_end = np.where( np.min(np.abs( (time2 + timeShift[k]) - time_vals[kj[i]][-1])) == np.abs( (time2 + timeShift[k]) - time_vals[kj[i]][-1]) )[0][0];
                    if( (k_strt.size != 0) & (k_end.size != 0) & (k_strt != k_end) & (x_vals[kj[i]].size == (k_end-k_strt+1)) ):
                        ax[1].plot((time2[k_strt.item():k_end.item()+1])/reportDivisor[0], sig2[k_strt.item():k_end.item()+1], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                    #END IF
                #END FOR i
                
                #--- Axis and Titles and Stuff ---
                if( FLG_fancyPlot_runner[figr] == 0 ): #only title non-fancy plot
                    if( intervalType == 'both' ):
                        intervalType = 'Pos. & Neg';
                    #END IF
                    string_Title = mode.capitalize()+' Method, '+intervalType.capitalize()+'. Intervals Longer than '+textNice(timeInterval/reportDivisor[0])+' '+reportDivisor[1]+' Highlighted\n2nd Plot has Time Shift of '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1];
                    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
                #END IF
                ax[1].set_xlabel('Time ['+reportDivisor[1]+']',fontproperties=FONT_axisLabelFM); #set the x axis label
                ax[0].set_ylabel(plotName[:plotName.find(' & ')],fontproperties=FONT_axisLabelFM); #set the y axis label
                ax[1].set_ylabel(plotName[plotName.find(' & ')+3:],fontproperties=FONT_axisLabelFM); #set the y axis label
                #END IF
                
                #nice axis ticks
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],FLG_removeLabels=True,FLG_tickDirIn=True);
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True);
                else:
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=True,FLG_tickDirIn=True);
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=False,FLG_tickDirIn=True);
                #END IF
    
                
                figFitter(fig); #fit the fig fast
                if( FLG_fancyPlot_runner[figr] != 0 ):
                    fig.savefig(os.path.join(path_corrFancyPlot,'corr_'+mode+plotName_fileName+'_highlite'+'.png')); #save the figure
                    plt.close(); #close figure b/c it lurks apparently
                    plt.ion(); #re-enable it for later stuff
                #END IF
            #END FOR figr
        elif( (FLG_plot == True) & (k.size == 0) ):
            print('WARNING in subfun_correlator: No plots produced because no max found (likely due to one of the data inputs being all NaNs or times not aligning at all.');
        #END IF
        
    elif( mode == 'interval manual' ):
        if( type(timeInterval) == type([]) ):
            timeInterval = np.int64(np.round(np.array(timeInterval))); #convert, round, and force to be an integer
        #END IF
        
        #--- ROLLING TIME shift TO FIND OPTIMAL CORR COEFF ---
        corr = np.zeros(timeShift.size); #preallocate
        
        #GET SIG1 ONLY WHERE MATTERS
        sig1_holder = []; #prep
        sig1_holderSize = []; #prep
        sig1_indexes = []; #prep
        for j in range(0,len(timeInterval)):
            k_strt = np.where( np.min(np.abs( time1 - timeInterval[j,0])) == np.abs( time1 - timeInterval[j,0]) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
            k_end = np.where( np.min(np.abs( time1 - timeInterval[j,1])) == np.abs( time1 - timeInterval[j,1]) )[0][0];
            sig1_holder.append(sig1[k_strt.item():k_end.item()+1]); #get the sig1 bits
            sig1_holderSize.append(k_end.item()-k_strt.item()+1); #get the sig1 bit size
            sig1_indexes.append([k_strt.item(),k_end.item()+1]); #get the indexes
        #END FOR j
        sig1_indexes = np.array(sig1_indexes); #convert to useful array
        sig1_only = np.zeros( np.sum(sig1_holderSize) ); #prep
        cntr = 0; #prep cntr
        for j in range(0,len(timeInterval)):
            sig1_only[cntr:cntr+sig1_holderSize[j]] = sig1_holder[j]; #get the sig1
            cntr += sig1_holderSize[j]; #increment
        #END IF
        
        #CALC PWRS
        pwr_sig1_only = np.sqrt(1/sig1_only.size*np.sum(sig1_only**2)); #estimate power of signal
        sig1_only_normd = 1/pwr_sig1_only*sig1_only; #norm it once b/c doesn't change
        for i in range(0,timeShift.size):
            #--- Adjust by the time offset ---
            sig2_adj = np.zeros( np.sum(sig1_holderSize) ); #preallocate
            cntr = 0; #prep cntr
            for j in range(0,len(timeInterval)):
                k_strt = np.where( np.min(np.abs( (time2+timeShift[i]) - timeInterval[j,0])) == np.abs( (time2+timeShift[i]) - timeInterval[j,0]) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                k_end = np.where( np.min(np.abs( (time2+timeShift[i]) - timeInterval[j,1])) == np.abs( (time2+timeShift[i]) - timeInterval[j,1]) )[0][0];
                if( (k_strt.size != 0) & (k_end.size != 0) & (k_strt != k_end) & (sig1_holderSize[j] == (k_end-k_strt+1)) ):
                    sig2_adj[cntr:cntr+sig1_holderSize[j]] = sig2[k_strt.item():k_end.item()+1]; #get the sig2 adjusted
                else:
                    sig2_adj[cntr:cntr+sig1_holderSize[j]] = np.nan; #set to NaN
                #END IF
                cntr += sig1_holderSize[j]; #increment
            #END FOR j
            jk = ~np.isnan(sig2_adj);
            if( jk.sum() > 0 ):
                sig2_adj = sig2_adj[jk]; #remove NaNs from shifts that were out of range
                
                if( jk.sum() == sig1_only_normd.size ):
                    sig1_only_normd_temp = sig1_only_normd; #keep it real
                else:
                    sig1_only_normd_temp = sig1_only_normd[jk]; #also remove this data to match size and time
                #END IF
                
                #CALC PWRS
                pwr_sig2_adj = np.sqrt(1/sig2_adj.size*np.sum(sig2_adj**2)); #estimate power of signal
                
                #Real quick side move to calc correlation coefficients
                corr[i] = np.corrcoef(sig1_only_normd_temp,1/pwr_sig2_adj*sig2_adj)[0,1];
            else:
                corr[i] = np.nan; #if timeshift is totally out of range and there's no data at all - return nan
            #END IF
        #END FOR i
        k = np.where( np.abs(corr) == np.nanmax(np.abs(corr)))[0]; #get maximum values
        corr_textResults = 'Corr coeff max: '+textNice(np.round(corr[k],reportRounder))+' at time shift: '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1]+' for time intervals '+str(timeInterval/reportDivisor[0])+' '+reportDivisor[1];
        if( FLG_enableText ):
            print(corr_textResults);#print resupts
        #END IF
        if( k.size > 0 ):
            corrRet['corr'] = corr[k]; #record
            corrRet['time shift'] = timeShift[k]; #record
        else:
            corrRet['corr'] = np.array((np.nan,));
            corrRet['time shift'] = np.array((np.nan,));
        #END IF
        corrRet['time interval'] = timeInterval; #record
        corrRet['data rate'] = dataRate; #record
        corrRet['corr vect'] = corr; #record
        corrRet['time shift vect'] = timeShift; #record
        corrRet['plot name'] = plotName; #record
        corrRet['text results'] = corr_textResults; #record
        
        if( (FLG_plot == True) & (k.size > 0) ):
            #--- Declare & Unpack ---
            if( (FLG_fancyPlot >= 1) & FLG_enableText ):
                print('MAKING FANCY PLOT: corr_'+mode+plotName_fileName+'_highlite'+' IN fancyPlot FOLDER'); #report since you won't see anything
            #END IF
            # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
            PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
            PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
            # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
            # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
            # PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
            # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
            FONT_titleFM = settings_plot['font title FM'];
            # FONT_axisTick = settings_plot['font axis tick'];
            FONT_axisLabelFM = settings_plot['font axis label FM'];
            
            if( FLG_fancyPlot == 0 ): #make a loop that can do fancy and non-fancy in one (correlation calcs are expensive when done en masse)
                FLG_fancyPlot_runner = np.array( (0,) );
            elif( FLG_fancyPlot == 1 ):
                FLG_fancyPlot_runner = np.array( (0,1) );
            elif( FLG_fancyPlot == 2 ):
                FLG_fancyPlot_runner = np.array( (1,) );
            #END IF
            
            #--- Prep Plot ---
            for figr in range(0, FLG_fancyPlot_runner.size):
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
                    figManager = fig.canvas.manager; #req to maximize
                    figManager.window.showMaximized(); #force maximized
                else:
                    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                    fig, ax = plt.subplots(nrows=2, ncols=1,figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
                #END IF
                fig.subplots_adjust(hspace=0);
                # divider = make_axes_locatable(ax); #prep to add an axis
                # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
                
                #Remove the aspect ratio from the basemap so it fills the screen better
                ax[0].set_aspect('auto');
                ax[1].set_aspect('auto');
                
                #--- Actual Plotting ---
                ax[0].plot(time1/reportDivisor[0], sig1, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                for j in range(0,len(timeInterval)):
                    k_strt = np.where( np.min(np.abs( time1 - timeInterval[j,0])) == np.abs( time1 - timeInterval[j,0]) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                    k_end = np.where( np.min(np.abs( time1 - timeInterval[j,1])) == np.abs( time1 - timeInterval[j,1]) )[0][0];
                    ax[0].plot((time1[k_strt.item():k_end.item()+1])/reportDivisor[0], sig1_holder[j], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                #END FOR j
                
                ax[1].plot(time2/reportDivisor[0], sig2, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean', zorder=1);
                for j in range(0,len(timeInterval)):
                    k_strt = np.where( np.min(np.abs( (time2+timeShift[k]) - timeInterval[j,0])) == np.abs( (time2+timeShift[k]) - timeInterval[j,0]) )[0][0]; #get where time2 matches up to the time1 times adjusted by timeShift for sig2
                    k_end = np.where( np.min(np.abs( (time2+timeShift[k]) - timeInterval[j,1])) == np.abs( (time2+timeShift[k]) - timeInterval[j,1]) )[0][0];
                    if( (k_strt.size != 0) & (k_end.size != 0) & (k_strt != k_end) & (sig1_holderSize[j] == (k_end-k_strt+1)) ):
                        ax[1].plot((time2[k_strt.item():k_end.item()+1])/reportDivisor[0], sig2[k_strt.item():k_end.item()+1], linewidth=PLOT_lineWidthDoublePlus, color='xkcd:brick red', zorder=2);
                    #END IF
                #END FOR j
                
                #--- Axis and Titles and Stuff ---
                if( FLG_fancyPlot_runner[figr] == 0 ): #only title non-fancy plot
                    string_Title = mode.capitalize()+' Method, Defined Intervals Highlighted (2nd Plot has Time Shift '+textNice(timeShift[k]/reportDivisor[0])+' '+reportDivisor[1]+')';
                    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
                #END IF
                ax[1].set_xlabel('Time ['+reportDivisor[1]+']',fontproperties=FONT_axisLabelFM); #set the x axis label
                ax[0].set_ylabel(plotName[:plotName.find(' & ')],fontproperties=FONT_axisLabelFM); #set the y axis label
                ax[1].set_ylabel(plotName[plotName.find(' & ')+3:],fontproperties=FONT_axisLabelFM); #set the y axis label
                
                #nice axis ticks
                if( FLG_fancyPlot_runner[figr] == 0 ):
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],FLG_removeLabels=True,FLG_tickDirIn=True);
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True);
                else:
                    GRITI_plotHelper_axisizerTime(time1/reportDivisor[0],ax=ax[0],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=True,FLG_tickDirIn=True);
                    GRITI_plotHelper_axisizerTime(time2/reportDivisor[0],ax=ax[1],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=False,FLG_tickDirIn=True);
                #END IF
                
                figFitter(fig); #fit the fig fast
                if( FLG_fancyPlot_runner[figr] != 0 ):
                    fig.savefig(os.path.join(path_corrFancyPlot,'corr_'+mode+plotName_fileName+'_highlite'+'.png')); #save the figure
                    plt.close(); #close figure b/c it lurks apparently
                    plt.ion(); #re-enable it for later stuff
                #END IF
            #END IF
        #END FOR figr
    elif( (FLG_plot == True) & (k.size == 0) ):
        print('WARNING in subfun_correlator: No plots produced because no max found (likely due to one of the data inputs being all NaNs or times not aligning at all.');
    #END IF
    
    if( (FLG_plot == True) & (mode != 'corr') & (k.size > 0) ):
        if( (FLG_fancyPlot >= 1) & FLG_enableText ):
            print('MAKING FANCY PLOT: corr_'+mode+plotName_fileName+' IN fancyPlot FOLDER'); #report since you won't see anything
        #END IF
        
        #--- Declare & Unpack ---
        # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
        # PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
        PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
        # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
        # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
        # PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
        # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
        FONT_titleFM = settings_plot['font title FM'];
        # FONT_axisTick = settings_plot['font axis tick'];
        FONT_axisLabelFM = settings_plot['font axis label FM'];
        
        if( FLG_fancyPlot == 0 ): #make a loop that can do fancy and non-fancy in one (correlation calcs are expensive when done en masse)
            FLG_fancyPlot_runner = np.array( (0,) );
        elif( FLG_fancyPlot == 1 ):
            FLG_fancyPlot_runner = np.array( (0,1) );
        elif( FLG_fancyPlot == 2 ):
            FLG_fancyPlot_runner = np.array( (1,) );
        #END IF
        
        #--- Prep Plot ---
        for figr in range(0, FLG_fancyPlot_runner.size):
            if( FLG_fancyPlot_runner[figr] == 0 ):
                fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
                figManager = fig.canvas.manager; #req to maximize
                figManager.window.showMaximized(); #force maximized
            else:
                plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
            #END IF
            # divider = make_axes_locatable(ax); #prep to add an axis
            # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax.set_aspect('auto');
            
            #--- Actual Plotting ---
            ax.plot(timeShift/reportDivisor[0], corr, linewidth=PLOT_lineWidthPlus, color='xkcd:cerulean');
            ax.scatter(timeShift[k]/reportDivisor[0], corr[k], marker='*', s=380, color='xkcd:brick red',zorder=5);
            
            #--- Axis and Titles and Stuff ---
            if( FLG_fancyPlot_runner[figr] == 0 ): #only title non-fancy plot
                if( mode == 'shift' ):
                    string_Title = mode.capitalize()+' Method Corr Coeffs';
                elif( mode == 'range' ):
                    if( ('min' in reportDivisor[1]) & (reportDivisor[0] == 60) ):
                        string_Title = mode.capitalize()+' Method Corr Coeffs for Time Range '+textNice(timeRange[0]/3600)+' to '+textNice(timeRange[1]/3600)+' hr';
                    else:
                        string_Title = mode.capitalize()+' Method Corr Coeffs for Time Range '+textNice(timeRange[0]/reportDivisor[0])+' to '+textNice(timeRange[1]/reportDivisor[0])+' '+reportDivisor[1];
                    #END IF
                elif( mode == 'interval' ):
                    if( intervalType == 'both' ):
                        intervalType = 'Pos. & Neg';
                    #END IF
                    string_Title = mode.capitalize()+' Method Corr Coeffs for '+intervalType.capitalize()+'. Intervals Longer than '+textNice(timeInterval/reportDivisor[0])+' '+reportDivisor[1];
                #END IF
                elif( mode == 'interval' ):
                    string_Title = mode.capitalize()+' Method Corr Coeffs for Defined Intervals ';
                #END IF
                ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            #END IF
            ax.set_xlabel('Time Shift ['+reportDivisor[1]+']',fontproperties=FONT_axisLabelFM); #set the x axis label
            ax.set_ylabel('Corr. Coeff. ('+plotName+')',fontproperties=FONT_axisLabelFM); #set the y axis label
            
            #nice axis ticks
            if( FLG_fancyPlot_runner[figr] == 0 ):
                GRITI_plotHelper_axisizerTime(timeShift/reportDivisor[0],ax=ax,unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True);
            else:
                GRITI_plotHelper_axisizerTime(timeShift/reportDivisor[0],ax=ax,unit=reportDivisor[1],tickNumGoal=19,FLG_removeLabels=False,FLG_tickDirIn=True);
            #END IF
            
            figFitter(fig); #fit the fig fast
            if( FLG_fancyPlot_runner[figr] != 0 ):
                fig.savefig(os.path.join(path_corrFancyPlot,'corr_'+mode+plotName_fileName+settings_plot['save file type'])); #save the figure
                plt.close(); #close figure b/c it lurks apparently
                plt.ion(); #re-enable it for later stuff
            #END IF
        #END IF
    #END FOR figr

    return corrRet #return dict that holds everything needed
#END DEF