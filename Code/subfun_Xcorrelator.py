"""
Finds the optimal correlation between two signals sig1 & sig2
sig1 & sig2 must be on the same cadence.
- modes:
    Xcorr : calcs correlation coeff for entirety of sig1, sig2 - time is not needed
    shift : calcs optimal correlation coeff by shifting sig2 around in time. sig1 & sig2 must be aligned in time and the same length.
    range : calcs optimal correlation coeff for a specific time range for sig1 by shifting sig2 around in time. sig1 & sig2 must be at the same cadence and sig2 needs to have enough data to be able to shift around the time range.
        timeRange - the time range in np.array((t1,t2)) or (t1,t2) or [t1,t2] format where t1 is less than t2
"""
#this was written in a feverdream (incl. subfun_peakFinder) driven by michelle branch's enitre discography on repeat
#(best song is a shootout between fault line and I'd rather be in love, I'd know)
import numpy as np
from Code.subfun_filter import subfun_filter
# from scipy.ndimage.filters import uniform_filter1d
from scipy import signal
from Code.subfun_textNice import textNice
import matplotlib.pyplot as plt
from Code.subfun_figFitter import figFitter
from Code.subfun_spectra import subfun_spectra
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings

def subfun_Xcorrelator(sig1, sig2, sig1_filt=None, sig2_filt=None, sig1_noise=None, sig2_noise=None, time1=None, time2=None, time2shift=None, \
                       dataRate=None, timeRange=None, settings_spectra=None, FLG_interpGaps=False, reportDivisor=[1,'time-units'], reportRounder=3, \
                       FLG_plot=True, FLG_plot_onlyXcorr=False, settings_plot=None, settings_paths=None, dates=None, plotName='', FLG_fancyPlot=0):
    #--- Determine mode automatically ---
    if( np.all(time2shift == None) & np.all(timeRange == None) ):
        mode = 'Xcorr';
    elif( np.all(time2shift != None) & np.all(timeRange == None) ):
        mode = 'shift';
    elif( np.all(time2shift != None) & np.all(timeRange != None) ):
        mode = 'range';
    #END IF
    
    #--- Make sure needed info is provided ---
    if( (mode != 'Xcorr') & np.all(time2 == None)):
        print('Error in subfun_Xcorrelator: Mode is "'+mode+'" but no time2 input was provided. Returning an error.'); #report error
        return 'Error in subfun_Xcorrelator' #very helpful
    #END IF
    if( (mode == 'range') & np.all(time1 == None) ):
        print('Error in subfun_Xcorrelator: Mode is "'+mode+'" but no time1 input was provided. Returning an error.'); #report error
        return 'Error in subfun_Xcorrelator' #very helpful
    #END IF
    if(  (mode == 'Xcorr') & (sig1.size != sig2.size) ):
        print('Error in subfun_Xcorrelator: sig1 and sig2 are not the same sizes ('+textNice(sig1.size)+' and '+textNice(sig2.size)+', respectively) and mode "'+mode+'" requires that. Returning an error.'); #report error
        return 'Error in subfun_Xcorrelator' #very helpful
    #END IF
    if( (mode == 'range') & np.all(timeRange == None) ):
        print('Error in subfun_Xcorrelator: Mode is "'+mode+'" but no timeRange was provided. Returning an error.'); #report error
        return 'Error in subfun_Xcorrelator' #very helpful
    #END IF
    
    #--- Say it's go time ---
    if( plotName == '' ):
        print('-- Correlation Calcing ---');
    else:
        print('-- On: '+plotName.replace('$\mathregular{','').replace('}$','')+' ---');
    #END IF
    
    #--- Prepare return variable ---
    XcorrRet = {}; #prep
    
    #--- Estimate dataRate if not avail ---
    if( dataRate == None ):
        dataRate = np.median(np.diff(time2)); #$TIME, delta of time between readings
        if( np.isclose(np.mod(dataRate,1),0.0) ):
            dataRate = np.int64(dataRate); #if it's an integer, make it one to improve accuracy/alignment w/o floating point errors
        #END IF
    #END IF
    if( settings_spectra == None ):
        winnow = None; #set to none
        nooverlap = None; #set to none
        nfft = 512; #set to 512 (6 min default value)
    else:
        if( np.isclose(360, dataRate) | (360/dataRate < 1) ):
            nfft = settings_spectra['nfft']['6min'];
            nooverlap = settings_spectra['noverlap']; #keep
            winnow = settings_spectra['window']; #keep
        else:
            #scale as needed
            nfft = settings_spectra['nfft']['6min']*360//dataRate;
            nooverlap = settings_spectra['noverlap']*360//dataRate; #adjust
            if( settings_spectra['window type'] == 'hamm' ):
                winnow = np.hamming(settings_spectra['windowLength']*360//dataRate); #adjust
            else:
                print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                import sys
                sys.crash();
            #END IF
        #END IF
    #END IF
    
    #--- Interpolate Gaps if Needed ---
    if( FLG_interpGaps == True ):
        if( np.any(np.diff(time1) > dataRate) ):
            print('Warning in subfun_Xcorrelator: There are gaps in time1 that will be interpolated over.');
            sig1 = subfun_filter( sig1, 'interp', dataTime = time1, dataRate = dataRate);
            time1 = np.arange(time1[0],time1[-1]+dataRate,dataRate); #make full time range
        #END IF
        if( np.any(np.diff(time2) > dataRate) ):
            print('Warning in subfun_Xcorrelator: There are gaps in time2 that will be interpolated over.');
            sig2 = subfun_filter( sig2, 'interp', dataTime = time2, dataRate = dataRate);
            time2 = np.arange(time2[0],time2[-1]+dataRate,dataRate); #make full time range
        #END IF
    #END IF
    
    #--- Filter if req'd ---
    if( sig1_filt != None ):
        sig1 = subfun_filter( sig1, sig1_filt, dataTime = time1, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
        if( np.all(sig1_noise != None) ):
            if( type(sig1_noise) is list ):
                for k in range(0,len(sig1_noise)):
                    sig1_noise[k] = subfun_filter( sig1_noise[k], sig1_filt, dataTime = time1, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
                #END FOR k
            elif( sig1_noise.ndim == 2 ):
                for k in range(0,sig1_noise.shape[0]):
                    sig1_noise[k,:] = subfun_filter( sig1_noise[k,:], sig1_filt, dataTime = time1, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
                #END FOR k
            else:
                sig1_noise = subfun_filter( sig1_noise, sig1_filt, dataTime = time1, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
            #END IF
        #END IF
    #END IF
    if( sig2_filt != None ):
        sig2 = subfun_filter( sig2, sig2_filt, dataTime = time2, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
        if( np.all(sig2_noise != None) ):
            if( type(sig2_noise) is list ):
                for k in range(0,len(sig2_noise)):
                    sig2_noise[k] = subfun_filter( sig2_noise[k], sig2_filt, dataTime = time2, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
                #END FOR k
            elif( sig2_noise.ndim == 2 ):
                for k in range(0,sig2_noise.shape[0]):
                    sig2_noise[k,:] = subfun_filter( sig2_noise[k,:], sig2_filt, dataTime = time2, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
                #END FOR k
            else:
                sig2_noise = subfun_filter( sig2_noise, sig2_filt, dataTime = time2, dataRate = dataRate, settings_spectra = settings_spectra, reduceWindow = 0, FLG_reportNaNs = False);
            #END IF
        #END IF
    #END IF
    
    #--- Calculate correlation coefficients for the various modes ---
    if( mode == 'Xcorr' ):
        #CALC PWRS
        pwr_sig1 = np.sqrt(1/sig1.size*np.sum(sig1**2)); #estimate power of signal
        pwr_sig2 = np.sqrt(1/sig2.size*np.sum(sig2**2)); #estimate power of signal
        
        #XCORR TIME
        [freqs_sig1vsig2,Cxy_sig1vsig2] = signal.csd(1/pwr_sig1*sig1,1/pwr_sig2*sig2,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
        # Axy_sig1vsig2 = np.angle(Cxy_sig1vsig2)*180/np.pi; 
        Pxy_sig1vsig2 = np.abs(Cxy_sig1vsig2);
        
        #SAVE IT
        XcorrRet['Xcorr'] = Pxy_sig1vsig2; #record
        XcorrRet['freqs'] = freqs_sig1vsig2; #record
        
        #noise
        if( np.all(sig1_noise != None) ):
            if( type(sig1_noise) is list ):
                Pxy_sig1_noisevsig2 = np.empty((len(sig1_noise),nfft//2+1));
                for k in range(0,len(sig1_noise)):
                    pwr_sig1_noise = np.sqrt(1/sig1_noise[k].size*np.sum(sig1_noise[k]**2)); #estimate power of signal
                    [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_noise*sig1_noise[k],1/pwr_sig2*sig2,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                    Pxy_sig1_noisevsig2[k,:] = np.abs(Cxy_sig1_noisevsig2);
                #END FOR k
                XcorrRet['Xcorr noise 1'] = np.mean(Pxy_sig1_noisevsig2,axis=0); #record
                #for plotting
                sig1_adj_noise = sig1_noise[-1];
                pwr_sig1_adj_noise = pwr_sig1_noise;
                sig1_noise_type = '('+str(len(sig1_noise))+' Iterations)'; #note it's many iterations
            elif( sig1_noise.ndim == 2 ):
                Pxy_sig1_noisevsig2 = np.empty((sig1_noise.shape[0],nfft//2+1));
                for k in range(0,sig1_noise.shape[0]):
                    pwr_sig1_noise = np.sqrt(1/sig1_noise[k,:].size*np.sum(sig1_noise[k,:]**2)); #estimate power of signal
                    [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_noise*sig1_noise[k,:],1/pwr_sig2*sig2,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                    Pxy_sig1_noisevsig2[k,:] = np.abs(Cxy_sig1_noisevsig2);
                #END FOR k
                XcorrRet['Xcorr noise 1'] = np.mean(Pxy_sig1_noisevsig2,axis=0); #record
                #for plotting
                sig1_adj_noise = sig1_noise[-1,:];
                pwr_sig1_adj_noise = pwr_sig1_noise;
                sig1_noise_type = '('+str(sig1_noise.shape[0])+' Iterations)'; #note it's many iterations
            else:
                pwr_sig1_noise = np.sqrt(1/sig1_noise.size*np.sum(sig1_noise**2)); #estimate power of signal
                [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_noise*sig1_noise,1/pwr_sig2*sig2,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                Pxy_sig1_noisevsig2 = np.abs(Cxy_sig1_noisevsig2);
                XcorrRet['Xcorr noise 1'] = Pxy_sig1_noisevsig2; #record
                #for plotting
                sig1_adj_noise = sig1_noise;
                pwr_sig1_adj_noise = pwr_sig1_noise;
                sig1_noise_type = '(Last Iteration)'; #note it's only 1 iteration
            #END IF
        #END IF
        if( np.all(sig2_noise != None) ):
            if( type(sig2_noise) is list ):
                Pxy_sig1vsig2_noise = np.empty((len(sig2_noise),nfft//2+1));
                for k in range(0,len(sig2_noise)):
                    pwr_sig2_noise = np.sqrt(1/sig2_noise[k].size*np.sum(sig2_noise[k]**2)); #estimate power of signal
                    [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1*sig1,1/pwr_sig2_noise*sig2_noise[k],window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                    Pxy_sig1vsig2_noise[k,:] = np.abs(Cxy_sig1vsig2_noise);
                #END FOR k
                XcorrRet['Xcorr noise 2'] = np.mean(Pxy_sig1vsig2_noise,axis=0); #record
                #for plotting
                sig2_adj_noise = sig2_noise[-1];
                pwr_sig2_adj_noise = pwr_sig2_noise;
                sig2_noise_type = '('+str(len(sig2_noise))+' Iterations)'; #note it's many iterations
            elif( sig2_noise.ndim == 2 ):
                Pxy_sig1vsig2_noise = np.empty((sig2_noise.shape[0],nfft//2+1));
                for k in range(0,sig2_noise.shape[0]):
                    pwr_sig2_noise = np.sqrt(1/sig2_noise[k,:].size*np.sum(sig2_noise[k,:]**2)); #estimate power of signal
                    [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1*sig1,1/pwr_sig2_noise*sig2_noise[k,:],window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                    Pxy_sig1vsig2_noise[k,:] = np.abs(Cxy_sig1vsig2_noise);
                #END FOR k
                XcorrRet['Xcorr noise 2'] = np.mean(Pxy_sig1vsig2_noise,axis=0); #record
                #for plotting
                sig2_adj_noise = sig2_noise[-1,:];
                pwr_sig2_adj_noise = pwr_sig2_noise;
                sig2_noise_type = '('+str(sig2_noise.shape[0])+' Iterations)'; #note it's many iterations
            else:
                pwr_sig2_noise = np.sqrt(1/sig2_noise.size*np.sum(sig2_noise**2)); #estimate power of signal
                [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1*sig1,1/pwr_sig2_noise*sig2_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                Pxy_sig1vsig2_noise = np.abs(Cxy_sig1vsig2_noise);
                XcorrRet['Xcorr noise 2'] = Pxy_sig1vsig2_noise; #record
                #for plotting
                sig2_adj_noise = sig2_noise;
                pwr_sig2_adj_noise = pwr_sig2_noise;
                sig2_noise_type = '(Last Iteration)'; #note it's only 1 iteration
            #END IF
        #END IF
        
        #for plotting
        time1_adj = time1; # no shift
        time2_adj = time2; # no shift
        sig1_adj = sig1;
        pwr_sig1_adj = pwr_sig1;
        sig2_adj = sig2;
        pwr_sig2_adj = pwr_sig2;        
    
    elif( mode == 'shift' ):
        #--- Adjust by the time offset ---
        if( time2shift != None ):
            time2_adj = time2 + time2shift; #adjust timeRange by a timeShift for sig2
        else:
            time2_adj = time2; # no shift
        #END IF
        
        #sig2 time adjustments
        time2_adj_indexes = np.array( ( np.where(np.min(np.abs( time2_adj - time2[0] )) == np.abs( time2_adj - time2[0] ) )[0][0] , \
            np.where(np.min(np.abs( time2_adj - time2[-1] )) == np.abs( time2_adj - time2[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
        sig2_adj = sig2[time2_adj_indexes[0]:time2_adj_indexes[1]+1];
        time2_adj = time2_adj[time2_adj_indexes[0]:time2_adj_indexes[1]+1];
        
        #sig1 time adjustments to match size of sig2
        if( time2shift != 0 ):
            sig1_adj = sig1[:sig2_adj.size]; #cut off the end of sig1 to match the size of sig2 (with pos/neg time shift the sig1[0] wants to match sig2_adj[0] and to keep size same so cut off end, future me - I am 95% confident in this)
            time1_adj = time1[:sig2_adj.size]; #cut off the end of sig1 to match the size of sig2 (with pos/neg time shift the sig1[0] wants to match sig2_adj[0] and to keep size same so cut off end, future me - I am 95% confident in this)
        else:
            sig1_adj = sig1; #timeShift is 0 and its regular
            time1_adj = time1; # no shift
        #END IF
        
        #CALC PWRS
        pwr_sig1_adj = np.sqrt(1/sig1_adj.size*np.sum(sig1_adj**2)); #estimate power of signal
        pwr_sig2_adj = np.sqrt(1/sig2_adj.size*np.sum(sig2_adj**2)); #estimate power of signal
        
        #MAKE WINDOW FIT
        if( settings_spectra['window'].size > sig2_adj.size ):
            if( settings_spectra['window type'] == 'hamm' ):
                winnow = np.hamming(sig2_adj.size - 1); #adjust
            else:
                print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                import sys
                sys.crash();
            #END IF
            nooverlap = winnow.size - 10; #adjust
        #END IF
        
        #XCORR TIME
        [freqs_sig1vsig2,Cxy_sig1vsig2] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
        # Axy_sig1vsig2 = np.angle(Cxy_sig1vsig2)*180/np.pi; 
        Pxy_sig1vsig2 = np.abs(Cxy_sig1vsig2);
        #SAVE IT
        XcorrRet['Xcorr'] = Pxy_sig1vsig2; #record
        XcorrRet['freqs'] = freqs_sig1vsig2; #record
        XcorrRet['time shift'] = time2shift; #record
        
        #noise
        if( np.all(sig1_noise != None) ):
            if( type(sig1_noise) is list ):
                Pxy_sig1_noisevsig2 = np.empty((len(sig1_noise),nfft//2+1));
                for k in range(0,len(sig1_noise)):
                    #sig1 time adjustments to match size of sig2
                    if( time2shift != 0 ):
                        sig1_adj_noise = sig1_noise[k][:sig2_adj.size]; #cut off the end of sig1 to match the size of sig2 (with pos/neg time shift the sig1[0] wants to match sig2_adj[0] and to keep size same so cut off end, future me - I am 95% confident in this)
                    else:
                        sig1_adj_noise = sig1_noise[k]; #timeShift is 0 and its regular
                    #END IF
                    pwr_sig1_adj_noise = np.sqrt(1/sig1_adj_noise.size*np.sum(sig1_adj_noise**2)); #estimate power of signal
                    [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_adj_noise*sig1_adj_noise,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                    Pxy_sig1_noisevsig2[k,:] = np.abs(Cxy_sig1_noisevsig2);
                #END FOR k
                XcorrRet['Xcorr noise 1'] = np.mean(Pxy_sig1_noisevsig2,axis=0); #record
                sig1_noise_type = '('+str(len(sig1_noise))+' Iterations)'; #note it's many iterations
            elif( sig1_noise.ndim == 2 ):
                Pxy_sig1_noisevsig2 = np.empty((sig1_noise.shape[0],nfft//2+1));
                for k in range(0,sig1_noise.shape[0]):
                    #sig1 time adjustments to match size of sig2
                    if( time2shift != 0 ):
                        sig1_adj_noise = sig1_noise[k,:sig2_adj.size]; #cut off the end of sig1 to match the size of sig2 (with pos/neg time shift the sig1[0] wants to match sig2_adj[0] and to keep size same so cut off end, future me - I am 95% confident in this)
                    else:
                        sig1_adj_noise = sig1_noise[k,:]; #timeShift is 0 and its regular
                    #END IF
                    pwr_sig1_adj_noise = np.sqrt(1/sig1_adj_noise.size*np.sum(sig1_adj_noise**2)); #estimate power of signal
                    [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_adj_noise*sig1_adj_noise,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                    Pxy_sig1_noisevsig2[k,:] = np.abs(Cxy_sig1_noisevsig2);
                #END FOR k
                XcorrRet['Xcorr noise 1'] = np.mean(Pxy_sig1_noisevsig2,axis=0); #record
                sig1_noise_type = '('+str(sig1_noise.shape[0])+' Iterations)'; #note it's many iterations
            else:
                #sig1 time adjustments to match size of sig2
                if( time2shift != 0 ):
                    sig1_adj_noise = sig1_noise[:sig2_adj.size]; #cut off the end of sig1 to match the size of sig2 (with pos/neg time shift the sig1[0] wants to match sig2_adj[0] and to keep size same so cut off end, future me - I am 95% confident in this)
                else:
                    sig1_adj_noise = sig1_noise; #timeShift is 0 and its regular
                #END IF
                pwr_sig1_adj_noise = np.sqrt(1/sig1_adj_noise.size*np.sum(sig1_adj_noise**2)); #estimate power of signal
                [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_adj_noise*sig1_adj_noise,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                Pxy_sig1_noisevsig2 = np.abs(Cxy_sig1_noisevsig2);
                XcorrRet['Xcorr noise 1'] = Pxy_sig1_noisevsig2; #record
                sig1_noise_type = '(Last Iteration)'; #note it's only 1 iteratio
            #END IF
        #END IF
        if( np.all(sig2_noise != None) ):
            if( type(sig2_noise) is list ):
                Pxy_sig1vsig2_noise = np.empty((len(sig2_noise),nfft//2+1));
                for k in range(0,len(sig2_noise)):
                    sig2_adj_noise = sig2_noise[k][time2_adj_indexes[0]:time2_adj_indexes[1]+1];
                    pwr_sig2_adj_noise = np.sqrt(1/sig2_adj_noise.size*np.sum(sig2_adj_noise**2)); #estimate power of signal
                    [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj_noise*sig2_adj_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                    Pxy_sig1vsig2_noise[k,:] = np.abs(Cxy_sig1vsig2_noise);
                #END FOR k
                XcorrRet['Xcorr noise 2'] = np.mean(Pxy_sig1vsig2_noise,axis=0); #record
                sig2_noise_type = '('+str(len(sig2_noise))+' Iterations)'; #note it's many iterations
            elif( sig2_noise.ndim == 2 ):
                Pxy_sig1vsig2_noise = np.empty((sig2_noise.shape[0],nfft//2+1));
                for k in range(0,sig2_noise.shape[0]):
                    sig2_adj_noise = sig2_noise[k,time2_adj_indexes[0]:time2_adj_indexes[1]+1];
                    pwr_sig2_adj_noise = np.sqrt(1/sig2_adj_noise.size*np.sum(sig2_adj_noise**2)); #estimate power of signal
                    [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj_noise*sig2_adj_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                    Pxy_sig1vsig2_noise[k,:] = np.abs(Cxy_sig1vsig2_noise);
                #END FOR k
                XcorrRet['Xcorr noise 2'] = np.mean(Pxy_sig1vsig2_noise,axis=0); #record
                sig2_noise_type = '('+str(sig2_noise.shape[0])+' Iterations)'; #note it's many iterations
            else:
                sig2_adj_noise = sig2_noise[time2_adj_indexes[0]:time2_adj_indexes[1]+1];
                pwr_sig2_adj_noise = np.sqrt(1/sig2_adj_noise.size*np.sum(sig2_adj_noise**2)); #estimate power of signal
                [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj_noise*sig2_adj_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                Pxy_sig1vsig2_noise = np.abs(Cxy_sig1vsig2_noise);
                XcorrRet['Xcorr noise 2'] = Pxy_sig1vsig2_noise; #record
                sig2_noise_type = '(Last Iteration)'; #note it's only 1 iteratio
            #END IF
        #END IF
            
    elif( mode == 'range' ):
        time1_adj_indexes = np.array( ( np.where(np.min(np.abs( time1 - timeRange[0] )) == np.abs( time1 - timeRange[0] ) )[0][0] , \
                np.where(np.min(np.abs( time1 - timeRange[-1] )) == np.abs( time1 - timeRange[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
        time1_adj = time1[time1_adj_indexes[0]:time1_adj_indexes[1]+1];
        sig1_adj = sig1[time1_adj_indexes[0]:time1_adj_indexes[1]+1];
        sig1_adj = subfun_filter( np.copy(sig1_adj), '0 mean'); #filter (or not)
        pwr_sig1_adj = np.sqrt(1/sig1_adj.size*np.sum(sig1_adj**2)); #estimate power of signal
        
        #--- Adjust by the time offset ---
        if( time2shift != None ):
            time2_adj = time2 + time2shift; #adjust timeRange by a timeShift for sig2
        else:
            time2_adj = time2; # no shift
        #END IF
        
        #sig2 time adjustments
        time2_adj_indexes = np.array( ( np.where(np.min(np.abs( time2_adj - timeRange[0] )) == np.abs( time2_adj - timeRange[0] ) )[0][0] , \
            np.where(np.min(np.abs( time2_adj - timeRange[-1] )) == np.abs( time2_adj - timeRange[-1] ) )[0][0] ) ); #get the indexes for that time cutout range
        sig2_adj = sig2[time2_adj_indexes[0]:time2_adj_indexes[1]+1];
        time2_adj = time2_adj[time2_adj_indexes[0]:time2_adj_indexes[1]+1];
        
        #CALC PWRS
        sig2_adj = subfun_filter( np.copy(sig2_adj), '0 mean'); #filter (or not)
        pwr_sig2_adj = np.sqrt(1/sig2_adj.size*np.sum(sig2_adj**2)); #estimate power of signal
        
        #MAKE WINDOW FIT
        if( winnow.size > sig2_adj.size ):
            if( settings_spectra['window type'] == 'hamm' ):
                winnow = np.hamming(sig2_adj.size - 1); #adjust
            else:
                print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                import sys
                sys.crash();
            #END IF
            nooverlap = winnow.size - 10; #adjust
        #END IF
        
        #XCORR TIME
        [freqs_sig1vsig2,Cxy_sig1vsig2] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
        # Axy_sig1vsig2 = np.angle(Cxy_sig1vsig2)*180/np.pi; 
        Pxy_sig1vsig2 = np.abs(Cxy_sig1vsig2);
        #SAVE IT
        XcorrRet['Xcorr'] = Pxy_sig1vsig2; #record
        XcorrRet['freqs'] = freqs_sig1vsig2; #record
        XcorrRet['time shift'] = time2shift; #record
        XcorrRet['time range'] = timeRange; #record
        
        #noise
        if( np.all(sig1_noise != None) ):
            if( type(sig1_noise) is list ):
                Pxy_sig1_noisevsig2 = np.empty((len(sig1_noise),nfft//2+1));
                for k in range(0,len(sig1_noise)):
                    #sig1 time adjustments to match size of sig2
                    sig1_adj_noise = sig1_noise[k][time1_adj_indexes[0]:time1_adj_indexes[1]+1];
                    pwr_sig1_adj_noise = np.sqrt(1/sig1_adj_noise.size*np.sum(sig1_adj_noise**2)); #estimate power of signal
                    [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_adj_noise*sig1_adj_noise,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                    Pxy_sig1_noisevsig2[k,:] = np.abs(Cxy_sig1_noisevsig2);
                #END FOR k
                XcorrRet['Xcorr noise 1'] = np.mean(Pxy_sig1_noisevsig2,axis=0); #record
                sig1_noise_type = '('+str(len(sig1_noise))+' Iterations)'; #note it's many iterations
            elif( sig1_noise.ndim == 2 ):
                Pxy_sig1_noisevsig2 = np.empty((sig1_noise.shape[0],nfft//2+1));
                for k in range(0,sig1_noise.shape[0]):
                    #sig1 time adjustments to match size of sig2
                    sig1_adj_noise = sig1_noise[k,time1_adj_indexes[0]:time1_adj_indexes[1]+1];
                    pwr_sig1_adj_noise = np.sqrt(1/sig1_adj_noise.size*np.sum(sig1_adj_noise**2)); #estimate power of signal
                    [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_adj_noise*sig1_adj_noise,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                    Pxy_sig1_noisevsig2[k,:] = np.abs(Cxy_sig1_noisevsig2);
                #END FOR k
                XcorrRet['Xcorr noise 1'] = np.mean(Pxy_sig1_noisevsig2,axis=0); #record
                sig1_noise_type = '('+str(sig1_noise.shape[0])+' Iterations)'; #note it's many iterations
            else:
                #sig1 time adjustments to match size of sig2
                sig1_adj_noise = sig1_noise[time1_adj_indexes[0]:time1_adj_indexes[1]+1];
                pwr_sig1_adj_noise = np.sqrt(1/sig1_adj_noise.size*np.sum(sig1_adj_noise**2)); #estimate power of signal
                [freqs_sig1_noisevsig2,Cxy_sig1_noisevsig2] = signal.csd(1/pwr_sig1_adj_noise*sig1_adj_noise,1/pwr_sig2_adj*sig2_adj,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                # Axy_sig1_noisevsig2 = np.angle(Cxy_sig1_noisevsig2)*180/np.pi; 
                Pxy_sig1_noisevsig2 = np.abs(Cxy_sig1_noisevsig2);
                XcorrRet['Xcorr noise 1'] = Pxy_sig1_noisevsig2; #record
                sig1_noise_type = '(Last Iteration)'; #note it's only 1 iteratio
            #END IF
        #END IF
        if( np.all(sig2_noise != None) ):
            if( type(sig2_noise) is list ):
                Pxy_sig1vsig2_noise = np.empty((len(sig2_noise),nfft//2+1));
                for k in range(0,len(sig2_noise)):
                    sig2_adj_noise = sig2_noise[k][time2_adj_indexes[0]:time2_adj_indexes[1]+1];
                    pwr_sig2_adj_noise = np.sqrt(1/sig2_adj_noise.size*np.sum(sig2_adj_noise**2)); #estimate power of signal
                    [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj_noise*sig2_adj_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                    Pxy_sig1vsig2_noise[k,:] = np.abs(Cxy_sig1vsig2_noise);
                #END FOR k
                XcorrRet['Xcorr noise 2'] = np.mean(Pxy_sig1vsig2_noise,axis=0); #record
                sig2_noise_type = '('+str(len(sig2_noise))+' Iterations)'; #note it's many iterations
            elif( sig2_noise.ndim == 2 ):
                Pxy_sig1vsig2_noise = np.empty((sig2_noise.shape[0],nfft//2+1));
                for k in range(0,sig2_noise.shape[0]):
                    sig2_adj_noise = sig2_noise[k,time2_adj_indexes[0]:time2_adj_indexes[1]+1];
                    pwr_sig2_adj_noise = np.sqrt(1/sig2_adj_noise.size*np.sum(sig2_adj_noise**2)); #estimate power of signal
                    [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj_noise*sig2_adj_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                    # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                    Pxy_sig1vsig2_noise[k,:] = np.abs(Cxy_sig1vsig2_noise);
                #END FOR k
                XcorrRet['Xcorr noise 2'] = np.mean(Pxy_sig1vsig2_noise,axis=0); #record
                sig2_noise_type = '('+str(sig2_noise.shape[0])+' Iterations)'; #note it's many iterations
            else:
                sig2_adj_noise = sig2_noise[time2_adj_indexes[0]:time2_adj_indexes[1]+1];
                pwr_sig2_adj_noise = np.sqrt(1/sig2_adj_noise.size*np.sum(sig2_adj_noise**2)); #estimate power of signal
                [freqs_sig1vsig2_noise,Cxy_sig1vsig2_noise] = signal.csd(1/pwr_sig1_adj*sig1_adj,1/pwr_sig2_adj_noise*sig2_adj_noise,window=winnow,noverlap=nooverlap,nfft=nfft,fs=1/dataRate);
                # Axy_sig1vsig2_noise = np.angle(Cxy_sig1vsig2_noise)*180/np.pi; 
                Pxy_sig1vsig2_noise = np.abs(Cxy_sig1vsig2_noise);
                XcorrRet['Xcorr noise 2'] = Pxy_sig1vsig2_noise; #record
                sig2_noise_type = '(Last Iteration)'; #note it's only 1 iteratio
            #END IF
        #END IF
    #END IF
    
    if( FLG_plot == True):
        #--- Declare & Unpack ---
        if( FLG_fancyPlot >= 1 ):
            journal_dpi = settings_plot['journal dpi'];
            if( plotName != '' ):
                plotName_fileName = '_'+plotName.lower().replace(' ','-');
            else:
                plotName_fileName = plotName.lower().replace(' ','-');
            #END IF
        #END IF
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
        
        if( FLG_plot_onlyXcorr == False ):
            #====== TIME SERIES PLOT ======
            #--- Prep Plot ---
            if( FLG_fancyPlot == 0 ):
                fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
                figManager = fig.canvas.manager; #req to maximize
                figManager.window.showMaximized(); #force maximized
            else:
                print('MAKING FANCY PLOT: Xcorr_'+mode+plotName_fileName+'_timeSeries IN fancyPlot FOLDER'); #report since you won't see anything
                plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
            #END IF
            # divider = make_axes_locatable(ax); #prep to add an axis
            # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax.set_aspect('auto');
            
            #--- Actual Plotting ---  
            cntr = 0; #prep incrementor
            pZ = []; #prep
            lZ = []; #prep
            
            pT, = ax.plot(time1_adj/reportDivisor[0],1/pwr_sig1_adj*sig1_adj,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            cntr += 1; #increment
            pZ.append(pT); #add onto the plot list
            lZ.append(plotName[:plotName.find(' & ')]); #add onto the legend list
            
            pT, = ax.plot(time2_adj/reportDivisor[0],1/pwr_sig2_adj*sig2_adj,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            cntr += 1; #increment
            pZ.append(pT); #add onto the plot list
            lZ.append(plotName[plotName.find(' & ')+3:]); #add onto the legend list
            
            # if( np.all(sig1_noise != None) ):
            #     pT, = ax.plot(time1_adj/reportDivisor[0],1/pwr_sig1_adj_noise*sig1_adj_noise,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            #     cntr += 1; #increment
            #     pZ.append(pT); #add onto the plot list
            #     lZ.append(plotName[:plotName.find(' & ')]+' Noise '+'(Last Iteration)'); #add onto the legend list
            # #END IF
            # if( np.all(sig2_noise != None) ):
            #     pT, = ax.plot(time2_adj/reportDivisor[0],1/pwr_sig2_adj_noise*sig2_adj_noise,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            #     cntr += 1; #increment
            #     pZ.append(pT); #add onto the plot list
            #     lZ.append(plotName[plotName.find(' & ')+3:]+' Noise '+'(Last Iteration)'); #add onto the legend list
            # #END IF
            
            ax.legend(pZ, lZ, loc='upper right', framealpha=0.5);
            
            #--- Axis and Titles and Stuff ---
            if( FLG_fancyPlot == 0 ): #only title non-fancy plot
                if( mode == 'Xcorr' ):
                    string_Title = plotName+' X-Corr Time Series';
                elif( mode == 'shift' ):
                    string_Title = plotName+' X-Corr Series with Time Shift of '+textNice(time2shift/reportDivisor[0])+' '+reportDivisor[1];
                elif( mode == 'range' ):
                    string_Title = plotName+' X-Corr Series for Time Range '+textNice(timeRange[0]/reportDivisor[0])+' to '+textNice(timeRange[1]/reportDivisor[0])+' '+reportDivisor[1]+' with Time Shift of '+textNice(time2shift/reportDivisor[0])+' '+reportDivisor[1];
                #END IF
                ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            #END IF
            ax.set_xlabel('Time [hr] - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
            ax.set_ylabel('Normalized Amplitude',fontproperties=FONT_axisLabelFM); #set the y axis label
            
            figFitter(fig); #fit the fig fast
            if( FLG_fancyPlot != 0 ):
                fig.savefig(os.path.join(settings_paths['fancyPlots'],'Xcorr_'+mode+plotName_fileName.replace(' ','+')+'_timeSeries'+'.png')); #save the figure
                plt.close(); #close figure b/c it lurks apparently
                plt.ion(); #re-enable it for later stuff
            #END IF
            
            
            #====== FFT PLOT ======
            XcorrRet['FFT'] = {}; #prep a sub-dict
            #--- Prep Plot ---
            if( FLG_fancyPlot == 0 ):
                fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
                figManager = fig.canvas.manager; #req to maximize
                figManager.window.showMaximized(); #force maximized
            else:
                print('MAKING FANCY PLOT: Xcorr_'+mode+plotName_fileName+'_FFT IN fancyPlot FOLDER'); #report since you won't see anything
                plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
                fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
            #END IF
            # divider = make_axes_locatable(ax); #prep to add an axis
            # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax.set_aspect('auto');
            
            #--- Actual Plotting ---  
            cntr = 0; #prep incrementor
            pZ = []; #prep
            lZ = []; #prep
            spectraMethod = 'fft'; #set 
            XcorrRet['FFT']['spectra method'] = spectraMethod;
            
            if( spectraMethod.lower() == 'fft' ):
                Cxx, period = subfun_spectra( sig1_adj, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
            elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                Cxx, period, _ = subfun_spectra( sig1_adj, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
            #END IF  
            pT, = ax.plot(period/reportDivisor[0],Cxx,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            cntr += 1; #increment
            pZ.append(pT); #add onto the plot list
            lZ.append(plotName[:plotName.find(' & ')]); #add onto the legend list
            XcorrRet['FFT']['pwr 1'] = Cxx; #record
            XcorrRet['FFT']['periods'] = period; #record
            
            if( spectraMethod.lower() == 'fft' ):
                Cxx, period = subfun_spectra( sig2_adj, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
            elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                Cxx, period, _ = subfun_spectra( sig2_adj, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
            #END IF  
            pT, = ax.plot(period/reportDivisor[0],Cxx,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            cntr += 1; #increment
            pZ.append(pT); #add onto the plot list
            lZ.append(plotName[plotName.find(' & ')+3:]); #add onto the legend list
            XcorrRet['FFT']['pwr 2'] = Cxx; #record
            
            if( np.all(sig1_noise != None) ):
                if( type(sig1_noise) is list ):
                    Cxx_avg = [None for k in range(0,len(sig1_noise))]; # preallocate
                    for k in range(0,len(sig1_noise)):
                        sig1_adj_noise = sig1_noise[k][time1_adj_indexes[0]:time1_adj_indexes[1]+1]; #get the noise portion
                        if( spectraMethod.lower() == 'fft' ):
                            Cxx_avg[k], period = subfun_spectra( sig1_adj_noise, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                            Cxx_avg[k], period, _ = subfun_spectra( sig1_adj_noise, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        #END IF
                    #END FOR k
                    pT, = ax.plot(period/reportDivisor[0],np.mean(np.asarray(Cxx_avg),axis=0),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
                    XcorrRet['FFT']['pwr noise 1'] = np.mean(np.asarray(Cxx_avg),axis=0); #record
                elif( sig1_noise.ndim == 2 ):
                    Cxx_avg = [None for k in range(0,len(sig1_noise))]; # preallocate
                    for k in range(0,sig1_noise.shape[0]):
                        sig1_adj_noise = sig1_noise[k,time1_adj_indexes[0]:time1_adj_indexes[1]+1]; #get the noise portion
                        if( spectraMethod.lower() == 'fft' ):
                            Cxx_avg[k], period = subfun_spectra( sig1_adj_noise, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                            Cxx_avg[k], period, _ = subfun_spectra( sig1_adj_noise, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        #END IF
                    #END FOR k
                    pT, = ax.plot(period/reportDivisor[0],np.mean(np.asarray(Cxx_avg),axis=0),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
                    XcorrRet['FFT']['pwr noise 1'] = np.mean(np.asarray(Cxx_avg),axis=0); #record
                else:
                    sig1_adj_noise = sig1_noise[time1_adj_indexes[0]:time1_adj_indexes[1]+1]; #get the noise portion
                    if( spectraMethod.lower() == 'fft' ):
                        Cxx_avg, period = subfun_spectra( sig1_adj_noise, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                    elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                        Cxx_avg, period, _ = subfun_spectra( sig1_adj_noise, time1_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                    #END IF
                    pT, = ax.plot(period/reportDivisor[0],Cxx_avg,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
                    XcorrRet['FFT']['pwr noise 1'] = Cxx_avg; #record
                #END IF
                cntr += 1; #increment
                pZ.append(pT); #add onto the plot list
                lZ.append(plotName[:plotName.find(' & ')]+' Noise '+sig1_noise_type); #add onto the legend list
            #END IF
            if( np.all(sig2_noise != None) ):
                if( type(sig2_noise) is list ):
                    Cxx_avg = [None for k in range(0,len(sig2_noise))]; # preallocate
                    for k in range(0,len(sig2_noise)):
                        sig2_adj_noise = sig2_noise[k][time2_adj_indexes[0]:time2_adj_indexes[1]+1]; #get the noise portion
                        if( spectraMethod.lower() == 'fft' ):
                            Cxx_avg[k], period = subfun_spectra( sig2_adj_noise, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                            Cxx_avg[k], period, _ = subfun_spectra( sig2_adj_noise, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        #END IF
                    #END FOR k
                    pT, = ax.plot(period/reportDivisor[0],np.mean(np.asarray(Cxx_avg),axis=0),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
                    XcorrRet['FFT']['pwr noise 2'] = np.mean(np.asarray(Cxx_avg),axis=0); #record
                elif( sig2_noise.ndim == 2 ):
                    Cxx_avg = [None for k in range(0,len(sig2_noise))]; # preallocate
                    for k in range(0,sig2_noise.shape[0]):
                        sig2_adj_noise = sig2_noise[k,time2_adj_indexes[0]:time2_adj_indexes[1]+1]; #get the noise portion
                        if( spectraMethod.lower() == 'fft' ):
                            Cxx_avg[k], period = subfun_spectra( sig2_adj_noise, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                            Cxx_avg[k], period, _ = subfun_spectra( sig2_adj_noise, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                        #END IF
                    #END FOR k
                    pT, = ax.plot(period/reportDivisor[0],np.mean(np.asarray(Cxx_avg),axis=0),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
                    XcorrRet['FFT']['pwr noise 2'] = np.mean(np.asarray(Cxx_avg),axis=0); #record
                else:
                    sig2_adj_noise = sig2_noise[time2_adj_indexes[0]:time2_adj_indexes[1]+1]; #get the noise portion
                    if( spectraMethod.lower() == 'fft' ):
                        Cxx_avg, period = subfun_spectra( sig2_adj_noise, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                    elif( spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle'):
                        Cxx_avg, period, _ = subfun_spectra( sig2_adj_noise, time2_adj, spectraMethod, settings_spectra, dataRate = dataRate, reduceWindow = 1); #get spectra
                    #END IF
                    pT, = ax.plot(period/reportDivisor[0],Cxx_avg,color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
                    XcorrRet['FFT']['pwr noise 2'] = Cxx_avg; #record
                #END IF
                cntr += 1; #increment
                pZ.append(pT); #add onto the plot list
                lZ.append(plotName[plotName.find(' & ')+3:]+' Noise '+sig2_noise_type); #add onto the legend list
            #END IF
            
            ax.legend(pZ, lZ, loc='upper right', framealpha=0.5);
            
            #--- Axis and Titles and Stuff ---
            if( FLG_fancyPlot == 0 ): #only title non-fancy plot
                if( mode == 'Xcorr' ):
                    string_Title = plotName+' X-Corr FFT';
                elif( mode == 'shift' ):
                    string_Title = plotName+' X-Corr FFT with Time Shift of '+textNice(time2shift/reportDivisor[0])+' '+reportDivisor[1];
                elif( mode == 'range' ):
                    string_Title = plotName+' X-Corr FFT for Time Range '+textNice(timeRange[0]/reportDivisor[0])+' to '+textNice(timeRange[1]/reportDivisor[0])+' '+reportDivisor[1]+' with Time Shift of '+textNice(time2shift/reportDivisor[0])+' '+reportDivisor[1];
                #END IF
                ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            #END IF
            ax.set_xlabel('Period ['+reportDivisor[1]+'] - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM);
            ax.set_ylabel('Arb. Power',fontproperties=FONT_axisLabelFM); #set the y axis label
            xAxisTicks = np.arange( 0, settings_spectra['period limit max']/reportDivisor[0]+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
            ax.set_xticks(xAxisTicks); #set x axis ticks
            ax.set_xlim( (settings_spectra['period limit min']/reportDivisor[0], settings_spectra['period limit max']/reportDivisor[0]) );
            
            figFitter(fig); #fit the fig fast
            if( FLG_fancyPlot != 0 ):
                fig.savefig(os.path.join(settings_paths['fancyPlots'],'Xcorr_'+mode+plotName_fileName.replace(' ','+')+'_FFT'+'.png')); #save the figure
                plt.close(); #close figure b/c it lurks apparently
                plt.ion(); #re-enable it for later stuff
            #END IF
        #END IF FLG_plot_onlyXcorr
        
        #====== XCORR PLOT ======
        warnings.filterwarnings("ignore", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
        #--- Prep Plot ---
        if( FLG_fancyPlot == 0 ):
            fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
            figManager = fig.canvas.manager; #req to maximize
            figManager.window.showMaximized(); #force maximized
        else:
            print('MAKING FANCY PLOT: Xcorr_'+mode+plotName_fileName+'_Xcorr IN fancyPlot FOLDER'); #report since you won't see anything
            plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
            fig, ax = plt.subplots(figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
        #END IF
        # divider = make_axes_locatable(ax); #prep to add an axis
        # cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax.set_aspect('auto');
        
        #--- Actual Plotting ---  
        cntr = 0; #prep incrementor
        pZ = []; #prep
        lZ = []; #prep
        
        pT, = ax.plot(1/freqs_sig1vsig2/reportDivisor[0],(Pxy_sig1vsig2),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
        cntr += 1; #increment
        pZ.append(pT); #add onto the plot list
        lZ.append(plotName); #add onto the legend list
        
        if( np.all(sig1_noise != None) ):
            pT, = ax.plot(1/freqs_sig1vsig2/reportDivisor[0],(XcorrRet['Xcorr noise 1']),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            cntr += 1; #increment
            pZ.append(pT); #add onto the plot list
            lZ.append(plotName[:plotName.find(' & ')]+' Noise & '+plotName[plotName.find(' & ')+3:]); #add onto the legend list
        #END IF
        if( np.all(sig2_noise != None) ):
            pT, = ax.plot(1/freqs_sig1vsig2/reportDivisor[0],(XcorrRet['Xcorr noise 2']),color=settings_plot['color'][cntr],linewidth=PLOT_lineWidthPlus, linestyle=settings_plot['line style'][cntr]);
            cntr += 1; #increment
            pZ.append(pT); #add onto the plot list
            lZ.append(plotName+' Noise'); #add onto the legend list
        #END IF    
        
        ax.set_ylabel('Arb. Power',fontproperties=FONT_axisLabelFM);
        
        if( mode == 'Xcorr' ):
            string_Title = plotName+' X-Corr CSD Series';
        elif( mode == 'shift' ):
            string_Title = plotName+' X-Corr CSD with Time Shift of '+textNice(time2shift/reportDivisor[0])+' '+reportDivisor[1];
        elif( mode == 'range' ):
            string_Title = plotName+' X-Corr CSD for Time Range '+textNice(timeRange[0]/reportDivisor[0])+' to '+textNice(timeRange[1]/reportDivisor[0])+' '+reportDivisor[1]+' with Time Shift of '+textNice(time2shift/reportDivisor[0])+' '+reportDivisor[1];
        #END IF
        
        ax.set_title(string_Title, fontproperties=FONT_titleFM);
        ax.legend(pZ, lZ, loc='upper right', framealpha=0.5);
        
        xAxisTicks = np.arange( 0, settings_spectra['period limit max']/reportDivisor[0]+10, 10); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
        ax.set_xticks(xAxisTicks); #set x axis ticks
        ax.set_xlim( (settings_spectra['period limit min']/reportDivisor[0], settings_spectra['period limit max']/reportDivisor[0]) );
        #END FOR i
        ax.set_xlabel('Period ['+reportDivisor[1]+'] - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM);
    
        figFitter(fig); #fit the fig fast
        if( FLG_fancyPlot != 0 ):
            fig.savefig(os.path.join(settings_paths['fancyPlots'],'Xcorr_'+mode+plotName_fileName.replace(' ','+')+'_Xcorr'+'.png')); #save the figure
            plt.close(); #close figure b/c it lurks apparently
            plt.ion(); #re-enable it for later stuff
        #END IF
        warnings.filterwarnings("default", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
    #END IF

    return XcorrRet #return dict that holds everything needed
#END DEF