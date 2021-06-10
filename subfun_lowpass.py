"""
#GOAL: High-pass filter something
#RD on 8/27/2018
#
#INPUT: 
time (unit doesn't matter - except must match filter_cutoffPeriod)
data - same size as time
filter_cutoffPeriod - same unit as time, default 2 hr
filter_order - default 42
windowType - default 'hann'
axisToUse - default 1 (only for multi-dimensional data inputs)
#OUTPUT:
data_hp - high-passed data
"""

import numpy as np
from scipy import signal
#from numba import jit
#
#@jit(nopython=True,nogil=False,cache=True,fastmath=True)
def subfun_lowpass(timeData, data, filter_cutoffPeriod = 0.5, filter_order = 42, windowType = 'hann', axisToUse = 1, reduceWindow = 0):
    # ===============High Pass Filtering to Excentuate Hourly Periods=========
    #These are unused, but accentuate the periods desired
    #Original file used highpass_fir.m, it has been merged
    # lp=1; # Lower period (in hrs)
    # hp=2; # Higher period (in hrs)
    # lf=(1/hp);  # Lowpass frequency corner (1/hr)
    # hf=(1/lp);  # Highpass frequency corner (1/hr)
    
    
    bs = (1/filter_cutoffPeriod);   # highpass cuttoff frequency (not sure what to make of it)
    #I think it is related to lf above (which is 1/2 )
    #ls = (1/(1/2)); #low-pass cutoff frequency
    
    time_delt = np.median(timeData[1:] - timeData[0:-1]); #Calculates a time delta based on avg of all the deltas
    
    # ===============Highpass filtering on original Zenith SNR================
    
#    n=42; # order of the Hamming window used in the custom function (uses 43 in reality)
    # c = 3.32*pi; #Hamming constant
    # M = n/2;
    # bandwidth = c/M;
    #The above was included to investigate the bandwidth - 1/2. Related to
    #bs/lf?
    
    fp = bs; # stop band stoppin freq (or pass band passing freq, depending on how you look at it)
    
    f= 1/(time_delt); #1/hr, the sampling frequency, based off of the time delta calc'd
    
    wp = 2*fp/f; # Normalizing the frequencies (Matlab is 0 to 1)
    #wl = 2*ls/f; #norm this
    #Calculation of filter coefficients
    # [b,a]=fir1(n,wp,'high'); #This takes the order and lower cut off
    #Uses the default hamming window
    # frequency ORIG
    #W = np.hanning(n+1);
    #W = signal.get_window('hann',n+1);
    # [b,a] = fir1(n,[wp,wl],W); #option for band-pass (30 min to 2 hr)
    #[b,a] = fir1(n,wp,'high',W); #just high-pass (2 hr and lower OK)
    b = signal.firwin(filter_order+1, wp, window=windowType, pass_zero=True); #just high-pass (2 hr and lower OK)
    a = 1; #for FIR a is 1
    #Applys Hanning Window for shiggles and giggles
    #NOTE: b is *very slightly* different from MATLAB output.
    
    #equiv to MATLAB's freqz(b,a,512) for checking
    #plt.figure();
    #w, h = signal.freqz(b, worN=512);
    #plt.plot( (w/np.pi), np.log10(np.abs(h))*20, linewidth=2);
    #plt.xlabel('Normalized Freqeuency (x pi rad/sample)');
    #plt.ylabel('Magnitude (dB)');
    #plt.title('Frequency Response');
    #plt.grid(True);
    #plt.show(); #req to make plot show up
    
    if( reduceWindow == 0 ):
        padlen=3*(b.size-1); #standard size to match Matlab
    else:
        padlen=(data.size-1); #make it work
    #END IF
    
    if(data.ndim == 1):
        data_hp = signal.filtfilt(b,a,data,padtype = 'odd', padlen=padlen); #Appies the filter
        #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    else:
        data_hp = signal.filtfilt(b,a,data,axis=1,padtype = 'odd', padlen=padlen); #Appies the filter
        #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    #END IF
    
    return data_hp