#Gets the spectrum of whatever comes in
import numpy as np

def subfun_spectra( data, dataTime, spectraMethod, settings_spectra, dataRate = None, reduceWindow = 0, returnFreqs = 0):
    
    #----- Normalized FFT -----
    if( spectraMethod.lower() == 'fft' ):
        from scipy import signal
        import warnings
        
        if( dataRate == None ):
            dataRate = np.median(np.diff(dataTime)); #estimate the data rate
        #END IF
        scaler = 360//dataRate; #it scales to the 6min==512 NFFT length
        if( scaler == 0 ):
            scaler = 1; #keep it together
        #END IF
        
        nfft = settings_spectra['nfft']['6min']*scaler;
        nooverlap = settings_spectra['noverlap']*scaler; #keep reg
        if( settings_spectra['window type'] == 'hamm' ):
            winnow = np.hamming(settings_spectra['windowLength']*scaler); #adjust
        else:
            print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
            import sys
            sys.crash();
        #END IF
        
        if( reduceWindow == 1 ):
            if( winnow.size > data.size ):
                if( settings_spectra['window type'] == 'hamm' ):
                    winnow = np.hamming(data.size - 1); #adjust
                else:
                    print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                    import sys
                    sys.crash();
                #END IF
                nooverlap = winnow.size - 10; #adjust
                # if( np.isclose(30.,dataRate) ):
                #    nfft = settings_spectra['nfft']['30sec'];
                # elif( np.isclose(60.,dataRate) ):
                #     nfft = settings_spectra['nfft']['1min'];
                # elif( np.isclose(360.,dataRate) ):
                #     nfft = settings_spectra['nfft']['6min'];
                # else:
                #     nfft = settings_spectra['nfft']['6min']*360//dataRate;
                # #END IF
                # nfft = settings_spectra['nfft']['6min']*scaler;
            else:
                pass; #obsolted
                
                # if( np.isclose(30.,dataRate) ):
                #     nfft = settings_spectra['nfft']['30sec'];
                #     nooverlap = np.int64(nooverlap*(360/30)); #adjust
                #     if( settings_spectra['window type'] == 'hamm' ):
                #         winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/30))); #adjust
                #     else:
                #         print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                #         import sys
                #         sys.crash();
                #     #END IF
                # elif( np.isclose(60.,dataRate) ):
                #     nfft = settings_spectra['nfft']['1min'];
                #     nooverlap = np.int64(nooverlap*(360/60)); #adjust
                #     if( settings_spectra['window type'] == 'hamm' ):
                #         winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/60))); #adjust
                #     else:
                #         print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                #         import sys
                #         sys.crash();
                #     #END IF
                # # elif( np.isclose(360.,dataRate) ):
                # #     nfft = settings_spectra['nfft']['6min'];
                # #     nooverlap = np.int64(nooverlap*(360/360)); #adjust
                # #     if( settings_spectra['window type'] == 'hamm' ):
                # #         winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/360))); #adjust
                # #     else:
                # #         print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                # #         import sys
                # #         sys.crash();
                # #     #END IF
                # else:
                #     nfft = settings_spectra['nfft']['default'];
                #     #all other stuff is set up for 6 min
                # #END IF
            #END IF
            
        else:
            pass; #obsoleted
            # winnow = settings_spectra['window']; #keep reg
            # nooverlap = settings_spectra['noverlap']; #keep reg
            
            # nfft = settings_spectra['nfft']['6min']*scaler;
            # nooverlap = nooverlap*scaler; #adjust
            # if( settings_spectra['window type'] == 'hamm' ):
            #     winnow = np.hamming(settings_spectra['windowLength']*scaler); #adjust
            # else:
            #     print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
            #     import sys
            #     sys.crash();
            # #END IF
            # if( np.isclose(30.,dataRate) ):
            #     nfft = settings_spectra['nfft']['30sec'];
            #     nooverlap = np.int64(nooverlap*(360/30)); #adjust
            #     if( settings_spectra['window type'] == 'hamm' ):
            #         winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/30))); #adjust
            #     else:
            #         print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
            #         import sys
            #         sys.crash();
            #     #END IF
            # elif( np.isclose(60.,dataRate) ):
            #     nfft = settings_spectra['nfft']['1min'];
            #     nooverlap = np.int64(nooverlap*(360/60)); #adjust
            #     if( settings_spectra['window type'] == 'hamm' ):
            #         winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/60))); #adjust
            #     else:
            #         print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
            #         import sys
            #         sys.crash();
            #     #END IF
            # else:
            #     nfft = settings_spectra['nfft']['6min'];
            #     #all other stuff is set up for 6 min
            # #END IF
        #END IF
        
        Fs = 1/dataRate; #get the freq of the data rate
        
        pwr = np.sqrt(1/data.size*np.sum(data**2)); #estimate power of signal
        [freqs,Cxx] = signal.welch( 1/pwr*data, window=winnow,noverlap=nooverlap,nfft=nfft,fs=Fs); #get power from welch's method
        
        if( returnFreqs == 0 ):
            warnings.filterwarnings("ignore", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
            period = 1/freqs; #calc the periods
            warnings.filterwarnings("default", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
        else:
            period = freqs; #keep the freqs as freqs
        #END IF
            
        return Cxx, period
    
    #----- FFT without Normalization -----
    elif( spectraMethod.lower() == 'fft no norm' ):
        from scipy import signal
        import warnings
        
        if( dataRate == None ):
            dataRate = np.median(np.diff(dataTime)); #estimate the data rate
        #END IF
        scaler = 360//dataRate; #it scales to the 6min==512 NFFT length
        # if( reduceWindow == 1 ):
        #     if( settings_spectra['window'].size > data.size ):
        #         if( settings_spectra['window type'] == 'hamm' ):
        #             winnow = np.hamming(data.size - 1); #adjust
        #         else:
        #             print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
        #             import sys
        #             sys.crash();
        #         #END IF
        #         nooverlap = winnow.size - 10; #adjust
        #         if( np.isclose(30.,dataRate) ):
        #            nfft = settings_spectra['nfft']['30sec'];
        #         elif( np.isclose(60.,dataRate) ):
        #             nfft = settings_spectra['nfft']['1min'];
        #         else:
        #             nfft = settings_spectra['nfft']['6min'];
        #         #END IF
        #     else:
        #         winnow = settings_spectra['window']; #keep reg
        #         nooverlap = settings_spectra['noverlap']; #keep reg
        #         if( np.isclose(30.,dataRate) ):
        #             nfft = settings_spectra['nfft']['30sec'];
        #             nooverlap = np.int64(nooverlap*(360/30)); #adjust
        #             if( settings_spectra['window type'] == 'hamm' ):
        #                 winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/30))); #adjust
        #             else:
        #                 print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
        #                 import sys
        #                 sys.crash();
        #             #END IF
        #         elif( np.isclose(60.,dataRate) ):
        #             nfft = settings_spectra['nfft']['1min'];
        #             nooverlap = np.int64(nooverlap*(360/60)); #adjust
        #             if( settings_spectra['window type'] == 'hamm' ):
        #                 winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/30))); #adjust
        #             else:
        #                 print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
        #                 import sys
        #                 sys.crash();
        #             #END IF
        #         else:
        #             nfft = settings_spectra['nfft']['6min'];
        #             #all other stuff is set up for 6 min
        #         #END IF
        #     #END IF
        # else:
        #     winnow = settings_spectra['window']; #keep reg
        #     nooverlap = settings_spectra['noverlap']; #keep reg
        #     if( np.isclose(30.,dataRate) ):
        #         nfft = settings_spectra['nfft']['30sec'];
        #         nooverlap = np.int64(nooverlap*(360/30)); #adjust
        #         if( settings_spectra['window type'] == 'hamm' ):
        #             winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/30))); #adjust
        #         else:
        #             print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
        #             import sys
        #             sys.crash();
        #         #END IF
        #     elif( np.isclose(60.,dataRate) ):
        #         nfft = settings_spectra['nfft']['1min'];
        #         nooverlap = np.int64(nooverlap*(360/60)); #adjust
        #         if( settings_spectra['window type'] == 'hamm' ):
        #             winnow = np.hamming(np.int64(settings_spectra['windowLength']*(360/30))); #adjust
        #         else:
        #             print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
        #             import sys
        #             sys.crash();
        #         #END IF
        #     else:
        #         nfft = settings_spectra['nfft']['6min'];
        #         #all other stuff is set up for 6 min
        #     #END IF
        # #END IF
        if( reduceWindow == 1 ):
            if( settings_spectra['window'].size > data.size ):
                if( settings_spectra['window type'] == 'hamm' ):
                    winnow = np.hamming(data.size - 1); #adjust
                else:
                    print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                    import sys
                    sys.crash();
                #END IF
                nooverlap = winnow.size - 10; #adjust
                nfft = settings_spectra['nfft']['6min']*scaler;
            else:
                winnow = settings_spectra['window']; #keep reg
                nooverlap = settings_spectra['noverlap']; #keep reg
 
                nfft = settings_spectra['nfft']['6min']*scaler;
                if( settings_spectra['window type'] == 'hamm' ):
                    winnow = np.hamming(settings_spectra['windowLength']*scaler); #adjust
                else:
                    print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                    import sys
                    sys.crash();
                #END IF
            #END IF
        else:
            winnow = settings_spectra['window']; #keep reg
            nooverlap = settings_spectra['noverlap']; #keep reg
            
            nfft = settings_spectra['nfft']['6min']*scaler;
            nooverlap = nooverlap*scaler; #adjust
            if( settings_spectra['window type'] == 'hamm' ):
                winnow = np.hamming(settings_spectra['windowLength']*scaler); #adjust
            else:
                print('ERROR: Unsupported window type \''+settings_spectra['window type']+'\', add support for it here.');
                import sys
                sys.crash();
            #END IF
        #END IF
        
        Fs = 1/dataRate; #get the freq of the data rate
        
        [freqs,Cxx] = signal.welch(data, \
            window=winnow,noverlap=nooverlap,nfft=nfft,fs=Fs); #get power from welch's method
        
        if( returnFreqs == 0 ):
            warnings.filterwarnings("ignore", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
            period = 1/freqs; #calc the periods
            warnings.filterwarnings("default", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
        else:
            period = freqs; #keep the freqs as freqs
        #EEND IF
            
        return Cxx, period
    
    #----- Lomb-Scargle Periodogram -----
    elif( (spectraMethod.lower().replace('-','').replace(' ','') == 'lombscargle')  | (spectraMethod.lower().replace('-','').replace(' ','') == 'scargle') ):
        from Code.subfun_lombscargle import subfun_lombscargle 
        period, powerNormalized, gf = subfun_lombscargle(dataTime, data); #get the power from lomb-scargle periodogram
        
        if( returnFreqs == 1 ):
            period = 1/period; #calc the freqs
        #END IF
        
        return powerNormalized, period, gf
        
    #END IF
#END DEF