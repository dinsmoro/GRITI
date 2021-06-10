#Filters whatever data comes in

def subfun_filter( dataFilt, dataTime, filtMethod, settings_spectra, dataRate = None, reduceWindow = 0):
    
    filtMethod = filtMethod.split('&'); #this allows multiple filtering methods to be chained together
    
    for i in range(0,len(filtMethod)):
        #----- No Filter Option (makes code easier) -----
        if( filtMethod[i].lower().replace('-','').replace(' ','') == 'none' ): #only original data, filtering for FFT
            
            pass; #nothing goes on ehre
        
        #----- High-pass option -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'highpass' ): #high-passed data for FFT
            from subfun_highpass import subfun_highpass
            
            dataFilt = subfun_highpass(dataTime, dataFilt, filter_cutoffPeriod=settings_spectra['filter cutoff period'], filter_order=settings_spectra['filter order'], windowType=settings_spectra['window type'], reduceWindow = reduceWindow); #high-pass that data
        
        #----- Low-pass option -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'lowpass' ): #low-passed data for FFT
            from subfun_lowpass import subfun_lowpass
            
            dataFilt = subfun_lowpass(dataTime, dataFilt, filter_cutoffPeriod=settings_spectra['filter cutoff period'], filter_order=settings_spectra['filter order'], windowType=settings_spectra['window type'], reduceWindow = reduceWindow); #low-pass that data
        
        #----- Sav-Gol option (effectively a high-pass) -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'savgol' ): #Sav-Gol filter on the data
            import numpy as np
            from scipy.signal import savgol_filter
    
            if( dataRate == None ):
                dataRate = np.median(np.diff(dataTime)); #estimate the data rate
            #END IF
            
            #prep the Sav-Gol filter for debiasing
            windowLen_savGol = np.int64(np.round(settings_spectra['savgol filter period']/dataRate)); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
            if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
                windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
            #END IF
            filtered = savgol_filter(dataFilt, windowLen_savGol, settings_spectra['savgol filter order'] ); #filter it up
            dataFilt = dataFilt - filtered; #remove the filtered data to get the "delta"
        
        #----- Simply 0 the mean -----
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == '0mean') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'mean') ): #only original data, filtering for FFT
            import numpy as np
            
            dataFilt = dataFilt - np.nanmean(dataFilt); #0 the mean
            
        #----- Simply take the log10 -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'log10' ): #only original data, filtering for FFT
            import numpy as np
            
            if( np.any(dataFilt <= 0) ):
                #print('WARNING: In subfun_filter log10 has values less than or equal to 0, adjusting them to 1 so the log works.');
                dataFilt[dataFilt <= 0] = 1; #adjust
            #END IF
            dataFilt = np.log10(dataFilt); #log10 of the data
            
        #----- Simply take the abs -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'abs' ): #only original data, filtering for FFT
            import numpy as np
            
            dataFilt = np.abs(dataFilt); #abs of the data
        #END IF
    #END FOR i
    
    return dataFilt
#END DEF