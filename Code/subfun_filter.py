"""
Filters whatever data comes in
Combine filters in order via & -> highpass & log10 will highpass the data then log10 it before returning the result.
Filters be named "high-pass" or "High Pass" and still properly "highpass"
none - does nothing
highpass - high-passes the data according to settings_spectra
lowpass - low-passes the data
savgol - Sav-Gol filters the data
0mean/mean - 0's the mean
log10 - takes log10 of the data
abs - takes absolute value of the data
neg - negates the data
mean+neg/0mean+neg - 0's the mean, negates, adds back in mean
nan/s - linearly interpolates over NaNs
zeronan/s - zeros NaNs
removenan/s - removes NaNs, do NOT filter after this b/c time steps will no longer be monotomically increasing (but filter won't be able to tell - only you can prevent invalid filtering!)
timematch - matches the data to the provided data rate
                                        
SPECIALS
"sav-gol denoise *winlen 15" signifies a filter-specific option of 15, multiple *options can be supplied
"sav-gol denoise|15" is shorthand for filters that only need 1 input to function
                                                                        
OPTIONS
dataRate - if not provided estimates it, not the end of the world to not provide
reduceWindow - if 0 will not reduce window automatically and error out, if 1 will redeuce window length to prevent error
FLG_reportNaNs - if True will report NaN stats for interpolated over or deleted when using the nan/removenan calls
"""
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from Code.subfun_highpass import subfun_highpass
from Code.subfun_lowpass import subfun_lowpass
from Code.subfun_timeMatch import subfun_timeMatch
from Code.subfun_strstr import strstr
from Code.subfun_strfind import strfind
from Code.subfun_textNice import textNice

def subfun_filter( dataFilt, filtMethod, dataTime = None, dataRate = None, settings_spectra = None, reduceWindow = 0, FLG_reportNaNs = False):
    
    FLG_returnTime = False; #usually not needed
    if( np.all(filtMethod == None) ):
        filtMethod = ''; #support None instead of empty string
    #END IF
    
    filtMethod = filtMethod.split('&'); #this allows multiple filtering methods to be chained together
    #python 3.10 introduced match/case which should be faster than if/elseif but it's so new I'm not upgrading since it'll break everything if I do rip
    for i in range(0,len(filtMethod)):
        filtMethod_curr = filtMethod[i].lower().replace('-','').replace(' ',''); #get the current method
        
        if( strfind(filtMethod[i],'|',opt=1) > 0 ): #check for filt info
            filtMethod_specifics = filtMethod[i][filtMethod[i].find('|')+1:].split(','); #get any filt specifics
            filtMethod[i] = filtMethod[i][:filtMethod[i].find('|')]; #remove the filt info
        else:
            filtMethod_specifics = None; #no specifics
        #END IF
        
        #----- No Filter Option (makes code easier) -----
        if( (filtMethod[i] == '') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'none') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'pass') ): #only original data, filtering for FFT
            
            pass; #nothing goes on ehre
        
        #----- High-pass option -----
        elif( 'highpass' in filtMethod_curr ): #high-passed data for FFT
        
            opts = filtMethod_curr.split('*')[1:];
            opt_spec = 'order'; #current option to check
            if( opt_spec in opts ):
                kj = [ij for ij, strang in enumerate(opts) if opt_spec in strang][0]; #silly strings
                filterOrder = int(opts[kj][len(opt_spec):]); #use the input
            else:
                filterOrder = settings_spectra['filter order']; #use the saved
            #END IF
        
            dataFilt = subfun_highpass(dataTime, dataFilt, filter_cutoffPeriod=settings_spectra['filter cutoff period'], filter_order=filterOrder, windowType=settings_spectra['window type'], reduceWindow = reduceWindow); #high-pass that data
        
        #----- Low-pass option -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'lowpass' ): #low-passed data for FFT
            
            dataFilt = subfun_lowpass(dataTime, dataFilt, filter_cutoffPeriod=settings_spectra['filter cutoff period'], filter_order=settings_spectra['filter order'], windowType=settings_spectra['window type'], reduceWindow = reduceWindow); #low-pass that data
        
        #----- Sav-Gol option (effectively a high-pass) -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'savgol' ): #Sav-Gol filter on the data
    
            if( dataRate == None ):
                dataRate = np.median(np.diff(dataTime)); #estimate the data rate
            #END IF
            if( np.isclose(np.mod(dataRate,1),0.0) ):
                dataRate = np.int64(dataRate); #if it's an integer, make it one to improve accuracy/alignment w/o floating point errors
            #END IF
            
            if( dataTime.size != dataFilt.size ): #only filters if data/time are same size (doesn't req time at all - just safety)
                print('WARNING in subfun_filter: In savgol filter data and time vectors are not the same size. Something is wrong or removenan was called before this filter method. Still filtering, may not be valid.');
            #END IF
            
            #prep the Sav-Gol filter for debiasing
            windowLen_savGol = np.int64(np.round(settings_spectra['savgol filter period']/dataRate)); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
            if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
                windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
            #END IF
            filtered = savgol_filter(dataFilt, windowLen_savGol, settings_spectra['savgol filter order'] ); #filter it up
            dataFilt = dataFilt - filtered; #remove the filtered data to get the "delta"
            
        #----- Sav-Gol smooth option (effectively a low-pass) -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'savgolsmooth' ): #Sav-Gol filter on the data
    
            if( dataRate == None ):
                dataRate = np.median(np.diff(dataTime)); #estimate the data rate
            #END IF
            if( np.isclose(np.mod(dataRate,1),0.0) ):
                dataRate = np.int64(dataRate); #if it's an integer, make it one to improve accuracy/alignment w/o floating point errors
            #END IF
            
            if( dataTime.size != dataFilt.size ): #only filters if data/time are same size (doesn't req time at all - just safety)
                print('WARNING in subfun_filter: In savgol filter data and time vectors are not the same size. Something is wrong or removenan was called before this filter method. Still filtering, may not be valid.');
            #END IF
            
            #prep the Sav-Gol filter for debiasing
            windowLen_savGol = np.int64(np.round(settings_spectra['savgol filter period']/dataRate)); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
            if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
                windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
            #END IF
            
            dataFilt = savgol_filter(dataFilt, windowLen_savGol, settings_spectra['savgol filter order'] ); #filter it up
            
        #----- Sav-Gol denoise option (combines a few pts to produce a smoothed signal so a lowpass but not wide band) -----
        elif( 'savgoldenoise' in filtMethod_curr ): #Sav-Gol filter on the data
    
            if( dataRate == None ):
                dataRate = np.median(np.diff(dataTime)); #estimate the data rate
            #END IF
            if( np.isclose(np.mod(dataRate,1),0.0) ):
                dataRate = np.int64(dataRate); #if it's an integer, make it one to improve accuracy/alignment w/o floating point errors
            #END IF
            
            if( dataTime.size != dataFilt.size ): #only filters if data/time are same size (doesn't req time at all - just safety)
                print('WARNING in subfun_filter: In savgol filter data and time vectors are not the same size. Something is wrong or removenan was called before this filter method. Still filtering, may not be valid.');
            #END IF
            
            #prep the Sav-Gol filter for debiasing
            # windowLen_savGol = np.int64(np.round(settings_spectra['savgol filter period']/dataRate)); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
            # if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
            #     windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
            # #END IF
            if( filtMethod_specifics == None ):
                windowLen_savGol = 11; #default
            else:
                windowLen_savGol = int(filtMethod_specifics[0]); #use what was provided
            #END IF
            
            opts = filtMethod_curr.split('*')[1:];
            opt_spec = 'winlen'; #current option to check
            if( opt_spec in opts ):
                kj = [ij for ij, strang in enumerate(opts) if opt_spec in strang][0]; #silly strings
                windowLen_savGol = int(opts[kj][len(opt_spec):]); #use the input
            else:
                windowLen_savGol = 11; #default
            #END IF
            
            dataFilt = savgol_filter(dataFilt, windowLen_savGol, settings_spectra['savgol filter order'] ); #filter it up
        
        #----- Simply 0 the mean -----
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == '0mean') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'mean') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'zeromean') ): #only original data, filtering for FFT
            
            dataFilt += -np.nanmean(dataFilt); #0 the mean
            
        #----- Simply subtract the median -----
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == '0median') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'median') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'zeromedian') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'submedian') ): #only original data, filtering for FFT
            
            dataFilt += -np.nanmedian(dataFilt); #0 the median
            
        #----- Simply take the log10 -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'log10' ): #only original data, filtering for FFT
            
            if( np.any(dataFilt <= 0) ):
                #print('WARNING: In subfun_filter log10 has values less than or equal to 0, adjusting them to 1 so the log works.');
                dataFilt[dataFilt <= 0] = 1; #adjust
            #END IF
            dataFilt = np.log10(dataFilt); #log10 of the data
            
        #----- Simply take the abs -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'abs' ): #only original data, filtering for FFT
            
            dataFilt = np.abs(dataFilt); #abs of the data
        
        #----- Simply negate -----
        elif( filtMethod[i].lower().replace('-','').replace(' ','') == 'neg'  ): #only original data, filtering for FFT
            
            dataFilt = -dataFilt; #negate the data
        
        #----- Remove mean, negate around 0, add back in the old mean -----
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == 'mean+neg') | (filtMethod[i].lower().replace('-','').replace(' ','') == '0mean+neg') ): #only original data, filtering for FFT
            
            tempMean = np.nanmean(dataFilt); #get the mean
            dataFilt += -tempMean; #remove the mean
            dataFilt = -dataFilt; #negate the data (should flip around 0 now so +->- and -->+)
            dataFilt += tempMean; #add back in the mean
            
        #--- Interpolate over Gaps ---
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == 'interp') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'interpolate') ):
            if( dataRate == None ):
                dataRate = np.median(np.diff(dataTime)); #estimate the data rate
            #END IF
            if( np.isclose(np.mod(dataRate,1),0.0) ):
                dataRate = np.int64(dataRate); #if it's an integer, make it one to improve accuracy/alignment w/o floating point errors
            #END IF
            if( np.sum(np.diff(dataTime) != dataRate) > 0 ): #only interp if req'd
                dataTime_full = np.arange(dataTime[0],dataTime[-1]+dataRate,dataRate); #make full time range
                jk = ~np.isin(dataTime_full,dataTime); #get missing data points
                interper = interp1d(dataTime,dataFilt,kind='linear',fill_value='extrapolate'); #make an interpolator out of known data
                dataFilt_temp = np.copy(dataFilt); #offload this var
                dataFilt = np.empty(dataTime_full.size); #preallocate
                dataFilt[~jk] = dataFilt_temp; #fill in the known data
                dataFilt[jk] = interper(dataTime_full[jk]); #make data for missing times
            
                if( FLG_reportNaNs == True ):
                    print('In subfun_filter: '+textNice(np.sum(jk))+'/'+textNice(dataFilt.size)+' ('+textNice(np.round(np.sum(jk)/dataFilt.size*100,2))+'%) of data were missing and interpolated over.');
                #END IF
            #END IF
            
        #--- Interpolate over NaNs ---
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == 'nan') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'nans') ):
            #interpolate over NaNs linearly
            jk = np.isnan(dataFilt); #find NaNs existing
            if( np.sum(jk) > 0 ):
                interper = interp1d(dataTime[~jk],dataFilt[~jk],kind='linear',fill_value='extrapolate'); #make an interpolator
                dataFilt[jk] = interper(dataTime[jk]); #make data for NaN times
                if( FLG_reportNaNs == True ):
                    print('In subfun_filter: '+textNice(np.sum(jk))+'/'+textNice(dataFilt.size)+' ('+textNice(np.round(np.sum(jk)/dataFilt.size*100,2))+'%) of data were NaNs and interpolated over.');
                #END IF
            #END IF
            
        #--- Zero NaNs ---
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == 'zeronan') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'zeronans') ):
            jk = np.isnan(dataFilt); #find NaNs existing
            # jl = np.where( np.isnan(dataFilt) != 1)[0]; #find NaNs not existing
            #purge NANs from data_data so it can be scargled w/o issue, also make a time var for it
            dataFilt[jk] = 0; #get rid of the NaNs
            # data_timeMatch_time = dataTime[jl]; #get the times that correspond
            if( FLG_reportNaNs == True ):
                print('In subfun_filter: '+textNice(np.sum(jk))+'/'+textNice(dataFilt.size)+' ('+textNice(np.round((np.sum(jk))/dataFilt.size*100,2))+'%) of time matched data were NaNs and zeroed.');
            #END IF
        
        #--- Remove NaNs ---
        elif( (filtMethod[i].lower().replace('-','').replace(' ','') == 'removenan') | (filtMethod[i].lower().replace('-','').replace(' ','') == 'removenans') ):
        #    jk = np.where( np.isnan(dataFilt) == 1)[0]; #find NaNs existing
            jl = ~np.isnan(dataFilt); #find NaNs not existing
            #purge NANs from data_data so it can be scargled w/o issue, also make a time var for it
            dataFilt = dataFilt[jl]; #get rid of the NaNs
            # data_timeMatch_time = dataTime[jl]; #get the times that correspond
            if( FLG_reportNaNs == True ):
                print('In subfun_filter: '+textNice(dataFilt.size-np.sum(jl))+'/'+textNice(dataFilt.size)+' ('+textNice(np.round((dataFilt.size-np.sum(jl))/dataFilt.size*100,2))+'%) of time matched data were NaNs and deleted.');
            #END IF
            
        #--- Change data rate ---
        elif( ('timematch' in filtMethod_curr) | ('dataratematch' in filtMethod_curr) ):
            if( filtMethod_specifics == None ):
                opts = filtMethod_curr.split('*')[1:];
                opt_spec = 'datarate'; #current option to check (#)
                if( strstr(opts, opt_spec).size > 0 ):
                    kj = [ij for ij, strang in enumerate(opts) if opt_spec in strang][0]; #silly strings
                    dataRate_goal = int(opts[kj][len(opt_spec):]); #use the input
                else:
                    print('WARNING in subfun_filter: '+filtMethod_curr+' needs a data rate to match to, but was not provided via the ".datarate 30" notation. Not doing anything.');
                #END IF
                opt_spec = 'usesum'; #current option to check (0 or 1)
                if( strstr(opts, opt_spec).size > 0 ):
                    kj = [ij for ij, strang in enumerate(opts) if opt_spec in strang][0]; #silly strings
                    useSum = int(opts[kj][len(opt_spec):]); #use the input
                else:
                    useSum = 0;
                #END IF
                opt_spec = 'removenans'; #current option to check (0 or 1)
                if( strstr(opts, opt_spec).size > 0 ):
                    kj = [ij for ij, strang in enumerate(opts) if opt_spec in strang][0]; #silly strings
                    removeNaNs = int(opts[kj][len(opt_spec):]); #use the input
                else:
                    removeNaNs = 0;
                #END IF
                opt_spec = 'reportnans'; #current option to check (0 or 1)
                if( strstr(opts, opt_spec).size > 0 ):
                    kj = [ij for ij, strang in enumerate(opts) if opt_spec in strang][0]; #silly strings
                    reportNaNs = int(opts[kj][len(opt_spec):]); #use the input
                else:
                    reportNaNs = 0;
                #END IF
            else:
                dataRate_goal = int(filtMethod_specifics[0]); #use what was provided
                useSum = 0; #set defaults
                removeNaNs = 0;
                reportNaNs = 0;
            #END IF
            timeMatch2 = np.arange(dataTime[0],dataTime[-1]+dataRate_goal,dataRate_goal); #get a new time range based on the goal data rate
            dataFilt, dataTime = subfun_timeMatch(dataFilt, dataTime, timeMatch2, timeMatch_delta=dataRate_goal, FLG_removeNaNs=removeNaNs, FLG_reportNaNs=reportNaNs, FLG_useSum=useSum); #yee
            dataRate = dataRate_goal; #update if filters later happen
            FLG_returnTime = True;
            
        #----- Catch all -----
        else:
            print('WARNING in subfun_filter: Unsupported filter method "'+filtMethod[i]+'". Ignoring and continuing on.');
        #END IF
    #END FOR i
    
    if( FLG_returnTime == False ):
        return dataFilt
    else:
        return dataFilt, dataTime
    #END IF
#END DEF