"""
Time matches to some time scale (faster data_time scale to a slower timeMatch scale)
data_data on the data_time time stamps is matched to the timeMatch time stamps through averaging
"""
import numpy as np
from scipy.interpolate import interp1d
from numba import jit#, prange
from Code.subfun_textNice import textNice

def subfun_timeMatch(data_data, data_time, timeMatch2, timeMatch_delta=None, FLG_removeNaNs=0, FLG_reportNaNs=False, FLG_useSum=0):
    timeMatch = np.copy(timeMatch2); #fixes a horrifying thing where timeMatch would be modified outside of this function (not even sure where it is modified in here but whatever)
    if( timeMatch_delta == None ):
        timeMatch_delta = np.median(np.diff(timeMatch)); #time unit, delta of time between readings
    #END IF
    
    if( FLG_useSum == 0 ):
        if( np.any(np.isnan(data_data)) ):
            if( data_data.ndim == 2 ):
                indx2use = np.where(np.asarray(data_data.shape) != data_time.size)[0][0]; #get index to shape on
                if( indx2use != 0 ):
                    data_data = data_data.T; #flip it, head hurts gotta finish
                #END IF
                data_timeMatch = np.empty((data_data.shape[0],timeMatch.size),dtype=data_data.dtype); #preallocate
                for k in range(0,data_data.shape[0]):
                    data_timeMatch[k,:] = subfun_timeMatch_nanmeanACCELERATOR(data_data[k,:],data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
                #EMD FOR k
            else:
                data_timeMatch = subfun_timeMatch_nanmeanACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
            #END IF
        else:
            if( data_data.ndim == 2 ):
                indx2use = np.where(np.asarray(data_data.shape) != data_time.size)[0][0]; #get index to shape on
                if( indx2use != 0 ):
                    data_data = data_data.T; #flip it, head hurts gotta finish
                #END IF
                data_timeMatch = np.empty((data_data.shape[0],timeMatch.size),dtype=data_data.dtype); #preallocate
                for k in range(0,data_data.shape[0]):
                    data_timeMatch[k,:] = subfun_timeMatch_meanACCELERATOR(data_data[k,:],data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
                #EMD FOR k
            else:
                data_timeMatch = subfun_timeMatch_meanACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
            #END IF
        #END IF
    else:
        if( np.any(np.isnan(data_data)) ):
            if( data_data.ndim == 2 ):
                indx2use = np.where(np.asarray(data_data.shape) != data_time.size)[0][0]; #get index to shape on
                if( indx2use != 0 ):
                    data_data = data_data.T; #flip it, head hurts gotta finish
                #END IF
                data_timeMatch = np.empty((data_data.shape[0],timeMatch.size),dtype=data_data.dtype); #preallocate
                for k in range(0,data_data.shape[0]):
                    data_timeMatch[k,:] = subfun_timeMatch_nansumACCELERATOR(data_data[k,:],data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
                #EMD FOR k
            else:
                data_timeMatch = subfun_timeMatch_nansumACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
            #END IF
        else:
            if( data_data.ndim == 2 ):
                indx2use = np.where(np.asarray(data_data.shape) != data_time.size)[0][0]; #get index to shape on
                if( indx2use != 0 ):
                    data_data = data_data.T; #flip it, head hurts gotta finish
                #END IF
                data_timeMatch = np.empty((data_data.shape[0],timeMatch.size),dtype=data_data.dtype); #preallocate
                for k in range(0,data_data.shape[0]):
                    data_timeMatch[k,:] = subfun_timeMatch_sumACCELERATOR(data_data[k,:],data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
                #EMD FOR k
            else:
                data_timeMatch = subfun_timeMatch_sumACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
            #END IF
        #END IF
    #END IF
    
    #--- Remove NaNs ---
    if( FLG_removeNaNs == 2 ):
        if( data_data.ndim == 2 ):
            for k in range(0,data_data.shape[0]):
                #interpolate over NaNs linearly
                jk = np.isnan(data_timeMatch[k,:]); #find NaNs existing
                if( np.any(jk) ):
                    interper = interp1d(timeMatch[~jk],data_timeMatch[k,~jk],kind='linear',fill_value='extrapolate'); #make an interpolator
                    data_timeMatch[k,jk] = interper(timeMatch[jk]); #make data for NaN times
                    if( FLG_reportNaNs == True ):
                        print('Fun Time Match: '+textNice(np.sum(jk))+'/'+textNice(data_timeMatch.size)+' ('+textNice(np.round(np.sum(jk)/data_timeMatch.size*100,2))+'%) of time matched data [indx '+str(k+1)+'(+1!) of '+str(data_data.shape[0])+'] were NaNs and interpolated over.');
                    #END IF
                #END IF
            #END FOR k
        else:
            #interpolate over NaNs linearly
            jk = np.isnan(data_timeMatch); #find NaNs existing
            if( np.any(jk) ):
                interper = interp1d(timeMatch[~jk],data_timeMatch[~jk],kind='linear',fill_value='extrapolate'); #make an interpolator
                data_timeMatch[jk] = interper(timeMatch[jk]); #make data for NaN times
                if( FLG_reportNaNs == True ):
                    print('Fun Time Match: '+textNice(np.sum(jk))+'/'+textNice(data_timeMatch.size)+' ('+textNice(np.round(np.sum(jk)/data_timeMatch.size*100,2))+'%) of time matched data were NaNs and interpolated over.');
                #END IF
            #END IF
        #END IF
        data_timeMatch_time = timeMatch; #set equal if since interpolated over NaNs
    elif( FLG_removeNaNs == 1 ):
        if( data_data.ndim == 2 ):
            print('Fun Time Match: WARNING! FLG_removeNaNs == 1 is a BAD choice for data that has 2 dims!! It will remove the entire time indexes even if only 1 value out of N at that time index is NaN!');
            jl = ~np.isnan(data_timeMatch); #find NaNs not existing
            if( np.any(jl) ):
                #purge NANs from data_data so it can be scargled w/o issue, also make a time var for it
                data_timeMatch = data_timeMatch[:,jl]; #get rid of the NaNs
                data_timeMatch_time = timeMatch[jl]; #get the times that correspond
                if( FLG_reportNaNs == True ):
                    print('Fun Time Match: '+textNice(np.sum(~jl))+'/'+textNice(jl.size)+' ('+textNice(np.round((np.sum(~jl))/jl.size*100,2))+'%) of time matched data were NaNs and deleted.');
                #END IF
            #END IF
        else:
            jl = ~np.isnan(data_timeMatch); #find NaNs not existing
            if( np.any(jl) ):
                #purge NANs from data_data so it can be scargled w/o issue, also make a time var for it
                data_timeMatch = data_timeMatch[jl]; #get rid of the NaNs
                data_timeMatch_time = timeMatch[jl]; #get the times that correspond
                if( FLG_reportNaNs == True ):
                    print('Fun Time Match: '+textNice(np.sum(~jl))+'/'+textNice(jl.size)+' ('+textNice(np.round((np.sum(~jl))/jl.size*100,2))+'%) of time matched data were NaNs and deleted.');
                #END IF
            #END IF
        #END IF
    else:
        data_timeMatch_time = timeMatch; #set equal if keeping NaNs
    #END IF
    if( data_data.ndim == 2 ):
        if( indx2use != 0 ):
            data_timeMatch = data_timeMatch.T; #flip it back, head hurts gotta finish
        #END IF
    #END IF
        
    return data_timeMatch, data_timeMatch_time
#END DEF

@jit(nopython=True,nogil=True,cache=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_timeMatch_meanACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta):
    data_timeMatch = np.empty(timeMatch.size, dtype=data_data.dtype)*np.nan; #preallocate
    
    #deals with the data before the timeMatch period occurs
    jInRange = ((timeMatch[0]-timeMatch_delta) < data_time) & (timeMatch[0] >= data_time);
    if( np.any(jInRange) ):
        data_timeMatch[0] = np.mean(data_data[jInRange]); # average over the timeMatch time period
    #END IF
    
    for i in range(1,timeMatch.size ): #no nan time table is different
            jInRange = (timeMatch[i-1] < data_time) & (timeMatch[i] >= data_time);
            if( np.any(jInRange) ):
                data_timeMatch[i] = np.mean(data_data[jInRange]); # average over the timeMatch time period
            #END IF
        #END IF
    #END FOR i
    
    #fill in nan gaps if the data rate allows for it
    if( np.median(np.diff(data_time)) > np.float64(timeMatch_delta) ):
        for i in range(1,timeMatch.size ): #no nan time table is different
            if( np.isnan(data_timeMatch[i]) & (not np.isnan(data_timeMatch[i-1])) & (not np.isnan(data_timeMatch[i+1])) ):
                data_timeMatch[i] = (data_timeMatch[i-1] + data_timeMatch[i+1])/2; #avg them
            #END IF
        #END FOR i
    #END IF
    
    return data_timeMatch
#END DEF

@jit(nopython=True,nogil=True,cache=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_timeMatch_nanmeanACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta):
    data_timeMatch = np.empty(timeMatch.size, dtype=data_data.dtype)*np.nan; #preallocate
    
    #deals with the data before the timeMatch period occurs
    jInRange = ((timeMatch[0]-timeMatch_delta) < data_time) & (timeMatch[0] >= data_time);
    if( np.any(jInRange) ):
        data_timeMatch[0] = np.nanmean(data_data[jInRange]); # average over the timeMatch time period
    #END IF
    
    for i in range(1,timeMatch.size ): #no nan time table is different
        jInRange = (timeMatch[i-1] < data_time) & (timeMatch[i] >= data_time);
        if( np.any(jInRange) ):
            data_timeMatch[i] = np.nanmean(data_data[jInRange]); # average over the timeMatch time period
        #END IF
    #END FOR i
    
    #fill in nan gaps if the data rate allows for it
    if( np.median(np.diff(data_time)) > np.float64(timeMatch_delta) ):
        for i in range(1,timeMatch.size ): #no nan time table is different
            if( np.isnan(data_timeMatch[i]) & (not np.isnan(data_timeMatch[i-1])) & (not np.isnan(data_timeMatch[i+1])) ):
                data_timeMatch[i] = (data_timeMatch[i-1] + data_timeMatch[i+1])/2; #avg them
            #END IF
        #END FOR i
    #END IF
    
    return data_timeMatch
#END DEF

@jit(nopython=True,nogil=True,cache=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_timeMatch_sumACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta):
    data_timeMatch = np.empty(timeMatch.size, dtype=data_data.dtype)*np.nan; #preallocate
    
    #deals with the data before the timeMatch period occurs
    jInRange = ((timeMatch[0]-timeMatch_delta) < data_time) & (timeMatch[0] >= data_time);
    if( np.any(jInRange) ):
        data_timeMatch[0] = np.sum(data_data[jInRange]); # average over the timeMatch time period
    #END IF
    
    for i in range(1,timeMatch.size ): #no nan time table is different
        jInRange = (timeMatch[i-1] < data_time) & (timeMatch[i] >= data_time);
        if( np.any(jInRange) ):
            data_timeMatch[i] = np.sum(data_data[jInRange]); # average over the timeMatch time period
        #END IF
    #END FOR i
    
    #fill in nan gaps if the data rate allows for it
    if( np.median(np.diff(data_time)) > np.float64(timeMatch_delta) ):
        for i in range(1,timeMatch.size ): #no nan time table is different
            if( np.isnan(data_timeMatch[i]) & (not np.isnan(data_timeMatch[i-1])) & (not np.isnan(data_timeMatch[i+1])) ):
                data_timeMatch[i] = (data_timeMatch[i-1] + data_timeMatch[i+1])/2; #avg them
            #END IF
        #END FOR i
    #END IF
    
    return data_timeMatch
#END DEF

@jit(nopython=True,nogil=True,cache=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_timeMatch_nansumACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta):
    data_timeMatch = np.empty(timeMatch.size, dtype=data_data.dtype)*np.nan; #preallocate
    
    #deals with the data before the timeMatch period occurs
    jInRange = ((timeMatch[0]-timeMatch_delta) < data_time) & (timeMatch[0] >= data_time);
    if( np.any(jInRange) ):
        data_timeMatch[0] = np.nansum(data_data[jInRange]); # average over the timeMatch time period
    #END IF
    
    for i in range(1,timeMatch.size ): #no nan time table is different
        jInRange = (timeMatch[i-1] < data_time) & (timeMatch[i] >= data_time);
        if( np.any(jInRange) ):
            data_timeMatch[i] = np.nansum(data_data[jInRange]); # average over the timeMatch time period
        #END IF
    #END FOR i
    
    #fill in nan gaps if the data rate allows for it
    if( np.median(np.diff(data_time)) > np.float64(timeMatch_delta) ):
        for i in range(1,timeMatch.size ): #no nan time table is different
            if( np.isnan(data_timeMatch[i]) & (not np.isnan(data_timeMatch[i-1])) & (not np.isnan(data_timeMatch[i+1])) ):
                data_timeMatch[i] = (data_timeMatch[i-1] + data_timeMatch[i+1])/2; #avg them
            #END IF
        #END FOR i
    #END IF
    
    return data_timeMatch
#END DEF