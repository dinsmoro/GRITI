"""
Time matches to some time scale (faster data_time scale to a slower timeMatch scale)
data_data on the data_time time stamps is matched to the timeMatch time stamps through averaging
"""
import numpy as np
from scipy.interpolate import interp1d
from numba import jit#, prange

def subfun_timeMatch(data_data, data_time, timeMatch, timeMatch_delta=None, FLG_removeNaNs=0, FLG_useSum=0):
    if( timeMatch_delta == None ):
        timeMatch_delta = np.median(np.diff(timeMatch)); #days, delta of time between readings
    #END IF
    
    if( FLG_useSum == 0 ):
        data_timeMatch = subfun_timeMatch_meanACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
    else:
        data_timeMatch = subfun_timeMatch_sumACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
    #END IF
    
    #--- Remove NaNs ---
    if( FLG_removeNaNs == 2 ):
        #interpolate over NaNs linearly
        jk = np.isnan(data_timeMatch); #find NaNs existing
        interper = interp1d(timeMatch[~jk],data_timeMatch[~jk],kind='linear',fill_value='extrapolate'); #make an interpolator
        data_timeMatch[jk] = interper(timeMatch[jk]); #make data for NaN times
        data_timeMatch_time = timeMatch; #set equal if since interpolated over NaNs
    elif( FLG_removeNaNs == 1 ):
    #    jk = np.where( np.isnan(data_timeMatch) == 1)[0]; #find NaNs existing
        jl = np.where( np.isnan(data_timeMatch) != 1)[0]; #find NaNs not existing
        #purge NANs from data_data so it can be scargled w/o issue, also make a time var for it
        data_timeMatch = data_timeMatch[jl]; #get rid of the NaNs
        data_timeMatch_time = timeMatch[jl]; #get the times that correspond
    else:
        data_timeMatch_time = timeMatch; #set equal if keeping NaNs
    #END IF
    
    return data_timeMatch, data_timeMatch_time
#END DEF

@jit(nopython=True,nogil=False,cache=True,fastmath=False) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_timeMatch_meanACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta):
    data_timeMatch = np.ones(timeMatch.size)*np.nan; #preallocate
    
    for i in range(0,timeMatch.size ): #no nan time table is different
        if(i == 0): #deals with the data before the timeMatch period occurs
            jInRange = ((timeMatch[i]-timeMatch_delta) < data_time) & (timeMatch[i] >= data_time)
            if( np.any(jInRange) ):
                data_timeMatch[i] = np.mean(data_data[jInRange]); # average over the timeMatch time period
            #END IF
            
            # jold = np.where( np.min(np.abs((timeMatch[i]-timeMatch_delta) - data_time)) == np.abs((timeMatch[i]-timeMatch_delta) - data_time) )[0][0]; #get the matching time timeMatch_delta prev. (since it's an integration up to the time stamp given)
            # jcurrent = np.where( np.min(np.abs((timeMatch[i]) - data_time)) == np.abs((timeMatch[i]) - data_time) )[0][0]+1; #get the matching time
            
            # data_timeMatch[i] = np.mean(data_data[jold:jcurrent]); # average over the timeMatch time period
        else:
            jInRange = (timeMatch[i-1] < data_time) & (timeMatch[i] >= data_time)
            if( np.any(jInRange) ):
                data_timeMatch[i] = np.mean(data_data[jInRange]); # average over the timeMatch time period
            #END IF
                
            # jold = jcurrent; #set the old time
            # jcurrent = np.where( np.min(np.abs((timeMatch[i]) - data_time)) == np.abs((timeMatch[i]) - data_time) )[0][0]+1; #get the matching time
            
            # if( jold != jcurrent ):
            #     data_timeMatch[i] = np.mean(data_data[jold:jcurrent]); # average over the timeMatch time period (increment jold by 1 so not using same data)
            # else:
            #     data_timeMatch[i] = np.nan; #set to NaN if it would be a mean of nothing
            #END IF
        #END IF
    #END FOR i
    
    return data_timeMatch
#END DEF

@jit(nopython=True,nogil=False,cache=True,fastmath=False) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_timeMatch_sumACCELERATOR(data_data,data_time,timeMatch,timeMatch_delta):
    data_timeMatch = np.ones(timeMatch.size)*np.nan; #preallocate
    
    for i in range(0,timeMatch.size ): #no nan time table is different
        if(i == 0): #deals with the data before the timeMatch period occurs
            jInRange = ((timeMatch[i]-timeMatch_delta) < data_time) & (timeMatch[i] >= data_time)
            if( np.any(jInRange) ):
                data_timeMatch[i] = np.sum(data_data[jInRange]); # average over the timeMatch time period
            #END IF
        
            # jold = np.where( np.min(np.abs((timeMatch[i]-timeMatch_delta) - data_time)) == np.abs((timeMatch[i]-timeMatch_delta) - data_time) )[0][0]; #get the matching time timeMatch_delta prev. (since it's an integration up to the time stamp given)
            # jcurrent = np.where( np.min(np.abs((timeMatch[i]) - data_time)) == np.abs((timeMatch[i]) - data_time) )[0][0]+1; #get the matching time
            
            # data_timeMatch[i] = np.sum(data_data[jold:jcurrent]); # sum over the timeMatch time period
        else:
            jInRange = (timeMatch[i-1] < data_time) & (timeMatch[i] >= data_time)
            if( np.any(jInRange) ):
                data_timeMatch[i] = np.sum(data_data[jInRange]); # average over the timeMatch time period
            #END IF
            
            # jold = jcurrent; #set the old time
            # jcurrent = np.where( np.min(np.abs((timeMatch[i]) - data_time)) == np.abs((timeMatch[i]) - data_time) )[0][0]+1; #get the matching time
            
            # if( jold != jcurrent ):
            #     data_timeMatch[i] = np.sum(data_data[jold:jcurrent]); # sum over the timeMatch time period (increment jold by 1 so not using same data)
            # else:
            #     data_timeMatch[i] = np.nan; #set to NaN if it would be a sum of nothing
            # #END IF
        #END IF
    #END FOR i
    
    return data_timeMatch
#END DEF