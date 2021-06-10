"""
Accelerates the time matching with NUMBA
"""
import numpy as np
from numba import jit, prange

@jit(nopython=True,nogil=False,cache=True,fastmath=False) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_TEC_avgPt_timeMatchACCELERATOR(avgPt_vTEC,avgPt_vTEC_time,Zenith_time,Zenith_time_delta):
    avgPt_vTEC_timeMatch = np.zeros(Zenith_time.size); #preallocate
    # avgPt_vTEC_timeMatch_HP = np.zeros(Zenith_time.size); #preallocate (non-HP'd data HP'd again after averaging - seems this method of averaging HP'd data shifts power to low-frequencies strongly which is opposite of goal)
    for i in range(0,Zenith_time.size ): #no nan time table is different
        
        if(i == 0): #deals with the fact there is TEC data before the ISR data started
            jold = np.where( np.min(np.abs((Zenith_time[i]-Zenith_time_delta) - avgPt_vTEC_time)) == np.abs((Zenith_time[i]-Zenith_time_delta) - avgPt_vTEC_time) )[0][0]; #get the matching time 5 min prev. (since it's an integration up to the time stamp given)
            jcurrent = np.where( np.min(np.abs((Zenith_time[i]) - avgPt_vTEC_time)) == np.abs((Zenith_time[i]) - avgPt_vTEC_time) )[0][0]+1; #get the matching time
            
            avgPt_vTEC_timeMatch[i] = np.mean(avgPt_vTEC[jold:jcurrent]); # average over the ISR time period
            # avgPt_vTEC_timeMatch_HP[i] = np.mean(avgPt_vTEC_HP[jold:jcurrent]); # average over the ISR time period
        else:
            jold = jcurrent; #set the old time
            jcurrent = np.where( np.min(np.abs((Zenith_time[i]) - avgPt_vTEC_time)) == np.abs((Zenith_time[i]) - avgPt_vTEC_time) )[0][0]+1; #get the matching time
            
            if( jold != jcurrent ):
                avgPt_vTEC_timeMatch[i] = np.mean(avgPt_vTEC[jold:jcurrent]); # average over the ISR time period (increment jold by 1 so not using same data)
                # avgPt_vTEC_timeMatch_HP[i] = np.mean(avgPt_vTEC_HP[jold:jcurrent]); # average over the ISR time period (increment jold by 1 so not using same data)
            else:
                avgPt_vTEC_timeMatch[i] = np.nan; #set to NaN if it would be a mean of nothing
                # avgPt_vTEC_timeMatch_HP[i] = np.nan; #set to NaN if it would be a mean of nothing
            #END IF
        #END IF
    #END FOR i
    
    return avgPt_vTEC_timeMatch