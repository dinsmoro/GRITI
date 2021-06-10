"""
Time matches to some time scale
"""
import numpy as np
from GRITI_TEC_avgPt_timeMatchACCELERATOR import GRITI_TEC_avgPt_timeMatchACCELERATOR
from subfun_highpass import subfun_highpass

def GRITI_TEC_avgPt_timeMatch(avgPt_vTEC,avgPt_vTEC_time,Zenith_time,dateRange_dayNum_zeroHr,filter_cutoffPeriod):

    Zenith_time_delta = np.median(np.diff(Zenith_time)); #days, delta of time between readings
    avgPt_vTEC_timeMatch = \
        GRITI_TEC_avgPt_timeMatchACCELERATOR(avgPt_vTEC,avgPt_vTEC_time,Zenith_time,Zenith_time_delta); #put the delta-vTEC point-averaged data on the same time cadence as the ISR data
    # avgPt_vTEC_timeMatch = np.zeros(Zenith_time.size); #preallocate
    # avgPt_vTEC_timeMatch_HP = np.zeros(Zenith_time.size); #preallocate
    # for i in range(0,Zenith_time.size ): #no nan time table is different
        
    #     if(i == 0): #deals with the fact there is TEC data before the ISR data started
    #         jold = np.where( np.min(np.abs((Zenith_time[i]-Zenith_time_delta) - avgPt_vTEC_time)) == np.abs((Zenith_time[i]-Zenith_time_delta) - avgPt_vTEC_time) )[0][0]; #get the matching time 5 min prev. (since it's an integration up to the time stamp given)
    #         jcurrent = np.where( np.min(np.abs((Zenith_time[i]) - avgPt_vTEC_time)) == np.abs((Zenith_time[i]) - avgPt_vTEC_time) )[0][0]+1; #get the matching time
            
    #         avgPt_vTEC_timeMatch[i] = np.mean(avgPt_vTEC[jold:jcurrent]); # average over the ISR time period
    #         avgPt_vTEC_timeMatch_HP[i] = np.mean(avgPt_vTEC_HP[jold:jcurrent]); # average over the ISR time period
    #     else:
    #         jold = jcurrent; #set the old time
    #         jcurrent = np.where( np.min(np.abs((Zenith_time[i]) - avgPt_vTEC_time)) == np.abs((Zenith_time[i]) - avgPt_vTEC_time) )[0][0]+1; #get the matching time
            
    #         avgPt_vTEC_timeMatch[i] = np.mean(avgPt_vTEC[jold:jcurrent]); # average over the ISR time period (increment jold by 1 so not using same data)
    #         avgPt_vTEC_timeMatch_HP[i] = np.mean(avgPt_vTEC_HP[jold:jcurrent]); # average over the ISR time period (increment jold by 1 so not using same data)
    #     #END IF
    # #END FOR i
    
    #now, make optional dataset that has the NaNs interpolated over
#    jk = np.where( np.isnan(avgPT_vTEC_timeMatch) == 1)[0]; #find NaNs existing
    jl = np.where( np.isnan(avgPt_vTEC_timeMatch) != 1)[0]; #find NaNs not existing
    #purge NANs from avgPt_vTEC so it can be scargled w/o issue, also make a time var for it
    avgPt_vTEC_timeMatch = avgPt_vTEC_timeMatch[jl]; #get rid of the NaNs
    avgPt_vTEC_timeMatch_time = Zenith_time[jl]; #get the times that correspond
    #This seems the better way to keep high-frequency stuff (which is the goal) versus averaging the previously high-passed stuff (peaks same spots, just larger here at normalized power, indicating power shifts to low-frequency when averaging high-passed data)
    avgPt_vTEC_timeMatch_HP = subfun_highpass( (avgPt_vTEC_timeMatch_time-dateRange_dayNum_zeroHr[1])*24 , avgPt_vTEC_timeMatch,filter_cutoffPeriod=filter_cutoffPeriod); #highpass that data to boot
    # avgPt_vTEC_timeMatch_HP = avgPt_vTEC_timeMatch_HP[jl]; #get rid of the NaNs [undeeded, already removed NaNs from orig array]
    
    return avgPt_vTEC_timeMatch, avgPt_vTEC_timeMatch_HP, avgPt_vTEC_timeMatch_time