"""
Averages around a point
"""
import numpy as np
from scipy.interpolate import interp1d
from Code.GRITI_TEC_avgPtACCELERATOR import GRITI_TEC_avgPtACCELERATOR
from Code.subfun_highpass import subfun_highpass
import warnings

def GRITI_TEC_avgPt(TEC_timeUnique,TEC_lat,TEC_long,TEC_time,TEC_dTEC, \
    avgPt_coords,avgPt_pointRadius,Re,dateRange_dayNum_zeroHr, \
    dataReject,dataRejectOrig,dataRejectLimit,dataRejectLimitOrig,dataRejectMax,FLG_report=1):
    # warnings.filterwarnings("ignore", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo

    pointRadiusAngular = (avgPt_pointRadius/Re)*180/np.pi; #deg, (based on s=R*theta) angular radius around the point allowed
    k = pointRadiusAngular**2 >= ((TEC_lat - avgPt_coords[0])**2 + (TEC_long - avgPt_coords[1])**2); #get all of the data pts within the radius
    tempTime = TEC_time[k];
    tempTEC = TEC_dTEC[k];
    
    avgPt_vTEC = GRITI_TEC_avgPtACCELERATOR(TEC_timeUnique,TEC_timeUnique.size,tempTime,tempTEC,
        dataReject,dataRejectOrig,dataRejectLimit,dataRejectLimitOrig,dataRejectMax); #call a numba-accelerated function
    
#     avgPt_vTEC = np.zeros( (TEC_timeUnique.size)); #preallocate
#     for t in range(0,TEC_timeUnique.size):
        
#         j = tempTime == TEC_timeUnique[t]; #find the time that connect to the time we are at
        
#         tempTempTEC = tempTEC[j]; #get the delta-vTEC involved, do some filtering
        
#         tempMean = np.mean( tempTempTEC ); #get the mean
        
#         tempVar = np.var( tempTempTEC, ddof=1 ); #get variance, ddof=1 to reduce bias (matlab does automatically)
        
#         #loop to prevent too much data being nixed
#         while( (tempTempTEC.size - np.sum(((tempMean+dataReject*tempVar) > tempTempTEC) & ((tempMean-dataReject*tempVar) < tempTempTEC) ) )*100/tempTempTEC.size > dataRejectLimit ):
#             dataReject = dataReject*1.05; #increase the data rejection ratio to include more data
#             if( dataReject > dataRejectMax ):
#                 dataRejectLimit = 100; #if this goes on for a bit then it is halted with this manuever
#             #END IF
# #             (length(tempTempTEC) - sum(((tempMean+dataReject*tempVar) > tempTempTEC) & ((tempMean-dataReject*tempVar) < tempTempTEC) ) )*100/length(tempTempTEC)
# #             dataReject
# #             dataRejectLimit
#         #END WHILE
        
#         if( dataReject <= dataRejectMax ): #if this DID NOT occur leave the data, it is too sparse or too varied to deal with effectively
#             tempTempTEC = tempTempTEC[ ((tempMean+dataReject*tempVar) > tempTempTEC) & ((tempMean-dataReject*tempVar) < tempTempTEC) ]; #delete some extraneous data
#             #this is normal operation (data isn't too sparse)
#         #END IF
        
#         dataReject = dataRejectOrig; #reset dataReject
#         dataRejectLimit = dataRejectLimitOrig; #reset dataRejectLimit
        
#         avgPt_vTEC[t] = np.mean(tempTempTEC); #TECU, average of the delta-vTEC points in that radius at the same time into one value
#     #END FOR t
    
    # warnings.filterwarnings("default", category=RuntimeWarning); #silences warnings about NaNs used in a logical comparisons #yolo
    
    #now, make optional dataset that has the NaNs interpolated over
    jk = np.where( np.isnan(avgPt_vTEC) == 1)[0]; #find NaNs existing
    jl = np.where( np.isnan(avgPt_vTEC) != 1)[0]; #find NaNs not existing
    #interped data unused since spectrum analysis done using Lomb-Scargle approach which doens't require consistent timing intervals
#    avgPt_vTEC_interped = avgPt_vTEC; #save the O.G. data
#    avgPt_vTEC_interped[jk] = np.interp(jk,jl,avgPt_vTEC[jl]); #fill in the gaps with the spline fill in
    if( FLG_report == 1 ):
        print('\navgPt_vTEC Calcs:\nAt lat '+str(np.round(avgPt_coords[0],3))+' deg/long '+str(np.round(avgPt_coords[1],3))+' deg w/ '+str(avgPt_pointRadius)+' km radius, '+str(np.round(jk.size/avgPt_vTEC.size*100,3))+'% of data was NaN and interpolated over.'); #report % data filled in
    #END IF
    
    #purge NANs from avgPt_vTEC so it can be scargled w/o issue, also make a time var for it
    avgPt_vTEC_noInterp = avgPt_vTEC[jl]; #get rid of the NaNs [will have missing time steps - bad for fft-based things] (OK FOR THE SCARGLE)
    avgPt_vTEC_time_noInterp = TEC_timeUnique[jl]; #get the times that correspond [will have missing time steps - bad for fft-based things] (OK FOR THE SCARGLE)
    
    try:
        avgPt_vTEC_HP_noInterp = subfun_highpass( avgPt_vTEC_time_noInterp , avgPt_vTEC_noInterp); #highpass that data to boot [this is inherently wrong b/c of missing data points - if low % should not be too wrong]
    except:
        avgPt_vTEC_HP_noInterp = np.ones(avgPt_vTEC_noInterp.shape)*np.nan; #just nan it if it fails
    #END TRY
    
    avgPt_vTEC_time = TEC_timeUnique; #set the time
    interpFun = interp1d(avgPt_vTEC_time_noInterp, avgPt_vTEC_noInterp, kind='linear', fill_value='extrapolate'); #get an interpolation function
    avgPt_vTEC[jk] = interpFun(avgPt_vTEC_time[jk]); #only interpolate at the points that were NaNs
    
    avgPt_vTEC_HP = subfun_highpass( avgPt_vTEC_time , avgPt_vTEC); #highpass that data to boot
    
    return avgPt_vTEC, avgPt_vTEC_HP, avgPt_vTEC_time, avgPt_vTEC_noInterp, avgPt_vTEC_HP_noInterp, avgPt_vTEC_time_noInterp