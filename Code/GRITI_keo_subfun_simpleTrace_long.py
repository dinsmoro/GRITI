#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff

import numpy as np
from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)
# @jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
# def GRITI_keo_subfun_simpleTrace_long(time_limd,TEC_timeUnique,pplat_limd,pplong_limd,splits,latLims,vTEC_limd,avg_anyAngle_N):
    
#     vTECChunked_anyAngleAvg = np.empty( (TEC_timeUnique.size,avg_anyAngle_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
#     # latLims_max = np.float32(np.max(latLims));
#     # latLims_min = np.float32(np.min(latLims));
    
#     longLims_logical = pplat_limd <= np.float32(np.max(latLims));
#     longLims_logical[longLims_logical] = pplat_limd[longLims_logical] >= np.float32(np.min(latLims)); #this double step thing is snazzy to reduce the work on successive checks
    
#     vTEC_limd_nonans = np.logical_not(np.isnan(vTEC_limd)); #get nan locations once so don't work on NaNs at all (saves work for NaN data)
    
#     time_limd_deux = time_limd[vTEC_limd_nonans & longLims_logical]; #pre-limit
#     vTEC_limd_deux = vTEC_limd[vTEC_limd_nonans & longLims_logical]; #pre-limit
#     pplong_limd_deux = pplong_limd[vTEC_limd_nonans & longLims_logical]; #pre-limit
#     for i in prange(0, TEC_timeUnique.size ): #run through each unique time period
#         #Corral the data to the right place  
#         # k = np.where((time_limd == TEC_timeUnique[i]) & vTEC_limd_nonans & longLims_logical)[0]; #gets during a time period
#         k = np.where(time_limd_deux == TEC_timeUnique[i])[0]; #gets during a time period
        
#         temp_vTEC = vTEC_limd_deux[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
#         # temp_pplat = pplat_limd[k];
#         temp_pplong = pplong_limd_deux[k];
    
#         # kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
        
#         for j in range(0,avg_anyAngle_N): #run through each segment
            
#             # kl = (temp_pplong < splits[j+1]) &  (temp_pplong >= splits[j]) & (temp_pplat <= latLims_max) & (temp_pplat >= latLims_min); #get data in the averaging zone
#             kl = (temp_pplong < splits[j+1]) &  (temp_pplong >= splits[j]); #get data in the averaging zone
#             # for l in range(0,temp_pplat.size): #turns out it's faster to do this in a loop with numba than in an array [not for this case, samesies]
#             #     #this loop runs through each point and checks to see if it is inside the given rectangle
#             #     kl[l] = (temp_pplong[l] < splits[j+1]) &  (temp_pplong[l] >= splits[j]) & (temp_pplat[l] <= np.max(latLims)) & (temp_pplat[l] >= np.min(latLims));
#             # #END FOR l
            
#             # if( np.sum(kl) != 0 ):
#             if( np.any(kl) ):
#                 vTECChunked_anyAngleAvg[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
#                 #vTECChunked_anyAngleAvg[j] = np.sum(temp_vTEC[kl])/np.sum(kl);'
#             else:
#                 vTECChunked_anyAngleAvg[i,j] = np.nan; #otherwise NaN it instead of 0 it
#             #END IF
#         #END FOR j
#     #END FOR i
        
#     return vTECChunked_anyAngleAvg; #donezo
# #END DEF
    
@jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_simpleTrace_long(time_limd,time_unique,pplat_limd,pplong_limd,splits,latLims,data_limd,splits_N, FLG_memSafe = 1):
    
    data_chunked = np.empty( (time_unique.size,splits_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    # longLims_max = np.float32(np.max(latLims));
    # longLims_min = np.float32(np.min(latLims));
    k = (pplat_limd <= np.float32(np.max(latLims))) & (pplat_limd >= np.float32(np.min(latLims))) & np.logical_not(np.isnan(data_limd)); #make sure in longitude range AND avoid working on NaN data     

    if( FLG_memSafe ):
        time_limd = time_limd[k].copy(); #pre-limit
        data_limd = data_limd[k].copy(); #pre-limit
        pplong_limd = pplong_limd[k].copy(); #pre-limit
    else:
        time_limd = time_limd[k]; #pre-limit
        data_limd = data_limd[k]; #pre-limit
        pplong_limd = pplong_limd[k]; #pre-limit
    #END IF
    
    #no where, no search, I realized sort is the way
    jk = time_limd.argsort(); #sort em
    time_limd = time_limd[jk]; #apply sort
    data_limd = data_limd[jk]; #apply sort
    pplong_limd = pplong_limd[jk]; #apply sort
    
    time_limd_splits = np.searchsorted(time_limd, time_unique); #get the split locs
    time_limd_splits_next = np.append(time_limd_splits[1:], time_limd.size); #get the next split locks to go to - this enables parallelization b/c parallel is auto-canceled if you i+1 (also adds the last value to go to too)
    
    for i in prange(0, time_unique.size ): #run through each unique time period
        #Corral the data to the right place  
        # k = np.where((time_limd == time_unique[i]) & data_limd_nonans & longLims_logical)[0]; #gets during a time period
        # k = np.where(time_limd == time_unique[i])[0]; #gets during a time period
        # k = time_limd[time_limd_splits[i]:time_limd_splits[i+1]]; #gets during a time period
        
        temp_data = data_limd[time_limd_splits[i]:time_limd_splits_next[i]]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplong = pplong_limd[time_limd_splits[i]:time_limd_splits_next[i]];
        
        jk = temp_pplong.argsort(); #sort em
        temp_pplong = temp_pplong[jk]; #apply sort
        temp_data = temp_data[jk]; #apply sort
        
        temp_pplong_splits = np.searchsorted(temp_pplong, splits); #get the split locs
    
        for j in range(0,splits_N): #run through each segment
            if( (temp_pplong_splits[j+1]-temp_pplong_splits[j]) > 0 ):
                data_chunked[i,j] = np.mean(temp_data[temp_pplong_splits[j]:temp_pplong_splits[j+1]]); #average the vTEC in this lat band with the given long range
                #data_chunked[j] = np.sum(temp_data[temp_pplong_splits[j]:temp_pplong_splits[j+1]]])/np.sum(temp_pplong_splits[j]:temp_pplong_splits[j+1]]);'
            else:
                data_chunked[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF
        #END FOR j
        
    return data_chunked; #donezo
#END DEF