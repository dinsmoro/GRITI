#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff

import numpy as np
from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)
# @jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
# def GRITI_keo_subfun_simpleTrace_lat(time_limd,time_unique,pplat_limd,pplong_limd,splits,longLims,data_limd,splits_N):
    
#     data_chunked = np.empty( (time_unique.size,splits_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
#     # longLims_max = np.float32(np.max(longLims));
#     # longLims_min = np.float32(np.min(longLims));
    
#     longLims_logical = pplong_limd <= np.float32(np.max(longLims));
#     longLims_logical[longLims_logical] = pplong_limd[longLims_logical] >= np.float32(np.min(longLims));
    
#     data_limd_nonans = np.logical_not(np.isnan(data_limd)); #get nan locations once so don't work on NaNs at all (saves work for NaN data)
    
#     k = data_limd_nonans & longLims_logical; #prelimits
#     time_limd_deux = time_limd[k]; #pre-limit
#     data_limd_deux = data_limd[k]; #pre-limit
#     pplat_limd_deux = pplat_limd[k]; #pre-limit
#     for i in prange(0, time_unique.size ): #run through each unique time period
#         #Corral the data to the right place  
#         # k = np.where((time_limd == time_unique[i]) & data_limd_nonans & longLims_logical)[0]; #gets during a time period
#         k = np.where(time_limd_deux == time_unique[i])[0]; #gets during a time period
        
#         temp_data = data_limd_deux[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
#         temp_pplat = pplat_limd_deux[k];
#         # temp_pplong = pplong_limd[k];

#         # kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
#         # kr = np.ones( temp_pplat.size , np.bool_); #preallocate
        
#         for j in range(0,splits_N): #run through each segment
            
#             # kl = (temp_pplong < longLims_max) &  (temp_pplong >= longLims_min) & (temp_pplat <= splits[j+1]) & (temp_pplat >= splits[j]); #get data in the averaging zone
#             kl = (temp_pplat < splits[j+1]) & (temp_pplat >= splits[j]); #get data in the averaging zone
#             # for l in range(0,temp_pplat.size): #turns out it's faster to do this in a loop with numba than in an array [not for this case, samesies]
#             #     #this loop runs through each point and checks to see if it is inside the given rectangle
#             #     kl[l] = (temp_pplong[l] < longLims_max) &  (temp_pplong[l] >= longLims_min) & (temp_pplat[l] <= splits[j+1]) & (temp_pplat[l] >= splits[j]);
#             # #END FOR l
#             #remove data that's already been used in the averaging zone
#             # kr[kr] = kr[kr] & ~kl; #yeet [did not help]
            
#             # if( np.sum(kl) != 0 ):
#             if( np.any(kl) ):
#                 data_chunked[i,j] = np.mean(temp_data[kl]); #average the vTEC in this lat band with the given long range
#                 #data_chunked[j] = np.sum(temp_data[kl])/np.sum(kl);'
#             else:
#                 data_chunked[i,j] = np.nan; #otherwise NaN it instead of 0 it
#             #END IF
#         #END FOR j
#     #END FOR i
        
#     return data_chunked; #donezo
# #END DEF

@jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_simpleTrace_lat(time_limd,time_unique,pplat_limd,pplong_limd,splits,longLims,data_limd,splits_N, FLG_memSafe = 1):
    
    data_chunked = np.empty( (time_unique.size,splits_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    # longLims_max = np.float32(np.max(longLims));
    # longLims_min = np.float32(np.min(longLims));
    k = (pplong_limd <= np.float32(np.max(longLims))) & (pplong_limd >= np.float32(np.min(longLims))) & np.logical_not(np.isnan(data_limd)); #make sure in longitude range AND avoid working on NaN data     

    if( FLG_memSafe ):
        time_limd = time_limd[k].copy(); #pre-limit
        data_limd = data_limd[k].copy(); #pre-limit
        pplat_limd = pplat_limd[k].copy(); #pre-limit
    else:
        time_limd = time_limd[k]; #pre-limit
        data_limd = data_limd[k]; #pre-limit
        pplat_limd = pplat_limd[k]; #pre-limit
    #END IF
    
    #no where, no search, I realized sort is the way
    jk = time_limd.argsort(); #sort em
    time_limd = time_limd[jk]; #apply sort
    data_limd = data_limd[jk]; #apply sort
    pplat_limd = pplat_limd[jk]; #apply sort
    
    time_limd_splits = np.searchsorted(time_limd, time_unique); #get the split locs
    time_limd_splits_next = np.append(time_limd_splits[1:], time_limd.size); #get the next split locks to go to - this enables parallelization b/c parallel is auto-canceled if you i+1 (also adds the last value to go to too)
    
    for i in prange(0, time_unique.size ): #run through each unique time period
        #Corral the data to the right place  
        # k = np.where((time_limd == time_unique[i]) & data_limd_nonans & longLims_logical)[0]; #gets during a time period
        # k = np.where(time_limd == time_unique[i])[0]; #gets during a time period
        # k = time_limd[time_limd_splits[i]:time_limd_splits[i+1]]; #gets during a time period
        
        temp_data = data_limd[time_limd_splits[i]:time_limd_splits_next[i]]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplat = pplat_limd[time_limd_splits[i]:time_limd_splits_next[i]];
        
        jk = temp_pplat.argsort(); #sort em
        temp_pplat = temp_pplat[jk]; #apply sort
        temp_data = temp_data[jk]; #apply sort
        
        temp_pplat_splits = np.searchsorted(temp_pplat, splits); #get the split locs
    
        for j in range(0,splits_N): #run through each segment
            if( (temp_pplat_splits[j+1]-temp_pplat_splits[j]) > 0 ):
                data_chunked[i,j] = np.mean(temp_data[temp_pplat_splits[j]:temp_pplat_splits[j+1]]); #average the vTEC in this lat band with the given long range
                #data_chunked[j] = np.sum(temp_data[temp_pplat_splits[j]:temp_pplat_splits[j+1]]])/np.sum(temp_pplat_splits[j]:temp_pplat_splits[j+1]]);'
            else:
                data_chunked[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF
        #END FOR j
        
    return data_chunked; #donezo
#END DEF