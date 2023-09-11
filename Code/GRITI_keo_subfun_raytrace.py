#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff

import numpy as np
from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)

@jit(nopython=True,nogil=True,parallel=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_raytrace_cosApprox(time_limd,time_unique,pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,data_limd,splits_N,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT, FLG_memSafe=1):
    
    data_chunked = np.empty( (time_unique.size,splits_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    temp_Lat_List_min = np.empty( splits_N );
    temp_Lat_List_max = np.empty( splits_N );
    temp_Long_List_min = np.empty( splits_N );
    temp_Long_List_max = np.empty( splits_N );
    for j in range(0,splits_N): #run through each segment
        temp_Lat_List_min[j] = np.min(temp_Lat_List[j,:]);
        temp_Lat_List_max[j] = np.max(temp_Lat_List[j,:]);
        temp_Long_List_min[j] = np.min(temp_Long_List[j,:]);
        temp_Long_List_max[j] = np.max(temp_Long_List[j,:]);
    #END FOR j
    
    k = np.logical_not(np.isnan(data_limd)); #get nan locations once so don't work on NaNs at all (saves work for NaN data)
    
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
        # k = np.where((time_limd == time_unique[i]) & k)[0]; #gets during a time period
        # k = np.where(time_limd == time_unique[i])[0];
        # temp_data = data_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        # temp_pplat = pplat_limd[k];
        # temp_pplong = pplong_limd[k];
        
        temp_data = data_limd[time_limd_splits[i]:time_limd_splits_next[i]]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplat = pplat_limd[time_limd_splits[i]:time_limd_splits_next[i]];
        temp_pplong = pplong_limd[time_limd_splits[i]:time_limd_splits_next[i]];
    
        # kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
        
        for j in range(0,splits_N): #run through each segment
            # kj = (temp_Lat_List_max[j] >= temp_pplat) & (temp_Lat_List_min[j] <= temp_pplat) & (temp_Long_List_max[j] >= temp_pplong) & (temp_Long_List_min[j] <= temp_pplong); #get the stuff that could be in the current box only
            kj = temp_Lat_List_max[j] >= temp_pplat; #step it through to do less work
            kj[kj] = temp_Lat_List_min[j] <= temp_pplat[kj];
            kj[kj] = temp_Long_List_max[j] >= temp_pplong[kj];
            kj[kj] = temp_Long_List_min[j] <= temp_pplong[kj]; #get the stuff that could be in the current box only
            temp_data_bit = temp_data[kj]; #do less work
            temp_pplat_bit = temp_pplat[kj]; #do less work
            temp_pplong_bit = temp_pplong[kj]; #do less work
            kl_bit = np.zeros( temp_pplat_bit.size , np.bool_); #preallocate
            for l in range(0,temp_pplat_bit.size): #turns out it's faster to do this in a loop with numba than in an array
                #this loop runs through each point and checks to see if it is inside the 4 sided polygon area
                P1PtVx = temp_pplong_bit[l] - temp_Long_List[j,0]; #effort to do own in polygon
                P1PtVy = temp_pplat_bit[l] - temp_Lat_List[j,0];
                P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVx),np.abs(P1PtVy)) + sqrtApproxBo*np.minimum(np.abs(P1PtVx),np.abs(P1PtVy));
                # P1PtVm = np.sqrt( (P1PtVx)**2 + (P1PtVy)**2 ); #no approx for magnitude used
                #P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat_bit - temp_Lat_List[j,2]]); #effort to do own in polygon
                P3PtVx = temp_pplong_bit[l] - temp_Long_List[j,2]; #effort to do own in polygon
                P3PtVy = temp_pplat_bit[l] - temp_Lat_List[j,2];
                P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVx),np.abs(P3PtVy)) + sqrtApproxBo*np.minimum(np.abs(P3PtVx),np.abs(P3PtVy));
                # P3PtVm = np.sqrt( (P3PtVx)**2 + (P3PtVy)**2 ); #no approx for magnitude used
            
                P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVx + P1P2Vxy[1,j]*P1PtVy)/(P1P2Vm[j]*P1PtVm);
                P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVx + P1P4Vxy[1,j]*P1PtVy)/(P1P4Vm[j]*P1PtVm);
                P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVx + P3P2Vxy[1,j]*P3PtVy)/(P3P2Vm[j]*P3PtVm);
                P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVx + P3P4Vxy[1,j]*P3PtVy)/(P3P4Vm[j]*P3PtVm);
                kl_bit[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
                # kl_bit[l] = P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]; #piece meal it to make it less work overall
                # kl_bit[l][kl_bit[l]] = P1P4P1PtdotCosT[kl_bit[l]] >= P1P2P1P4dotCosT[j]; #each successive one is less checks
                # kl_bit[l][kl_bit[l]] = P3P2P3PtdotCosT[kl_bit[l]] >= P3P2P3P4dotCosT[j]; #it worked elsewhere, maybe here too
                # kl_bit[l][kl_bit[l]] = P3P4P3PtdotCosT[kl_bit[l]] >= P3P2P3P4dotCosT[j]; #didn't test, I see this is working on a single #
            #END FOR l
            
            # if( np.sum(kl_bit) != 0 ):
            if( np.any(kl_bit) ):
                data_chunked[i,j] = np.mean(temp_data_bit[kl_bit]); #average the vTEC in this lat band with the given long range
                #data_chunked[j] = np.sum(temp_data[kl])/np.nansum(kl);
            else:
                data_chunked[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF

        #END FOR j
    #END FOR i
        
    return data_chunked; #donezo
#END DEF

@jit(nopython=True,nogil=True,parallel=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_raytrace_cosExact(time_limd,time_unique,pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,data_limd,splits_N,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT, FLG_memSafe=1):
    
    data_chunked = np.empty( (time_unique.size,splits_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    temp_Lat_List_min = np.empty( splits_N );
    temp_Lat_List_max = np.empty( splits_N );
    temp_Long_List_min = np.empty( splits_N );
    temp_Long_List_max = np.empty( splits_N );
    for j in range(0,splits_N): #run through each segment
        temp_Lat_List_min[j] = np.min(temp_Lat_List[j,:]);
        temp_Lat_List_max[j] = np.max(temp_Lat_List[j,:]);
        temp_Long_List_min[j] = np.min(temp_Long_List[j,:]);
        temp_Long_List_max[j] = np.max(temp_Long_List[j,:]);
    #END FOR j
    
    k = np.logical_not(np.isnan(data_limd)); #get nan locations once so don't work on NaNs at all (saves work for NaN data)
    
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
        # k = np.where((time_limd == time_unique[i]) & k)[0]; #gets during a time period
        # k = (time_limd == time_unique[i]);
        # temp_data = data_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        # temp_pplat = pplat_limd[k];
        # temp_pplong = pplong_limd[k];
        
        temp_data = data_limd[time_limd_splits[i]:time_limd_splits_next[i]]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplat = pplat_limd[time_limd_splits[i]:time_limd_splits_next[i]];
        temp_pplong = pplong_limd[time_limd_splits[i]:time_limd_splits_next[i]];
    
        # kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
        
        for j in range(0,splits_N): #run through each segment
            # kj = (temp_Lat_List_max[j] >= temp_pplat) & (temp_Lat_List_min[j] <= temp_pplat) & (temp_Long_List_max[j] >= temp_pplong) & (temp_Long_List_min[j] <= temp_pplong); #get the stuff that could be in the current box only
            kj = temp_Lat_List_max[j] >= temp_pplat; #step it through to do less work
            kj[kj] = temp_Lat_List_min[j] <= temp_pplat[kj];
            kj[kj] = temp_Long_List_max[j] >= temp_pplong[kj];
            kj[kj] = temp_Long_List_min[j] <= temp_pplong[kj]; #get the stuff that could be in the current box only
            temp_data_bit = temp_data[kj]; #do less work
            temp_pplat_bit = temp_pplat[kj]; #do less work
            temp_pplong_bit = temp_pplong[kj]; #do less work
            kl_bit = np.zeros( temp_pplat_bit.size , np.bool_); #preallocate
        
            for l in range(0,temp_pplat_bit.size): #turns out it's faster to do this in a loop with numba than in an array
                #this loop runs through each point and checks to see if it is inside the 4 sided polygon area
                P1PtVx = temp_pplong_bit[l] - temp_Long_List[j,0]; #effort to do own in polygon
                P1PtVy = temp_pplat_bit[l] - temp_Lat_List[j,0];
                # P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVx),np.abs(P1PtVy)) + sqrtApproxBo*np.minimum(np.abs(P1PtVx),np.abs(P1PtVy));
                P1PtVm = np.sqrt( (P1PtVx)**2 + (P1PtVy)**2 ); #no approx for magnitude used
                #P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat_bit - temp_Lat_List[j,2]]); #effort to do own in polygon
                P3PtVx = temp_pplong_bit[l] - temp_Long_List[j,2]; #effort to do own in polygon
                P3PtVy = temp_pplat_bit[l] - temp_Lat_List[j,2];
                # P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVx),np.abs(P3PtVy)) + sqrtApproxBo*np.minimum(np.abs(P3PtVx),np.abs(P3PtVy));
                P3PtVm = np.sqrt( (P3PtVx)**2 + (P3PtVy)**2 ); #no approx for magnitude used
            
                P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVx + P1P2Vxy[1,j]*P1PtVy)/(P1P2Vm[j]*P1PtVm);
                P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVx + P1P4Vxy[1,j]*P1PtVy)/(P1P4Vm[j]*P1PtVm);
                P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVx + P3P2Vxy[1,j]*P3PtVy)/(P3P2Vm[j]*P3PtVm);
                P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVx + P3P4Vxy[1,j]*P3PtVy)/(P3P4Vm[j]*P3PtVm);
                kl_bit[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
                # kl_bit[l] = P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]; #piece meal it to make it less work overall
                # kl_bit[l][kl_bit[l]] = P1P4P1PtdotCosT[kl_bit[l]] >= P1P2P1P4dotCosT[j]; #each successive one is less checks
                # kl_bit[l][kl_bit[l]] = P3P2P3PtdotCosT[kl_bit[l]] >= P3P2P3P4dotCosT[j]; #it worked elsewhere, maybe here too
                # kl_bit[l][kl_bit[l]] = P3P4P3PtdotCosT[kl_bit[l]] >= P3P2P3P4dotCosT[j]; #didn't test, I see this is working on a single #
            #END FOR l
            
            if( np.sum(kl_bit) != 0 ):
                data_chunked[i,j] = np.nanmean(temp_data_bit[kl_bit]); #average the vTEC in this lat band with the given long range
                #data_chunked[j] = np.nansum(temp_data[kl])/np.nansum(kl);
            else:
                data_chunked[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF

        #END FOR j
    #END FOR i
        
    return data_chunked; #donezo
#END DEF
    
@jit(nopython=True,nogil=True,parallel=True,fastmath={'ninf','nsz','arcp','contract','afn','reassoc'}) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_raytrace_tan(time_limd,time_unique,pplat_limd,pplong_limd,temp_Lat_List_max,temp_Lat_List_min,temp_Long_List_max,temp_Long_List_min,data_limd,keo_N,P1Pttan_limd,P3Pttan_limd,P1Pmaxtan,P1Pmintan,P3Pmaxtan,P3Pmintan, FLG_memSafe):
    #!! this mode isn't fully baked ¯\_(ツ)_/¯ roudning issues it seems, works mostly I think !!
    keo = np.empty( (time_unique.size,keo_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
        
    k = np.logical_not(np.isnan(data_limd)); #get nan locations once so don't work on NaNs at all (saves work for NaN data)
    
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
    
    for i in prange(0,time_unique.size):
        # k = np.where((time_limd == time_unique[i]) & k)[0]; #gets during a time period
        # # k = (time_limd == time_unique[i]);
        # temp_data = data_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        # # temp_pplat = pplat_limd[k];
        # # temp_pplong = pplong_limd[k];
        # temp_P1Pttan = P1Pttan_limd[k];
        # temp_P3Pttan = P3Pttan_limd[k];
        
        temp_data = data_limd[time_limd_splits[i]:time_limd_splits_next[i]]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_P1Pttan = P1Pttan_limd[time_limd_splits[i]:time_limd_splits_next[i]];
        temp_P3Pttan = P3Pttan_limd[time_limd_splits[i]:time_limd_splits_next[i]];
        
        for j in range(0,keo_N):
            kj = (temp_Lat_List_max[j] >= pplat_limd[k]) & (temp_Lat_List_min[j] <= pplat_limd[k]) & (temp_Long_List_max[j] >= pplong_limd[k]) & (temp_Long_List_min[j] <= pplong_limd[k]); #get the stuff that could be in the current box only
            temp_data_bit = temp_data[kj]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
            # temp_P1Pttan_bit = temp_P1Pttan[kj];
            # temp_P3Pttan_bit = temp_P3Pttan[kj];
            kl_bit = (P1Pmaxtan[j] >= temp_P1Pttan[kj]) & (P1Pmintan[j] <= temp_P1Pttan[kj]) & (P3Pmaxtan[j] >= temp_P3Pttan[kj]) & (P3Pmintan[j] <= temp_P3Pttan[kj]); #keep the data within the averaging zone
            if( np.sum(kl_bit) != 0 ):
                keo[i,j] = np.nanmean(temp_data_bit[kl_bit]); #average the vTEC in this lat band with the given long range
                #keo[j] = np.nansum(temp_data_bit[kl_bit])/np.nansum(kl_bit);
            else:
                keo[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF
        #END FOR j
    #END FOR i
        
    return keo; #donezo
#END DEF