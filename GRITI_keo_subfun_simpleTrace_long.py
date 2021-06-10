#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff

import numpy as np
from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)
@jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_simpleTrace_long(time_limd,TEC_timeUnique,pplat_limd,pplong_limd,splits,latLims,vTEC_limd,avg_anyAngle_N):
    
    vTECChunked_anyAngleAvg = np.empty( (TEC_timeUnique.size,avg_anyAngle_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    for i in prange(0, TEC_timeUnique.size ): #run through each unique time period
        #Corral the data to the right place  
        k = np.where(time_limd == TEC_timeUnique[i]); #gets during a time period
        
        temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplat = pplat_limd[k];
        temp_pplong = pplong_limd[k];
    
        # kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
        
        for j in range(0,avg_anyAngle_N): #run through each segment
            
            kl = (temp_pplong <= splits[j+1]) &  (temp_pplong >= splits[j]) & (temp_pplat <= np.max(latLims)) & (temp_pplat >= np.min(latLims)); #get data in the averaging zone
            # for l in range(0,temp_pplat.size): #turns out it's faster to do this in a loop with numba than in an array [not for this case, samesies]
            #     #this loop runs through each point and checks to see if it is inside the given rectangle
            #     kl[l] = (temp_pplong[l] <= splits[j+1]) &  (temp_pplong[l] >= splits[j]) & (temp_pplat[l] <= np.max(latLims)) & (temp_pplat[l] >= np.min(latLims));
            # #END FOR l
            
            if( np.sum(kl) != 0 ):
                vTECChunked_anyAngleAvg[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
                #vTECChunked_anyAngleAvg[j] = np.sum(temp_vTEC[kl])/np.sum(kl);'
            else:
                vTECChunked_anyAngleAvg[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF
        #END FOR j
    #END FOR i
        
    return vTECChunked_anyAngleAvg; #donezo
    
    
    