#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff

import numpy as np
from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)
@jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_keo_subfun_raytrace(time_limd,TEC_timeUnique,pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,vTEC_limd,avg_anyAngle_N,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT):
    
    vTECChunked_anyAngleAvg = np.empty( (TEC_timeUnique.size,avg_anyAngle_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    for i in prange(0, TEC_timeUnique.size ): #run through each unique time period
        #Corral the data to the right place  
        k = np.where(time_limd == TEC_timeUnique[i]); #gets during a time period
        
        temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplat = pplat_limd[k];
        temp_pplong = pplong_limd[k];
    
        kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
        
        for j in range(0,avg_anyAngle_N): #run through each segment
            
            for l in range(0,temp_pplat.size): #turns out it's faster to do this in a loop with numba than in an array
                #this loop runs through each point and checks to see if it is inside the 4 sided polygon area
                P1PtVx = temp_pplong[l] - temp_Long_List[j,0]; #effort to do own in polygon
                P1PtVy = temp_pplat[l] - temp_Lat_List[j,0];
                P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVx),np.abs(P1PtVy)) + sqrtApproxBo*np.minimum(np.abs(P1PtVx),np.abs(P1PtVy));
                # P1PtVm = np.sqrt( (P1PtVx)**2 + (P1PtVy)**2 ); #no approx for magnitude used
                #P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat - temp_Lat_List[j,2]]); #effort to do own in polygon
                P3PtVx = temp_pplong[l] - temp_Long_List[j,2]; #effort to do own in polygon
                P3PtVy = temp_pplat[l] - temp_Lat_List[j,2];
                P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVx),np.abs(P3PtVy)) + sqrtApproxBo*np.minimum(np.abs(P3PtVx),np.abs(P3PtVy));
                # P3PtVm = np.sqrt( (P3PtVx)**2 + (P3PtVy)**2 ); #no approx for magnitude used
            
                P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVx + P1P2Vxy[1,j]*P1PtVy)/(P1P2Vm[j]*P1PtVm);
                P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVx + P1P4Vxy[1,j]*P1PtVy)/(P1P4Vm[j]*P1PtVm);
                P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVx + P3P2Vxy[1,j]*P3PtVy)/(P3P2Vm[j]*P3PtVm);
                P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVx + P3P4Vxy[1,j]*P3PtVy)/(P3P4Vm[j]*P3PtVm);
                kl[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
            #END FOR l
            
            if( np.sum(kl) != 0 ):
                vTECChunked_anyAngleAvg[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
                #vTECChunked_anyAngleAvg[j] = np.sum(temp_vTEC[kl])/np.sum(kl);'
            else:
                vTECChunked_anyAngleAvg[i,j] = np.nan; #otherwise NaN it instead of 0 it
            #END IF

        #END FOR j
    #END FOR i
        
    return vTECChunked_anyAngleAvg; #donezo
    
    
    