#Function to ray trace real fast
#RD on 5/15/2020

#Ray tracing is happening because of stuff

import numpy as np
from numba import jit #, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)
@jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def subfun_raytrace(pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,data_limd,averageSplitNum,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT):
    
    data_averaged = np.empty( (averageSplitNum,) ,dtype=np.float32 ); #preallocate, don't care what's in it
    
    kl = np.zeros( pplat_limd.size , np.bool_); #preallocate
    for j in range(0,averageSplitNum): #run through each segment
        
        for l in range(0,pplat_limd.size): #turns out it's faster to do this in a loop with numba than in an array
            #this loop runs through each point and checks to see if it is inside the 4 sided polygon area
            P1PtVx = pplong_limd[l] - temp_Long_List[j,0]; #effort to do own in polygon
            P1PtVy = pplat_limd[l] - temp_Lat_List[j,0];
            P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVx),np.abs(P1PtVy)) + sqrtApproxBo*np.minimum(np.abs(P1PtVx),np.abs(P1PtVy));
            #P3PtVxy = np.array([pplong_limd - temp_Long_List[j,2] , pplat_limd - temp_Lat_List[j,2]]); #effort to do own in polygon
            P3PtVx = pplong_limd[l] - temp_Long_List[j,2]; #effort to do own in polygon
            P3PtVy = pplat_limd[l] - temp_Lat_List[j,2];
            P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVx),np.abs(P3PtVy)) + sqrtApproxBo*np.minimum(np.abs(P3PtVx),np.abs(P3PtVy));
        
            P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVx + P1P2Vxy[1,j]*P1PtVy)/(P1P2Vm[j]*P1PtVm);
            P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVx + P1P4Vxy[1,j]*P1PtVy)/(P1P4Vm[j]*P1PtVm);
            P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVx + P3P2Vxy[1,j]*P3PtVy)/(P3P2Vm[j]*P3PtVm);
            P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVx + P3P4Vxy[1,j]*P3PtVy)/(P3P4Vm[j]*P3PtVm);
            kl[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
        #END FOR l
        
        if( np.sum(kl) != 0 ):
            data_averaged[j] = np.mean(data_limd[kl]); #average the sTEC in this lat band with the given long range
            #sTECChunked_anyAngleAvg[j] = np.sum(temp_sTEC[kl])/np.sum(kl);'
        else:
            data_averaged[j] = np.nan; #otherwise NaN it instead of 0 it
        #END IF
    #END FOR j
        
    return data_averaged; #donezo