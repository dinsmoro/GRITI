#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff

#inpoly = []; #prep a list
#for j in range( 0, avg_anyAngle_N):
#    inpoly.append([(temp_Long_List[j,0],temp_Lat_List[j,1]), (temp_Long_List[j,0], temp_Lat_List[j,2]), \
#            (temp_Long_List[j,1], temp_Lat_List[j,2]), (temp_Long_List[j,1], temp_Lat_List[j,2])]);  # square with legs length 1 and bottom left corner at the origin
##END FOR j
#    
#inpolygon_latMin = np.min(temp_Lat_List[j][1:3]);
#inpolygon_latMax = np.max(temp_Lat_List[j][1:3]);
#inpolygon_longMin = np.min(temp_Long_List[j][0:2]);
#inpolygon_longMax = np.max(temp_Long_List[j][0:2]);
#
#
#
#for i in range(0,len(temp_pplatpplong)):
#    
#    kLatLong = (temp_pplatpplong[:,1] < inpolygon_latMax) & (temp_pplatpplong[:,1] > inpolygon_latMin) & (temp_pplatpplong[:,0] < inpolygon_longMax) & (temp_pplatpplong[:,0] > inpolygon_longMin) 
#    latlongs = temp_pplatpplong[kLatLong,:]
#
#    kl2 = 

import numpy as np
from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)
@jit(nopython=True,nogil=True,parallel=False) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_TEC_avgAnyAngle_subfun_Raytrace(temp_pplat,temp_pplong,temp_Lat_List,temp_Long_List,temp_vTEC,avg_anyAngle_N,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT):
    
    vTECChunked_anyAngleAvg = np.zeros( avg_anyAngle_N , dtype=np.float32 ); #preallocate #,dtype=np.float32
    #vTECChunked_anyAngleAvg = np.empty( avg_anyAngle_N  , dtype=np.float32); #faster at least
    kl = np.zeros( temp_pplat.size , np.bool_); #preallocate
    
    for j in prange(0,avg_anyAngle_N): #avg_anyAngle_N
        #average vTEC for a range chunk on an angle
        #kl = np.where(inpolgygon[j].contains_points( temp_pplatpplong ))[0];  #get pts inside the area defined
        #MATPLOTLIB'S INPOLYGON IS INCORRECT!
        #Gets pts inside the range
        
#        #P1PtVxy = np.array( (temp_pplong - temp_Long_List[j,0] , temp_pplat - temp_Lat_List[j,0]) ); #effort to do own in polygon
#        P1PtVx = temp_pplong - temp_Long_List[j,0]; #effort to do own in polygon
#        P1PtVy = temp_pplat - temp_Lat_List[j,0];
#        P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVx),np.abs(P1PtVy)) + sqrtApproxBo*np.minimum(np.abs(P1PtVx),np.abs(P1PtVy));
#        #P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat - temp_Lat_List[j,2]]); #effort to do own in polygon
#        P3PtVx = temp_pplong - temp_Long_List[j,2]; #effort to do own in polygon
#        P3PtVy = temp_pplat - temp_Lat_List[j,2];
#        P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVx),np.abs(P3PtVy)) + sqrtApproxBo*np.minimum(np.abs(P3PtVx),np.abs(P3PtVy));
#  
#        for l in range(0,temp_pplat.size):
#            P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVx[l] + P1P2Vxy[1,j]*P1PtVy[l])/(P1P2Vm[j]*P1PtVm[l]);
#            P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVx[l] + P1P4Vxy[1,j]*P1PtVy[l])/(P1P4Vm[j]*P1PtVm[l]);
#            P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVx[l] + P3P2Vxy[1,j]*P3PtVy[l])/(P3P2Vm[j]*P3PtVm[l]);
#            P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVx[l] + P3P4Vxy[1,j]*P3PtVy[l])/(P3P4Vm[j]*P3PtVm[l]);
#            kl[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
#        #END FOR l
        
        #one shot requires big matrix stuff to happen, a loop might help the time
#        P1P2P1PtdotCosT = np.einsum('ij,ij->j',np.tile(P1P2Vxy[:,j],[P1PtVx.size,1]).T, np.array((P1PtVx,P1PtVy)) )/(np.tile(P1P2Vm[j],[P1PtVx.size,])*P1PtVm);
#        P1P4P1PtdotCosT = np.einsum('ij,ij->j',np.tile(P1P4Vxy[:,j],[P1PtVx.size,1]).T, np.array((P1PtVx,P1PtVy)) )/(np.tile(P1P4Vm[j],[P1PtVx.size,])*P1PtVm);
#        P3P2P3PtdotCosT = np.einsum('ij,ij->j',np.tile(P3P2Vxy[:,j],[P1PtVx.size,1]).T, np.array((P3PtVx,P3PtVy)) )/(np.tile(P3P2Vm[j],[P1PtVx.size,])*P3PtVm);
#        P3P4P3PtdotCosT = np.einsum('ij,ij->j',np.tile(P3P4Vxy[:,j],[P1PtVx.size,1]).T, np.array((P3PtVx,P3PtVy)) )/(np.tile(P3P4Vm[j],[P1PtVx.size,])*P3PtVm);
#        kl = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);    
#        P1P2P1PtdotCosT = np.sum(np.tile(P1P2Vxy[:,j],[P1PtVx.size,1]).T.conj()*np.array((P1PtVx,P1PtVy)) ,axis=0)/(np.tile(P1P2Vm[j],[P1PtVx.size,])*P1PtVm);
#        P1P4P1PtdotCosT = np.sum(np.tile(P1P4Vxy[:,j],[P1PtVx.size,1]).T.conj()*np.array((P1PtVx,P1PtVy)) ,axis=0)/(np.tile(P1P4Vm[j],[P1PtVx.size,])*P1PtVm);
#        P3P2P3PtdotCosT = np.sum(np.tile(P3P2Vxy[:,j],[P1PtVx.size,1]).T.conj()*np.array((P3PtVx,P3PtVy)) ,axis=0)/(np.tile(P3P2Vm[j],[P1PtVx.size,])*P3PtVm);
#        P3P4P3PtdotCosT = np.sum(np.tile(P3P4Vxy[:,j],[P1PtVx.size,1]).T.conj()*np.array((P3PtVx,P3PtVy)) ,axis=0)/(np.tile(P3P4Vm[j],[P1PtVx.size,])*P3PtVm);
#        kl = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]); 
        
        for l in range(0,temp_pplat.size):
            P1PtVx = temp_pplong[l] - temp_Long_List[j,0]; #effort to do own in polygon
            P1PtVy = temp_pplat[l] - temp_Lat_List[j,0];
            P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVx),np.abs(P1PtVy)) + sqrtApproxBo*np.minimum(np.abs(P1PtVx),np.abs(P1PtVy));
            #P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat - temp_Lat_List[j,2]]); #effort to do own in polygon
            P3PtVx = temp_pplong[l] - temp_Long_List[j,2]; #effort to do own in polygon
            P3PtVy = temp_pplat[l] - temp_Lat_List[j,2];
            P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVx),np.abs(P3PtVy)) + sqrtApproxBo*np.minimum(np.abs(P3PtVx),np.abs(P3PtVy));
        
            P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVx + P1P2Vxy[1,j]*P1PtVy)/(P1P2Vm[j]*P1PtVm);
            P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVx + P1P4Vxy[1,j]*P1PtVy)/(P1P4Vm[j]*P1PtVm);
            P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVx + P3P2Vxy[1,j]*P3PtVy)/(P3P2Vm[j]*P3PtVm);
            P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVx + P3P4Vxy[1,j]*P3PtVy)/(P3P4Vm[j]*P3PtVm);
            kl[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
        #END FOR l
        
        if( np.sum(kl) != 0 ):
            vTECChunked_anyAngleAvg[j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
            #vTECChunked_anyAngleAvg[j] = np.sum(temp_vTEC[kl])/np.sum(kl);
        #END IF
        #vTECChunked_anyAngleAvg[j] = np.sum(kl); #check what the kl sums are (or number of points in the box found)
    #END FOR j
        
    return vTECChunked_anyAngleAvg; #donezo
    
    
    