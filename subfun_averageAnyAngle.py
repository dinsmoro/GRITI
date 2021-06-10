"""
Averages an area

data is array of data points to average (scattered) [numData , 3] sized - [0 is value, 1 is x location, 2 is y location]
angle is angle to average on
xRange is the x range of the data, same units as the x location values in data

globalMode means the xRange and yRange are actually longitudeRange and latitudeRange, respectively. That means latitudeRange is bounded [-90 to 90] and longitudeRange wraps at the -180 and 180 boundaries
timeToTime - 0 for no time, 1 for time data in data[numData,3] column
"""
import numpy as np
from subfun_raytrace import subfun_raytrace
from GRITI_TEC_avgAnyAngle_subfun_Raytrace import GRITI_TEC_avgAnyAngle_subfun_Raytrace
# from numba import jit

# @jit(nopython=True,nogil=True,parallel=True,fastmath=True)
def subfun_averageAnyAngle(data,angle,xRange,yRange,averageWidth,averageSplitNum,angleUnit=None,globalMode=0,timeToTime=0):    
    
    if( (angleUnit == 'degrees') | (angleUnit == 'Degrees') ):
        #AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
        if( (angle == 0) | (angle == 270) ):
            angle = angle + 0.0001; #deg, adjust so not exactly 0 or 270 (to be symmetrically consistent)
        elif( (angle == 90) | (angle == 180) ):
            angle = angle - 0.0001; #deg, adjust so not exactly 90 or 180 (to be symmetrically consistent)
        #END IF    
        angle = angle*np.pi/180; #rad, convert to radians
    else:
        #ironically it compares an infinite number as equal, don't worry
        #AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
        if( (angle == 0) | (angle == np.pi+np.pi/2) ):
            angle = angle + 0.0001*np.pi/180; #rad, adjust so not exactly 0 or 270 (to be symmetrically consistent)
        elif( (angle == np.pi/2) | (angle == np.pi) ):
            angle = angle - 0.0001*np.pi/180; #rad, adjust so not exactly 90 or 180 (to be symmetrically consistent)
        #END IF 
    #END IF
    
    angle_slope = np.tan(angle); #get slope of line required
    #this conversion is for y=LATITUDE x=LONGITUDE line action
    
    #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) 90 deg (real angle) to the req
    #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
    angle_upLine_int = averageWidth/2*np.sin(angle + np.pi/2) + np.mean(yRange)  \
        - angle_slope*(averageWidth/2*np.cos(angle + np.pi/2) + np.mean(xRange)); #get intercept of upper line
    #upper and lower lines are parallel
    #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) -90 deg (real angle) to the req
    #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
    angle_loLine_int = averageWidth/2*np.sin(angle - np.pi/2) + np.mean(yRange)  \
        - angle_slope*(averageWidth/2*np.cos(angle - np.pi/2) + np.mean(xRange)); #get intercept of lower line
    
    angle_LatLim_upLine = angle_slope*np.sort(xRange) + angle_upLine_int; #arcdeg,
    #latitude range from largest possible longitude range
    angle_LatLim_upLine[ angle_LatLim_upLine > np.max(yRange)] = np.max(yRange); #arcdeg, limit lat to current range
    angle_LatLim_upLine[ angle_LatLim_upLine < np.min(yRange)] = np.min(yRange); #arcdeg, limit lat to current range
    
    #Project upper line pts to lower pts
    angle_LatLim_loLine = np.array( [np.min(angle_LatLim_upLine) + averageWidth*np.sin(angle - np.pi/2) , np.max(angle_LatLim_upLine) + averageWidth*np.sin(angle - np.pi/2)] ); #arcdeg, calc lower pts based on upper
    angle_LatLim_loLine[ angle_LatLim_loLine > np.max(yRange)] = np.max(yRange); #arcdeg, limit lat to current range
    angle_LatLim_loLine[ angle_LatLim_loLine < np.min(yRange)] = np.min(yRange); #arcdeg, limit lat to current range
    #redefine upper pts based on possibly adjusted lower pts
    angle_LatLim_upLine = np.array( [np.min(angle_LatLim_loLine) + averageWidth*np.sin(angle + np.pi/2) , np.max(angle_LatLim_loLine) + averageWidth*np.sin(angle + np.pi/2)] ); #arcdeg, calc upper pts based on lower
    angle_LongLim_upLine = (angle_LatLim_upLine - angle_upLine_int)/angle_slope; #arcdeg, get longitudes that match
    angle_LongLim_loLine = (angle_LatLim_loLine - angle_loLine_int)/angle_slope; #arcdeg, get longitudes that match
    
    angle_Range = np.concatenate( [np.concatenate(  [angle_LatLim_upLine.T , angle_LatLim_loLine.T]  ).reshape(-1,1) , np.concatenate( [angle_LongLim_upLine , angle_LongLim_loLine] ).reshape(-1,1)] , axis = 1); #arcdeg, record pts for use
    
    if( globalMode == 1 ):
        angle_Range[np.where(angle_Range[:,1] > 180),1]  = 180.; #fix over 180 so it flips over to the other side
        angle_Range[np.where(angle_Range[:,1] < -180),1]  = -180.; #fix less than -180 so it flips over to the other side
        angle_Range[np.where(angle_Range[:,0] > 90),0]  = 90.; #fix over 90 so it flips over to the other side
        angle_Range[np.where(angle_Range[:,0] < -90),0]  = -90.; #fix less than -90 so it flips over to the other side
    #END IF
    
    angle_Range_Chunks_Long_up = np.linspace( np.min(angle_Range[0:2,1]) , np.max(angle_Range[0:2,1]) , averageSplitNum + 1); #arcdeg, chunks in longitude to go between
    #chose longitude because 0 deg will be stable - 90 deg would never be with
    #my hella math *wasn't stable at 0 anyway lol*
    angle_Range_Chunks_Long_lo = np.linspace( np.min(angle_Range[2:4,1]) , np.max(angle_Range[2:4,1]) , averageSplitNum + 1); #arcdeg, chunks in longitude to go between
    # angle_Range_Chunks_Lat = linspace( min(angle_Range(1:2,1)) , max(angle_Range(1:2,1)) , averageSplitNum + 1)'; %arcdeg, chunks in latitude to go between
    
    temp_Longs_up = np.zeros( [averageSplitNum,2] ); #preallocate
    temp_Lats_up = np.zeros( [averageSplitNum,2] ); #preallocate
    temp_Longs_lo = np.zeros( [averageSplitNum,2] ); #preallocate
    temp_Lats_lo = np.zeros( [averageSplitNum,2] ); #preallocate
    for j in range(0,averageSplitNum): #preallocate and fill
        temp_Longs_up[j,:] = [angle_Range_Chunks_Long_up[j] , angle_Range_Chunks_Long_up[j+1]]; #arcdeg, get longitudes needed upper line
        temp_Longs_lo[j,:] = np.flip( np.array( [angle_Range_Chunks_Long_lo[j] , angle_Range_Chunks_Long_lo[j+1]] ) ); #arcdeg, get longitudes needed low
        temp_Lats_up[j,:] = angle_slope*temp_Longs_up[j,:] + angle_upLine_int; #arcdeg, get latitudes needed up
        temp_Lats_lo[j,:] = angle_slope*temp_Longs_lo[j,:] + angle_loLine_int; #arcdeg, get latitudes needed lower line
    #END IF
    
    temp_Long_List = np.zeros( [averageSplitNum,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
    temp_Lat_List = np.zeros( [averageSplitNum,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
    for j in range(0,averageSplitNum):
        temp_Long_List[j,:] = np.hstack( [temp_Longs_up[j,:],temp_Longs_lo[j,:],temp_Longs_up[j,0]] ); #so this isn't done dynamtically
        temp_Lat_List[j,:] = np.hstack( [temp_Lats_up[j,:],temp_Lats_lo[j,:],temp_Lats_up[j,0]] );
    #END IF
    
    
    #limit the work needed to calc the stuff by making sure only using data inside the grid
    sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
    sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
    P1P2Vxy = np.array([angle_Range[1,1] - angle_Range[0,1] , angle_Range[1,0] - angle_Range[0,0]]); #effort to do own in polygon
    P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1]));
    #P1P2Vm = np.sqrt( (P1P2Vxy[0])**2 + (P1P2Vxy[1])**2 ); #no approx for magnitude used
    P1P4Vxy = np.array([angle_Range[2,1] - angle_Range[0,1] , angle_Range[2,0] - angle_Range[0,0]]); #effort to do own in polygon
    P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1]));
    #P1P4Vm = np.sqrt( (P1P4Vxy[0])**2 + (P1P4Vxy[1])**2 ); #no approx for magnitude used
    P3P2Vxy = np.array([angle_Range[1,1] - angle_Range[3,1] , angle_Range[1,0] - angle_Range[3,0]]); #effort to do own in polygon
    P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1]));
    #P3P2Vm = np.sqrt( (P3P2Vxy[0])**2 + (P3P2Vxy[1])**2 ); #no approx for magnitude used
    P3P4Vxy = np.array([angle_Range[2,1] - angle_Range[3,1] , angle_Range[2,0] - angle_Range[3,0]]); #effort to do own in polygon
    P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1]));
    #P3P4Vm = np.sqrt( (P3P4Vxy[0])**2 + (P3P4Vxy[1])**2 ); #no approx for magnitude used
    P1P2P1P4dotCosT = np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
    P3P2P3P4dotCosT = np.sum(P3P2Vxy.conj()*P3P4Vxy,axis=0)/(P3P2Vm*P3P4Vm);
    
    P1PtVxy = np.array([data[:,1] - angle_Range[0,1] , data[:,2] - angle_Range[0,0]]); #effort to do own in polygon
    P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
    #P1PtVm = np.sqrt( (P1PtVxy[0,:])**2 + (P1PtVxy[1,:])**2 ); #no approx for magnitude used
    P3PtVxy = np.array([data[:,1] - angle_Range[3,1] , data[:,2] - angle_Range[3,0]]); #effort to do own in polygon
    P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:]));
    #P3PtVm = np.sqrt( (P3PtVxy[0,:])**2 + (P3PtVxy[1,:])**2 ); #no approx for magnitude used
    
    #not enough RAM for this! uh oh
    P1P2P1PtdotCosT = (P1P2Vxy[0]*P1PtVxy[0,:] + P1P2Vxy[1]*P1PtVxy[1,:])/(P1P2Vm*P1PtVm);
    P1P4P1PtdotCosT = (P1P4Vxy[0]*P1PtVxy[0,:] + P1P4Vxy[1]*P1PtVxy[1,:])/(P1P4Vm*P1PtVm);
    P3P2P3PtdotCosT = (P3P2Vxy[0]*P3PtVxy[0,:] + P3P2Vxy[1]*P3PtVxy[1,:])/(P3P2Vm*P3PtVm);
    P3P4P3PtdotCosT = (P3P4Vxy[0]*P3PtVxy[0,:] + P3P4Vxy[1]*P3PtVxy[1,:])/(P3P4Vm*P3PtVm);
    keepr = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT);    
    #sys.exit()
    
#    keepr = np.zeros((P3PtVxy[1,:].size,),dtype=np.uint8); #preallocate    
#    for j in range(0,P3PtVxy[1,:].size):
#        P1P2P1PtdotCosT = (P1P2Vxy[0]*P1PtVxy[0,j] + P1P2Vxy[1]*P1PtVxy[1,j])/(P1P2Vm*P1PtVm[j]);
#        P1P4P1PtdotCosT = (P1P4Vxy[0]*P1PtVxy[0,j] + P1P4Vxy[1]*P1PtVxy[1,j])/(P1P4Vm*P1PtVm[j]);
#        P3P2P3PtdotCosT = (P3P2Vxy[0]*P3PtVxy[0,j] + P3P2Vxy[1]*P3PtVxy[1,j])/(P3P2Vm*P3PtVm[j]);
#        P3P4P3PtdotCosT = (P3P4Vxy[0]*P3PtVxy[0,j] + P3P4Vxy[1]*P3PtVxy[1,j])/(P3P4Vm*P3PtVm[j]);
#        
#        keepr[j] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT);    
#    #END FOR j
                
    #keepr = pplong < max(max(temp_Long_List)) &  pplong > min(min(temp_Long_List)) & pplat < max(max(temp_Lat_List)) & pplat > min(min(temp_Lat_List)); %limit memory sent to parallel workers
    #keepr = inpolygon(pplong,pplat, [angle_Range(1:2,2) ; flipud(angle_Range(3:4,2)) ; angle_Range(1,2)] , [angle_Range(1:2,1) ; flipud(angle_Range(3:4,1)) ; angle_Range(1,1)]);
    if( timeToTime == 1 ):
        time_limd =  data[keepr,3];
    #END IF
    data_limd =  data[keepr,0].T;
    pplat_limd = data[keepr,2];
    pplong_limd = data[keepr,1];
    
    #sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
    #sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
    P1P2Vxy = np.array([temp_Long_List[:,1] - temp_Long_List[:,0] , temp_Lat_List[:,1] - temp_Lat_List[:,0]]); #effort to do own in polygon
    P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0,:]),np.abs(P1P2Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0,:]),np.abs(P1P2Vxy[1,:]));
    P1P4Vxy = np.array([temp_Long_List[:,3] - temp_Long_List[:,0] , temp_Lat_List[:,3] - temp_Lat_List[:,0]]); #effort to do own in polygon
    P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0,:]),np.abs(P1P4Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0,:]),np.abs(P1P4Vxy[1,:]));
    P3P2Vxy = np.array([temp_Long_List[:,1] - temp_Long_List[:,2] , temp_Lat_List[:,1] - temp_Lat_List[:,2]]); #effort to do own in polygon
    P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0,:]),np.abs(P3P2Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0,:]),np.abs(P3P2Vxy[1,:]));
    P3P4Vxy = np.array([temp_Long_List[:,3] - temp_Long_List[:,2] , temp_Lat_List[:,3] - temp_Lat_List[:,2]]); #effort to do own in polygon
    P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0,:]),np.abs(P3P4Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0,:]),np.abs(P3P4Vxy[1,:]));
    P1P2P1P4dotCosT = np.einsum('ij,ij->j',P1P2Vxy,P1P4Vxy)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
    P3P2P3P4dotCosT = np.einsum('ij,ij->j',P3P2Vxy,P3P4Vxy)/(P3P2Vm*P3P4Vm);
     
    if( timeToTime == 0 ):
        data_averaged = subfun_raytrace(pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,data_limd,averageSplitNum,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT); #avg without moving in time - e.g. only avg 1 image
    else:
        #this one supports time
        data_averaged = GRITI_TEC_avgAnyAngle_subfun_Raytrace(time_limd,np.unique(data[:,3]),pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,data_limd,averageSplitNum,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT); #call a function to avg all the TEC together into bands
    #END IF
    
    delta = np.median(np.sqrt(np.diff(temp_Longs_up)**2 + np.diff(temp_Lats_up)**2)); #degc, get the delta between averaging points
    
    return data_averaged, delta