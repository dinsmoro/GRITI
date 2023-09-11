#Function to ray trace real fast
#RD on 12/19/2018

#Ray tracing is happening because of stuff
import cupy as np

def GRITI_keo_subfun_simpleTrace_lat_GPU(time_limd,TEC_timeUnique,pplat_limd,pplong_limd,splits,longLims,vTEC_limd,avg_anyAngle_N):
    vTECChunked_anyAngleAvg = np.empty( (TEC_timeUnique.size,avg_anyAngle_N) ,dtype=np.float32 )*np.nan; #preallocate, don't care what's in it and nan it to be sure
    
    # longLims_max = np.float32(np.max(longLims));
    # longLims_min = np.float32(np.min(longLims));
    
    longLims_logical = pplong_limd <= np.float32(np.max(longLims));
    longLims_logical[longLims_logical] = pplong_limd[longLims_logical] >= np.float32(np.min(longLims));
    
    vTEC_limd_nonans = np.logical_not(np.isnan(vTEC_limd)); #get nan locations once so don't work on NaNs at all (saves work for NaN data)
    
    time_limd_deux = time_limd[vTEC_limd_nonans & longLims_logical]; #pre-limit
    vTEC_limd_deux = vTEC_limd[vTEC_limd_nonans & longLims_logical]; #pre-limit
    pplat_limd_deux = pplat_limd[vTEC_limd_nonans & longLims_logical]; #pre-limit
    for i in range(0, TEC_timeUnique.size ): #run through each unique time period
        #Corral the data to the right place  
        # k = np.where((time_limd == TEC_timeUnique[i]) & vTEC_limd_nonans & longLims_logical)[0]; #gets during a time period
        # k = np.where(time_limd_deux == TEC_timeUnique[i])[0]; #gets during a time period
        
        temp_vTEC = vTEC_limd_deux[time_limd_deux == TEC_timeUnique[i]]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        temp_pplat = pplat_limd_deux[time_limd_deux == TEC_timeUnique[i]];
        # temp_pplong = pplong_limd[time_limd_deux == TEC_timeUnique[i]];
        
        for j in range(0,avg_anyAngle_N): #run through each segment
            vTECChunked_anyAngleAvg[i,j] = np.mean(temp_vTEC[(temp_pplat <= splits[j+1]) & (temp_pplat >= splits[j])]); #average the vTEC in this lat band with the given long range
        #END FOR j
    #END FOR i
        
    return vTECChunked_anyAngleAvg; #donezo
#END DEF