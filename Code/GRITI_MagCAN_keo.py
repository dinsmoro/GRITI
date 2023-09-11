#GOAL: Run any angle averaging strips on TEC data in requested area, also plot visualization of averaging area
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: vTECChunked_anyAngleAvg, keo_angle, keo_width,  keo_range_chunks_long_plot, keo_range_chunks_long_plot_name

import numpy as np #import in here I dunno
# import os
import time
from scipy.signal import savgol_filter
from Code.subfun_filter import subfun_filter
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.ticker import FormatStrFormatter

from Code.GRITI_MagCAN_keo_area import GRITI_MagCAN_keo_area
# from Code.GRITI_MagCAN_keo_fancyPlot_area import GRITI_MagCAN_keo_fancyPlot_area
from Code.GRITI_keo_subfun_raytrace import GRITI_keo_subfun_raytrace_cosApprox


def GRITI_MagCAN_keo(data,dates,settings,FLG_disablePlot=0):
    
    #==============UNPACK SETTINGS==============
    keo_angle = np.copy(settings['MagCAN']['keo angle']).item(); #make a copy just in case python memory stuff
    keo_width = np.copy(settings['MagCAN']['keo width orig']).item(); #make a copy just in case python memory stuff
    keo_N = settings['MagCAN']['keo N'];
    keo_45vsLatLong = settings['MagCAN']['keo 45 lat/long'];
    map_plotLatRange = settings['MagCAN']['lat range'];
    map_plotLongRange = settings['MagCAN']['long range'];
    FLG_fancyPlot = settings['plot']['fancy plot'];
    
    #==============PREP DATA TO BE IN RIGHT FORM==============
    if( settings['MagCAN']['keo set stations'] != 1 ):
        siteNames = data['MagCAN']['site names']; #get the site names
    else:
        siteNames = settings['MagCAN']['keo set stations names']; #set the site names to whatever the user had
    #END IF
    removeList = []; #prep
    for j in range(0,len(siteNames)):
        if( ~((data['MagCAN'][siteNames[j]]['lat'] >= np.min(settings['MagCAN']['lat range'])) & \
                (data['MagCAN'][siteNames[j]]['lat'] <= np.max(settings['MagCAN']['lat range'])) & \
                (data['MagCAN'][siteNames[j]]['long'] >= np.min(settings['MagCAN']['long range'])) & \
                (data['MagCAN'][siteNames[j]]['long'] <= np.max(settings['MagCAN']['long range']))) ):
            removeList.append(siteNames[j]); #prep for removal b/c not within the range we want
        #END IF
    #END FOR j
    #solving lists takes a lot of slowww lists ohw ell
    for j in range(0,len(removeList)):
        siteNames.remove(removeList[j]); #remove the stuff we don't need
    #END FOR j
    
    cntr = 0; #prep cntr
    for j in range(0,len(siteNames)):
        cntr += data['MagCAN'][siteNames[j]]['MagCANF'].size; #get the total size
    #END FOR j
    MagCAN_time = np.zeros( cntr, dtype=np.float64); #preallocate
    MagCAN_deltaMagCANF = np.zeros( cntr, dtype=np.float64); #preallocate
    MagCAN_lat = np.zeros( cntr, dtype=np.float64); #preallocate
    MagCAN_long = np.zeros( cntr, dtype=np.float64); #preallocate
    cntr = 0; #Prep cntr
    maxRec = 0; #prep max recorder
    for j in range(0,len(siteNames)):
        tempSize = data['MagCAN'][siteNames[j]]['MagCANF'].size; #prep
        MagCAN_time[cntr:cntr+tempSize] = np.int64(data['MagCAN'][siteNames[j]]['dayNumF'])*86400 + np.int64(data['MagCAN'][siteNames[j]]['secF']); #sec, get unique times (v useful)
        #filter function here
        temp_deltaMagCANF = subfun_filter( data['MagCAN'][siteNames[j]]['MagCANF'], data['MagCAN'][siteNames[j]]['secF'], settings['MagCAN']['delta method'], settings['spectra'], dataRate = data['MagCAN'][siteNames[j]]['dataRateF'], reduceWindow = 0); #filter using the all-powerful filter function
        
        MagCAN_deltaMagCANF[cntr:cntr+tempSize] = temp_deltaMagCANF; #write data
        MagCAN_lat[cntr:cntr+tempSize] = np.tile(data['MagCAN'][siteNames[j]]['lat'], tempSize); #tile data (lat doesn't change)
        MagCAN_long[cntr:cntr+tempSize] = np.tile(data['MagCAN'][siteNames[j]]['long'], tempSize); #tile data (long doesn't change)
        cntr += tempSize; #increment
        if( np.nanmax(np.abs(temp_deltaMagCANF)) > maxRec ):
            maxRec = np.nanmax(np.abs(temp_deltaMagCANF)); #set it
        #END IF
    #END FOR j
    if( settings['MagCAN']['keo normalize'] == 1 ):
        cntr = 0; #Prep cntr
        for j in range(0,len(siteNames)):
            tempSize = data['MagCAN'][siteNames[j]]['MagCANF'].size; #prep
            MagCAN_deltaMagCANF[cntr:cntr+tempSize] = MagCAN_deltaMagCANF[cntr:cntr+tempSize]*(maxRec/np.nanmax(np.abs(MagCAN_deltaMagCANF[cntr:cntr+tempSize]))); #normalize
            #above is better than power normalizing
            cntr += tempSize; #increment
        #END FOR j
    #END IF
    
    #==============Analysis: Any Angle AVG (Keogram)==============
    print('Analyis: MagCANnetometer Keogram Beginning:'); 
    tic = time.time();
    
    #AVG any angle also can't do exactly the full width (so a width that takes the whole plot area)
    if( ( (np.round(keo_angle) == 0) | (np.round(keo_angle) == 180) ) & (keo_width >= np.abs(map_plotLatRange[0] - map_plotLatRange[1])) ):
        if( keo_width > np.abs(map_plotLatRange[0] - map_plotLatRange[1]) ):
            print("\n==============~Warning~==============");
            print('Keogram width of '+str(keo_width)+' arcdeg is larger than the latitudinal plot area of '+str(np.abs(map_plotLatRange[1] - map_plotLatRange[0]))+' arcdeg. Reducing to be size of latitudinal plot area.');
            keo_width = np.abs(map_plotLatRange[0] - map_plotLatRange[1]) - 0.001; #arcdeg, adjust so not exactly plot width, but basically since overstepped the range already
        else:
            keo_width = keo_width - 0.001; #arcdeg, adjust so not exactly plot width
        #END IF
    elif( ( (np.round(keo_angle) == 90) | (np.round(keo_angle) == 270) ) & (keo_width >= np.abs(map_plotLongRange[0] - map_plotLongRange[1])) ):
        if( keo_width > np.abs(map_plotLongRange[0] - map_plotLongRange[1]) ):
            print("\n==============~Warning~==============");
            print('Keogram of '+str(keo_width)+' arcdeg is larger than the longitudinal plot area of '+str(np.abs(map_plotLongRange[0] - map_plotLongRange[1]))+' arcdeg. Reducing to be size of longitudinal plot area.');
            keo_width = np.abs(map_plotLongRange[0] - map_plotLongRange[1]) - 0.001; #arcdeg, adjust so not exactly plot width, but just under it since overstepped the range already
        else:
            keo_width = keo_width - 0.001; #arcdeg, adjust so not exactly plot width
        #END IF
    #END IF
        
    #AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
    if( (keo_angle == 0) | (keo_angle == 270) ):
        keo_angle = keo_angle + 0.0001; #deg, adjust so not exactly 0 or 270 (to be symmetrically consistent)
    elif( (keo_angle == 90) | (keo_angle == 180) ):
        keo_angle = keo_angle - 0.0001; #deg, adjust so not exactly 90 or 180 (to be symmetrically consistent)
    #END IF    
    
    keo_angleRad = keo_angle*np.pi/180; #rad, convert to radians
    
    keo_slope = np.tan(keo_angleRad); #get slope of line required
    #this conversion is for y=LATITUDE x=LONGITUDE line action
    
    #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) 90 deg (real angle) to the req
    #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
    keo_upLine_int = keo_width/2*np.sin(keo_angleRad + np.pi/2) + np.mean(map_plotLatRange)  \
        - keo_slope*(keo_width/2*np.cos(keo_angleRad + np.pi/2) + np.mean(map_plotLongRange)); #get intercept of upper line
    #upper and lower lines are parallel
    #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) -90 deg (real angle) to the req
    #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
    keo_loLine_int = keo_width/2*np.sin(keo_angleRad - np.pi/2) + np.mean(map_plotLatRange)  \
        - keo_slope*(keo_width/2*np.cos(keo_angleRad - np.pi/2) + np.mean(map_plotLongRange)); #get intercept of lower line
    
    keo_LatLim_upLine = keo_slope*np.sort(map_plotLongRange) + keo_upLine_int; #arcdeg,
    #latitude range from largest possible longitude range
    keo_LatLim_upLine[ keo_LatLim_upLine > np.max(map_plotLatRange)] = np.max(map_plotLatRange); #arcdeg, limit lat to current range
    keo_LatLim_upLine[ keo_LatLim_upLine < np.min(map_plotLatRange)] = np.min(map_plotLatRange); #arcdeg, limit lat to current range
    
    #Project upper line pts to lower pts
    keo_LatLim_loLine = np.array( [np.min(keo_LatLim_upLine) + keo_width*np.sin(keo_angleRad - np.pi/2) , np.max(keo_LatLim_upLine) + keo_width*np.sin(keo_angleRad - np.pi/2)] ); #arcdeg, calc lower pts based on upper
    keo_LatLim_loLine[ keo_LatLim_loLine > np.max(map_plotLatRange)] = np.max(map_plotLatRange); #arcdeg, limit lat to current range
    keo_LatLim_loLine[ keo_LatLim_loLine < np.min(map_plotLatRange)] = np.min(map_plotLatRange); #arcdeg, limit lat to current range
    #redefine upper pts based on possibly adjusted lower pts
    keo_LatLim_upLine = np.array( [np.min(keo_LatLim_loLine) + keo_width*np.sin(keo_angleRad + np.pi/2) , np.max(keo_LatLim_loLine) + keo_width*np.sin(keo_angleRad + np.pi/2)] ); #arcdeg, calc upper pts based on lower
    keo_LongLim_upLine = (keo_LatLim_upLine - keo_upLine_int)/keo_slope; #arcdeg, get longitudes that match
    keo_LongLim_loLine = (keo_LatLim_loLine - keo_loLine_int)/keo_slope; #arcdeg, get longitudes that match
    
    keo_range = np.concatenate( [np.concatenate(  [keo_LatLim_upLine.T , keo_LatLim_loLine.T]  ).reshape(-1,1) , np.concatenate( [keo_LongLim_upLine , keo_LongLim_loLine] ).reshape(-1,1)] , axis = 1); #arcdeg, record pts for use
    
    keo_range[np.where(keo_range[:,1] > 180),1]  = 180.; #fix over 180 so it flips over to the other side
    keo_range[np.where(keo_range[:,1] < -180),1]  = -180.; #fix less than -180 so it flips over to the other side
    keo_range[np.where(keo_range[:,0] > 90),0]  = 90.; #fix over 90 so it flips over to the other side
    keo_range[np.where(keo_range[:,0] < -90),0]  = -90.; #fix less than -90 so it flips over to the other side
    
    keo_range_chunks_long_up = np.linspace( np.min(keo_range[0:2,1]) , np.max(keo_range[0:2,1]) , keo_N + 1); #arcdeg, chunks in longitude to go between
    #chose longitude because 0 deg will be stable - 90 deg would never be with
    #my hella math *wasn't stable at 0 anyway lol*
    keo_range_chunks_long_lo = np.linspace( np.min(keo_range[2:4,1]) , np.max(keo_range[2:4,1]) , keo_N + 1); #arcdeg, chunks in longitude to go between
    # keo_range_Chunks_Lat = linspace( min(keo_range(1:2,1)) , max(keo_range(1:2,1)) , keo_N + 1)'; %arcdeg, chunks in latitude to go between
    
    keo_range_chunks_long_plot = np.linspace( (keo_range_chunks_long_up[0] + keo_range_chunks_long_lo[0])/2 ,\
        (keo_range_chunks_long_up[-1] + keo_range_chunks_long_lo[-1])/2 , keo_N+1); #for plotting, the mean between each point
    #for plotting using up
    
    if( (keo_45vsLatLong == 1) & ((keo_angleRad == 45*np.pi/180) | (keo_angleRad == 135*np.pi/180)) ): #override mechanism to allow LATITUDE when normally longitude is defaulted to
        
        keo_range_chunks_long_plot_name = 'Latitude'; #name of Latitude
        
        keo_range_chunks_long_plot_int = np.mean(map_plotLatRange)  \
            - keo_slope*np.mean(map_plotLongRange); #get intercept of upper line
        
        keo_range_chunks_long_plot = keo_slope*keo_range_chunks_long_plot + keo_range_chunks_long_plot_int;
        #for plotting rocking it up
        
    elif( (keo_angleRad <= 45*np.pi/180) | (keo_angleRad >= 135*np.pi/180) ): #actually LONGITUDE on the axis
    
        keo_range_chunks_long_plot_name = 'Longitude'; #name of longitude
    
    else: #otherwise LATITUDE on the axis
        keo_range_chunks_long_plot_name = 'Latitude'; #name of Latitude
        
        keo_range_chunks_long_plot_int = np.mean(map_plotLatRange)  \
            - keo_slope*np.mean(map_plotLongRange); #get intercept of upper line
        
        keo_range_chunks_long_plot = keo_slope*keo_range_chunks_long_plot + keo_range_chunks_long_plot_int;
        #for plotting rocking it up
    #END IF
    #keo_range_chunks_long_plot = np.round(keo_range_chunks_long_plot,decimals=5); #keeps things from getting out of hand
    
    temp_longs_up = np.zeros( [keo_N,2] ); #preallocate
    temp_lats_up = np.zeros( [keo_N,2] ); #preallocate
    temp_Longs_lo = np.zeros( [keo_N,2] ); #preallocate
    temp_lats_lo = np.zeros( [keo_N,2] ); #preallocate
    for j in range(0,keo_N): #preallocate and fill
        temp_longs_up[j,:] = [keo_range_chunks_long_up[j] , keo_range_chunks_long_up[j+1]]; #arcdeg, get longitudes needed upper line
        temp_Longs_lo[j,:] = np.flip( np.array( [keo_range_chunks_long_lo[j] , keo_range_chunks_long_lo[j+1]] ) ); #arcdeg, get longitudes needed low
        temp_lats_up[j,:] = keo_slope*temp_longs_up[j,:] + keo_upLine_int; #arcdeg, get latitudes needed up
        temp_lats_lo[j,:] = keo_slope*temp_Longs_lo[j,:] + keo_loLine_int; #arcdeg, get latitudes needed lower line
    #END IF
    
    temp_long_list = np.zeros( [keo_N,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
    temp_lat_list = np.zeros( [keo_N,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
    for j in range(0,keo_N):
        temp_long_list[j,:] = np.hstack( [temp_longs_up[j,:],temp_Longs_lo[j,:],temp_longs_up[j,0]] ); #so this isn't done dynamtically
        temp_lat_list[j,:] = np.hstack( [temp_lats_up[j,:],temp_lats_lo[j,:],temp_lats_up[j,0]] );
    #END IF
    
    # #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    # if( np.isscalar(TEC_plotLimValu) == 1 ):
    #     TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    # #END IF
    
    if( FLG_disablePlot == 0 ):
        #THIS IS TO VIEW DATA BEING AVERAGED
        GRITI_MagCAN_keo_area(data, dates, settings, keo_width, keo_range, temp_lats_up, temp_longs_up, temp_lat_list, temp_long_list); #plot the area
    #END IF
    if( FLG_fancyPlot == 1 ):
        #THIS IS TO VIEW DATA BEING AVERAGED ~fancily~
        pass;
        # GRITI_MagCAN_keo_fancyPlot_area(data, dates, settings, keo_range, temp_lats_up, temp_longs_up, temp_Lat_List, temp_Long_List); #plot the area
    #END IF

    #limit the work needed to calc the stuff by making sure only using data inside the grid
    sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
    sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
    P1P2Vxy = np.array([keo_range[1,1] - keo_range[0,1] , keo_range[1,0] - keo_range[0,0]]); #effort to do own in polygon
    P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1]));
    #P1P2Vm = np.sqrt( (P1P2Vxy[0])**2 + (P1P2Vxy[1])**2 ); #no approx for MagCANnitude used
    P1P4Vxy = np.array([keo_range[2,1] - keo_range[0,1] , keo_range[2,0] - keo_range[0,0]]); #effort to do own in polygon
    P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1]));
    #P1P4Vm = np.sqrt( (P1P4Vxy[0])**2 + (P1P4Vxy[1])**2 ); #no approx for MagCANnitude used
    P3P2Vxy = np.array([keo_range[1,1] - keo_range[3,1] , keo_range[1,0] - keo_range[3,0]]); #effort to do own in polygon
    P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1]));
    #P3P2Vm = np.sqrt( (P3P2Vxy[0])**2 + (P3P2Vxy[1])**2 ); #no approx for MagCANnitude used
    P3P4Vxy = np.array([keo_range[2,1] - keo_range[3,1] , keo_range[2,0] - keo_range[3,0]]); #effort to do own in polygon
    P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1]));
    #P3P4Vm = np.sqrt( (P3P4Vxy[0])**2 + (P3P4Vxy[1])**2 ); #no approx for MagCANnitude used
    P1P2P1P4dotCosT = np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
    P3P2P3P4dotCosT = np.sum(P3P2Vxy.conj()*P3P4Vxy,axis=0)/(P3P2Vm*P3P4Vm);
    
    P1PtVxy = np.array([MagCAN_long - keo_range[0,1] , MagCAN_lat - keo_range[0,0]]); #effort to do own in polygon
    P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
    #P1PtVm = np.sqrt( (P1PtVxy[0,:])**2 + (P1PtVxy[1,:])**2 ); #no approx for MagCANnitude used
    P3PtVxy = np.array([MagCAN_long - keo_range[3,1] , MagCAN_lat - keo_range[3,0]]); #effort to do own in polygon
    P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:]));
    #P3PtVm = np.sqrt( (P3PtVxy[0,:])**2 + (P3PtVxy[1,:])**2 ); #no approx for MagCANnitude used
    
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
            
    tic = time.time(); #get current time (it's like from 1970 or whatever time)
    
    #keepr = pplong < max(max(temp_long_list)) &  pplong > min(min(temp_long_list)) & pplat < max(max(temp_lat_list)) & pplat > min(min(temp_lat_list)); %limit memory sent to parallel workers
    #keepr = inpolygon(pplong,pplat, [keo_range(1:2,2) ; flipud(keo_range(3:4,2)) ; keo_range(1,2)] , [keo_range(1:2,1) ; flipud(keo_range(3:4,1)) ; keo_range(1,1)]);
    MagCAN_time =  MagCAN_time[keepr];
    MagCAN_deltaMagCANF =  MagCAN_deltaMagCANF[keepr];
    MagCAN_lat = MagCAN_lat[keepr];
    MagCAN_long = MagCAN_long[keepr];
    MagCAN_timeUnique = np.unique(MagCAN_time); #days, get the unique times for the conglomeration
    
    #sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
    #sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
    P1P2Vxy = np.array([temp_long_list[:,1] - temp_long_list[:,0] , temp_lat_list[:,1] - temp_lat_list[:,0]]); #effort to do own in polygon
    P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0,:]),np.abs(P1P2Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0,:]),np.abs(P1P2Vxy[1,:]));
    P1P4Vxy = np.array([temp_long_list[:,3] - temp_long_list[:,0] , temp_lat_list[:,3] - temp_lat_list[:,0]]); #effort to do own in polygon
    P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0,:]),np.abs(P1P4Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0,:]),np.abs(P1P4Vxy[1,:]));
    P3P2Vxy = np.array([temp_long_list[:,1] - temp_long_list[:,2] , temp_lat_list[:,1] - temp_lat_list[:,2]]); #effort to do own in polygon
    P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0,:]),np.abs(P3P2Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0,:]),np.abs(P3P2Vxy[1,:]));
    P3P4Vxy = np.array([temp_long_list[:,3] - temp_long_list[:,2] , temp_lat_list[:,3] - temp_lat_list[:,2]]); #effort to do own in polygon
    P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0,:]),np.abs(P3P4Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0,:]),np.abs(P3P4Vxy[1,:]));
    P1P2P1P4dotCosT = np.einsum('ij,ij->j',P1P2Vxy,P1P4Vxy)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
    P3P2P3P4dotCosT = np.einsum('ij,ij->j',P3P2Vxy,P3P4Vxy)/(P3P2Vm*P3P4Vm);
      
    MagCAN_keo = GRITI_keo_subfun_raytrace_cosApprox(MagCAN_time,MagCAN_timeUnique,MagCAN_lat,MagCAN_long,temp_lat_list,temp_long_list,MagCAN_deltaMagCANF,keo_N,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT); #call a function to avg all the TEC together into bands
    
#    #original, in-fun way (slower w/o jit etc.)
#        #pre-calc the inpolygon path stuff
#    inpolgygon = []; #prep a list
#    for j in range( 0, keo_N):
#        inpolgygon.append(Path([(temp_long_list[j,0],temp_lat_list[j,1]), (temp_long_list[j,0], temp_lat_list[j,2]), \
#                (temp_long_list[j,1], temp_lat_list[j,2]), (temp_long_list[j,1], temp_lat_list[j,2])]));  # square with legs length 1 and bottom left corner at the origin
#    #END FOR j
#    
#    finAnnounce_percentToUpdateAt = 0.5; # every % to update info at
#    finAnnounce_div = np.round(100/finAnnounce_percentToUpdateAt); #calc divisor to use
#    finAnnounce_modRef = np.round(len(TEC_timeUnique)/finAnnounce_div); #calc mod reference to be used
#    Ref1 = np.zeros( [len(TEC_timeUnique),keo_N] ,dtype=np.float32 ); #preallocate
#    for i in range(0, len(TEC_timeUnique) ): #87
#        #Corral the data to the right place  
#        k = np.where(time_limd == TEC_timeUnique[i]); #gets during a time period
#        
#        temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
#        #temp_pplatpplong = np.vstack( [pplong_limd[k], pplat_limd[k]] ).T; #record the pierce-point lat/long for the time in a way the Path function can use
#        temp_pplat = pplat_limd[k];
#        temp_pplong = pplong_limd[k];
#        
#        #put in a function for jit power
#        for j in range(0,keo_N): #keo_N
#            #average vTEC for a range chunk on an angle
#            #kl = np.where(inpolgygon[j].contains_points( temp_pplatpplong ))[0];  #get pts inside the area defined
#            #MATPLOTLIB'S INPOLYGON IS INCORRECT!
#            #Gets pts inside the range
#            
#            P1PtVxy = np.array([temp_pplong - temp_long_list[j,0] , temp_pplat - temp_lat_list[j,0]]); #effort to do own in polygon
#            P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
#            P3PtVxy = np.array([temp_pplong - temp_long_list[j,2] , temp_pplat - temp_lat_list[j,2]]); #effort to do own in polygon
#            P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:]));
#            
##            #one shot requires big matrix stuff to happen, a loop might help the time
##            P1P2P1PtdotCosT = np.einsum('ij,ij->j',np.tile(P1P2Vxy[:,j],[P1PtVxy.shape[1],1]).T,P1PtVxy)/(np.tile(P1P2Vm[j],[P1PtVxy.shape[1],])*P1PtVm);
##            P1P4P1PtdotCosT = np.einsum('ij,ij->j',np.tile(P1P4Vxy[:,j],[P1PtVxy.shape[1],1]).T,P1PtVxy)/(np.tile(P1P4Vm[j],[P1PtVxy.shape[1],])*P1PtVm);
##            P3P2P3PtdotCosT = np.einsum('ij,ij->j',np.tile(P3P2Vxy[:,j],[P3PtVxy.shape[1],1]).T,P3PtVxy)/(np.tile(P3P2Vm[j],[P3PtVxy.shape[1],])*P3PtVm);
##            P3P4P3PtdotCosT = np.einsum('ij,ij->j',np.tile(P3P4Vxy[:,j],[P3PtVxy.shape[1],1]).T,P3PtVxy)/(np.tile(P3P4Vm[j],[P3PtVxy.shape[1],])*P3PtVm);
##            kl = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);    
#            
#            kl = np.zeros((P3PtVxy[1,:].size,),dtype=np.bool_); #preallocate    
#            for l in range(0,P3PtVxy[1,:].size):
#                P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVxy[0,l] + P1P2Vxy[1,j]*P1PtVxy[1,l])/(P1P2Vm[j]*P1PtVm[l]);
#                P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVxy[0,l] + P1P4Vxy[1,j]*P1PtVxy[1,l])/(P1P4Vm[j]*P1PtVm[l]);
#                P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVxy[0,l] + P3P2Vxy[1,j]*P3PtVxy[1,l])/(P3P2Vm[j]*P3PtVm[l]);
#                P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVxy[0,l] + P3P4Vxy[1,j]*P3PtVxy[1,l])/(P3P4Vm[j]*P3PtVm[l]);
#                kl[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
#            #END FOR l
#            
#            #if( np.sum(kl) != 0):
#            Ref1[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
#        #END FOR j
#        
#        #Generic progress reporter
#        if( np.mod(i,finAnnounce_modRef) == 0 ):
#            toc = time.time() - tic; #sec, calc current toc
#            eta = toc/(finAnnounce_percentToUpdateAt*i/finAnnounce_modRef)*100-toc; #sec, estimate time till finished
#            print('PROGRESS: Averaging TEC\t %Complete: '+str(np.round(finAnnounce_percentToUpdateAt*i/finAnnounce_modRef,1))+'%\tRuntime: '+str(np.round(toc,2))+' sec/'+str(np.round(toc/60,2))+' min\tETA: '+str(np.round(eta,2))+' sec/'+str(np.round(eta/60,2))+' min'); #announce job, % complete, and time so far
#        #END IF 
#    #END FOR i
    
    toc = time.time() - tic; #sec, calc current toc
    print("Time to run: "+str(np.round(toc,2))+" sec/ "+str(np.round(toc/60,2))+" min");
    
    return MagCAN_keo, MagCAN_timeUnique, (MagCAN_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600, keo_angle, keo_width, keo_range_chunks_long_plot, keo_range_chunks_long_plot_name