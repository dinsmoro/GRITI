#GOAL: Run any angle averaging strips on TEC data in requested area, also plot visualization of averaging area
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: vTECChunked_anyAngleAvg, avg_anyAngle, avg_anyAngle_Width,  avg_anyAngle_Range_Chunks_Long_Plot, avg_anyAngle_Range_Chunks_Long_Plot_Name

import numpy as np #import in here I dunno
import os
import time
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.ticker import FormatStrFormatter

from Code.GRITI_TEC_keo_plot_area import GRITI_TEC_keo_plot_area
from Code.GRITI_TEC_keo_fancyPlot_area import GRITI_TEC_keo_fancyPlot_area
from Code.GRITI_keo_subfun_raytrace import GRITI_keo_subfun_raytrace_cosApprox

#this version exists to use noise instead of TEC_float, as modifying TEC_float within a function would use a lot of RAM and I don't ahve enough for that
def GRITI_TEC_keo_noise(plotLatRange,plotLongRange,TEC_noise,TEC_float,TEC_timeUnique,TEC_plotLimValu,
        colorMap,locFloat_dTEC,locFloat_time,locFloat_lat,locFloat_long, time_Ref, avg_anyAngle,avg_anyAngle_N,avg_anyAngle_Width,
        avg_anyAngle_45vsLatLong,avgPt_coords,geoMap_projectionStyle,dateRange_dayNum_zeroHr,plotLatRange_autoTick,
        plotLongRange_autoTick,plotLongRange_autoTick_Crunched, gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size,
        FONT_titleFM,FONT_axisTick,FONT_axisTickFM,FONT_axisLabelFM,BasemapFixDir,
        avg_anyAngle_dataType,avg_anyAngle_plotLabel,FLG_fancyPlot,PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi,
        avg_anyAngle_polarMode=0,FLG_disablePlot=0):
    
    #==============Analysis: Any Angle AVG (Keogram)==============
    print('Analyis: Any Angle AVG (Keogram) Beginning:'); 
    tic = time.time();
    
    #AVG any angle also can't do exactly the full width (so a width that takes the whole plot area)
    if( ( (np.round(avg_anyAngle) == 0) | (np.round(avg_anyAngle) == 180) ) & (avg_anyAngle_Width >= np.abs(plotLatRange[0] - plotLatRange[1])) ):
        if( avg_anyAngle_Width > np.abs(plotLatRange[0] - plotLatRange[1]) ):
            print("\n==============~Warning~==============");
            print('Avg any angle width of '+str(avg_anyAngle_Width)+' arcdeg is larger than the latitudinal plot area of '+str(np.abs(plotLatRange[1] - plotLatRange[0]))+' arcdeg. Reducing to be size of latitudinal plot area.');
            avg_anyAngle_Width = np.abs(plotLatRange[0] - plotLatRange[1]) - 0.001; #arcdeg, adjust so not exactly plot width, but basically since overstepped the range already
        else:
            avg_anyAngle_Width = avg_anyAngle_Width - 0.001; #arcdeg, adjust so not exactly plot width
        #END IF
    elif( ( (np.round(avg_anyAngle) == 90) | (np.round(avg_anyAngle) == 270) ) & (avg_anyAngle_Width >= np.abs(plotLongRange[0] - plotLongRange[1])) ):
        if( avg_anyAngle_Width > np.abs(plotLongRange[0] - plotLongRange[1]) ):
            print("\n==============~Warning~==============");
            print('Avg any angle width of '+str(avg_anyAngle_Width)+' arcdeg is larger than the longitudinal plot area of '+str(np.abs(plotLongRange[0] - plotLongRange[1]))+' arcdeg. Reducing to be size of longitudinal plot area.');
            avg_anyAngle_Width = np.abs(plotLongRange[0] - plotLongRange[1]) - 0.001; #arcdeg, adjust so not exactly plot width, but just under it since overstepped the range already
        else:
            avg_anyAngle_Width = avg_anyAngle_Width - 0.001; #arcdeg, adjust so not exactly plot width
        #END IF
    #END IF
        
    #AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
    if( (avg_anyAngle == 0) | (avg_anyAngle == 270) ):
        avg_anyAngle = avg_anyAngle + 0.0001; #deg, adjust so not exactly 0 or 270 (to be symmetrically consistent)
    elif( (avg_anyAngle == 90) | (avg_anyAngle == 180) ):
        avg_anyAngle = avg_anyAngle - 0.0001; #deg, adjust so not exactly 90 or 180 (to be symmetrically consistent)
    #END IF    
    
    avg_anyAngle = avg_anyAngle*np.pi/180; #rad, convert to radians
    
    avg_anyAngle_slope = np.tan(avg_anyAngle); #get slope of line required
    #this conversion is for y=LATITUDE x=LONGITUDE line action
    
    #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) 90 deg (real angle) to the req
    #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
    avg_anyAngle_upLine_int = avg_anyAngle_Width/2*np.sin(avg_anyAngle + np.pi/2) + np.mean(plotLatRange)  \
        - avg_anyAngle_slope*(avg_anyAngle_Width/2*np.cos(avg_anyAngle + np.pi/2) + np.mean(plotLongRange)); #get intercept of upper line
    #upper and lower lines are parallel
    #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) -90 deg (real angle) to the req
    #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
    avg_anyAngle_loLine_int = avg_anyAngle_Width/2*np.sin(avg_anyAngle - np.pi/2) + np.mean(plotLatRange)  \
        - avg_anyAngle_slope*(avg_anyAngle_Width/2*np.cos(avg_anyAngle - np.pi/2) + np.mean(plotLongRange)); #get intercept of lower line
    
    avg_anyAngle_LatLim_upLine = avg_anyAngle_slope*np.sort(plotLongRange) + avg_anyAngle_upLine_int; #arcdeg,
    #latitude range from largest possible longitude range
    avg_anyAngle_LatLim_upLine[ avg_anyAngle_LatLim_upLine > np.max(plotLatRange)] = np.max(plotLatRange); #arcdeg, limit lat to current range
    avg_anyAngle_LatLim_upLine[ avg_anyAngle_LatLim_upLine < np.min(plotLatRange)] = np.min(plotLatRange); #arcdeg, limit lat to current range
    
    #Project upper line pts to lower pts
    avg_anyAngle_LatLim_loLine = np.array( [np.min(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*np.sin(avg_anyAngle - np.pi/2) , np.max(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*np.sin(avg_anyAngle - np.pi/2)] ); #arcdeg, calc lower pts based on upper
    avg_anyAngle_LatLim_loLine[ avg_anyAngle_LatLim_loLine > np.max(plotLatRange)] = np.max(plotLatRange); #arcdeg, limit lat to current range
    avg_anyAngle_LatLim_loLine[ avg_anyAngle_LatLim_loLine < np.min(plotLatRange)] = np.min(plotLatRange); #arcdeg, limit lat to current range
    #redefine upper pts based on possibly adjusted lower pts
    avg_anyAngle_LatLim_upLine = np.array( [np.min(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*np.sin(avg_anyAngle + np.pi/2) , np.max(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*np.sin(avg_anyAngle + np.pi/2)] ); #arcdeg, calc upper pts based on lower
    avg_anyAngle_LongLim_upLine = (avg_anyAngle_LatLim_upLine - avg_anyAngle_upLine_int)/avg_anyAngle_slope; #arcdeg, get longitudes that match
    avg_anyAngle_LongLim_loLine = (avg_anyAngle_LatLim_loLine - avg_anyAngle_loLine_int)/avg_anyAngle_slope; #arcdeg, get longitudes that match
    
    avg_anyAngle_Range = np.concatenate( [np.concatenate(  [avg_anyAngle_LatLim_upLine.T , avg_anyAngle_LatLim_loLine.T]  ).reshape(-1,1) , np.concatenate( [avg_anyAngle_LongLim_upLine , avg_anyAngle_LongLim_loLine] ).reshape(-1,1)] , axis = 1); #arcdeg, record pts for use
    
    avg_anyAngle_Range[np.where(avg_anyAngle_Range[:,1] > 180),1]  = 180.; #fix over 180 so it flips over to the other side
    avg_anyAngle_Range[np.where(avg_anyAngle_Range[:,1] < -180),1]  = -180.; #fix less than -180 so it flips over to the other side
    avg_anyAngle_Range[np.where(avg_anyAngle_Range[:,0] > 90),0]  = 90.; #fix over 90 so it flips over to the other side
    avg_anyAngle_Range[np.where(avg_anyAngle_Range[:,0] < -90),0]  = -90.; #fix less than -90 so it flips over to the other side
    
    avg_anyAngle_Range_Chunks_Long_up = np.linspace( np.min(avg_anyAngle_Range[0:2,1]) , np.max(avg_anyAngle_Range[0:2,1]) , avg_anyAngle_N + 1); #arcdeg, chunks in longitude to go between
    #chose longitude because 0 deg will be stable - 90 deg would never be with
    #my hella math *wasn't stable at 0 anyway lol*
    avg_anyAngle_Range_Chunks_Long_lo = np.linspace( np.min(avg_anyAngle_Range[2:4,1]) , np.max(avg_anyAngle_Range[2:4,1]) , avg_anyAngle_N + 1); #arcdeg, chunks in longitude to go between
    # avg_anyAngle_Range_Chunks_Lat = linspace( min(avg_anyAngle_Range(1:2,1)) , max(avg_anyAngle_Range(1:2,1)) , avg_anyAngle_N + 1)'; %arcdeg, chunks in latitude to go between
    
    avg_anyAngle_Range_Chunks_Long_Plot = np.linspace( (avg_anyAngle_Range_Chunks_Long_up[0] + avg_anyAngle_Range_Chunks_Long_lo[0])/2 ,\
        (avg_anyAngle_Range_Chunks_Long_up[-1] + avg_anyAngle_Range_Chunks_Long_lo[-1])/2 , avg_anyAngle_N); #for plotting, the mean between each point
    #for plotting using up
    
    if( (avg_anyAngle_45vsLatLong == 1) & ((avg_anyAngle == 45*np.pi/180) | (avg_anyAngle == 135*np.pi/180)) ): #override mechanism to allow LATITUDE when normally longitude is defaulted to
        
        avg_anyAngle_Range_Chunks_Long_Plot_Name = 'Latitude'; #name of Latitude
        
        avg_anyAngle_Range_Chunks_Long_Plot_int = np.mean(plotLatRange)  \
            - avg_anyAngle_slope*np.mean(plotLongRange); #get intercept of upper line
        
        avg_anyAngle_Range_Chunks_Long_Plot = avg_anyAngle_slope*avg_anyAngle_Range_Chunks_Long_Plot + avg_anyAngle_Range_Chunks_Long_Plot_int;
        #for plotting rocking it up
        
    elif( (avg_anyAngle <= 45*np.pi/180) | (avg_anyAngle >= 135*np.pi/180) ): #actually LONGITUDE on the axis
    
        avg_anyAngle_Range_Chunks_Long_Plot_Name = 'Longitude'; #name of longitude
    
    else: #otherwise LATITUDE on the axis
        avg_anyAngle_Range_Chunks_Long_Plot_Name = 'Latitude'; #name of Latitude
        
        avg_anyAngle_Range_Chunks_Long_Plot_int = np.mean(plotLatRange)  \
            - avg_anyAngle_slope*np.mean(plotLongRange); #get intercept of upper line
        
        avg_anyAngle_Range_Chunks_Long_Plot = avg_anyAngle_slope*avg_anyAngle_Range_Chunks_Long_Plot + avg_anyAngle_Range_Chunks_Long_Plot_int;
        #for plotting rocking it up
    #END IF
    #avg_anyAngle_Range_Chunks_Long_Plot = np.round(avg_anyAngle_Range_Chunks_Long_Plot,decimals=5); #keeps things from getting out of hand
    
    temp_Longs_up = np.zeros( [avg_anyAngle_N,2] ); #preallocate
    temp_Lats_up = np.zeros( [avg_anyAngle_N,2] ); #preallocate
    temp_Longs_lo = np.zeros( [avg_anyAngle_N,2] ); #preallocate
    temp_Lats_lo = np.zeros( [avg_anyAngle_N,2] ); #preallocate
    for j in range(0,avg_anyAngle_N): #preallocate and fill
        temp_Longs_up[j,:] = [avg_anyAngle_Range_Chunks_Long_up[j] , avg_anyAngle_Range_Chunks_Long_up[j+1]]; #arcdeg, get longitudes needed upper line
        temp_Longs_lo[j,:] = np.flip( np.array( [avg_anyAngle_Range_Chunks_Long_lo[j] , avg_anyAngle_Range_Chunks_Long_lo[j+1]] ) ); #arcdeg, get longitudes needed low
        temp_Lats_up[j,:] = avg_anyAngle_slope*temp_Longs_up[j,:] + avg_anyAngle_upLine_int; #arcdeg, get latitudes needed up
        temp_Lats_lo[j,:] = avg_anyAngle_slope*temp_Longs_lo[j,:] + avg_anyAngle_loLine_int; #arcdeg, get latitudes needed lower line
    #END IF
    
    temp_Long_List = np.zeros( [avg_anyAngle_N,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
    temp_Lat_List = np.zeros( [avg_anyAngle_N,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
    for j in range(0,avg_anyAngle_N):
        temp_Long_List[j,:] = np.hstack( [temp_Longs_up[j,:],temp_Longs_lo[j,:],temp_Longs_up[j,0]] ); #so this isn't done dynamtically
        temp_Lat_List[j,:] = np.hstack( [temp_Lats_up[j,:],temp_Lats_lo[j,:],temp_Lats_up[j,0]] );
    #END IF
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(TEC_plotLimValu) == 1 ):
        TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    #END IF
    
    if( FLG_disablePlot == 0 ):
        #THIS IS TO VIEW DATA BEING AVERAGED
        GRITI_TEC_keo_plot_area(TEC_float,locFloat_time,locFloat_dTEC,\
            locFloat_lat,locFloat_long,TEC_timeUnique, TEC_plotLimValu, time_Ref, \
            temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, avgPt_coords, avg_anyAngle_Range, \
            avg_anyAngle, avg_anyAngle_Width, avg_anyAngle_N, avg_anyAngle_dataType, avg_anyAngle_plotLabel, \
            plotLatRange,plotLongRange,plotLatRange_autoTick,plotLongRange_autoTick, \
            dateRange_dayNum_zeroHr, gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size, \
            colorMap, avg_anyAngle_polarMode,geoMap_projectionStyle,BasemapFixDir, \
            FONT_titleFM, FONT_axisLabelFM, FONT_axisTickFM, FONT_axisTick); #plot the area
    if( FLG_fancyPlot == 1 ):
        #THIS IS TO VIEW DATA BEING AVERAGED ~fancily~
        GRITI_TEC_keo_fancyPlot_area(TEC_float,locFloat_time,locFloat_dTEC,\
            locFloat_lat,locFloat_long,TEC_timeUnique, TEC_plotLimValu, time_Ref, \
            temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, avgPt_coords, avg_anyAngle_Range, \
            avg_anyAngle, avg_anyAngle_Width, avg_anyAngle_N, avg_anyAngle_dataType, avg_anyAngle_plotLabel, \
            plotLatRange,plotLongRange,plotLatRange_autoTick,plotLongRange_autoTick_Crunched, \
            dateRange_dayNum_zeroHr, gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size, \
            colorMap, avg_anyAngle_polarMode,geoMap_projectionStyle,BasemapFixDir, \
            FONT_titleFM, FONT_axisLabelFM, FONT_axisTickFM, FONT_axisTick, \
            PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi); #plot the area
    #END IF

    #limit the work needed to calc the stuff by making sure only using data inside the grid
    sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
    sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
    P1P2Vxy = np.array([avg_anyAngle_Range[1,1] - avg_anyAngle_Range[0,1] , avg_anyAngle_Range[1,0] - avg_anyAngle_Range[0,0]]); #effort to do own in polygon
    P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1]));
    #P1P2Vm = np.sqrt( (P1P2Vxy[0])**2 + (P1P2Vxy[1])**2 ); #no approx for magnitude used
    P1P4Vxy = np.array([avg_anyAngle_Range[2,1] - avg_anyAngle_Range[0,1] , avg_anyAngle_Range[2,0] - avg_anyAngle_Range[0,0]]); #effort to do own in polygon
    P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1]));
    #P1P4Vm = np.sqrt( (P1P4Vxy[0])**2 + (P1P4Vxy[1])**2 ); #no approx for magnitude used
    P3P2Vxy = np.array([avg_anyAngle_Range[1,1] - avg_anyAngle_Range[3,1] , avg_anyAngle_Range[1,0] - avg_anyAngle_Range[3,0]]); #effort to do own in polygon
    P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1]));
    #P3P2Vm = np.sqrt( (P3P2Vxy[0])**2 + (P3P2Vxy[1])**2 ); #no approx for magnitude used
    P3P4Vxy = np.array([avg_anyAngle_Range[2,1] - avg_anyAngle_Range[3,1] , avg_anyAngle_Range[2,0] - avg_anyAngle_Range[3,0]]); #effort to do own in polygon
    P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1]));
    #P3P4Vm = np.sqrt( (P3P4Vxy[0])**2 + (P3P4Vxy[1])**2 ); #no approx for magnitude used
    P1P2P1P4dotCosT = np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
    P3P2P3P4dotCosT = np.sum(P3P2Vxy.conj()*P3P4Vxy,axis=0)/(P3P2Vm*P3P4Vm);
    
    P1PtVxy = np.array([TEC_float[:,locFloat_long] - avg_anyAngle_Range[0,1] , TEC_float[:,locFloat_lat] - avg_anyAngle_Range[0,0]]); #effort to do own in polygon
    P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
    #P1PtVm = np.sqrt( (P1PtVxy[0,:])**2 + (P1PtVxy[1,:])**2 ); #no approx for magnitude used
    P3PtVxy = np.array([TEC_float[:,locFloat_long] - avg_anyAngle_Range[3,1] , TEC_float[:,locFloat_lat] - avg_anyAngle_Range[3,0]]); #effort to do own in polygon
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
            
    tic = time.time(); #get current time (it's like from 1970 or whatever time)
    
    #keepr = pplong < max(max(temp_Long_List)) &  pplong > min(min(temp_Long_List)) & pplat < max(max(temp_Lat_List)) & pplat > min(min(temp_Lat_List)); %limit memory sent to parallel workers
    #keepr = inpolygon(pplong,pplat, [avg_anyAngle_Range(1:2,2) ; flipud(avg_anyAngle_Range(3:4,2)) ; avg_anyAngle_Range(1,2)] , [avg_anyAngle_Range(1:2,1) ; flipud(avg_anyAngle_Range(3:4,1)) ; avg_anyAngle_Range(1,1)]);
    time_limd =  TEC_float[keepr,locFloat_time];
    vTEC_limd =  TEC_noise[keepr].T;
    pplat_limd = TEC_float[keepr,locFloat_lat];
    pplong_limd = TEC_float[keepr,locFloat_long];
    
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
      
    vTECChunked_anyAngleAvg = GRITI_keo_subfun_raytrace(time_limd,TEC_timeUnique,pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,vTEC_limd,avg_anyAngle_N,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT); #call a function to avg all the TEC together into bands
    
#    #original, in-fun way (slower w/o jit etc.)
#        #pre-calc the inpolygon path stuff
#    inpolgygon = []; #prep a list
#    for j in range( 0, avg_anyAngle_N):
#        inpolgygon.append(Path([(temp_Long_List[j,0],temp_Lat_List[j,1]), (temp_Long_List[j,0], temp_Lat_List[j,2]), \
#                (temp_Long_List[j,1], temp_Lat_List[j,2]), (temp_Long_List[j,1], temp_Lat_List[j,2])]));  # square with legs length 1 and bottom left corner at the origin
#    #END FOR j
#    
#    finAnnounce_percentToUpdateAt = 0.5; # every % to update info at
#    finAnnounce_div = np.round(100/finAnnounce_percentToUpdateAt); #calc divisor to use
#    finAnnounce_modRef = np.round(len(TEC_timeUnique)/finAnnounce_div); #calc mod reference to be used
#    Ref1 = np.zeros( [len(TEC_timeUnique),avg_anyAngle_N] ,dtype=np.float32 ); #preallocate
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
#        for j in range(0,avg_anyAngle_N): #avg_anyAngle_N
#            #average vTEC for a range chunk on an angle
#            #kl = np.where(inpolgygon[j].contains_points( temp_pplatpplong ))[0];  #get pts inside the area defined
#            #MATPLOTLIB'S INPOLYGON IS INCORRECT!
#            #Gets pts inside the range
#            
#            P1PtVxy = np.array([temp_pplong - temp_Long_List[j,0] , temp_pplat - temp_Lat_List[j,0]]); #effort to do own in polygon
#            P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
#            P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat - temp_Lat_List[j,2]]); #effort to do own in polygon
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
    
    return vTECChunked_anyAngleAvg, avg_anyAngle, avg_anyAngle_Width,  avg_anyAngle_Range_Chunks_Long_Plot, avg_anyAngle_Range_Chunks_Long_Plot_Name