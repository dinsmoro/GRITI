#GOAL: Run any angle averaging strips on TEC data in requested area, also plot visualization of averaging area
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: keo, keo_angle, keo_width,  keo_range_plotLatLong_chunks, keo_range_plotLatLong_name

import numpy as np #import in here I dunno
import os
import time
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.ticker import FormatStrFormatter

from Code.GRITI_keo_plot_area import GRITI_keo_plot_area
from Code.GRITI_keo_plot_dataDensity import GRITI_keo_plot_dataDensity
# from Code.GRITI_keo_plot_dataReceiverLocation import GRITI_keo_plot_dataReceiverLocation
from Code.GRITI_keo_subfun_raytrace import GRITI_keo_subfun_raytrace_cosApprox, GRITI_keo_subfun_raytrace_cosExact, GRITI_keo_subfun_raytrace_tan
from Code.GRITI_keo_subfun_simpleTrace_lat import GRITI_keo_subfun_simpleTrace_lat
from Code.GRITI_keo_subfun_simpleTrace_long import GRITI_keo_subfun_simpleTrace_long
from Code.subfun_textNice import textNice


def GRITI_keo_keogrammer(data_data ,data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates,
        settings_dataSpecific, settings_paths, settings_map, settings_plot,
        FLG_fancyPlot=0,FLG_disablePlot=0,FLG_dataDensity=0,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True):
    #==============Analysis: Any Angle AVG (Keogram)==============
    if( FLG_disableText == 0 ):
        print('Analyis: Keogram Beginning:'); 
    #END IF
    # tic = time.time();
    #----- Unpack -----
    plotLatRange = settings_map['lat range']; #unpack
    plotLongRange = settings_map['long range']; #unpack
    plotLimValu = settings_dataSpecific['keo plot lim']; #unpack
    keo_angle = settings_dataSpecific['keo angle'];
    keo_width = settings_dataSpecific['keo width orig'];
    keo_N = settings_dataSpecific['keo N'];
    keo_45vsLatLong = settings_dataSpecific['keo 45 lat or long'];
    keo_name = settings_dataSpecific['keo labels'];
    if( FLG_disableCache == 0 ):
        keo_dataVersion = settings_dataSpecific['version']; #get the data version to prevent cached old data versions from being used
    #END IF
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(plotLimValu) == 1 ):
        plotLimValu = np.array( (-plotLimValu,plotLimValu) ); #make it a vector
    #END IF
    
    # #AVG any angle also can't do exactly the full width (so a width that takes the whole plot area)
    # if( ( (np.round(keo_angle) == 0) | (np.round(keo_angle) == 180) ) & (keo_width >= np.abs(plotLatRange[0] - plotLatRange[1])) ):
    #     if( keo_width > np.abs(plotLatRange[0] - plotLatRange[1]) ):
    #         print("\n==============~Warning~==============");
    #         print('Avg any angle width of '+str(keo_width)+' arcdeg is larger than the latitudinal plot area of '+str(np.abs(plotLatRange[1] - plotLatRange[0]))+' arcdeg. Reducing to be size of latitudinal plot area.');
    #         keo_width = np.abs(plotLatRange[0] - plotLatRange[1])- 0.001; #arcdeg, adjust so not exactly plot width, but basically since overstepped the range already
    #     else:
    #         keo_width = keo_width- 0.001; #arcdeg, adjust so not exactly plot width
    #     #END IF
    # elif( ( (np.round(keo_angle) == 90) | (np.round(keo_angle) == 270) ) & (keo_width >= np.abs(plotLongRange[0] - plotLongRange[1])) ):
    #     if( keo_width > np.abs(plotLongRange[0] - plotLongRange[1]) ):
    #         print("\n==============~Warning~==============");
    #         print('Avg any angle width of '+str(keo_width)+' arcdeg is larger than the longitudinal plot area of '+str(np.abs(plotLongRange[0] - plotLongRange[1]))+' arcdeg. Reducing to be size of longitudinal plot area.');
    #         keo_width = np.abs(plotLongRange[0] - plotLongRange[1])- 0.001; #arcdeg, adjust so not exactly plot width, but just under it since overstepped the range already
    #     else:
    #         keo_width = keo_width- 0.001; #arcdeg, adjust so not exactly plot width
    #     #END IF
    # #END IF
    # 
    # #AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
    # if( (keo_angle == 0) | (keo_angle == 270) ):
    #     keo_angle = keo_angle + 0.0001; #deg, adjust so not exactly 0 or 270 (to be symmetrically consistent)
    # elif( (keo_angle == 90) | (keo_angle == 180) ):
    #     keo_angle = keo_angle - 0.0001; #deg, adjust so not exactly 90 or 180 (to be symmetrically consistent)
    # #END IF
        
    if( (FLG_useRightExact == 1) & ((keo_angle == 0) | (keo_angle == 270) | (keo_angle == 90) | (keo_angle == 180)) ): #don't use this on the reg, ray trace is x3 faster
        #raytracing is actually faster somehow
        #for angles on the axes, raytracing isn't needed
        
        #clip keo_width here because it's based on the plot area directly
        if( ( (np.round(keo_angle) == 0) | (np.round(keo_angle) == 180) ) & (keo_width > np.abs(plotLatRange[0] - plotLatRange[1])) ):
            if( FLG_disableText == 0 ):
                print("\n==============~Warning~==============");
                print('Avg any angle width of '+str(keo_width)+' arcdeg is larger than the latitudinal plot area of '+str(np.abs(plotLatRange[1] - plotLatRange[0]))+' arcdeg. Reducing to be size of latitudinal plot area.');
            #END IF
            keo_width = np.abs(plotLatRange[0] - plotLatRange[1]); #arcdeg, adjust so not exactly plot width, but basically since overstepped the range already
        elif( ( (np.round(keo_angle) == 90) | (np.round(keo_angle) == 270) ) & (keo_width > np.abs(plotLongRange[0] - plotLongRange[1])) ):
            if( FLG_disableText == 0 ):
                print("\n==============~Warning~==============");
                print('Avg any angle width of '+str(keo_width)+' arcdeg is larger than the longitudinal plot area of '+str(np.abs(plotLongRange[0] - plotLongRange[1]))+' arcdeg. Reducing to be size of longitudinal plot area.');
            #END IF
            keo_width = np.abs(plotLongRange[0] - plotLongRange[1]); #arcdeg, adjust so not exactly plot width, but just under it since overstepped the range already
        #END IF
        
        settings_dataSpecific['keo width'] = keo_width; #record new keo width
                        
        if( (keo_angle == 0) | (keo_angle == 180) ):
            latLims = np.array( (np.mean(plotLatRange)-keo_width/2, np.mean(plotLatRange)+keo_width/2) ); #degc, lat range is based on the keo_width
            longLims = np.array( (np.min(plotLongRange), np.max(plotLongRange)) ); #degc, long range is the plot range
            
            keepr = (data_long <= np.max(longLims)) &  (data_long >= np.min(longLims)) & (data_lat <= np.max(latLims)) & (data_lat >= np.min(latLims)); #limit work done
            if( np.all(keepr) == 0 ):
                time_limd =  data_time[keepr];
                vTEC_limd =  data_data[keepr];
                pplat_limd = data_lat[keepr];
                pplong_limd = data_long[keepr];
            else:
                time_limd =  data_time;
                vTEC_limd =  data_data;
                pplat_limd = data_lat;
                pplong_limd = data_long;
            #END IF
            
            #Check for cached version
            if( FLG_disableCache == 0 ):
                import pickle
                cachedName = keo_name+str(dates['date range dayNum']).replace('\n','')+str(plotLatRange)+\
                            str(plotLongRange)+textNice(np.round(keo_angle,2))+'_'+textNice(np.round(keo_width,2))+'_'+\
                            str(keo_N)+'_'+str(keo_45vsLatLong)+'_'+str(pplat_limd.size)+\
                            str(settings_dataSpecific['source to use'])+str(keo_dataVersion)+settings_map['coord type']+'_'+\
                            'rightexact-'+str(FLG_useRightExact); #make a big string to uniquely identify this keogram
                if( os.path.isfile(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl')) == 1 ):
                    #if the data already exists, load it in
                    with open(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl'), 'rb') as keo_cacher:
                        keo = pickle.load(keo_cacher); #load in the data
                    #END WITH
                    FLG_cachedAvail = 1; #set the flag to cached data available
                else:
                    FLG_cachedAvail = 0; #set the flag to no cached data avail
                    print('---CATCHED NOT FOUND!---');
                #END IF
            else:
                FLG_cachedAvail = 0; #set the flag to no cached data avail
            #END IF
            
            splits = np.linspace( np.min(longLims) , np.max(longLims) , keo_N + 1); #arcdeg, chunks in latitude to go between
            
            #prep making variables needed for the plotting
            keo_range = np.concatenate( [np.concatenate(  [latLims.T , latLims.T]  ).reshape(-1,1) , np.sort(np.concatenate( [longLims , longLims] ).reshape(-1,1),axis=0)] , axis = 1); #arcdeg, record pts for use
            
            temp_Lats_up = np.tile( np.min(latLims), (keo_N,2) ); #replicate other code output here
            temp_Longs_up = np.concatenate( (splits[0:-1].reshape(-1,1),splits[1:].reshape(-1,1)) , axis = 1 ); #make a sandwich of pt-to-pt
            
            temp_Lat_List = np.concatenate( (np.tile(np.min(latLims), (keo_N,1)) , np.tile(latLims, (keo_N,1)) , np.tile(np.flip(latLims), (keo_N,1))) , axis = 1); #make a sandwich of pt-to-pt
            temp_Long_List = np.concatenate( (splits[0:-1].reshape(-1,1),splits[1:].reshape(-1,1),splits[1:].reshape(-1,1),splits[0:-1].reshape(-1,1),splits[0:-1].reshape(-1,1)) , axis = 1 ); #make a sandwich of pt-to-pt
            
            keo_range_plotLatLong_chunks = splits; #for plotting, the mean between each point?
            keo_range_plotLatLong_name = 'Longitude'; #name of Longitude
            
            #plot here just so it doesn't end up in the timer
            if( FLG_disablePlot == 0 ):
                #THIS IS TO VIEW DATA BEING AVERAGED
                GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                    settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                    temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                    FLG_fancyPlot = 0); #plot the area
                if( FLG_dataDensity == 1 ):
                    GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                        FLG_fancyPlot = 0); #plot the area with data density
                #END IF
            #END IF
            if( (FLG_fancyPlot == 1) & (FLG_disablePlot != 2) ):
                #THIS IS TO VIEW DATA BEING AVERAGED ~fancily~
                GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                    settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                    temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                    FLG_fancyPlot = 1); #plot the area
                if( FLG_dataDensity == 1 ):
                    GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                        FLG_fancyPlot = 1); #plot the area with data density
                #END IF
            #END IF
            
            tic = time.time(); #get current time (it's like from 1970 or whatever time)
            if( FLG_cachedAvail == 0 ):     
                keepr = (data_lat <= np.max(latLims)) &  (data_lat >= np.min(latLims)) & (data_lat <= np.max(latLims)) & (data_lat >= np.min(latLims)); #limit work done
                # if( np.all(keepr) == 0 ):
                #     time_limd =  data_time[keepr];
                #     vTEC_limd =  data_data[keepr];
                #     pplat_limd = data_lat[keepr];
                #     pplong_limd = data_long[keepr];
                # else:
                #     time_limd =  data_time;
                #     vTEC_limd =  data_data;
                #     pplat_limd = data_lat;
                #     pplong_limd = data_long;
                # #END IF
                if( np.all(keepr) == 0 ):
                    keo = GRITI_keo_subfun_simpleTrace_long(data_time[keepr],data_timeUnique,data_lat[keepr],data_long[keepr],splits,latLims,data_data[keepr],keo_N, FLG_memSafe=0); #do the below alg faster
                else:
                    keo = GRITI_keo_subfun_simpleTrace_long(data_time,data_timeUnique,data_lat,data_long,splits,latLims,data_data,keo_N, FLG_memSafe=0); #do the below alg faster
                #END IF
                
                #faster jit function here
                # keo = GRITI_keo_subfun_simpleTrace_long(time_limd,data_timeUnique,pplat_limd,pplong_limd,splits,latLims,vTEC_limd,keo_N, FLG_memSafe=0); #do the below alg faster
                # #time to crunch keo
                # keo = np.zeros( [len(data_timeUnique),keo_N] ,dtype=np.float32 ); #preallocate
                # for i in range(0, len(data_timeUnique) ): #87
                #     #Corral the data to the right place  
                #     k = np.where(time_limd == data_timeUnique[i]); #gets during a time period
                    
                #     temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
                #     temp_pplat = pplat_limd[k];
                #     temp_pplong = pplong_limd[k];
                    
                #     for j in range(0,keo_N): #keo_N
                #         kl = (temp_pplong <= np.max(longLims)) &  (temp_pplong >= np.min(longLims)) & (temp_pplat <= splits[j+1]) & (temp_pplat >= splits[j]); #get data in the averaging zone
                #         keo[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
                #     #END FOR j
                # #END FOR i
                
                #save the data
                if( FLG_disableCache == 0 ):
                    with open(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl'), 'wb') as keo_cacher:
                        pickle.dump(keo, keo_cacher); #save that data
                    #END WITH
                #END IF
            #END IF
        else:
            #90 or 270 then
            latLims = np.array( (np.min(plotLatRange), np.max(plotLatRange)) ); #degc, lat range is the plot range
            longLims = np.array( (np.mean(plotLongRange)-keo_width/2, np.mean(plotLongRange)+keo_width/2) ); #degc, long range is based on the keo_width
            
            keepr = (data_long <= np.max(longLims)) &  (data_long >= np.min(longLims)) & (data_lat <= np.max(latLims)) & (data_lat >= np.min(latLims)); #limit work done
            if( np.all(keepr) == 0 ):
                time_limd =  data_time[keepr];
                vTEC_limd =  data_data[keepr];
                pplat_limd = data_lat[keepr];
                pplong_limd = data_long[keepr];
            else:
                time_limd =  data_time;
                vTEC_limd =  data_data;
                pplat_limd = data_lat;
                pplong_limd = data_long;
            #END IF
            
            #Check for cached version
            if( FLG_disableCache == 0 ):
                import pickle
                cachedName = keo_name+str(dates['date range dayNum']).replace('\n','')+str(plotLatRange)+\
                            str(plotLongRange)+textNice(np.round(keo_angle,2))+'_'+textNice(np.round(keo_width,2))+'_'+\
                            str(keo_N)+'_'+str(keo_45vsLatLong)+'_'+str(pplat_limd.size)+\
                            str(settings_dataSpecific['source to use'])+str(keo_dataVersion)+settings_map['coord type']+'_'+\
                            'rightexact-'+str(FLG_useRightExact); #make a big string to uniquely identify this keogram
                if( os.path.isfile(os.path.join(settings_paths['cache'], 'keo_'+cachedName+'.pkl') ) == 1 ):
                    #if the data already exists, load it in
                    with open(os.path.join(settings_paths['cache'], 'keo_'+cachedName+'.pkl'), 'rb') as keo_cacher:
                        keo = pickle.load(keo_cacher); #load in the data
                    #END WITH
                    FLG_cachedAvail = 1; #set the flag to cached data available
                else:
                    FLG_cachedAvail = 0; #set the flag to no cached data avail
                    print('---CATCHED NOT FOUND!---');
                #END IF
            else:
                FLG_cachedAvail = 0; #set the flag to no cached data avail
            #END IF
            
            splits = np.linspace( np.min(latLims) , np.max(latLims) , keo_N + 1); #arcdeg, chunks in latitude to go between
            
            #prep making variables needed for the plotting
            keo_range = np.concatenate( [np.concatenate(  [latLims.T , latLims.T]  ).reshape(-1,1) , np.sort(np.concatenate( [longLims , longLims] ).reshape(-1,1),axis=0)] , axis = 1); #arcdeg, record pts for use
            
            temp_Lats_up = np.concatenate( (splits[0:-1].reshape(-1,1),splits[1:].reshape(-1,1)) , axis = 1 ); #make a sandwich of pt-to-pt
            temp_Longs_up = np.tile( np.min(longLims), (keo_N,2) ); #replicate other code output here
            
            temp_Lat_List = np.concatenate( (splits[0:-1].reshape(-1,1),splits[1:].reshape(-1,1),splits[1:].reshape(-1,1),splits[0:-1].reshape(-1,1),splits[0:-1].reshape(-1,1)) , axis = 1 ); #make a sandwich of pt-to-pt
            temp_Long_List = np.concatenate( (np.tile(np.min(longLims), (keo_N,1)) , np.tile(longLims, (keo_N,1)) , np.tile(np.flip(longLims), (keo_N,1))) , axis = 1); #make a sandwich of pt-to-pt
            
            keo_range_plotLatLong_chunks = splits; #for plotting, the mean between each point?
            keo_range_plotLatLong_name = 'Latitude'; #name of Latitude
            
            #plot here just so it doesn't end up in the timer
            if( FLG_disablePlot == 0 ):
                #THIS IS TO VIEW DATA BEING AVERAGED
                GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                    settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                    temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                    FLG_fancyPlot = 0); #plot the area
                if( FLG_dataDensity == 1 ):
                    GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                        FLG_fancyPlot = 0); #plot the area with data density
                # elif( FLG_dataDensity == 2 ):
                #     GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                #         FLG_fancyPlot = 0); #plot the area with data density
                # elif( FLG_dataDensity == 3 ):
                #     GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                #         FLG_fancyPlot = 0); #plot the area with data density
                #     GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                #         FLG_fancyPlot = 0); #plot the area with data density
                #END IF
            #END IF
            if( (FLG_fancyPlot == 1) & (FLG_disablePlot != 2) ):
                #THIS IS TO VIEW DATA BEING AVERAGED ~fancily~
                GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                    settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                    temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                    FLG_fancyPlot = 1); #plot the area
                if( FLG_dataDensity == 1 ):
                    GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                        settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                        temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                        FLG_fancyPlot = 1); #plot the area with data density
                # elif( FLG_dataDensity == 2 ):
                #     GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                #         FLG_fancyPlot = 1); #plot the area with data density
                # elif( FLG_dataDensity == 3 ):
                #     GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                #         FLG_fancyPlot = 1); #plot the area with data density
                #     GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                #         FLG_fancyPlot = 1); #plot the area with data density
                #END IF
            #END IF
            
            tic = time.time(); #get current time (it's like from 1970 or whatever time)
            #Check for a cached version
            if( FLG_cachedAvail == 0 ):
                keepr = (data_long <= np.max(longLims)) &  (data_long >= np.min(longLims)) & (data_lat <= np.max(latLims)) & (data_lat >= np.min(latLims)); #limit work done
                # if( np.all(keepr) == 0 ):
                #     time_limd =  data_time[keepr];
                #     vTEC_limd =  data_data[keepr];
                #     pplat_limd = data_lat[keepr];
                #     pplong_limd = data_long[keepr];
                # else:
                #     time_limd =  data_time;
                #     vTEC_limd =  data_data;
                #     pplat_limd = data_lat;
                #     pplong_limd = data_long;
                # #END IF
                if( np.all(keepr) == 0 ):
                    keo = GRITI_keo_subfun_simpleTrace_lat(data_time[keepr],data_timeUnique,data_lat[keepr],data_long[keepr],splits,longLims,data_data[keepr],keo_N, FLG_memSafe=0); #do the below alg faster
                else:
                    keo = GRITI_keo_subfun_simpleTrace_lat(data_time,data_timeUnique,data_lat,data_long,splits,longLims,data_data,keo_N, FLG_memSafe=0); #do the below alg faster
                #END IF
                
                #faster jit function here
                # keo = GRITI_TEC_avgAnyAngle_subfun_Simpletrace_lat(time_limd,data_timeUnique,pplat_limd,pplong_limd,splits,longLims,vTEC_limd,keo_N); #do the below alg faster
                # #time to crunch keo
                # keo = np.zeros( [len(data_timeUnique),keo_N] ,dtype=np.float32 ); #preallocate
                # for i in range(0, len(data_timeUnique) ): #87
                #     #Corral the data to the right place  
                #     k = np.where(time_limd == data_timeUnique[i]); #gets during a time period
                    
                #     temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
                #     temp_pplat = pplat_limd[k];
                #     temp_pplong = pplong_limd[k];
                    
                #     for j in range(0,keo_N): #keo_N
                #         kl = (temp_pplong <= np.max(longLims)) &  (temp_pplong >= np.min(longLims)) & (temp_pplat <= splits[j+1]) & (temp_pplat>= splits[j]); #get data in the averaging zone
                #         keo[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
                #     #END FOR j
                # #END FOR i   
                
                #save the data
                if( FLG_disableCache == 0 ):
                    with open(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl'), 'wb') as keo_cacher:
                        pickle.dump(keo, keo_cacher); #save that data
                    #END WITH
                #END IF
            #END IF
        #END IF
    else:
        #AVG any angle also can't do exactly the full width (so a width that takes the whole plot area)
        if( ( (np.round(keo_angle) == 0) | (np.round(keo_angle) == 180) ) & (keo_width >= np.abs(plotLatRange[0] - plotLatRange[1])) ):
            if( keo_width > np.abs(plotLatRange[0] - plotLatRange[1]) ):
                if( FLG_disableText == 0 ):
                    print("\n==============~Warning~==============");
                    print('Avg any angle width of '+str(keo_width)+' arcdeg is larger than the latitudinal plot area of '+str(np.abs(plotLatRange[1] - plotLatRange[0]))+' arcdeg. Reducing to be size of latitudinal plot area.');
                #END IF
                keo_width = np.abs(plotLatRange[0] - plotLatRange[1])- 0.001; #arcdeg, adjust so not exactly plot width, but basically since overstepped the range already
            else:
                keo_width = keo_width- 0.001; #arcdeg, adjust so not exactly plot width
            #END IF
        elif( ( (np.round(keo_angle) == 90) | (np.round(keo_angle) == 270) ) & (keo_width >= np.abs(plotLongRange[0] - plotLongRange[1])) ):
            if( keo_width > np.abs(plotLongRange[0] - plotLongRange[1]) ):
                if( FLG_disableText == 0 ):
                    print("\n==============~Warning~==============");
                    print('Avg any angle width of '+str(keo_width)+' arcdeg is larger than the longitudinal plot area of '+str(np.abs(plotLongRange[0] - plotLongRange[1]))+' arcdeg. Reducing to be size of longitudinal plot area.');
                #END IF
                keo_width = np.abs(plotLongRange[0] - plotLongRange[1])- 0.001; #arcdeg, adjust so not exactly plot width, but just under it since overstepped the range already
            else:
                keo_width = keo_width- 0.001; #arcdeg, adjust so not exactly plot width
            #END IF
        #END IF
        
        #AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
        if( (keo_angle == 0) | (keo_angle == 270) ):
            keo_angle = keo_angle + 0.0001; #deg, adjust so not exactly 0 or 270 (to be symmetrically consistent)
        elif( (keo_angle == 90) | (keo_angle == 180) ):
            keo_angle = keo_angle - 0.0001; #deg, adjust so not exactly 90 or 180 (to be symmetrically consistent)
        #END IF    
        
        # settings_dataSpecific['keo angle'] = keo_angle; #record new keo angle
        settings_dataSpecific['keo width'] = keo_width; #record new keo width
        
        #Check for cached version
        if( FLG_disableCache == 0 ):
            import pickle
            cachedName = keo_name+str(dates['date range dayNum']).replace('\n','')+str(plotLatRange)+\
                str(plotLongRange)+textNice(np.round(keo_angle,2))+'_'+textNice(np.round(keo_width,2))+'_'+\
                str(keo_N)+'_'+str(keo_45vsLatLong)+'_'+str(data_lat.size)+\
                str(settings_dataSpecific['source to use'])+str(keo_dataVersion)+settings_map['coord type']+'_'+\
                'raytrace-'+str(FLG_raytraceMethod)+'_'+str(FLG_raytraceExactSqrt); #make a big string to uniquely identify this keogram
            if( os.path.isfile(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl')) == 1 ):
                #if the data already exists, load it in
                with open(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl'), 'rb') as keo_cacher:
                    keo = pickle.load(keo_cacher); #load in the data
                #END WITH
                FLG_cachedAvail = 1; #set the flag to cached data available
            else:
                FLG_cachedAvail = 0; #set the flag to no cached data avail
                print('---CATCHED NOT FOUND!---');
            #END IF
        else:
            FLG_cachedAvail = 0; #set the flag to no cached data avail
        #END IF
        
        #any angle average needs ray tracing
        keo_rad = keo_angle*np.pi/180; #rad, convert to radians
        
        keo_slope = np.tan(keo_rad); #get slope of line required
        #this conversion is for y=LATITUDE x=LONGITUDE line action
        
        #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) 90 deg (real angle) to the req
        #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
        keo_upLine_int = keo_width/2*np.sin(keo_rad + np.pi/2) + np.mean(plotLatRange)  \
            - keo_slope*(keo_width/2*np.cos(keo_rad + np.pi/2) + np.mean(plotLongRange)); #get intercept of upper line
        #upper and lower lines are parallel
        #idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) -90 deg (real angle) to the req
        #angle - y=Latitude x=Longitude solving for intercept with 2 points and slope
        keo_loLine_int = keo_width/2*np.sin(keo_rad - np.pi/2) + np.mean(plotLatRange)  \
            - keo_slope*(keo_width/2*np.cos(keo_rad - np.pi/2) + np.mean(plotLongRange)); #get intercept of lower line
        
        keo_LatLim_upLine = keo_slope*np.sort(plotLongRange) + keo_upLine_int; #arcdeg,
        #latitude range from largest possible longitude range
        keo_LatLim_upLine[ keo_LatLim_upLine > np.max(plotLatRange)] = np.max(plotLatRange); #arcdeg, limit lat to current range
        keo_LatLim_upLine[ keo_LatLim_upLine < np.min(plotLatRange)] = np.min(plotLatRange); #arcdeg, limit lat to current range
        
        #Project upper line pts to lower pts
        keo_LatLim_loLine = np.array( [np.min(keo_LatLim_upLine) + keo_width*np.sin(keo_rad - np.pi/2) , np.max(keo_LatLim_upLine) + keo_width*np.sin(keo_rad - np.pi/2)] ); #arcdeg, calc lower pts based on upper
        keo_LatLim_loLine[ keo_LatLim_loLine > np.max(plotLatRange)] = np.max(plotLatRange); #arcdeg, limit lat to current range
        keo_LatLim_loLine[ keo_LatLim_loLine < np.min(plotLatRange)] = np.min(plotLatRange); #arcdeg, limit lat to current range
        #redefine upper pts based on possibly adjusted lower pts
        keo_LatLim_upLine = np.array( [np.min(keo_LatLim_loLine) + keo_width*np.sin(keo_rad + np.pi/2) , np.max(keo_LatLim_loLine) + keo_width*np.sin(keo_rad + np.pi/2)] ); #arcdeg, calc upper pts based on lower
        keo_LongLim_upLine = (keo_LatLim_upLine - keo_upLine_int)/keo_slope; #arcdeg, get longitudes that match
        keo_LongLim_loLine = (keo_LatLim_loLine - keo_loLine_int)/keo_slope; #arcdeg, get longitudes that match
        
        keo_range = np.concatenate( [np.concatenate(  [keo_LatLim_upLine.T , keo_LatLim_loLine.T]  ).reshape(-1,1) , np.concatenate( [keo_LongLim_upLine , keo_LongLim_loLine] ).reshape(-1,1)] , axis = 1); #arcdeg, record pts for use
        
        #not super sure this actually fixes it, rather just makes it go through truncation idk been a while
        keo_range[np.where(keo_range[:,1] > 180),1]  = 180.; #fix over 180 so it flips over to the other side
        keo_range[np.where(keo_range[:,1] < -180),1]  = -180.; #fix less than -180 so it flips over to the other side
        keo_range[np.where(keo_range[:,0] > 90),0]  = 90.; #fix over 90 so it flips over to the other side
        keo_range[np.where(keo_range[:,0] < -90),0]  = -90.; #fix less than -90 so it flips over to the other side
        
        keo_range_chunksLongUp = np.linspace( np.min(keo_range[0:2,1]) , np.max(keo_range[0:2,1]) , keo_N + 1); #arcdeg, chunks in longitude to go between
        #chose longitude because 0 deg will be stable - 90 deg would never be with
        #my hella math *wasn't stable at 0 anyway lol*
        keo_range_chunksLongLo = np.linspace( np.min(keo_range[2:4,1]) , np.max(keo_range[2:4,1]) , keo_N + 1); #arcdeg, chunks in longitude to go between
        # keo_Range_Chunks_Lat = linspace( min(keo_range(1:2,1)) , max(keo_range(1:2,1)) , keo_N + 1)'; %arcdeg, chunks in latitude to go between
        
        keo_range_plotLatLong_chunks = np.linspace( (keo_range_chunksLongUp[0] + keo_range_chunksLongLo[0])/2 ,\
            (keo_range_chunksLongUp[-1] + keo_range_chunksLongLo[-1])/2 , keo_N+1); #for plotting, the mean between each point
        #for plotting using up
        
        if( (keo_45vsLatLong == 1) & ((keo_angle == 45) | (keo_angle == 135)) ): #override mechanism to allow LATITUDE when normally longitude is defaulted to
            
            keo_range_plotLatLong_name = 'Latitude'; #name of Latitude
            
            keo_range_plotLatLong_int = np.mean(plotLatRange)  \
                - keo_slope*np.mean(plotLongRange); #get intercept of upper line
            
            keo_range_plotLatLong_chunks = keo_slope*keo_range_plotLatLong_chunks + keo_range_plotLatLong_int;
            #for plotting rocking it up
            
        elif( (keo_angle <= 45) | (keo_angle >= 135) ): #actually LONGITUDE on the axis
        
            keo_range_plotLatLong_name = 'Longitude'; #name of longitude
        
        else: #otherwise LATITUDE on the axis
            keo_range_plotLatLong_name = 'Latitude'; #name of Latitude
            
            keo_range_plotLatLong_int = np.mean(plotLatRange)  \
                - keo_slope*np.mean(plotLongRange); #get intercept of upper line
            
            keo_range_plotLatLong_chunks = keo_slope*keo_range_plotLatLong_chunks + keo_range_plotLatLong_int;
            #for plotting rocking it up
        #END IF
        #keo_range_plotLatLong_chunks = np.round(keo_range_plotLatLong_chunks,decimals=5); #keeps things from getting out of hand
        
        temp_Longs_up = np.zeros( [keo_N,2] ); #preallocate
        temp_Lats_up = np.zeros( [keo_N,2] ); #preallocate
        temp_Longs_lo = np.zeros( [keo_N,2] ); #preallocate
        temp_Lats_lo = np.zeros( [keo_N,2] ); #preallocate
        for j in range(0,keo_N): #preallocate and fill
            temp_Longs_up[j,:] = [keo_range_chunksLongUp[j] , keo_range_chunksLongUp[j+1]]; #arcdeg, get longitudes needed upper line
            temp_Longs_lo[j,:] = np.flip( np.array( [keo_range_chunksLongLo[j] , keo_range_chunksLongLo[j+1]] ) ); #arcdeg, get longitudes needed low
            temp_Lats_up[j,:] = keo_slope*temp_Longs_up[j,:] + keo_upLine_int; #arcdeg, get latitudes needed up
            temp_Lats_lo[j,:] = keo_slope*temp_Longs_lo[j,:] + keo_loLine_int; #arcdeg, get latitudes needed lower line
        #END IF
        
        temp_Long_List = np.zeros( [keo_N,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
        temp_Lat_List = np.zeros( [keo_N,5] ); #preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
        for j in range(0,keo_N):
            temp_Long_List[j,:] = np.hstack( [temp_Longs_up[j,:],temp_Longs_lo[j,:],temp_Longs_up[j,0]] ); #so this isn't done dynamtically
            temp_Lat_List[j,:] = np.hstack( [temp_Lats_up[j,:],temp_Lats_lo[j,:],temp_Lats_up[j,0]] );
        #END IF
        
        #plot here just so it doesn't end up in the timer
        if( FLG_disablePlot == 0 ):
            #THIS IS TO VIEW DATA BEING AVERAGED
            GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                FLG_fancyPlot = 0); #plot the area
            if( FLG_dataDensity == 1 ):
                GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                    settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                    temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                    FLG_fancyPlot = 0); #plot the area with data density
            # elif( FLG_dataDensity == 2 ):
            #     GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
            #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
            #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
            #         FLG_fancyPlot = 0); #plot the area with data density
            # elif( FLG_dataDensity == 3 ):
            #     GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
            #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
            #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
            #         FLG_fancyPlot = 0); #plot the area with data density
            #     GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
            #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
            #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
            #         FLG_fancyPlot = 0); #plot the area with data density
            #END IF
        if( (FLG_fancyPlot == 1) & (FLG_disablePlot != 2) ):
            #THIS IS TO VIEW DATA BEING AVERAGED ~fancily~
            GRITI_keo_plot_area(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                FLG_fancyPlot = 1); #plot the area
            if( FLG_dataDensity == 1 ):
                GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
                    settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
                    temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
                    FLG_fancyPlot = 1); #plot the area with data density
            # elif( FLG_dataDensity == 2 ):
            #     GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
            #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
            #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
            #         FLG_fancyPlot = 1); #plot the area with data density
            # elif( FLG_dataDensity == 3 ):
            #     GRITI_keo_plot_dataDensity(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
            #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
            #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
            #         FLG_fancyPlot = 1); #plot the area with data density
            #     GRITI_keo_plot_dataReceiverLocation(data_data, data_time, data_lat, data_long, data_timeUnique, data_timeRef, dates, \
            #         settings_dataSpecific ,settings_plot, settings_paths, settings_map, \
            #         temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, keo_range, \
            #         FLG_fancyPlot = 1); #plot the area with data density
            #END IF
            
            # from GRITI_TEC_avgAnyAngle_fancyPlot_area_explainer import GRITI_TEC_avgAnyAngle_fancyPlot_area_explainer
            # GRITI_TEC_avgAnyAngle_fancyPlot_area_explainer(TEC_float,locFloat_time,locFloat_dTEC,\
            #     locFloat_lat,locFloat_long,data_timeUnique, plotLimValu, time_Ref, \
            #     temp_Lats_up, temp_Longs_up, temp_Lat_List, temp_Long_List, avgPt_coords, keo_range, \
            #     keo_angle, keo_width, keo_N, keo_dataType, keo_plotLabel, \
            #     plotLatRange,plotLongRange,plotLatRange_autoTick,plotLongRange_autoTick_Crunched, \
            #     dateRange_dayNum_zeroHr, gif_Millstone_Marker, gif_Millstone_Marker_Color, gif_Millstone_Marker_Size, \
            #     colorMap, keo_polarMode,geoMap_projectionStyle,BasemapFixDir, \
            #     FONT_titleFM, FONT_axisLabelFM, FONT_axisTickFM, FONT_axisTick, \
            #     PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi); #for methods paper only
        #END IF
        
        tic = time.time(); #get current time (it's like from 1970 or whatever time)
        if( FLG_cachedAvail == 0 ):
            if( FLG_raytraceMethod == 'cos' ):
                #--- prep for raytracing ---
                #limit the work needed to calc the stuff by making sure only using data inside the grid
                sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
                sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
                P1P2Vxy = np.array([keo_range[1,1] - keo_range[0,1] , keo_range[1,0] - keo_range[0,0]]); #effort to do own in polygon
                P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0]),np.abs(P1P2Vxy[1]));
                # P1P2Vm = np.sqrt( (P1P2Vxy[0])**2 + (P1P2Vxy[1])**2 ); #no approx for magnitude used
                P1P4Vxy = np.array([keo_range[2,1] - keo_range[0,1] , keo_range[2,0] - keo_range[0,0]]); #effort to do own in polygon
                P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0]),np.abs(P1P4Vxy[1]));
                # P1P4Vm = np.sqrt( (P1P4Vxy[0])**2 + (P1P4Vxy[1])**2 ); #no approx for magnitude used
                P3P2Vxy = np.array([keo_range[1,1] - keo_range[3,1] , keo_range[1,0] - keo_range[3,0]]); #effort to do own in polygon
                P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0]),np.abs(P3P2Vxy[1]));
                # P3P2Vm = np.sqrt( (P3P2Vxy[0])**2 + (P3P2Vxy[1])**2 ); #no approx for magnitude used
                P3P4Vxy = np.array([keo_range[2,1] - keo_range[3,1] , keo_range[2,0] - keo_range[3,0]]); #effort to do own in polygon
                P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0]),np.abs(P3P4Vxy[1]));
                # P3P4Vm = np.sqrt( (P3P4Vxy[0])**2 + (P3P4Vxy[1])**2 ); #no approx for magnitude used
                P1P2P1P4dotCosT = np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
                P3P2P3P4dotCosT = np.sum(P3P2Vxy.conj()*P3P4Vxy,axis=0)/(P3P2Vm*P3P4Vm);
                
                P1PtVxy = np.array([data_long - keo_range[0,1] , data_lat - keo_range[0,0]]); #effort to do own in polygon
                # P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
                P1PtVm = np.sqrt( (P1PtVxy[0,:])**2 + (P1PtVxy[1,:])**2 ); #no approx for magnitude used
                P3PtVxy = np.array([data_long - keo_range[3,1] , data_lat - keo_range[3,0]]); #effort to do own in polygon
                # P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:]));
                P3PtVm = np.sqrt( (P3PtVxy[0,:])**2 + (P3PtVxy[1,:])**2 ); #no approx for magnitude used
                
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
                #keepr = inpolygon(pplong,pplat, [keo_range(1:2,2) ; flipud(keo_range(3:4,2)) ; keo_range(1,2)] , [keo_range(1:2,1) ; flipud(keo_range(3:4,1)) ; keo_range(1,1)]);
                time_limd =  data_time[keepr];
                vTEC_limd =  data_data[keepr].T;
                pplat_limd = data_lat[keepr];
                pplong_limd = data_long[keepr];
                
                #sqrtApproxAo = 2*np.cos(np.pi/8)/(1+np.cos(np.pi/8)); #https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm
                #sqrtApproxBo = 2*np.sin(np.pi/8)/(1+np.cos(np.pi/8)); #uses this alg for approx sqrt stuff, works v well for comparison that has same approx applied
                P1P2Vxy = np.array([temp_Long_List[:,1] - temp_Long_List[:,0] , temp_Lat_List[:,1] - temp_Lat_List[:,0]]); #effort to do own in polygon
                P1P4Vxy = np.array([temp_Long_List[:,3] - temp_Long_List[:,0] , temp_Lat_List[:,3] - temp_Lat_List[:,0]]); #effort to do own in polygon
                P3P2Vxy = np.array([temp_Long_List[:,1] - temp_Long_List[:,2] , temp_Lat_List[:,1] - temp_Lat_List[:,2]]); #effort to do own in polygon
                P3P4Vxy = np.array([temp_Long_List[:,3] - temp_Long_List[:,2] , temp_Lat_List[:,3] - temp_Lat_List[:,2]]); #effort to do own in polygon
                if( FLG_raytraceExactSqrt == False ):
                    P1P2Vm = sqrtApproxAo*np.maximum(np.abs(P1P2Vxy[0,:]),np.abs(P1P2Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1P2Vxy[0,:]),np.abs(P1P2Vxy[1,:]));
                    P1P4Vm = sqrtApproxAo*np.maximum(np.abs(P1P4Vxy[0,:]),np.abs(P1P4Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1P4Vxy[0,:]),np.abs(P1P4Vxy[1,:]));
                    P3P2Vm = sqrtApproxAo*np.maximum(np.abs(P3P2Vxy[0,:]),np.abs(P3P2Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3P2Vxy[0,:]),np.abs(P3P2Vxy[1,:]));
                    P3P4Vm = sqrtApproxAo*np.maximum(np.abs(P3P4Vxy[0,:]),np.abs(P3P4Vxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3P4Vxy[0,:]),np.abs(P3P4Vxy[1,:]));
                else:
                    P1P2Vm = np.sqrt( (P1P2Vxy[0,:])**2 + (P1P2Vxy[1,:])**2 ); #no approx for magnitude used
                    P1P4Vm = np.sqrt( (P1P4Vxy[0,:])**2 + (P1P4Vxy[1,:])**2 ); #no approx for magnitude used
                    P3P2Vm = np.sqrt( (P3P2Vxy[0,:])**2 + (P3P2Vxy[1,:])**2 ); #no approx for magnitude used
                    P3P4Vm = np.sqrt( (P3P4Vxy[0,:])**2 + (P3P4Vxy[1,:])**2 ); #no approx for magnitude used
                #END IF
                P1P2P1P4dotCosT = np.einsum('ij,ij->j',P1P2Vxy,P1P4Vxy)/(P1P2Vm*P1P4Vm); #np.sum(P1P2Vxy.conj()*P1P4Vxy,axis=0)/(P1P2Vm*P1P4Vm) alt
                P3P2P3P4dotCosT = np.einsum('ij,ij->j',P3P2Vxy,P3P4Vxy)/(P3P2Vm*P3P4Vm);
                
                if( FLG_raytraceExactSqrt == False ):
                    keo = GRITI_keo_subfun_raytrace_cosApprox(time_limd,data_timeUnique,pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,vTEC_limd,keo_N,sqrtApproxAo,sqrtApproxBo,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT); #call a function to avg all the TEC together into bands
                else:
                    keo = GRITI_keo_subfun_raytrace_cosExact(time_limd,data_timeUnique,pplat_limd,pplong_limd,temp_Lat_List,temp_Long_List,vTEC_limd,keo_N,P1P2Vxy,P1P4Vxy,P3P2Vxy,P3P4Vxy,P1P2Vm,P1P4Vm,P3P2Vm,P3P4Vm,P1P2P1P4dotCosT,P3P2P3P4dotCosT); #call a function to avg all the TEC together into bands
                #END IF
            elif( FLG_raytraceMethod == 'tan' ):
                #!! this mode isn't fully baked ¯\_(ツ)_/¯ roudning issues it seems, works mostly I think !!
                #--- only work on data within the outer rectangle ---
                keo_rangeRnded = np.round(keo_range,2); #round it otherwise it lose some points, not good sign I guess
                P1P2tan = np.arctan2(keo_rangeRnded[1,0] - keo_rangeRnded[0,0], keo_rangeRnded[1,1] - keo_rangeRnded[0,1]); #rad
                P1P4tan = np.arctan2(keo_rangeRnded[2,0] - keo_rangeRnded[0,0], keo_rangeRnded[2,1] - keo_rangeRnded[0,1]); #rad
                P3P2tan = np.arctan2(keo_rangeRnded[1,0] - keo_rangeRnded[3,0], -(keo_rangeRnded[1,1] - keo_rangeRnded[3,1])); #rad, P3 is oppsoite of P1
                P3P4tan = np.arctan2(keo_rangeRnded[2,0] - keo_rangeRnded[3,0], -(keo_rangeRnded[2,1] - keo_rangeRnded[3,1])); #rad, x is negative for stability, was having issues w/ 90 deg + -180 deg instead of 90 deg + 180 deg. this gets 90 deg + -0 deg which won't rely on the weird trig 180/-180 flip to math
                
                P1Pttan = np.arctan2(data_lat - keo_rangeRnded[0,0], data_long - keo_rangeRnded[0,1]); #effort to do own in polygon
                P3Pttan = np.arctan2(data_lat - keo_rangeRnded[3,0], -(data_long - keo_rangeRnded[3,1])); #effort to do own in polygon
                
                keepr = (np.max((P1P2tan,P1P4tan)) >= P1Pttan) & (np.min((P1P2tan,P1P4tan)) <= P1Pttan) & (np.max((P3P2tan,P3P4tan)) >= P3Pttan) & (np.min((P3P2tan,P3P4tan)) <= P3Pttan); #keep the data within the averaging zone
                #keep just the data within the outer rectangle
                time_limd =  data_time[keepr];
                vTEC_limd =  data_data[keepr].T;
                pplat_limd = data_lat[keepr];
                pplong_limd = data_long[keepr];
                P1Pttan_limd = P1Pttan[keepr];
                P3Pttan_limd = P3Pttan[keepr];
                
                #--- reuse vars to be the individual smaller rectangles
                P1P2tan = np.arctan2(temp_Lat_List[:,1] - temp_Lat_List[:,0], temp_Long_List[:,1] - temp_Long_List[:,0]); #effort to do own in polygon
                P1P4tan = np.arctan2(temp_Lat_List[:,3] - temp_Lat_List[:,0], temp_Long_List[:,3] - temp_Long_List[:,0]); #effort to do own in polygon
                P3P2tan = np.arctan2(temp_Lat_List[:,1] - temp_Lat_List[:,2], -(temp_Long_List[:,1] - temp_Long_List[:,2])); #effort to do own in polygon
                P3P4tan = np.arctan2(temp_Lat_List[:,3] - temp_Lat_List[:,2], -(temp_Long_List[:,3] - temp_Long_List[:,2])); #effort to do own in polygon
                P1Pmaxtan = np.max((P1P2tan,P1P4tan),axis=0);
                P1Pmintan = np.min((P1P2tan,P1P4tan),axis=0);
                P3Pmaxtan = np.max((P3P2tan,P3P4tan),axis=0);
                P3Pmintan = np.min((P3P2tan,P3P4tan),axis=0);
                
                temp_Lat_List_min = np.min(temp_Lat_List,axis=1);
                temp_Lat_List_max = np.max(temp_Lat_List,axis=1);
                temp_Long_List_min = np.min(temp_Long_List,axis=1);
                temp_Long_List_max = np.max(temp_Long_List,axis=1);
                
                #!! this mode isn't fully baked ¯\_(ツ)_/¯ roudning issues it seems, works mostly I think !!
                #--- do the ray tracing ---
                # keo = np.empty( (data_timeUnique.size,keo_N) ,dtype=np.float32 ); #preallocate, don't care what's in it
                # for i in range(0,data_timeUnique.size):
                #     k = np.where(time_limd == data_timeUnique[i])[0]; #gets during a time period
                #     # k = (time_limd == data_timeUnique[i]);
                    
                #     temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
                #     # temp_pplat = pplat_limd[k];
                #     # temp_pplong = pplong_limd[k];
                #     temp_P1Pttan = P1Pttan_limd[k];
                #     temp_P3Pttan = P3Pttan_limd[k];
                #     for j in range(0,keo_N):
                #         kj = (temp_Lat_List_max[j] >= pplat_limd[k]) & (temp_Lat_List_min[j] <= pplat_limd[k]) & (temp_Long_List_max[j] >= pplong_limd[k]) & (temp_Long_List_min[j] <= pplong_limd[k]); #get the stuff that could be in the current box only
                #         temp_vTEC_bit = temp_vTEC[kj]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
                #         # temp_P1Pttan_bit = temp_P1Pttan[kj];
                #         # temp_P3Pttan_bit = temp_P3Pttan[kj];
                #         kl_bit = (P1Pmaxtan[j] >= temp_P1Pttan[kj]) & (P1Pmintan[j] <= temp_P1Pttan[kj]) & (P3Pmaxtan[j] >= temp_P3Pttan[kj]) & (P3Pmintan[j] <= temp_P3Pttan[kj]); #keep the data within the averaging zone
                #         if( np.sum(kl_bit) != 0 ):
                #             keo[i,j] = np.mean(temp_vTEC_bit[kl_bit]); #average the vTEC in this lat band with the given long range
                #             #keo[j] = np.sum(temp_vTEC_bit[kl_bit])/np.sum(kl_bit);'
                #         else:
                #             keo[i,j] = np.nan; #otherwise NaN it instead of 0 it
                #         #END IF
                #     #END FOR j
                # #END FOR i
                keo = GRITI_keo_subfun_raytrace_tan(time_limd,data_timeUnique,pplat_limd,pplong_limd,temp_Lat_List_max,temp_Lat_List_min,temp_Long_List_max,temp_Long_List_min,vTEC_limd,keo_N,P1Pttan,P3Pttan,P1Pmaxtan,P1Pmintan,P3Pmaxtan,P3Pmintan); #call a function to avg all the TEC together into bands
            #END IF
            #save the data
            if( FLG_disableCache == 0 ):
                with open(os.path.join(settings_paths['cache'],'keo_'+cachedName+'.pkl'), 'wb') as keo_cacher:
                    pickle.dump(keo, keo_cacher); #save that data
                #END WITH
            #END IF
        #END IF
        
        # #original, in-fun way (slower w/o jit etc.)
        #     #pre-calc the inpolygon path stuff
        # # inpolgygon = []; #prep a list
        # # for j in range( 0, keo_N):
        # #     inpolgygon.append(Path([(temp_Long_List[j,0],temp_Lat_List[j,1]), (temp_Long_List[j,0], temp_Lat_List[j,2]), \
        # #             (temp_Long_List[j,1], temp_Lat_List[j,2]), (temp_Long_List[j,1], temp_Lat_List[j,2])]));  # square with legs length 1 and bottom left corner at the origin
        # # #END FOR j
        
        # finAnnounce_percentToUpdateAt = 0.5; # every % to update info at
        # finAnnounce_div = np.round(100/finAnnounce_percentToUpdateAt); #calc divisor to use
        # finAnnounce_modRef = np.round(len(data_timeUnique)/finAnnounce_div); #calc mod reference to be used
        # keo = np.zeros( [len(data_timeUnique),keo_N] ,dtype=np.float32 ); #preallocate
        # for i in range(0, len(data_timeUnique) ): #87
        #     #Corral the data to the right place  
        #     k = np.where(time_limd == data_timeUnique[i]); #gets during a time period
            
        #     temp_vTEC = vTEC_limd[k]; #record the vTEC for the time, this is faster than TEC_float[k,0][kl] by farrr
        #     #temp_pplatpplong = np.vstack( [pplong_limd[k], pplat_limd[k]] ).T; #record the pierce-point lat/long for the time in a way the Path function can use
        #     temp_pplat = pplat_limd[k];
        #     temp_pplong = pplong_limd[k];
            
        #     #put in a function for jit power
        #     for j in range(0,keo_N): #keo_N
        #         #average vTEC for a range chunk on an angle
        #         #kl = np.where(inpolgygon[j].contains_points( temp_pplatpplong ))[0];  #get pts inside the area defined
        #         #MATPLOTLIB'S INPOLYGON IS INCORRECT!
        #         #Gets pts inside the range
                
        #         P1PtVxy = np.array([temp_pplong - temp_Long_List[j,0] , temp_pplat - temp_Lat_List[j,0]]); #effort to do own in polygon
        #         P1PtVm = sqrtApproxAo*np.maximum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P1PtVxy[0,:]),np.abs(P1PtVxy[1,:]));
        #         P3PtVxy = np.array([temp_pplong - temp_Long_List[j,2] , temp_pplat - temp_Lat_List[j,2]]); #effort to do own in polygon
        #         P3PtVm = sqrtApproxAo*np.maximum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:])) + sqrtApproxBo*np.minimum(np.abs(P3PtVxy[0,:]),np.abs(P3PtVxy[1,:]));
                
        #         #one shot requires big matrix stuff to happen, a loop might help the time
        #         P1P2P1PtdotCosT = np.einsum('ij,ij->j',np.tile(P1P2Vxy[:,j],[P1PtVxy.shape[1],1]).T,P1PtVxy)/(np.tile(P1P2Vm[j],[P1PtVxy.shape[1],])*P1PtVm);
        #         P1P4P1PtdotCosT = np.einsum('ij,ij->j',np.tile(P1P4Vxy[:,j],[P1PtVxy.shape[1],1]).T,P1PtVxy)/(np.tile(P1P4Vm[j],[P1PtVxy.shape[1],])*P1PtVm);
        #         P3P2P3PtdotCosT = np.einsum('ij,ij->j',np.tile(P3P2Vxy[:,j],[P3PtVxy.shape[1],1]).T,P3PtVxy)/(np.tile(P3P2Vm[j],[P3PtVxy.shape[1],])*P3PtVm);
        #         P3P4P3PtdotCosT = np.einsum('ij,ij->j',np.tile(P3P4Vxy[:,j],[P3PtVxy.shape[1],1]).T,P3PtVxy)/(np.tile(P3P4Vm[j],[P3PtVxy.shape[1],])*P3PtVm);
        #         kl = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);    
                
        #         # kl = np.zeros((P3PtVxy[1,:].size,),dtype=np.bool_); #preallocate    
        #         # for l in range(0,P3PtVxy[1,:].size):
        #         #     P1P2P1PtdotCosT = (P1P2Vxy[0,j]*P1PtVxy[0,l] + P1P2Vxy[1,j]*P1PtVxy[1,l])/(P1P2Vm[j]*P1PtVm[l]);
        #         #     P1P4P1PtdotCosT = (P1P4Vxy[0,j]*P1PtVxy[0,l] + P1P4Vxy[1,j]*P1PtVxy[1,l])/(P1P4Vm[j]*P1PtVm[l]);
        #         #     P3P2P3PtdotCosT = (P3P2Vxy[0,j]*P3PtVxy[0,l] + P3P2Vxy[1,j]*P3PtVxy[1,l])/(P3P2Vm[j]*P3PtVm[l]);
        #         #     P3P4P3PtdotCosT = (P3P4Vxy[0,j]*P3PtVxy[0,l] + P3P4Vxy[1,j]*P3PtVxy[1,l])/(P3P4Vm[j]*P3PtVm[l]);
        #         #     kl[l] = (P1P2P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P1P4P1PtdotCosT >= P1P2P1P4dotCosT[j]) & (P3P2P3PtdotCosT >= P3P2P3P4dotCosT[j]) & (P3P4P3PtdotCosT >= P3P2P3P4dotCosT[j]);
        #         # #END FOR l
                
        #         #if( np.sum(kl) != 0):
        #         keo[i,j] = np.mean(temp_vTEC[kl]); #average the vTEC in this lat band with the given long range
        #     #END FOR j
            
        #     #Generic progress reporter
        #     if( np.mod(i,finAnnounce_modRef) == 0 ):
        #         toc = time.time() - tic; #sec, calc current toc
        #         eta = toc/(finAnnounce_percentToUpdateAt*i/finAnnounce_modRef)*100-toc; #sec, estimate time till finished
        #         print('PROGRESS: Averaging TEC\t %Complete: '+str(np.round(finAnnounce_percentToUpdateAt*i/finAnnounce_modRef,1))+'%\tRuntime: '+str(np.round(toc,2))+' sec/'+str(np.round(toc/60,2))+' min\tETA: '+str(np.round(eta,2))+' sec/'+str(np.round(eta/60,2))+' min'); #announce job, % complete, and time so far
        #     #END IF 
        # #END FOR i
    #END IF
    
    settings_dataSpecific['keo plot latlong name'] = keo_range_plotLatLong_name; #set lat or longitude name
    settings_dataSpecific['keo plot latlong chunks'] = keo_range_plotLatLong_chunks; #set lat or longitude name
    
    toc = time.time() - tic; #sec, calc current toc
    if( FLG_disableText == 0 ):
        print("Time to run: "+str(np.round(toc,2))+" sec/ "+str(np.round(toc/60,2))+" min");
    #END IF
    
    return keo, settings_dataSpecific