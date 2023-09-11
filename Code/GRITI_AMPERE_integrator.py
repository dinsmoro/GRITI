#AMPERE integrator
# 0 = Integrate within plot lat/long range
# 1 = Integrate within plot long range with max latitude of 90
# 2 = Integrate within plot long range with user-defined lat max value
# 3 = Integrate entire hemisphere
# 4 = Integrate within plot long range with user-defined lat min value and max of 90
# 5 = Integrate entire hemisphere within user-defined lat min value and max of 90 (so all long)
# 6 = Integrate entire hemisphere less than user-defined lat value to equator (and all long values)
# 7 = Integrate at a radius around a point (it's back!)
import numpy as np
from Code.subfun_time_to_dateRange import subfun_time_to_dateRange

def GRITI_AMPERE_integrator(AMPERE_data, dates, settings_AMPERE, plotLatRange, plotLongRange, AMPERE_integrateMethod, AMPERE_integrateMethod_val, 
                            AMPERE_integrateMethod_coordType='reg', AMPERE_integrateMethod_coordType_global=None, GRITI_import_AMPERE=None, AMPERE_desired_latLongSteps=[1,15],
                            AMPERE_import_AMPERE_hemi='north', settings_config=None, settings_paths=None,
                            AMPERE_integrateMethod_log=0, AMPERE_integrateMethod_radiusLoc=None, AMPERE_integrateMethod_radius=None):
    #-----Unpack-----
    AMPERE_timeUnique = AMPERE_data['time unique'];
    # AMPERE_integrateMethod = settings_AMPERE['integrate method'];
    # AMPERE_integrateMethod_val = settings_AMPERE['integrate method lat val'];
    # AMPERE_integrateMethod_log = settings_AMPERE['integrate method log'];
    # plotLatRange = settings['map']['lat range'];
    # plotLongRange = settings['map']['long range'];
    AMPERE_dataType = settings_AMPERE['data type']; #get the current AMPERE data type
    
    FLG_paddedDays = ~np.all(dates['date range full dayNum'] == subfun_time_to_dateRange(AMPERE_data['time unique'], dates_zeroHr = dates['date range zero hr dayNum'], FLG_timesWRTzeroHr = False, options = 0)); #see if the time indices exceed the expected date range or not
    
    if( (AMPERE_integrateMethod_coordType != 'reg') & (AMPERE_integrateMethod_coordType != AMPERE_integrateMethod_coordType_global)  ):
        if( GRITI_import_AMPERE == None ):
            #default to a method
            from Code.GRITI_import_AMPERE_direct import GRITI_import_AMPERE_direct as GRITI_import_AMPERE
        #END IF
        AMPERE_data = GRITI_import_AMPERE(dates, settings_paths, AMPERE_import_AMPERE_hemi,
                                          settings_config,AMPERE_coordType=AMPERE_integrateMethod_coordType,AMPERE_desired_latLongSteps=AMPERE_desired_latLongSteps,
                                          FLG_paddedDays=FLG_paddedDays,FLG_dataMix=1,FLG_dataPreference=0,FLG_float64=0); #import pre-processed AMPERE data
        # plotLatRange = 
        # plotLatRange = 
    #END IF
    
    #prep to plot
    if( AMPERE_integrateMethod == 0 ):
        #regular
        kInRange = (AMPERE_data['lat'] >= np.min(plotLatRange)) & (AMPERE_data['lat'] <= np.max(plotLatRange)) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 1 ):
        #max is 90 now
        kInRange = (AMPERE_data['lat'] >= np.min(plotLatRange)) & (AMPERE_data['lat'] <= 90) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 2 ):
        #max is user-defined now
        kInRange = (AMPERE_data['lat'] >= np.min(plotLatRange)) & (AMPERE_data['lat'] <= AMPERE_integrateMethod_val) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 3 ):
        if( (np.min(plotLatRange) <= 0) & (np.max(plotLatRange) >= 0) ):
            kInRange = np.ones(AMPERE_data['lat'].shape,dtype=np.bool_); #all are good to go, both hemisphere integration yolo
        else:
            if( (np.min(plotLatRange) >= 0) & (np.max(plotLatRange) >= 0) ):
                #northern hemisphere
                kInRange = AMPERE_data['lat'] >= 0; #get data in the range
            else:
                #southern hemisphere
                kInRange = AMPERE_data['lat'] <= 0; #get data in the range
            #END IF
        #END IF
    elif( AMPERE_integrateMethod == 4 ):
        kInRange = (AMPERE_data['lat'] <= 90) & (AMPERE_data['lat'] >= AMPERE_integrateMethod_val) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 5 ):
        kInRange = (AMPERE_data['lat'] <= 90) & (AMPERE_data['lat'] >= AMPERE_integrateMethod_val); #get data in the range
    elif( AMPERE_integrateMethod == 6 ):
        kInRange = (AMPERE_data['lat'] <= AMPERE_integrateMethod_val); #get data in the range
    elif( AMPERE_integrateMethod == 7 ):
        if( (AMPERE_integrateMethod_radius == None) | (AMPERE_integrateMethod_radiusLoc == None) ):
            print('ERROR IN GRITI_AMPERE_integrator: AMPERE_integrateMethod 7 chosen but AMPERE_integrateMethod_radius OR AMPERE_integrateMethod_radiusLoc NOT provided. Crashing.');
            import sys
            sys.crash()
        #END IF
        # geopy was too slow for this
        # try:
        #     import geopy.distance as geodist
        #     import time
        #     tic = time.time(); #start
        #     AMPERE_integrateMethod_radius_dist = np.zeros(AMPERE_data['lat'].size); 
        #     for i in range(0,AMPERE_data['lat'].size):
        #         AMPERE_integrateMethod_radius_dist[i] = geodist.distance(AMPERE_integrateMethod_radiusLoc,(AMPERE_data['lat'][i],AMPERE_data['long'][i])).km;
        #         if( np.mod(i+1,AMPERE_data['lat'].size//10) == 0 ):
        #             print('10% more done at '+str(np.round(time.time()-tic,2)));
        #         #END IF
        #     #END FOR i
        #     kInRange = (AMPERE_data['lat'] <= AMPERE_integrateMethod_val); #get data in the range
        # except:
            # print('WARNING in GRITI_AMPERE_integrator: Geopy not installed. Using less advanced distance formula.');
        from Code.subfun_dist import distWGS84_notBad as geodist
        AMPERE_integrateMethod_radius_dist = geodist(AMPERE_integrateMethod_radiusLoc[0], AMPERE_integrateMethod_radiusLoc[1], AMPERE_data['lat'], AMPERE_data['long']);
        kInRange = (AMPERE_integrateMethod_radius_dist <= AMPERE_integrateMethod_radius); #get data in the range
        #END TRY
    #END IF
    
    kNotNan = ~np.isnan(AMPERE_data[AMPERE_dataType]); #get where NaNs are NOT at
    #--- this took much longer to run than I expected it to, so I tried ---
    import importlib
    if( (importlib.find_loader('cupy') != None) & settings_config['acceleration']['use_GPU'] ):
        #cupy was the fastest by far (impressive for "simple" math), needs an nvidia GPU [single threaded]
        import cupy
        AMPERE_time_inRange = cupy.array(AMPERE_data['time'][kInRange&kNotNan]); #limit it here once
        AMPERE_data_inRange = cupy.array(AMPERE_data[AMPERE_dataType][kInRange&kNotNan]); #limit here once
        AMPERE_timeUnique_GPU = cupy.array(AMPERE_timeUnique);
        AMPERE_integrate = cupy.empty( AMPERE_timeUnique_GPU.size , dtype=AMPERE_data_inRange.dtype); #prep integrated joule heating
        for i in range(AMPERE_timeUnique_GPU.size):
            AMPERE_integrate[i] = cupy.sum(AMPERE_data_inRange[cupy.equal(AMPERE_timeUnique_GPU[i], AMPERE_time_inRange)]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
        #END FOR i
        AMPERE_integrate = cupy.asnumpy(AMPERE_integrate); #convert to a numpy array
        del AMPERE_time_inRange, AMPERE_data_inRange, AMPERE_timeUnique_GPU
        cupy._default_memory_pool.free_all_blocks(); #clean GPU memory I hope
        
    elif( importlib.find_loader('numba') != None ):
        #numba and joblib are (very oddly) neck-and-neck, but numba has less moving parts (on my end) [parallel]
        from numba import jit, prange
        
        @jit(nopython=True,nogil=True,parallel=True,fastmath=True)
        def GRITI_AMPERE_integrator_summer(AMPERE_dataNow, AMPERE_time, AMPERE_timeUnique):
            AMPERE_integrate = np.empty( AMPERE_timeUnique.size , dtype=AMPERE_dataNow.dtype); #prep integrated joule heating
            for i in prange(AMPERE_timeUnique.size):
                AMPERE_integrate[i] = np.sum(AMPERE_dataNow[AMPERE_timeUnique[i] == AMPERE_time]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
            #END FOR i
            return AMPERE_integrate
        #END DEF
        
        AMPERE_integrate = GRITI_AMPERE_integrator_summer(AMPERE_data[AMPERE_dataType][kInRange&kNotNan], AMPERE_data['time'][kInRange&kNotNan], AMPERE_timeUnique); #numba loop to jooce it
        
    elif( importlib.find_loader('joblib') != None ):
        #joblib was actually performant using an "implied" variable system (e.g. non-explicitly-defined "global" variables) [parallel]
        import joblib
        from multiprocessing import cpu_count
        parallel_numThreads = cpu_count(); #use multiprocess to get # of CPU threads (can't differentiate between cores/threads so will give # of threads only) - built-in and ez

        AMPERE_time_inRange = AMPERE_data['time'][kInRange&kNotNan]; #limit it here once
        AMPERE_data_inRange = AMPERE_data[AMPERE_dataType][kInRange&kNotNan]; #limit here once
        parallel_splitterIndexes = np.int64(np.round(np.linspace(0,AMPERE_timeUnique.size,parallel_numThreads+1,endpoint=True))); #split up the indexes to be parallelized
        parallel_list = [None for j in range(0,parallel_splitterIndexes.size-1)]; #Prep
        for j in range(0,parallel_splitterIndexes.size-1):
            parallel_list[j] = AMPERE_timeUnique[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]]; #load in the variables needed
        #END FOR j
        
        def GRITI_AMPERE_integrator_summerJL(AMPERE_timeUnique):
            AMPERE_integrate = np.empty( AMPERE_timeUnique.size , dtype=AMPERE_data_inRange.dtype); #prep integrated joule heating
            for i in range(AMPERE_timeUnique.size):
                AMPERE_integrate[i] = np.sum(AMPERE_data_inRange[AMPERE_timeUnique[i] == AMPERE_time_inRange]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
            #END FOR i
            return AMPERE_integrate
        #END DEF
        
        with joblib.parallel_backend('loky'):
            parallel_return = joblib.Parallel(n_jobs=parallel_numThreads)(joblib.delayed(GRITI_AMPERE_integrator_summerJL)(k) for k in parallel_list); #will this not destroy the world?
        #END WITH
        
        AMPERE_integrate = np.empty( AMPERE_timeUnique.size , dtype=AMPERE_data_inRange.dtype); #prep integrated joule heating
        for j in range(0,parallel_splitterIndexes.size-1):
            AMPERE_integrate[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_return[j]; #unpack
        #END FOR j
    else:
        #otherwise straight numpy [single threaded]
        print('WARNING in GRITI_AMPERE_integrator: cupy (NVIDIA GPU only!)/numba/joblib all aren\'t available and the fallback code is v slow. Get one of them!');
        AMPERE_time_inRange = AMPERE_data['time'][kInRange&kNotNan]; #limit it here once
        AMPERE_data_inRange = AMPERE_data[AMPERE_dataType][kInRange&kNotNan]; #limit here once
        AMPERE_integrate = np.empty( AMPERE_timeUnique.size , dtype=AMPERE_data_inRange.dtype); #prep integrated joule heating
        for i in range(AMPERE_timeUnique.size):
            AMPERE_integrate[i] = np.sum(AMPERE_data_inRange[AMPERE_timeUnique[i] == AMPERE_time_inRange]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
        #END FOR i
    #END IF
            
    if( AMPERE_integrateMethod_log == 1 ):
        #Can't log 0's, interpolate over them
        k = np.where( AMPERE_integrate == 0 )[0]; #get where 0's are
        AMPERE_integrate[k] = 1; #set to 1 to prevent error
        AMPERE_integrate = np.log10(AMPERE_integrate); #log it
    #END IF
    
    return AMPERE_integrate
#END DEF

# @jit(nopython=True,nogil=True,parallel=True,fastmath=True) #{'nnan','ninf','nsz','contract','afn'}
# def GRITI_AMPERE_integrator_summer(AMPERE_dataNow, AMPERE_time, AMPERE_timeUnique, kInRange):
#     AMPERE_integrate = np.empty( AMPERE_timeUnique.size , dtype=AMPERE_dataNow.dtype); #prep integrated joule heating
#     for i in prange(AMPERE_timeUnique.size):
#         k = AMPERE_timeUnique[i] == AMPERE_time; #get the right time
#         AMPERE_integrate[i] = np.sum(AMPERE_dataNow[k&kInRange]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
#     #END FOR i
#     return AMPERE_integrate
# #END DEF



# @jit(nopython=True,nogil=True,parallel=False,fastmath=True)
# def GRITI_AMPERE_integrator_summer2S(AMPERE_dataNow, AMPERE_time, AMPERE_timeUnique):
#     AMPERE_integrate = np.empty( AMPERE_timeUnique.size , dtype=AMPERE_dataNow.dtype); #prep integrated joule heating
#     for i in prange(AMPERE_timeUnique.size):
#         AMPERE_integrate[i] = np.sum(AMPERE_dataNow[AMPERE_timeUnique[i] == AMPERE_time]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
#     #END FOR i
#     return AMPERE_integrate
# #END DEF