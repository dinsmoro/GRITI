# supporting functions for GRITI_import_TEC_xxx
import numpy as np
from Code.subfun_textNice import textNice
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import h5py
from time import time
from Code.subfun_strfind import strfind
import sys, os
import psutil, joblib

def GRITI_import_TEC_support_filterVersion(): #sends the verison how illegal is this lolo I'm not gonna use a class
    version_filt = 4.0; #filtered algorithm version
    #1 10/8/2019 - initial algorithm
    #1.1 9/11/2020 - fixed non-0/30 time step handling (29/59 were big ones)    
    #1.2 9/16/2020 - removed highpass filter, only dampens data
    #1.3 9/17/2020 - outlier control introduced
    #1.4 9/17/2020 - linear interpolation to make savgol filter work correctly for gappy data
    #2.0 9/17/2020 - interpolated data used (large groups of receivers would not have data at the same time, leaving vertical "streaks" of no-data)
    #3.0 4/9/2021 - moved from monolithic data bloc to slices of data, moved to dict royale
    #3.1 11/30/2021 - fixed missing pierceAlt info & added geomagnetic coordinate support
    #4.0 2/8/2022 - moved to interpolation for most time keeping controls during a pass & supports different satellite types
    return version_filt;
#END DEF

def GRITI_import_TEC_support_filter(settings_paths,settings_config,paths_TEC,TEC_dataAvail,TEC_dataFileName,TEC_dataFileNameUnfilt,dateRange_full,dateRange_dayNum_full,dateRange_dayNum_full_orig, \
                                    TEC_dataRate,TEC_timeTolerance,TEC_maxAmpAllowed,filter_savGolPeriod,order_savGol,minElevation,minimumTimeGap,dataType_meth, \
                                    TEC_dataFilePathUnfilt, TEC_dataAgg_timeAdditionLimit, TEC_dataAgg_distToPts_degcSq, deltaTEC_compareValue, \
                                    FLG_deleteUnfilt,FLG_reqPaddedDays,FLG_dataAggregation):
    #==============Filter Data after download & conversion==============
    version_unfiltReq = 3.0; #requires this unfilt version to function
    version_filt = GRITI_import_TEC_support_filterVersion(); #get function version from elsewhere
    
    # establish which days need to to be converted
    instances_toDo = np.zeros(TEC_dataAvail.size, dtype=np.bool_); #prep
    for i in range(0,len(dateRange_dayNum_full[:,0])): #just unpack needed data
        if( TEC_dataAvail[i] == 4 ): #4 means data was downloaded and converted but needs to be filtered (this is where 7 jumps off)
            instances_toDo[i] = True; # gotta do it
        #END IF
    #END FOR i
    instances_toDo_where = np.where(instances_toDo)[0]; #get where we gotta do
    
    if( instances_toDo_where.size > 0 ):
        #--- estimate the memory needed to filter 1 day ---
        i = instances_toDo_where[0]; #start off with i=first one
        j = 0; #start off with j=0
        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i+j], 'r', rdcc_nbytes =500*1024*1024) as TEC_fileUnfilt:
            #----- Read in the unfiltered data, size it up -----
            unfilt_dictTemp_elev = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
            TEC_goodVals = unfilt_dictTemp_elev >= minElevation; #get the locations where elevation is above the min elevation
            # TEC_goodVals = np.ones(TEC_fileData_float_temp[:,locFloatUnfilt_elev].shape,dtype=bool); #set all to good [elevation clipped later]
        #END WITH
        TEC_dataSize = np.int64(TEC_goodVals.sum()); #get the data estimated data size
        TEC_dataSize += np.int64(TEC_dataSize*1/2); #add on half to represent the two side days for calc'n the edges good
        del unfilt_dictTemp_elev, TEC_goodVals #save memory
        floatBytes = np.dtype(dataType_meth).itemsize; #get the bytes used by the floats, which can be 64
        TEC_dataMem = TEC_dataSize*2*6 + TEC_dataSize*floatBytes*6 + TEC_dataSize*8*1 + TEC_dataSize*1*1 + TEC_dataSize*4*1 + \
            TEC_dataSize*2*6 + TEC_dataSize*floatBytes*5 + TEC_dataSize*4*1 + TEC_dataSize*1*1 + TEC_dataSize*4*1; #bytes, valid for unfilt 3.0 and filt 4.0 (characters are 1 byte per character in numpy)
        
        comp_availMem = np.int64(psutil.virtual_memory().total*.9) #bytes, use 90% of available RAM (.available is also an option to use what's available - but I feel we gotta go hard to do the big calcs)
        
        numInstances_memLim = np.uint32(comp_availMem//TEC_dataMem); # number of instances that can be fit
        
        if( numInstances_memLim < 1 ):
            numInstances_toDo = 1; #just one, no parallel possible due to memory limit
        else:
            if( numInstances_memLim > settings_config['parallel num threads'] ):
                numInstances_toDo = settings_config['parallel num threads']; #enough memory to use all threads
            else:
                numInstances_toDo = numInstances_memLim; #just use as much memory as possible, threads won't be filled out
            #END IF
        #END IF
        if( instances_toDo_where.size < numInstances_toDo ):
            numInstances_toDo = instances_toDo_where.size; #limit to number of days to work on if that's under the computation limit
        #END IF
        
        if( numInstances_toDo > 1 ):
            #--- pack up for parallel ---
            print('WARNING in GRITI_import_TEC_support_filter: Parallelization possible for converting unfiltered data to filtered data, which will speed it up. Downside is that no updates will be printed. See you on the other side!');
            tic = time();
            parallel_list = []; #Prep
            for i in range(0, instances_toDo_where.size): # Every iteration appends a list of inputs to parallel_list, each index of parallel_list will be independently run
                parallel_list.append([ instances_toDo_where[i], dateRange_dayNum_full, dateRange_dayNum_full_orig, dateRange_full, \
                                       TEC_dataFileNameUnfilt, TEC_dataFileName, \
                                       TEC_dataAvail,  TEC_dataRate, TEC_timeTolerance, TEC_maxAmpAllowed, \
                                       minElevation, minimumTimeGap, dataType_meth, \
                                       filter_savGolPeriod, order_savGol, version_filt, version_unfiltReq, \
                                       settings_paths, paths_TEC, TEC_dataFilePathUnfilt, \
                                       TEC_dataAgg_timeAdditionLimit, TEC_dataAgg_distToPts_degcSq, deltaTEC_compareValue, \
                                       FLG_reqPaddedDays, FLG_dataAggregation, FLG_deleteUnfilt, True ]);
            #END FOR i
            
            #--- apply parallel calc on function ---
            with joblib.parallel_backend('loky'):
                with joblib.Parallel(n_jobs=numInstances_toDo,pre_dispatch=numInstances_toDo,batch_size=1) as parallel_arbiter:    
                    parallel_results = parallel_arbiter(joblib.delayed(GRITI_import_TEC_support_filter_perDay)(*sublist) for sublist in parallel_list); #annoyingly I need to unpack every variable input manually
                #END WITH
            #END WITH
            del parallel_list #save some mem
            #--- unpack parallel results ---
            for i in range(0,len(parallel_results)):
                TEC_dataAvail[instances_toDo_where[i]] = parallel_results[i];
            #END FOR i
            del parallel_results #save some mem
            print('\nTime to filter '+str(instances_toDo_where.size)+' days in parallel took: '+str(np.round((time()-tic)/60,2))+' min\n'); #extra space at end  
        else:
            # no need to parallelize if just 1 instance can run
            for i in range(0, instances_toDo_where.size):
                TEC_dataAvail[instances_toDo_where[i]] = GRITI_import_TEC_support_filter_perDay(instances_toDo_where[i], dateRange_dayNum_full, dateRange_dayNum_full_orig, dateRange_full, \
                                    TEC_dataFileNameUnfilt, TEC_dataFileName, \
                                    TEC_dataAvail,  TEC_dataRate, TEC_timeTolerance, TEC_maxAmpAllowed, \
                                    minElevation, minimumTimeGap, dataType_meth, \
                                    filter_savGolPeriod, order_savGol, version_filt, version_unfiltReq, \
                                    settings_paths, paths_TEC, TEC_dataFilePathUnfilt, \
                                    TEC_dataAgg_timeAdditionLimit, TEC_dataAgg_distToPts_degcSq, deltaTEC_compareValue, \
                                    FLG_reqPaddedDays, FLG_dataAggregation, FLG_deleteUnfilt, False);
            #END FOR i
        #END IF
    #END IF
    
    return TEC_dataAvail
#END DEF

# this as a function enables it to be parallelized!
def GRITI_import_TEC_support_filter_perDay(i, dateRange_dayNum_full, dateRange_dayNum_full_orig, dateRange_full, \
                                           TEC_dataFileNameUnfilt, TEC_dataFileName, \
                                           TEC_dataAvail,  TEC_dataRate, TEC_timeTolerance, TEC_maxAmpAllowed, \
                                           minElevation, minimumTimeGap, dataType_meth, \
                                           filter_savGolPeriod, order_savGol, version_filt, version_unfiltReq, \
                                           settings_paths, paths_TEC, TEC_dataFilePathUnfilt, \
                                           TEC_dataAgg_timeAdditionLimit, TEC_dataAgg_distToPts_degcSq, deltaTEC_compareValue, \
                                           FLG_reqPaddedDays, FLG_dataAggregation, FLG_deleteUnfilt, FLG_parallel):
            
    #-----FILTER THE DATA-----
    if( TEC_dataAvail[i] == 4 ): #4 means data was downloaded and converted but needs to be filtered (this is where 7 jumps off)
        
        #prep the Sav-Gol filter for debiasing
        windowLen_savGol = np.int64(np.round(filter_savGolPeriod/TEC_dataRate)); #window length, 60 minutes in seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
        #from conversations with AC ^
        if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
            windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
        #END IF
        
        #-----Declare some constants used to force the data to the required data rate-----
        TEC_dataRate_allowedStamps = np.arange(0,60,TEC_dataRate); #sec, make an array of allowed second timestamps
        TEC_dataRate_allowedStamps60_sameOrder = np.copy(TEC_dataRate_allowedStamps); #copy
        TEC_dataRate_allowedStamps60_sameOrder[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
        # TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60_sameOrder); #sort it 
        
        #TIME TO IMPORT DATA TO FILTER!
        if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == True ): #if statement to find original days requested
            if( FLG_parallel == False ):
                print("\nFiltering {} to {}.".format(TEC_dataFileNameUnfilt[i],TEC_dataFileName[i]) );
                tic = time(); #for time testing
            #END IF
            
            TEC_fileData_paddingWarning = []; #prime this list
            unfilt_dict = {}; #prep this
            for j in range(-1,2): #run through the days to get the main day and the days around it
                if( (TEC_dataAvail[i+j] == 4) | (TEC_dataAvail[i+j] == 7) | (TEC_dataAvail[i+j] == -3) ): #make sure data is there before reading it (4 or 7 mean it is there) and -3 also means it is there, just a padded day
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i+j], 'r', rdcc_nbytes =500*1024*1024) as TEC_fileUnfilt:
                        #----- Read in the unfiltered data, keep only the good data -----
                        keyz = list(unfilt_dict.keys()); #get the saved keys
                        keyzNew = list(TEC_fileUnfilt.keys()); #get the saved keys
                        if( j == -1 ):
                            #j = -1 is guaranteed to be 1st, so can assume unfilt_dict is unfilled yet
                            unfilt_dict['time'] = TEC_fileUnfilt.get('time')[()]; #get that dataset out
                            unfilt_dict['elev'] = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
                            TEC_goodVals = unfilt_dict['time']/3600 > 18; #get UTC hours 18 and greater (6 hour buffer) to cut down on the data imported (would have to think more about the max theoretical observable GPS satellite time period - 6 is safe.)
                            TEC_goodVals = (unfilt_dict['elev'] >= minElevation) & TEC_goodVals; #get the locations where elevation is above the min elevation
                            unfilt_dict['time'] = [unfilt_dict['time'][TEC_goodVals]] #keep the good stuff
                            unfilt_dict['elev'] = [unfilt_dict['elev'][TEC_goodVals]]; #keep the good stuff
                            keyzNew.remove('time'); #remove from the list, manually got it for data corralling
                            keyzNew.remove('elev'); #remove from the list, manually got it for data corralling
                        elif( j == 1 ):
                            unfilt_dictTemp_time = TEC_fileUnfilt.get('time')[()]; #get that dataset out
                            unfilt_dictTemp_elev = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
                            TEC_goodVals = unfilt_dictTemp_time/3600 < 6; #get UTC hours 6 and less (6 hour buffer) to cut down on the data imported (would have to think more about the max theoretical observable GPS satellite time period - 6 is safe.)
                            TEC_goodVals = (unfilt_dictTemp_elev >= minElevation) & TEC_goodVals; #get the locations where elevation is above the min elevation
                            unfilt_dictTemp_time = unfilt_dictTemp_time[TEC_goodVals]; #keep the good stuff
                            unfilt_dictTemp_elev = unfilt_dictTemp_elev[TEC_goodVals]; #keep the good stuff
                            keyzNew.remove('time'); #remove from the list, manually got it for data corralling
                            keyzNew.remove('elev'); #remove from the list, manually got it for data corralling
                            if( strfind(keyz,'time',opt=1) > 0 ):
                                unfilt_dict['time'].append(unfilt_dictTemp_time); #tack that dataset on
                                del unfilt_dictTemp_time; #clean up the memory
                            else:
                                #otherwise it's a new data type to add in
                                unfilt_dict['time'] = [np.copy(unfilt_dictTemp_time)]; #get that dataset out
                                del unfilt_dictTemp_time; #clean up the memory
                            #END IF
                            if( strfind(keyz,'elev',opt=1) > 0 ):
                                unfilt_dict['elev'].append(unfilt_dictTemp_elev); #tack that dataset on
                                del unfilt_dictTemp_elev; #clean up the memory
                            else:
                                #otherwise it's a new data type to add in
                                unfilt_dict['elev'] = [np.copy(unfilt_dictTemp_elev)]; #get that dataset out
                                del unfilt_dictTemp_elev; #clean up the memory
                            #END IF
                        else:
                            unfilt_dictTemp_elev = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
                            TEC_goodVals = unfilt_dictTemp_elev >= minElevation; #get the locations where elevation is above the min elevation
                            # TEC_goodVals = np.ones(TEC_fileData_float_temp[:,locFloatUnfilt_elev].shape,dtype=bool); #set all to good [elevation clipped later]
                            unfilt_dictTemp_elev = unfilt_dictTemp_elev[TEC_goodVals]; #keep the good stuff
                            keyzNew.remove('elev'); #remove from the list, manually got it for data corralling
                            if( strfind(keyz,'elev',opt=1) > 0 ):
                                unfilt_dict['elev'].append(unfilt_dictTemp_elev); #tack that dataset on
                                del unfilt_dictTemp_elev; #clean up the memory
                            else:
                                #otherwise it's a new data type to add in
                                unfilt_dict['elev'] = [np.copy(unfilt_dictTemp_elev)]; #get that dataset out
                                del unfilt_dictTemp_elev; #clean up the memory
                            #END IF
                        #END IF
                        for k in range(0,len(keyzNew)):
                            if( strfind(keyz,keyzNew[k],opt=1) > 0 ):
                                unfilt_dict[keyzNew[k]].append(TEC_fileUnfilt.get(keyzNew[k])[()][TEC_goodVals]); #tack that dataset on
                            else:
                                #otherwise it's a new data type to add in
                                unfilt_dict[keyzNew[k]] = [TEC_fileUnfilt.get(keyzNew[k])[()][TEC_goodVals]]; #get that dataset out, only keep the good stuff
                            #END IF
                        #END FOR k
                        #--- Read the attributes in ---
                        keyzNew = list(TEC_fileUnfilt.attrs.keys()); #get the attribute keys
                        # keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                        for k in range(0,len(keyzNew)):
                            if( strfind(keyz,keyzNew[k],opt=1) > 0 ):
                                if( np.isclose(unfilt_dict[keyzNew[k]],TEC_fileUnfilt.attrs[keyzNew[k]]) == False ):
                                    #only worry is if the attribute isn't consistent
                                    print('-----Warning-----');
                                    print('Attribute '+keyzNew[k]+' isn\'t the same as the previously recorded value from another file of '+ \
                                        str(unfilt_dict[keyzNew[k]])+' and this file\'s value of '+str(TEC_fileUnfilt.attrs[keyzNew[k]])+ \
                                        '.\n NaN\'ing it and try to sort that out.');
                                    unfilt_dict[keyzNew[k]] = np.nan; #nan that attribute, figure it out later
                                #END IF
                            else:
                                unfilt_dict[keyzNew[k]] = TEC_fileUnfilt.attrs[keyzNew[k]]; #get that attribute out
                            #END IF
                        #END FOR k
                        del TEC_goodVals #clean up that memory
                    #END WITH
                                        
                elif(FLG_reqPaddedDays == 0): #otherwise data isn't there - because strict padding is off and data wasn't there
                    if( j == -1 ): #before padding not there
                        TEC_fileData_paddingWarning.append("Before"); #append to var
                    elif( j == 1): #after padding not there
                        TEC_fileData_paddingWarning.append("After"); #append to var
                    else: #something bad
                        print("\n==============ERROR==============");
                        print("There's no data on a main day requested {}/{} (Y/#D) - which shouldn't have happened. Help. Exiting.\n".format(dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
                        sys.crash(); #yolo outa here
                    #END TRY
                else:
                    print("\n==============ERROR==============");
                    print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are required via passed flag.\nPrinting availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i+j,0],dateRange_full[i+j,1],dateRange_full[i+j,2]) );
                    print("{}".format(TEC_dataAvail)); #print for error - lets user know available days
                    print("{}".format(i+j)); #print for error - lets user know available days
                    # print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                    TEC_dataAvail[i+j] = -1; #note the error
                    sys.crash(); #yolo out
                #END IF              
            #END FOR j
            #--- Convert to numpy arrays from the dynamically added lists ---
            keyz = list(unfilt_dict.keys()); #get the current keys
            for j in range(0,len(keyz)):
                if( np.isscalar(unfilt_dict[keyz[j]]) == False ):
                    #if not a scalar, apply the logical mask
                    unfilt_dict[keyz[j]] = np.hstack(unfilt_dict[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
                #END IF
            #END FOR j
            #--- Make sure req'd unfilt version is here ---
            if( version_unfiltReq > unfilt_dict['version'] ):
                print("\n==============ERROR==============");
                print('Unfiltered version is '+str(unfilt_dict['version'])+' but required unfiltered verison is '+str(version_unfiltReq)+'. Mitigate in import_TEC function and re-run.');
                sys.crash(); #yolo out
            #END IF
            
            if( FLG_parallel == False ):
                toc = time() - tic; #for time testing
                print("Time to import all needed data into memory: {} min".format(round(toc/60,2))); #extra space at end
                tic = time(); #for time testing
            #END IF
            
            #==============THIS WILL GO IN A FUNCTION FOR MAX SPEED==============
            #TIME TO FILTER!
            
            #-----prep dynamic list to hold data as it is accepted-----
            filt_dict = {
                'sat':[],
                'year':[],
                'dayNum':[],
                'hour':[],
                'min':[],
                'sec':[],
                'time':[],
                'lat':[],
                'long':[],
                'elev':[],
                'dTEC':[],
                'dTECerr':[],
                'site':[],
                'satType':[],
                }; #prep a dict to hold the data as it's made
            
#            unfilt_site_unique = np.unique(unfilt_dict['site']); #get unique site names
            unfilt_site_unique, unfilt_site_unique_indexes = np.unique(unfilt_dict['site'] , return_inverse=True); #get unique site names and the indexes to speed up searching for where to look
            #inspired by https://stackoverflow.com/questions/33281957/faster-alternative-to-numpy-where
#            unfilt_site_unique_sortedIndex = np.argsort(unfilt_dict['site'], kind='mergesort'); #returns the indexes IF the array unfilt_dict['site'] was sorted
            unfilt_site_unique_currentSiteArray = np.split(np.argsort(unfilt_dict['site'], kind='mergesort'), np.cumsum(np.bincount(unfilt_site_unique_indexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
            #faster than [for j in range... currentSite_loc = np.where(j == unfilt_site_unique_indexes)[0]; #get the data locations for 1 site per]
            
            if( FLG_dataAggregation == 1 ):
                TEC_timeUnique, TEC_timeUniqueIndexes , TEC_timeUniqueCount = np.unique( unfilt_dict['time'] , return_inverse=True , return_counts=True);
                #Cut off time stamps with very little data in the range selected
    #            TEC_dataAvgNum = np.sum(TEC_timeUniqueCount)/TEC_timeUnique.size; #average data per time
    #            TEC_dataLim = np.round(TEC_dataLimPercent*TEC_dataAvgNum); #min number before eject time set
    #            TEC_timeUnique = TEC_timeUnique[ TEC_timeUniqueCount > TEC_dataLim ]; #remove the low data # stuff
                #disabled for now
                TEC_timeUnique_currentTimeArray = np.split(np.argsort(unfilt_dict['time'], kind='mergesort'), np.cumsum(np.bincount(TEC_timeUniqueIndexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
            #END IF
            if( FLG_parallel == False ):
                print("Time to prep for site calcs: {} min".format(round((time() - tic)/60,2))); #extra space at end
                tic2 = time(); #get a timer for just calcs
            #END IF
            
            # cntr = 0; #for controlling the number of plots during debugging
            if( FLG_dataAggregation == 0 ): #with this off, just use what is directly there
                #cruise through every site (from the unique of unfilt_dict['site'])
                for j in range(0,unfilt_site_unique.size):
                    # print('unfilt_site_unique.size = '+str(unfilt_site_unique.size))
                    currentSite_loc = unfilt_site_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                    _, currentSat_uniqueInv = np.unique(unfilt_dict['sat'][currentSite_loc], return_inverse=True); #get the unique sats in at that current site
                    _, currentSatType_uniqueInv = np.unique(unfilt_dict['satType'][currentSite_loc], return_inverse=True); #get the unique sats in at that current site
                    currentSatNSatType = np.vstack((currentSat_uniqueInv,currentSatType_uniqueInv)).T; #combined, use indexes as a proxy since can't make int16 and S1 types touch ez
                    currentSatNSatType_unique = np.unique( currentSatNSatType, axis=0); #unique 
                    
                    #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])             
                    for k in range(0,currentSatNSatType_unique.shape[0]):
                        # print('currentSat_unique.size = '+str(currentSat_unique.size))
                        currentSat_loc = np.where( np.all(currentSatNSatType_unique[k,:] == currentSatNSatType, axis=1) )[0]; #get the data locations for 1 sat at that 1 site
                        #-----GET NEEDED LOCAL VARS-----
                        #---INTS---
                        currentSat = unfilt_dict['sat'][currentSite_loc][currentSat_loc]; #sat#, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentYear = unfilt_dict['year'][currentSite_loc][currentSat_loc]; #year, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentDayNum = unfilt_dict['dayNum'][currentSite_loc][currentSat_loc]; #day#, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentHour = unfilt_dict['hour'][currentSite_loc][currentSat_loc]; #hr, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentMin = unfilt_dict['min'][currentSite_loc][currentSat_loc]; #min, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentSec = unfilt_dict['sec'][currentSite_loc][currentSat_loc]; #sec, get the day data for 1 sat at that 1 site - uses abs illegal double indexing                        
                        #---FLOATS---
                        currentLat = unfilt_dict['lat'][currentSite_loc][currentSat_loc];
                        currentLong = unfilt_dict['long'][currentSite_loc][currentSat_loc];
                        currentElv = unfilt_dict['elev'][currentSite_loc][currentSat_loc];
                        #---STRINGS---
                        currentSite = unfilt_dict['site'][currentSite_loc][currentSat_loc];
                        currentSatType = unfilt_dict['satType'][currentSite_loc][currentSat_loc];
                        #---UNFILT DATA---
                        currentvTEC = unfilt_dict['vTEC'][currentSite_loc][currentSat_loc]; #TECU, get the vTEC data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentsTECerror = unfilt_dict['sTECerr'][currentSite_loc][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentSecTotal = np.int64(unfilt_dict['dayNum'][currentSite_loc][currentSat_loc])*86400 + unfilt_dict['time'][currentSite_loc][currentSat_loc]; #sec, calculate days/hour/min/sec into total seconds but do it at 64 bit for extra good precision in the seconds [no overflow]
                        # currentsTEC = TEC_fileData_float[currentSite_loc,locFloatUnfilt_sTEC][currentSat_loc]; #TECU, get the sTEC data for 1 sat at that 1 site - uses abs illegal double indexing [only for debug]
                        
                        currentTimeSplits_loc = np.append( np.insert( np.where(np.diff(currentSecTotal) > minimumTimeGap)[0]+1 , 0, 0), currentSecTotal.shape ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                        
                        for l in range(0,len(currentTimeSplits_loc)-1):
                            # print('len(currentTimeSplits_loc)-1 = '+str(len(currentTimeSplits_loc)-1))
                            #-----Identify current TEC_dataRate & prep the filters-----
                            if( currentDayNum[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]].size > 1 ):
                                current_TEC_dataRate = np.diff(currentSecTotal[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #sec, get the current data rate for the site
                                uni, cnts = np.unique(current_TEC_dataRate,return_counts=True); #get unique rates and their counts
                                current_TEC_dataRate = uni[cnts >= np.mean(cnts)]; #choose the counts that are larger than the mean (heavily influenced by outlier uniques with big counts (most likely the correct data rate))
                                if( current_TEC_dataRate.size == 1 ):
                                    current_TEC_dataRate = np.int64(current_TEC_dataRate.item()); #get the current data rate as an integer
                                else:
                                    current_TEC_dataRate = np.int64(np.min(current_TEC_dataRate)); #choose the minimum one, assuming 15's and 30's mean that the real time is 15 and there's just a lot of skips
                                #END IF
                            else:
                                current_TEC_dataRate = 0; #set current data rate to 0 so code doesn't run on it 
                            #END IF
                            
                            #make sure do work only on the day we want
                            #only get the data for times ranges thare are two times the filter period or longer
                            if( ( (currentSecTotal[currentTimeSplits_loc[l+1]-1]-currentSecTotal[currentTimeSplits_loc[l]]) > filter_savGolPeriod*2) & #make sure the data time period is 2*filterPeriod
                            (np.any(currentDayNum[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]] == dateRange_dayNum_full[i,1]) == True) &
                            (current_TEC_dataRate != 0) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                #-----GET NEEDED LOCAL-LOCAL VARS-----
                                #---INTS---
                                currentSat_singleSatLine = currentSat[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentYear_singleSatLine = currentYear[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentDayNum_singleSatLine = currentDayNum[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentHour_singleSatLine = currentHour[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentMin_singleSatLine = currentMin[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentSec_singleSatLine = currentSec[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                #---FLOATS---
                                currentLat_singleSatLine = currentLat[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentLong_singleSatLine = currentLong[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentElv_singleSatLine = currentElv[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                #---STRINGS---
                                currentSite_singleSatLine = currentSite[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentSatType_singleSatLine = currentSatType[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                #---UNFILT DATA---
                                currentvTEC_singleSatLine = np.float64(currentvTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #this will increase readability
                                currentsTECerror_singleSatLine = currentsTECerror[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentSecTotal_singleSatLine = currentSecTotal[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                # currentsTEC_singleSatLine = np.float64(currentsTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #this will increase readability [only for debug]
                                
                                #-----This code bit exists to catch if some weird duplicate timestamp thing happens-----
                                #it happens
                                if( np.unique(currentSecTotal_singleSatLine).size != currentSecTotal_singleSatLine.size ):
                                    # print(str(unfilt_site_unique[j])+' sat# '+str(currentSatNSatType_unique[k,:])+' date '+'['+str(currentDayNum_singleSatLine[0])+','+str(currentDayNum_singleSatLine[-1])+']'); #report duplicate sites
                                    checkr_uniq, checkr_uniqCnts = np.unique(currentSecTotal_singleSatLine, return_counts=True); #get the uniques and the counts
                                    jk = np.isin(currentSecTotal_singleSatLine,checkr_uniq[checkr_uniqCnts>1]); #use uniqs and cnts to get where the repeated data are at
                                    jkCnts = checkr_uniqCnts[checkr_uniqCnts>1]; #get the amt of repeats for each
                                    currentsTEC_singleSatLine = np.float64(unfilt_dict['sTEC'][currentSite_loc][currentSat_loc][currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #activate abs illegal quadruple indexing
                                    jkWhere = np.where(jk)[0]; #get where
                                    jkWhere_oop = np.copy(jkWhere); #make where copy for some index work for sTEC
                                    jkWhere_oopSkips = np.where(np.diff(jkWhere_oop) > 1)[0]; #get skips
                                    for op in range(0,jkWhere_oopSkips.size):
                                        jkWhere_oop[jkWhere_oopSkips[op]+1] = jkWhere_oop[jkWhere_oopSkips[op]]+1; #replace value since sTEC is out of phase so diff needs to diff end not start
                                    #END FOR op
                                    if( np.all(np.isclose(np.diff(currentElv_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentLat_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentLong_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentsTEC_singleSatLine[jkWhere_oop])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0)) & \
                                        np.all(np.diff(currentSec_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]] == 0) ):
                                        #in this instance time/lat/long/elv are duplicated and sTEC (col11/col2fromRight) is too - but it is out of phase
                                        #Note vTEC (col6) is calc'd based on sTEC & elev (col8) which are out of phase so it looks unique every time - it's not
                                        #e.g.:
                                        # 327 4875   7  20  30    8.769    3.110   53.223   17.425  -63.475   10.613  302.234 [REAL]
                                        # 328 4875   7  20  30    8.685    3.110   53.223   17.425  -63.475   10.511  302.234 [FAKE]
                                        # 329 4875   7  21   0    8.683    3.119   53.199   17.440  -63.469   10.511  302.641 [REAL]
                                        # 330 4875   7  21   0    8.627    3.119   53.199   17.440  -63.469   10.444  302.641 [FAKE]
                                        checkr_keep = jkWhere[np.cumsum(jkCnts)-jkCnts[0]]; # keep 1st, delete the 2nd/3rd/4th duplicate
                                        checkr_del = jkWhere[~np.isin(jkWhere,checkr_keep)]; #delete index array        
                                    elif( np.all(np.isclose(np.diff(currentElv_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentLat_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentLong_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        ((np.all(np.isclose(np.diff(currentsTEC_singleSatLine[jkWhere_oop])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0)) == False) & \
                                        np.all(np.isclose(np.diff(currentsTEC_singleSatLine[jkWhere_oop])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0) | ~np.isclose(np.diff(currentDayNum_singleSatLine[jkWhere_oop])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0))) & \
                                        np.all(np.diff(currentSec_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]] == 0) ):
                                        #for the case where the sTEC catch doesn't quite work as hoped - occured during a day transition
                                        # 859 4875  23  58  30    5.384   19.449   66.177  -36.272  -67.897    5.821  112.627 [FAKE]
                                        # 860 4875  23  59   0    5.382   19.458   66.113  -36.261  -67.886    5.821  112.033 [REAL]
                                        # 861 4875  23  59   0    5.413   19.458   66.113  -36.261  -67.886    5.854  112.033 [FAKE] (technically could make another data pt at 23 59 30 w/ this sTEC but would need to properly propagate sat's lat/long/elv to shake out vTEC)
                                        # 3 4876   0   0  30    5.127   -4.513   65.331  -36.266  -67.825    5.575  111.455 [REAL]
                                        # 4 4876   0   0  30    5.155   -4.513   65.331  -36.266  -67.825    5.605  111.455 [FAKE]
                                        # 5 4876   0   1   0    5.151   -4.504   65.217  -36.257  -67.811    5.605  110.975 [REAL]
                                        checkr_keep = jkWhere[np.cumsum(jkCnts)-jkCnts[0]]; # keep 1st, delete the 2nd/3rd/4th duplicate
                                        checkr_del = jkWhere[~np.isin(jkWhere,checkr_keep)]; #delete index array     
                                    elif( np.all(np.isclose(np.diff(currentElv_singleSatLine[jk])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0)) & \
                                          np.all(np.isclose(np.diff(currentLat_singleSatLine[jk])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0)) & \
                                          np.all(np.isclose(np.diff(currentLong_singleSatLine[jk])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0)) & \
                                          np.all(np.isclose(np.diff(currentsTEC_singleSatLine[jkWhere_oop])[np.insert(np.cumsum(jkCnts)[:-2]+1,0,1)], 0)) & \
                                          np.all(np.diff(currentSec_singleSatLine[jk])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1] == 0) ):
                                        #for the case where it's the opposite - real is 2nd value not 1st
                                        # 909 4877   2  21  30   13.488    -.169   29.962   -8.478  -37.906   23.286  173.318 [FAKE] (not incl. in dataset due to elevation)
                                        # 910 4877   2  21  30   13.431    -.169   29.962   -8.478  -37.906   23.188  173.318 [REAL] (not incl. in dataset due to elevation)
                                        # 911 4877   2  21  30   13.431    -.169   29.962   -8.478  -37.906   23.188  173.318 [FAKE] (not incl. in dataset due to elevation)
                                        # 912 4877   2  21  45   13.444    -.164   30.012   -8.468  -37.898   23.188  173.209 [FAKE]
                                        # 913 4877   2  21  45   13.404    -.164   30.012   -8.468  -37.898   23.118  173.209 [REAL] !note the real one isn't 1st!
                                        # 914 4877   2  21  45   13.404    -.164   30.012   -8.468  -37.898   23.118  173.209 [FAKE]
                                        # 915 4877   2  22   0   13.417    -.159   30.062   -8.459  -37.890   23.118  173.100 [FAKE]
                                        # 916 4877   2  22   0   13.366    -.159   30.062   -8.459  -37.890   23.031  173.100 [REAL]
                                        # 917 4877   2  22  15   13.380    -.155   30.112   -8.449  -37.882   23.031  172.991 [FAKE]
                                        # 918 4877   2  22  15   13.320    -.155   30.112   -8.449  -37.882   22.927  172.991 [REAL]
                                        # 919 4877   2  22  15   13.320    -.155   30.112   -8.449  -37.882   22.927  172.991 [FAKE]
                                        # 920 4877   2  22  30   13.333    -.150   30.162   -8.440  -37.874   22.927  172.883 [FAKE]
                                        # 921 4877   2  22  30   13.272    -.150   30.162   -8.440  -37.874   22.823  172.883 [REAL]
                                        # 922 4877   2  22  30   13.272    -.150   30.162   -8.440  -37.874   22.823  172.883 [FAKE]
                                        # 923 4877   2  22  45   13.285    -.145   30.212   -8.430  -37.866   22.823  172.774 [FAKE]
                                        # 924 4877   2  22  45   13.239    -.145   30.212   -8.430  -37.866   22.743  172.774 [REAL]
                                        # 925 4877   2  22  45   13.239    -.145   30.212   -8.430  -37.866   22.743  172.774 [FAKE]
                                        # 926 4877   2  23   0   13.252    -.141   30.262   -8.421  -37.858   22.743  172.665 [FAKE]
                                        # 927 4877   2  23   0   13.193    -.141   30.262   -8.421  -37.858   22.643  172.665 [REAL]
                                        # 928 4877   2  23   0   13.193    -.141   30.262   -8.421  -37.858   22.643  172.665 [FAKE]
                                        # 929 4877   2  23  15   13.206    -.136   30.312   -8.411  -37.850   22.643  172.556 [FAKE]
                                        # 930 4877   2  23  15   13.143    -.136   30.312   -8.411  -37.850   22.534  172.556 [REAL]
                                        # 931 4877   2  23  15   13.143    -.136   30.312   -8.411  -37.850   22.534  172.556 [FAKE]
                                        # 932 4877   2  23  30   13.156    -.131   30.362   -8.402  -37.842   22.534  172.447 [FAKE]
                                        checkr_keep = jkWhere[np.cumsum(jkCnts)-jkCnts[0]+1]; # keep 2nd, delete the 1st/3rd/4th duplicate
                                        checkr_del = jkWhere[~np.isin(jkWhere,checkr_keep)]; #delete index array     
                                    elif( np.all(np.isclose(np.diff(currentElv_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentLat_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        np.all(np.isclose(np.diff(currentLong_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]], 0)) & \
                                        (np.all(np.isclose(np.diff(currentsTEC_singleSatLine[jkWhere_oop])[np.cumsum(jkCnts)[:-1]-jkCnts[0]+1], 0)) == False) & \
                                        np.all(np.diff(currentSec_singleSatLine[jk])[np.cumsum(jkCnts)-jkCnts[0]] == 0) ):
                                        #for the case where sTEC is never duplicated, guess avg the 2 values together?
                                        # 4983 6283  23  57  15   11.039   19.292   53.294  -38.029  -69.927   13.351   96.955 [REAL, no duplicate]
                                        # 4987 6283  23  57  30   11.043   19.297   53.220  -38.024  -69.919   13.366   96.796 [keep, for now]
                                        # 4988 6283  23  57  30   11.035   19.297   53.220  -38.024  -69.919   13.356   96.796
                                        # 4993 6283  23  57  45   11.045   19.302   53.146  -38.018  -69.911   13.379   96.636 [keep, for now]
                                        # 4994 6283  23  57  45   11.037   19.302   53.146  -38.018  -69.911   13.370   96.636
                                        # 4998 6283  23  58   0   11.043   19.306   53.073  -38.013  -69.903   13.387   96.477 [keep, for now]
                                        # 4999 6283  23  58   0   11.050   19.306   53.073  -38.013  -69.903   13.396   96.477
                                        # 5004 6283  23  58  15   11.039   19.311   52.999  -38.007  -69.895   13.393   96.319 [keep, for now]
                                        # 5005 6283  23  58  15   11.053   19.311   52.999  -38.007  -69.895   13.410   96.319
                                        # 5008 6283  23  58  30   11.025   19.316   52.925  -38.001  -69.887   13.387   96.161 [REAL, no duplicate]
                                        checkr_keep = jkWhere[np.cumsum(jkCnts)-jkCnts[0]]; # keep 1st, delete the 2nd/3rd/4th duplicate
                                        #for now just keep 1st, below code includes vTEC mapping from an avg'd sTEC pt but there's an issue w/ altFactor that doesn't match what is used for LISN for the expected alt of 350km (confirmed LISN uses 400 km for mapping due to heritage)
                                        # if( np.isscalar(unfilt_dict['pierceAlt']) == True ):
                                        #     altFactor = (6371+0)/(6371+unfilt_dict['pierceAlt']); #calcs altitude factor assuming receiver height is 0, ionosphere is at pierceAlt, and earth is perfect sphere
                                        # else:
                                        #     altFactor = (6371+0)/(6371+unfilt_dict['pierceAlt'][0]); #calcs altitude factor assuming receiver height is 0, ionosphere is at pierceAlt, and earth is perfect sphere
                                        # #END IF
                                        # for jj in range(0,checkr_keep.size):
                                        #     currentsTEC_singleSatLine[checkr_keep[jj]] =  np.mean(currentsTEC_singleSatLine[checkr_keep[jj]:checkr_del[jj]+1]); #avg the two values together
                                        #     currentvTEC_singleSatLine[checkr_keep[jj]] = currentsTEC_singleSatLine[checkr_keep[jj]]*np.cos(np.arcsin(altFactor*np.cos(currentElv_singleSatLine[checkr_keep[jj]]*np.pi/180))); #TECU, calc new vTEC based off of newly avg'd sTEC
                                        #     #eqvtec(iprn,m) = tec(nk)*cos(asin(0.94092*cos(xelev(nk)*degrad)))
                                        #     #The factor = 0.94092 includes the altitude of 350.0 km (which is the sub-ionospheric point)
                                        #END FOR jj
                                        checkr_del = jkWhere[~np.isin(jkWhere,checkr_keep)]; #delete index array  
                                    else:
                                        print('Duplicate timestamps occured that don\'t have a pattern that has been coded to be fixed (or can be fixed dundun). Code something to fix this pls.');
                                        sys.crash();
                                    #END IF
                                    #Remove whatever is deemed needed to be removed
                                    #---INTS---
                                    currentSat_singleSatLine = np.delete(currentSat_singleSatLine,checkr_del); #remove duplicates that fit the above pattern
                                    currentYear_singleSatLine = np.delete(currentYear_singleSatLine,checkr_del);
                                    currentDayNum_singleSatLine = np.delete(currentDayNum_singleSatLine,checkr_del);
                                    currentHour_singleSatLine = np.delete(currentHour_singleSatLine,checkr_del);
                                    currentMin_singleSatLine = np.delete(currentMin_singleSatLine,checkr_del);
                                    currentSec_singleSatLine = np.delete(currentSec_singleSatLine,checkr_del);
                                    #---FLOATS---
                                    currentLat_singleSatLine = np.delete(currentLat_singleSatLine,checkr_del);
                                    currentLong_singleSatLine = np.delete(currentLong_singleSatLine,checkr_del);
                                    currentElv_singleSatLine = np.delete(currentElv_singleSatLine,checkr_del);
                                    #---STRINGS---
                                    currentSite_singleSatLine = np.delete(currentSite_singleSatLine,checkr_del);
                                    currentSatType_singleSatLine = np.delete(currentSatType_singleSatLine,checkr_del);
                                    #---UNFILT DATA---
                                    currentvTEC_singleSatLine = np.delete(currentvTEC_singleSatLine,checkr_del);
                                    currentsTECerror_singleSatLine = np.delete(currentsTECerror_singleSatLine,checkr_del);
                                    currentSecTotal_singleSatLine = np.delete(currentSecTotal_singleSatLine,checkr_del);
                                    # currentsTEC_singleSatLine = np.delete(currentsTEC_singleSatLine,checkr_del); #[only for debug]
                                    del currentsTEC_singleSatLine; #clean the mem  
                                    if( np.unique(currentSecTotal_singleSatLine).size != currentSecTotal_singleSatLine.size ):
                                        print('ERROR uniques not successfully removed. Failure. Crashing. Code this better pls.');
                                        sys.crash();
                                    #END IF
                                #END IF
                                
                                #-----This code bit exists to deal with a day change b/c the TEC bias offset for a day is diff from one to the next-----
                                if( currentDayNum_singleSatLine[0] != currentDayNum_singleSatLine[-1] ): #if this is true, we are looking at different days' data being used
                                # if( np.any(currentDayNum_singleSatLine != dateRange_dayNum_full[i,1]) ): #if this is true, we are looking at different days' data being used [ponder this implementation]
                                    #the problem with that is that it seems the site can have different TEC fit values for different days
                                    #so the lines are continuous - but have a big jump where the different days are fit to different means(?) I'm not sure
                                    #but there's def a big jump we're gonna fix here
                                    currentDayNum_singleSatLine_firstDayLocations = np.where( currentDayNum_singleSatLine[0] == currentDayNum_singleSatLine )[0]; #get index where 
                                    currentDayNum_singleSatLine_secondDayLocations = np.where( currentDayNum_singleSatLine[-1] == currentDayNum_singleSatLine )[0]; #get index where 
                                    currentTEC_singleSatLine_diff = np.diff(currentvTEC_singleSatLine); #get the deltas between the TEC values
                                    currentTEC_singleSatLine_firstDayLast = currentvTEC_singleSatLine[currentDayNum_singleSatLine_firstDayLocations[-1]]; #get the first day's last TEC value
                                    currentTEC_singleSatLine_secondDay = currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations]; #get the second day's TEC values
                                    if( currentDayNum_singleSatLine_firstDayLocations.size > 1 ): #make sure currentDayNum_singleSatLine_firstDayLocations[-2] exists
                                        currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_firstDayLocations[-2]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                        #sorry about these confusing explanations but it def makes sense to me rightnow
                                    else: #if it doesn't, time to use currentDayNum_singleSatLine_secondDayLocations[0] instead of currentDayNum_singleSatLine_firstDayLocations[-2] (so projecting from reverse instead of forward)
                                        currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_secondDayLocations[0]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                        #sorry about these confusing explanations but it def makes sense to me rightnow
                                    #END IF
                                    currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations] = np.float64(currentTEC_singleSatLine_secondDay); #put em in
                                #END IF

                                if( np.all(currentSecTotal_singleSatLine[:-1] <= currentSecTotal_singleSatLine[1:]) == False ):
                                    print('NOT SORTED AS MONOTONICALLY INCREASING TIME CORRECTLY, CHECK IT OUT');
                                    sys.crash(); #not a real call, actually crashes
                                #END IF
                                
                                #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                
                                #get the time steps right is the first step!
                                #~~~COMPACT THE TIME STAMPS IN THE DATA~~~
                                #basically, it needs to meet the lowest common demoninator (30 sec per satellite reading)
                                #10 sec data will be compacted by averaging 1-10, 11-20, and 21-30 (shown as 10, 20, 30 in the data) into just 30 at the
                                #position of the 30 sec data pt
                                #now it could be also an average of the position datas, but I'm going with the last position for now since I don't know how the position is chosen to start with so the easier way is the way to go!
                                
                                # if( current_TEC_dataRate > 66 ):
                                #     # if( np.all(np.diff(currentSecTotal_singleSatLine)==np.diff(currentSecTotal_singleSatLine)[0]) == False ):
                                #     print('yo check this')
                                #     sys.crash();
                                #     #END IF
                                # #EN DIF
                                                                        
                                #-----CONVERT NEARLY 0/30 TIMESTAMPS INTO 0/30 TIMESTAMPS-----
                                if( np.all(np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps) == False) ): #only run if its not already perfect 0/30 time stamps
                                    #this saves interpolation efforts on 29/30/31 time stamps that vary just alittle bit, w/o this those would be interpolated to 30 while it doesn't really matter if it's 1 sec off
                                    for m in range(0,TEC_dataRate_allowedStamps.size):
                                        allowedStamps_range = np.array( (TEC_dataRate_allowedStamps[m]-TEC_dataRate*TEC_timeTolerance,TEC_dataRate_allowedStamps[m]+TEC_dataRate*TEC_timeTolerance) ); #get the current allowed time stamp tolerance range
                                        if( np.any(allowedStamps_range < 0) ):
                                            #Deals with a -3 to 3 time frame around 0 or something
                                            allowedStamps_range[allowedStamps_range < 0] = allowedStamps_range[allowedStamps_range < 0]+60; #keep in the 0-59 sec range
                                            current_timeInTolerance = ((currentSec_singleSatLine >= allowedStamps_range[0]) | (currentSec_singleSatLine <= allowedStamps_range[1])) & (currentSec_singleSatLine != TEC_dataRate_allowedStamps[m]); #get the times within the time stamp tolerance range (so for 0 sec, 57 to 3 sec is taken as 0 sec)
                                            #---CHECK SECONDS---
                                            if( np.any(currentSec_singleSatLine >= allowedStamps_range[0]) ): #catches a 58 and rolls it over to 0 (+1 minute as well)
                                                kj = currentSec_singleSatLine >= allowedStamps_range[0]; #get where the minute rolls over
                                                currentMin_singleSatLine[kj] += 1; #min, the minutes are incremented as well if the seconds rolled over     
                                                #---CHECK MINUTES---
                                                kj = currentMin_singleSatLine >= 60; #find incorrect time keeping
                                                if( np.sum(kj) > 0 ):
                                                    #increment the main hour variable where minutes are 60 or more
                                                    currentHour_singleSatLine[kj] += 1; #hour, increment time by 1
                                                    currentMin_singleSatLine[kj] = 0; #min, remove 60 time units from the time keeping
                                                    #---CHECK HOURS---
                                                    kj = currentHour_singleSatLine >= 24; #find incorrect time keeping
                                                    if( np.sum(kj) > 0 ):
                                                        #increment the main day variable where hours are 24 or more
                                                        currentDayNum_singleSatLine[kj] += 1; #day number, increment time by 1 [integer version]
                                                        currentHour_singleSatLine[kj] = 0; #hour, remove 24 time units from the time keeping
                                                        #---CHECK DAYS---
                                                        #deal with day limit is based on leap year or not
                                                        dayLim = np.ones(currentDayNum_singleSatLine.shape,dtype=np.int16)*365; #day number, get the day number limits as 365
                                                        #adjust leap years to 366
                                                        leapYears = (np.mod(currentYear_singleSatLine,4) == 0) & (np.mod(currentYear_singleSatLine,100) != 0) & (np.mod(currentYear_singleSatLine,400) == 0)
                                                        dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
                                                        kj = currentDayNum_singleSatLine >= dayLim; #find incorrect time keeping
                                                        if( np.sum(kj) > 0 ):
                                                            #increment the main year variable where day number is equal to the day number limit or higher than it
                                                            currentYear_singleSatLine += 1; #year, increment time by 1
                                                            currentDayNum_singleSatLine[kj] = 1; #day num, remove the day number limit time units from the time keeping
                                                        #END IF
                                                    #END IF
                                                #END IF
                                            #END IF
                                        else:
                                            current_timeInTolerance = (currentSec_singleSatLine >= allowedStamps_range[0]) & (currentSec_singleSatLine <= allowedStamps_range[1]) & (currentSec_singleSatLine != TEC_dataRate_allowedStamps[m]); #get the times within the time stamp tolerance range (so for 30 sec, 27 to 33 sec is taken as 30 sec)
                                        #END IF
                                        currentSec_singleSatLine[current_timeInTolerance] = TEC_dataRate_allowedStamps[m]; #set the times within the tolerance range to the allowed time stamp for that range
                                    #END FOR m
                                    currentSecTotal_singleSatLine = np.int64(currentDayNum_singleSatLine)*86400 + np.int64(currentHour_singleSatLine)*3600 + np.int64(currentMin_singleSatLine)*60 + np.int64(currentSec_singleSatLine); #sec, calculate days/hour/min/sec into total seconds but do it at 64 bit for extra good precision in the seconds [no overflow]
                                    
                                    #-----DEAL WITH NEWLY MINTED REPEATS-----
                                    #!note this is not needed if above code is ditched and only interpolation is used!
                                    if( np.unique(currentSecTotal_singleSatLine).size != currentSecTotal_singleSatLine.size ):
                                        # try:
                                        #     cntr_activater += 1;
                                        # except:
                                        #     cntr_activater = 1; #start
                                        # print('WARNING NEW REPEATS CODE WAS ACTIVATED '+str(cntr_activater));
                                        # print('size start: '+str(currentSecTotal_singleSatLine.size)+' & uniq size: '+str(np.unique(currentSecTotal_singleSatLine).size));
                                        #this is for a weird thing where one had a 0 second time stamp and one had a 1 second time stamp w/ same day/hr/min so it was unique at the check but condensed to be 0/0 after above code. gonna just avg and call it a day
                                        _, current_repeatsIndex, current_repeatsCounts = np.unique(currentSecTotal_singleSatLine, return_index=True, return_counts=True); #get indexes and counts
                                        current_repeatsCounts = np.where(current_repeatsCounts >= 2)[0]; #get where > 2
                                        for kl in range(0,current_repeatsCounts.size): #to support more than 1 > 2 at once (worst case)
                                            current_repeatsValue = currentSecTotal_singleSatLine[current_repeatsIndex[current_repeatsCounts[kl]]]; #get the repeated value
                                            track_repeats = currentSecTotal_singleSatLine == current_repeatsValue; #get where the repeat value is
                                            track_repeatsFirst = np.where(track_repeats)[0][0]; #get the 1st one to edit
                                            
                                            #-----AVG DUPLICATE DATA INTO 1ST INSTANCE-----
                                            #---FLOATS---
                                            currentLat_singleSatLine[track_repeatsFirst] = np.mean(currentLat_singleSatLine[track_repeats]); #keep the good stuff
                                            currentLong_singleSatLine[track_repeatsFirst] = np.mean(currentLong_singleSatLine[track_repeats]); #keep the good stuff
                                            currentElv_singleSatLine[track_repeatsFirst] = np.mean(currentElv_singleSatLine[track_repeats]); #keep the good stuff
                                            #---UNFILT DATA---
                                            currentvTEC_singleSatLine[track_repeatsFirst] = np.mean(currentvTEC_singleSatLine[track_repeats]); #keep the good stuff
                                            currentsTECerror_singleSatLine[track_repeatsFirst] = np.mean(currentsTECerror_singleSatLine[track_repeats]); #keep the good stuff
                                            currentSecTotal_singleSatLine[track_repeatsFirst] = np.mean(currentSecTotal_singleSatLine[track_repeats]); #keep the good stuff
                                            
                                            #-----REMOVE DUPLICATE TIME STAMP DATA-----
                                            track_repeats[track_repeatsFirst] = False; #set 1st to false so don't delete it
                                            track_repeats = ~track_repeats; #logical flip for keeping good stuff
                                            #---INTS---
                                            currentSat_singleSatLine = currentSat_singleSatLine[track_repeats]; #keep the good stuff
                                            currentYear_singleSatLine = currentYear_singleSatLine[track_repeats]; #keep the good stuff
                                            currentDayNum_singleSatLine = currentDayNum_singleSatLine[track_repeats]; #keep the good stuff
                                            currentHour_singleSatLine = currentHour_singleSatLine[track_repeats]; #keep the good stuff
                                            currentMin_singleSatLine = currentMin_singleSatLine[track_repeats]; #keep the good stuff
                                            currentSec_singleSatLine = currentSec_singleSatLine[track_repeats]; #keep the good stuff
                                            #---FLOATS---
                                            # currentDayNumF_singleSatLine = currentDayNumF_singleSatLine[track_repeats]; #keep the good stuff
                                            currentLat_singleSatLine = currentLat_singleSatLine[track_repeats]; #keep the good stuff
                                            currentLong_singleSatLine = currentLong_singleSatLine[track_repeats]; #keep the good stuff
                                            currentElv_singleSatLine = currentElv_singleSatLine[track_repeats]; #keep the good stuff
                                            #---STRINGS---
                                            currentSite_singleSatLine = currentSite_singleSatLine[track_repeats]; #keep the good stuff
                                            currentSatType_singleSatLine = currentSatType_singleSatLine[track_repeats]; #keep the good stuff
                                            #---UNFILT DATA---
                                            currentvTEC_singleSatLine = currentvTEC_singleSatLine[track_repeats]; #keep the good stuff
                                            currentsTECerror_singleSatLine = currentsTECerror_singleSatLine[track_repeats]; #keep the good stuff
                                            currentSecTotal_singleSatLine = currentSecTotal_singleSatLine[track_repeats]; #keep the good stuff
                                        #END FOR kl
                                        # print('size end: '+str(currentSecTotal_singleSatLine.size)+' & uniq size: '+str(np.unique(currentSecTotal_singleSatLine).size));
                                    #END IF
                                #END IF
                                
                                #-----IDENTIFY TIME GAPS-----
                                current_timeDiff = np.diff(currentSecTotal_singleSatLine); #get the diff between the times
                                if( current_TEC_dataRate <= TEC_dataRate*(1+TEC_timeTolerance) ):
                                    track_gaps = np.where( current_timeDiff > TEC_dataRate*2+TEC_dataRate*(TEC_timeTolerance) )[0]; #catch time gaps that can't be interpolated (missing 1 is reliably interpolated) [sub-30 second time gaps aren't an issue]
                                    track_gaps_fill = np.where( current_timeDiff > TEC_dataRate+TEC_dataRate*(TEC_timeTolerance) )[0]; #catch time gaps that CAN be interpolated (missing 1 is reliably interpolated) [sub-30 second time gaps aren't an issue]
                                else:
                                    track_gaps = np.where( current_timeDiff > current_TEC_dataRate*(1+TEC_timeTolerance) )[0]; #catch time gaps [for 60 second time rates or more, can't interpolate time gaps b/c there's already implied time gaps]
                                    track_gaps_fill = track_gaps; #set same as track_gaps b/c the current_TEC_dataRate is higher than the desired data rate - no chance for fillable gaps
                                #END IF
                                track_gapsTimes = np.empty((track_gaps.size,2),dtype=np.int64); #preallocate
                                for kj in range(0,track_gaps.size):
                                    track_gapsTimes[kj,:] = currentSecTotal_singleSatLine[track_gaps[kj]:track_gaps[kj]+2]; #fill in the gap times, will use later to enforce gap times
                                #END FOR kj
                                
                                #-----INTERPOLATE TO REQUIRED DATA RATE TIME STAMPS (only needed if gaps or not already on 0/30 allowed time stamps)-----
                                if( np.any(np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps) == False) | (track_gaps_fill.size > 0) | (current_TEC_dataRate > TEC_dataRate*(1+TEC_timeTolerance)) ):
                                    #---GET INTERP'D TIMES---
                                    if( np.any(np.isin(TEC_dataRate_allowedStamps,currentSec_singleSatLine[0])) ):
                                        ideal_timeStart = currentSecTotal_singleSatLine[0]; #set the ideal start time to the real value
                                    else:
                                        ideal_timeStartDelta = TEC_dataRate_allowedStamps60_sameOrder-currentSec_singleSatLine[0]; #get the deltas possible
                                        ideal_timeStartDelta = ideal_timeStartDelta[ideal_timeStartDelta > 0]; #get the positive deltas
                                        ideal_timeStartDelta = ideal_timeStartDelta[np.min(ideal_timeStartDelta) == ideal_timeStartDelta].item(); #get the delta to use
                                        ideal_timeStart = currentSecTotal_singleSatLine[0] + ideal_timeStartDelta; #tack on the required delta to move up to the next 0/30 allowed time step
                                    #END IF
                                    if( np.any(np.isin(TEC_dataRate_allowedStamps,currentSec_singleSatLine[-1])) ):
                                        ideal_timeEnd = currentSecTotal_singleSatLine[-1]; #set the ideal start time to the real value
                                    else:
                                        ideal_timeEndDelta = TEC_dataRate_allowedStamps60_sameOrder-currentSec_singleSatLine[-1]; #get the deltas possible
                                        ideal_timeEndDelta = ideal_timeEndDelta[ideal_timeEndDelta > 0]; #get the positive deltas
                                        ideal_timeEndDelta = ideal_timeEndDelta[np.min(ideal_timeEndDelta) == ideal_timeEndDelta].item(); #get the delta to use
                                        ideal_timeEnd = currentSecTotal_singleSatLine[-1] + ideal_timeEndDelta; #tack on the required delta to move up to the next 0/30 allowed time step
                                    #END IF
                                    ideal_time = np.arange(ideal_timeStart,ideal_timeEnd+TEC_dataRate,TEC_dataRate); #sec, get the ideal time
                                    # print('ideal time size: '+str(ideal_time.size)+' & orig time size: '+str(currentSecTotal_singleSatLine.size)+' & detected data rate: '+str(current_TEC_dataRate))
                                    real_times = np.isin(ideal_time,currentSecTotal_singleSatLine); #get a logical array of the real times
                                    real_timesRev = np.isin(currentSecTotal_singleSatLine,ideal_time); #get a logical array of the real times
                                    #---GET TEC INTERPER---
                                    TEC_interper = interp1d(currentSecTotal_singleSatLine,currentvTEC_singleSatLine,kind='linear',fill_value='extrapolate'); #make an interpolator
                                    #---INTERP TEC---
                                    ideal_TEC = TEC_interper(ideal_time); #make TEC for all the times
                                    ideal_TEC[real_times] = currentvTEC_singleSatLine[real_timesRev]; #make sure the real times are the real data 
                                    currentvTEC_singleSatLine = np.copy(ideal_TEC); #overwrite-filter later
                                    #FILTER LATER TO MAKE CODE JIVE
                                    #---INTERP TEC ERR(???)---
                                    #!!need to deal with mapping sTEC->vTEC on error, vTEC->delta-vTEC on error!!
                                    #note I'm not super sure how to handle TEC error, this is just to get the ball rolling (figure out good way later)
                                    TECerr_interper = interp1d(currentSecTotal_singleSatLine,currentsTECerror_singleSatLine,kind='linear',fill_value='extrapolate'); #make an interpolator
                                    ideal_TECerr = TECerr_interper(ideal_time); #make TEC for all the times
                                    ideal_TECerr[real_times] = currentsTECerror_singleSatLine[real_timesRev]; #make sure the real times are the real data
                                    currentsTECerror_singleSatLine = np.copy(ideal_TECerr); #for now keep it
                                    #---INTERP LAT---
                                    #!!THIS MAY BE A BAD FIT!!
                                    lat_interper = interp1d(currentSecTotal_singleSatLine,currentLat_singleSatLine,kind='linear',fill_value='extrapolate'); #make an interpolator
                                    ideal_lat = lat_interper(ideal_time); #make lat for all the times
                                    ideal_lat[real_times] = currentLat_singleSatLine[real_timesRev]; #make sure the real times are the real data
                                    currentLat_singleSatLine = np.copy(ideal_lat); #keep it now
                                    #---INTERP LONG---
                                    #!!THIS MAY BE A BAD FIT!!
                                    long_interper = interp1d(currentSecTotal_singleSatLine,currentLong_singleSatLine,kind='linear',fill_value='extrapolate'); #make an interpolator
                                    ideal_long = long_interper(ideal_time); #make long for all the times
                                    ideal_long[real_times] = currentLong_singleSatLine[real_timesRev]; #make sure the real times are the real data
                                    currentLong_singleSatLine = np.copy(ideal_long); #keep it now
                                    #---INTERP ELV---
                                    #!!THIS MAY BE A BAD FIT!!
                                    elv_interper = interp1d(currentSecTotal_singleSatLine,currentElv_singleSatLine,kind='linear',fill_value='extrapolate'); #make an interpolator
                                    ideal_elv = elv_interper(ideal_time); #make lat for all the times
                                    ideal_elv[real_times] = currentElv_singleSatLine[real_timesRev]; #make sure the real times are the real data
                                    currentElv_singleSatLine = np.copy(ideal_elv); #keep it now
                                    #---CALC DATA THAT IS NOT CONSTANT---
                                    currentDayNum_singleSatLine = np.int16(ideal_time//86400); #get day number
                                    midTime = np.mod(ideal_time,86400); #get mid time w/o days
                                    currentHour_singleSatLine = np.int16(midTime//3600); #get hour
                                    midTime = np.mod(midTime,3600); #get mid time w/o days
                                    currentMin_singleSatLine = np.int16(midTime//60); #get min
                                    midTime = np.mod(midTime,60); #get mid time w/o days
                                    currentSec_singleSatLine = np.int16(midTime); #get sec
                                    #old way that replied on floating point error for subsequent calcs
                                    # midTime = ideal_time - np.int64(currentDayNum_singleSatLine)*86400; #get mid time w/o days
                                    # currentHour_singleSatLine = np.int16(midTime//3600); #get hour
                                    # midTime = midTime - np.int64(currentHour_singleSatLine)*3600; #get mid time w/o days
                                    # currentMin_singleSatLine = np.int16(midTime//60); #get min
                                    # midTime = midTime - np.int64(currentMin_singleSatLine)*60; #get mid time w/o days
                                    # currentSec_singleSatLine = np.int16(midTime); #get sec
                                    #Calc day num float w/ all the recovered times above
                                    # currentDayNumF_singleSatLine = np.float64(currentDayNum_singleSatLine) + currentHour_singleSatLine/24 + currentMin_singleSatLine/1440 + currentSec_singleSatLine/86400; #days, calculate hour/min/sec into days and add to the current day but do it at 64 bit for extra good precision in the seconds
                                    #Calc sec int w/ all the recovered times above
                                    # currentSecTotal_singleSatLine = np.int64(currentDayNum_singleSatLine)*86400 + np.int64(currentHour_singleSatLine)*3600 + np.int64(currentMin_singleSatLine)*60 + np.int64(currentSec_singleSatLine); #sec, calculate hour/min/sec into sec and add to the current day in sec (no overflow)
                                    currentSecTotal_singleSatLine = np.copy(ideal_time); #it's the same so just copy over
                                    #---REPEAT DATA THAT IS CONSTANT---
                                    # if( np.all(currentSat_singleSatLine == currentSat_singleSatLine[0]) == False ):
                                    #     print('yo bad news'); #error check
                                    #     sys.crash()
                                    # #END IF
                                    # if( np.all(currentYear_singleSatLine == currentYear_singleSatLine[0]) == False ):
                                    #     print('yo bad news'); #error check
                                    #     sys.crash()
                                    # #END IF
                                    # if( np.all(currentSite_singleSatLine == currentSite_singleSatLine[0]) == False ):
                                    #     print('yo bad news'); #error check
                                    #     sys.crash()
                                    # #END IF
                                    currentSat_singleSatLine = np.tile(currentSat_singleSatLine[0],(ideal_time.shape)); #copy sat ID (doesn't change)
                                    currentYear_singleSatLine  = np.tile(currentYear_singleSatLine[0],(ideal_time.shape)); #copy year (doesn't change FOR NOW)
                                    currentSite_singleSatLine = np.tile(currentSite_singleSatLine[0],(ideal_time.shape)); #copy site name (doesn't change)
                                    currentSatType_singleSatLine = np.tile(currentSatType_singleSatLine[0],(ideal_time.shape)); #copy site name (doesn't change)
                                #END IF
                                
                                #-----CALC delta-vTEC FOR "STANDARDIZED" TIME STEP DATA-----
                                currentPolyYvals = savgol_filter(currentvTEC_singleSatLine, windowLen_savGol, order_savGol ); #filter it up
                                current_deltavTEC = currentvTEC_singleSatLine - currentPolyYvals;
                                #recalc here so don't have to do it twice
                                # currentDayNumF_singleSatLine = np.float64(currentDayNum_singleSatLine) + currentHour_singleSatLine/24 + currentMin_singleSatLine/1440 + currentSec_singleSatLine/86400; #days, calculate hour/min/sec into days and add to the current day but do it at 64 bit for extra good precision in the seconds
                                
                                #This is for debug viewing/methods paper
                                # from importlib import reload
                                # import GRITI_import_TEC_Madrigal_viewSatellitePass
                                # reload(GRITI_import_TEC_Madrigal_viewSatellitePass)
                                # from GRITI_import_TEC_Madrigal_viewSatellitePass import GRITI_import_TEC_Madrigal_viewSatellitePass
                                # GRITI_import_TEC_Madrigal_viewSatellitePass(currentsTEC_singleSatLine, currentvTEC_singleSatLine, currentPolyYvals, current_deltavTEC, \
                                #     currentElv_singleSatLine, currentSat_singleSatLine, currentSite_singleSatLine, currentYear_singleSatLine, currentDayNum_singleSatLine, \
                                #     currentHour_singleSatLine, currentMin_singleSatLine, currentSec_singleSatLine, currentLat_singleSatLine, currentLong_singleSatLine, \
                                #     0.5, minElevation, folder, savedVersion = True); #savedVersion also saves a nicer looking version of the plot
                                
                                #-----Shore up time gaps for greater-than-desired-data-rate-----
                                # if( current_TEC_dataRate/2 > TEC_dataRate*(1+TEC_timeTolerance) ):
                                #     #this logic is to deal with a 120 sec data rate where every 60 seconds can be interpolated - but every 30 sec were just interpolated.. gotta figure it out!
                                #     interpolatableTime = np.arange(ideal_timeStart,ideal_timeEnd+current_TEC_dataRate/2,current_TEC_dataRate/2); #sec, get overloaded gaps
                                #     track_gaps = np.append(track_gaps,interpolatableTime.size); #tack on the size
                                #     track_gapsTimes = np.append(track_gapsTimes,np.repeat(interpolatableTime,2).reshape(interpolatableTime.size//2,2)); #tack on the gaps that need to be removed due to too slow a data rate (e.g. we can only keep 1/2 the current_TEC_dataRate through interpolation, 1/4 is too far for interpolation)
                                # #END IF
                                
                                #-----ENFORCE TIME GAPS-----
                                if( track_gaps.size > 0 ):
                                    #-----CREATE KEEP MASK-----
                                    track_matches = np.ones(ideal_time.size,dtype=np.bool_); #preallocate true
                                    for kj in range(0,track_gaps.size):
                                        track_matches[ (ideal_time > track_gapsTimes[kj,0]) & (ideal_time < track_gapsTimes[kj,1]) ] = False; #set to false to be deleted
                                    #END kj
                                    # print('gap size pre: '+str(track_matches.size)+' & size post: '+str(track_matches.sum()));
                                    # print('max gap time: '+str(np.max(np.abs(np.diff(track_gapsTimes,axis=1)))));
                                    
                                    #-----REMOVE UNDEEDED DATA-----
                                    #---INTS---
                                    currentSat_singleSatLine = currentSat_singleSatLine[track_matches]; #keep the good stuff
                                    currentYear_singleSatLine = currentYear_singleSatLine[track_matches]; #keep the good stuff
                                    currentDayNum_singleSatLine = currentDayNum_singleSatLine[track_matches]; #keep the good stuff
                                    currentHour_singleSatLine = currentHour_singleSatLine[track_matches]; #keep the good stuff
                                    currentMin_singleSatLine = currentMin_singleSatLine[track_matches]; #keep the good stuff
                                    currentSec_singleSatLine = currentSec_singleSatLine[track_matches]; #keep the good stuff
                                    #---FLOATS---
                                    # currentDayNumF_singleSatLine = currentDayNumF_singleSatLine[track_matches]; #keep the good stuff
                                    currentLat_singleSatLine = currentLat_singleSatLine[track_matches]; #keep the good stuff
                                    currentLong_singleSatLine = currentLong_singleSatLine[track_matches]; #keep the good stuff
                                    currentElv_singleSatLine = currentElv_singleSatLine[track_matches]; #keep the good stuff
                                    #---STRINGS---
                                    currentSite_singleSatLine = currentSite_singleSatLine[track_matches]; #keep the good stuff
                                    currentSatType_singleSatLine = currentSatType_singleSatLine[track_matches]; #keep the good stuff
                                    #---UNFILT DATA---
                                    currentvTEC_singleSatLine = currentvTEC_singleSatLine[track_matches]; #keep the good stuff
                                    currentsTECerror_singleSatLine = currentsTECerror_singleSatLine[track_matches]; #keep the good stuff
                                    currentSecTotal_singleSatLine = currentSecTotal_singleSatLine[track_matches]; #keep the good stuff
                                    #---FILT DATA---
                                    current_deltavTEC = current_deltavTEC[track_matches]; #keep the good stuff
                                #END IF
                                
                                #-----RECORD DATA IN THE BIG LISTS-----
                                filt_dict['sat'].append(np.int16(currentSat_singleSatLine)); #tack that data on
                                filt_dict['year'].append(np.int16(currentYear_singleSatLine)); #tack that data on
                                filt_dict['dayNum'].append(np.int16(currentDayNum_singleSatLine)); #tack that data on
                                filt_dict['hour'].append(np.int16(currentHour_singleSatLine)); #tack that data on
                                filt_dict['min'].append(np.int16(currentMin_singleSatLine)); #tack that data on
                                filt_dict['sec'].append(np.int16(currentSec_singleSatLine)); #tack that data on
                                filt_dict['time'].append(np.int32(currentSecTotal_singleSatLine)); #tack that data on
                                filt_dict['lat'].append(dataType_meth(currentLat_singleSatLine)); #tack that data on
                                filt_dict['long'].append(dataType_meth(currentLong_singleSatLine)); #tack that data on
                                filt_dict['elev'].append(dataType_meth(currentElv_singleSatLine)); #tack that data on
                                filt_dict['dTEC'].append(dataType_meth(current_deltavTEC)); #tack that data on
                                filt_dict['dTECerr'].append(dataType_meth(currentsTECerror_singleSatLine)); #tack that data on
                                filt_dict['site'].append(currentSite_singleSatLine); #tack that data on
                                filt_dict['satType'].append(currentSatType_singleSatLine); #tack that data on
                            #END IF
                        #END FOR l
                    #END FOR k
                    if( (FLG_parallel == False) & ((np.mod(j+1,100) == 0) | (j == 0) | (j == unfilt_site_unique.size-1)) ): #report only every 100 sites, makes it go quicker but still often enough to keep you up to the minute updated
                        toc = time() - tic2; #for time testing
                        sys.stdout.write('\rSite {} of {} | {} min so far, {} sec per site, ETA {} min\t\t\t\t'.format(j+1,unfilt_site_unique.size,textNice(np.round(toc/60,2)),textNice(np.round(toc/(j+1),4)),textNice(np.round(toc/(j+1)*(unfilt_site_unique.size-j-1)/60,2)))); #report
                        sys.stdout.flush();   
                    #END IF
                #END FOR j
            else: #with dataAggregation on, time to try to include other data to make longer satellite tracks
                #cruise through every site (from the unique of unfilt_dict['site'])
#                    cntr = 0; #for testing
                for j in range(0,unfilt_site_unique.size): 
                    
    #                currentSite_loc = np.where(j == unfilt_site_unique_indexes)[0]; #get the data locations for 1 site
                    currentSite_loc = unfilt_site_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                                    
                    currentSat_unique = np.unique(unfilt_dict['sat'][currentSite_loc]); #get the unique sats in at that current site
                    #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])
                    
#                    currentSite_latRange = np.array( (np.min(TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat]), np.max(TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat]) ) ); #get the lat range
#                    currentSite_longRange = np.array( (np.min(TEC_fileData_float[currentSite_loc,locFloatUnfilt_long]), np.max(TEC_fileData_float[currentSite_loc,locFloatUnfilt_long]) ) ); #get the long range
#                    
#                    currentSite_nearbyLat = TEC_fileData_float[ (TEC_fileData_float[:,locFloatUnfilt_lat] >= currentSite_latRange[0]) & (TEC_fileData_float[:,locFloatUnfilt_lat] <= currentSite_latRange[0])  ,locFloatUnfilt_lat]
                    
                    for k in range(0,currentSat_unique.size):
                        currentSat_loc = np.where( currentSat_unique[k] == unfilt_dict['sat'][currentSite_loc] )[0]; #get the data locations for 1 sat at that 1 site
                        currentvTEC = unfilt_dict['vTEC'][currentSite_loc][currentSat_loc]; #TECU, get the vTEC data for 1 sat at that 1 site - uses abs illegal double indexing
    #                        currentsTECerror = TEC_fileData_float[currentSite_loc,locFloatUnfilt_sTECerr][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentDayNumF = np.float64(unfilt_dict['dayNum'][currentSite_loc][currentSat_loc]) + np.float64(unfilt_dict['time'][currentSite_loc][currentSat_loc])/86400; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
#                        currentElv = TEC_fileData_float[currentSite_loc,locFloatUnfilt_elev][currentSat_loc];
                        
                        currentTimeSplits_loc = np.append( np.insert( np.where(np.diff(currentDayNumF)*86400 > minimumTimeGap)[0]+1 , 0, 0), currentDayNumF.shape ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                        
                        currentLat = unfilt_dict['lat'][currentSite_loc][currentSat_loc]; #degc, get the lat data for 1 sat at that 1 site - uses abs illegal double indexing
                        currentLong = unfilt_dict['long'][currentSite_loc][currentSat_loc]; #degc, get the TEC long for 1 sat at that 1 site - uses abs illegal double indexing
                        
                        for l in range(0,len(currentTimeSplits_loc)-1):
                            currentDayNumF_singleSatLine = currentDayNumF[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                            currentDayNum_singleSatLine_min = (currentDayNumF_singleSatLine-np.round(np.mean(currentDayNumF_singleSatLine)))*144; #min/10, convert to minutes and /10 so it's more stable when fitting a curve
                            currentvTEC_singleSatLine = currentvTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability

                            currentLat_singleSatLine = currentLat[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                            currentLong_singleSatLine = currentLong[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                            
#                            currentElv_singleSatLine = currentElv[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]];

#                            plt.figure();
#                            plt.scatter( currentDayNumF_singleSatLine , currentLat_singleSatLine );
#                            plt.title("Latitude Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                            plt.xlabel('Time [days]');
#                            plt.ylabel('Lat [degc]');
#                            plt.show();
#                            
#                            plt.figure();
#                            plt.scatter( currentDayNumF_singleSatLine , currentLong_singleSatLine );
#                            plt.title("Longitude Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                            plt.xlabel('Time [days]');
#                            plt.ylabel('Long [degc]');
#                            plt.show();
                            
#                            plt.figure();
#                            plt.scatter( currentLong_singleSatLine , currentLat_singleSatLine );
#                            plt.title("Lat/Long Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                            plt.xlabel('Long [degc]');
#                            plt.ylabel('Lat [degc]');
#                            plt.show();
                            
#                            plt.figure();
#                            plt.scatter( currentDayNumF_singleSatLine , currentElv_singleSatLine );
#                            plt.title("Elevation Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                            plt.xlabel('Time [days]');
#                            plt.ylabel('Elevation [deg]');
#                            plt.show();
                                                                                
                            currentDataAgg_stepsFwdBack= np.int64(TEC_dataAgg_timeAdditionLimit/TEC_dataRate); #get the number of steps to take forward and back in time
                            currentDataAgg_dayNumSSLBefore = TEC_timeUnique[ np.where( np.min(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]-currentDataAgg_stepsFwdBack:np.where( np.min(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]  ]; #days, get the times before the observation period
                            currentDataAgg_dayNumSSLAfter = TEC_timeUnique[ (np.where( np.max(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]+1):np.where( np.max(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]+currentDataAgg_stepsFwdBack+1  ]; #days, get the times before the observation period
                            
                            #only do this if any of the days are in the correct day
                            if( (np.any(np.int16(currentDayNumF[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                
                                currentDataAgg_dayNumSSLBefore_min = (currentDataAgg_dayNumSSLBefore - np.round(np.mean(currentDayNumF_singleSatLine)))*144; #min/10, get the times before the observation period
                                currentDataAgg_dayNumSSLAfter_min = (currentDataAgg_dayNumSSLAfter - np.round(np.mean(currentDayNumF_singleSatLine)))*144; #min/10, get the times after the observation period
    
                                currentDataAgg_coefs = np.polynomial.polynomial.polyfit(currentDayNum_singleSatLine_min, currentLat_singleSatLine, 3); #gets the coefficients that fit an x^2 polynomial - what a function to call am I right
                                currentDataAgg_latGuessBefore = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLBefore_min, currentDataAgg_coefs); #calc the guessed values 
                                currentDataAgg_latGuessAfter = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLAfter_min, currentDataAgg_coefs); #calc the guessed values 
#                                    currentDataAgg_latGuess = np.hstack( (currentDataAgg_latGuessBefore, currentLat_singleSatLine , currentDataAgg_latGuessAfter) ); #degc, combine lat guess into one variable
                                                                                                 
                                currentDataAgg_coefs = np.polynomial.polynomial.polyfit(currentDayNum_singleSatLine_min, currentLong_singleSatLine, 3); #gets the coefficients that fit an x^2 polynomial - what a function to call am I right
                                currentDataAgg_longGuessBefore = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLBefore_min, currentDataAgg_coefs); #calc the guessed values 
                                currentDataAgg_longGuessAfter = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLAfter_min, currentDataAgg_coefs); #calc the guessed values 
#                                    currentDataAgg_longGuess = np.hstack( (currentDataAgg_longGuessBefore, currentLong_singleSatLine , currentDataAgg_longGuessAfter) ); #degc, combine lat guess into one variable
                                
#                                plt.figure();
#                                plt.scatter( currentDataAgg_longGuess , currentDataAgg_latGuess );
#                                plt.title("Predicted Lat/Long Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                                plt.xlabel('Long [degc]');
#                                plt.ylabel('Lat [degc]');
#                                plt.show();
                                
                                currentDataAgg_TECBefore = np.ones(currentDataAgg_stepsFwdBack , dtype=dataType_meth)*np.nan; #preallocate
                                cntrSkip = 0; #start up a counter
                                cntrTime = 0; #start up another counter
                                #run backwards
                                currentDataAgg_dayNumSSLBeforeRev = np.flip(currentDataAgg_dayNumSSLBefore); #flip it to make it easy
                                currentDataAgg_dayNumSSLBeforeRevStartIndex = np.where(TEC_timeUnique == currentDataAgg_dayNumSSLBeforeRev[0])[0][0]; #where to start it off
                                while( (cntrSkip < np.int64(minimumTimeGap/TEC_dataRate)) & (cntrTime < currentDataAgg_stepsFwdBack) ):
                                    
                                    currentDataAgg_currentTimeIndexes = TEC_timeUnique_currentTimeArray[currentDataAgg_dayNumSSLBeforeRevStartIndex+cntrTime]; #increment it up as we go to save processing
                                    
                                    currentDataAgg_closeIndexes = np.where( TEC_dataAgg_distToPts_degcSq >= ((unfilt_dict['lat'][currentDataAgg_currentTimeIndexes] - currentDataAgg_latGuessBefore[cntrTime])**2 + (unfilt_dict['long'][currentDataAgg_currentTimeIndexes] - currentDataAgg_longGuessBefore[cntrTime])**2) )[0];
                                    
                                    if( currentDataAgg_closeIndexes.size != 0 ): #record some data
                                        currentDataAgg_TECBefore[cntrTime] = np.mean(unfilt_dict['vTEC'][currentDataAgg_currentTimeIndexes[currentDataAgg_closeIndexes]]); #get the mean of the vTEC around the point we guessed
                                        cntrSkip = 0; #reset
                                    else: #else, note data skipped
                                        cntrSkip = cntrSkip + 1; #increment skip counter
                                    #END IF
                                    cntrTime = cntrTime + 1; #increment total time counter
                                #END WHILE
                                currentDataAgg_TECBefore = np.flip(currentDataAgg_TECBefore); #reverse back to correct time sequence
                                
                                currentDataAgg_TECAfter = np.ones(currentDataAgg_stepsFwdBack , dtype=dataType_meth)*np.nan; #preallocate
                                cntrSkip = 0; #start up a counter
                                cntrTime = 0; #start up another counter
                                #run forwards
                                currentDataAgg_dayNumSSLAfterStartIndex = np.where(TEC_timeUnique == currentDataAgg_dayNumSSLAfter[0])[0][0]; #where to start it off
                                while( (cntrSkip < np.int64(minimumTimeGap/TEC_dataRate)) & (cntrTime < currentDataAgg_stepsFwdBack) ):
                                    
                                    currentDataAgg_currentTimeIndexes = TEC_timeUnique_currentTimeArray[currentDataAgg_dayNumSSLAfterStartIndex + cntrTime]; #increment it up as we go to save processing
                                    
                                    currentDataAgg_closeIndexes = np.where( TEC_dataAgg_distToPts_degcSq >= ((unfilt_dict['lat'][currentDataAgg_currentTimeIndexes] - currentDataAgg_latGuessAfter[cntrTime])**2 + (unfilt_dict['long'][currentDataAgg_currentTimeIndexes] - currentDataAgg_longGuessAfter[cntrTime])**2) )[0];
                                    
                                    if( currentDataAgg_closeIndexes.size != 0 ): #record some data
                                        currentDataAgg_TECAfter[cntrTime] = np.mean(unfilt_dict['vTEC'][currentDataAgg_currentTimeIndexes[currentDataAgg_closeIndexes]]); #get the mean of the vTEC around the point we guessed
                                        cntrSkip = 0; #reset
                                    else: #else, note data skipped
                                        cntrSkip = cntrSkip + 1; #increment skip counter
                                    #END IF
                                    cntrTime = cntrTime + 1; #increment total time counter
                                #END WHILE
                                
                                currentDataAgg_nansBefore = np.logical_not(np.isnan(currentDataAgg_TECBefore)); #find the NaNs
                                currentDataAgg_dayNumSSLBefore = currentDataAgg_dayNumSSLBefore[currentDataAgg_nansBefore]; #delete the NaNs
                                currentDataAgg_TECBefore = currentDataAgg_TECBefore[currentDataAgg_nansBefore]; #delete the NaNs
                                                                
                                currentDataAgg_nansAfter = np.logical_not(np.isnan(currentDataAgg_TECAfter)); #find the NaNs
                                currentDataAgg_dayNumSSLAfter = currentDataAgg_dayNumSSLAfter[currentDataAgg_nansAfter]; #delete the NaNs
                                currentDataAgg_TECAfter = currentDataAgg_TECAfter[currentDataAgg_nansAfter]; #delete the NaNs
                                
                                currentDataAgg_numberBefore = currentDataAgg_dayNumSSLBefore.size; #the size
                                currentDataAgg_numberAfter = currentDataAgg_dayNumSSLAfter.size; #the size
                                currentDataAgg_number = currentDataAgg_numberBefore + currentDataAgg_numberAfter; #sum em for ez
                                currentDataAgg_dayNumSSLExtended = np.hstack( (currentDataAgg_dayNumSSLBefore , currentDayNumF_singleSatLine , currentDataAgg_dayNumSSLAfter ) ); #days, stack em
                                                                
                            else:
                                currentDataAgg_numberBefore = 0; #set to 0
                                currentDataAgg_numberAfter = 0; #set to 0
                                currentDataAgg_number = 0; #set to 0
                            #END IF
                        
                            #make sure do work only on the day we want
                            #only get the data for times ranges thare are two times the filter period or longer
                            if( ( ((currentDayNumF[currentTimeSplits_loc[l+1]-1]-currentDayNumF[currentTimeSplits_loc[l]])*86400+currentDataAgg_number*TEC_dataRate) > filter_savGolPeriod*2) & #make sure the data time period is 2*filterPeriod
                               ( (currentTimeSplits_loc[l+1]-currentTimeSplits_loc[l]+currentDataAgg_number) >= windowLen_savGol) & #and that the data number is greater than the winodw length (in the event data is skipped and it's not perfectly 30 sec)
                               # ( (currentTimeSplits_loc[l+1]-currentTimeSplits_loc[l]+currentDataAgg_number) > (3*(filter_b.size-1)) ) & #and that the data number is greater than the filter padding length requirement (cannot be equal)
                               (np.any(np.int16(currentDayNumF[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                               
                                if( np.int16(currentDayNumF_singleSatLine[0]) != np.int16(currentDayNumF_singleSatLine[-1]) ): #if this is true, we are looking at different days' data being used
                                    #the problem with that is that it seems the site can have different TEC fit values for different days
                                    #so the lines are continuous - but have a big jump where the different days are fit to different means(?) I'm not sure
                                    #but there's def a big jump we're gonna fix here
                                    currentDayNum_singleSatLine_firstDayLocations = np.where( np.int16(currentDayNumF_singleSatLine[0]) == np.int16(currentDayNumF_singleSatLine) )[0]; #get index where 
                                    currentDayNum_singleSatLine_secondDayLocations = np.where( np.int16(currentDayNumF_singleSatLine[-1]) == np.int16(currentDayNumF_singleSatLine) )[0]; #get index where 
                                    currentTEC_singleSatLine_diff = np.diff(currentvTEC_singleSatLine); #get the deltas between the TEC values
                                    currentTEC_singleSatLine_firstDayLast = currentvTEC_singleSatLine[currentDayNum_singleSatLine_firstDayLocations[-1]]; #get the first day's last TEC value
                                    currentTEC_singleSatLine_secondDay = currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations]; #get the second day's TEC values
                                    if( currentDayNum_singleSatLine_firstDayLocations.size > 1 ): #make sure currentDayNum_singleSatLine_firstDayLocations[-2] exists
                                        currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_firstDayLocations[-2]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                        #sorry about these confusing explanations but it def makes sense to me rightnow
                                    else: #if it doesn't, time to use currentDayNum_singleSatLine_secondDayLocations[0] instead of currentDayNum_singleSatLine_firstDayLocations[-2] (so projecting from reverse instead of forward)
                                        currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_secondDayLocations[0]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                        #sorry about these confusing explanations but it def makes sense to me rightnow
                                    #END IF
                                    currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations] = currentTEC_singleSatLine_secondDay; #put em in
                                #END IF
                                
                                #now make sure the TEC data before and after connects to the current TEC data (cause of means and day differences and whatever)                                                        
                                if( currentDataAgg_numberBefore > 0 ): #make sure there's some before data
                                    currentDataAgg_TECBefore = currentDataAgg_TECBefore*( 2*currentvTEC_singleSatLine[0] - currentvTEC_singleSatLine[1] )/currentDataAgg_TECBefore[-1]; ##scale it so that the last Before data pt is the same as the first real data pt plus the delta between the first and second points
                                #END IF
                                if( currentDataAgg_numberAfter > 0 ): #make sure there's some after data
                                    currentDataAgg_TECAfter = currentDataAgg_TECAfter*( 2*currentvTEC_singleSatLine[-1] - currentvTEC_singleSatLine[-2] )/currentDataAgg_TECAfter[0]; ##scale it so that the first After data pt is the same as the last real data pt plus the delta between the last and second-to-last points
                                #END IF
                                currentDataAgg_TEC = np.hstack( ( currentDataAgg_TECBefore, currentvTEC_singleSatLine, currentDataAgg_TECAfter ) ); #stack em into one
                                
                                
                                #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                
                                #windowLen_savGol = np.int64(np.round(len(currentvTEC_singleSatLine)/4)) - (np.mod(np.int64(np.round(len(currentvTEC_singleSatLine)/4)),2) - 1); #unused
                                
                                currentPolyYvals = savgol_filter(currentDataAgg_TEC,windowLen_savGol,order_savGol ); #filter it up
                                current_deltavTEC = currentDataAgg_TEC - currentPolyYvals;
                                
                                #-----This code bit exists to catch where the fit greatly deviates from the start OR end, which seems to happen often-----
                                #use median/median distance instead of mean/stdev for stability apparently
                                #inspired by https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
                                current_deltaTEC_diff = np.diff(current_deltavTEC); #get the difference between each value and its following neighbor value
                                current_deltaTEC_dist = np.abs(current_deltaTEC_diff - np.median(current_deltaTEC_diff)); #get the absolute distance between the difference and the median of the difference
                                current_deltaTEC_comparitor = current_deltaTEC_dist/np.median(current_deltaTEC_dist); #scale the dist from the median by the median of the dist fromt he median
                                current_deltaTEC_diffWhere = np.where(current_deltaTEC_comparitor > deltaTEC_compareValue)[0]; #get where the scaled distance is greater than the hard-coded comparator value (set to 3.5 when I did this - is a variable)
                                
                                if( current_deltaTEC_diffWhere.size != 0 ): #makes sure the culling only occurs if it is needed (if current_deltaTEC_diffWhere is empty then everything is good)                   
                                    current_deltaTEC_diffWhere_diff = np.diff(current_deltaTEC_diffWhere); #the deltas of the indexes of the deltas that are too great am I right
                                    current_deltaTEC_diffWhere_diffWhere = np.where(current_deltaTEC_diffWhere_diff >= np.int64(deltaTEC_compareValue*2))[0]; #where the index difference is greater than deltaTEC_compareValue*2
                                    if( current_deltaTEC_diffWhere[0] <= np.int64(deltaTEC_compareValue*2) ): #if true, at the start, the actual values and the fit line greatly diverge
                                        if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the last current_deltaTEC_diffWhere since we're coming from the start
                                            current_deltaTEC_startCut = current_deltaTEC_diffWhere[ -1 ] + 1; #get the last difference where the delta was out of the delta-TEC-compareValue range and use it.
                                            #add the +1 cause of that slicing weirdness in python
                                        else:
                                            current_deltaTEC_startCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[0] ] + 1; #get the first big difference (where the start divergence ends - before the first big difference the index deltas are small till we hit the start!)
                                            #add the +1 cause of that slicing weirdness in python
                                        #END IF
                                    else:
                                        current_deltaTEC_startCut = 0; #set no values to be cut off
                                    #END IF
                                    if( current_deltaTEC_diffWhere[-1] >= (current_deltaTEC_diff.size - np.int64(deltaTEC_compareValue*2)) ): #if true, at the end, the actual values and the fit line greatly diverge
                                        if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the first current_deltaTEC_diffWhere since we're coming from the end
                                            current_deltaTEC_endCut = current_deltaTEC_diffWhere[0] + 1; #get the first diff where the delta was out of the deltaTEC_compareValue range and use it.
                                            #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
                                        else:
                                            current_deltaTEC_endCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[-1]+1 ] + 1; #get the last big difference (where the end divergence ends - after the last big difference the index deltas are small till we run hit the end!)
                                            #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
                                        #END IF
                                    else:
                                        current_deltaTEC_endCut = current_deltavTEC.size; #set no values to be cut off
                                    #END IF  
                                    
                                    current_deltavTEC[0:current_deltaTEC_startCut] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
                                    current_deltavTEC[current_deltaTEC_endCut:current_deltavTEC.size] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
                                #END IF
                                
                                #DEBUG PLOTS
                                #if( np.any( np.abs(current_deltavTEC) > 1.5) ):
#                                plt.figure();
#                                plt.scatter( currentDayNumF_singleSatLine , currentvTEC_singleSatLine , 20 , "r" );
#                                plt.scatter( currentDayNumF_singleSatLine , currentPolyYvals , 20 );
#                                plt.title("Fit Sav-Gol: Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                                plt.xlabel('Time [days]');
#                                plt.ylabel('los TEC [TECU]');
#                                plt.show();
    #                            
    #                            plt.figure();
    #                            plt.scatter( currentDayNumF_singleSatLine , currentElv_singleSatLine );
    #                            plt.title("Elevation Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Elevation [deg]');
    #                            plt.show();
    #                           
#                                if(  currentDataAgg_number > 0 ):
#                                    plt.figure();
#                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentDataAgg_TEC , 20 , "r" );
#                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentPolyYvals , 20 );
#                                    plt.title("Fit Sav-Gol: Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                                    plt.xlabel('Time [days]');
#                                    plt.ylabel('los TEC [TECU]');
#                                    plt.show();
#                                    
#                                    plt.figure();
#                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentDataAgg_TEC , 20 , "r" );
#                                    plt.scatter( currentDayNumF_singleSatLine , currentvTEC_singleSatLine , 20 );
#                                    plt.title("OG TEC B, Extra TEC R: Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                                    plt.xlabel('Time [days]');
#                                    plt.ylabel('los TEC [TECU]');
#                                    plt.show();
#                                    
#                                    plt.figure();
#                                    plt.scatter( currentDataAgg_dayNumSSLExtended , current_deltavTEC );
#                                    plt.title("Projected Delta TEC Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                                    plt.xlabel('Time [days]');
#                                    plt.ylabel('delta los TEC [TECU]');
#                                    plt.show();
#                                #END IF
                                
                                
    
                                
                                                            
                                # ===============Highpass filtering================
                                # current_deltavTEC = signal.filtfilt(filter_b,filter_a,current_deltavTEC,padtype='odd',padlen=3*(filter_b.size-1)); #Appies the filter
                                #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
                                
                                # ===============Cut to the correct size================
                                if( (currentDataAgg_numberBefore >= 0) & (currentDataAgg_numberAfter > 0) ):
                                    current_deltavTEC = current_deltavTEC[currentDataAgg_numberBefore:-currentDataAgg_numberAfter]; #remove any of the added TEC values - they were just used to increase the possible "visible" time of a satellite
                                elif( (currentDataAgg_numberBefore > 0) & (currentDataAgg_numberAfter == 0) ): #this needs some lil help
                                    current_deltavTEC = current_deltavTEC[currentDataAgg_numberBefore:currentDataAgg_dayNumSSLExtended.size+1]; #remove any of the added TEC values - they were just used to increase the possible "visible" time of a satellite
                                #END IF
                                
#                                if(  currentDataAgg_number > 0 ):
#                                    plt.figure();
#                                    plt.scatter( currentDayNumF_singleSatLine , current_deltavTEC );
#                                    plt.title("Trimmed Delta TEC Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
#                                    plt.xlabel('Time [days]');
#                                    plt.ylabel('delta los TEC [TECU]');
#                                    plt.show();
#                                    cntr = cntr + 1;
#                                    if( cntr == 11 ):
#                                        sys.crash();
#                                    #END IF
#                                #END IF
                                
                                # ===============Force mean to 0================
                                #I decided this is wrong
                                # current_deltavTEC = current_deltavTEC - np.nanmean(current_deltavTEC); #subtract the mean to make it a mean of 0
                                if( np.sum(np.isnan(current_deltavTEC)) == current_deltavTEC.size):
                                    sys.crash();
                                #END IF
                                                        
                            else:
    #                            currentvTEC_singleSatLine = currentvTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                current_deltavTEC = np.ones(currentTimeSplits_loc[l+1]-currentTimeSplits_loc[l],dtype=dataType_meth)*np.nan; #no filter, just filler (NaN to easy identify)
                            #END IF
                            
                            #no matter what, write it in!
    #                        TEC_fileData_float[currentSite_loc,locFloatUnfilt_dTEC][currentSat_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]]  = current_deltavTEC; #save that deltaTEC into wherever it should go - assigning with double indexing DOES NOT WORK
#                                TEC_fileData_float[currentSite_loc[currentSat_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]],locFloatUnfilt_dTEC] = current_deltavTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                            # TEC_fileData_dTEC[currentSite_loc[currentSat_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]]] = current_deltavTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                            #WRONG BELOW
    #                        TEC_fileData_float[currentSite_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]],locFloatUnfilt_dTEC] = current_deltavTEC; #save that deltaTEC into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                        #END FOR l
                    #END FOR k
                    if( (FLG_parallel == False) & ((np.mod(j+1,100) == 0) | (j == 0) | (j == unfilt_site_unique.size-1)) ): #report only every 100 sites, makes it go quicker but still often enough to keep you up to the minute updated
                        toc = time() - tic2; #for time testing
                        sys.stdout.write('\rSite {} of {} | {} min so far, {} sec per site, ETA {} min\t\t\t\t'.format(j+1,unfilt_site_unique.size,textNice(np.round(toc/60,2)),textNice(np.round(toc/(j+1),4)),textNice(np.round(toc/(j+1)*(unfilt_site_unique.size-j-1)/60,2)))); #report
                        sys.stdout.flush();   
                    #END IF
                    #END FOR j
            #END IF dataAggregation
            
            #-----COMBINE LISTS OF GOOD DATA INTO MECHA ARRAYS-----
            keyz = list(filt_dict.keys()); #get the current keys
            for j in range(0,len(keyz)):
                if( np.isscalar(filt_dict[keyz[j]]) == False ):
                    #if not a scalar, apply the logical mask
                    filt_dict[keyz[j]] = np.hstack(filt_dict[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
                #END IF
            #END FOR j
            
            if( FLG_parallel == False ):
                toc = time() - tic; #for time testing
                print("\nTime to filter all data: {} min".format(np.round(toc/60,2))); #extra space at end 
                tic = time(); #start a new timer
            #END IF
            
            #-----Get where the good data is-----
            #---NAN CONTROL---
            #Time to remove NaNs that were accumulated for TEC_fileData_float[:,0] (they mean data was bad)
#            TEC_fileData_logical_TECnotnans = np.where(np.logical_not(np.isnan(TEC_fileData_float[:,locFloatUnfilt_dTEC])) == True)[0]; #find the not NaNs (index seems to be a bit faster maybe)
            TEC_logical_TECnotnans = np.logical_not(np.isnan(filt_dict['dTEC'])); #find the not NaNs    
            #---NON-DAY CONTROL---
            TEC_logical_onDay = (filt_dict['dayNum'] == dateRange_dayNum_full[i,1]) & (filt_dict['year'] == dateRange_dayNum_full[i,0]); #find when the day reported is the day we want and the year
            #---OUTLIER CONTROL---
            #adapted from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list b/c median and median dist are more stable
            distFromMedian = np.abs(filt_dict['dTEC'] - np.nanmedian(filt_dict['dTEC'])); #get distance from median
            medianDistFromMedian = np.nanmedian(distFromMedian); #get median distance from median
            distFromMedian = distFromMedian/medianDistFromMedian if medianDistFromMedian else 0.0; #get normalized distance from median [reuse var to save memory]
            refDistFromMedian = np.abs(TEC_maxAmpAllowed - np.nanmedian(filt_dict['dTEC'])); #distance from median for 6 TECU [for dTEC that's a lot]
            refDistFromMedian = refDistFromMedian/medianDistFromMedian; #norm it so it can be used to compare
            TEC_logical_withinBounds = distFromMedian < refDistFromMedian; #get where the normalized distance from median is under the reference ceiling normalized value
            del distFromMedian
            #---MIN ELEVATION CONTROL---
            TEC_logical_withinElevation = filt_dict['elev'] >= minElevation; #get the locations where elevation is above the min elevation
            #---COMBINE ALL CONTROLS---
            TEC_logical_combined = np.where(TEC_logical_TECnotnans & TEC_logical_onDay & TEC_logical_withinBounds & TEC_logical_withinElevation)[0]; #combine them, get the index (seems slightly faster maybe)
            del TEC_logical_TECnotnans, TEC_logical_onDay, TEC_logical_withinBounds, TEC_logical_withinElevation; #save RAM
            #---CLEAR OUT STUFF---
            keyz = list(filt_dict.keys()); #get the current keys
            for j in range(0,len(keyz)):
                if( np.isscalar(filt_dict[keyz[j]]) == False ):
                    #if not a scalar, apply the logical mask
                    filt_dict[keyz[j]] = filt_dict[keyz[j]][TEC_logical_combined]; #keep only the good stuff
                #END IF
            #END FOR j
            del TEC_logical_combined; #save RAM
                            
            #-----TIME TO SAVE FILTERED DATA!-----
            #with h5py.File(TEC_dataFilePath[i], 'w') as TEC_file:
            h5pyChunkShape = filt_dict['dTEC'].shape; #get the shape of one of the vectors and use it as a chunk (read only whole chunks)
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'w', rdcc_nbytes =500*1024*1024) as TEC_file:
                for j in range(0,len(keyz)):
                    if( np.isscalar(filt_dict[keyz[j]]) == False ):
                        TEC_file.create_dataset(keyz[j], data=filt_dict[keyz[j]], chunks=h5pyChunkShape, compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #write that data
                    else:
                        #if size 1, add it as an attribute
                        TEC_file.attrs[keyz[j]] = filt_dict[keyz[j]]; #save the attribute
                    #END IF
                    del filt_dict[keyz[j]]; #clean up that memory
                    if( FLG_parallel == False ):
                        sys.stdout.write("\rWriting {} to file & {} min | {} out of {}\t\t\t\t\t".format(keyz[j],np.round((time()-tic)/60,2),j+1,len(keyz)));
                        sys.stdout.flush();
                    #END IF
                #END FOR j
                #add on non-data-related attributes
                if( len(TEC_fileData_paddingWarning) > 0 ):
                    TEC_file.attrs['paddedDayMissing'] = " and ".join(TEC_fileData_paddingWarning) + " missing"; #record the attribute
                else:
                    TEC_file.attrs['paddedDayMissing'] = "None, before and after padded"; #record the attribute
                #END IF
                TEC_file.attrs['TECmaxAmplitudeAllowed'] = TEC_maxAmpAllowed; #record the attribute
                TEC_file.attrs['pierceAlt'] = unfilt_dict['pierceAlt']; #record the attribute
                TEC_file.attrs['forcedTECdataRateSec'] = TEC_dataRate; #record the attribute
                # TEC_file.attrs['medianRejectRatio'] = deltaTEC_compareValue; #record the attribute
                TEC_file.attrs['savgolFiltPeriodSec'] = filter_savGolPeriod; #record the attribute
                # TEC_file.attrs['highpassFiltPeriodHr'] = filter_cutoffPeriod; #record the attribute
                # TEC_file.attrs['highpassFiltOrd'] = filter_n+1; #record the attribute
                # TEC_file.attrs['highpassFiltWindow'] = "FIR Hanning"; #record the attribute
                # TEC_file.attrs['highpassFiltType'] = "signal.filtfilt, padtype='odd',padlen=3*(b.size-1)"; #record the attribute
                TEC_file.attrs['version'] = version_filt; #record the filtered algorithm version
                if( FLG_parallel == False ):
                    print('\nDone writing filtered file!');
                #END IF
            #END WITH
            
                            
            #----DELETE UNNEEDED VARS----
            #Save memory, done with it
            del unfilt_dict;
            del filt_dict; #clean memory
            
            TEC_dataAvail[i] = 7; #data filtered and done! setting to 7 is good, setting to 1 is bad because 7 implies 4 while 1 doesn't guarantee that at alllllll
            #this is an artform of index juggling buddi
            
            if( FLG_parallel == False ):
                toc = time() - tic; #for time testing
                print('\nTime to save filtered data: {} min'.format(np.round(toc/60,2))); #extra space at end 
            #END IF
            
        else:
            if( TEC_dataAvail[i] == 4 ): #only do this if there is data there
                TEC_dataAvail[i] = -3; #data is a padded day and the filtering is done for non-padded days, so these are donezo
            #END IF
        #END IF
        
    #END IF  
    
    if( (FLG_deleteUnfilt == 1) & ((TEC_dataAvail[i] == 4) | (TEC_dataAvail[i] == 7) | (TEC_dataAvail[i] == 1) | (TEC_dataAvail[i] == -2)) ): #if delete unfiltered data is on - delete it
        os.remove(TEC_dataFilePathUnfilt[i]); #delete the unfiltered file instead of keeping it for more filtering stuff (saves hard drive space)
    #END IF
    
    return TEC_dataAvail[i]
#END DEF