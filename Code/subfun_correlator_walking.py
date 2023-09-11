#  data1 stays static and data2 slides if FLG_correlator == 1
#  data2 stays static and data1 slides if FLG_correlator == 2
# -NEG lags mean the sliding data was lagged (backward) in time, +POS lags mean the sliding data was pushed forward in time

import numpy as np
# import joblib
from copy import deepcopy
from Code.subfun_timeMatch import subfun_timeMatch
# from Code.subfun_correlator import subfun_correlator
from Code.subfun_strfind import strfind
# from Code.subfun_strstr import strstr
# from Code.subfun_textNice import textNice
from Code.subfun_correlator_corraler import subfun_correlator_corraler
import time

#time2bound should be an array vector/tuple/list of size 2 like [inclusive start, inclusive stop] that's in ABSOLUTE time for the year (not relative zero hour'd)
# None default works right at least

#FLG_clipData2 if True will allow for time ranges to be tested that cause data2 to run out of data for some time offsets, 
#   False prevents that and only allows for time ranges where data2 has full data availability (may clip beginning/end of time2bound or data1 if time2bound=None)

def subfun_correlator_walking(
        data1, settings1, name1, data1_setNames, \
        data2, settings2, name2, data2_setNames, \
        dates, settings_plot, settings_paths, settings_config, \
        FLG_correlator, FLG_correlator_options_filler, FLG_correlator_plot = False, 
        FLG_correlator_tabulator = False,  FLG_enableText = False, \
        FLG_correlator_shiftDir = 'both', FLG_correlator_timeLim = 14400, \
        filt1 = None,  filt2 = None, settings_spectra=None, reportDivisor=[60,'min'],
        time2span=6*3600, time2step=1*3600, time2bound=None, FLG_clipData2=False, FLG_nanLimitData2=False):    
        
        if( not isinstance(data1_setNames,list) and not isinstance(data1_setNames,tuple) ):
            data1_setNames = [data1_setNames]; #make it a list
        #END IF
        if( not isinstance(data2_setNames,list) and not isinstance(data2_setNames,tuple) ):
            data2_setNames = [data2_setNames]; #make it a list
        #END IF
        data1_setNames = data1_setNames.copy(); #protect it from possible edits
        data2_setNames = data2_setNames.copy(); #protect it from possible edits
        
        if( settings1 == '!DIRECT' ):
            #---FOR THIS CASE data2 MUST BE A LIST OF [data2, timeUnique2] so things still work without the whole dict construct---
            #---name1 and data1_setNames should represent the <name1=name of overarching data type> and <data1_setNames=sub-set of data>--- 
            #(for real data it would be name1='OMNI', data1_setNames='SYM/H') 
            #(like for double keo data, name1='Delta-vTEC E USA' and data1_setNames='Delta-vTEC E USA' too so the name and tabulation files are named rightish)
            #denotes direct data/label were sent in, not a "standardized" dict
            data1 = {data1_setNames[0]:data1[0].copy(), 'time unique':data1[1].copy(), 'data rate':np.median(np.diff(data1[1]))}; #wrap it so it works
            settings1 = {'labels':{data1_setNames[0]:data1_setNames[0]},'data type':data1_setNames[0]}; #wrap it so it works
        else:
            data1 = deepcopy(data1); #make sure outside is fine
        #END IF
        
        if( settings2 == '!DIRECT' ):
            #---FOR THIS CASE data2 MUST BE A LIST OF [data2, timeUnique2] so things still work without the whole dict construct---
            #---name1 and data1_setNames should represent the <name1=name of overarching data type> and <data1_setNames=sub-set of data>--- 
            #(for real data it would be name1='OMNI', data1_setNames='SYM/H') 
            #(like for double keo data, name1='Delta-vTEC E USA' and data1_setNames='Delta-vTEC E USA' too so the name and tabulation files are named rightish)
            #denotes direct data/label were sent in, not a "standardized" dict
            data2 = {data2_setNames[0]:data2[0].copy(), 'time unique':data2[1].copy(), 'data rate':np.median(np.diff(data2[1]))}; #wrap it so it works
            settings2 = {'labels':{data2_setNames[0]:data2_setNames[0]},'data type':data2_setNames[0]}; #wrap it so it works
        else:
            data2 = deepcopy(data2); #make sure outside is fine
        #END IF
        
        #--- create walk range ---
        if( np.any(time2bound == None) ):
            time2bound_range = np.array( (data1['time unique'][0], data1['time unique'][-1]+data1['data rate']) );
            time2bound = data1['time unique'][-1]+data1['data rate'] - data1['time unique'][0]; #make the var if not provided based on the data range
            if( np.isclose(time2bound, np.int64(time2bound)) ):
                time2bound = np.int64(time2bound); #ensure its an integer
            #END IF
            limitr = np.ceil((time2bound)/time2step).astype(np.int64)+1; #how many times to go
        else:
            if( isinstance(time2bound,list) | isinstance(time2bound,tuple) ):
                time2bound = np.asarray(time2bound); #convert to numpy
            #END IF
            time2bound_range = time2bound.copy(); #copy the range
            time2bound = time2bound_range[-1] - time2bound_range[0]; #get it
            if( np.isclose(time2bound, np.int64(time2bound)) ):
                time2bound = np.int64(time2bound); #ensure its an integer
            #END IF
            limitr = np.ceil((time2bound)/time2step).astype(np.int64)+1; #how many times to go
        #END IF
        time_cutout_range_walking = np.tile(np.array( (0,time2span) ), (limitr,1) )+np.tile(np.arange(0,limitr)*time2step, (2,1)).T; #yeet
        
        #--- restrict walk range based on if front/back of data1 input is NaNs (happens if time match to data2 with a wider time range occurs - or bad data I guess) ---
        for d1 in range(0,len(data1_setNames)):
            jk =  np.isnan(data1[data1_setNames[d1]]);
            jk_where = np.where(jk)[0];
            jk_where_diff_not1 = np.where(np.diff(jk_where) != 1 )[0];
            if( jk[0] == True ):
                if( jk_where_diff_not1.size > 0 ):
                    time2ditch_ditcher = jk_where[jk_where_diff_not1[0]]+1; #spot to ditch at
                else:
                    time2ditch_ditcher = jk_where[-1]; #spot to ditch at (entire front is NaNs apparently)
                #END IF
                if( time2ditch_ditcher > 180/data1['data rate'] ): #make sure just some minor bits at the front aren't missing
                    time2ditch = np.where(time_cutout_range_walking < data1['time unique'][time2ditch_ditcher] - time2bound_range[0])[0];
                    time_cutout_range_walking = np.delete(time_cutout_range_walking,time2ditch,axis=0); #ensure nothing goes over the data1 data range
                #END IF
            #END IF
            if( jk[-1] == True ):
                if( jk_where_diff_not1.size > 0 ):
                    time2ditch_ditcher = jk_where[jk_where_diff_not1[0]+1]-1; #spot to ditch at
                else:
                    time2ditch_ditcher = jk_where[0]; #spot to ditch at (entire front is NaNs apparently)
                #END IF
                if( (jk.size - time2ditch_ditcher) > 180/data1['data rate'] ): #make sure just some minor bits at the end aren't missing
                    time2ditch = np.where(time_cutout_range_walking > data1['time unique'][time2ditch_ditcher] - time2bound_range[0])[0];
                    time_cutout_range_walking = np.delete(time_cutout_range_walking,time2ditch,axis=0); #ensure nothing goes over the data1 data range
                #END IF
            #END IF
        #END FOR d1
        
        #--- restrict walk range based on if front/back of data2 input is NaNs (happens if time match to data2 with a wider time range occurs - or bad data I guess) ---
        if( FLG_nanLimitData2 == True ):
            #turned off normally b/c for large sets there can be a single station outage that limits the entire group
            for d2 in range(0,len(data2_setNames)):
                
                #&| means sub-dict, $| means sub-array
                if( (strfind(data2_setNames[d2],'&|',opt=1) > 0) | (strfind(data2_setNames[d2],'$|',opt=1) > 0) ):
                    # data2_adder = '$\mathregular{_F}$'; #default to F
                    data2_setNames_subSegs = data2_setNames[d2].split('&|'); #get the sub-dicts
                    # data2_setNames_base = data2_setNames_subSegs[0]; #set the 1st thing to be the thing for the label to work
                    data2_alias = data2; #get an alias going
                    for d2s in range(0,len(data2_setNames_subSegs)):
                        if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) == 0 ): #check for sub-array for the sub-dict
                            data2_alias = data2_alias[data2_setNames_subSegs[d2s]]; #move down in the data structure
                        else:
                            data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                            data2_alias = data2_alias[data2_setNames_subSegs_subVects[0]]; #move down in the data structure
                            if( data2_alias.shape[0] > data2_alias.shape[1] ):
                                data2_alias = data2_alias[:,int(data2_setNames_subSegs_subVects[1])]; #choose this array based on smaller array dim being the component
                            else:
                                data2_alias = data2_alias[int(data2_setNames_subSegs_subVects[1]),:]; #choose this array based on smaller array dim being the component
                            #END IF
                            # if( data2_setNames_subSegs_subVects[1] == '0' ):
                            #     data2_adder = '$\mathregular{_X}$';
                            # elif( data2_setNames_subSegs_subVects[1] == '1' ):
                            #     data2_adder = '$\mathregular{_Y}$';
                            # elif( data2_setNames_subSegs_subVects[1] == '2' ):
                            #     data2_adder = '$\mathregular{_Z}$';
                            # elif( data2_setNames_subSegs_subVects[1] == '3' ):
                            #     data2_adder = '$\mathregular{_T}$'; #(By^2 + Bz^2)^(1/2) - not implemented yet
                            # #END IF
                        #END IF
                    #END FOR d2s
                    jk = np.isnan(data2_alias); #get it from the aliased data
                else:
                    #it's straightforward - no shennanigans needed
                    jk =  np.isnan(data2[data2_setNames[d2]]);
                #END IF
                
                jk_where = np.where(jk)[0];
                jk_where_diff_not1 = np.where(np.diff(jk_where) != 1 )[0];
                if( jk[0] == True ):
                    if( jk_where_diff_not1.size > 0 ):
                        time2ditch_ditcher = jk_where[jk_where_diff_not1[0]]+1; #spot to ditch at
                    else:
                        time2ditch_ditcher = jk_where[-1]; #spot to ditch at (entire front is NaNs apparently)
                    #END IF
                    if( time2ditch_ditcher > 180/data2['data rate'] ): #make sure just some minor bits at the front aren't missing (basically allow for 180 sec gaps, which is 3 60s data pts)
                        time2ditch = np.where(time_cutout_range_walking < data2['time unique'][time2ditch_ditcher] - time2bound_range[0])[0];
                        time_cutout_range_walking = np.delete(time_cutout_range_walking,time2ditch,axis=0); #ensure nothing goes over the data2 data range
                    #END IF
                #END IF
                if( jk[-1] == True ):
                    if( jk_where_diff_not1.size > 0 ):
                        time2ditch_ditcher = jk_where[jk_where_diff_not1[0]+1]-1; #spot to ditch at
                    else:
                        time2ditch_ditcher = jk_where[0]; #spot to ditch at (entire back is NaNs apparently)
                    #END IF
                    if( (jk.size - time2ditch_ditcher) > 180/data2['data rate'] ): #make sure just some minor bits at the end aren't missing (basically allow for 180 sec gaps, which is 3 60s data pts)
                        time2ditch = np.where(time_cutout_range_walking > data2['time unique'][time2ditch_ditcher] - time2bound_range[0])[0];
                        time_cutout_range_walking = np.delete(time_cutout_range_walking,time2ditch,axis=0); #ensure nothing goes over the data2 data range
                    #END IF
                #END IF
            #END FOR d2
        #END IF
        
        #--- restrict walk range to data1 availability on data1 edges ---
        #now limit time_cutout_range_walking based on if there's data1 data to support that walk (subfun_correlator_corraler shouldn't error out, but it won't give a full picture - so just don't take the picture)
        #check neg range
        time2cut1 = np.where(time_cutout_range_walking < (data1['time unique'][0]-time2bound_range[0]))[0]; #get where to cut
        time_cutout_range_walking = np.delete(time_cutout_range_walking,time2cut1,axis=0); #ensure nothing goes over the data1 data range
        #check pos range
        time2cut1 = np.where(time_cutout_range_walking > (data1['time unique'][-1]+data1['data rate']-time2bound_range[0]))[0]; #get where to cut
        time_cutout_range_walking = np.delete(time_cutout_range_walking,time2cut1,axis=0); #ensure nothing goes over the data1 data range
        
        #--- restrict walk range to data2 availability on data2 edges (for the time lag check, otherwise the start(+lag)/end(-lag) 
        #   will be clipped and those respective lags won't be analyzed b/c no data 
        #   (but the opposite lags will be analyzed) if FLG_clipData2 is True) ---
        if( FLG_clipData2 == False ):
            #now limit time_cutout_range_walking based on if there's data2 data to support that walk (subfun_correlator_corraler shouldn't error out, but it won't give a full picture - so just don't take the picture)
            #check neg range
            time2cut2 = np.where(time_cutout_range_walking < (FLG_correlator_timeLim+data2['time unique'][0]-time2bound_range[0]))[0]; #get where to cut
            time_cutout_range_walking = np.delete(time_cutout_range_walking,time2cut2,axis=0); #ensure nothing goes over the data1 data range
            #check pos range
            time2cut2 = np.where(time_cutout_range_walking > (data2['time unique'][-1]+data2['data rate']-time2bound_range[0]-FLG_correlator_timeLim))[0]; #get where to cut
            time_cutout_range_walking = np.delete(time_cutout_range_walking,time2cut2,axis=0); #ensure nothing goes over the data1 data range
        #END IF
        
        time_cutout_range_walking = time_cutout_range_walking - (dates['date range zero hr dayNum'][1]*86400 - time2bound_range[0]); #remove zero hour time offset (absolute time reference has been time2bound_range[0] so align to that
        #--- restrict based on time2bound_range (added it later so just a bandaid, could be included above implicitly) ---
        time_cutout_range_walking = np.delete(time_cutout_range_walking,np.where((time2bound_range[0]-dates['date range zero hr dayNum'][1]*86400) > time_cutout_range_walking[:,0])[0],axis=0); #ensure nothing goes over the time2bound range
        time_cutout_range_walking = np.delete(time_cutout_range_walking,np.where((time2bound_range[-1]-dates['date range zero hr dayNum'][1]*86400) < time_cutout_range_walking[:,0])[0],axis=0); #ensure nothing goes over the time2bound range
        
        
        #note: this was quite a big hassle (specifically the MagCAN support), wooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooof
        #----- do time matching once so it doesn't have to be done repeatedly in the subfun_correlator_corraler (esp. brutal for MagCAN's 1 sec data rate) -----
        def dictDiver_insert(dictDiver, dictDiver_subdicts, dictDiver_insert, createMissing=False):
            #inspired by https://stackoverflow.com/a/49290758
            dictDiver_alias = dictDiver; #make an alias in here
            for dic in range(0,len(dictDiver_subdicts)-1): #dive through all but the last key (which is set)
                if( dictDiver_subdicts[dic] in dictDiver_alias ):
                    dictDiver_alias = dictDiver_alias[dictDiver_subdicts[dic]]
                elif( createMissing ):
                    dictDiver_alias = dictDiver_alias.setdefault(dictDiver_subdicts[dic], {}); #create the key as a dictionary
                else:
                    print('WARNING in dictDiver_insert: Dict Key '+str(dictDiver_subdicts[dic])+' not in the supplied dict and createMissing=False. Value was NOT inserted.');
                    return dictDiver
                #END IF
            #END FOR dic
            if( (dictDiver_subdicts[-1] in dictDiver_alias) or createMissing ):
                dictDiver_alias[dictDiver_subdicts[-1]] = dictDiver_insert;
            elif( (dictDiver_subdicts[-1] not in dictDiver_alias) and (not createMissing) ):
                print('WARNING in dictDiver_insert: Final Dict Key '+str(dictDiver_subdicts[-1])+' not in the supplied dict and createMissing=False. Value was NOT inserted.');
            #END IF
            
            return dictDiver #return the original that's been updated
        #END DEF
        if( data1['data rate'] > data2['data rate'] ): #to cover extended time snipping that affects whatever ends up as sig2 (which is why they're different)
            kr_dataB_extended = (data2['time unique'] >= (dates['date range full dayNum'][0,1]-np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400) & (data2['time unique'] < (dates['date range full dayNum'][-1,1]+1+np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400); #get stuff outside the date range
        elif( data1['data rate'] < data2['data rate'] ):
            kr_dataB_extended = (data1['time unique'] >= (dates['date range full dayNum'][0,1]-np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400) & (data1['time unique'] < (dates['date range full dayNum'][-1,1]+1+np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400); #get stuff outside the date range
        #END IF
        data1_dataRate_all = np.ones(len(data1_setNames),dtype=type(data1['data rate']))*data1['data rate']; #get an array of data rates to go through (to make sure all are changed)
        data2_dataRate_all = np.ones(len(data2_setNames),dtype=type(data2['data rate']))*data2['data rate']; #get an array of data rates to go through (to make sure all are changed)
        decorator_supporter = []; #prep
        for d1 in range(0,len(data1_setNames)):
            for d2 in range(0,len(data2_setNames)):
                # &| means sub-dict, $| means sub-array
                if( (strfind(data2_setNames[d2],'&|',opt=1) > 0) | (strfind(data2_setNames[d2],'$|',opt=1) > 0) ):
                    data2_setNames_subSegs = data2_setNames[d2].split('&|'); #get the sub-dicts
                    data2_alias = data2; #get an alias going
                    for d2s in range(0,len(data2_setNames_subSegs)):
                        if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) == 0 ): #check for sub-array for the sub-dict
                            data2_alias = data2_alias[data2_setNames_subSegs[d2s]]; #move down in the data structure
                        else:
                            data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                            data2_alias = data2_alias[data2_setNames_subSegs_subVects[0]]; #move down in the data structure
                            #here we are time matching the entire array, subfun_timeMatch is cool with it
                            # if( data2_alias.shape[0] > data2_alias.shape[1] ):
                            #     data2_alias = data2_alias[:,int(data2_setNames_subSegs_subVects[1])]; #choose this array based on smaller array dim being the component
                            # else:
                            #     data2_alias = data2_alias[int(data2_setNames_subSegs_subVects[1]),:]; #choose this array based on smaller array dim being the component
                            # #END IF
                        #END IF
                    #END FOR d2s
                    if( (np.isclose(data1_dataRate_all[d1],data2_dataRate_all[d2]) == False) & (data2_setNames_subSegs[0] not in decorator_supporter) ):
                        if( data1_dataRate_all[d1] > data2_dataRate_all[d2] ):
                            if( (data1['time unique'].size < data2['time unique'][kr_dataB_extended].size) & (FLG_correlator == 1) ):
                                #attempt to keep all data2_alias data available by extending data1['time unique'] so that data2_alias can be time matched properly but extend past the real data1[data2_setNames[d2]]/data1['time unique']
                                time1p = data1['time unique'].copy(); #copy
                                if( data1['time unique'][-1] < data2['time unique'][kr_dataB_extended][-1] ):
                                    time1p = np.concatenate( (time1p, np.arange(data1['time unique'][-1]+data1_dataRate_all[d1], data2['time unique'][kr_dataB_extended][-1]+data1_dataRate_all[d1]+FLG_correlator_timeLim, data1_dataRate_all[d1])) ); #tack it on
                                #END IF
                                if( data1['time unique'][0] > data2['time unique'][kr_dataB_extended][0] ):
                                    time1p = np.concatenate( (np.flip(np.arange(data1['time unique'][0]-data1_dataRate_all[d1], data2['time unique'][kr_dataB_extended][0]-data1_dataRate_all[d1]-FLG_correlator_timeLim, -data1_dataRate_all[d1])), time1p) ); #tack it on
                                #END IF
                                if( data2_alias.shape[0] == data2['time unique'].size ):
                                    data2_alias, data2_timeUnique = subfun_timeMatch(data2_alias[kr_dataB_extended, :], data2['time unique'][kr_dataB_extended], time1p, timeMatch_delta=data1_dataRate_all[d1], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                                else:
                                    data2_alias, data2_timeUnique = subfun_timeMatch(data2_alias[:, kr_dataB_extended], data2['time unique'][kr_dataB_extended], time1p, timeMatch_delta=data1_dataRate_all[d1], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                                #END IF
                            else:
                                if( data2_alias.shape[0] == data2['time unique'].size ):
                                    data2_alias, data2_timeUnique = subfun_timeMatch(data2_alias[kr_dataB_extended, :], data2['time unique'][kr_dataB_extended], data1['time unique'], timeMatch_delta=data1_dataRate_all[d1], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                                else:
                                    data2_alias, data2_timeUnique = subfun_timeMatch(data2_alias[:, kr_dataB_extended], data2['time unique'][kr_dataB_extended], data1['time unique'], timeMatch_delta=data1_dataRate_all[d1], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                                #END IF
                            #END IF
                            data2_dataRate_all[d2] = data1_dataRate_all[d1]; #set it
                            
                            #prep to set this value since it needs more work than others (data2 isn't linked in memory to data2_alias after the vector shennanigans I guess)
                            dictDiver_subdicts = []; #prep a list
                            for d2s in range(0,len(data2_setNames_subSegs)):
                                if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) == 0 ): #check for sub-array for the sub-dict
                                    dictDiver_subdicts.append(data2_setNames_subSegs[d2s]); #move down in the data structure
                                else:
                                    data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                                    dictDiver_subdicts.append(data2_setNames_subSegs_subVects[0]); #move down in the data structure
                                #END IF
                            #END FOR d2s
                                                        
                            data2 = dictDiver_insert(data2, dictDiver_subdicts, data2_alias, createMissing=False); #insert the data2_alias value
                        else:
                            if( (data2['time unique'].size < data1['time unique'][kr_dataB_extended].size) & (FLG_correlator == 2) ):
                                #attempt to keep all data2_alias data available by extending data1['time unique'] so that data2_alias can be time matched properly but extend past the real data1[data2_setNames[d2]]/data1['time unique']
                                time2p = data2['time unique'].copy(); #copy
                                if( data2['time unique'][-1] < data1['time unique'][kr_dataB_extended][-1] ):
                                    time2p = np.concatenate( (time2p, np.arange(data2['time unique'][-1]+data2_dataRate_all[d2], data1['time unique'][kr_dataB_extended][-1]+data2_dataRate_all[d2]+FLG_correlator_timeLim, data2_dataRate_all[d2])) ); #tack it on
                                #END IF
                                if( data2['time unique'][0] > data1['time unique'][kr_dataB_extended][0] ):
                                    time2p = np.concatenate( (np.flip(np.arange(data2['time unique'][0]-data2_dataRate_all[d2], data1['time unique'][kr_dataB_extended][0]-data2_dataRate_all[d2]-FLG_correlator_timeLim, -data2_dataRate_all[d2])), time2p) ); #tack it on
                                #END IF                            
                                data1[data1_setNames[d1]], data1_timeUnique = subfun_timeMatch(data1[data1_setNames[d1]][kr_dataB_extended], data1['time unique'][kr_dataB_extended], time2p, timeMatch_delta=data2_dataRate_all[d2], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            else:
                                data1[data1_setNames[d1]], data1_timeUnique = subfun_timeMatch(data1[data1_setNames[d1]][kr_dataB_extended], data1['time unique'][kr_dataB_extended], data2['time unique'], timeMatch_delta=data2_dataRate_all[d2], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            #END IF
                            data1_dataRate_all[d1] = data2_dataRate_all[d2]; #set it
                        #END IF
                        decorator_supporter.append(data2_setNames_subSegs[0]); #finished with that array (even if other indexes are called, the entire array was worked on
                    elif( (np.isclose(data1_dataRate_all[d1],data2_dataRate_all[d2]) == False) & (data2_setNames_subSegs[0] in decorator_supporter) ):
                        if( data1_dataRate_all[d1] > data2_dataRate_all[d2] ):
                            data2_dataRate_all[d2] = data1_dataRate_all[d1]; #set it
                        elif( data2_dataRate_all[d2] > data1_dataRate_all[d1] ):
                            data1_dataRate_all[d1] = data2_dataRate_all[d2]; #set it
                        #END IF
                    #END IF
                else:
                    #it's straightforward - no shennanigans needed
                    if( np.isclose(data1_dataRate_all[d1],data2_dataRate_all[d2]) == False ):
                        if( data1_dataRate_all[d1] > data2_dataRate_all[d2] ):
                            if( (data1['time unique'].size < data2['time unique'][kr_dataB_extended].size) & (FLG_correlator == 1) ):
                                #attempt to keep all data2[data2_setNames[d2]] data available by extending data1['time unique'] so that data2[data2_setNames[d2]] can be time matched properly but extend past the real data1[data2_setNames[d2]]/data1['time unique']
                                time1p = data1['time unique'].copy(); #copy
                                if( data1['time unique'][-1] < data2['time unique'][kr_dataB_extended][-1] ):
                                    time1p = np.concatenate( (time1p, np.arange(data1['time unique'][-1]+data1_dataRate_all[d1], data2['time unique'][kr_dataB_extended][-1]+data1_dataRate_all[d1]+FLG_correlator_timeLim, data1_dataRate_all[d1])) ); #tack it on
                                #END IF
                                if( data1['time unique'][0] > data2['time unique'][kr_dataB_extended][0] ):
                                    time1p = np.concatenate( (np.flip(np.arange(data1['time unique'][0]-data1_dataRate_all[d1], data2['time unique'][kr_dataB_extended][0]-data1_dataRate_all[d1]-FLG_correlator_timeLim, -data1_dataRate_all[d1])), time1p) ); #tack it on
                                #END IF
                                data2[data2_setNames[d2]], data2_timeUnique = subfun_timeMatch(data2[data2_setNames[d2]][kr_dataB_extended], data2['time unique'][kr_dataB_extended], time1p, timeMatch_delta=data1_dataRate_all[d1], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            else:
                                data2[data2_setNames[d2]], data2_timeUnique = subfun_timeMatch(data2[data2_setNames[d2]][kr_dataB_extended], data2['time unique'][kr_dataB_extended], data1['time unique'], timeMatch_delta=data1_dataRate_all[d1], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            #END IF
                            data2_dataRate_all[d2] = data1_dataRate_all[d1]; #set it
                        else:
                            if( (data2['time unique'].size < data1['time unique'][kr_dataB_extended].size) & (FLG_correlator == 2) ):
                                #attempt to keep all data2[data2_setNames[d2]] data available by extending data1['time unique'] so that data2[data2_setNames[d2]] can be time matched properly but extend past the real data1[data2_setNames[d2]]/data1['time unique']
                                time2p = data2['time unique'].copy(); #copy
                                if( data2['time unique'][-1] < data1['time unique'][kr_dataB_extended][-1] ):
                                    time2p = np.concatenate( (time2p, np.arange(data2['time unique'][-1]+data2_dataRate_all[d2], data1['time unique'][kr_dataB_extended][-1]+data2_dataRate_all[d2]+FLG_correlator_timeLim, data2_dataRate_all[d2])) ); #tack it on
                                #END IF
                                if( data2['time unique'][0] > data1['time unique'][kr_dataB_extended][0] ):
                                    time2p = np.concatenate( (np.flip(np.arange(data2['time unique'][0]-data2_dataRate_all[d2], data1['time unique'][kr_dataB_extended][0]-data2_dataRate_all[d2]-FLG_correlator_timeLim, -data2_dataRate_all[d2])), time2p) ); #tack it on
                                #END IF                            
                                data1[data1_setNames[d1]], data1_timeUnique = subfun_timeMatch(data1[data1_setNames[d1]][kr_dataB_extended], data1['time unique'][kr_dataB_extended], time2p, timeMatch_delta=data2_dataRate_all[d2], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            else:
                                data1[data1_setNames[d1]], data1_timeUnique = subfun_timeMatch(data1[data1_setNames[d1]][kr_dataB_extended], data1['time unique'][kr_dataB_extended], data2['time unique'], timeMatch_delta=data2_dataRate_all[d2], FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            #END IF
                            data1_dataRate_all[d1] = data2_dataRate_all[d2]+data2_setNames[d2]; #set it
                        #END IF
                    #END IF
                #END IF
            #END FOR d2
        #END FOR d1
        #finalize it b/c just did it for each data type above but couldn't edit the global data rates yet
        if( np.isclose(data1['data rate'],data2['data rate']) == False ):
            if( data1['data rate'] > data2['data rate'] ):
                data2['time unique'] = data2_timeUnique; #set this too
                data2['data rate'] = data1['data rate']; #set it
            else:
                data1['time unique'] = data1_timeUnique; #set this too
                data1['data rate'] = data2['data rate']; #set it
            #END IF
        #END IF
            
        time2bound_range_boost = time2bound_range.copy();
        time2bound_range_boost[-1] = time2bound_range_boost[-1] + time2span + FLG_correlator_timeLim; #this allows for the final correlation coeff to be calc'd (extend to negative if I ever make it so the walk represents middle of range not start)
        corrRet = [None for i in range(0,time_cutout_range_walking.shape[0])]; #hold it
        for i in range(0,time_cutout_range_walking.shape[0]):
            FLG_correlator_options = [{'mode':'range','time range':time_cutout_range_walking[i,:]} for _ in range(len(data2_setNames))]; #make the options
            corrRet[i] = subfun_correlator_corraler(
                    data1, settings1, name1, data1_setNames, \
                    data2, settings2, name2, data2_setNames, \
                    dates, settings_plot, settings_paths, settings_config, \
                    FLG_correlator, FLG_correlator_options, FLG_correlator_plot = FLG_correlator_plot, 
                    FLG_correlator_tabulator = FLG_correlator_tabulator, FLG_tabulator_extraBit = 'walking', FLG_enableText = FLG_enableText, \
                    FLG_correlator_shiftDir = FLG_correlator_shiftDir, FLG_correlator_timeLim = FLG_correlator_timeLim, FLG_correlator_timeBound = time2bound_range_boost, \
                    filt1 = filt1,  filt2 = filt2, settings_spectra=settings_spectra, reportDivisor=reportDivisor);
            # corrRet[i][0]['time range walking'] = time_cutout_range_walking; #also tack the walking on
        #END FOR i
        
        return corrRet, time_cutout_range_walking
#END DEF