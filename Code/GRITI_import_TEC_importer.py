import numpy as np
from time import time
import gc, sys
from Code.GRITI_import_TEC_Madrigal import GRITI_import_TEC_Madrigal
from Code.GRITI_import_TEC_otherSources import GRITI_import_TEC_otherSources
from Code.GRITI_import_TEC_LISN import GRITI_import_TEC_LISN
from Code.GRITI_TEC_randomSynth import GRITI_TEC_randomSynth
from Code.subfun_isin_row import isin_row
from Code.subfun_satConv import satConv_toInt_fromStr
from Code.subfun_textNice import textNice

def GRITI_import_TEC_importer(data, dates, settings, FLG_timeToCheck=False):
    FLG_dataTypes = settings['data']['data types req']; #unpack
    FLG_TECloc = settings['data']['data types'].index('TEC');
    
    FLG_justChecking = settings['TEC']['justChecking'];
    if( FLG_timeToCheck == True ):
        FLG_justChecking = False; #override, it's time to check for realsies
    #END IF
    
    if( (not 'TEC' in list(data.keys())) | (FLG_timeToCheck == True) ):
        #~~~IMPORT PRE-PROCESSED DATA~~~
        if( np.any(np.asarray(settings['TEC']['source to use']) == 0) | np.any(np.asarray(settings['TEC']['source to use']) == 1) ):
            if( FLG_dataTypes[FLG_TECloc] == 1):
                print("Importing TEC");
                # (TEC_int, TEC_float, TEC_str, FLG_dataTypes[FLG_TECloc]) = \
                (data['TEC'], FLG_dataTypes[FLG_TECloc]) = \
                    GRITI_import_TEC_otherSources(dates['date range full dayNum'],dates['date range zero hr dayNum'],settings['paths'],settings['TEC import']['lat range'],settings['TEC import']['long range'], \
                        TEC_dataRate=settings['TEC import']['TEC_dataRate'],TEC_minimumTimeGap=settings['TEC import']['TEC_minimumTimeGap'],
                        TEC_timeTolerance=settings['TEC import']['TEC_timeTolerance']); #import pre-processed TEC data from other sources
            #END IF
        #END IF
            #~~~IMPORT MADRIGAL DATA~~~
        if( np.any(np.asarray(settings['TEC']['source to use']) == 0) | np.any(np.asarray(settings['TEC']['source to use']) == 2) ):    
            if( ((FLG_dataTypes[FLG_TECloc] == 2) | np.any(np.asarray(settings['TEC']['source to use']) == 2)) & (FLG_dataTypes[FLG_TECloc] != 0) ):
                
                data['TEC'] = GRITI_import_TEC_Madrigal( dates, settings, FLG_justChecking=FLG_justChecking ); #import TEC data for the date range
                # if( 'justChecking' not in data['TEC'] ):
                #     TEC_dataRate = data['TEC']['data rate']; #record the data rate for functions and stuff
                # #END IF
                
            #END IF    
        #END IF
            #~~~IMPORT LISN DATA~~~
        if( np.any(np.asarray(settings['TEC']['source to use']) == 0) | np.any(np.asarray(settings['TEC']['source to use']) == 3) ):    
            if( FLG_dataTypes[FLG_TECloc] >= 1 ):
                TEC_dataTemp = GRITI_import_TEC_LISN( dates, settings, FLG_justChecking=FLG_justChecking ); #import TEC data for the date range
                                    
                if( 'justChecking' not in data['TEC'] ):
                    keyz_new = list(TEC_dataTemp.keys()); #get the keys
                    if( 'TEC' in list(data.keys()) ):
                        tic_LISN = time(); #time it
                        print('\nCombining LISN and current data together.');
                        keyz_curr = list(data['TEC'].keys()); #get the keys
                        TECcurr_site_unique, TECcurr_site_unique_indexes = np.unique(data['TEC']['site'] , return_inverse=True); #get unique site names and the indexes to speed up searching for where to look
                        #inspired by https://stackoverflow.com/questions/33281957/faster-alternative-to-numpy-where
                        TECcurr_site_unique_currentSiteArray = np.split(np.argsort(data['TEC']['site'], kind='mergesort'), np.cumsum(np.bincount(TECcurr_site_unique_indexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
                        TECtemp_site_unique, TECtemp_site_unique_indexes = np.unique(TEC_dataTemp['site'] , return_inverse=True); #get unique site names and the indexes to speed up searching for where to look
                        #inspired by https://stackoverflow.com/questions/33281957/faster-alternative-to-numpy-where
                        TECtemp_site_unique_currentSiteArray = np.split(np.argsort(TEC_dataTemp['site'], kind='mergesort'), np.cumsum(np.bincount(TECtemp_site_unique_indexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
                        TECtemp_site_isIn = np.isin(TECtemp_site_unique,TECcurr_site_unique); #precalc if it is in b/c ez
                        TECtemp_site_maskr = np.ones(TEC_dataTemp['site'].size,dtype=np.bool_); #preallocate
                        for j in range(0,TECtemp_site_unique.size):
                            if( TECtemp_site_isIn[j] ): #use the isin mask (default mask value is true, so just need to do work if site is in)
                                TECcurr_site_loc = TECcurr_site_unique_currentSiteArray[np.where(TECtemp_site_unique[j] == TECcurr_site_unique)[0].item()]; #pull it out of the pre-calc'd list of data locations for 1 site
                                TECtemp_site_loc = TECtemp_site_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                                if( TECcurr_site_loc.size == TECtemp_site_loc.size ):
                                    TECtemp_site_maskr[TECtemp_site_loc] = False; #don't keep it b/c data is already all there
                                else:
                                    #otherwise check each row in the temp data to see if it is duplicated in the current data - if duplicated do not keep
                                    #deep Q is of course how can 1 site have 2 completely different measurement lat/long at the same time in 2 diff repos - but for another day (maybe site names reused but diff actual sites??)
                                    # TECtemp_isin = isin_row(np.vstack((TEC_dataTemp['lat'][TECtemp_site_loc],TEC_dataTemp['long'][TECtemp_site_loc],TEC_dataTemp['time'][TECtemp_site_loc])).T, np.vstack((data['TEC']['lat'][TECcurr_site_loc],data['TEC']['long'][TECcurr_site_loc],data['TEC']['time'][TECcurr_site_loc])).T); #check row-wise if temp data is in current data - keep anything that's not duplicated (the most permissive mode)
                                    #moved to sat/satType/time instead of lat/long/time b/c assumed ionosphere altitude can differ lat/long but theoretically only 1 measurement for sat/satType/time/site should occur
                                    TECtemp_isin = isin_row(np.vstack((TEC_dataTemp['sat'][TECtemp_site_loc],satConv_toInt_fromStr(TEC_dataTemp['satType'][TECtemp_site_loc]),TEC_dataTemp['time'][TECtemp_site_loc])).T, np.vstack((data['TEC']['sat'][TECcurr_site_loc],satConv_toInt_fromStr(data['TEC']['satType'][TECcurr_site_loc]),data['TEC']['time'][TECcurr_site_loc])).T); #check row-wise if temp data is in current data - keep anything that's not duplicated (the most permissive mode)
                                    TECtemp_site_maskr[TECtemp_site_loc[TECtemp_isin]] = False; #don't keep it b/c data is already all there
                                #END IF
                            #END IF
                        #END FOR j
                        #now load in the temp data into the main array w/ the mask array
                        for j in range(0,len(keyz_new)):
                            if( keyz_new[j] in keyz_curr ):
                                if( isinstance(data['TEC'][keyz_new[j]],(np.ndarray)) ):
                                    data['TEC'][keyz_new[j]] = np.hstack( (data['TEC'][keyz_new[j]], TEC_dataTemp[keyz_new[j]][TECtemp_site_maskr]) ); #smack them together
                                elif( isinstance(data['TEC'][keyz_new[j]],(list)) ):
                                    if( TECtemp_site_maskr.size == TECtemp_site_maskr.sum() ):
                                        data['TEC'][keyz_new[j]] += TEC_dataTemp[keyz_new[j]]; #append if list
                                    else:
                                        data['TEC'][keyz_new[j]] += [b for a, b in zip(TECtemp_site_maskr.tolist(), TEC_dataTemp[keyz_new[j]]) if a]; #append the data to keep w/ some list shennanigans to use the bool array
                                    #END IF
                                elif( np.isscalar(data['TEC'][keyz_new[j]]) ):
                                    if( data['TEC'][keyz_new[j]] != TEC_dataTemp[keyz_new[j]]):
                                        print('WARNING: TEC scalar numbers don\'t match but not a big deal, deal with this later\n'+\
                                              'Orig: '+str(data['TEC'][keyz_new[j]])+'\t New: '+str(TEC_dataTemp[keyz_new[j]])); #tell it's not quite right
                                    #END IF
                                #END IF
                            else:
                                data['TEC'][keyz_new[j]] = TEC_dataTemp[keyz_new[j]]; #new addition
                                keyz_curr = list(data['TEC'].keys()); #get the keys again
                            #END IF
                        #END FOR j
                        print('Combining LISN and current data together took '+textNice(np.round(time()-tic_LISN,2))+' sec.');
                    else:
                        data['TEC'] = {}; #prep it as a dict sicne it doesn't exist
                        for j in range(0,len(keyz_new)):
                            data['TEC'][keyz_new[j]] = TEC_dataTemp[keyz_new[j]]; #create the dict
                        #END FOR j
                    #END IF
                else:
                    keyz_new = list(TEC_dataTemp.keys()); #get the keys
                    if( 'TEC' not in list(data.keys()) ):
                        data['TEC'] = {}; #prep it as a dict sicne it doesn't exist
                        for j in range(0,len(keyz_new)):
                            data['TEC'][keyz_new[j]] = TEC_dataTemp[keyz_new[j]]; #create the dict
                        #END FOR j
                    #END IF
                #END IF
                del TEC_dataTemp
                gc.collect(); #clean the garbage
            #END IF
        #END IF
        
        if( 'justChecking' not in data['TEC'] ):
            if( FLG_dataTypes[FLG_TECloc] >= 1 ):
                # TEC_timeUnique = np.unique(data['TEC']['time']); #days, get unique times (v useful)
                tic = time()
                data['TEC']['time unique'] = np.unique(data['TEC']['time']); #sec, get unique times (v useful)
                # TEC_timeUnique = data['TEC']['time unique']; #sec set it for reuse in main
                data['TEC']['time unique aligned'] = np.unique(data['TEC']['time aligned']); #sec, get unique times (v useful) [aligned has year support and is set around the defined 0 hour]
                # TEC_timeUniqueAligned = data['TEC']['time unique aligned']; #sec set it for reuse in main
                print('Time to unique: '+str(time()-tic))
            #END IF
                
            if( (settings['TEC']['noise settings']['noise mode'] >= 1) & (FLG_dataTypes[FLG_TECloc] >= 1) ): #replace the delta-vTEC data with random data OR random data with synth waves embedded
                data['TEC']['dTEC'] = GRITI_TEC_randomSynth(data['TEC']['dTEC'].size,data['TEC']['lat'],data['TEC']['long'],data['TEC']['time'], \
                    settings['TEC']['noise settings']['noise mean'],settings['TEC']['noise settings']['noise stdev'],settings['map']['Re'],dates['date range zero hr'], \
                    settings['map']['lat range'],settings['map']['long range'],settings['map']['lat autotick'],settings['map']['long autotick'], \
                    settings['TEC']['snyth wave settings']['lat range'],settings['TEC']['snyth wave settings']['long range'], \
                    settings['TEC']['snyth wave settings']['N'],settings['TEC']['snyth wave settings']['angle'],settings['TEC']['snyth wave settings']['phase'], \
                    settings['TEC']['snyth wave settings']['wave length'],settings['TEC']['snyth wave settings']['period'],settings['TEC']['snyth wave settings']['amplitude'], \
                    settings, np.max(np.abs(settings['TEC']['plot lim'])),settings['TEC']['noise settings']['noise mode'],FLG_plotStuff=1); #replace the delta-vTEC data with random data OR random data with synth waves embedded
            #END IF
        #END IF
    else:
        print('Importing TEC\nTEC data still in memory from previous run. No import needed.\n');
    #END IF
    if( FLG_dataTypes[FLG_TECloc] >= 1 ):
        if( 'justChecking' not in data['TEC'] ):
            if( np.all(np.asarray(data['TEC']['version']) == data['TEC']['version'][0]) ):
                settings['TEC']['version'] = data['TEC']['version'][0]; #settings['TEC'] will update with these automagically since mem is shared
                settings['TEC']['keo']['version'] = data['TEC']['version'][0];
            else:
                print('ERROR in MAIN: TEC Data Versions ('+str(data['TEC']['version'])+') are not all the same, write some code to deal with this');
                sys.crash();
            #END IF
            if( np.all(np.asarray(data['TEC']['pierceAlt']) == data['TEC']['pierceAlt'][0]) ):
                settings['TEC']['pierceAlt'] = data['TEC']['pierceAlt'][0]; #settings['TEC'] will update with these automagically since mem is shared
                settings['TEC']['keo']['pierceAlt'] = data['TEC']['pierceAlt'][0];
            else:
                print('ERROR in MAIN: TEC Data Versions ('+str(data['TEC']['pierceAlt'])+') are not all the same, write some code to deal with this');
                sys.crash();
            #END IF
        #END IF
    #END IF
    
    return data, settings
#END DEF