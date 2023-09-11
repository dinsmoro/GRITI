#Function to import SuperMAG data from the pre-downloaded data
#RD on 12/10/2018

'''
GOAL: Get SuperMAG data
expecting: 
 dateRange_dayNum_full - [year,day#min;year,day#inbetween;year,day#max] numerical format. So, say [2013,126;2013,127;2013,128] for 2013 Day 126 to Day 128. No chars pls
 folder - list with folder[0] being the running folder (where the code is) and folder[1] being where the data is kept (can be the same, or not)
NOTE: SuperMAG is every 1 minute.
NOTE: Download year-long all stations from https://supermag.jhuapl.edu/mag/ and weep 10G files

'''

import numpy as np
import h5py
import os
import copy
from urllib.request import urlopen
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date

def GRITI_import_SuperMAG_stations(dateRange_dayNum_full, SuperMAGstations_names, settings_paths, settings_config, coordType='geo', FLG_overwrite=False, FLG_detrended=False):
    #--- Prep housekeeping ---
    version_alg = 1.0; #algorithm version
    #1.0 5/26/2022 - initial algorithm
    
    #--- Prepare web interface if needed ---
    site_inventory_base = 'https://supermag.jhuapl.edu/services/inventory.php?fmt=json&nohead&logon='+settings_config['login SuperMAG']['user']; #+'&start=2019-10-15T10:40&extent=3600'
    site_station_base = 'https://supermag.jhuapl.edu/services/data-api.php?fmt=json&nohead&logon='+settings_config['login SuperMAG']['user']; #+'&start=2019-10-15T10:40&extent=3600&all&station=NCK'
    #site_indices_base =  https://supermag.jhuapl.edu/services/indices.php?fmt=json&nohead&logon=YOURNAME&start=2019-10-15T10:40&extent=3600&all

    #--- Prep vars needed ---
    if( coordType == 'mag' ):
        coordType = 'magC'; #special for here b/c I name stuff dumdum
    #END IF
    if( FLG_detrended == True ):
        detrendedString = ' delta'; #to access SuperMAG detrended data (good for comparisons with SMR, SMU, etc. b/c they're based on detrended data
    else:
        detrendedString = ''; #to access SuperMAG un-detrended data (default, live your own life)
    #END IF
    # dateRange_yearRange = np.arange(dateRange_dayNum_full[0,0],dateRange_dayNum_full[-1,0]+1,1,dtype=np.int16); #get the full year range from min to max
    path_base = os.path.join(settings_paths['data'],'SuperMAG','Stations'); #prep path base
    SuperMAGstations_data = {}; #prep dict
    if( (not isinstance(SuperMAGstations_names,list)) & (not isinstance(SuperMAGstations_names,tuple)) ):
        SuperMAGstations_names = [SuperMAGstations_names]; #make sure list so can iterate over it
    #END IF
    #Read in the data
    for i in range(0,dateRange_dayNum_full.shape[0]): #steps through each day
        FLG_getStation = False; #prime every time
        path_yeared = os.path.join(path_base,str(dateRange_dayNum_full[i,0])); #add in day
        if( os.path.isdir(path_yeared) == False ):
            os.makedirs(path_yeared); #make the path
        #END IF
        for stat in range(0,len(SuperMAGstations_names)):
            SuperMAGstations_names[stat] = SuperMAGstations_names[stat].upper(); #make sure all caps (safety 1st)
            # jk = np.where(dateRange_yearRange[i] == dateRange_dayNum_full[:,0])[0]; #get where dateRange is the current year
            #check to see if data is already processed
            SuperMAG_dataFilePath =  os.path.join(path_yeared,'SuperMAGstations_' + str(dateRange_dayNum_full[i,0]) + '_' + str(dateRange_dayNum_full[i,1]) + '.h5'); #create the expected data path
            if( os.path.isfile(SuperMAG_dataFilePath) & (FLG_overwrite == False) ):
                SuperMAGstations_TEMP = {}; #read it in
                # keyz = list(SuperMAGstations_data.keys()); #get the current keys in AMPERE_temp     
                try:
                    with h5py.File(SuperMAG_dataFilePath, 'r') as h5file:
                        #version check first!
                        h5file_version = h5file.attrs['version']; #get that version
                        if( h5file_version != version_alg ):
                            print('WARNING in GRITI_import_SuperMAG_stations: File '+SuperMAG_dataFilePath+
                                  'for '+str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+
                                  ' has OLD version '+str(h5file_version)+', current version is '+str(version_alg)+
                                  '. Renaming and rebuilding!');
                            FLG_getStation = True; #needa get that station data
                            os.rename(SuperMAG_dataFilePath, \
                                      SuperMAG_dataFilePath+'_oldV'+str(h5file_version).replace('.','p')); #rename
                        else:
                            keyzNew = list(h5file.keys()); #get the saved keys
                            if( SuperMAGstations_names[stat] not in keyzNew ):
                                #make sure station is in the saved file
                                FLG_getStation = True; #needa get that station data
                            #END IF
                        #END IF
                    #END WITH
                except:
                    #looks like file is bjorked
                    print('WARNING in GRITI_import_SuperMAG_stations: File '+SuperMAG_dataFilePath+
                          'for '+str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+
                          ' is unreadable. Renaming and rebuilding!');
                    FLG_getStation = True; #needa get that station data
                    os.rename(SuperMAG_dataFilePath, \
                              SuperMAG_dataFilePath+'_BAD'); #rename
                #END TRY
                
                if( FLG_getStation == False ):
                    with h5py.File(SuperMAG_dataFilePath, 'r') as h5file:
                        #--- Read the keyed data in ---
                        # keyzStation = list(h5file.keys()); #get the saved keys
                        # jk = keyzStation.index(SuperMAGstations_names[stat]); #get index where station name is
                        SuperMAGstations_stationDict = {}; #prep a dict
                        keyzNew = list(h5file.get(SuperMAGstations_names[stat]).keys()); #get the group keys
                        for j in range(0,len(keyzNew)):
                            if( type(h5file.get(SuperMAGstations_names[stat]).get(keyzNew[j])) is h5py.Group ):
                                sub_keyz = list(h5file.get(SuperMAGstations_names[stat]).get(keyzNew[j]).keys()); #get the group keys
                                SuperMAGstations_stationDict[keyzNew[j]] = {}; #declare it a dict
                                for k in range(0,len(sub_keyz)):
                                    SuperMAGstations_stationDict[keyzNew[j]][sub_keyz[k]] = h5file[SuperMAGstations_names[stat]][keyzNew[j]].get(sub_keyz[k])[()]; #get that dataset out
                                #END FOR k
                                #--- Read the attributes in ---
                                sub_keyz = list(h5file.get(SuperMAGstations_names[stat]).get(keyzNew[j]).attrs.keys()); #get the attribute keys
                                # keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                                for k in range(0,len(sub_keyz)):
                                    SuperMAGstations_stationDict[keyzNew[j]][sub_keyz[k]] = h5file[SuperMAGstations_names[stat]][keyzNew[j]].attrs[sub_keyz[k]]; #get that attribute out
                                #END FOR k
                            else:
                                #simplified b/c attributes are read differently and not in the keyz
                                SuperMAGstations_stationDict[keyzNew[j]] = h5file.get(SuperMAGstations_names[stat]).get(keyzNew[j])[()]; #get that dataset out
                            #END IF
                        #END FOR j
                        #--- Read the attributes in ---
                        keyzNew = list(h5file.get(SuperMAGstations_names[stat]).attrs.keys()); #get the attribute keys
                        # keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                        for j in range(0,len(keyzNew)):
                            SuperMAGstations_stationDict[keyzNew[j]] = h5file.get(SuperMAGstations_names[stat]).attrs[keyzNew[j]]; #get that attribute out
                        #END FOR j
                    #END WITH
                    
                    if( ('mag'+detrendedString in SuperMAGstations_stationDict['geo'].keys()) & ('mag'+detrendedString in SuperMAGstations_stationDict['magC'].keys()) ):
                        FLG_getStation = False; #no need to get station data, it had it!
                    else:
                        FLG_getStation = True; #needs whatever sub-mag type was requested
                    #END IF
                #END IF
            else:
                FLG_getStation = True;
            #END IF
            
            if( FLG_getStation == True ):
                dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #can't use date dict here because date range may be artificially extended
                #--- See if station needed is there ---
                site_inventory = site_inventory_base + '&start='+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1])+'-'+str(dateRange_full[i,2])+'T00:00&extent='+str(86400); #create full URL
                try:
                    SuperMAGstations_raw = urlopen(site_inventory); #download the data needed
                    SuperMAGstations_rawRead = SuperMAGstations_raw.read(); #read it into a much more useful format
                    charset = SuperMAGstations_raw.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                    if( charset == None ):
                        charset = 'utf-8'; #set it to this
                    #END IF
                    SuperMAGstations_stations = SuperMAGstations_rawRead.decode(charset).split('\n')[:-1]; #"decode" the HTML content so it's more legible
                except:
                    import requests #requests doesn't seem to have cert errors, use it if urlopen is complainy
                    SuperMAGstations_raw = requests.get(site_inventory, stream=True); #get the data, stream is key I think
                    # SuperMAGstations_raw.encoding = SuperMAGstations_raw.apparent_encoding; #not needed, UTF-8 detected correctly
                    SuperMAGstations_stations = SuperMAGstations_raw.text.split('\n')[:-1]; #decode to text
                #END TRY
                
                if( int(SuperMAGstations_stations[0]) == len(SuperMAGstations_stations)-1 ):
                    SuperMAGstations_stations.pop(0); #get rid of the #
                else:
                    print('ERROR in GRITI_import_SuperMAG_stations: '+site_inventory+' returned incorrect # of sites actually listed! Check this out and fix it in the code. Crashing!');
                    import sys
                    sys.crash(); #yeet
                #END
                                
                if( SuperMAGstations_names[stat] in SuperMAGstations_stations ):
                    import json #for reading json dicts
                    import datetime #for reading time since 1970
                    if( FLG_detrended == False ):
                        site_station = site_station_base + '&start='+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1])+'-'+str(dateRange_full[i,2])+'T00:00&extent='+str(86400)+'&mag&geo&mlt&decl&sza&delta=none&baseline=none&station='+SuperMAGstations_names[stat]; #create full URL
                    else:
                        site_station = site_station_base + '&start='+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1])+'-'+str(dateRange_full[i,2])+'T00:00&extent='+str(86400)+'&mag&geo&mlt&decl&sza&baseline=all&station='+SuperMAGstations_names[stat]; #create full URL [not super clear what delta is so leaving it "default" by not mentioning it]
                    #END IF
                    print('In GRITI_import_SuperMAG_stations: Downloading data for station '+SuperMAGstations_names[stat]+' on '+str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+', it will take a few seconds prob.' );
                    try:
                        SuperMAGstations_web_raw = urlopen(site_station); #download the data needed
                        SuperMAGstations_web_rawRead = SuperMAGstations_web_raw.read(); #read it into a much more useful format
                        charset = SuperMAGstations_web_raw.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                        if( charset == None ):
                            charset = 'utf-8'; #set it to this
                        #END IF
                        SuperMAGstations_web_data = SuperMAGstations_web_rawRead.decode(charset).split('\n')[:-1]; #"decode" the HTML content so it's more legible
                    except:
                        import requests #requests doesn't seem to have cert errors, use it if urlopen is complainy
                        SuperMAGstations_web_raw = requests.get(site_station, stream=True); #get the data, stream is key I think
                        # SuperMAGstations_raw.encoding = SuperMAGstations_raw.apparent_encoding; #not needed, UTF-8 detected correctly
                        SuperMAGstations_web_data = SuperMAGstations_web_raw.text.split('\n')[:-1]; #decode to text
                    #END TRY
                    
                    SuperMAGstations_stationDict = {'geo':{},'magC':{}}; #prep it up, will keep stuff separate
                    SuperMAGstations_dataLen = len(SuperMAGstations_web_data); #precalc
                    #NOTE: the #s provided by the web api have 6 decimal places - float32 is easily good enough and will 1/2 memory
                    SuperMAGstations_stationDict['geo']['mag'+detrendedString] = np.empty( (SuperMAGstations_dataLen,3), dtype=np.float32); #nT, called 'geo' by SuperMAG [will be saved in order of N/E/Z (X/Y/Z)
                    SuperMAGstations_stationDict['magC']['mag'+detrendedString] = np.empty( (SuperMAGstations_dataLen,3), dtype=np.float32); #nT, called 'nez' by SuperMAG
                    SuperMAGstations_stationDict['geo']['lat'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #arcdeg, magnetic latitude is not as absolute as geo lat
                    SuperMAGstations_stationDict['geo']['long'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #arcdeg, magnetic latitude is not as absolute as geo lat
                    SuperMAGstations_stationDict['magC']['lat'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #arcdeg, magnetic latitude is not as absolute as geo lat
                    SuperMAGstations_stationDict['magC']['long'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #arcdeg, magnetic latitude is not as absolute as geo lat
                    SuperMAGstations_stationDict['magC']['MLT'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #"hours", MLT is hella not absolute
                    SuperMAGstations_stationDict['decl'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #sun stuff idk "declination from IGRF model" just for completeness
                    SuperMAGstations_stationDict['sza'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.float32); #sun stuff idk "solar zenith angle"
                    SuperMAGstations_stationDict['year'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.int16); #year
                    SuperMAGstations_stationDict['dayNum'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.int16); #daynum
                    SuperMAGstations_stationDict['sec'] = np.empty( (SuperMAGstations_dataLen,1), dtype=np.int32); #sec (wraps on a day), cuts down on int16 hr/min/sec to just int32 sec
                    for j in range(0,SuperMAGstations_dataLen): #gotta parse individual json disasters
                        if( j != 0 ):
                            tempDict = json.loads(SuperMAGstations_web_data[j][0:-1]);
                            if( tempDict['ext'] != SuperMAGstations_stationDict['data rate'] ):
                                #sanity check
                                print('ERROR in GRITI_import_SuperMAG_stations: Station '+SuperMAGstations_names[stat]+' on '+
                                      str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+
                                      ' returned a dict data rate '+tempDict['ext']+' but expected '+
                                      str(SuperMAGstations_stationDict['data rate'])+'. j = '+str(j)+'. Printing dict and Crashing!');
                                print(tempDict)
                                import sys
                                sys.crash(); #yeet
                            #END IF
                        else:
                            #start has extra annoying [ for some reason (chaos?)
                            tempDict = json.loads(SuperMAGstations_web_data[j][1:-1]);
                            SuperMAGstations_stationDict['data rate'] = tempDict['ext']; #record now, only need to record once
                        #END IF
                        if( tempDict['iaga'] != SuperMAGstations_names[stat] ):
                            #sanity check
                            print('ERROR in GRITI_import_SuperMAG_stations: Station '+SuperMAGstations_names[stat]+' on'+
                                  str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+
                                  ' returned a dict with station name '+tempDict['iaga']+'. j = '+str(j)+'. Printing dict and Crashing!');
                            print(tempDict)
                            import sys
                            sys.crash(); #yeet
                        #END IF
                        #save time
                        timeObj = datetime.datetime.utcfromtimestamp(tempDict['tval']); #time obj to deal with the timestamp (UTC)
                        SuperMAGstations_stationDict['year'][j] = np.int16(timeObj.strftime('%Y')); #year
                        SuperMAGstations_stationDict['dayNum'][j] = np.int16(timeObj.strftime('%j')); #dayNum
                        SuperMAGstations_stationDict['sec'][j] = np.int32(timeObj.strftime('%H'))*3600 + np.int32(timeObj.strftime('%M'))*60 + np.int32(timeObj.strftime('%S')); #sec
                        #save mag
                        SuperMAGstations_stationDict['geo']['mag'+detrendedString][j,0] = np.float32(tempDict['N']['geo']); #nT N/X
                        SuperMAGstations_stationDict['geo']['mag'+detrendedString][j,1] = np.float32(tempDict['E']['geo']); #nT E/Y
                        SuperMAGstations_stationDict['geo']['mag'+detrendedString][j,2] = np.float32(tempDict['Z']['geo']); #nT Z/Z (into earth)
                        #I'm in the gremlin zone for naming this 'mag''mag'
                        SuperMAGstations_stationDict['magC']['mag'+detrendedString][j,0] = np.float32(tempDict['N']['nez']); #nT N/X
                        SuperMAGstations_stationDict['magC']['mag'+detrendedString][j,1] = np.float32(tempDict['E']['nez']); #nT E/Y
                        SuperMAGstations_stationDict['magC']['mag'+detrendedString][j,2] = np.float32(tempDict['Z']['nez']); #nT Z/Z (into earth)
                        #save coords
                        SuperMAGstations_stationDict['geo']['lat'][j] = np.float32(tempDict['glat']); #arcdeg
                        SuperMAGstations_stationDict['geo']['long'][j] = np.float32(tempDict['glon']); #arcdeg
                        SuperMAGstations_stationDict['magC']['lat'][j] = np.float32(tempDict['mlat']); #arcdeg
                        SuperMAGstations_stationDict['magC']['long'][j] = np.float32(tempDict['mlon']); #arcdeg
                        SuperMAGstations_stationDict['magC']['MLT'][j] = np.float32(tempDict['mlt']); #"hours", may save a conversion
                        SuperMAGstations_stationDict['decl'][j] = np.float32(tempDict['decl']); #declination from IGRF model
                        SuperMAGstations_stationDict['sza'][j] = np.float32(tempDict['sza']); #solar zenith angle
                    #END FOR j
                    
                    SuperMAGstations_stationDict_keyz = list(SuperMAGstations_stationDict.keys()); #get em
                    # h5pyChunkShape = 1440; #a day's worth of data I guess            
                    if( os.path.isfile(SuperMAG_dataFilePath) ):
                        #append to the file if it does exist
                        with h5py.File(SuperMAG_dataFilePath, 'a') as h5file:
                            keyzThere = list(h5file.keys()); #get the group keys
                            if( SuperMAGstations_names[stat] in keyzThere ):
                                #if already there just need to append new mag type
                                #write in the geo data
                                siteGroup_sub = h5file[SuperMAGstations_names[stat]]['geo']; #alias it
                                dataHolder = siteGroup_sub.create_dataset('mag'+detrendedString, \
                                    SuperMAGstations_stationDict['geo']['mag'+detrendedString].shape, \
                                    dtype=SuperMAGstations_stationDict['geo']['mag'+detrendedString].dtype, \
                                    compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset  chunks=h5pyChunkShape, 
                                dataHolder[...] = SuperMAGstations_stationDict['geo']['mag'+detrendedString]; #write that data
                                #write in the mag data
                                siteGroup_sub = h5file[SuperMAGstations_names[stat]]['magC']; #alias it
                                dataHolder = siteGroup_sub.create_dataset('mag'+detrendedString, \
                                    SuperMAGstations_stationDict['magC']['mag'+detrendedString].shape, \
                                    dtype=SuperMAGstations_stationDict['magC']['mag'+detrendedString].dtype, \
                                    compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset  chunks=h5pyChunkShape, 
                                dataHolder[...] = SuperMAGstations_stationDict['magC']['mag'+detrendedString]; #write that data
                            else:
                                #if not there need to add entire new station group
                                siteGroup = h5file.create_group(SuperMAGstations_names[stat]); #create a group for each site
                                for j in range(0,len(SuperMAGstations_stationDict_keyz)):
                                    if( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]],dict) ):
                                        siteGroup_sub = siteGroup.create_group(SuperMAGstations_stationDict_keyz[j]); #create a group for sub-dicts
                                        SuperMAGstations_stationDict_subKeyz = list(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].keys()); #get em
                                        for k in range(0,len(SuperMAGstations_stationDict_subKeyz)):
                                            if( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]],np.ndarray) ):
                                                dataHolder = siteGroup_sub.create_dataset(SuperMAGstations_stationDict_subKeyz[k], \
                                                    SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]].shape, \
                                                    dtype=SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]].dtype, \
                                                    compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset  chunks=h5pyChunkShape, 
                                                dataHolder[...] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]; #write that data
                                            elif( np.isscalar(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]) ):
                                                siteGroup_sub.attrs[SuperMAGstations_stationDict_subKeyz[k]] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]; #record the attribute
                                            #END IF
                                        #END FOR k
                                    elif( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]],np.ndarray) ):
                                        dataHolder = siteGroup.create_dataset(SuperMAGstations_stationDict_keyz[j], \
                                            SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].shape, \
                                            dtype=SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].dtype, \
                                            compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset  chunks=h5pyChunkShape, 
                                        dataHolder[...] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]; #write that data
                                    elif( np.isscalar(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]) ):
                                        siteGroup.attrs[SuperMAGstations_stationDict_keyz[j]] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]; #record the attribute
                                    #END IF
                                #END FOR j
                            #END IF
                            # appendNum = h5file.attrs['append#']; #get append num
                            # appendNum += 1; #increment
                            # h5file.attrs['append#'] = appendNum; #reapply
                        #END WITH
                        # if( appendNum > 5 ):
                        #     pass; #nothing yet, future possible plan is to resave entire hdf5 file if it seems appending is shoddy
                        # #END IF
                    else:
                        #create the file if it doesn't exist
                        with h5py.File(SuperMAG_dataFilePath, 'w') as h5file:
                            siteGroup = h5file.create_group(SuperMAGstations_names[stat]); #create a group for each site
                            for j in range(0,len(SuperMAGstations_stationDict_keyz)):
                                if( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]],dict) ):
                                    siteGroup_sub = siteGroup.create_group(SuperMAGstations_stationDict_keyz[j]); #create a group for sub-dicts
                                    SuperMAGstations_stationDict_subKeyz = list(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].keys()); #get em
                                    for k in range(0,len(SuperMAGstations_stationDict_subKeyz)):
                                        if( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]],np.ndarray) ):
                                            dataHolder = siteGroup_sub.create_dataset(SuperMAGstations_stationDict_subKeyz[k], \
                                                SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]].shape, \
                                                dtype=SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]].dtype, \
                                                compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset  chunks=h5pyChunkShape, 
                                            dataHolder[...] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]; #write that data
                                        elif( np.isscalar(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]) ):
                                            siteGroup_sub.attrs[SuperMAGstations_stationDict_subKeyz[k]] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]; #record the attribute
                                        #END IF
                                    #END FOR k
                                elif( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]],np.ndarray) ):
                                    dataHolder = siteGroup.create_dataset(SuperMAGstations_stationDict_keyz[j], \
                                        SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].shape, \
                                        dtype=SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].dtype, \
                                        compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset  chunks=h5pyChunkShape, 
                                    dataHolder[...] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]; #write that data
                                elif( np.isscalar(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]) ):
                                    siteGroup.attrs[SuperMAGstations_stationDict_keyz[j]] = SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]; #record the attribute
                                #END IF
                            #END FOR j
                            #file-wide values
                            h5file.attrs['version'] = version_alg; #record the attribute
                            # h5file.attrs['append#'] = 0; #record the attribute [exists to enforce housekeeping]
                        #END WITH
                    #END IF
                    SuperMAGstations_stationDict['version'] = version_alg; #record into the dict no matter if append or write
                else:
                    print('ERROR in GRITI_import_SuperMAG_stations: Station '+SuperMAGstations_names[stat]+
                          ' is NOT available for '+str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+
                          '. Printing available stations, also see https://supermag.jhuapl.edu/mag/ for a map of stations. Crashing!');
                    print(SuperMAGstations_stations);
                    import sys
                    sys.crash(); #yeet
                #END IF
            #END IF
            
            #--- MECHA append code to append whatever was made or read off of the cached HDF5 files ---
            if( SuperMAGstations_names[stat] in SuperMAGstations_data ):
                SuperMAGstations_stationDict_keyz = list(SuperMAGstations_stationDict.keys()); #get em
                for j in range(0,len(SuperMAGstations_stationDict_keyz)):
                    if( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]],dict) ):
                        SuperMAGstations_stationDict_subKeyz = list(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].keys()); #get em
                        for k in range(0,len(SuperMAGstations_stationDict_subKeyz)):
                            if( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]],np.ndarray) ):
                                if( SuperMAGstations_stationDict_subKeyz[k] in SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]] ):
                                    #do array appending
                                    SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]] = \
                                        np.concatenate((SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]],SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]),axis=0); #append on, no copy needed b/c def changed a lot
                                else:
                                    #load it in new
                                    SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]] = \
                                        SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]].copy(); #load it in and copy it
                                #END IF
                            elif( np.isscalar(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]) ):
                                #check if same
                                if( SuperMAGstations_stationDict_subKeyz[k] in SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]] ):
                                    if( SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]] != \
                                       SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]] ):
                                        #if not same, bad times :(
                                        print('WARNING in GRITI_import_SuperMAG_stations: Scalar value '+SuperMAGstations_stationDict_subKeyz[k]+\
                                              ' in '+SuperMAGstations_stationDict_keyz[j]+' dict '+\
                                              ' was '+SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]+\
                                              ' but now is '+SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]+\
                                              ' on '+str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+'. Doing nothing, just warning.');
                                    #END IF
                                else:
                                    #load it in new
                                    SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]] = \
                                        copy.deepcopy(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]][SuperMAGstations_stationDict_subKeyz[k]]); #load it in and copy it (may not be required for python scalars but w/e small fry)
                                #END IF
                            #END IF
                        #END FOR k
                    elif( isinstance(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]],np.ndarray) ):
                        if( SuperMAGstations_stationDict_keyz[j] in SuperMAGstations_data[SuperMAGstations_names[stat]] ):
                            #do array appending
                            SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]] = \
                                np.concatenate((SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]],SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]),axis=0); #append on, no copy needed b/c def changed a lot
                        else:
                            #load it in new
                            SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]] = \
                                SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]].copy(); #load it in and copy it
                        #END IF
                    elif( np.isscalar(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]) ):
                        #check if same
                        if( SuperMAGstations_stationDict_keyz[j] in SuperMAGstations_data[SuperMAGstations_names[stat]] ):
                            if( SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]] != SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]] ):
                                #if not same, bad times :(
                                print('WARNING in GRITI_import_SuperMAG_stations: Scalar value '+SuperMAGstations_stationDict_keyz[j]+\
                                      ' was '+SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]]+\
                                      ' but now is '+SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]+\
                                      ' on '+str(dateRange_dayNum_full[i,0])+'/'+str(dateRange_dayNum_full[i,1])+'. Doing nothing, just warning.');
                            #END IF
                        else:
                            #load it in new
                            SuperMAGstations_data[SuperMAGstations_names[stat]][SuperMAGstations_stationDict_keyz[j]] = \
                                copy.deepcopy(SuperMAGstations_stationDict[SuperMAGstations_stationDict_keyz[j]]); #load it in and copy it (may not be required for python scalars but w/e small fry)
                        #END IF
                    #END IF
                #END FOR j
            else:
                SuperMAGstations_data[SuperMAGstations_names[stat]] = copy.deepcopy(SuperMAGstations_stationDict); #load it in! so easy
            #END IF
        #END FOR i
    #END FOR stat
    
    #get out the needed type & condense lat/long
    SuperMAG_station_keys = list(SuperMAGstations_data.keys()); #get the keys
    for j in range(0,len(SuperMAG_station_keys)):
        SuperMAG_coord_keys = list(SuperMAGstations_data[SuperMAG_station_keys[j]][coordType].keys()); #get the coord keys
        for k in range(0,len(SuperMAG_coord_keys)):
            SuperMAGstations_data[SuperMAG_station_keys[j]][SuperMAG_coord_keys[k]] = SuperMAGstations_data[SuperMAG_station_keys[j]][coordType][SuperMAG_coord_keys[k]].copy(); #bring out everything in the sub-dict
        #END FOR k
        SuperMAGstations_data[SuperMAG_station_keys[j]].pop('geo', None); #remove geo and mag sub-dicts b/c pulled out
        SuperMAGstations_data[SuperMAG_station_keys[j]].pop('magC', None); #this is automatically overrun
        if( FLG_detrended == True ):
            SuperMAGstations_data[SuperMAG_station_keys[j]]['mag'] = copy.deepcopy(SuperMAGstations_data[SuperMAG_station_keys[j]]['mag'+detrendedString]); #overwrite b/c want detrended data
            SuperMAGstations_data[SuperMAG_station_keys[j]].pop('mag'+detrendedString, None); #remove excess
        #END IF
        
        SuperMAGstations_data[SuperMAG_station_keys[j]]['lat'] = np.median(SuperMAGstations_data[SuperMAG_station_keys[j]]['lat']); #get the median instead of a big array [like this b/c mag lat/long can move]
        SuperMAGstations_data[SuperMAG_station_keys[j]]['long'] = np.median(SuperMAGstations_data[SuperMAG_station_keys[j]]['lat']); #get the median instead of a big array
    #END FOR j
    
    #Do some minor time keeping stuff - just use last site I guess?
    SuperMAGstations_data['data rate'] = SuperMAGstations_data[SuperMAGstations_names[stat]]['data rate']; #get the time step
    SuperMAGstations_data['time'] = np.int64(dateRange_dayNum_full[0,1])*86400 + np.arange(0,dateRange_dayNum_full.shape[0]*86400,SuperMAGstations_data['data rate']); #time in sec, create in case there are data gaps I guess
    SuperMAGstations_data['time unique'] = SuperMAGstations_data['time']; #alias
    
    if( np.isclose(np.mod(SuperMAGstations_data['data rate'],1),0) == True ):
        SuperMAGstations_data['data rate'] = np.int64(SuperMAGstations_data['data rate']); #make it an integer if its an integer
    #END IF
    
    #NaN 999999 values
    for j in range(0,len(SuperMAG_station_keys)):
        SuperMAG_data_keys = list(SuperMAGstations_data[SuperMAG_station_keys[j]].keys());
        for k in range(0,len(SuperMAG_data_keys)):
            kj = np.isclose(SuperMAGstations_data[SuperMAG_station_keys[j]][SuperMAG_data_keys[k]],999999); #get where 999999 are at
            if( np.any(kj) ):
                if( SuperMAGstations_data[SuperMAG_station_keys[j]][SuperMAG_data_keys[k]].dtype == np.int64() ): #can't nan int64
                    SuperMAGstations_data[SuperMAG_station_keys[j]][SuperMAG_data_keys[k]] = np.float64(SuperMAGstations_data[SuperMAG_data_keys[j]][SuperMAG_data_keys[k]]); #put into float64 so we can NaN it
                #END IF
                SuperMAGstations_data[SuperMAG_station_keys[j]][SuperMAG_data_keys[k]][kj] = np.nan; #nan it out
            #END IF
        #END FOR k
    #END FOR j
    
    return SuperMAGstations_data #finish it out