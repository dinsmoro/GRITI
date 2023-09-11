#Function to import SuperMAG data from the pre-downloaded data
#RD on 12/10/2018

'''
GOAL: Get SuperMAG data
expecting: 
 dateRange_dayNum_full - [year,day#min;year,day#inbetween;year,day#max] numerical format. So, say [2013,126;2013,127;2013,128] for 2013 Day 126 to Day 128. No chars pls
 folder - list with folder[0] being the running folder (where the code is) and folder[1] being where the data is kept (can be the same, or not)
NOTE: SuperMAG is every 1 minute.
NOTE: Download year-long indices from https://supermag.jhuapl.edu/indices/ & read descriptions of indices

'''

import numpy as np
import h5py
import os
from Code.subfun_strfind import strfind

def GRITI_import_SuperMAG(dateRange_dayNum_full, settings_paths, FLG_deleteOrig = False):
    
    
    dateRange_yearRange = np.arange(dateRange_dayNum_full[0,0],dateRange_dayNum_full[-1,0]+1,1,dtype=np.int16); #get the full year range from min to max
    
    SuperMAG_data = {}; #prep dict
    
    #Read in the data
    for i in range(0,dateRange_yearRange.size): #steps through each year
        #check to see if data is already processed
        SuperMAG_dataFilePath = os.path.join(settings_paths['data'], 'SuperMAG', 'Indices', 'supermag-indices-' + str(dateRange_yearRange[i]) + '.h5'); #create the expected data path
        if( os.path.isfile(SuperMAG_dataFilePath) ): #check if SuperMAG data file exists 
            with h5py.File(SuperMAG_dataFilePath, 'r') as h5file:
                #----- Read the data types in -----
                #--- First prep the mask and import needed vars for it ---
                SuperMAG_data_raw_keys = list(h5file.keys()); #get the saved keys
                SuperMAG_data_raw = {}; #prep a dict
                for j in range(0,len(SuperMAG_data_raw_keys)):
                    SuperMAG_data_raw[SuperMAG_data_raw_keys[j]] = h5file.get(SuperMAG_data_raw_keys[j])[()]; #get that dataset out
                #END FOR j
            #END WITH
        else:
            print('WARNING in GRITI_import_SuperMAG: .h5 file for year '+ str(dateRange_yearRange[i]) +' isn\'t available, converting from CSV. This will take a bit (~15 sec) b/c of a string date thing.');
            import time
            tic = time.time(); #time it
            SuperMAG_dataFilePath_CSV = os.path.join(settings_paths['data'], 'SuperMAG', 'Indices', 'supermag-indices-' + str(dateRange_yearRange[i]) + '.csv'); #create the expected data path
            if( os.path.isfile(SuperMAG_dataFilePath_CSV) ): #check if SuperMAG data file exists   
                import csv
                from Code.subfun_strstr import strstrNB as strstr
                from Code.subfun_date_to_dayNum import subfun_date_to_dayNum #these are only needed to convert CSV -> hdf5
                try: #gonna try to read the file - if we fail, it's borked or something
                    # SuperMAG_data_raw = pandas.read_csv(SuperMAG_dataFilePath_CSV, delimiter=','); #fast but new dependency not worth it b/c of time conversion reqs
                    #more-or-less equivalent with built-in python
                    SuperMAG_data_raw = {}; #prep a dict
                    with open(SuperMAG_dataFilePath_CSV,'r') as csvFile:
                        csvReader = csv.reader(csvFile); #make a file reader object b/c object oriented
                        header = next(csvReader); #read the header
                        data_raw = [row for row in csvReader]; #read the data
                        for j in range(0,len(header)):
                            SuperMAG_data_raw[header[j]] = [None for k in range(0,len(data_raw))]; #preallocate
                            for k in range(0,len(data_raw)):
                                SuperMAG_data_raw[header[j]][k] = data_raw[k][j]; #load in the data nice and slow b/c list
                            #END FOR k
                        #END FOR j
                        del data_raw, header
                    #END WITH
                except OSError:
                    #if this happens then the data is corrupted or something [happens implicitly]
                    print('ERROR in GRITI_import_SuperMAG: "'+SuperMAG_dataFilePath_CSV+'" is broken. Redownload it from SuperMAG for the year '+str(dateRange_yearRange[i])+'. Returning gibberish ¯\_(ツ)_/¯');
                    return 'no'
                #END TRY
            else:
                print('ERROR in GRITI_import_SuperMAG: "'+SuperMAG_dataFilePath_CSV+'" is not downloaded. Download it from SuperMAG for the year '+str(dateRange_yearRange[i])+'. Returning gibberish ¯\_(ツ)_/¯');
                return 'no'
            #END IF
            
            SuperMAG_data_raw_keys = list(SuperMAG_data_raw.keys()); #get the keys
            #Convert lists to numpy arrays
            for j in range(0,len(SuperMAG_data_raw_keys)):
                if('Date_UTC' != SuperMAG_data_raw_keys[j]):
                    if( strfind(SuperMAG_data_raw[SuperMAG_data_raw_keys[j]],'.',1) > 0 ):
                        SuperMAG_data_raw[SuperMAG_data_raw_keys[j]] = np.asarray(SuperMAG_data_raw[SuperMAG_data_raw_keys[j]],dtype=np.float64); #convert to numpy array
                    else:
                        SuperMAG_data_raw[SuperMAG_data_raw_keys[j]] = np.asarray(SuperMAG_data_raw[SuperMAG_data_raw_keys[j]],dtype=np.int64); #convert to numpy array
                    #END IF
                #END IF
            #END FOR j
            
            #Deal with Date_UTC being a useless string
            SuperMAG_data_raw['time'] = np.empty(len(SuperMAG_data_raw['Date_UTC']), dtype=np.int64); #preallocate
            SuperMAG_data_raw['hr'] = np.empty(len(SuperMAG_data_raw['Date_UTC']), dtype=np.int16); #preallocate
            SuperMAG_data_raw['min'] = np.empty(len(SuperMAG_data_raw['Date_UTC']), dtype=np.int16); #preallocate
            SuperMAG_data_raw['sec'] = np.empty(len(SuperMAG_data_raw['Date_UTC']), dtype=np.int16); #preallocate
            SuperMAG_data_raw_dateHolder = np.empty((len(SuperMAG_data_raw['Date_UTC']),3), dtype=np.int16); #preallocate
            loc_hyphens = strstr(SuperMAG_data_raw['Date_UTC'][0],'-'); #constant spacing lets us do this once
            loc_space = SuperMAG_data_raw['Date_UTC'][0].find(' ');
            loc_colons = strstr(SuperMAG_data_raw['Date_UTC'][0],':');
            for j in range(0,len(SuperMAG_data_raw['Date_UTC'])):
                SuperMAG_data_raw_dateHolder[j,:] = (int(SuperMAG_data_raw['Date_UTC'][j][:loc_hyphens[0]]),int(SuperMAG_data_raw['Date_UTC'][j][loc_hyphens[0]+1:loc_hyphens[1]]),int(SuperMAG_data_raw['Date_UTC'][j][loc_hyphens[1]+1:loc_space])); #load it in
            #END FOR j
            SuperMAG_data_raw_dateHolder = subfun_date_to_dayNum( SuperMAG_data_raw_dateHolder ); #convert to dayNum for this date
            SuperMAG_data_raw['year'] = SuperMAG_data_raw_dateHolder[:,0]; #record year
            SuperMAG_data_raw['dayNum'] = SuperMAG_data_raw_dateHolder[:,1]; #record year
            SuperMAG_data_raw_dateHolder = np.int64(SuperMAG_data_raw_dateHolder[:,1])*86400; #convert dayNum to seconds, ignore year yolo
            for j in range(0,len(SuperMAG_data_raw['Date_UTC'])):
                SuperMAG_data_raw['time'][j] = SuperMAG_data_raw_dateHolder[j] + np.int64(SuperMAG_data_raw['Date_UTC'][j][loc_space+1:loc_colons[0]])*3600 + np.int64(SuperMAG_data_raw['Date_UTC'][j][loc_colons[0]+1:loc_colons[1]])*60 + np.int64(SuperMAG_data_raw['Date_UTC'][j][loc_colons[1]+1:]); #get the time in seconds
                SuperMAG_data_raw['hr'][j] = SuperMAG_data_raw['Date_UTC'][j][loc_space+1:loc_colons[0]];
                SuperMAG_data_raw['min'][j] = SuperMAG_data_raw['Date_UTC'][j][loc_colons[0]+1:loc_colons[1]];
                SuperMAG_data_raw['sec'][j] = SuperMAG_data_raw['Date_UTC'][j][loc_colons[1]+1:];
            #END FOR j
            del SuperMAG_data_raw['Date_UTC']; #done with Date_UTC now
            
            #Save as .h5 file
            SuperMAG_data_raw_keys = list(SuperMAG_data_raw.keys()); #get the keys
            h5pyChunkShape = SuperMAG_data_raw['time'].shape; #get the shape of one of the vectors and use it as a chunk (read only whole chunks)
            with h5py.File(SuperMAG_dataFilePath, 'w') as h5file: #, rdcc_nbytes =500*1024*1024
                for j in range(0,len(SuperMAG_data_raw_keys)):
                    h5file.create_dataset(SuperMAG_data_raw_keys[j], data=SuperMAG_data_raw[SuperMAG_data_raw_keys[j]], chunks=h5pyChunkShape, compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #write that data
                #END FOR j
            #END WITH
            
            if( FLG_deleteOrig == True ):
                os.remove(SuperMAG_dataFilePath_CSV); #only delete if req'd
            #END IF
            
            toc = time.time() - tic; #time it
            print('WARNING in GRITI_import_SuperMAG: Conversion to .h5 from .csv for year '+ str(dateRange_yearRange[i])+' took '+str(np.round(toc,2))+' sec.');
        #END IF
        
        #Load in and append on the data
        SuperMAG_data_raw_keys = list(SuperMAG_data_raw.keys()); #get the keys
        for j in range(0,len(SuperMAG_data_raw_keys)):
            if( SuperMAG_data_raw_keys[j] in SuperMAG_data):
                SuperMAG_data[SuperMAG_data_raw_keys[j]] = np.append(SuperMAG_data[SuperMAG_data_raw_keys[j]],SuperMAG_data_raw[SuperMAG_data_raw_keys[j]]); #add on data
            else:
                SuperMAG_data[SuperMAG_data_raw_keys[j]] = SuperMAG_data_raw[SuperMAG_data_raw_keys[j]]; #create key
            #END IF                
        #END FOR j
    #END FOR i
    
    #Do some minor time keeping stuff
    SuperMAG_data['time unique'] = SuperMAG_data['time']; #alias
    SuperMAG_data['data rate'] = np.median(np.diff(SuperMAG_data['time unique'])); #get the time step
    if( np.isclose(np.mod(SuperMAG_data['data rate'],1),0) == True ):
        SuperMAG_data['data rate'] = np.int64(SuperMAG_data['data rate']); #make it an integer if its an integer
    #END IF
    
    #NaN 999999 values
    SuperMAG_data_keys = list(SuperMAG_data.keys()); #get the keys
    for j in range(0,len(SuperMAG_data_keys)):
        kj = np.isclose(SuperMAG_data[SuperMAG_data_keys[j]],999999); #get where 999999 are at
        if( np.any(kj) ):
            if( SuperMAG_data[SuperMAG_data_keys[j]].dtype == np.int64() ): #can't nan int64
                SuperMAG_data[SuperMAG_data_keys[j]] = np.float64(SuperMAG_data[SuperMAG_data_keys[j]]); #put into float64 so we can NaN it
            #END IF
            SuperMAG_data[SuperMAG_data_keys[j]][kj] = np.nan; #nan it out
        #END IF
    #END FOR j
    
    #Make SMLrXX time stamps from SMErXX and SMUrXX
    # SME = SMU-SML, SML = SMU-SME
    try:
        for i in range(0,24):
            SuperMAG_data['SMLr'+str(i).zfill(2)] = SuperMAG_data['SMUr'+str(i).zfill(2)] - SuperMAG_data['SMEr'+str(i).zfill(2)]; #calculate it
        #END FOR i
    except:
        print('WARNING: SMUrXX or SMErXX (MLT zone SMU/SME) are missing. SMLrXX will not be calculated and also those won\'t exist. If used code will crash. Redownload from site with everything.');
    #END TRY
    
    return SuperMAG_data #finish it out