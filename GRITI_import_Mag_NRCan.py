#Function to import Mag NRCan data from the internet
#RD on 12/10/2018

'''
GOAL: Get Magnetometer data from NRCan
#https://www.geomag.nrcan.gc.ca/data-donnee/sd-en.php
expecting: 
 dates dictionary
 settings dictonary
 FLG_overwrite [optional, default 0] - 1 downloads data again and overwrites old data, 0 just reads previously saved local data. Has no effect on new data.
NOTE: Mag is every 1 sec
NOTE: Mag data is the 'sec' mag data in the 'variation' format
NOTE: Processed Mag data goes in a .h5 file in Data/Magnetometer/[Year]/
NOTE: Raw Mag data goes in Data/Magnetometer/[Year]/NRCan/
NOTE: Raw Mag data can be set to be deleted, default is it is kept
'''

import numpy as np
import h5py
import os
import glob
import gzip #gzip is built into default python
from subfun_strfind import strfind

def GRITI_import_Mag_NRCan(dates, settings, FLG_deleteRaw = 0, FLG_overwrite = 0):
    version = 1.0; #import Mag NRCan version
    
    dataPath = settings['paths']['data']+'\\Magnetometer\\'; #prep the data path
    dateRange_full = dates['date range full']; #get the full date range in daynum format out
    dateRange_dayNum_full = dates['date range full dayNum']; #get the full date range in daynum format out
    
    #make sure folder exists
    if( os.path.isdir(dataPath) == False ):
        os.makedirs(dataPath); #make the path
    #END IF
    
    FLG_dataAlreadyProcessed = 0; #flag to show if data is already downloaded or not
    for i in range(0,dateRange_dayNum_full.shape[0]): #steps through each day
        
        #check to see if data is already processed
        dataPath_currFile = dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan_'+str(dateRange_dayNum_full[i,1])+'.h5'; #create the expected data path 
        #make sure folder exists
        if( os.path.isdir(dataPath + str(dateRange_dayNum_full[i,0])) == False ):
            os.makedirs(dataPath + str(dateRange_dayNum_full[i,0])); #make the path
        #END IF
        
        if( os.path.isfile(dataPath_currFile) == 1 ): #check if omni data file exists            
            try: #gonna try to read the file - if we fail, it's borked or something
                with h5py.File(dataPath_currFile, 'r') as filez:
                    Mag_dataTemp_test = filez["/siteNames"][:].astype("S3"); #try to read some data
                    file_version = filez.attrs['version']; #read the attribute
    #                Mag_folders = list(filez.keys()); #lists the data folders in the HDF5 file
                #END WITH
                if( file_version == version ):
                    FLG_dataAlreadyProcessed = 1; #set the flag to not process
                elif( file_version < version):
                    print('WARNING: MAGNETOMETER NRCAN FILE VERSION IS OLD '+str(file_version)+', CURRENT ALG VERSION IS '+str(version)+'. OVERWRITING CURRENT DATA. ON '+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1]).zfill(2)+'-'+str(dateRange_full[i,2]).zfill(2));
                    FLG_dataAlreadyProcessed = 0; #set the flag to process [old file version]
                else:
                    print('WARNING: MAGNETOMETER NRCAN FILE VERSION '+str(file_version)+' IS NEWER THAN THE ALG VERSION '+str(version)+'. CRASHING AND QUITTING. ON '+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1]).zfill(2)+'-'+str(dateRange_full[i,2]).zfill(2));
                    import sys
                    sys.crash(); #crash it
                #END IF
            except:
                #if this happens then the data is corrupted or something
                FLG_dataAlreadyProcessed = 0; #set the flag to process
            #END TRY
        else:
            #file not there, process it
            FLG_dataAlreadyProcessed = 0; #set the flag to process
        #END IF
        
        if( (FLG_dataAlreadyProcessed != 1) | (FLG_overwrite == 1) ):
            if( FLG_overwrite == 1 ):
                print('WARNING: MAGNETOMETER NRCAN OVERWRITE IS ON. OVERWRITING CURRENT DATA. ON '+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1]).zfill(2)+'-'+str(dateRange_full[i,2]).zfill(2));
            #END IF
            print('=== Processing date '+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1]).zfill(2)+'-'+str(dateRange_full[i,2]).zfill(2)+' ('+str(dateRange_dayNum_full[i,0])+'-'+str(dateRange_dayNum_full[i,1])+')');
            import pandas
            from subfun_strstr import strstr #this one is aware of \n's
            from subfun_strstrNB import strstrNB #this one ignores \n's
            #Process the data into a nice HDF5 file        
            siteNames = np.array(['BLC','BRD','CBB','EUA','FCC','IQA','MEA','OTT','RES','SNK','STJ','VIC','YKC']); #names of the magnetometer stations
            siteData = glob.glob(dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan\\*'+str(dateRange_full[i,0])+str(dateRange_full[i,1]).zfill(2)+str(dateRange_full[i,2]).zfill(2)+'vsec.sec.gz')
            # siteNames = np.array( np.zeros(len(siteData)), dtype='<U3' ); #preallocate
            for j in range(0,len(siteData)): 
                siteNames[j] = siteData[j][siteData[j].find('NRCan')+6:siteData[j].rfind(str(dateRange_full[i,0]))].upper(); #get the site names from the data
            #END FOR j
            siteNamesFound = np.zeros(len(siteNames),dtype=bool); #prep array to record which were missing
            for j in range(0,len(siteNames)):
                siteFound_temp = strfind(siteData,siteNames[j].lower()); #search the site data for each site name
                if( siteFound_temp.sum() == 1 ):
                    siteNamesFound = siteNamesFound | siteFound_temp.astype(bool); #increment a found site match
                elif( siteFound_temp.sum() > 1 ):
                    import sys
                    print('ERROR: MULTIPLE SITE MATCHES FOUND FOR \"'+siteNames[j].lower()+'\". PRINTING LIST OF SITE FILES:\n');
                    print(siteData);
                    print('QUITTING. FIX THIS ISSUE. SHOULDN\'T HAPPEN EVERR');
                    sys.crash(); #crash the program
                #END IF
            #END FOR j
            if( siteNamesFound.sum() != len(siteNames) ):
                #Sites missing, go to the web to fill in the data
                #USE THIS https://stackoverflow.com/questions/3451111/unzipping-files-in-python
                #ALSO THIS https://stackoverflow.com/questions/30905198/python-to-select-a-dropdown-menu-on-a-website-and-submit
                print('HAVEN\'T DONE THIS YET, GOTTA READ A COMPLICATED BROWSER THINGY. QUITTING GOOD LUCK. CAN ALSO MANUALLY DOWNLOAD FROM\n'+ \
                      'https://www.geomag.nrcan.gc.ca/data-donnee/dl/dl-en.php\nDAY IN QUESTION IS '+ \
                      str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1]).zfill(2)+'-'+str(dateRange_full[i,2]).zfill(2)+' AND FULL DATE RANGE:\n');
                print(dateRange_full);
                sys.crash(); #crash the program
                #make sure folder exists
                if( os.path.isdir(dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan') == False ):
                    os.makedirs(dataPath + str(dateRange_dayNum_full[i,0]))+ '\\NRCan'; #make the path
                #END IF
            #END IF
            
            #Prep mag data
            siteNames = siteNames[siteNamesFound]; #make sure only to record the available site names
            Mag_dataTemp = {'site names':siteNames}; #use a dictionary to keep all the data
            #Sites expected are there, go for it
            for j in range(0,len(siteNames)):
                Mag_dataTemp[siteNames[j]] = {}; #prep a dict for the site
                siteFound_temp = np.where(strfind(siteData,siteNames[j].lower()))[0].item(); #search the site data for each site name
                
                with gzip.open(siteData[siteFound_temp],'rt') as filez:
                    Mag_dataTemp_raw = filez.read(); #read that file
                #END WITH
                with open(dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan\\temp.sec', 'w') as filez:
                    filez.write(Mag_dataTemp_raw); #save as a sec?
                #END WITH
                
                #FIND HEADER LINE START
                headerLineDataStart = Mag_dataTemp_raw.rfind('|'); #gonna use it to speed up searches
                headerLine = strstr(Mag_dataTemp_raw[0:headerLineDataStart],r'\n').size; #counts how many /n's occur before the data starts
                #Find geodedic latitude
                tempStart = strstrNB(Mag_dataTemp_raw[0:headerLineDataStart],'Latitude').item(); #get where geodedic latitude is delcated
                tempLine = Mag_dataTemp_raw[tempStart:tempStart+Mag_dataTemp_raw[tempStart:headerLineDataStart].find('|')]; #get the line (each line ends with a |)
                Mag_dataTemp[siteNames[j]]['lat'] = float(tempLine.split()[1]); #record it
                #Find geodedic longitude
                tempStart = strstrNB(Mag_dataTemp_raw[0:headerLineDataStart],'Longitude').item(); #get where geodedic latitude is delcated
                tempLine = Mag_dataTemp_raw[tempStart:tempStart+Mag_dataTemp_raw[tempStart:headerLineDataStart].find('|')]; #get the line (each line ends with a |)
                Mag_dataTemp[siteNames[j]]['long'] = float(tempLine.split()[1]); #record it
                #Find elevation
                tempStart = strstrNB(Mag_dataTemp_raw[0:headerLineDataStart],'Elevation').item(); #get where geodedic latitude is delcated
                tempLine = Mag_dataTemp_raw[tempStart:tempStart+Mag_dataTemp_raw[tempStart:headerLineDataStart].find('|')]; #get the line (each line ends with a |)
                Mag_dataTemp[siteNames[j]]['elv'] = float(tempLine.split()[1]); #record it
                #Find data rate
                tempStart = strstr(Mag_dataTemp_raw[0:headerLineDataStart],'1-second').size
                if( tempStart == 1 ):
                    Mag_dataTemp[siteNames[j]]['dataRate'] = 1; #seconds, record
                    Mag_dataTemp[siteNames[j]]['dataRateF'] = 5; #seconds, record
                else:
                    #if the 1-second phrase wasn't found it's 1 min
                    Mag_dataTemp[siteNames[j]]['dataRate'] = 60; #seconds, record
                    Mag_dataTemp[siteNames[j]]['dataRateF'] = 60; #seconds, record
                #END IF
                
                #Read the actual data
                Mag_dataTemp_raw = pandas.read_csv(dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan\\temp.sec',delim_whitespace=True,header=headerLine); #Read the data quick
                
                temp_time = Mag_dataTemp_raw.iloc[:,1].to_numpy(dtype=object, copy=True); #get the time
                temp_time = np.array([jk.split(':') for jk in temp_time.ravel()],dtype=object); #split by : [from https://stackoverflow.com/questions/50459578/split-2d-numpy-array-of-strings-on-character]
                Mag_timeUnique_sec = temp_time[:,0].astype(np.int32)*3600 + temp_time[:,1].astype(np.int32)*60 + np.int32(np.round(temp_time[:,2].astype(np.float64))); #get the total seconds in the day
                Mag_mag = Mag_dataTemp_raw.iloc[:,3:6].to_numpy(dtype=np.float32, copy=True); #mag data, [X, Y, Z]
                Mag_magF = Mag_dataTemp_raw.iloc[:,6].to_numpy(dtype=np.float32, copy=True); #mag Fdata, [F is every 5 sec]
                k = np.isclose(Mag_magF,88888); #get where it's 88888 [means no data]
                Mag_magF = np.delete(Mag_magF,k); #remove the useless data
                Mag_timeUnique_secF = np.delete(Mag_timeUnique_sec,k); #get the time unique to match
                
                Mag_timeUnique_year = np.tile(dateRange_dayNum_full[i,0], (Mag_dataTemp_raw.shape[0])); #get year (it's always the same)
                Mag_timeUnique_dayNum = np.tile(dateRange_dayNum_full[i,1], (Mag_dataTemp_raw.shape[0])); #get daynum (it's always the same)
    
                #Write it into the dictionary
                Mag_dataTemp[siteNames[j]]['mag'] = Mag_mag; #the X, Y, Z magnetometer data
                Mag_dataTemp[siteNames[j]]['magF'] = Mag_magF; #the magnitude 'F' magnetometer data
                Mag_dataTemp[siteNames[j]]['sec'] = Mag_timeUnique_sec; #storing the seconds is more accurate than a float
                Mag_dataTemp[siteNames[j]]['secF'] = Mag_timeUnique_secF; #storing the seconds is more accurate than a float
                Mag_dataTemp[siteNames[j]]['year'] = Mag_timeUnique_year; #storing the seconds is more accurate than a float
                Mag_dataTemp[siteNames[j]]['yearF'] = np.delete(Mag_timeUnique_year,k); #storing the seconds is more accurate than a float
                Mag_dataTemp[siteNames[j]]['dayNum'] = Mag_timeUnique_dayNum; #storing the seconds is more accurate than a float
                Mag_dataTemp[siteNames[j]]['dayNumF'] = np.delete(Mag_timeUnique_dayNum,k); #storing the seconds is more accurate than a float
                
                #Remove major outliers
                if( np.all(np.isclose(np.diff(Mag_dataTemp[siteNames[j]]['magF']),0)) == True ):
                    #this case is for no magF data, but there's mag data. Can calculate averages for this
                    kk = np.concatenate( (np.where(~k)[0],np.array( (k.size+1,)) ) ); #get indexes where the avg pts are in the big mag vector
                    for k in range(0,kk.size-1):
                        Mag_dataTemp[siteNames[j]]['magF'][k] = np.sqrt( \
                            np.nanmean(Mag_dataTemp[siteNames[j]]['mag'][kk[k]:kk[k+1],0])**2 + \
                            np.nanmean(Mag_dataTemp[siteNames[j]]['mag'][kk[k]:kk[k+1],1])**2 + \
                            np.nanmean(Mag_dataTemp[siteNames[j]]['mag'][kk[k]:kk[k+1],2])**2 ); #this isn't quite what they do, but it's pretty close (basic magnitude w/ average in time)
                    #END FOR k
                #END IF
                k = np.abs(Mag_dataTemp[siteNames[j]]['magF'] - np.nanmedian(Mag_dataTemp[siteNames[j]]['magF']))/np.nanmedian(np.abs(Mag_dataTemp[siteNames[j]]['magF'] - np.nanmedian(Mag_dataTemp[siteNames[j]]['magF']))) > 100; #gets extreme outliers [which are bad data points]
                Mag_dataTemp[siteNames[j]]['magF'][k] = np.interp(Mag_dataTemp[siteNames[j]]['secF'][k], Mag_dataTemp[siteNames[j]]['secF'][~k], Mag_dataTemp[siteNames[j]]['magF'][~k]); #linearly interpolate between that point
            #END FOR j
            
            #write all this stuff into a nice HDF5 file
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(dataPath_currFile, 'w') as filez:
                siteList = filez.create_dataset('siteNames', (len(Mag_dataTemp['site names']),), dtype='S3' ,compression="gzip"); #create dataset
                siteList[...] = Mag_dataTemp['site names'].astype('S3'); #write that data
                for j in range(0,len(Mag_dataTemp['site names'])):
                    siteGroup = filez.create_group(Mag_dataTemp['site names'][j]); #create a group for each site
                    dataHolder = siteGroup.create_dataset('mag', Mag_dataTemp[Mag_dataTemp['site names'][j]]['mag'].shape, dtype='float32' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['mag']; #write that data
                    dataHolder = siteGroup.create_dataset('magF', Mag_dataTemp[Mag_dataTemp['site names'][j]]['magF'].shape, dtype='float32' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['magF']; #write that data
                    dataHolder = siteGroup.create_dataset('sec', Mag_dataTemp[Mag_dataTemp['site names'][j]]['sec'].shape, dtype='int32' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['sec']; #write that data
                    dataHolder = siteGroup.create_dataset('secF', Mag_dataTemp[Mag_dataTemp['site names'][j]]['secF'].shape, dtype='int32' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['secF']; #write that data
                    dataHolder = siteGroup.create_dataset('year', Mag_dataTemp[Mag_dataTemp['site names'][j]]['year'].shape, dtype='int16' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['year']; #write that data
                    dataHolder = siteGroup.create_dataset('yearF', Mag_dataTemp[Mag_dataTemp['site names'][j]]['yearF'].shape, dtype='int16' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['yearF']; #write that data
                    dataHolder = siteGroup.create_dataset('dayNum', Mag_dataTemp[Mag_dataTemp['site names'][j]]['dayNum'].shape, dtype='int16' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['dayNum']; #write that data
                    dataHolder = siteGroup.create_dataset('dayNumF', Mag_dataTemp[Mag_dataTemp['site names'][j]]['dayNumF'].shape, dtype='int16' ,compression="gzip"); #create dataset
                    dataHolder[...] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['dayNumF']; #write that data
                    #Write in attributes
                    siteGroup.attrs['lat'] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['lat']; #record the attribute
                    siteGroup.attrs['long'] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['long']-360; #record the attribute [transform to -180 to 180 longitude from 0 to 360 longitude]
                    siteGroup.attrs['elv'] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['elv']; #record the attribute
                    siteGroup.attrs['dataRate'] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['dataRate']; #record the attribute
                    siteGroup.attrs['dataRateF'] = Mag_dataTemp[Mag_dataTemp['site names'][j]]['dataRateF']; #record the attribute
                #END FOR j
                #Write in global attributes
                filez.attrs['version'] = version; #record the attribute [saves to the main file, not per data type]
            #END WITH
            
            os.remove(dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan\\temp.sec'); #delete the temp file
        #END IF
        
        #DELETE RAW DATA IF IT'S ON
        if( FLG_deleteRaw == 1 ):
            print('WARNING: MAGNETOMETER NRCAN DELETE RAW DATA IS ON. DELETING CURRENT RAW DATA [PROCESSING IS DONE]. ON '+str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1]).zfill(2)+'-'+str(dateRange_full[i,2]).zfill(2));
            siteData = glob.glob(dataPath + str(dateRange_dayNum_full[i,0]) + '\\NRCan\\*'+str(dateRange_full[i,0])+str(dateRange_full[i,1]).zfill(2)+str(dateRange_full[i,2]).zfill(2)+'vsec.sec.gz')
            for j in range(0,len(siteData)):
                os.remove(siteData[j]); #delete the temp file
            #END FOR j
        #END IF
        
        #TIME TO READ IN THAT DATA
        with h5py.File(dataPath_currFile, 'r') as filez:
            if( i == 0 ): #init the vars
                siteNames = filez['/siteNames'][:].astype('S3'); #import site names
                siteNames = np.char.decode(siteNames).astype(object); #convert back to unconstrained strings
                Mag_data = {'site names':siteNames}; #use a dictionary to keep all the data
                for j in range(0,len(siteNames)):
                    Mag_data[siteNames[j]] = {}; #prep a dict for the site
                    Mag_data[siteNames[j]]['mag'] = filez['/'+siteNames[j]+'/mag'][:].astype('float32'); #import data
                    Mag_data[siteNames[j]]['magF'] = filez['/'+siteNames[j]+'/magF'][:].astype('float32'); #import data
                    Mag_data[siteNames[j]]['sec'] = filez['/'+siteNames[j]+'/sec'][:].astype('int32'); #import data
                    Mag_data[siteNames[j]]['secF'] = filez['/'+siteNames[j]+'/secF'][:].astype('int32'); #import data
                    Mag_data[siteNames[j]]['year'] = filez['/'+siteNames[j]+'/year'][:].astype('int16'); #import data
                    Mag_data[siteNames[j]]['yearF'] = filez['/'+siteNames[j]+'/yearF'][:].astype('int16'); #import data
                    Mag_data[siteNames[j]]['dayNum'] = filez['/'+siteNames[j]+'/dayNum'][:].astype('int16'); #import data
                    Mag_data[siteNames[j]]['dayNumF'] = filez['/'+siteNames[j]+'/dayNumF'][:].astype('int16'); #import data
                    #read the attributes
                    Mag_data[siteNames[j]]['lat'] = filez[siteNames[j]].attrs['lat']; #read the attribute
                    Mag_data[siteNames[j]]['long'] = filez[siteNames[j]].attrs['long']; #read the attribute
                    Mag_data[siteNames[j]]['elv'] = filez[siteNames[j]].attrs['elv']; #read the attribute
                    Mag_data[siteNames[j]]['dataRate'] = filez[siteNames[j]].attrs['dataRate']; #read the attribute
                    Mag_data[siteNames[j]]['dataRateF'] = filez[siteNames[j]].attrs['dataRateF']; #read the attribute
                #END FOR j
            else: #add to the vars
                siteNamesTemp = filez['/siteNames'][:].astype('S3'); #import site names
                siteNamesTemp = np.char.decode(siteNamesTemp).astype(object); #convert back to unconstrained strings
                if( len(siteNamesTemp) > len(siteNames) ):
                    Mag_data['site names'] = siteNames; #update site names with the larger amount of site names
                #END IF
                for j in range(0,len(siteNames)):
                    if( (siteNames[j] in Mag_data) == False ):
                        Mag_data[siteNames[j]] = {}; #prep a dict for the site b/c is new
                    #END IF
                    Mag_data[siteNames[j]]['mag'] = np.concatenate( (Mag_data[siteNames[j]]['mag'], filez['/'+siteNames[j]+'/mag'][:].astype('float32')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['magF'] = np.concatenate( (Mag_data[siteNames[j]]['magF'], filez['/'+siteNames[j]+'/magF'][:].astype('float32')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['sec'] = np.concatenate( (Mag_data[siteNames[j]]['sec'], filez['/'+siteNames[j]+'/sec'][:].astype('int32')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['secF'] = np.concatenate( (Mag_data[siteNames[j]]['secF'], filez['/'+siteNames[j]+'/secF'][:].astype('int32')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['year'] = np.concatenate( (Mag_data[siteNames[j]]['year'], filez['/'+siteNames[j]+'/year'][:].astype('int16')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['yearF'] = np.concatenate( (Mag_data[siteNames[j]]['yearF'], filez['/'+siteNames[j]+'/yearF'][:].astype('int16')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['dayNum'] = np.concatenate( (Mag_data[siteNames[j]]['dayNum'], filez['/'+siteNames[j]+'/dayNum'][:].astype('int16')) , axis = 0); #import data
                    Mag_data[siteNames[j]]['dayNumF'] = np.concatenate( (Mag_data[siteNames[j]]['dayNumF'], filez['/'+siteNames[j]+'/dayNumF'][:].astype('int16')) , axis = 0); #import data
                    #read the attributes [they shouldn't change]
                    # Mag_data[siteNames[j]]['lat'] = filez[siteNames[j]].attrs['lat']; #read the attribute
                    # Mag_data[siteNames[j]]['long'] = filez[siteNames[j]].attrs['long']; #read the attribute
                    # Mag_data[siteNames[j]]['elv'] = filez[siteNames[j]].attrs['elv']; #read the attribute
                    # Mag_data[siteNames[j]]['dataRate'] = filez[siteNames[j]].attrs['dataRate']; #read the attribute
                    # Mag_data[siteNames[j]]['dataRateF'] = filez[siteNames[j]].attrs['dataRateF']; #read the attribute
                #END FOR j
            #END IF
        #END WITH
    #END FOR i
    
    return Mag_data #return Mag data dict with everything in it
