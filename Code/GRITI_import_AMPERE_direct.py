#Function to import AMPERE data
#RD on 3/5/2021

#To properly use, place pre-proccessed sources in their respective year folders in the data folder you've specified

#Info on stuff you can send:
# AMPERE_hemi - either 'north' or 'south'
    
#FLG_dataMix = 0; #0 prevents data sources to mixing to fill a time span, 1 allows data sources to mix to fill a time span
#FLG_dataPreference = 0; #preffered data type (by order as appended below at folder_fileNameFormat)

import numpy as np
import os
from urllib.request import urlretrieve
from Code.subfun_downloadProgress import downloadProgress
import h5py
import netCDF4 as netCDF4 #h5py won't vibe with AMPERE netCDF4 (even w/ the translation layer package)
from copy import deepcopy
import datetime
import aacgmv2 #install with: pip install aacgmv2 [need buildtools what a pain]
from scipy import interpolate
from Code.subfun_h5_reader import subfun_h5_reader
from Code.subfun_h5_writer import subfun_h5_writer
from Code.subfun_dictDelver import dictDelver_saver
from Code.subfun_strfind import strfind
from Code.subfun_daysInAYear import subfun_daysInAYear
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.subfun_textNice import textNice
# import time
import joblib
import subprocess
from Code.GRITI_import_AMPERE_direct_subfun_adelphi import adelphi_modeler
from Code.subfun_rbf_interp import rbf_interp, distance_calc_identical
from scipy.linalg import lu_factor, lu_solve

##-----Testing variables-----
##Date range goes Month-Day-Year
##dateRange = np.array([[2013,5,8],[2013,5,10]],dtype="int16"); #dates are in int16 because they can be
#dateRange = np.array([[2013,5,6],[2013,5,8]],dtype="int16"); #dates are in int16 because they can be
##dateRange = np.array([[2014,12,31],[2015,1,1]],dtype="int16"); #for debug, check year success
##dates better go earlier -> later
##print("{}".format(dateRange))
#folder = [os.getcwd()]; #current working directory, leave it as this call usually
#folder.append('E:\Big Data'); #place to save data files to
##folder var structure: 0 = running folder, 1 = data folder
#from subfun_date_to_dayNum import subfun_date_to_dayNum
#from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units
#FLG_dataMix = 0; #0 prevents data sources to mixing to fill a time span, 1 allows data sources to mix to fill a time span
#FLG_dataPreference = 0; #preffered data type (by order as appended below at folder_fileNameFormat)


def GRITI_import_AMPERE_direct(dates,settings_paths,AMPERE_hemi,settings_config,
                               AMPERE_coordType='geo',AMPERE_desired_latLongSteps=[1,15],FLG_paddedDays=True,FLG_dataMix=0,FLG_dataPreference=0,FLG_float64=0):
    version_alg = 2.4; #algorithm version (temporarially reduced)
    #1.0 11/15/2021 - initial algorithm
    #1.1 5/8/2022 - improved MLT handling
    #1.2 5/8/2022 - improved MLT handling more [last version before JH added, currently closed source so version based on if JH code there]
    #2.0 7/25/2022 - added JH-calc'n code named adelphi, added support for AMPERE NEXT [data still has teething issues]
    #2.1 10/14/2022 - let 1D adelphi parameters & basic data info stay with all resolutions
    #2.2 10/24/2022 - use custom spherical RBF interpolator instead of Cartesian RBF interpolator
    #2.3 12/20/2022 - add back in MLT being saved, support pos_geo's 3 dimensions (XYZ) properly, support less-than-720 times in the data, adelphi 1D params have time
    #2.4 6/27/2023 - fix interpolated offset (my b)
    
    parallel_numCores = settings_config['parallel num cores']; #num cores
    parallel_numThreads = settings_config['parallel num threads']; #num threads
    parallel_threadsPerProcess = settings_config['parallel threads per process']; #threads per process (ensures calcs are packed)
    
    #----- Unpack -----
    if( FLG_paddedDays == True ):
        dateRange_full = dates['date range full padded']; #use full padded range for analysis
        dateRange_dayNum_full = dates['date range full padded dayNum']; #use full padded range for analysis
    else:
        dateRange_full = dates['date range full']; #unpack it
        dateRange_dayNum_full = dates['date range full dayNum']; #unpack it
    #END IF    
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum']; #unpack it
    
    if( FLG_float64 == 0 ):
        dataAccuracy = np.float32; #set the data accuracy to 32 bit floats (1/2 the memory)
    else:
        dataAccuracy = np.float64; #set the data accuracy to 64 bit floats
    #END IF
    
    print("AMPERE Direct - Date range requested (yr/day num format): {}/{} to {}/{}.\n".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )

    #==============Constants Needed==============
    folder_AMPERE = 'AMPERE'; #name for the AMPERE folder
    folder_fileNameFormat = []; #prep a list
    # #YR for year
    # #DN for day num
    # #MO for month
    # #DY for day
    #.append more if more arise
    folder_fileNameFormat.append('AMPERED_#YR_#DN_#HEMI.h5'); #0 - 1st in the list is the compiled version, after that are possible names that need compiling
    folder_fileNameFormat.append('ampere_#YR#MO#DY_#HEMI_grd.nc'); #1 - expected file name style
    folder_fileNameFormat.append('ampere_#YR#MO#DY_#HEMI_grd.nc.7z'); #2 - expected file name style
    folder_fileNameFormat.append('ampere_#YR#MO#DY_#HEMI_grd.nc.gz'); #3 - expected file name style
    folder_fileNameFormat.append('ampere.#YR#MO#DY.k060_m08.#HEMI.grd.nc'); #4 - expected file name style
    folder_fileNameFormat.append('ampere.#YR#MO#DY.k060_m08.#HEMI.grd.nc.7z'); #5 - expected file name style
    folder_fileNameFormat.append('ampere.#YR#MO#DY.k060_m08.#HEMI.grd.nc.gz'); #6 - expected file name style
    folder_fileNameFormat.append('#YR#MO#DY.0000.86400.120.#HEMI.grd.ncdf'); #7 - expected file name style [old version of AMPERE data files]
    # support not quite yet folder_fileNameFormat.append('#YR#MO#DY.0000.86400.120.south.grd.ncdf'); #2 - expected file name style
    # folder_fileNameFormat.append('adelphi 2D Output #YR#MO#DY North Geographic'); #2 - expected file name style 2nd version
    #folder_fileNameFormat.append('otherFormat_#DY_#MO.h5'); #1 - tack on another (example)
    # AMPERE_dataRate = 120; #declare data rate, it's 2 min and it's in the file name (not anymore)
    
    #these are for writing a compact HDF5 file, these are the names that will be used for the dict names
    # headerNames = ['lat','JR','dJr','dBnorth1','dBeast1','dBnorth2','dBeast2','ddBnorth1','ddBeast1','ddBnorth2','ddBeast2','dBr','ddBr','dBprt']; #note the nice names [names not here will be tacked on as they are in the header]
    
                                 
    #==============Prep to look for the data==============       
    if( os.path.isdir( os.path.join(settings_paths['data'], folder_AMPERE) ) == 0 ): #check if AMPERE folder exists
        #if not, make it
        os.makedirs( os.path.join(settings_paths['data'], folder_AMPERE) );
        print("NOTA BENE: Importing AMPERE Func - Created AMPERE directory: {}\n".format( os.path.join(settings_paths['data'], folder_AMPERE) ) );
    #END IF
    
    dateRange_uniqueYears = np.unique(dateRange_dayNum_full[:,0]); #get the unique years involved
    AMPERE_dataAmnt = len(dateRange_dayNum_full[:,1]); #get number of days needed to investigate
    AMPERE_dataAvail_perSource = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for finding where data is available
    AMPERE_dataAvail_toUse = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for recording which piece of data to use on the day
    #0 = no data available on required days, will quit
    #1 = note data is there, ready to use
    
    for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
        if( os.path.isdir( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_uniqueYears[i])) ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_uniqueYears[i])) );
            print("NOTA BENE: Importing AMPERE Func - Created AMPERE subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i], os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_uniqueYears[i])) ));
        #END IF
    #END FOR
    
    #==============Look for data in expected naming formats==============
    for i in range(0,len(folder_fileNameFormat)): #loop through the different name formats    
        
        for j in range(0,AMPERE_dataAmnt): #loop through the different years needed
            
            AMPERE_fileName = folder_fileNameFormat[i].replace('#DN', str(dateRange_dayNum_full[j,1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName = AMPERE_fileName.replace('#YR', str(dateRange_dayNum_full[j,0]) ); #replace any #YR with the current year
            AMPERE_fileName = AMPERE_fileName.replace('#MO', str(dateRange_full[j,1]).zfill(2) ); #replace any #MO with the current month (padded w/ 0's so always 2 long)
            AMPERE_fileName = AMPERE_fileName.replace('#DY', str(dateRange_full[j,2]).zfill(2) ); #replace any #DY with current day (padded w/ 0's so always 2 long)
            AMPERE_fileName = AMPERE_fileName.replace('#HEMI', AMPERE_hemi ); #replace any #HEMI with current hemisphere
            
            AMPERE_dataAvail_perSource[j,i] = os.path.isfile( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[j,0]), AMPERE_fileName) ); #record if data is there or not
        #END FOR j
    #END FOR i
        
    FLG_dataAvail_entireSpan = 0; #flag that needs to be set to 1 with this data finding alg
    if( FLG_dataMix == 1): #this allows for the data sources to mix to cover more availability
        AMPERE_dataAvail_entireSpan = np.any( AMPERE_dataAvail_perSource , axis=1 ); #see if each day has some availability
        
        if( np.all(AMPERE_dataAvail_entireSpan) ): #if all are 1 this is true
            #only do the work if we're good to go
            AMPERE_dataAvail_toUse = np.copy(AMPERE_dataAvail_perSource); #copy this over
            for i in range(0,AMPERE_dataAmnt): #check each entry for multiples, etc.
                
                if( np.sum(AMPERE_dataAvail_toUse[i,:]) > 1 ): #if the sum is greater than 1, choose one
                    
                    if( AMPERE_dataAvail_toUse[i,FLG_dataPreference] == 1 ): #easy, set it
                        AMPERE_dataAvail_toUse[i,:] = False; #0 it out
                        AMPERE_dataAvail_toUse[i,FLG_dataPreference] = 1; #set it to 1
                    else:
                        AMPERE_dataAvail_toUse[i,np.where(AMPERE_dataAvail_perSource[i,:])[0][0]] = 1; #use the closest to 0
                        # print("\n==============ERROR in GRITI_import_AMPERE_direct==============");
                        # print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                        # #return("No"); #return something that will def crash things
                        # import sys #yolo import
                        # sys.crash(); #more def will crash
                        # #don't feel like writing a disp, but def I didn't get around to writing alternate stuff
                    #END IF
                else:
                    AMPERE_dataAvail_toUse[i,np.where(AMPERE_dataAvail_perSource[i,:])[0][0]] = 1; #use the closest to 0
                #END IF
            #END FOR i    
            
            FLG_dataAvail_entireSpan = 1; #set the flag to good
        #END IF    
    
    else: #otherwise one data source has to cover the entire span
        
        AMPERE_dataAvail_entireSpan = np.all( AMPERE_dataAvail_perSource , axis=0 ); #see if each source has complete availability
        
        if( np.any(AMPERE_dataAvail_entireSpan) ): #if any source has data for all dates, this will be true
            #don't do the work if it's not true
        
            if( np.sum(AMPERE_dataAvail_entireSpan) > 1): #make sure to choose one data source
                            
                if( AMPERE_dataAvail_entireSpan[FLG_dataPreference] == 1 ): #easy, data type preferred is there
                    # AMPERE_dataAvail_entireSpan[:] = False; #0 it out
                    # AMPERE_dataAvail_entireSpan[FLG_dataPreference]  = 1; #put a 1 where the preffered data is there
                    AMPERE_dataAvail_toUse[:,FLG_dataPreference] = 1; #use the data preference
                else: 
                    # #don't feel like writing a disp, but def I didn't get around to writing alternate stuff
                    # print("\n==============ERROR in GRITI_import_AMPERE_direct==============");
                    # print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                    # #return("No"); #return something that will def crash things
                    # import sys #yolo import
                    # sys.exit(); #more def will crash
                    AMPERE_dataAvail_toUse[:,np.where(AMPERE_dataAvail_entireSpan)[0][0]] = 1; #use the closest to 0
                #END IF
            #END IF
                
            # AMPERE_dataAvail_toUse[:,AMPERE_dataAvail_entireSpan] = 1; #set 
            
            FLG_dataAvail_entireSpan = 1; #set the flag to good
        else:
            print('WARNING in GRITI_import_AMPERE_direct: FLG_dataMix = 0 and there wasn\'t one type of data to cover the entire span. This will likely cause all new raw data to be downloaded that will be the same type - it\'ll take a bit. And if that wasn\'t your goal (custom raw data or something), go fix that tiger');
        #END IF
    #END IF
        
    if( FLG_dataAvail_entireSpan == 0 ): #make sure there's data
        # print("\n==============ERROR in GRITI_import_AMPERE_direct==============");
        # print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
        # print("Printing file name formats supported:");
        # print("{}".format(folder_fileNameFormat)); #print for error
        # print("Printing available data matrix (made of dates and file name formats):");
        # print("{}".format(AMPERE_dataAvail_perSource)); #print for error - lets user know available days
        # print("Will exit via returning no");
        # #return("No"); #return something that will def crash things
        # import sys #yolo import
        # sys.exit(); #more def will crash
        
        # curl "https://ampere.jhuapl.edu/services/data-grd.php?logon=rld5204&start=2013-05-06T00:00&extent=86400&pole=north" > north.nc
        
        AMPERE_dataAvail_toUse = np.copy(AMPERE_dataAvail_perSource); #copy this over
        
        AMPERE_dataAvail_missing = np.where(AMPERE_dataAvail_entireSpan == False)[0]; #there's a way to download via easy links now!
        for i in range(0,AMPERE_dataAvail_missing.size):
            AMPERE_fileName = folder_fileNameFormat[1].replace('#DN', str(dateRange_dayNum_full[AMPERE_dataAvail_missing[i],1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName = AMPERE_fileName.replace('#YR', str(dateRange_dayNum_full[AMPERE_dataAvail_missing[i],0]) ); #replace any #YR with the current year
            AMPERE_fileName = AMPERE_fileName.replace('#MO', str(dateRange_full[AMPERE_dataAvail_missing[i],1]).zfill(2) ); #replace any #MO with the current month (padded w/ 0's so always 2 long)
            AMPERE_fileName = AMPERE_fileName.replace('#DY', str(dateRange_full[AMPERE_dataAvail_missing[i],2]).zfill(2) ); #replace any #DY with current day (padded w/ 0's so always 2 long)
            AMPERE_fileName = AMPERE_fileName.replace('#HEMI', AMPERE_hemi ); #replace any #HEMI with current hemisphere
            
            AMPERE_dateMissing = str(dateRange_full[AMPERE_dataAvail_missing[i],0])+'-'+ \
                str(dateRange_full[AMPERE_dataAvail_missing[i],1]).zfill(2)+'-'+ \
                str(dateRange_full[AMPERE_dataAvail_missing[i],2]).zfill(2);
                
            print('WARNING in GRITI_import_AMPERE_direct: Day '+AMPERE_dateMissing+' was missing any processed or raw data. Downloading raw now. (note that file doesn\'t declare how big it is so the download % is actually just the MB downloaded so far (files usually end around 179.74 MB [~171 MiB in Windows])');
            tryLim = 3;
            tryNum = 0;
            while(tryNum < tryLim):
                try:
                    urlretrieve('https://ampere.jhuapl.edu/services/data-grd.php?logon='+settings_config['login AMPERE']['user']+ \
                                '&start='+AMPERE_dateMissing+ \
                                'T00:00&extent=86400&pole='+AMPERE_hemi.lower(), \
                                os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[AMPERE_dataAvail_missing[i],0]), AMPERE_fileName), \
                                reporthook=downloadProgress); #download that file
                    tryNum = tryLim+1; #free of the loop once download succeeds
                except:
                    tryNum += 1; #increment
                    print('WARNING in GRITI_import_AMPERE_direct: Download failed. Trying again, try #'+str(tryNum+1)+'/'+str(tryLim));
                #END TRY
            #END WHILE
            if( tryNum > tryLim ):
                AMPERE_dataAvail_entireSpan[AMPERE_dataAvail_missing[i]] = True; #after download set to true, we do ahve the data
                AMPERE_dataAvail_toUse[AMPERE_dataAvail_missing[i],1] = 1; #set since we do have the data
                print(''); #extra line needed at end
            else:
                print('ERROR in GRITI_import_AMPERE_direct: Downloading '+AMPERE_fileName+' failed '+str(tryLim)+' times (the coded limit). Crashing b/c data won\'t be there for the code. Try downloading it yourself from "https://ampere.jhuapl.edu/download/?page=dataTab" - download from "Daily NetCDF Files" tab and the "Filtered" data for your day(s) & hemisphere. (the files can stay zipped)');
                os.remove(os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[AMPERE_dataAvail_missing[i],0]), AMPERE_fileName)); #cleanup since partial file will be there
                import sys #yolo import
                sys.exit(); #more def will crash
            #END IF
            
            if( os.path.getsize(os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[AMPERE_dataAvail_missing[i],0]), AMPERE_fileName)) < 100000 ):
                # a really small file shows up if the data wasn't there apparently
                print('ERROR in GRITI_import_AMPERE_direct: '+AMPERE_fileName+' is not available, DL\'d a smol placeholder file. There is no data for this day, change date range or code a way to bypass a day with NaNs or something. Check out "https://ampere.jhuapl.edu/download/?page=dataTab" at "Daily NetCDF Files" tab and the "Filtered" data for your day(s) & hemisphere - but it ain\'t there. ¯\_(ツ)_/¯');
                os.remove(os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[AMPERE_dataAvail_missing[i],0]), AMPERE_fileName)); #cleanup since partial file will be there
                import sys #yolo import
                sys.exit(); #more def will crash
            #END IF
        #END FOR i
        FLG_dataAvail_entireSpan = 1; #set the flag to good
    #END IF
    
    #Code to read year from file name, decided unneeded
    #AMPERE_date = np.zeros(size(fileName_AMPERE,1),2,'int16'); #good till yr 32,767
    #for( i = 1:size(fileName_AMPERE,1) )
    #    jk = strfind(fileName_AMPERE(i,:),'_'); #get all the _'s
    #    tempYr = str2double(fileName_AMPERE(i, (jk(end)+1):((jk(end)+1)+3) )); #yr, get the year
    #    tempMon = str2double(fileName_AMPERE(i, (jk(end)+1+4):((jk(end)+1)+5) )); #mon, get the month
    #    tempDay = str2double(fileName_AMPERE(i, (jk(end)+1+6):end )); #day, get the day
    #    AMPERE_date(i,1) = tempYr; #yr, save year
    #    AMPERE_date(i,2) = sFUN_dateToDayNum([tempYr,tempMon,tempDay]); #dayNum, save day number (what we use for ez-er math)
    #end
    
    #==============Read in the data==============
    AMPERE_timeRuns = 0; #counts the runs (for calcin/recalcin data)
    AMPERE_data = {}; #prep data dict to hold the data
    FLG_firstFileMade = True; #flag for first file made
    for i in range(0,AMPERE_dataAmnt ):
        FLG_recalc = False; #set recalc flag to false every time
        FLG_append = False; #set append flag to false every time
        AMPERE_sourceIndex = np.where(AMPERE_dataAvail_toUse[i,:] == 1)[0]; #get the location of the index (corresponds to which source)
        if( np.any(AMPERE_sourceIndex == FLG_dataPreference) ): #if any of the ones that have data are the chosen source, then use it
            AMPERE_sourceIndex = np.int64(FLG_dataPreference); #just set it n forget it
        else: #otherwise just use the next one in line (don't have a hierarchy yet)
            AMPERE_sourceIndex = AMPERE_sourceIndex[0]; #just get the first one in line
        #END IF
        AMPERE_fileName = folder_fileNameFormat[ AMPERE_sourceIndex ].replace('#DN', str(dateRange_dayNum_full[i,1]).zfill(3) ); #replace any #DN with the current day number
        AMPERE_fileName = AMPERE_fileName.replace('#YR', str(dateRange_dayNum_full[i,0]) ); #replace any #YR with the current year
        AMPERE_fileName = AMPERE_fileName.replace('#MO', str(dateRange_full[i,1]).zfill(2) ); #replace any #MO with the current month
        AMPERE_fileName = AMPERE_fileName.replace('#DY', str(dateRange_full[i,2]).zfill(2) ); #replace any #DY with current day
        AMPERE_fileName = AMPERE_fileName.replace('#HEMI', AMPERE_hemi ); #replace any #HEMI with current hemisphere
        #prep the name needed
        
        if( np.all(AMPERE_desired_latLongSteps == None) ):
            AMPERE_desiredResStr = '1&15'; #the default data type AMPERE comes with, hardcoded b/c I don't wanna bootstrap my self outa this
            # AMPERE_desired_latDelta = 1; #so it's also accessible
            # AMPERE_desired_longDelta = 15; #so it's also accessible
        else:
            AMPERE_desired_latDelta = AMPERE_desired_latLongSteps[0]; #get lat steps
            AMPERE_desired_longDelta = AMPERE_desired_latLongSteps[1]; #get long steps
            AMPERE_desiredResStr = textNice(AMPERE_desired_latDelta)+'&'+textNice(AMPERE_desired_longDelta); #desired res string for naming
        #END IF
        
        if( AMPERE_sourceIndex == 0 ): #this means it is the first type that has been refactored into an HDF5 file
            #load saved AMPERE to speed up
            keyz = list(AMPERE_data.keys()); #get the current keys in AMPERE_data
            with h5py.File( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName), 'r') as AMPERE_file:
                try:
                    if( AMPERE_file.attrs['version'] >= version_alg ):
                        #--- Read the data types in ---
                        keyzNew = list(AMPERE_file.keys()); #get the saved keys
                        if( AMPERE_coordType+AMPERE_desiredResStr in AMPERE_file.keys() ):
                            AMPERE_data = subfun_h5_reader(AMPERE_file, h5_path_prefix=AMPERE_coordType+AMPERE_desiredResStr, diction=AMPERE_data); #call h5 reader to convert to dict
                            
                            #--- Read the top level attributes in ---
                            for key, item in AMPERE_file.attrs.items():
                                if( key in AMPERE_data.keys() ):
                                    if( np.isclose(AMPERE_data[key], item) == False ):
                                        #only worry is if the attribute isn't consistent
                                        print('-----Warning in subfun_h5_reader-----');
                                        print('Attribute '+key+' isn\'t the same as the previously recorded value from another file of '+ \
                                            str(AMPERE_data[key])+' and this ('+str(AMPERE_file.filename)+') file\'s value of '+str(item)+ \
                                            '.\n NaN\'ing it and try to sort that out.');
                                        AMPERE_data[key] = np.nan; #nan that attribute, figure it out later
                                    #END IF
                                else:
                                    AMPERE_data[key] = item; #get that attribute out
                                #END IF
                            #END FOR key, item
                        else:
                            FLG_append = True; #set flag to true, data append for different gridding is needed
                            AMPERE_fileName_append = AMPERE_fileName; #keep the file name around
                        #END IF
                    else:
                        FLG_recalc = True; #set flag to true, recalc is needed
                        oldVersion = AMPERE_file.attrs['version'];
                    #END IF
                except KeyError:
                    FLG_recalc = True; #file is broken, recalc it
                    oldVersion = 0; #set to 0 b/c it's broken but need the thing to work
                    print('---WARNING in GRITI_import_AMPERE_direct---');
                    print( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName)+ \
                        ' file is broken, KeyError occured. File version set to 0 to make code go. Recalcing over it.');
                #END TRY
            #END WITH
        #END IF
        if( (FLG_recalc == True) ):
            print('---WARNING in GRITI_import_AMPERE_direct---');
            print( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName)+ \
                ' file version is '+str(oldVersion)+' which is less than the current file version of '+str(version_alg)+'!'+\
                '\nThat\'s bad! Renaming file to add a \'_oldV'+str(oldVersion).replace('.','p')+'\' at the end, then gonna calculate a new version of the file!');
            os.rename( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName), \
                os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName+'_oldV'+str(oldVersion).replace('.','p')) ); #rename           
        
            #rerun the sourceIndexer b/c needed to import the data but use _perSource b/c its req'd
            AMPERE_dataAvail_toUse_tmp = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for recording which piece of data to use on the day
            if( FLG_dataMix == 0 ):
                AMPERE_dataAvail_entireSpan_tmp = np.all( AMPERE_dataAvail_perSource[:,1:] , axis=0 ); #see if each source has complete availability
                if( np.any(AMPERE_dataAvail_entireSpan_tmp) ):
                    AMPERE_dataAvail_entireSpan_where = np.where(AMPERE_dataAvail_entireSpan_tmp)[0]+1;
                    if( np.any(AMPERE_dataAvail_entireSpan_where == FLG_dataPreference) ):
                        AMPERE_dataAvail_toUse_tmp[:,FLG_dataPreference] = 1; #use the data preference
                    else:
                        AMPERE_dataAvail_toUse_tmp[:,AMPERE_dataAvail_entireSpan_where[0]] = 1; #use the closest to 0
                    #END IF
                else:
                    print("\n==============ERROR in GRITI_import_AMPERE_direct==============");
                    print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
                    print("Printing file name formats supported:");
                    print("{}".format(folder_fileNameFormat)); #print for error
                    print("Printing available data matrix (made of dates and file name formats):");
                    print("{}".format(AMPERE_dataAvail_perSource)); #print for error - lets user know available days
                    print("Will exit via returning no");
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.crash(); #more def will crash
                #END IF
            else:
                if( np.sum(AMPERE_dataAvail_perSource[i,1:]) > 0 ): #if the sum is greater than 0, choose one
                    if( (AMPERE_dataAvail_perSource[i,FLG_dataPreference] == 1) & (FLG_dataPreference != 0) ): #easy, set it
                        AMPERE_dataAvail_toUse_tmp[i,FLG_dataPreference] = 1; #set it to 1
                    else:
                        AMPERE_dataAvail_toUse_tmp[i,np.where(AMPERE_dataAvail_perSource[i,1:])[0][0]+1] = 1; #use the closest to 0
                    #END IF
                else:
                    print("\n==============ERROR in GRITI_import_AMPERE_direct==============");
                    print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
                    print("Printing file name formats supported:");
                    print("{}".format(folder_fileNameFormat)); #print for error
                    print("Printing available data matrix (made of dates and file name formats):");
                    print("{}".format(AMPERE_dataAvail_perSource)); #print for error - lets user know available days
                    print("Will exit via returning no");
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.crash(); #more def will crash
                #END IF
            #END IF
            AMPERE_sourceIndex = np.where(AMPERE_dataAvail_toUse_tmp[i,1:] == 1)[0]+1; #get the location of the index (corresponds to which source)
            if( np.any(AMPERE_sourceIndex == FLG_dataPreference) ): #if any of the ones that have data are the chosen source, then use it
                AMPERE_sourceIndex = np.int64(FLG_dataPreference); #just set it n forget it
            else: #otherwise just use the next one in line (don't have a hierarchy yet)
                AMPERE_sourceIndex = AMPERE_sourceIndex[0]; #just get the first one in line
            #END IF
            AMPERE_fileName = folder_fileNameFormat[ AMPERE_sourceIndex ].replace('#DN', str(dateRange_dayNum_full[i,1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName = AMPERE_fileName.replace('#YR', str(dateRange_dayNum_full[i,0]) ); #replace any #YR with the current year
            AMPERE_fileName = AMPERE_fileName.replace('#MO', str(dateRange_full[i,1]).zfill(2) ); #replace any #MO with the current month
            AMPERE_fileName = AMPERE_fileName.replace('#DY', str(dateRange_full[i,2]).zfill(2) ); #replace any #DY with current day
            AMPERE_fileName = AMPERE_fileName.replace('#HEMI', AMPERE_hemi ); #replace any #HEMI with current hemisphere
        elif( FLG_append == True ):
            print('---WARNING in GRITI_import_AMPERE_direct---');
            print( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName)+ \
                ' is lacking requested lat&long ('+AMPERE_desiredResStr+') so calculating it and appending to current file.');
        #END IF
        
        if( AMPERE_timeRuns == 0 ):
            if( FLG_append == False ):
                AMPERE_timeEst = 9.5*(((50/AMPERE_desired_latLongSteps[0])*(360/AMPERE_desired_latLongSteps[1]))/(1200))**1.7; #estimates the times for newly calc'd data, scales based on resolution roughly
            else:
                AMPERE_timeEst = 8.5*(((50/AMPERE_desired_latLongSteps[0])*(360/AMPERE_desired_latLongSteps[1]))/(1200))**1.7; #estimates the times for appended data
            #END IF
        #END IF
        if( (AMPERE_sourceIndex >= 1) | (FLG_append == True) ): #this means it is the first type that has been refactored into an HDF5 file 
            #else read AMPERE files to get the data needed
            if( FLG_append == False ):
                print('Some or none of "cached" pre-read AMPERE files found; beginning to read AMPERE data and save as HDF5 for future use - est. '+textNice(np.round(AMPERE_timeEst,2))+' min per ('+textNice(np.round(AMPERE_timeEst*(AMPERE_dataAmnt-np.sum(np.int64(AMPERE_dataAvail_perSource[:,0]))-AMPERE_timeRuns),2))+' min total) on fast comp');
            else:
                print('Appending to "cached" pre-read AMPERE file - est. '+textNice(np.round(AMPERE_timeEst,2))+' min per on fast comp');
            #ENDIF
            import time
            # from subfun_strstr import strstr
            tic = time.time(); #start timing
            
            FLG_directRead = True; #idea here is appending requires 1&15 gridding
            if( FLG_append == True ):
                #prime the looking
                AMPERE_save = {}; #prep
                
                with h5py.File( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName_append), 'r') as AMPERE_file:
                    #--- Read the data types in ---
                    keyzNew = list(AMPERE_file.keys()); #get the saved keys
                    if( ('mag orig' in keyzNew) ):
                        keyzNew = ['mag orig']; #only read the needed keys
                        for j in range(0,len(keyzNew)):
                            sub_keyz = list(AMPERE_file.get(keyzNew[j]).keys()); #get the group keys [assume group here is safe]
                            AMPERE_save[keyzNew[j]] = {}; #prep the dict to hold things
                            for k in range(0,len(sub_keyz)):
                                if( not type(AMPERE_file[keyzNew[j]].get(sub_keyz[k])) is h5py.Group ):
                                    AMPERE_save[keyzNew[j]][sub_keyz[k]] = AMPERE_file[keyzNew[j]].get(sub_keyz[k])[()]; #get that dataset out
                                else:
                                    AMPERE_save[keyzNew[j]][sub_keyz[k]] = {}; #prep a sub-group
                                    subsub_keyz = list(AMPERE_file[keyzNew[j]][sub_keyz[k]].keys()); #get the group keys
                                    for kj in range(0,len(subsub_keyz)):
                                        AMPERE_save[keyzNew[j]][sub_keyz[k]][subsub_keyz[kj]] = AMPERE_file[keyzNew[j]][sub_keyz[k]].get(subsub_keyz[kj])[()]; #get that dataset out
                                    #END FOR kj
                                    subsub_keyz = list(AMPERE_file[keyzNew[j]][sub_keyz[k]].attrs.keys()); #get the attribute keys
                                    for kj in range(0,len(subsub_keyz)):
                                        AMPERE_save[keyzNew[j]][sub_keyz[k]][subsub_keyz[kj]] = AMPERE_file[keyzNew[j]][sub_keyz[k]].attrs[subsub_keyz[kj]]; #get that attribute out
                                    #END FOR kj
                                #END IF
                            #END FOR k
                            #--- Read the attributes in ---
                            sub_keyzCurr = list(AMPERE_save[keyzNew[j]]); #get the keys that are now in there
                            sub_keyz = list(AMPERE_file[keyzNew[j]].attrs.keys()); #get the attribute keys
                            sub_keyz = np.delete(np.asarray(sub_keyz),np.where(strfind(sub_keyz,'version'))[0]).tolist(); #remove version from the list, only needed here
                            for k in range(0,len(sub_keyz)):
                                if( strfind(sub_keyzCurr,sub_keyz[k],opt=1) > 0 ):
                                    if( np.isclose(AMPERE_file[keyzNew[j]].attrs[sub_keyz[k]], AMPERE_save[keyzNew[j]][sub_keyz[k]]) == False ):
                                        #only worry is if the attribute isn't consistent
                                        print('-----Warning in GRITI_import_AMPERE_direct-----');
                                        print('Attribute '+sub_keyz[k]+' isn\'t the same as the previously recorded value from another file of '+ \
                                            str(AMPERE_save[keyzNew[j]][sub_keyz[k]])+' and this file\'s value of '+str(AMPERE_file[keyzNew[j]].attrs[sub_keyz[k]])+ \
                                            '.\n NaN\'ing it and try to sort that out.');
                                        AMPERE_save[keyzNew[j]][sub_keyz[k]] = np.nan; #nan that attribute, figure it out later
                                    #END IF
                                else:
                                    AMPERE_save[keyzNew[j]][sub_keyz[k]] = AMPERE_file[keyzNew[j]].attrs[sub_keyz[k]]; #get that attribute out
                                #END IF
                            #END FOR k
                        #END FOR j
                        #--- Read the attributes in ---
                        keyzNew = list(AMPERE_file.attrs.keys()); #get the attribute keys
                        keyzNewCurr = list(AMPERE_save); #get the keys that are now in there
                        keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                        for j in range(0,len(keyzNew)):
                            if( strfind(keyzNewCurr,keyzNew[j],opt=1) > 0 ):
                                if( np.isclose(AMPERE_file.attrs[keyzNew[j]], AMPERE_save[keyzNew[j]]) == False ):
                                    #only worry is if the attribute isn't consistent
                                    print('-----Warning in GRITI_import_AMPERE_direct-----');
                                    print('Attribute '+keyzNew[j]+' isn\'t the same as the previously recorded value from another file of '+ \
                                        str(AMPERE_save[keyzNew[j]])+' and this file\'s value of '+str(AMPERE_file.attrs[keyzNew[j]])+ \
                                        '.\n NaN\'ing it and try to sort that out.');
                                    AMPERE_save[keyzNew[j]] = np.nan; #nan that attribute, figure it out later
                                #END IF
                            else:
                                AMPERE_save[keyzNew[j]] = AMPERE_file.attrs[keyzNew[j]]; #get that attribute out
                            #END IF
                        #END FOR j
                        FLG_directRead = False; #no need for a direct read
                        #assemble some vars needed from the saved data
                        AMPERE_temp_timeData = AMPERE_save['mag orig']['data info']['num time steps']; #number of time step
                        AMPERE_temp_numData = AMPERE_save['mag orig']['data info']['num data per time step']; #number of data per time step
                        AMPERE_temp_latNum = AMPERE_save['mag orig']['data info']['lat num']; #number of lat num
                        AMPERE_temp_longNum = AMPERE_save['mag orig']['data info']['long num']; #number of long num
                        AMPERE_temp_latDelta = AMPERE_save['mag orig']['data info']['lat delta']; #deg, delta of lat pts provided
                        AMPERE_temp_longDelta = AMPERE_save['mag orig']['data info']['long delta']; #deg, delta of long pts provided (assume 360 coverage is safe)
                    #END IF
                #END WITH
            #END IF
            if( FLG_directRead == True ):
                #load netcdf4 file from AMPERE
                AMPERE_temp = {}; #hold the raw info from the netcdf4 file
                # keyz = list(AMPERE_temp.keys()); #get the current keys in AMPERE_temp    
                if( (AMPERE_fileName[AMPERE_fileName.rfind('.'):] == '.gz') | \
                    (AMPERE_fileName[AMPERE_fileName.rfind('.'):] == '.7z') ):
                    #must unzip if so
                    subprocess.run('"'+settings_paths['7zip']+'"'+' e '+'"'+ \
                                   os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]),AMPERE_fileName)+'"' \
                                   ' -o"'+os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]))+'"'); #call 7z to extract
                    AMPERE_fileName_origName = AMPERE_fileName; #record orig name
                    AMPERE_fileName = AMPERE_fileName[:AMPERE_fileName.rfind('.')]; #remove compression ending
                else:
                    AMPERE_fileName_origName = AMPERE_fileName; #samesies
                #END IF
                
                with netCDF4.Dataset(os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]),AMPERE_fileName), 'r') as AMPERE_file:
                    #--- Read the data types in ---
                    keyzNew = list(AMPERE_file.variables.keys()); #get the saved keys               
                    for j in range(0,len(keyzNew)):
                        if( AMPERE_file[keyzNew[j]][:].data.ndim < 3 ):
                            if( np.issubdtype(AMPERE_file[keyzNew[j]][:].dtype, np.integer) ):
                                AMPERE_temp[keyzNew[j]] = AMPERE_file[keyzNew[j]][:].data.flatten(); #get that dataset out, keep integers integers
                            else:
                                AMPERE_temp[keyzNew[j]] = dataAccuracy(AMPERE_file[keyzNew[j]][:].data.flatten()); #get that dataset out, makesure float data accuracy is right
                        else:
                            if( np.issubdtype(AMPERE_file[keyzNew[j]][:].dtype, np.integer) ):
                                AMPERE_temp[keyzNew[j]] = AMPERE_file[keyzNew[j]][:].data.reshape(-1, AMPERE_file[keyzNew[j]][:].data.shape[-1]); #get that dataset out, keep integers integers
                            else:
                                AMPERE_temp[keyzNew[j]] = dataAccuracy(AMPERE_file[keyzNew[j]][:].data.reshape(-1, AMPERE_file[keyzNew[j]][:].data.shape[-1])); #get that dataset out
                            #END IF
                        #END IF
                    #END FOR j
                #END WITH
                
                if( AMPERE_fileName_origName[AMPERE_fileName_origName.rfind('.'):] == '.gz' ):
                    os.remove(os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]),AMPERE_fileName_origName)); #remove original compressed file, gonna do better b/c don't hit this often
                #END IF
                subprocess.run('"'+settings_paths['7zip']+'"'+' a '+ \
                               '-o"'+os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]))+'"'+' -mx9 '+ \
                               '"'+os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]),AMPERE_fileName)+'.7z" '+ \
                               '"'+os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]),AMPERE_fileName)+'"'); #call 7z to compress MAXIMUM POWER
                os.remove(os.path.join(settings_paths['data'],folder_AMPERE,str(dateRange_full[i,0]),AMPERE_fileName)); #remove uncompressed file b/c why keep it
                
                #Push V2 (AMPERE NEXT) format into V1 naming so that I don't have to change the later code (committed to no regartts here)
                if( 'colat' not in list(AMPERE_temp.keys()) ):
                    AMPERE_temp['nlat'] = AMPERE_temp['nLatGrid']; #simple alias
                    AMPERE_temp['nlon'] = AMPERE_temp['nLonGrid']; #simple alias
                    AMPERE_temp['colat'] = AMPERE_temp['cLat_deg'].ravel(); #simple alias + flatten
                    AMPERE_temp['mlt'] = AMPERE_temp['mlt_hr'].ravel(); #simple alias + flatten
                    AMPERE_temp['start_yr'] = AMPERE_temp['year']; #simple alias
                    AMPERE_temp['start_dyN'] = AMPERE_temp['doy']; #simple alias
                    AMPERE_tempDates = subfun_dayNum_to_date(np.vstack((AMPERE_temp['start_yr'],AMPERE_temp['start_dyN'])).T); #get it back to yr/month/day stuff for compatibility
                    AMPERE_temp['start_mo'] = AMPERE_tempDates[:,1]; #pull em out
                    AMPERE_temp['start_dy'] = AMPERE_tempDates[:,2]; #pull em out
                    AMPERE_temp['start_hr'] = np.int64(AMPERE_temp['time']);
                    AMPERE_temp['start_mt'] = np.int64(np.round((AMPERE_temp['time']-AMPERE_temp['start_hr'])*60,3));
                    AMPERE_temp['start_sc'] = np.int64(np.round(((AMPERE_temp['time']-AMPERE_temp['start_hr'])*60-AMPERE_temp['start_mt'])*60,3));                 
                    AMPERE_temp['Jr'] = AMPERE_temp['jPar']; #simple alias
                    AMPERE_temp['del_Jr'] = AMPERE_temp['del_jPar']; #simple alias
                    
                    AMPERE_dataRate = np.int64(np.median(np.diff(AMPERE_temp['time']*3600))); #sec, data rate
                    
                    #early data had an error where mlt_hr was the same as cLatDeg - not very useful
                    if( np.all(np.isclose(AMPERE_temp['mlt'],AMPERE_temp['colat'])) ): #catches error in the data
                        AMPERE_temp['mlt'] = np.tile(np.repeat(np.arange(0,24),AMPERE_temp['nlat'][0]),(AMPERE_temp['time'].size,1)).ravel(); #make it myself I guess
                    #END IF
                else:
                    #V1 naming to V2 naming for limited things (dBnorth1 -> db_Th_Th which I assume is better representative of what it is)
                    AMPERE_temp['db_Th_Th'] = AMPERE_temp['dBnorth1']; #simple alias
                    AMPERE_temp['db_Ph_Th'] = AMPERE_temp['dBeast1']; #simple alias
                    AMPERE_temp['db_Th_Ph'] = AMPERE_temp['dBnorth2']; #simple alias
                    AMPERE_temp['db_Ph_Ph'] = AMPERE_temp['dBeast2']; #simple alias
                    AMPERE_temp['db_R'] = AMPERE_temp['dBr']; #simple alias
                    AMPERE_temp['del_db_Th_Th'] = AMPERE_temp['ddBnorth1']; #simple alias
                    AMPERE_temp['del_db_Ph_Th'] = AMPERE_temp['ddBeast1']; #simple alias
                    AMPERE_temp['del_db_Th_Ph'] = AMPERE_temp['ddBnorth2']; #simple alias
                    AMPERE_temp['del_db_Ph_Ph'] = AMPERE_temp['ddBeast2']; #simple alias
                    AMPERE_temp['del_db_R'] = AMPERE_temp['ddBr']; #simple alias
                    AMPERE_temp['del_Jr'] = AMPERE_temp['dJr']; #simple alias
                    #fillers
                    AMPERE_temp['R'] = np.empty(AMPERE_temp['Jr'].size,dtype=AMPERE_temp['Jr'].dtype); #fill in for R, if V1 data and R are ever needed gotta write some stuff
                    AMPERE_temp['pos_geo'] = np.empty(AMPERE_temp['Jr'].size,dtype=AMPERE_temp['Jr'].dtype); #fill in for pos_geo, if V1 data and R are ever needed gotta write some stuff
                #END IF
                                
                #!! turns out the data served can be bad and will have 0's if an error occured with the data serving system!!
                #catch bad data and remove it (all 0's for every input)
                AMPERE_temp_badData = np.where(np.diff(AMPERE_temp['colat']) == 0)[0]; #get where bad data are
                if( AMPERE_temp_badData.size > 0 ):
                    print('-----ERROR in GRITI_import_AMPERE_direct-----');
                    print('The file '+os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName)+' was downloaded BAD! There\'s data that\'s just 0\'s - but real data gaps are ommitted, those 0\'s aren\'t real data gaps. Redownload the file from source! CRASHIN\'');
                    sys.crash();
                    # AMPERE_temp_badData = np.append(AMPERE_temp_badData,AMPERE_temp_badData[-1]+1); #add last value b/c diff misses it
                    # AMPERE_temp_badData_refSize = AMPERE_temp['colat'].size;
                    # keyz = list(AMPERE_temp.keys()); #get the current keys in AMPERE_temp
                    # for j in range(0,len(keyz)):
                    #     if( AMPERE_temp[keyz[j]].size == AMPERE_temp_badData_refSize ): #only remove naughty data from full-sized arrays
                    #         AMPERE_temp[keyz[j]] = np.delete(AMPERE_temp[keyz[j]],AMPERE_temp_badData); #remove the naughty data
                    #     #END IF
                    # #END FOR j
                #END IF
                
                #Calc number of datas per time
                AMPERE_temp_goodIndex = np.where(AMPERE_temp['nlon'] > 0)[0][0]; #avoids all 0's data at the beginning of a day, hasn't happened but it could you never know
                AMPERE_temp_numData = np.int64(AMPERE_temp['nlon'][AMPERE_temp_goodIndex]*AMPERE_temp['nlat'][AMPERE_temp_goodIndex]); #number of data per time step
                AMPERE_temp_timeData = AMPERE_temp['colat'].size//AMPERE_temp_numData; #number of time steps
                AMPERE_temp_longDelta = 360/AMPERE_temp['nlon'][AMPERE_temp_goodIndex]; #deg, delta of long pts provided (assume 360 coverage is safe)
                AMPERE_temp_latDelta = np.median(np.diff(np.unique(AMPERE_temp['colat']))); #deg, delta of lat pts provided
                AMPERE_temp_latNum = np.int64(AMPERE_temp['nlat'][AMPERE_temp_goodIndex]);
                AMPERE_temp_longNum = np.int64(AMPERE_temp['nlon'][AMPERE_temp_goodIndex]);
                
                #divine the data rate
                AMPERE_dataRate = np.int64(np.median(np.diff(AMPERE_temp['start_mt']))*60); #sec, data rate
                
                #Fire up a save dict
                AMPERE_save = {
                    'data rate':AMPERE_dataRate, #sec, record data rate (it's in the file name) - not any more
                    'mag orig':
                        {'year':np.ones(AMPERE_temp['colat'].size,dtype=np.int16)*np.int16(dateRange_dayNum_full[i,0]), #yr, record the year b/c it won't change here
                         'dayNum':np.ones(AMPERE_temp['colat'].size,dtype=np.int16)*np.int16(dateRange_dayNum_full[i,1]), #dayNum, record the dayNum b/c it won't change here
                         'data rate':AMPERE_dataRate, #sec, record data rate (it's in the file name) - not any more
                        },
                    }; #hold mag info (its how AMPERE data comes in)
                    
                #list of AMPERE_temp keys to skip when saving to AMPERE_save
                listToSkip = ['nlon', 'nlat', 'nLatGrid', 'nLonGrid', 'start_yr', 'start_mo', 'start_dy', 'start_hr', 'start_mt', 'start_sc',
                                'end_yr', 'end_mo', 'end_dy', 'end_hr', 'end_mt', 'end_sc', 'npnt', 'year', 'doy', 'start_dyN', 'avgint', 
                                'kmax', 'mmax', 'res_deg', 'npnt', 'dBprt', 'db_T', 'db_P', 'db_geo', 'del_db_T', 'del_db_P', 'del_db_geo',
                                'dBnorth1', 'dBnorth2', 'dBeast1', 'dBeast2', 'ddBnorth1', 'ddBnorth2', 'ddBeast1', 'ddBeast2', 'dBr', 'ddBr',
                                'colat', 'mlt_hr', 'cLat_deg', 'geo_cLat_deg', 'geo_lon_deg', 'jPar', 'del_jPar', 'dJr', 'time'];
                
                #Deal with the time - the used functions build time from the data, should deal with a random data outage correctly (if that ever happened)
                AMPERE_save['mag orig']['hour'] = np.repeat(AMPERE_temp['start_hr'],AMPERE_temp_numData).astype(np.int16); #there's AMPERE_temp_numData # of each time step w/ AMPERE_temp_numData # of 0's in a row, then next time step etc..
                AMPERE_save['mag orig']['min'] = np.repeat(AMPERE_temp['start_mt'],AMPERE_temp_numData).astype(np.int16); 
                AMPERE_save['mag orig']['sec'] = np.repeat(AMPERE_temp['start_sc'],AMPERE_temp_numData).astype(np.int16); 
                AMPERE_save['mag orig']['time'] = np.int64(AMPERE_save['mag orig']['hour'])*3600 + np.int64(AMPERE_save['mag orig']['min'])*60 + np.int64(AMPERE_save['mag orig']['sec']); #construct time from time
                # #========== Ensure AMPERE data has all expected times (should only be working on a single day) ==========
                # time_expected = np.arange(0, 86400, AMPERE_save['mag orig']['data rate']); #get the expected times
                # time_got = np.unique(AMPERE_save['mag orig']['time']); #get the times provided
                # if( time_expected.size != time_got.size ):
                #     time_missing_colat = np.tile(np.arange(1,AMPERE_temp_latNum+1,1).astype(AMPERE_temp['colat'].dtype),AMPERE_temp_longNum); #recreate the colat needed to insert
                #     time_missing_mlt = np.tile(np.arange(1,AMPERE_temp_longNum+1,1).astype(AMPERE_temp['mlt'].dtype),AMPERE_temp_latNum); #recreate the mlt needed to insert
                #     time_missing_nan = np.ones(AMPERE_temp_numData,dtype=dataAccuracy)*np.nan; #premake a nan array
                #     time_missing = np.where(np.in1d(time_expected, time_got) == False)[0]; #get the missing times
                #     for j in range(0,time_missing.size):
                #         time_missing_curr = time_expected[time_missing[j]]; #get the time missing
                #         if( time_missing[j] != 0 ):
                #             # time_missing_prev = time_expected[time_missing[j]-1]; #get the time before the missing time
                #             time_missing_insertAt = np.where(time_expected[time_missing[j]-1] == AMPERE_save['mag orig']['time'])[0][-1]+1; #get index to insert at
                #         else:
                #             time_missing_insertAt = 0; #get index to insert at (it's 0 b/c 0 is missing)
                #         #END IF

                #         #assemble times to insert
                #         time_missing_curr_hr = time_missing_curr//3600; #get the hour
                #         time_missing_curr_min = (time_missing_curr - time_missing_curr_hr*3600)//60; #get the min
                #         time_missing_curr_sec = (time_missing_curr - time_missing_curr_hr*3600) - time_missing_curr_min*60; #get the sec
                        
                #         #fix up the time and colat/mlt which are needed for later
                #         AMPERE_save['mag orig']['hour'] = np.insert(AMPERE_save['mag orig']['hour'], time_missing_insertAt, np.ones(AMPERE_temp_numData ,dtype=AMPERE_save['mag orig']['hour'].dtype)*time_missing_curr_hr); #insert
                #         AMPERE_save['mag orig']['min'] = np.insert(AMPERE_save['mag orig']['min'], time_missing_insertAt, np.ones(AMPERE_temp_numData ,dtype=AMPERE_save['mag orig']['min'].dtype)*time_missing_curr_min); #insert
                #         AMPERE_save['mag orig']['sec'] = np.insert(AMPERE_save['mag orig']['sec'], time_missing_insertAt, np.ones(AMPERE_temp_numData ,dtype=AMPERE_save['mag orig']['sec'].dtype)*time_missing_curr_sec); #insert
                #         AMPERE_save['mag orig']['time'] = np.insert(AMPERE_save['mag orig']['time'], time_missing_insertAt, np.ones(AMPERE_temp_numData ,dtype=AMPERE_save['mag orig']['time'].dtype)*time_missing_curr); #insert
                #         AMPERE_save['mag orig']['dayNum'] = np.insert(AMPERE_save['mag orig']['dayNum'], time_missing_insertAt, np.ones(AMPERE_temp_numData ,dtype=AMPERE_save['mag orig']['dayNum'].dtype)*np.int16(dateRange_dayNum_full[i,1])); #insert
                #         AMPERE_save['mag orig']['year'] = np.insert(AMPERE_save['mag orig']['year'], time_missing_insertAt, np.ones(AMPERE_temp_numData ,dtype=AMPERE_save['mag orig']['year'].dtype)*np.int16(dateRange_dayNum_full[i,0])); #insert
                #         AMPERE_temp['colat'] = np.insert(AMPERE_temp['colat'], time_missing_insertAt, time_missing_colat); #insert colat needed
                #         AMPERE_temp['mlt'] = np.insert(AMPERE_temp['mlt'], time_missing_insertAt, time_missing_mlt); #insert MLT needed
                        
                #         keyzFix = list(AMPERE_temp.keys()); #get the current keys in AMPERE_temp    
                #         for jk in range(0,len(keyzFix)):
                #             if( ( AMPERE_temp[keyzFix[jk]].size > AMPERE_temp_timeData ) & (not keyzFix[jk] in listToSkip) ):
                #                 if( AMPERE_temp[keyzFix[jk]].ndim == 1 ):
                #                     AMPERE_temp[keyzFix[jk]] = np.insert(AMPERE_temp[keyzFix[jk]], time_missing_insertAt, time_missing_nan); #insert as needed
                #                 else:
                #                     AMPERE_temp[keyzFix[jk]] = np.insert(AMPERE_temp[keyzFix[jk]], time_missing_insertAt, np.tile(time_missing_nan,3).reshape(AMPERE_temp_numData,3), axis=0); #insert as needed
                #                 #END IF
                #             #END IF
                #         #END FOR jk
                #     #END FOR j
                # #END IF
                
                # AMPERE_save['mag orig']['time'] = np.repeat(np.arange(0,AMPERE_temp_timeData*AMPERE_dataRate,AMPERE_dataRate,dtype=np.int32),AMPERE_temp_numData); #file name ensures 120 is the data rate, there's AMPERE_temp_numData # of each time step w/ AMPERE_temp_numData # of 0's in a row, then next time step etc..
                # AMPERE_save['mag orig']['hour'] = np.int32(AMPERE_save['mag orig']['time']//3600); #get hours
                # AMPERE_save['mag orig']['min'] = np.int32(AMPERE_save['mag orig']['time']//60-AMPERE_save['mag orig']['hour']*60); #get the minutes
                # AMPERE_save['mag orig']['sec'] = np.int32(AMPERE_save['mag orig']['time']-AMPERE_save['mag orig']['min']*60-AMPERE_save['mag orig']['hour']*3600); #get the seconds
                # Calc total sec for the entire time period aligned to the zero hr [debuting year support! I hope it works]
                AMPERE_save['mag orig']['time aligned'] = np.int32(AMPERE_save['mag orig']['dayNum']-dateRange_dayNum_zeroHr[1]-(dateRange_dayNum_zeroHr[0]-AMPERE_save['mag orig']['year'])*subfun_daysInAYear(AMPERE_save['mag orig']['year']))*86400 + \
                    AMPERE_save['mag orig']['time']; #sec, calc total sec for the day range
                AMPERE_save['mag orig']['time'] = np.int32(AMPERE_save['mag orig']['dayNum'])*86400 + AMPERE_save['mag orig']['time']; #sec, calc total sec for the day range
                
                #Load in the mag data that's native to AMPERE raw
                AMPERE_save['mag orig']['lat'] = 90 - AMPERE_temp['colat'];
                #Calculate the geomag data from the mag data
                timeThing = [None for j in range(0,AMPERE_save['mag orig']['time'].size)]; #preallocate
                for j in range(0,AMPERE_save['mag orig']['time'].size):
                    timeThing[j] = datetime.datetime(dateRange_full[i,0], dateRange_full[i,1], dateRange_full[i,2],\
                        hour = AMPERE_save['mag orig']['hour'][j], minute = AMPERE_save['mag orig']['min'][j], second = AMPERE_save['mag orig']['sec'][j]); #date time object for aacgmv2    
                #END FOR j
                AMPERE_save['mag orig']['long'] = dataAccuracy(np.asarray(aacgmv2.convert_mlt(AMPERE_temp['mlt'],timeThing,m2a=True))); #converts from MLT to AACGMV2 coords, whatever they are - this one seems to be vectorized (checked and is correct & vectorized)
                # THIS LINE IS NOT IT, MLT HAS 12 ALWAYS AT SUN! AMPERE_save['mag orig']['long'] = AMPERE_temp['mlt']*15; #file gives mlt, convert to longitude 0 to 360 (east longitude only, if you will)
                
                AMPERE_save['mag orig']['long'][AMPERE_save['mag orig']['long'] >= 180] = AMPERE_save['mag orig']['long'][AMPERE_save['mag orig']['long'] >= 180] - 360; #convert to longitude -180 to 180 (this is OK, AMPERE ref code does this, this is geomagnetic coords tho)
                keyzKeep = list(AMPERE_temp.keys()); #get the current keys in AMPERE_temp    
                for j in range(0,len(keyzKeep)):
                    if( not keyzKeep[j] in listToSkip ):
                        AMPERE_save['mag orig'][keyzKeep[j]] = AMPERE_temp[keyzKeep[j]].copy(); #copy it over if it's not some of the useless headers
                    #END IF
                #END FOR j
                AMPERE_save['mag orig']['data info'] = {
                    'num data per time step':AMPERE_temp_numData,
                    'num time steps':AMPERE_temp_timeData,
                    'lat delta':AMPERE_temp_latDelta,
                    'long delta':AMPERE_temp_longDelta,
                    'lat num':AMPERE_temp_latNum,
                    'long num':AMPERE_temp_longNum,
                    }; #record some essential numbers in case need to backout stuff
                
                AMPERE_save['mag orig']['lat geo'] = np.empty(AMPERE_temp['colat'].size,dtype=dataAccuracy); #preallocate
                AMPERE_save['mag orig']['long geo'] = np.empty(AMPERE_temp['colat'].size,dtype=dataAccuracy);
                for j in range(0,AMPERE_save['mag orig']['time'].size):
                    [AMPERE_save['mag orig']['lat geo'][j], AMPERE_save['mag orig']['long geo'][j], _] = aacgmv2.convert_latlon(AMPERE_save['mag orig']['lat'][j], AMPERE_save['mag orig']['long'][j], 120., timeThing[j], method_code='A2G'); #converts from AACGMV2 to geographic (120. is altitude) - this one is not vectorized big rip
                #END FOR j
                                
                AMPERE_save['mag orig'] = adelphi_modeler(AMPERE_save['mag orig'], AMPERE_hemi, dateRange_full[i,:], dateRange_dayNum_full[i,:]); #calls adelphi code to return data calculated by it (uses default smoothing)
            #END IF
            
            #----- Interpolate values here to new values (will also interpolate geo/mag to different steps if requested) -----
            if( np.all(AMPERE_desired_latLongSteps == None) ):
                # AMPERE_desired_latDelta = AMPERE_temp_latDelta; #no changes, keep orig [defined earlier]
                # AMPERE_desired_longDelta = AMPERE_temp_longDelta; #no changes, keep orig
                # AMPERE_desiredResStr = textNice(AMPERE_desired_latDelta)+'&'+textNice(AMPERE_desired_longDelta); #desired res string for naming
                # AMPERE_save['mag'+AMPERE_desiredResStr] = AMPERE_save['mag orig']; #no changes (not true b/c of MLT making longs wonky, want on a stable grid!)
                AMPERE_save['mag'+AMPERE_desiredResStr] = {}; #prep a dict
                AMPERE_save['geo'+AMPERE_desiredResStr] = {}; #prep a dict
                AMPERE_desired_numData = AMPERE_temp_numData; #set the same
                
                #no need to recalc geo lat/long & time values b/c same as mag orig
                AMPERE_save['geo'+AMPERE_desiredResStr]['lat'] = AMPERE_save['mag orig']['lat']; #copy over these shared lat grid & time values
                AMPERE_desired_latDelta = AMPERE_save['mag orig']['data info']['lat delta']; #desired lat delta (default)
                AMPERE_desired_longDelta = AMPERE_save['mag orig']['data info']['long delta']; #desired long delta (default) long grid is not shared b/c of MLT changing with time
                AMPERE_desired_latNum = AMPERE_save['mag orig']['data info']['lat num']; #number of desired lat data pts
                AMPERE_desired_longNum = AMPERE_save['mag orig']['data info']['long num']; #number of desired long data pts
                AMPERE_save['geo'+AMPERE_desiredResStr]['long'] = np.tile(np.repeat(np.arange(0,360,AMPERE_desired_longDelta),AMPERE_desired_latNum),AMPERE_temp_timeData); #doing it piece-wise will deal with things that wouldn't divide cleanly into 180
                AMPERE_save['geo'+AMPERE_desiredResStr]['long'][AMPERE_save['geo'+AMPERE_desiredResStr]['long'] >= 180] = AMPERE_save['geo'+AMPERE_desiredResStr]['long'][AMPERE_save['geo'+AMPERE_desiredResStr]['long'] >= 180] - 360; #convert to longitude -180 to 180 (this is OK, AMPERE ref code does this, this is geomagnetic coords tho)
                AMPERE_save['geo'+AMPERE_desiredResStr]['year'] = AMPERE_save['mag orig']['year'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['dayNum'] = AMPERE_save['mag orig']['dayNum'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['time'] = AMPERE_save['mag orig']['time'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['hour'] = AMPERE_save['mag orig']['hour'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['min'] = AMPERE_save['mag orig']['min'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['sec'] = AMPERE_save['mag orig']['sec'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['time aligned'] = AMPERE_save['mag orig']['time aligned'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['1D'] = AMPERE_save['mag orig']['1D'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info'] = deepcopy(AMPERE_save['mag orig']['data info']);
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['lat delta'] = AMPERE_desired_latDelta;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['long delta'] = AMPERE_desired_longDelta;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['lat num'] = AMPERE_desired_latNum;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['long num'] = AMPERE_desired_longNum;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['num data per time step'] = np.int64(AMPERE_desired_latNum*AMPERE_desired_longNum);
                
                #gotta interpolate mag too b/c not orig coords (at least this is on a grid!)
                AMPERE_save['mag'+AMPERE_desiredResStr]['lat'] = AMPERE_save['geo'+AMPERE_desiredResStr]['lat']; #copy over these shared lat/long grid & time values
                AMPERE_save['mag'+AMPERE_desiredResStr]['long'] = AMPERE_save['geo'+AMPERE_desiredResStr]['long'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['year'] = AMPERE_save['geo'+AMPERE_desiredResStr]['year'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['dayNum'] = AMPERE_save['geo'+AMPERE_desiredResStr]['dayNum'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['time'] = AMPERE_save['geo'+AMPERE_desiredResStr]['time'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['hour'] = AMPERE_save['geo'+AMPERE_desiredResStr]['hour'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['min'] = AMPERE_save['geo'+AMPERE_desiredResStr]['min'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['sec'] = AMPERE_save['geo'+AMPERE_desiredResStr]['sec'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['time aligned'] = AMPERE_save['geo'+AMPERE_desiredResStr]['time aligned'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['1D'] = AMPERE_save['mag orig']['1D'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info'] = deepcopy(AMPERE_save['mag orig']['data info']);
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['lat delta'] = AMPERE_desired_latDelta;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['long delta'] = AMPERE_desired_longDelta;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['lat num'] = AMPERE_desired_latNum;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['long num'] = AMPERE_desired_longNum;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['num data per time step'] = np.int64(AMPERE_desired_latNum*AMPERE_desired_longNum);
                
                timeThing = [None for j in range(0,AMPERE_save['mag'+AMPERE_desiredResStr]['time'].size)]; #preallocate
                for j in range(0,AMPERE_save['mag'+AMPERE_desiredResStr]['time'].size):
                    timeThing[j] = datetime.datetime(dateRange_full[i,0], dateRange_full[i,1], dateRange_full[i,2],\
                        hour = AMPERE_save['mag'+AMPERE_desiredResStr]['hour'][j], minute = AMPERE_save['mag'+AMPERE_desiredResStr]['min'][j], second = AMPERE_save['mag'+AMPERE_desiredResStr]['sec'][j]); #date time object for aacgmv2    
                #END FOR j
                
                # #Rebase to x y z b/c spherical no go good
                # AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig = geo_to_xyz(AMPERE_save['mag orig']['lat geo'], AMPERE_save['mag orig']['long geo'], np.ones(AMPERE_save['mag orig']['lat geo'].size)*120.); #calc
                # AMPERE_x, AMPERE_y, AMPERE_z = geo_to_xyz(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'], AMPERE_save['geo'+AMPERE_desiredResStr]['long'], np.ones(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'].size)*120.); #calc
            else:
                AMPERE_desired_latDelta = AMPERE_desired_latLongSteps[0]; #desired lat delta
                AMPERE_desired_longDelta = AMPERE_desired_latLongSteps[1]; #desired long delta
                # AMPERE_desiredResStr = textNice(AMPERE_desired_latDelta)+'&'+textNice(AMPERE_desired_longDelta); #desired res string for naming [defined earlier]
                AMPERE_save['mag'+AMPERE_desiredResStr] = {}; #prep a dict
                AMPERE_save['geo'+AMPERE_desiredResStr] = {}; #prep a dict
                
                AMPERE_desired_latNum = np.int64((90-AMPERE_save['mag orig']['lat'].min())/AMPERE_desired_latDelta); #number of desired lat data pts
                AMPERE_desired_longNum = np.int64(360/AMPERE_desired_longDelta); #number of desired long data pts
                AMPERE_desired_numData = np.int64(AMPERE_desired_latNum*AMPERE_desired_longNum); #number of data per time step
                
                AMPERE_save['geo'+AMPERE_desiredResStr]['lat'] = np.tile(np.tile(np.arange(90-AMPERE_desired_latDelta,AMPERE_save['mag orig']['lat'].min()-AMPERE_desired_latDelta,-AMPERE_desired_latDelta),AMPERE_desired_longNum),AMPERE_temp_timeData);
                AMPERE_save['geo'+AMPERE_desiredResStr]['long'] = np.tile(np.repeat(np.arange(0,360,AMPERE_desired_longDelta),AMPERE_desired_latNum),AMPERE_temp_timeData); #doing it piece-wise will deal with things that wouldn't divide cleanly into 180
                AMPERE_save['geo'+AMPERE_desiredResStr]['long'][AMPERE_save['geo'+AMPERE_desiredResStr]['long'] >= 180] = AMPERE_save['geo'+AMPERE_desiredResStr]['long'][AMPERE_save['geo'+AMPERE_desiredResStr]['long'] >= 180] - 360; #convert to longitude -180 to 180 (this is OK, AMPERE ref code does this, this is geomagnetic coords tho)
                # longNew = np.tile(np.repeat(np.vstack((np.arange(0,180,AMPERE_desired_longDelta),np.arange(-180,0,AMPERE_desired_longDelta))),AMPERE_desired_latNum),AMPERE_temp_timeData);
                
                #update the time keeping to match the extra interpolated data pts
                AMPERE_save['geo'+AMPERE_desiredResStr]['year'] = np.ones(AMPERE_save['geo'+AMPERE_desiredResStr]['long'].size,dtype=np.int16)*np.int16(dateRange_dayNum_full[i,0]); #yr, record the year b/c it won't change here
                AMPERE_save['geo'+AMPERE_desiredResStr]['dayNum'] = np.ones(AMPERE_save['geo'+AMPERE_desiredResStr]['long'].size,dtype=np.int16)*np.int16(dateRange_dayNum_full[i,1]); #dayNum, record the dayNum b/c it won't change here
                AMPERE_save['geo'+AMPERE_desiredResStr]['time'] = np.repeat(np.unique(np.int32(AMPERE_save['mag orig']['hour'])*3600+np.int32(AMPERE_save['mag orig']['min'])*60+np.int32(AMPERE_save['mag orig']['sec'])),AMPERE_desired_numData); #file name ensures 120 is the data rate, there's AMPERE_temp_numData # of each time step w/ AMPERE_temp_numData # of 0's in a row, then next time step etc..
                AMPERE_save['geo'+AMPERE_desiredResStr]['hour'] = np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['time']//3600); #get hours
                AMPERE_save['geo'+AMPERE_desiredResStr]['min'] = np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['time']//60-np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['hour'])*60); #get the minutes
                AMPERE_save['geo'+AMPERE_desiredResStr]['sec'] = np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['time']-np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['min'])*60-np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['hour'])*3600); #get the seconds
                # Calc total sec for the entire time period aligned to the zero hr [debuting year support! I hope it works]
                AMPERE_save['geo'+AMPERE_desiredResStr]['time aligned'] = np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['dayNum']-dateRange_dayNum_zeroHr[1]-(dateRange_dayNum_zeroHr[0]-AMPERE_save['geo'+AMPERE_desiredResStr]['year'])*subfun_daysInAYear(AMPERE_save['geo'+AMPERE_desiredResStr]['year']))*86400 + \
                    AMPERE_save['geo'+AMPERE_desiredResStr]['time']; #sec, calc total sec for the day range
                AMPERE_save['geo'+AMPERE_desiredResStr]['time'] = np.int32(AMPERE_save['geo'+AMPERE_desiredResStr]['dayNum'])*86400 + AMPERE_save['geo'+AMPERE_desiredResStr]['time']; #sec, calc total sec for the day range
                AMPERE_save['geo'+AMPERE_desiredResStr]['1D'] = AMPERE_save['mag orig']['1D'];
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info'] = deepcopy(AMPERE_save['mag orig']['data info']);
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['lat delta'] = AMPERE_desired_latDelta;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['long delta'] = AMPERE_desired_longDelta;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['lat num'] = AMPERE_desired_latNum;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['long num'] = AMPERE_desired_longNum;
                AMPERE_save['geo'+AMPERE_desiredResStr]['data info']['num data per time step'] = np.int64(AMPERE_desired_latNum*AMPERE_desired_longNum);
                
                #gotta interpolate mag too b/c not orig coords (at least this is on a grid!)
                AMPERE_save['mag'+AMPERE_desiredResStr]['lat'] = AMPERE_save['geo'+AMPERE_desiredResStr]['lat']; #copy over these shared lat/long grid & time values
                AMPERE_save['mag'+AMPERE_desiredResStr]['long'] = AMPERE_save['geo'+AMPERE_desiredResStr]['long'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['year'] = AMPERE_save['geo'+AMPERE_desiredResStr]['year'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['dayNum'] = AMPERE_save['geo'+AMPERE_desiredResStr]['dayNum'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['time'] = AMPERE_save['geo'+AMPERE_desiredResStr]['time'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['hour'] = AMPERE_save['geo'+AMPERE_desiredResStr]['hour'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['min'] = AMPERE_save['geo'+AMPERE_desiredResStr]['min'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['sec'] = AMPERE_save['geo'+AMPERE_desiredResStr]['sec'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['time aligned'] = AMPERE_save['geo'+AMPERE_desiredResStr]['time aligned'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['1D'] = AMPERE_save['mag orig']['1D'];
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info'] = deepcopy(AMPERE_save['mag orig']['data info']);
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['lat delta'] = AMPERE_desired_latDelta;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['long delta'] = AMPERE_desired_longDelta;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['lat num'] = AMPERE_desired_latNum;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['long num'] = AMPERE_desired_longNum;
                AMPERE_save['mag'+AMPERE_desiredResStr]['data info']['num data per time step'] = np.int64(AMPERE_desired_latNum*AMPERE_desired_longNum);
                
                
                #Calculate the geomag data from the mag data, do it once b/c its big n it hurts
                def parallel_helper(temp_lat, temp_aacgmv2_long, temp_timeThing):
                    temp_AMPERE_interper_geo_latOrig = np.empty(temp_lat.size,dtype=dataAccuracy); #preallocate
                    temp_AMPERE_interper_geo_longOrig = np.empty(temp_lat.size,dtype=dataAccuracy); #preallocate
                    for j in range(0,temp_lat.size):
                        [temp_AMPERE_interper_geo_latOrig[j], temp_AMPERE_interper_geo_longOrig[j], _] = aacgmv2.convert_latlon(temp_lat[j], temp_aacgmv2_long[j], 120., temp_timeThing[j], method_code='A2G'); #converts from AACGMV2 to geographic (120. is altitude) - this one is not vectorized big rip
                    #END FOR j
                    return [temp_AMPERE_interper_geo_latOrig, temp_AMPERE_interper_geo_longOrig]
                #END DEF
                AMPERE_interper_geo_latOrig = np.empty(AMPERE_save['mag'+AMPERE_desiredResStr]['lat'].size,dtype=dataAccuracy); #preallocate
                AMPERE_interper_geo_longOrig = np.empty(AMPERE_save['mag'+AMPERE_desiredResStr]['lat'].size,dtype=dataAccuracy);
                aacgmv2_long = np.empty(AMPERE_save['mag'+AMPERE_desiredResStr]['lat'].size, dtype=dataAccuracy); #preallocate
                timeThing = [None for j in range(0,AMPERE_save['mag'+AMPERE_desiredResStr]['time'].size)]; #preallocate
                for j in range(0,AMPERE_save['mag'+AMPERE_desiredResStr]['time'].size):
                    timeThing[j] = datetime.datetime(dateRange_full[i,0], dateRange_full[i,1], dateRange_full[i,2],\
                        hour = AMPERE_save['mag'+AMPERE_desiredResStr]['hour'][j], minute = AMPERE_save['mag'+AMPERE_desiredResStr]['min'][j], second = AMPERE_save['mag'+AMPERE_desiredResStr]['sec'][j]); #date time object for aacgmv2    
                #END FOR j
                
                #NOTE-this is wronk, can't just do that even if it's nice and easy sadly
                # AMPERE_interper_MLT = AMPERE_save['mag'+AMPERE_desiredResStr]['long'].copy(); #gotta rebuild MLT
                # AMPERE_interper_MLT[AMPERE_interper_MLT < 0] = AMPERE_interper_MLT[AMPERE_interper_MLT < 0] + 360; #convert to 0 to 360
                # AMPERE_interper_MLT = AMPERE_interper_MLT/15; #convert to MLT [0 to 24), never touches 24 (like never touches 360 or 180)
                # ditch this too AMPERE_interper_MLT = aacgmv2.convert_mlt(AMPERE_save['mag'+AMPERE_desiredResStr]['long'],timeThing,m2a=False); #converts from MLT to AACGMV2 coords, whatever they are - this one seems to be vectorized (checked and is correct & vectorized)
                # for j in range(0,AMPERE_save['mag'+AMPERE_desiredResStr]['lat'].size):
                #     [AMPERE_interper_geo_latOrig[j], AMPERE_interper_geo_longOrig[j], _] = aacgmv2.convert_latlon(AMPERE_save['mag'+AMPERE_desiredResStr]['lat'][j], aacgmv2_long[j], 120., timeThing[j], method_code='A2G'); #converts from AACGMV2 to geographic (120. is altitude) - this one is not vectorized big rip
                # #END FOR j
                # tic = time.time()
                with joblib.parallel_backend('loky'):
                    with joblib.Parallel(n_jobs=parallel_numThreads,pre_dispatch=parallel_numThreads,batch_size=1) as parallel_arbiter:    
                        # parallel_list = []; #Prep
                        # parallel_splitterIndexes = np.int64(np.round(np.linspace(0,AMPERE_interper_MLT.size,parallel_numThreads+1,endpoint=True))); #split up the indexes to be parallelized
                        # for j in range(0,parallel_splitterIndexes.size-1):
                        #     parallel_list.append([AMPERE_interper_MLT[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]],timeThing[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]]]);
                        # #END FOR j
                        # aacgmv2_long_list = parallel_arbiter(joblib.delayed(aacgmv2.convert_mlt)(j, k, m2a=True) for j, k in parallel_list); #will this not destroy the world?
                        # for j in range(0,parallel_splitterIndexes.size-1):
                        #     aacgmv2_long[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = aacgmv2_long_list[j]; #load it in
                        # #END FOR j
                        aacgmv2_long = AMPERE_save['mag'+AMPERE_desiredResStr]['long']; #above was useless, it just converted the arbitrary, ideal longitude grid (in aacgmv2 mag coords) to MLT and then back to aacgmv2. nothing of value
                        parallel_splitterIndexes = np.int64(np.round(np.linspace(0,AMPERE_save['mag'+AMPERE_desiredResStr]['long'].size,parallel_numThreads+1,endpoint=True))); #split up the indexes to be parallelized
                        parallel_list = []; #Prep
                        for j in range(0,parallel_splitterIndexes.size-1):
                            parallel_list.append([AMPERE_save['mag'+AMPERE_desiredResStr]['lat'][parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]], aacgmv2_long[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]], timeThing[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]]]);
                        #END FOR j
                        parallel_geo_list = parallel_arbiter(joblib.delayed(parallel_helper)(*sublist) for sublist in parallel_list); #will this not destroy the world?
                        for j in range(0,parallel_splitterIndexes.size-1):
                            AMPERE_interper_geo_latOrig[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_geo_list[j][0]; #load it in
                            AMPERE_interper_geo_longOrig[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_geo_list[j][1]; #load it in
                        #END FOR j
                    #END WITH
                #END WITH
                del parallel_list, parallel_geo_list #save some mem aacgmv2_long_list

                # #calc these all in 1 shot to reduce CPU load later on
                # AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig = geo_to_xyz(AMPERE_interper_geo_latOrig, AMPERE_interper_geo_longOrig, np.ones(AMPERE_interper_geo_longOrig.size)*120.); #calc
                # AMPERE_x, AMPERE_y, AMPERE_z = geo_to_xyz(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'], AMPERE_save['geo'+AMPERE_desiredResStr]['long'], np.ones(AMPERE_save['geo'+AMPERE_desiredResStr]['long'].size)*120.); #calc
               
                # # this doesn't smooth enough, but its perfectly exact interpolation if that's needed ever
                # from RegridderZenMats import RegridderZenMats
                # from RegridderZenAccel import RegridderZenAccel
                # yMatrix, xMatrix = RegridderZenMats((AMPERE_temp_longNum,AMPERE_temp_latNum), (AMPERE_desired_longNum,AMPERE_desired_latNum));
            #END IF
            AMPERE_save['mag'+AMPERE_desiredResStr]['mlt'] = dataAccuracy(np.asarray(aacgmv2.convert_mlt(AMPERE_save['mag'+AMPERE_desiredResStr]['long'],timeThing,m2a=False))); #converts from  AACGMV2 coords to MLT, whatever they are - this one seems to be vectorized (checked and is correct & vectorized)
            
            #----- Prepare a list of dicts ------
            keyzCalc = list(AMPERE_save['mag orig'].keys()); #get the keys available
            parallel_list = []; #prep a parallel list
            for k in range(0,len(keyzCalc)):
                if( not keyzCalc[k] in ['year', 'dayNum', 'data rate', 'time', 'hour', 'min', 'sec', 'time aligned', 'lat', 'long', 'mlt', 'lat geo', 'long geo', 'data info', '1D'] ):
                    AMPERE_saveTemp = {
                        'mag'+AMPERE_desiredResStr:{
                            'lat':AMPERE_save['mag'+AMPERE_desiredResStr]['lat'],
                            'long':AMPERE_save['mag'+AMPERE_desiredResStr]['long'],
                            },
                        'AMPERE_temp_timeData':AMPERE_temp_timeData,
                        'AMPERE_temp_numData':AMPERE_temp_numData,
                        'AMPERE_desired_numData':AMPERE_desired_numData,
                        'AMPERE_desired_latNum':AMPERE_desired_latNum,
                        'AMPERE_desired_longNum':AMPERE_desired_longNum,
                        'AMPERE_desiredResStr':AMPERE_desiredResStr,
                        # 'AMPERE_xOrig':AMPERE_xOrig,
                        # 'AMPERE_yOrig':AMPERE_yOrig,
                        # 'AMPERE_zOrig':AMPERE_zOrig,
                        # 'AMPERE_x':AMPERE_x,
                        # 'AMPERE_y':AMPERE_y,
                        # 'AMPERE_z':AMPERE_z,
                        'AMPERE lat orig rad':AMPERE_interper_geo_latOrig*np.pi/180,
                        'AMPERE long orig rad':AMPERE_interper_geo_longOrig*np.pi/180,
                        'AMPERE lat rad':AMPERE_save['geo'+AMPERE_desiredResStr]['lat']*np.pi/180,
                        'AMPERE long rad':AMPERE_save['geo'+AMPERE_desiredResStr]['long']*np.pi/180,
                        'dataKey':keyzCalc[k],
                        'dataAccuracy':dataAccuracy,
                        'mag orig':{
                            keyzCalc[k]:AMPERE_save['mag orig'][keyzCalc[k]],
                            'lat':AMPERE_save['mag orig']['lat'],
                            'long':AMPERE_save['mag orig']['long'],
                            }, #needed for new mag calcs
                        }; #prime this up
                    
                    parallel_list.append(AMPERE_saveTemp); #tack on this dict of data
                #END IF
            #END FOR k
            
            #----- calc in parallel -----
            # # # AMPERE_desired_results = pool.map(regrid_parallel, parallel_list); #parallel calc the results [this is awfully coded on Windows]
            # # # AMPERE_desired_results = joblib.Parallel(n_jobs=parallel_numCores)(joblib.delayed(regrid_parallel)(i) for i in parallel_list); #will this not destroy the world?
            # # # testing: regrid_parallel(parallel_list[0]);
            
            os.environ['OPENBLAS_NUM_THREADS'] = str(parallel_threadsPerProcess); #set thread override for BLAS and MKL
            os.environ['MKL_NUM_THREADS'] = str(parallel_threadsPerProcess); #set thread override for BLAS and MKL
            with joblib.parallel_backend('loky', inner_max_num_threads=parallel_threadsPerProcess): #should only use the real cores
                AMPERE_desired_results = joblib.Parallel(n_jobs=parallel_numCores)(joblib.delayed(regrid_parallel)(i) for i in parallel_list); #will this not destroy the world?
            #END WITH
            os.environ.pop('OPENBLAS_NUM_THREADS',None); #undo the override
            os.environ.pop('MKL_NUM_THREADS',None); #undo the override
            
            # #non-parallel for debugging
            # AMPERE_desired_results = [];
            # for i in range(0,len(parallel_list)):
            #     i = 12;
            #     AMPERE_desired_results.append(regrid_parallel(parallel_list[i])); #will this not destroy the world?
            # #END WITH
            
            #-----unpack the AMPERE_desired_results -----
            for j in range(0,len(AMPERE_desired_results)):
                AMPERE_save['geo'+AMPERE_desiredResStr][parallel_list[j]['dataKey']] = AMPERE_desired_results[j]['geo'].copy(); #load in geo results
                AMPERE_save['mag'+AMPERE_desiredResStr][parallel_list[j]['dataKey']] = AMPERE_desired_results[j]['mag'].copy(); #load in mag results
            #END FOR j
            del parallel_list, AMPERE_desired_results #clean up mem
            #add in MLT
            
            # #-----this is the loop that's parallelized above-----
            # for k in range(0,len(keyzCalc)):
            #     if( not keyzCalc[k] in ['year', 'dayNum', 'data rate', 'time', 'hour', 'min', 'sec', 'time aligned', 'lat', 'long', 'lat geo', 'long geo', 'data info'] ):
            #         #Did NOT work: so many things
            #         AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]] = np.empty(AMPERE_save['geo'+AMPERE_desiredResStr]['long'].size,dtype=dataAccuracy); #preallocate
            #         AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]] = np.empty(AMPERE_save['geo'+AMPERE_desiredResStr]['long'].size,dtype=dataAccuracy); #preallocate
            #         # try:
            #         #     import geopy.distance as distance
            #         #     distr = distance.distance((lat1,long1),(lat2,long2)).km
            #         # except ImportError:
            #         #     if( FLG_geoWarn == True ):
            #         #         print('---WARNING in GRITI_import_AMPERE_direct---');
            #         #         print('The package geopy isn\'t installed. It has the best alg for geometric distance, so get it if you really want some top notch calcs. Using a minisculely worse alg (and if it fails, somehwat worse[spherical fallback]).');
            #         #         FLG_geoWarn = False; #turn off the first run warning flag
            #         #     #END IF
            #         # from subfun_dist import distWGS84_notBad as distr
            #         # from subfun_dist import dist_haversine as distr #SPEED
            #         # from GRITI_movieMaker_subfun_dataGridder import GRITI_movieMaker_subfun_dataGridder
            #         for j in range(0,AMPERE_temp_timeData):
            #             # for k in range(0,len(keyzKeep)):
            #             #     if( not keyzKeep[k] in ['nlon', 'nlat', 'start_yr', 'start_mo', 'start_dy', 'start_hr', 'start_mt', 'start_sc', 'end_yr', 'end_mo', 'end_dy', 'end_hr', 'end_mt', 'end_sc', 'colat', 'mlt'] ):
            #             #         AMPERE_save['mag'][keyzKeep[k]] = AMPERE_temp[keyzKeep[k]].copy(); #copy it over if it's not some of the useless headers
            #             #     #END IF
            #             # #END FOR k
                        
            #             # AMPERE_gridLatDelta = 1; #hardcoded
            #             # AMPERE_gridLongDelta = 15; #hardcoded
            #             # AMPERE_gridLat = np.arange(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'].min(),AMPERE_save['geo'+AMPERE_desiredResStr]['lat'].max()+AMPERE_gridLatDelta,AMPERE_gridLatDelta);
            #             # AMPERE_gridLong = np.arange(AMPERE_save['geo'+AMPERE_desiredResStr]['long'].min(),AMPERE_save['geo'+AMPERE_desiredResStr]['long'].max()+AMPERE_gridLatDelta,AMPERE_gridLatDelta);
                        
            #             # AMPERE_gridded = GRITI_movieMaker_subfun_dataGridder(AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData],AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData],AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData],AMPERE_gridLat,AMPERE_gridLong,AMPERE_gridLat.size-1,AMPERE_gridLong.size-1,AMPERE_gridLatDelta,AMPERE_gridLongDelta,100,101,0).T                    
            #             # AMPERE_griddedLong, AMPERE_griddedLat = np.meshgrid( AMPERE_gridLong, AMPERE_gridLat); #helps the pcolor work
            #             # kk = ~np.isnan(AMPERE_gridded.flatten()); #keep non-NaNs
            #             # AMPERE_gridded_flat = AMPERE_gridded.flatten();
                        
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = interpolate.griddata(np.vstack((AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], np.vstack((AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T, \
            #             #                                                                                method='linear');
                        
            #             # AMPERE_interper_geo = interpolate.SmoothSphereBivariateSpline((AMPERE_griddedLat.flatten()[kk]+90)*np.pi/180, (AMPERE_griddedLong.flatten()[kk]+180)*np.pi/180, AMPERE_gridded.flatten()[kk], s=1); #s=1 is req for stability, default s=0 does not work
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo((AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+90)*np.pi/180, (AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+180)*np.pi/180, grid=False); #interpolate!
        
            #             # AMPERE_interper_geo = interpolate.LinearNDInterpolator( np.vstack((AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] );
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo( np.vstack((AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T );
                        
            #             #Rebase to x y z b/c spherical I hate
            #             # AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig = geo_to_xyz(AMPERE_save['mag orig']['lat geo'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_save['mag orig']['long geo'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], np.ones(AMPERE_desired_numData)*120.); #calc
            #             # AMPERE_x, AMPERE_y, AMPERE_z = geo_to_xyz(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], np.ones(AMPERE_desired_numData)*120.); #calc
            #             # AMPERE_interper_geo = interpolate.LinearNDInterpolator( np.vstack((AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig)).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] );
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo( np.vstack((AMPERE_x, AMPERE_y, AMPERE_z)).T );
                        
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = interpolate.griddata(np.vstack((AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig)).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], np.vstack((AMPERE_x, AMPERE_y, AMPERE_z)).T, \
            #             #                                                                                 method='linear');
                            
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = interpolate.interp2d(AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig, kind='linear');
                        
            #                 #found that interpolating to smaller grid values w/ original data and not assuming gridding didn't work - just made hot spots where the pts were. using most powerful regridder to get the job done
            #                 # AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = RegridderZenAccel(AMPERE_save['mag orig'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].reshape(AMPERE_temp_longNum,AMPERE_temp_latNum), yMatrix, xMatrix).ravel()*(xMatrix.shape[1]/xMatrix.shape[0])*(yMatrix.shape[0]/yMatrix.shape[1]); #use the gridded values to spread the data
                            
            #                 #RegridderZen interpolates OK but does not smooth enough (it is perfectly - and I mean perfectly - exact tho if that's ever important)
            #                 #this interpolates just right, but it needs a lot of help and padding to prevent NaNs - longitude axis is rolled by 2 to align it correctly, idk why it gets unaligned and idk why its 2 to fix it but idk anymore
            #                 #!!DO NOT KNOW IF 2 KEEPS WORKING FOR OTHER VALUES (other than 3)!!
            #                 AMPERE_interper_dataOrig = AMPERE_save['mag orig'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            #                 AMPERE_interper_latOrig = AMPERE_save['mag orig']['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            #                 AMPERE_interper_longOrig = AMPERE_save['mag orig']['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            #                 AMPERE_interper_longOrig[AMPERE_interper_longOrig<0] = AMPERE_interper_longOrig[AMPERE_interper_longOrig<0]+360; #make it roll around, idk if its req'd but doin it
            #                 kr = AMPERE_interper_longOrig == 0; #get the edge
            #                 AMPERE_interper_longOrig = np.append(AMPERE_interper_longOrig,np.repeat(360,kr.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
            #                 AMPERE_interper_latOrig = np.append(AMPERE_interper_latOrig,AMPERE_interper_latOrig[kr]); #append on lat values (making it cyclical)
            #                 AMPERE_interper_dataOrig = np.append(AMPERE_interper_dataOrig,AMPERE_interper_dataOrig[kr]); #append on data values (making it cyclical)
            #                 kl = AMPERE_interper_longOrig == 345; #get the edge
            #                 AMPERE_interper_longOrig = np.insert(AMPERE_interper_longOrig,0,np.repeat(-15,kl.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
            #                 AMPERE_interper_latOrig = np.insert(AMPERE_interper_latOrig,0,AMPERE_interper_latOrig[kl]); #append on lat values (making it cyclical)
            #                 AMPERE_interper_dataOrig = np.insert(AMPERE_interper_dataOrig,0,AMPERE_interper_dataOrig[kl]); #append on data values (making it cyclical)
            #                 AMPERE_interper_long = AMPERE_save['mag'+AMPERE_desiredResStr]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy();
            #                 AMPERE_interper_long[AMPERE_interper_long<0] = AMPERE_interper_long[AMPERE_interper_long<0]+360; #make it roll around too
            #                 AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = \
            #                     np.roll(interpolate.griddata(np.vstack((AMPERE_interper_longOrig, AMPERE_interper_latOrig)).T, AMPERE_interper_dataOrig, \
            #                         np.vstack((AMPERE_interper_long, AMPERE_save['mag'+AMPERE_desiredResStr]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T,\
            #                         method='linear').reshape(AMPERE_desired_longNum,AMPERE_desired_latNum),2,axis=0).ravel(); #one helluva call to interpolate
                            
            #                 #gotta interpolate mag too b/c not orig coords (at least this is on a grid!) - turns out that's a bad thing kinda.
            #                 #note can't use gridded functions because 2d sphere grids suck nards and its a spherical surface number so 3d grids (x/y/z) don't work b/c no 3rd dim on the data
            #                 #Note that AMPERE_x/y/z are the same b/c geo lat == mag lat b/c it's all on the same arbitrary grid, so no recalc here. xOrig etc need change tho b/c they're coming from a grid instead of the shifted thing for geo
            #                 # AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig = geo_to_xyz(AMPERE_save['mag orig']['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['mag orig']['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], np.ones(AMPERE_temp_numData)*120.); #calc
            #                 # AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig)).T, AMPERE_save['mag orig'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=(AMPERE_desired_numData/AMPERE_temp_numData)**3.5-1, kernel='linear');
            #                 # AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x, AMPERE_y, AMPERE_z)).T); #interpolate!  
                            
            #                 #this does work on spherical data but spline interpolation overshoots a lot - no good here (esp. when regridding to a smaller grid value which should be very constrained)
            #                 # AMPERE_interper_latOrig = (np.arange(AMPERE_save['mag orig']['lat'].min(),90,AMPERE_temp_latDelta)+90)*np.pi/180; #this is flipped, but the data is also flipped lr to compensate (function demands increasing)
            #                 # AMPERE_interper_longOrig = np.arange(0,360,AMPERE_temp_longDelta)*np.pi/180; #this is not quite correct but for interpolation it shouldn't matter - its got the same continuity
            #                 # AMPERE_interper_geo = interpolate.RectSphereBivariateSpline(AMPERE_interper_latOrig, AMPERE_interper_longOrig, np.fliplr(AMPERE_save['mag orig'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].reshape(AMPERE_temp_longNum,AMPERE_temp_latNum)), s=0.0, pole_continuity=False);
            #                 # AMPERE_interper_lat = (np.arange(AMPERE_save['mag orig']['lat'].min(), 90, AMPERE_desired_latDelta)+90)*np.pi/180;
            #                 # AMPERE_interper_long = np.arange(AMPERE_desired_longDelta,360+AMPERE_desired_longDelta,AMPERE_desired_longDelta)*np.pi/180; #this is not quite correct but for interpolation it shouldn't matter - its got the same continuity
            #                 # AMPERE_interper_lat, AMPERE_interper_long = np.meshgrid(AMPERE_interper_lat, AMPERE_interper_long)
            #                 # AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = np.fliplr(AMPERE_interper_geo.ev(AMPERE_interper_lat.ravel(),AMPERE_interper_long.ravel()).reshape(AMPERE_desired_latNum, AMPERE_desired_longNum)).ravel(); #undo the flip and flatten it out
            #             #END IF
            #             if( np.all(AMPERE_desired_latLongSteps == None) ):
            #                 #Rebase to x y z b/c spherical no go good
            #                 # AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig = geo_to_xyz(AMPERE_save['mag orig']['lat geo'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['mag orig']['long geo'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], np.ones(AMPERE_temp_numData)*120.); #calc
            #                 # AMPERE_x, AMPERE_y, AMPERE_z = geo_to_xyz(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], np.ones(AMPERE_desired_numData)*120.); #calc
            #                 AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_yOrig[j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_zOrig[j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T, AMPERE_save['mag orig'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=0, kernel='linear');
            #                 AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
            #             else:
            #                 #found that interpolating to smaller grid values w/ original data and not assuming gridding didn't work - just made hot spots where the pts were. using most powerful regridder to get the job done
            #                 # AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = RegridderZenAccel(AMPERE_interper_dataOrig, yMatrix, xMatrix).ravel()*(xMatrix.shape[1]/xMatrix.shape[0])*(yMatrix.shape[0]/yMatrix.shape[1]); #use the gridded values to spread the data, multiply to renormalize so values are the same as original (this isn't finite flux spreading, its interpolating measurements across an area w/ no measurements)
                            
            #                 # AMPERE_interper_latOrig = np.arange(90-AMPERE_temp_latDelta,AMPERE_save['mag orig']['lat'].min()-AMPERE_temp_latDelta,-AMPERE_temp_latDelta); #this is flipped, but the data is also flipped lr to compensate (function demands increasing)
            #                 # AMPERE_interper_longOrig = np.arange(0,360,AMPERE_temp_longDelta); #this is not quite correct but for interpolation it shouldn't matter - its got the same continuity
            #                 # AMPERE_interper_longOrig, AMPERE_interper_latOrig = np.meshgrid(AMPERE_interper_longOrig, AMPERE_interper_latOrig); #mesh grid to get all pts possible
            #                 # AMPERE_interper_lat = np.arange(90-AMPERE_desired_latDelta,AMPERE_save['mag orig']['lat'].min()-AMPERE_desired_latDelta,-AMPERE_desired_latDelta);
            #                 # AMPERE_interper_long = np.arange(0,360,AMPERE_desired_longDelta); #this is not quite correct but for interpolation it shouldn't matter - its got the same continuity
            #                 # AMPERE_interper_long, AMPERE_interper_lat = np.meshgrid(AMPERE_interper_long, AMPERE_interper_lat); #mesh grid to get all pts possible
                            
            #                 # #RegridderZen interpolates OK but does not smooth enough (it is perfectly - and I mean perfectly - exact tho if that's ever important)
            #                 # #this interpolates just right, but it needs a lot of help and padding to prevent NaNs - longitude axis is rolled by 2 to align it correctly, idk why it gets unaligned and idk why its 2 to fix it but idk anymore
            #                 # #!!DO NOT KNOW IF 2 KEEPS WORKING FOR OTHER VALUES (other than 3)!!
            #                 # AMPERE_interper_dataOrig = AMPERE_save['geo1&15'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            #                 # AMPERE_interper_latOrig = AMPERE_save['geo1&15']['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            #                 # AMPERE_interper_longOrig = AMPERE_save['geo1&15']['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            #                 # AMPERE_interper_longOrig[AMPERE_interper_longOrig<0] = AMPERE_interper_longOrig[AMPERE_interper_longOrig<0]+360; #make it roll around, idk if its req'd but doin it
            #                 # kr = AMPERE_interper_longOrig == 0; #get the edge
            #                 # AMPERE_interper_longOrig = np.append(AMPERE_interper_longOrig,np.repeat(360,kr.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
            #                 # AMPERE_interper_latOrig = np.append(AMPERE_interper_latOrig,AMPERE_interper_latOrig[kr]); #append on lat values (making it cyclical)
            #                 # AMPERE_interper_dataOrig = np.append(AMPERE_interper_dataOrig,AMPERE_interper_dataOrig[kr]); #append on data values (making it cyclical)
            #                 # kl = AMPERE_interper_longOrig == 345; #get the edge
            #                 # AMPERE_interper_longOrig = np.insert(AMPERE_interper_longOrig,0,np.repeat(-15,kl.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
            #                 # AMPERE_interper_latOrig = np.insert(AMPERE_interper_latOrig,0,AMPERE_interper_latOrig[kl]); #append on lat values (making it cyclical)
            #                 # AMPERE_interper_dataOrig = np.insert(AMPERE_interper_dataOrig,0,AMPERE_interper_dataOrig[kl]); #append on data values (making it cyclical)
            #                 # AMPERE_interper_long = AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy();
            #                 # AMPERE_interper_long[AMPERE_interper_long<0] = AMPERE_interper_long[AMPERE_interper_long<0]+360; #make it roll around too
            #                 # AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = \
            #                 #     np.roll(interpolate.griddata(np.vstack((AMPERE_interper_longOrig, AMPERE_interper_latOrig)).T, AMPERE_interper_dataOrig, \
            #                 #         np.vstack((AMPERE_interper_long, AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T,\
            #                 #         method='linear').reshape(AMPERE_desired_longNum,AMPERE_desired_latNum),2,axis=0).ravel(); #one helluva call to interpolate
                                    
            #                 #this special mode is based off of the mag just calc'd and uses the same gridding function as previous
            #                 #Rebase to x y z b/c spherical no go good
            #                 # AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], neighbors=AMPERE_temp_numData*AMPERE_temp_timeData, smoothing=0, kernel='linear');
            #                 # AMPERE_interper_splittr = np.append(np.arange(0,AMPERE_desired_numData,AMPERE_temp_numData*AMPERE_temp_timeData),AMPERE_desired_numData); #split up the call to prevent sys crashes
            #                 # for jk in range(0,AMPERE_interper_splittr.size-1):
            #                 #     AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData+AMPERE_interper_splittr[jk]:j*AMPERE_desired_numData+AMPERE_interper_splittr[jk+1]] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData+AMPERE_interper_splittr[jk]:j*AMPERE_desired_numData+AMPERE_interper_splittr[jk+1]], AMPERE_y[j*AMPERE_desired_numData+AMPERE_interper_splittr[jk]:j*AMPERE_desired_numData+AMPERE_interper_splittr[jk+1]], AMPERE_z[j*AMPERE_desired_numData+AMPERE_interper_splittr[jk]:j*AMPERE_desired_numData+AMPERE_interper_splittr[jk+1]])).T); #interpolate!  
            #                 # #END FOR jk
            #                 #no limits
            #                 AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_save['mag'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], smoothing=0, kernel='linear');
            #                 AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
            #                 #******WARNING****** THIS CODE MAY CRASH YOUR COMPUTER LIKE STRAIGHT SHUT OFF AND TURN BACK ON (BC COMP NOT 100% STABLE) BUT ITS DONE IT ON AN i5-4690K SYSTEM AND FX-8530 SYSTEM (i7-7700K was fine)******************************************************
                            
            #                 #!! going with just straight regrid, but it is not as smooth as could be - need to determine a dynamic smoothing factor to smooth out !!
            #                 #turns out that there should be more smoothing than straight regridding can handle (no matter how exact it is)
            #                 # AMPERE_x, AMPERE_y, AMPERE_z = geo_to_xyz(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], np.ones(AMPERE_desired_numData)*120.); #calc
            #                 # AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_x, AMPERE_y, AMPERE_z)).T, AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], smoothing=AMPERE_desired_numData/AMPERE_temp_numData**2, kernel='thin_plate_spline');
            #                 # AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x, AMPERE_y, AMPERE_z)).T); #interpolate!  
                            
            #                 # AMPERE_interper_latOrig = np.arange(90-AMPERE_temp_latDelta,AMPERE_save['mag orig']['lat'].min()-AMPERE_temp_latDelta,-AMPERE_temp_latDelta); #this is flipped, but the data is also flipped lr to compensate (function demands increasing)
            #                 # AMPERE_interper_longOrig = np.arange(0,360,AMPERE_temp_longDelta); #this is not quite correct but for interpolation it shouldn't matter - its got the same continuity
            #                 # AMPERE_interper_lat = np.arange(90-AMPERE_desired_latDelta,AMPERE_save['mag orig']['lat'].min()-AMPERE_desired_latDelta,-AMPERE_desired_latDelta);
            #                 # AMPERE_interper_long = np.arange(0,360,AMPERE_desired_longDelta); #this is not quite correct but for interpolation it shouldn't matter - its got the same continuity
            #                 # AMPERE_interper_lat, AMPERE_interper_long = np.meshgrid(AMPERE_interper_lat, AMPERE_interper_long); #mesh grid to get all pts possible
            #                 # AMPERE_interper_geo = interpolate.RectBivariateSpline(AMPERE_interper_longOrig, AMPERE_interper_latOrig, AMPERE_save['mag orig'][keyzCalc[k]][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].reshape(AMPERE_temp_longNum,AMPERE_temp_latNum), kx=1, ky=1, s=0); #interpolate the grid (linearly!)
            #                 # AMPERE_save['geo'+AMPERE_desiredResStr][keyzCalc[k]][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo.ev(AMPERE_interper_lat.ravel(),AMPERE_interper_long.ravel())#.reshape(AMPERE_desired_longNum, AMPERE_desired_latNum, ); #get all the pts needed
            #             #END IF
                        
            #             #Rebase (this idd not work, but personalized rebasing did work)
            #             # center_lat = 90;
            #             # center_long = 0;
            #             # dist_sphere = distr(np.array(center_lat),np.array(center_long),AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData],AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]);
            #             # azimuth_sphere = np.arctan2(np.cos(center_lat*np.pi/180) * np.sin((center_long - AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])*np.pi/180),np.cos(AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) * np.sin(center_lat*np.pi/180) - np.sin(AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) * np.cos(center_lat*np.pi/180) * np.cos((center_long - AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])*np.pi/180));
            #             # x_centeredOrig = dist_sphere*np.cos(azimuth_sphere);
            #             # y_centeredOrig = dist_sphere*np.sin(azimuth_sphere);
            #             # # AMPERE_interper_geo = interpolate.interp2d(x_centeredOrig, y_centeredOrig, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], kind='linear');
            #             # #Rebase
            #             # dist_sphere = distr(np.array(center_lat),np.array(center_long),AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData],AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]);
            #             # azimuth_sphere = np.arctan2(np.cos(center_lat*np.pi/180) * np.sin((center_long - AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])*np.pi/180),np.cos(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) * np.sin(center_lat*np.pi/180) - np.sin(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) * np.cos(center_lat*np.pi/180) * np.cos((center_long - AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])*np.pi/180));
            #             # x_centered = dist_sphere*np.cos(azimuth_sphere);
            #             # y_centered = dist_sphere*np.sin(azimuth_sphere);
            #             # # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo(x_centered, y_centered);
            #             # # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = interpolate.griddata(np.vstack((x_centeredOrig, y_centeredOrig)).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], np.vstack((x_centered, y_centered)).T, \
            #             # #                                                                                 method='linear');
                        
            #             # AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((x_centeredOrig, y_centeredOrig)).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=0, kernel='linear');
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo(np.vstack((x_centered, y_centered)).T); #interpolate!
                            
            #             # # this worked but slowwww
            #             # for pp in range(0,AMPERE_temp_numData):                   
            #             #     dist_sphere = distr(np.array(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData+pp]),np.array(AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData+pp]),AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData],AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])
            #             #     azimuth_sphere = np.arctan2(np.cos(AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) * np.sin((AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] - AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData+pp])*np.pi/180),np.cos(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData+pp]*np.pi/180) * np.sin(AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) - np.sin(AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData+pp]*np.pi/180) * np.cos(AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]*np.pi/180) * np.cos((AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] - AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData+pp])*np.pi/180));
            #             #     x_centered = dist_sphere*np.cos(azimuth_sphere);
            #             #     y_centered = dist_sphere*np.sin(azimuth_sphere);
                            
            #             #     # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData+pp] = interpolate.griddata( (x_centered,y_centered), AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], (0,0), \
            #             #     #                                                                                 method='linear');
                                
            #             #     AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((x_centered, y_centered)).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=0, kernel='linear');
            #             #     AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData+pp] = AMPERE_interper_geo(np.vstack((0, 0)).T); #interpolate!
                            
            #             #     # AMPERE_interper_geo = interpolate.interp1d(dist_sphere, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], kind='linear', fill_value='extrapolate');
            #             #     # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData+pp] = AMPERE_interper_geo(0); #calculate the extrapolated value at 0
            #             # # #END FOR pp
                        
            #             # AMPERE_interper_geo = interpolate.SmoothSphereBivariateSpline((AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+90)*np.pi/180, (AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+180)*np.pi/180, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], s=150); #s=1 is req for stability, default s=0 does not work
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo((AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+90)*np.pi/180, (AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+180)*np.pi/180, grid=False); #interpolate!
                        
            #             # AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T, AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=0, kernel='linear');
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo(np.vstack((AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData])).T); #interpolate!
            
            #             # AMPERE_interper_geo = interpolate.Rbf((AMPERE_save['geo'+AMPERE_desiredResStr]['lat orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+90)*np.pi/180, (AMPERE_save['geo'+AMPERE_desiredResStr]['long orig'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+180)*np.pi/180, np.ones(AMPERE_temp_numData,dtype=dataAccuracy), AMPERE_save['mag']['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]);
            #             # AMPERE_save['geo'+AMPERE_desiredResStr]['Jr'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData] = AMPERE_interper_geo((AMPERE_save['geo'+AMPERE_desiredResStr]['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+90)*np.pi/180, (AMPERE_save['geo'+AMPERE_desiredResStr]['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData]+180)*np.pi/180, np.ones(AMPERE_temp_numData,dtype=dataAccuracy));
            #         #END FOR j
            #     #END IF
            # #END FOR k
            
            #--- Save the data as HDF5 ---
            #0 is the hdf5 format goal
            AMPERE_fileName_write = folder_fileNameFormat[ 0 ].replace('#DN', str(dateRange_dayNum_full[i,1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#YR', str(dateRange_dayNum_full[i,0]) ); #replace any #YR with the current year
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#MO', str(dateRange_full[i,1]).zfill(2) ); #replace any #MO with the current month
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#DY', str(dateRange_full[i,2]).zfill(2) ); #replace any #DY with current day
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#HEMI', AMPERE_hemi ); #replace any #HEMI with current hemisphere
            
            if( not os.path.isfile( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName_write) ) ):
                #now save the data to use again
                # keyz = list(AMPERE_save.keys()); #keys to the dict
                # h5pyChunkShape = AMPERE_save['hour'].shape; #get the shape of one of the vectors and use it as a chunk (read only whole chunks)
                with h5py.File( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName_write), 'w') as AMPERE_file:
                    subfun_h5_writer(AMPERE_file, AMPERE_save, h5_path_prefix=None, keyz_ignore=None); #does it all
                    # for j in range(0,len(keyz)):
                    #     if( np.isscalar(AMPERE_save[keyz[j]]) == False ):
                    #         if( not type(AMPERE_save[keyz[j]]) is dict ):
                    #             AMPERE_file.create_dataset(keyz[j], data=AMPERE_save[keyz[j]], chunks=AMPERE_save[keyz[j]].shape, compression="gzip", compression_opts=9, shuffle=True, fletcher32=True); #write that data
                    #         else:
                    #             sub_AMPERE_file = AMPERE_file.create_group(keyz[j]); #create a handle for a subgroup
                    #             sub_keyz = list(AMPERE_save[keyz[j]].keys()); #keys to the dict
                    #             for k in range(0,len(sub_keyz)):
                    #                 if( np.isscalar(AMPERE_save[keyz[j]][sub_keyz[k]]) == False ):
                    #                     if( not type(AMPERE_save[keyz[j]][sub_keyz[k]]) is dict ):
                    #                         sub_AMPERE_file.create_dataset(sub_keyz[k], data=AMPERE_save[keyz[j]][sub_keyz[k]], chunks=AMPERE_save[keyz[j]][sub_keyz[k]].shape, compression="gzip", compression_opts=9, shuffle=True, fletcher32=True); #write that data
                    #                     else:
                    #                         subsub_AMPERE_file = sub_AMPERE_file.create_group(sub_keyz[k]); #create a handle for a subgroup
                    #                         subsub_keyz = list(AMPERE_save[keyz[j]][sub_keyz[k]].keys()); #keys to the dict
                    #                         for kj in range(0,len(subsub_keyz)):
                    #                             if( np.isscalar(AMPERE_save[keyz[j]][sub_keyz[k]][subsub_keyz[kj]]) == False ):
                    #                                 subsub_AMPERE_file.create_dataset(subsub_keyz[kj], data=AMPERE_save[keyz[j]][sub_keyz[k]][subsub_keyz[kj]], chunks=AMPERE_save[keyz[j]][sub_keyz[k]][subsub_keyz[kj]].shape, compression="gzip", compression_opts=9, shuffle=True, fletcher32=True); #write that data
                    #                             else:
                    #                                 subsub_AMPERE_file.attrs[subsub_keyz[kj]] = AMPERE_save[keyz[j]][sub_keyz[k]][subsub_keyz[kj]]; #save the attribute
                    #                             #END IF
                    #                         #END FOR k
                    #                     #END IF
                    #                 else:
                    #                     sub_AMPERE_file.attrs[sub_keyz[k]] = AMPERE_save[keyz[j]][sub_keyz[k]]; #save the attribute
                    #                 #END IF
                    #             #END FOR k
                    #         #END IF
                    #     else:
                    #         #if size 1, add it as an attribute
                    #         AMPERE_file.attrs[keyz[j]] = AMPERE_save[keyz[j]]; #save the attribute
                    #     #END IF
                    # #END FOR j
                    #add on non-data-related attributes
                    AMPERE_file.attrs['version'] = version_alg; #save the attribute
                #END WITH
            else:
                #now save the data to use again (updated this time!)
                with h5py.File( os.path.join(settings_paths['data'], folder_AMPERE, str(dateRange_full[i,0]), AMPERE_fileName_write), 'a') as AMPERE_file:
                    subfun_h5_writer(AMPERE_file, AMPERE_save, h5_path_prefix=None, keyz_ignore=['mag orig', 'data rate']); #does it all
                #END WITH
            #END IF
            
            #--- Copy to main data var as well (that gets returned) ---
            AMPERE_data = dictDelver_saver(AMPERE_save, dict2save_subDict=AMPERE_coordType+AMPERE_desiredResStr, diction=AMPERE_data, FLG_saveTopLevelNonDicts=True);
            
            # keyz = list(AMPERE_save.keys()); #keys to the dict
            # if( FLG_firstFileMade == True ):
            #     #Pick out the coord type we want   
            #     for j in range(0,len(keyz)):
            #         if( not type(AMPERE_save[keyz[j]]) is dict ):
            #             if( np.isscalar(AMPERE_save[keyz[j]]) == False ):
            #                 AMPERE_data[keyz[j]] = [AMPERE_save[keyz[j]]]; #copy it over
            #             else:
            #                 AMPERE_data[keyz[j]] = AMPERE_save[keyz[j]]; #copy it over, don't wrap in list
            #             #END IF
            #         elif( keyz[j] == AMPERE_coordType+AMPERE_desiredResStr ):
            #             sub_keyz = list(AMPERE_save[keyz[j]].keys()); #keys to the dict
            #             for k in range(0,len(sub_keyz)):
            #                 if( type(AMPERE_save[keyz[j]][sub_keyz[k]]) is dict ):
            #                     AMPERE_data[sub_keyz[k]] = AMPERE_save[keyz[j]][sub_keyz[k]]; #copy it over, don't wrap in list
            #                 else:
            #                     AMPERE_data[sub_keyz[k]] = [AMPERE_save[keyz[j]][sub_keyz[k]]]; #copy it over
            #                 #END IF
            #             #END FOR k
            #         #END IF
            #     #END FOR j
            #     FLG_firstFileMade = False; #turn off first file flag
            # else:
            #     for j in range(0,len(keyz)):
            #         if( not type(AMPERE_save[keyz[j]]) is dict ):
            #             if( np.sum(strfind(list(AMPERE_data.keys()),keyz[j])) > 0 ): #make sure key exists
            #                 if( np.isscalar(AMPERE_save[keyz[j]]) == False ):
            #                     AMPERE_data[keyz[j]].append(AMPERE_save[keyz[j]]); #tack on
            #                 else:
            #                     if( np.isclose(AMPERE_save[keyz[j]], AMPERE_data[keyz[j]]) == False ):
            #                         #only worry is if the attribute isn't consistent
            #                         print('-----Warning-----');
            #                         print('Attribute '+keyz[j]+' isn\'t the same as the previously recorded value from another file of '+ \
            #                             str(AMPERE_data[keyz[j]])+' and this file\'s value of '+str(AMPERE_save[keyz[j]])+ \
            #                             '.\n NaN\'ing it and try to sort that out.');
            #                         AMPERE_data[keyz[j]] = np.nan; #nan that attribute, figure it out later
            #                     #END IF
            #                 #END IF
            #             else:
            #                 if( np.isscalar(AMPERE_save[keyz[j]]) == False ):
            #                     AMPERE_data[keyz[j]] = [AMPERE_save[keyz[j]]]; #record if doesn't exist
            #                 else:
            #                     AMPERE_data[keyz[j]] = AMPERE_save[keyz[j]]; #copy it over
            #                 #END IF
            #             #END IF
            #         elif( keyz[j] == AMPERE_coordType+AMPERE_desiredResStr ):
            #             sub_keyz = list(AMPERE_save[keyz[j]].keys()); #keys to the dict
            #             for k in range(0,len(sub_keyz)):
            #                 if( np.sum(strfind(list(AMPERE_data.keys()),sub_keyz[k])) > 0 ): #make sure key exists
            #                     if( np.isscalar(AMPERE_save[keyz[j]][sub_keyz[k]]) == False ):
            #                         AMPERE_data[sub_keyz[k]].append(AMPERE_save[keyz[j]][sub_keyz[k]]); #tack on
            #                     else:
            #                         if( np.isclose(AMPERE_save[keyz[j]][sub_keyz[k]], AMPERE_data[sub_keyz[k]]) == False ):
            #                             #only worry is if the attribute isn't consistent
            #                             print('-----Warning-----');
            #                             print('Attribute '+sub_keyz[k]+' isn\'t the same as the previously recorded value from another file of '+ \
            #                                 str(AMPERE_data[sub_keyz[k]])+' and this file\'s value of '+str(AMPERE_save[keyz[j]][sub_keyz[k]])+ \
            #                                 '.\n NaN\'ing it and try to sort that out.');
            #                             AMPERE_data[sub_keyz[k]] = np.nan; #nan that attribute, figure it out later
            #                         #END IF
            #                     #END IF
            #                 else:
            #                     if( np.isscalar(AMPERE_save[keyz[j]][sub_keyz[k]]) == False ):
            #                         AMPERE_data[sub_keyz[k]] = [AMPERE_save[keyz[j]][sub_keyz[k]]]; #record if doesn't exist
            #                     else:
            #                         if( np.isclose(AMPERE_save[keyz[j]][sub_keyz[k]], AMPERE_data[sub_keyz[k]]) == False ):
            #                             #only worry is if the attribute isn't consistent
            #                             print('-----Warning-----');
            #                             print('Attribute '+sub_keyz[k]+' isn\'t the same as the previously recorded value from another file of '+ \
            #                                 str(AMPERE_data[sub_keyz[k]])+' and this file\'s value of '+str(AMPERE_save[keyz[j]][sub_keyz[k]])+ \
            #                                 '.\n NaN\'ing it and try to sort that out.');
            #                             AMPERE_data[sub_keyz[k]] = np.nan; #nan that attribute, figure it out later
            #                         #END IF
            #                     #END IF
            #                 #END IF
            #             #END FOR k
            #         #END IF                
            #     #END FOR j
            # #END IF
    
            #clear AMPERE_data_par AMPERE_data_temp
            toc = time.time() - tic; #end timing
            print("AMPERE data parsing and re-saving as HDF5 took: "+str(np.round(toc,2))+" sec / "+str(np.round(toc/60,2))+" min");
            AMPERE_timeEst = toc/60; #record
            AMPERE_timeRuns += 1; #increment the runs
        #END IF
    #END FOR i
    #--- Convert to numpy arrays from the dynamically added lists ---
    # keyz = list(AMPERE_data.keys()); #get the current keys
    # for j in range(0,len(keyz)):
    #     if( isinstance(AMPERE_data[keyz[j]], dict) ):
    #         keyz_sub = list(AMPERE_data[keyz[j]].keys()); #can go one down, make afunction if more is needed and recurse it (see subfun_h5_reader)
    #         for kj in range(0,len(keyz_sub)):
    #             if( np.isscalar(AMPERE_data[keyz[j]][keyz_sub[kj]]) == False ):
    #                 if( AMPERE_data[keyz[j]][keyz_sub[kj]][0].ndim == 1 ):
    #                     AMPERE_data[keyz[j]][keyz_sub[kj]] = np.hstack(AMPERE_data[keyz[j]][keyz_sub[kj]]); #stack em if needed
    #                 else:
    #                     AMPERE_data[keyz[j]][keyz_sub[kj]] = np.vstack(AMPERE_data[keyz[j]][keyz_sub[kj]]); #stack em if needed
    #                 #END IF
    #             #END IF
    #         #END FOR kj
    #     elif( np.isscalar(AMPERE_data[keyz[j]]) == False ):
    #         #if not a scalar, apply the logical mask
    #         if( AMPERE_data[keyz[j]][0].ndim == 1 ):
    #             AMPERE_data[keyz[j]] = np.hstack(AMPERE_data[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)   
    #         else:
    #             AMPERE_data[keyz[j]] = np.vstack(AMPERE_data[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
    #         #END IF
    #     #END IF
    # #END FOR j
    
    def dictDelver_convert_custom(dictDelver):
        # dictDelver_alias = dictDelver; #make an alias in here
        for dic in dictDelver: #delve through all
            if( isinstance(dictDelver[dic], dict) ):
                dictDelver[dic] = dictDelver_convert_custom(dictDelver[dic]); #delve deeper
            elif( np.isscalar(dictDelver[dic]) == False ):
                if( isinstance(dictDelver[dic][0], dict) ):
                    #this is hacky but I'm tired
                    dictDelver[dic+'tmp'] = deepcopy(dictDelver[dic]); #make a copy
                    keyz = list(dictDelver[dic][0].keys()); #yoink
                    dictDelver[dic] = {}; #reset
                    for j in range(0,len(keyz)):
                        for jk in range(0,len(dictDelver[dic+'tmp'])):
                            if( keyz[j] in dictDelver[dic] ):
                                if( dictDelver[dic+'tmp'][jk][keyz[j]].ndim == 1 ):
                                    dictDelver[dic][keyz[j]] = np.hstack( (dictDelver[dic][keyz[j]], dictDelver[dic+'tmp'][jk][keyz[j]]) ); #convert from list to array (cause list can be stacked fast but is awful for using)   
                                else:
                                    dictDelver[dic][keyz[j]] = np.vstack( (dictDelver[dic][keyz[j]], dictDelver[dic+'tmp'][jk][keyz[j]]) ); #convert from list to array (cause list can be stacked fast but is awful for using)
                                #END IF
                            else:
                                dictDelver[dic][keyz[j]] = dictDelver[dic+'tmp'][jk][keyz[j]]; #load in
                            #END IF
                        #END FOR jk
                    #END FOR j
                else:
                    #this is the custom thing to convert a bunch of lists back into useful arrays
                    if( dictDelver[dic][0].ndim == 1 ):
                        dictDelver[dic] = np.hstack(dictDelver[dic]); #convert from list to array (cause list can be stacked fast but is awful for using)
                    else:
                        dictDelver[dic] = np.vstack(dictDelver[dic]); #convert from list to array (cause list can be stacked fast but is awful for using)
                    #END IF
                #END IF
            #END IF
        #END FOR dic
        return dictDelver
    #END DEF
    AMPERE_data = dictDelver_convert_custom(AMPERE_data); #badabing badaboom
    
    
    # #Clear out AMPERE data less than the minimum
    # k = np.where(np.min(AMPERE_jouleHeating_plotLimValu) > AMPERE_data[:,AMPERE_jouleHeating_pos])[0]; #find entries less than the min plotting number (clear it up)
    # AMPERE_data = np.delete(AMPERE_data,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    # #AMPERE_time = np.delete(AMPERE_time,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    # #AMPERE_lat = np.delete(AMPERE_lat,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    # #AMPERE_long = np.delete(AMPERE_long,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    
    #Pick out the coord type we want
    # AMPERE_data = {}; #prep a list
    # keyzKeep = list(AMPERE_save.keys()); #get the current keys in AMPERE_temp    
    # for j in range(0,len(keyzKeep)):
    #     if( not type(AMPERE_save[keyzKeep[j]]) is dict ):
    #         AMPERE_data[keyzKeep[j]] = AMPERE_save[keyzKeep[j]].copy(); #copy it over
    #     #END IF
    # #END FOR j
    # keyzKeep = list(AMPERE_save[AMPERE_coordType].keys()); #get the current keys in AMPERE_temp    
    # for j in range(0,len(keyzKeep)):
    #     AMPERE_data[keyzKeep[j]] = AMPERE_save[AMPERE_coordType][keyzKeep[j]].copy(); #copy it over
    # #END FOR j
        
    #Alias JR in to field-aligned current
    AMPERE_data['FAC'] = AMPERE_data['Jr']; #it's the same
    AMPERE_data['FAC +'] = np.copy(AMPERE_data['FAC']); #copy it over
    AMPERE_data['FAC +'][AMPERE_data['FAC'] < 0] = np.nan; #nan stuff that's negative
    AMPERE_data['FAC -'] = np.copy(AMPERE_data['FAC']); #copy it over
    AMPERE_data['FAC -'][AMPERE_data['FAC'] > 0] = np.nan; #nan stuff that's positive
    AMPERE_data['FAC Abs'] = np.abs(AMPERE_data['FAC']); #abs it
    
    AMPERE_data['Jr_out +'] = np.copy(AMPERE_data['Jr_out']); #copy it over
    AMPERE_data['Jr_out +'][AMPERE_data['Jr_out'] < 0] = np.nan; #nan stuff that's negative
    AMPERE_data['Jr_out -'] = np.copy(AMPERE_data['FAC']); #copy it over
    AMPERE_data['Jr_out -'][AMPERE_data['Jr_out'] > 0] = np.nan; #nan stuff that's positive
    AMPERE_data['Jr_out Abs'] = np.abs(AMPERE_data['Jr_out']); #abs it
    
    AMPERE_data['Jr_diff'] = AMPERE_data['Jr'] - AMPERE_data['Jr_out']; #the diff
    AMPERE_data['Jr_diff Abs'] = np.abs(AMPERE_data['Jr_diff']); #the abs of the diff
    #Get time ready
    AMPERE_data['time unique'] =  np.unique(AMPERE_data['time']); #sec, get unique times (v useful)
    #Make sure JH is always positive
    AMPERE_data['JH'][AMPERE_data['JH'] < 0] = 0; #0 negative (negative comes from regridding near-0 values me thinks)
    
    # #quick adjustment to dictionary holder method
    # AMPERE = {
    #     'time':AMPERE_data[:,5],
    #     'lat':AMPERE_data[:,6],
    #     'long':AMPERE_data[:,7],
    #     'Pedersen':AMPERE_data[:,0],
    #     'Hall':AMPERE_data[:,1],
    #     'JH':AMPERE_data[:,2],
    #     'elec potenl':AMPERE_data[:,3],
    #     'field-aligned current':AMPERE_data[:,4],
    #     'data rate':AMPERE_dataRate,
    #     }; #make a dict

    return AMPERE_data #return the success
#END DEF

#Rebase to xyz - spherical has made a bb enemy
def geo_to_xyz(lat, lon, alt): #inspired by https://gis.stackexchange.com/a/261230
    lat_rad = lat*np.pi/180.; #rad, convert to rad
    long_rad = lon*np.pi/180.; #rad, conver to rad

    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    f = 1/298.257223563; #unitless, flattening factor
    lat_radSin = np.sin(lat_rad); #reused
    lat_radCos = np.cos(lat_rad); #reused
    e2 = 1 - (1 - f) * (1 - f);
    v = a / np.sqrt(1 - e2 * lat_radSin * lat_radSin);

    x = (v + alt) * lat_radCos * np.cos(long_rad)/1000;
    y = (v + alt) * lat_radCos * np.sin(long_rad)/1000;
    z = (v * (1 - e2) + alt) * lat_radSin/1000;

    return x, y, z
#END DEF

#for parallel work
def regrid_parallel(AMPERE_saveTemp):
    #-----unpack-----
    AMPERE_temp_timeData = AMPERE_saveTemp['AMPERE_temp_timeData'];
    AMPERE_temp_numData = AMPERE_saveTemp['AMPERE_temp_numData'];
    AMPERE_desired_numData = AMPERE_saveTemp['AMPERE_desired_numData'];
    AMPERE_desired_latNum = AMPERE_saveTemp['AMPERE_desired_latNum'];
    AMPERE_desired_longNum = AMPERE_saveTemp['AMPERE_desired_longNum'];
    AMPERE_desiredResStr = AMPERE_saveTemp['AMPERE_desiredResStr'];
    # AMPERE_xOrig = AMPERE_saveTemp['AMPERE_xOrig'];
    # AMPERE_yOrig = AMPERE_saveTemp['AMPERE_yOrig'];
    # AMPERE_zOrig = AMPERE_saveTemp['AMPERE_zOrig'];
    # AMPERE_x = AMPERE_saveTemp['AMPERE_x'];
    # AMPERE_y = AMPERE_saveTemp['AMPERE_y'];
    # AMPERE_z = AMPERE_saveTemp['AMPERE_z'];
    dataKey = AMPERE_saveTemp['dataKey'];
    dataAccuracy = AMPERE_saveTemp['dataAccuracy'];
    
    #precalc what can be precalcd
    mag_string = 'mag'+AMPERE_desiredResStr;
    
    if( AMPERE_saveTemp['mag orig'][dataKey].ndim == 1 ):
        AMPERE_return_geo = np.empty(AMPERE_saveTemp[mag_string]['long'].size,dtype=dataAccuracy); #preallocate
        AMPERE_return_mag = np.empty(AMPERE_saveTemp[mag_string]['long'].size,dtype=dataAccuracy); #preallocate
        # AMPERE_return_mag = AMPERE_saveTemp[mag_string]; #use it for Geo, None it later b/c already have it
    
        for j in range(0,AMPERE_temp_timeData):
            #RegridderZen interpolates OK but does not smooth enough (it is perfectly - and I mean perfectly - exact tho if that's ever important)
            AMPERE_interper_dataOrig = AMPERE_saveTemp['mag orig'][dataKey][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            AMPERE_interper_latOrig = AMPERE_saveTemp['mag orig']['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            AMPERE_interper_longOrig = AMPERE_saveTemp['mag orig']['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
            AMPERE_interper_longOrig[AMPERE_interper_longOrig<0] = AMPERE_interper_longOrig[AMPERE_interper_longOrig<0]+360; #make it roll around, idk if its req'd but doin it
            AMPERE_interper_longOrig_unique = np.unique(AMPERE_interper_longOrig); #since MLT cruises nothing's static but it is always gonna be evenly spaced
            kr = AMPERE_interper_longOrig == AMPERE_interper_longOrig_unique[0]; #get the edge
            AMPERE_interper_longOrig = np.append(AMPERE_interper_longOrig,np.repeat(AMPERE_interper_longOrig_unique[-1]+15,kr.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
            AMPERE_interper_latOrig = np.append(AMPERE_interper_latOrig,AMPERE_interper_latOrig[kr]); #append on lat values (making it cyclical)
            AMPERE_interper_dataOrig = np.append(AMPERE_interper_dataOrig,AMPERE_interper_dataOrig[kr]); #append on data values (making it cyclical)
            kl = AMPERE_interper_longOrig == AMPERE_interper_longOrig_unique[-1]; #get the edge
            AMPERE_interper_longOrig = np.insert(AMPERE_interper_longOrig,0,np.repeat(AMPERE_interper_longOrig_unique[0]-15,kl.sum())); #append on -15 to make extra data for the interper to work with (to make it cyclical in a way)
            AMPERE_interper_latOrig = np.insert(AMPERE_interper_latOrig,0,AMPERE_interper_latOrig[kl]); #append on lat values (making it cyclical)
            AMPERE_interper_dataOrig = np.insert(AMPERE_interper_dataOrig,0,AMPERE_interper_dataOrig[kl]); #append on data values (making it cyclical)
            AMPERE_interper_long = AMPERE_saveTemp[mag_string]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy();
            AMPERE_interper_long[AMPERE_interper_long<0] = AMPERE_interper_long[AMPERE_interper_long<0]+360; #make it roll around too
            AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = \
                interpolate.griddata(np.vstack((AMPERE_interper_longOrig, AMPERE_interper_latOrig)).T, AMPERE_interper_dataOrig, \
                    np.vstack((AMPERE_interper_long, AMPERE_saveTemp[mag_string]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T,\
                    method='linear').reshape(AMPERE_desired_longNum,AMPERE_desired_latNum).ravel(); #one helluva call to interpolate
            # if( np.all(AMPERE_desired_latLongSteps == None) ): #note the below works b/c for default stepping 'mag'+AMPERE_desiredResStr is mapped to 'mag orig' so its the same deal
            #     AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_saveTemp['mag orig'][dataKey][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=0, kernel='linear');
            #     AMPERE_saveTemp[geo_string][dataKey][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
            # else:
            #this special mode is based off of the mag just calc'd and uses the same gridding function as previous
            #Rebase to x y z b/c spherical no go good
            # tic = time.time()
            # AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], smoothing=0, kernel='linear');
            # AMPERE_return_geo[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
            # # #******WARNING****** THIS CODE MAY CRASH YOUR COMPUTER LIKE STRAIGHT SHUT OFF AND TURN BACK ON (BC COMP NOT 100% STABLE) BUT ITS DONE IT ON AN i5-4690K SYSTEM AND FX-8530 SYSTEM (i7-7700K was fine) - def AVX stability (FX just...bad) -******************************************************
            # # #END IF
            # toc = time.time() - tic;
            # print('runo time: '+str(toc))
            
            # tic = time.time()
            # AMPERE_interper_geo2 = RBF_custom.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], smoothing=0, kernel='linear');
            # gg6 = AMPERE_interper_geo2(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
            # # #******WARNING****** THIS CODE MAY CRASH YOUR COMPUTER LIKE STRAIGHT SHUT OFF AND TURN BACK ON (BC COMP NOT 100% STABLE) BUT ITS DONE IT ON AN i5-4690K SYSTEM AND FX-8530 SYSTEM (i7-7700K was fine) - def AVX stability -******************************************************
            # # #END IF
            # toc = time.time() - tic;
            # print('run time: '+str(toc))
            
            #-----attempt to make a faster truly spherical RBF function, it worked (dark, dark magicks abound)-----
            # tic = time.time()
            # weightz = -np.linalg.solve(subfun_rbf_interp.distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
            # weightz = lu_solve(lu_factor(subfun_rbf_interp.distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]),overwrite_a=True, check_finite=False), -AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy(), trans=0, overwrite_b=True, check_finite=False); #fastest with the checks off and overwrites on ayyyy
            # weightz = solve(subfun_rbf_interp.distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), -AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy(), sym_pos=False, lower=False, overwrite_a=True, overwrite_b=True, check_finite=False, assume_a='sym', transposed=False); #slower than above oh well
            # _, _, weightz, _ = dsysv(subfun_rbf_interp.distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), -AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy(), overwrite_a=True, overwrite_b=True); #nope, super slow despite being for symmetric matricies - something up there      
            # weightz = subfun_rbf_interp.rbf_interp_weights(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])
            AMPERE_return_geo[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = rbf_interp(np.float32(AMPERE_saveTemp['AMPERE lat rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), np.float32(AMPERE_saveTemp['AMPERE long rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], \
                                               lu_solve(lu_factor(distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]),overwrite_a=True, check_finite=False), -AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy(), trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
            # toc = time.time() - tic;
            # print('runn time: '+str(toc))
        #END FOR j
        # else:
        #     #-----attempt to make a faster truly spherical RBF function, it worked (dark, dark magicks abound)-----
        #     AMPERE_return_geo[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = rbf_interp(np.float32(AMPERE_saveTemp['AMPERE lat rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), np.float32(AMPERE_saveTemp['AMPERE long rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], \
        #                                            lu_solve(lu_factor(distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]),overwrite_a=True, check_finite=False), -AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy(), trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
        # #END IF
    else:
        AMPERE_return_geo = np.empty( (AMPERE_saveTemp[mag_string]['long'].size,AMPERE_saveTemp['mag orig'][dataKey].shape[1]),dtype=dataAccuracy); #preallocate
        AMPERE_return_mag = np.empty( (AMPERE_saveTemp[mag_string]['long'].size,AMPERE_saveTemp['mag orig'][dataKey].shape[1]),dtype=dataAccuracy); #preallocate
        # AMPERE_return_mag = AMPERE_saveTemp[mag_string]; #use it for Geo, None it later b/c already have it
        
        for kj in range(0,AMPERE_saveTemp['mag orig'][dataKey].shape[1]):
            for j in range(0,AMPERE_temp_timeData):
                #RegridderZen interpolates OK but does not smooth enough (it is perfectly - and I mean perfectly - exact tho if that's ever important)
                AMPERE_interper_dataOrig = AMPERE_saveTemp['mag orig'][dataKey][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData,kj].copy();
                AMPERE_interper_latOrig = AMPERE_saveTemp['mag orig']['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
                AMPERE_interper_longOrig = AMPERE_saveTemp['mag orig']['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
                AMPERE_interper_longOrig[AMPERE_interper_longOrig<0] = AMPERE_interper_longOrig[AMPERE_interper_longOrig<0]+360; #make it roll around, idk if its req'd but doin it
                AMPERE_interper_longOrig_unique = np.unique(AMPERE_interper_longOrig); #since MLT cruises nothing's static but it is always gonna be evenly spaced
                kr = AMPERE_interper_longOrig == AMPERE_interper_longOrig_unique[0]; #get the edge
                AMPERE_interper_longOrig = np.append(AMPERE_interper_longOrig,np.repeat(AMPERE_interper_longOrig_unique[-1]+15,kr.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
                AMPERE_interper_latOrig = np.append(AMPERE_interper_latOrig,AMPERE_interper_latOrig[kr]); #append on lat values (making it cyclical)
                AMPERE_interper_dataOrig = np.append(AMPERE_interper_dataOrig,AMPERE_interper_dataOrig[kr]); #append on data values (making it cyclical)
                kl = AMPERE_interper_longOrig == AMPERE_interper_longOrig_unique[-1]; #get the edge
                AMPERE_interper_longOrig = np.insert(AMPERE_interper_longOrig,0,np.repeat(AMPERE_interper_longOrig_unique[0]-15,kl.sum())); #append on -15 to make extra data for the interper to work with (to make it cyclical in a way)
                AMPERE_interper_latOrig = np.insert(AMPERE_interper_latOrig,0,AMPERE_interper_latOrig[kl]); #append on lat values (making it cyclical)
                AMPERE_interper_dataOrig = np.insert(AMPERE_interper_dataOrig,0,AMPERE_interper_dataOrig[kl]); #append on data values (making it cyclical)
                AMPERE_interper_long = AMPERE_saveTemp[mag_string]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy();
                AMPERE_interper_long[AMPERE_interper_long<0] = AMPERE_interper_long[AMPERE_interper_long<0]+360; #make it roll around too
                AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData,kj] = \
                    interpolate.griddata(np.vstack((AMPERE_interper_longOrig, AMPERE_interper_latOrig)).T, AMPERE_interper_dataOrig, \
                        np.vstack((AMPERE_interper_long, AMPERE_saveTemp[mag_string]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T,\
                        method='linear').reshape(AMPERE_desired_longNum,AMPERE_desired_latNum).ravel(); #one helluva call to interpolate
                
                #-----attempt to make a faster truly spherical RBF function, it worked (dark, dark magicks abound)-----
                AMPERE_return_geo[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData,kj] = rbf_interp(np.float32(AMPERE_saveTemp['AMPERE lat rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), np.float32(AMPERE_saveTemp['AMPERE long rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]), AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], \
                                                   lu_solve(lu_factor(distance_calc_identical(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]),overwrite_a=True, check_finite=False), -AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData, kj].copy(), trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
            #END FOR j
        #END FOR kj
    #END IF
    
    # if( FLG_newMag == False ):
    #     AMPERE_return_mag = None; #nothing to return
    # #END IF
    
    return {'geo':AMPERE_return_geo, 'mag':AMPERE_return_mag}
#END DEF

#not positive definite, no cholesky
# from scipy.linalg import lu_factor, lu_solve #best so far
# from scipy.linalg import solve #nope
# from scipy.linalg.lapack import dsysv #way way slower
# from scipy.linalg import ldl #idk how to use solution (lu has lu_solve for lu_factor)

# import matplotlib.pyplot as plt

# gg_orig = np.flip(AMPERE_saveTemp['mag orig'][dataKey][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].reshape(24,50).astype('float64'),axis=1)
# gg_magUp = np.flip(AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(120,50).astype('float64'),axis=1)
# gg_magUp2 = np.flip(interpolate.griddata(np.vstack((AMPERE_interper_longOrig, AMPERE_interper_latOrig)).T, AMPERE_interper_dataOrig, \
#     np.vstack((AMPERE_interper_long, AMPERE_saveTemp[mag_string]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T,\
#     method='linear').reshape(AMPERE_desired_longNum,AMPERE_desired_latNum).ravel().reshape(120,50).astype('float64'),axis=1)
# gg2 = np.flip(AMPERE_return_geo[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(120,50).astype('float64'),axis=1)
# gg77 = np.flip(gg7.reshape(120,50).astype('float64'),axis=1)

# plt.figure()
# plt.pcolormesh(np.meshgrid(np.arange(0,360+15,15),np.arange(0,50+1,1))[0],np.meshgrid(np.arange(0,360+15,15),np.arange(0,50+1,1))[1],gg_orig.T,shading='flat')
# plt.figure()
# plt.pcolormesh(np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[0],np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[1],gg_magUp.T,shading='flat')
# plt.figure()
# plt.pcolormesh(np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[0],np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[1],gg_magUp2.T,shading='flat')
# plt.figure()
# plt.pcolormesh(np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[0],np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[1],gg2.T,shading='flat')
# plt.figure()
# plt.pcolormesh(np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[0],np.meshgrid(np.arange(0,360+3,3),np.arange(0,50+1,1))[1],gg77.T,shading='flat')

# import Code.subfun_custom_rbfinterp as RBF_custom
# reload(RBF_custom)

# # # from subfun_rbf_interp import rbf_interp, rbf_interp_weights
# import subfun_rbf_interp

# from importlib import reload
# reload(subfun_rbf_interp)



# dist_sphere = np.sin((np.float64(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(1,-1))-np.float64(AMPERE_saveTemp['AMPERE lat orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(-1,1)))*np.pi/360)**2 + np.cos(np.float64(AMPERE_saveTemp['AMPERE lat orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(-1,1))*np.pi/180) * np.cos(np.float64(AMPERE_saveTemp['AMPERE lat orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(1,-1))*np.pi/180) * np.sin((np.float64(AMPERE_saveTemp['AMPERE long orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(1,-1))-np.float64(AMPERE_saveTemp['AMPERE long orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].reshape(-1,1)))*np.pi/360)**2; #modified haversine formula
# haversine2 = -np.arctan2(np.sqrt(dist_sphere), np.sqrt(1 - dist_sphere)) #modified haversine formula



# AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] - AMPERE_saveTemp['AMPERE lat orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]

# np.subtract.outer(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_saveTemp['AMPERE lat orig'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])


# lat1 = np.float64(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
# lat2 = np.float64(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
# long1 = np.float64(AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
# long2 = np.float64(AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);

# haversine3 = np.zeros((lat1.size,lat2.size));
# for jj in range(0,lat1.size-1):
#     #this reduces the effort for identical vector inputs by 1/2 or something
#     dist_sphere = np.sin((lat1[jj]-lat2[jj+1:])*np.pi/360)**2 + np.cos(lat1[jj]*np.pi/180) * np.cos(lat2[jj+1:]*np.pi/180) * np.sin((long1[jj]-long2[jj+1:])*np.pi/360)**2; #modified haversine formula
#     haversine3[jj+1:,jj] =  -np.arctan2(np.sqrt(dist_sphere), np.sqrt(1 - dist_sphere)) #modified haversine formula
# #END FOR jj
# haversine3 += haversine3.T; #add on the other half


# pt_lat = np.float64(AMPERE_saveTemp['AMPERE lat rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
# data_lat = np.float64(AMPERE_saveTemp['AMPERE lat orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
# pt_long = np.float64(AMPERE_saveTemp['AMPERE long rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);
# data_long = np.float64(AMPERE_saveTemp['AMPERE long orig rad'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData]);

# data_out = np.empty(pt_lat.size,dtype=pt_lat.dtype);
# dist_combo = np.zeros((pt_lat.size,pt_lat.size));
# for jj in range(0,pt_lat.size):
#     dist_combo[:,jj] = subfun_rbf_interp.distance_calc(pt_lat[jj], pt_long[jj], data_lat, data_long)
#     data_out[jj] = (-subfun_rbf_interp.distance_calc(pt_lat[jj], pt_long[jj], data_lat, data_long)@weightz).item();
# #END FOR jj

# dist_combo2 = subfun_rbf_interp.distance_calc(pt_lat.reshape(-1,1), pt_long.reshape(-1,1), data_lat.reshape(1,-1), data_long.reshape(1,-1))



