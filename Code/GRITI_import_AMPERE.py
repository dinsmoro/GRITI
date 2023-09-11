#Function to import AMPERE data
#RD on 3/5/2021

#To properly use, place pre-proccessed sources in their respective year folders in the data folder you've specified

#Info on stuff you can send:
#FLG_dataMix = 0; #0 prevents data sources to mixing to fill a time span, 1 allows data sources to mix to fill a time span
#FLG_dataPreference = 0; #preffered data type (by order as appended below at folder_fileNameFormat)

import numpy as np
import os
import h5py
from subfun_strfind import strfind
from subfun_daysInAYear import subfun_daysInAYear
from subfun_dayNum_to_date import subfun_dayNum_to_date

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


def GRITI_import_AMPERE(dates,settings_paths,hemi,settings_loginAMPERE,
                        AMPERE_coordType='geo',FLG_dataMix=0,FLG_dataPreference=0,FLG_float64=0):
    
    version_alg = 1.3; #algorithm version
    #1.0 10/8/2019 - initial algorithm
    #1.1 11/10/2021 - update to deal with attributes when combining previously converted data & data newly being converted, also include chunk shape for faster reading (not big issue now, never know tho)
    #1.2 11/12/2021 - fix longitude distributing error that was causing longtidues to get aliased instead of in their right spots
    #1.3 11/15/2021 - adjust longitude to match old style where -180 is used instead of 180 (makes gridding algs work better)
    
    #----- Unpack -----
    dateRange_full = dates['date range full']; #unpack it
    dateRange_dayNum_full = dates['date range full dayNum']; #unpack it
    # dateRange_dayNum_full = dateRange_dayNum_full_adj; #set it [don't have extra days yet, code in an override so can use this]
    # dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #get it
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum']; #unpack it
    
    if( FLG_float64 == 0 ):
        dataAccuracy = np.float32; #set the data accuracy to 32 bit floats (1/2 the memory)
    else:
        dataAccuracy = np.float64; #set the data accuracy to 64 bit floats
    #END IF
    
    print("AMPERE - Date range requested (yr/day num format): {}/{} to {}/{}.\n".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )

    #==============Constants Needed==============
    folder_AMPERE = 'AMPERE'; #name for the AMPERE folder
    folder_fileNameFormat = []; #prep a list
    # #YR for year
    # #DN for day num
    # #MO for month
    # #DY for day
    #.append more if more arise
    folder_fileNameFormat.append('AMPEREM_#YR_#DN.h5'); #0 - 1st in the list is the compiled version, after that are possible names that need compiling
    folder_fileNameFormat.append('AuroraPHILE 2D Output #YR#MO#DY North Geographic.txt'); #1 - expected file name style
    folder_fileNameFormat.append('adelphi 2D Output #YR#MO#DY North Geographic'); #2 - expected file name style 2nd version
    #folder_fileNameFormat.append('otherFormat_#DY_#MO.h5'); #1 - tack on another (example)
    
    #these are for writing a compact HDF5 file, these are the names that will be used for the dict names
    headerNames = ['lat','Ped','Hall','Phi','E East','E North','J East','J North','E Flux','JH','JR In','JR Out']; #note the nice names [names not here will be tacked on as they are in the header]
    
    #AMPERE_jouleHeating_plotLimValu = [1,5]; #erg/(cm^2*sec), low and high limits for plotting joule heating - taken from example video (seems Joule Heating low end limit is 1 - nothing less than 1)
    #AMPERE_jouleHeating_pos = 2; #set the position of the joule heating data
    AMPERE_data_num = 5; #number of data entries per line (+3 for time/lat/long, added later on)
     #0 Pedersen Conductance / 1 Hall Conductance / 2 Joule Heat (ergs/(cm^2*sec)) / 3 Electric Potential / 4 Field-Algined Current / 5 time (hr) / 6 latitude / 7 longitude
    
                                 
    #==============Prep to look for the data==============                                 
    if( os.path.isdir(settings_paths['data'] + '\\' + folder_AMPERE) == 0 ): #check if AMPERE folder exists
        #if not, make it
        os.makedirs(settings_paths['data'] + '\\' + folder_AMPERE);
        print("NOTA BENE: Importing AMPERE Func - Created AMPERE directory: {}\n".format(settings_paths['data'] + '\\' + folder_AMPERE) );
    #END IF
    
    dateRange_uniqueYears = np.unique(dateRange_dayNum_full[:,0]); #get the unique years involved
    AMPERE_dataAmnt = len(dateRange_dayNum_full[:,1]); #get number of days needed to investigate
    AMPERE_dataAvail_perSource = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for finding where data is available
    AMPERE_dataAvail_toUse = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for recording which piece of data to use on the day
    #0 = no data available on required days, will quit
    #1 = note data is there, ready to use
    
    for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
        if( os.path.isdir(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_uniqueYears[i]) ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_uniqueYears[i]) );
            print("NOTA BENE: Importing AMPERE Func - Created AMPERE subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_uniqueYears[i]) ));
        #END IF
    #END FOR
    
    #==============Look for data in expected naming formats==============
    for i in range(0,len(folder_fileNameFormat)): #loop through the different name formats    
        
        for j in range(0,AMPERE_dataAmnt): #loop through the different years needed
            
            AMPERE_fileName = folder_fileNameFormat[i].replace('#DN', str(dateRange_dayNum_full[j,1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName = AMPERE_fileName.replace('#YR', str(dateRange_dayNum_full[j,0]) ); #replace any #YR with the current year
            AMPERE_fileName = AMPERE_fileName.replace('#MO', str(dateRange_full[j,1]).zfill(2) ); #replace any #MO with the current month (padded w/ 0's so always 2 long)
            AMPERE_fileName = AMPERE_fileName.replace('#DY', str(dateRange_full[j,2]).zfill(2) ); #replace any #DY with current day (padded w/ 0's so always 2 long)
            
            AMPERE_dataAvail_perSource[j,i] = os.path.isfile(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[j,0]) + '\\' + AMPERE_fileName); #record if data is there or not
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
                        # print("\n==============ERROR in GRITI_import_AMPERE==============");
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
                    # print("\n==============ERROR in GRITI_import_AMPERE==============");
                    # print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                    # #return("No"); #return something that will def crash things
                    # import sys #yolo import
                    # sys.exit(); #more def will crash
                    AMPERE_dataAvail_toUse[:,np.where(AMPERE_dataAvail_entireSpan)[0][0]] = 1; #use the closest to 0
                #END IF
            #END IF
                
            # AMPERE_dataAvail_toUse[:,AMPERE_dataAvail_entireSpan] = 1; #set 
            
            FLG_dataAvail_entireSpan = 1; #set the flag to good
        #END IF
    #END IF
        
    if( FLG_dataAvail_entireSpan == 0 ): #make sure there's data
        print("\n==============ERROR in GRITI_import_AMPERE==============");
        print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
        print("Printing file name formats supported:");
        print("{}".format(folder_fileNameFormat)); #print for error
        print("Printing available data matrix (made of dates and file name formats):");
        print("{}".format(AMPERE_dataAvail_perSource)); #print for error - lets user know available days
        print("Will exit via returning no");
        #return("No"); #return something that will def crash things
        import sys #yolo import
        sys.exit(); #more def will crash
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
    AMPERE_data = {}; #prep data dict to hold the data
    for i in range(0,AMPERE_dataAmnt ):
        FLG_recalc = False; #set recalc flag to false every time
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
        #prep the name needed
        
        if( AMPERE_sourceIndex == 0 ): #this means it is the first type that has been refactored into an HDF5 file
            #load saved AMPERE to speed up
            keyz = list(AMPERE_data.keys()); #get the current keys in AMPERE_data
            with h5py.File(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, 'r') as AMPERE_file:
                if( AMPERE_file.attrs['version'] >= version_alg ):
                    #--- Read the data types in ---
                    keyzNew = list(AMPERE_file.keys()); #get the saved keys
                    for j in range(0,len(keyzNew)):
                        if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                            AMPERE_data[keyzNew[j]].append(AMPERE_file.get(keyzNew[j])[()]); #tack that dataset on
                        else:
                            #otherwise it's a new data type to add in
                            AMPERE_data[keyzNew[j]] = [AMPERE_file.get(keyzNew[j])[()]]; #get that dataset out
                        #END IF
                    #END FOR j
                    #--- Read the attributes in ---
                    keyzNew = list(AMPERE_file.attrs.keys()); #get the attribute keys
                    keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                    for j in range(0,len(keyzNew)):
                        if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                            if( np.isclose(AMPERE_file.attrs[keyzNew[j]], AMPERE_data[keyzNew[j]]) == False ):
                                #only worry is if the attribute isn't consistent
                                print('-----Warning-----');
                                print('Attribute '+keyzNew[j]+' isn\'t the same as the previously recorded value from another file of '+ \
                                    str(AMPERE_data[keyzNew[j]])+' and this file\'s value of '+str(AMPERE_file.attrs[keyzNew[j]])+ \
                                    '.\n NaN\'ing it and try to sort that out.');
                                AMPERE_data[keyzNew[j]] = np.nan; #nan that attribute, figure it out later
                            #END IF
                        else:
                            AMPERE_data[keyzNew[j]] = AMPERE_file.attrs[keyzNew[j]]; #get that attribute out
                        #END IF
                    #END FOR j
                else:
                    FLG_recalc = True; #set flag to true, recalc is needed
                    oldVersion = AMPERE_file.attrs['version'];
                #END IF
            #END WITH
        #END IF
        if( FLG_recalc == True ):
            print('---WARNING in GRITI_import_AMPERE---');
            print(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+ \
                ' file version is '+str(oldVersion)+' which is less than the current file version of '+str(version_alg)+'!'+\
                '\nThat\'s bad! Renaming file to add a \'_oldV'+str(oldVersion).replace('.','p')+'\' at the end, then gonna calculate a new version of the file!');
            os.rename(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, \
                settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+'_oldV'+str(oldVersion).replace('.','p')); #rename
            
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
                    print("\n==============ERROR in GRITI_import_AMPERE==============");
                    print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
                    print("Printing file name formats supported:");
                    print("{}".format(folder_fileNameFormat)); #print for error
                    print("Printing available data matrix (made of dates and file name formats):");
                    print("{}".format(AMPERE_dataAvail_perSource)); #print for error - lets user know available days
                    print("Will exit via returning no");
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.exit(); #more def will crash
                #END IF
            else:
                if( np.sum(AMPERE_dataAvail_perSource[i,1:]) > 0 ): #if the sum is greater than 0, choose one
                    if( (AMPERE_dataAvail_perSource[i,FLG_dataPreference] == 1) & (FLG_dataPreference != 0) ): #easy, set it
                        AMPERE_dataAvail_toUse_tmp[i,FLG_dataPreference] = 1; #set it to 1
                    else:
                        AMPERE_dataAvail_toUse_tmp[i,np.where(AMPERE_dataAvail_perSource[i,1:])[0][0]+1] = 1; #use the closest to 0
                    #END IF
                else:
                    print("\n==============ERROR in GRITI_import_AMPERE==============");
                    print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
                    print("Printing file name formats supported:");
                    print("{}".format(folder_fileNameFormat)); #print for error
                    print("Printing available data matrix (made of dates and file name formats):");
                    print("{}".format(AMPERE_dataAvail_perSource)); #print for error - lets user know available days
                    print("Will exit via returning no");
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.exit(); #more def will crash
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
        #END IF
                    
        if( AMPERE_sourceIndex >= 1): #this means that its the second original format that needs to be parsed
    
            #else read AMPERE files to get the data needed
            print('Some or none of "cached" pre-read AMPERE files found; beginning to read AMPERE data and save as HDF5 for future use - est. 1 min per ('+str(np.round(1*(AMPERE_dataAmnt-np.sum(np.int64(AMPERE_dataAvail_perSource[:,0]))),2)).strip('0').strip('.')+' min total) on fast comp');
            import time
            from subfun_strstr import strstr
            import re
            tic = time.time(); #start timing
            
            #--- Read in the raw AMPERE file ---
            with open(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, 'r') as AMPERE_file: #open file pointer (with for safety cause life is hard apparently)
                #AMPERE_raws = textscan(AMPERE_file, '%s','delimiter','\n'); #multiple headers of varying sizes make reading this real annoying
                AMPERE_raws = AMPERE_file.readlines(); #multiple headers of varying sizes make reading this real annoying
            #END with
            AMPERE_raws_len = len(AMPERE_raws); #get the length for limit checks
            
            #--- Look for the start line (in case it isn't 0 in future iterations) ---
            startLine = np.array([]); #prep
            cntr = 0; #cntr for the while
            while( (startLine.size == 0) & (cntr != AMPERE_raws_len) ): #go till it is there
                startLine = strstr(AMPERE_raws[cntr],str(dateRange_dayNum_full[i,0])+str(dateRange_full[i,1]).zfill(2)+str(dateRange_full[i,2]).zfill(2)); #go from the top down looking for the start line
                cntr += 1; #increment
            #END WHILE
            if( cntr == AMPERE_raws_len ):
                print("\n==============ERROR in GRITI_import_AMPERE==============");
                print('File identified as the date '+str(dateRange_dayNum_full[i,0])+str(dateRange_full[i,1]).zfill(2)+str(dateRange_full[i,2]).zfill(2)+
                      ' does not contain that date!\nDelete the file (and preferably replace it with a correct one): '+
                      settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+'\nGonna crash now.');
                #return("No"); #return something that will def crash things
                import sys #yolo import
                sys.crash(); #more def will crash
            #END IF
            startLine = startLine.item(); #get the start line number
            
            #--- Analyze the file for data size and other "header" values ---
            # timeList_yr = []; #list of the times
            # timeList_dn = []; #list of the times
            timeList_loc = []; #list of the time locations
            timeList_hr = []; #list of the times
            timeList_min = []; #list of the times
            timeList_sec = []; #list of the times
            longList = []; #list of the longitudes
            dataList = []; #list of where the data is
            headerList = []; #list of the header info
            for j in range(startLine,AMPERE_raws_len):
                if( re.search('[a-zA-Z]', AMPERE_raws[j].strip('\n')) == None ):
                    #it's a data line
                    dataList.append(j); #add on the index
                else:
                    if( strstr(AMPERE_raws[j],str(dateRange_dayNum_full[i,0])+str(dateRange_full[i,1]).zfill(2)+str(dateRange_full[i,2]).zfill(2)).size != 0 ):
                        #it's a time line
                        timeTemp = float(AMPERE_raws[j][8:strstr(AMPERE_raws[j],' UT')[0]].strip(' ')); #hr, get only the time part of the string, if this code is used until years aren't 4 digits I'll be very disappointed in everyone who didn't make better code between now and then
                        # timeList_yr.append(np.int16(dateRange_dayNum_full[i,0])); #record
                        # timeList_dn.append(np.int16(dateRange_dayNum_full[i,1])); #record
                        timeTemp_hr = np.uint8(timeTemp); #hr, get the hours
                        timeTemp_min = (timeTemp - timeTemp_hr)*60; #min, get the minutes (still float)
                        if( np.isclose( timeTemp_min, np.round(timeTemp_min) , rtol=1.e-2) ):
                            timeTemp_min = np.round(timeTemp_min); #round it, deals with the fractional float
                        #END IF
                        timeTemp_min = np.uint8(timeTemp_min); #min, get the minutes
                        timeTemp_sec = np.uint8((timeTemp_min - (timeTemp - timeTemp_hr)*60)*60); #sec, get the seconds [no microseconds so this is where it ends]
                        timeList_hr.append(timeTemp_hr); #record
                        timeList_min.append(timeTemp_min); #record
                        timeList_sec.append(timeTemp_sec); #record
                        timeList_loc.append(j); #record the index
                    elif( strstr(AMPERE_raws[j].lower(),'longitude').size != 0 ):
                        #it's a longitude line
                        if( strstr(AMPERE_raws[j].lower(), 'east').size != 0 ):
                            #it's deg east units
                           longTemp = float(AMPERE_raws[j][strstr(AMPERE_raws[j].lower(), '=')[0]+1:strstr(AMPERE_raws[j].lower(), 'degrees')[0]].strip(' ')); #get the longitude
                           if( longTemp >= 180 ):
                               longTemp -= 360; #deg, keep it from going 0 to 360 to be -180 to 180
                           #END IF
                           longList.append(dataAccuracy(longTemp)); #append
                        else:
                            #it used west
                            longTemp = float(AMPERE_raws[j][strstr(AMPERE_raws[j].lower(), '=')[0]+1:strstr(AMPERE_raws[j].lower(), 'degrees')[0]].strip(' ')); #get the longitude
                            if( longTemp >= 180 ):
                                longTemp -= 360; #deg, keep it from going 0 to 360 to be -180 to 180
                            #END IF
                            longList.append(dataAccuracy(-longTemp)); #append [west is just negative of east I think, I didn't think hard this prob won't happen]
                        #END IF
                    else:
                        #it's a header line
                        if( len(headerList) == 0 ):
                            headerList = AMPERE_raws[j].split(); #fill in the header list
                        #END IF
                    #END IF
                #END IF
            #END FOR j
            # timeList_yr = np.asarray(timeList_yr); #convert to numpy array
            # timeList_dn = np.asarray(timeList_dn); #convert to numpy array
            timeList_loc.append(AMPERE_raws_len); #tack on max value to make later loops happy
            timeList_loc = np.asarray(timeList_loc); #covert to numpy array
            timeList_hr = np.asarray(timeList_hr); #convert to numpy array
            timeList_min = np.asarray(timeList_min); #convert to numpy array
            timeList_sec = np.asarray(timeList_sec); #convert to numpy array
            longList = np.asarray(longList); #convert to numpy array
            dataList = np.asarray(dataList); #convert to numpy array
            
            #--- Calculate data rate ---
            AMPERE_dataRate = np.median(np.diff(np.uint32(timeList_hr)*3600+np.uint32(timeList_min)*60+np.uint32(timeList_sec))); #sec, calculate the data rate
            
            #--- Preallocate time & long arrays ---
            AMPERE_save = {
                'year':np.ones(dataList.size,dtype=np.int16)*np.int16(dateRange_dayNum_full[i,0]), #yr, record the year b/c it won't change here
                'dayNum':np.ones(dataList.size,dtype=np.int16)*np.int16(dateRange_dayNum_full[i,1]), #dayNum, record the dayNum b/c it won't change here
                'hour':np.empty(dataList.size,dtype=np.uint8), #preallocate
                'min':np.empty(dataList.size,dtype=np.uint8), #preallocate
                'sec':np.empty(dataList.size,dtype=np.uint8), #preallocate
                'long':np.empty(dataList.size,dtype=dataAccuracy), #preallocate
                'data rate':AMPERE_dataRate, #sec, record data rate
                };
            
            #--- Preallocate data arrays ---
            headerLocs = -np.ones(len(headerNames),np.int64); #preallocate 1st is headerNames loc, 2nd is headerList loc
            headerList_aligned = [zz.lower().replace(' ','') for zz in headerList]; #had to use list loop oof
            for j in range(0,len(headerNames)):
                headerFound = strfind(headerList_aligned , headerNames[j].lower().replace(' ','') ); #find match
                if( np.sum(headerFound) > 0 ): #deal with unfound data types
                    headerLocs[j] = np.where(headerFound == 1)[0].item(); #get the location of the match
                else:
                    AMPERE_save[headerNames[j]] = np.ones(2,dtype=dataAccuracy)*np.nan; #record that the data isn't there [size 2 to avoid the attribute detecting heuristic]
                #END IF
            #END FOR j
            k = np.where(headerLocs < 0)[0]; #get where there are no data type matches [already recorded as NaN]
            headerLocs = np.delete(headerLocs,k); #delete location entries that are -1 with no matches [already recorded as NaN]
            headerNames = np.delete(np.asarray(headerNames),k).tolist(); #delete from headerName list as well, numpy flexes on lists
            if( len(headerList) > len(headerNames) ): #check for extra data types
                missingIndexes = np.setdiff1d(np.arange(0,len(headerList),step=1,dtype=np.int64), headerLocs, assume_unique=True); #get the missing indexes
                for j in range(0,missingIndexes.size):
                    headerLocs = np.append(headerLocs,missingIndexes[j]); #tack on the location of the missing data type
                    headerNames.append(headerList[missingIndexes[j]]); #tack on the missing data type
                #END FOR j
            #END IF
            for j in range(0,len(headerNames)): #preallocate all of the data types
                AMPERE_save[headerNames[j]] = np.empty(dataList.size,dtype=dataAccuracy); #preallocate
            #END FOR j
            
            #--- Record data into those arrays ---
            dataList_deltas = np.append(np.insert(np.where(np.diff(dataList)>1)[0],0,0),dataList.size-1); #get the line skips where there are the headers occur
            cntr_time = 0; #cntr for time
            cntr_long = 0; #cntr for longitude
            for j in range(0,dataList.size):
                indxr = dataList[j]; #get the current index
                # np.sum(indxr > timeList_loc)
                if( indxr > timeList_loc[cntr_time+1] ): #check to increment time cntr
                    cntr_time += 1; #increment
                #END IF
                if( indxr > dataList[dataList_deltas[cntr_long+1]] ): #check to increment long cntr
                    cntr_long += 1; #increment
                #END IF
                AMPERE_save['hour'][j] = timeList_hr[cntr_time]; #record
                AMPERE_save['min'][j] = timeList_min[cntr_time]; #record
                AMPERE_save['sec'][j] = timeList_sec[cntr_time]; #record
                AMPERE_save['long'][j] = longList[cntr_long]; #record
                #now read the data line itself
                dataLine = list(filter(None,AMPERE_raws[indxr].strip().split(sep = ' '))); #get the line, strip off spaces and \n, split by spaces and remove the empty strings from multiple spaces next to eachother
                #scan the line for overlapping values
                jk = np.where( strfind(dataLine,'.') > 1 )[0]; #get where there are more than 1 .'s (indicative of a line smooshed together)
                if( jk.size > 0 ):
                    # print(dataLine); #report
                    for k in range(jk.size-1,0-1,-1): #count backwards so inserting indexes doesn't mess up earlier indexes
                        splitr = strstr(dataLine[jk[k]],'-')[-1]; #split on a -, otherwise no hope to guess
                        line1 = dataLine[jk[k]][0:splitr]; #get the 1st line
                        line2 = dataLine[jk[k]][splitr:]; #get the 2nd line
                        dataLine[jk[k]] = line1; #set the 1st line to its right value
                        dataLine.insert(jk[k]+1,line2); #insert the 2nd line
                    #END FOR k
                #END IF
                #now actually record the data
                for k in range(0,headerLocs.size):
                    AMPERE_save[headerNames[k]][j] = dataAccuracy(dataLine[headerLocs[k]]); #record data with selected data accuracy
                #END FOR k
            #END FOR j            
            
            #--- Save the data as HDF5 ---
            #0 is the hdf5 format goal
            AMPERE_fileName_write = folder_fileNameFormat[ 0 ].replace('#DN', str(dateRange_dayNum_full[i,1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#YR', str(dateRange_dayNum_full[i,0]) ); #replace any #YR with the current year
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#MO', str(dateRange_full[i,1]).zfill(2) ); #replace any #MO with the current month
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#DY', str(dateRange_full[i,2]).zfill(2) ); #replace any #DY with current day
            
            #now save the data to use again
            keyz = list(AMPERE_save.keys()); #keys to the dict
            h5pyChunkShape = AMPERE_save['hour'].shape; #get the shape of one of the vectors and use it as a chunk (read only whole chunks)
            with h5py.File(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName_write, 'w') as AMPERE_file:
                for j in range(0,len(keyz)):
                    if( AMPERE_save[keyz[j]].size > 1 ):
                        AMPERE_file.create_dataset(keyz[j], data=AMPERE_save[keyz[j]], chunks=h5pyChunkShape, compression="gzip"); #write that data
                    else:
                        #if size 1, add it as an attribute
                        AMPERE_file.attrs[keyz[j]] = AMPERE_save[keyz[j]]; #save the attribute
                    #END IF
                #END FOR j
                #add on non-data-related attributes
                AMPERE_file.attrs['version'] = version_alg; #save the attribute
            #END WITH
            
            #--- Copy to main data var as well (that gets returned) ---
            if( i == 0 ):
                for j in range(0,len(keyz)):
                    AMPERE_data[keyz[j]] = [AMPERE_save[keyz[j]]]; #record
                #END FOR j
            else:
                for j in range(0,len(keyz)):
                    if( np.sum(strfind(list(AMPERE_data.keys()),keyz[j])) > 0 ): #make sure key exists
                        if( AMPERE_save[keyz[j]].size > 1 ):
                            AMPERE_data[keyz[j]].append(AMPERE_save[keyz[j]]); #tack on
                        else:
                            if( np.isclose(AMPERE_save[keyz[j]], AMPERE_data[keyz[j]]) == False ):
                                #only worry is if the attribute isn't consistent
                                print('-----Warning-----');
                                print('Attribute '+keyz[j]+' isn\'t the same as the previously recorded value from another file of '+ \
                                    str(AMPERE_data[keyz[j]])+' and this file\'s value of '+str(AMPERE_save[keyz[j]])+ \
                                    '.\n NaN\'ing it and try to sort that out.');
                                AMPERE_data[keyz[j]] = np.nan; #nan that attribute, figure it out later
                            #END IF
                        #END IF
                    else:
                        if( AMPERE_save[keyz[j]].size > 1 ):
                            AMPERE_data[keyz[j]] = [AMPERE_save[keyz[j]]]; #record if doesn't exist
                        else:
                            if( np.isclose(AMPERE_save[keyz[j]], AMPERE_data[keyz[j]]) == False ):
                                #only worry is if the attribute isn't consistent
                                print('-----Warning-----');
                                print('Attribute '+keyz[j]+' isn\'t the same as the previously recorded value from another file of '+ \
                                    str(AMPERE_data[keyz[j]])+' and this file\'s value of '+str(AMPERE_save[keyz[j]])+ \
                                    '.\n NaN\'ing it and try to sort that out.');
                                AMPERE_data[keyz[j]] = np.nan; #nan that attribute, figure it out later
                            #END IF
                        #END IF
                    #END IF
                #END FOR j
            #END IF
    
            #clear AMPERE_data_par AMPERE_data_temp
            toc = time.time() - tic; #end timing
            print("AMPERE data parsing and re-saving as HDF5 took: "+str(np.round(toc,2))+" sec / "+str(np.round(toc/60,2))+" min");
        
        #END IF
                
    #END FOR i
    #--- Convert to numpy arrays from the dynamically added lists ---
    keyz = list(AMPERE_data.keys()); #get the current keys
    for j in range(0,len(keyz)):
        if( np.isscalar(AMPERE_data[keyz[j]]) == False ):
            #if not a scalar, apply the logical mask
            AMPERE_data[keyz[j]] = np.hstack(AMPERE_data[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
        #END IF
    #END FOR j
    
    # #Clear out AMPERE data less than the minimum
    # k = np.where(np.min(AMPERE_jouleHeating_plotLimValu) > AMPERE_data[:,AMPERE_jouleHeating_pos])[0]; #find entries less than the min plotting number (clear it up)
    # AMPERE_data = np.delete(AMPERE_data,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    # #AMPERE_time = np.delete(AMPERE_time,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    # #AMPERE_lat = np.delete(AMPERE_lat,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    # #AMPERE_long = np.delete(AMPERE_long,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    
    # Calc total sec for the entire time period aligned to the zero hr [debuting year support! I hope it works]
    AMPERE_data['time aligned'] = np.int32(AMPERE_data['dayNum']-dateRange_dayNum_zeroHr[1]-(dateRange_dayNum_zeroHr[0]-AMPERE_data['year'])*subfun_daysInAYear(AMPERE_data['year']))*86400 + \
        np.int32(AMPERE_data['hour'])*3600 + np.int32(AMPERE_data['min'])*60 + np.int32(AMPERE_data['sec']); #sec, calc total sec for the day range
    AMPERE_data['time'] = np.int32(AMPERE_data['dayNum'])*86400 + \
        np.int32(AMPERE_data['hour'])*3600 + np.int32(AMPERE_data['min'])*60 + np.int32(AMPERE_data['sec']); #sec, calc total sec for the day range
    
    #Alias JR in to field-aligned current
    AMPERE_data['field-aligned current'] = AMPERE_data['JR In']; #it's the same
    
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