#Function to import AMPERE data
#RD on 3/19/2019
#info on how files were created:
#        for zt=0,239 do begin
#         printf,5,utout(zt)
#         for zlat=0,49 do begin
#           printf,5,clat(zlat)
#           for zlong=0,119 do begin
#             printf,5,qsigp(zt,zlong,zlat),qsigh(zt,zlong,zlat),qjh(zt,zlong,zlat),qpot(zt,zlong,zlat),qjp(zt,zlong,zlat)
#           endfor
#         endfor
#       endear
# In other words, each block of data starts with the UT, then a geographic latitude beginning with 40 degrees, and then a block of 120 lines, each containing six quantities:
# Pedersen Conductance, Hall Conductance, Joule Heat in ergs/cm2-sec, electric potential, field-aligned current
# The 120 lines correspond to geographic longitudes at 3-degree intervals beginning with 0 degrees.

#-----AMPERE File Layout----- (also data imported layout)
    #Float 32 layout
    #0 = Pedersen Conductance [?]
    #1 = Hall Conductance [?]
    #2 = Joule Heat [ergs/(cm^2*sec)]
    #3 = Electric Potential [?]
    #4 = Field-Algined Current [?]
    #5 = time [days] - does not support years 
    #6 = latitude [arcdeg]
    #7 = longitude [arcdeg]

#To properly use, place pre-proccessed sources in their respective year folders in the data folder you've specified

#Info on stuff you can send:
#FLG_dataMix = 0; #0 prevents data sources to mixing to fill a time span, 1 allows data sources to mix to fill a time span
#FLG_dataPreference = 0; #preffered data type (by order as appended below at folder_fileNameFormat)

import numpy as np
import os
import h5py
import aacgmv2 #install with: pip install aacgmv2 [need buildtools what a pain]
import gc
from multiprocessing import cpu_count #only for cpu_count
import joblib #lets multiprocess happen w/o an insane reloading of GRITI_main
from scipy import interpolate
from subfun_strfind import strfind
from subfun_dayNum_to_date import subfun_dayNum_to_date
from subfun_daysInAYear import subfun_daysInAYear
from subfun_textNice import textNice

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


def GRITI_import_AMPERE_preprepared(dates,settings_paths,hemi,settings_loginAMPERE,
                                    AMPERE_coordType='geo',FLG_dataMix=0,FLG_dataPreference=0,FLG_float64=0):
    
    version_alg = 1.0; #algorithm version
    #1 10/8/2019 - initial algorithm
    
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
    folder_fileNameFormat.append('AMPEREO_#YR_#DN.h5'); #0 - tack on another (example)
    folder_fileNameFormat.append('AMPERE_model_data_out_#YR#MO#DY'); #1 - expected file name style
    #folder_fileNameFormat.append('otherFormat_#DY_#MO.h5'); #1 - tack on another (example)
    tEst_mag = 2222/3*dateRange_dayNum_full.shape[0]; #sec, initial estimate for converison to mag
    
    #AMPERE_jouleHeating_plotLimValu = [1,5]; #erg/(cm^2*sec), low and high limits for plotting joule heating - taken from example video (seems Joule Heating low end limit is 1 - nothing less than 1)
    #AMPERE_jouleHeating_pos = 2; #set the position of the joule heating data
    # AMPERE_data_num = 5; #number of data entries per line (+3 for time/lat/long, added later on)
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
                        AMPERE_dataAvail_toUse[i,:] = AMPERE_dataAvail_toUse[i,:] & 0; #0 it out
                        AMPERE_dataAvail_toUse[i,FLG_dataPreference] = 1; #set it to 1
                    else:
                        print("\n==============ERROR==============");
                        print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                        #return("No"); #return something that will def crash things
                        import sys #yolo import
                        sys.crash(); #more def will crash
                        #don't feel like writing a disp, but def I didn't get around to writing alternate stuff
                    #END IF
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
                    AMPERE_dataAvail_entireSpan = AMPERE_dataAvail_entireSpan & 0; #0 it out
                    AMPERE_dataAvail_entireSpan[FLG_dataPreference]  = 1; #put a 1 where the preffered data is there
                else: 
                    #don't feel like writing a disp, but def I didn't get around to writing alternate stuff
                    print("\n==============ERROR==============");
                    print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.exit(); #more def will crash
                #END IF
            #END IF
                
            AMPERE_dataAvail_toUse[:,AMPERE_dataAvail_entireSpan] = 1; #set 
            
            FLG_dataAvail_entireSpan = 1; #set the flag to good
        #END IF
        
    #END IF
        
    if( FLG_dataAvail_entireSpan == 0 ): #make sure there's data
        print("\n==============ERROR==============");
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
    
    #==============Prep AMPERE-specific variables since file format is odd==============
    AMPERE_data = {}; #prep data dict to hold the data
    
    AMPERE_dataRate = 360; #sec, set the time step to split time into
    AMPERE_time_rangeSingle = np.arange(0,86400,AMPERE_dataRate); #get the number of steps needed
    AMPERE_time_rangeSingleNum = AMPERE_time_rangeSingle.size; #get the number of steps needed
    AMPERE_time_range = np.arange(dates['date range zero hr hour bounds'][0]*3600,dates['date range zero hr hour bounds'][1]*3600,AMPERE_dataRate); #sec, AMPERE data is every 6 min
    # AMPERE_time_range = AMPERE_time_range/24 - (np.where( np.where(dateRange_dayNum_zeroHr[0] == dateRange_dayNum_full[:,0])[0] == np.where(dateRange_dayNum_zeroHr[1] == dateRange_dayNum_full[:,1])[0] )[0][0]) + dateRange_dayNum_zeroHr[1]; #days, AMPERE data in day format
    #AMPERE_time_range = AMPERE_time_range - (np.where( np.where(dateRange_dayNum_zeroHr[0] == dateRange_dayNum_full[:,0])[0] == np.where(dateRange_dayNum_zeroHr[1] == dateRange_dayNum_full[:,1])[0] )[0][0])*24; #hr, AMPERE data 0 hr on the zero hour date given above (e.g. May 7th, 2013)
    AMPERE_lat_range = dataAccuracy(np.arange(40,90,1)); #degc, lat, ampere's lat range
    AMPERE_long_range = dataAccuracy(np.roll(np.arange(0,360,3) - 180, np.int64(np.round(np.arange(0,360,3).size/2)) )); #degc, long, ampere's long range (starts at 0 long, increases, wraps around to -180, etc)
    AMPERE_time_len = AMPERE_time_range.size; #length of time range (how many unique times there are)
    AMPERE_lat_len = AMPERE_lat_range.size; #length of lat range (how many unique lats there are)
    AMPERE_long_len = AMPERE_long_range.size; #length of long range (how many unique longs there are)
    AMPERE_data_len = AMPERE_time_len*AMPERE_lat_len*AMPERE_long_len; #length of data vector (time/lat/long repeated in varying ways to match data)
    AMPERE_time = np.zeros(AMPERE_data_len,dtype=np.int32); #sec, time that corresponds to data points (50*120=6,000 hr -24 then 6,000 hr -23.9 etc.)
    for i in range(0,len(AMPERE_time_range) ): #loop to make time data
        AMPERE_time[ i*AMPERE_lat_len*AMPERE_long_len:(i+1)*AMPERE_lat_len*AMPERE_long_len ] = np.tile(AMPERE_time_range[i],AMPERE_lat_len*AMPERE_long_len); #hr, replicate each AMPERE_time_range entry 6,000 times (50lat*120long) to match data cadence
    #END FOR i
    #AMPERE_lat_local = np.tile(AMPERE_lat_range,AMPERE_time_len); #degc, AMPERE's lat range that matches file readout for checking
    AMPERE_lat = dataAccuracy(np.zeros(AMPERE_lat_len*AMPERE_long_len)); #degc, lat that corresponds to data points prep (40 120 times then 41 120 times etc.)
    for i in range(0,len(AMPERE_lat_range) ): #loop to make lat data
       	AMPERE_lat[i*AMPERE_long_len:(i+1)*AMPERE_long_len] = np.tile(AMPERE_lat_range[i],AMPERE_long_len); #degc, put in (40 120 times then 41 120 times etc.)
    #END FOR i
    AMPERE_lat = np.tile(AMPERE_lat,AMPERE_time_len); #degc, replicate to cover all of the time steps
    AMPERE_long = np.tile(AMPERE_long_range,AMPERE_time_len*AMPERE_lat_len); #degc, long that corresponds to data points
    
    #key AMPERE data descrip here           
    AMPERE_data['year'] = np.sort(np.tile(np.int16(dateRange_dayNum_full[:,0]),AMPERE_time_rangeSingleNum*AMPERE_lat_len*AMPERE_long_len)); #yr, record the year b/c it won't change here
    AMPERE_data['dayNum'] = np.sort(np.tile(np.int16(dateRange_dayNum_full[:,1]),AMPERE_time_rangeSingleNum*AMPERE_lat_len*AMPERE_long_len)); #dayNum, record the dayNum b/c it won't change here
    AMPERE_data['data rate'] = AMPERE_dataRate; #sec, record data rate
    AMPERE_data['time'] = AMPERE_time; #sec, record
    AMPERE_data['lat'] = AMPERE_lat; #degc, latitude
    AMPERE_data['long'] = AMPERE_long; #degc, longitude
    # AMPERE_data['Ped'] = np.zeros(AMPERE_data_len,dtype=dataAccuracy); #0 Pedersen Conductance 
    # AMPERE_data['Hall'] = np.zeros(AMPERE_data_len,dtype=dataAccuracy); #/ 1 Hall Conductance
    # AMPERE_data['JH'] = np.zeros(AMPERE_data_len,dtype=dataAccuracy); #2 Joule Heat (ergs/(cm^2*sec))
    # AMPERE_data['elec poten'] = np.zeros(AMPERE_data_len,dtype=dataAccuracy); #3 Electric Potential
    # AMPERE_data['field-aligned current'] = np.zeros(AMPERE_data_len,dtype=dataAccuracy); #4 Field-Algined Current

    
    
    #==============Read in the data==============
    for i in range(0,AMPERE_dataAmnt ):
        
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
                if( AMPERE_file.attrs['version'] <= version_alg ):
                    #--- Read the data types in ---
                    keyzNew = list(AMPERE_file.keys()); #get the saved keys
                    for j in range(0,len(keyzNew)):
                        if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                            AMPERE_data[keyzNew[j]] = np.hstack( (AMPERE_data[keyzNew[j]],AMPERE_file.get(keyzNew[j])[()]) ); #tack that dataset on
                        else:
                            #otherwise it's a new data type to add in
                            AMPERE_data[keyzNew[j]] = AMPERE_file.get(keyzNew[j])[()]; #get that dataset out
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
                    print('\n==============ERROR==============');
                    print(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+ \
                        ' file version is '+str(AMPERE_file.attrs['version'])+' which is less than the current file version of '+str(version_alg)+'!'+\
                        '\nThat\'s bad! Renaming file to add a \'_old\' at the end, then gonna crash. Re-run and code will remake the file with the latest alg.\nI didn\'t code in a way to recover from this. Srry </3.');
                    os.rename(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, \
                        settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+'_old'); #rename
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.crash(); #more def will crash
                #END IF
            #END WITH                    
        elif( AMPERE_sourceIndex == 1): #this means that its the second original format that needs to be parsed
    
            #else read AMPERE files to get the data needed
            print('Some or none of "cached" pre-read AMPERE files found; beginning to read AMPERE data and save as HDF5 for future use - est. 1.5 min per ('+str(np.round(1.5*(AMPERE_dataAmnt-np.sum(np.int64(AMPERE_dataAvail_perSource[:,0]))),2)).strip('0').strip('.')+' min total) on fast comp');
            import time
            tic = time.time(); #start timing
            #create local time variable
            AMPERE_time_rangeSingle = np.arange(0,86400,AMPERE_dataRate); #get the number of steps needed
            AMPERE_time_rangeSingleNum = AMPERE_time_rangeSingle.size; #get the number of steps needed
            AMPERE_save_dataLen = AMPERE_time_rangeSingleNum*AMPERE_lat_len*AMPERE_long_len; #get the total data len for a single day
            
            AMPERE_save = {
                'Ped':np.zeros(AMPERE_save_dataLen,dtype=dataAccuracy), #0 Pedersen Conductance 
                'Hall':np.zeros(AMPERE_save_dataLen,dtype=dataAccuracy), #/ 1 Hall Conductance
                'JH': np.zeros(AMPERE_save_dataLen,dtype=dataAccuracy), #2 Joule Heat (ergs/(cm^2*sec))
                'elec poten':np.zeros(AMPERE_save_dataLen,dtype=dataAccuracy), #3 Electric Potential
                'field-aligned current':np.zeros(AMPERE_save_dataLen,dtype=dataAccuracy), #4 Field-Algined Current
                };
            
            #Import raw AMPERE data
            with open(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, 'r') as AMPERE_file: #open file pointer (with for safety cause life is hard apparently)
                #AMPERE_raws = textscan(AMPERE_file, '%s','delimiter','\n'); #multiple headers of varying sizes make reading this real annoying
                AMPERE_raws = AMPERE_file.readlines(); #multiple headers of varying sizes make reading this real annoying
            #END with
            #AMPERE_raws = AMPERE_raws{1,1}; #pull out the data to make it easier to work with
    
            cntr_time = -1; #prep counter for time (240*# of files)
            cntr_lat = AMPERE_lat_len+1; #prep counter for lats (50 per time) - starts off 50+1 for logic
            cntr_long = 0; #prep counter for longs (120 per lat)
            for j in range(0,len(AMPERE_raws)): #run through each line
                #python is weak and the strfind equivalents can't do more than 1 or are garbage to read and understand
                #found solution https://stackoverflow.com/questions/38974168/finding-entries-containing-a-substring-in-a-numpy-array needed list() so it'd do them all though. magical
                #k = strfind(AMPERE_raws{j,1},'.'); #get place of each .
                k = np.flatnonzero(np.core.defchararray.find(list(AMPERE_raws[j]),'.') != -1); #get place of each .
                #Time & Lat header info have 1 . per line and time is before lat when it occurs
                #long inputs have 5 . per line
                if( (len(k) == 1) & (cntr_lat >= AMPERE_lat_len) ): #if this is true it is a time entry
                    cntr_time = cntr_time + 1; #increment
                    cntr_lat = 0; #set to 0 to begin the fun
                    # str2double(AMPERE_raws{j,1}) ~= AMPERE_time_local(cntr_time) had issues with differing levels of precision
        #             if(  ismembertol(str2double(AMPERE_raws{j,1}),AMPERE_time_local(cntr_time)) ~= 1 )
        #                 error(['Expected time of ',num2str(AMPERE_time_local(cntr_time)),' is not same as read line of ',AMPERE_raws{j,1},' on line #',num2str(j),' of file ',fileName_AMPERE(i,:)]); #Report an error in data mismatch expected vs read
        #             end
        #         DISABLED ALL CHECKS TO IMPROVE SPEED FARTHER
    
                elif( (len(k) == 1) & (cntr_lat < AMPERE_lat_len) ): #if this is true it is a lat entry
                    cntr_lat = cntr_lat + 1; #increment lat cntr
                    # str2double(AMPERE_raws{j,1}) ~= AMPERE_lat_local(cntr_lat) might have same issue as time had with differing levels of precision
        #             if( ismembertol(str2double(AMPERE_raws{j,1}),AMPERE_lat_local(cntr_lat)) ~= 1 )
        #                 error(['Expected lat of ',num2str(AMPERE_lat_local(cntr_time)),' is not same as read line of ',AMPERE_raws{j,1},' on line #',num2str(j),' of file ',fileName_AMPERE(i,:)]); #Report an error in data mismatch expected vs read
        #             end
    
                else: #otherwise it is a long entry
                    
        #             C = strsplit(AMPERE_raws{j,1},' '); #split up the 5 numbers
        #             if( length(C) ~= 5 )
        #                 error(['Expected number of data pts for a data line is ',length(C),' not 5. The read line is ',AMPERE_raws{j,1},' on line #',num2str(j),' of file ',fileName_AMPERE(i,:)]); #Report an error in data mismatch expected vs read
        #             end
                    #C = np.core.defchararray.split(AMPERE_raws[j], sep=' '); #same but extra garbled
                    #C = AMPERE_raws[j].split(sep = ' '); #seperate by spaces (condensed below)
                    #C = list(filter(None,AMPERE_raws[j].split(sep = ' '))); #I had to use this awful thing to remove the empty stuff (condensed below)
                    
                    # AMPERE_data_save[cntr_long,0:AMPERE_data_num] = [dataAccuracy(x) for x in list(filter(None,AMPERE_raws[j].split(sep = ' ')))]; #place the 5 data points in, had to use awful list things    
                    
                    [ AMPERE_save['Ped'][cntr_long], AMPERE_save['Hall'][cntr_long], \
                        AMPERE_save['JH'][cntr_long], AMPERE_save['elec poten'][cntr_long], \
                        AMPERE_save['field-aligned current'][cntr_long] ] = \
                        [dataAccuracy(x) for x in list(filter(None,AMPERE_raws[j].split(sep = ' ')))]; #place the 5 data points in, had to use awful list things    
                    
                    cntr_long = cntr_long + 1; #increment
                #END IF
            #END FOR j
            
            #0 is the hdf5 format goal
            AMPERE_fileName_write = folder_fileNameFormat[ 0 ].replace('#DN', str(dateRange_dayNum_full[i,1]).zfill(3) ); #replace any #DN with the current day number
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#YR', str(dateRange_dayNum_full[i,0]) ); #replace any #YR with the current year
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#MO', str(dateRange_full[i,1]).zfill(2) ); #replace any #MO with the current month
            AMPERE_fileName_write = AMPERE_fileName_write.replace('#DY', str(dateRange_full[i,2]).zfill(2) ); #replace any #DY with current day
            
            #now save the data to use again
            keyz = list(AMPERE_save.keys()); #keys to the dict
            with h5py.File(settings_paths['data'] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName_write, 'w') as AMPERE_file:
                for j in range(0,len(keyz)):
                    if( AMPERE_save[keyz[j]].size > 1 ):
                        AMPERE_file.create_dataset(keyz[j], data=AMPERE_save[keyz[j]], compression="gzip"); #write that data
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
                    AMPERE_data[keyz[j]] = AMPERE_save[keyz[j]]; #record
                #END FOR j
            else:
                for j in range(0,len(keyz)):
                    if( np.sum(strfind(list(AMPERE_data.keys()),keyz[j])) > 0 ): #make sure key exists
                        AMPERE_data[keyz[j]] = np.hstack( (AMPERE_data[keyz[j]], AMPERE_save[keyz[j]]) ); #tack on
                    else:
                        AMPERE_data[keyz[j]] = AMPERE_save[keyz[j]]; #record if doesn't exist
                    #END IF
                #END FOR j
            #END IF
            #clear AMPERE_data_par AMPERE_data_temp
            toc = time.time() - tic; #end timing
            print("AMPERE data parsing and re-saving as HDF5 took: "+str(np.round(toc,2))+" sec / "+str(np.round(toc/60,2))+" min");
        #END IF                
    #END FOR i

    #Clear out AMPERE data less than the minimum
    # k = np.where(np.min(AMPERE_jouleHeating_plotLimValu) > AMPERE_data[:,AMPERE_jouleHeating_pos])[0]; #find entries less than the min plotting number (clear it up)
    # AMPERE_data = np.delete(AMPERE_data,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    #AMPERE_time = np.delete(AMPERE_time,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    #AMPERE_lat = np.delete(AMPERE_lat,k,axis=0); #delete entries that are less than the min plotting number (clear it up)
    #AMPERE_long = np.delete(AMPERE_long,k,axis=0); #delete entries that are less than the min plotting number (clear it up)

    # Calc total sec for the entire time period aligned to the zero hr [debuting year support! I hope it works]
    # AMPERE_data['time aligned'] = AMPERE_data['time']; #it's aligned time
    AMPERE_data['time aligned'] = np.copy(AMPERE_data['time']); #aligned time to zero hr
    AMPERE_data['time'] = AMPERE_data['time'] + dates['date range zero hr dayNum'][1]*86400; #sec, calc total sec for the day range not aligned
    
    if( AMPERE_coordType == 'mag' ):
        #!note picle is used instead of modifying the .h5 files because everything is jammed together and it would be annoying figuring out how to split it for a dead-end file format!
        import pickle as pkl
        magHash = 'AMPEREO_mag_'+str(dateRange_dayNum_full[0,:])+str(dateRange_dayNum_full[-1,:])+'_'+str(version_alg)+'.pkl'; #make the file name to use
        if( os.path.isfile(settings_paths['cache']+'\\'+magHash) == 1 ):
            #if the data already exists, load it in
            with open(settings_paths['cache']+'\\'+magHash, 'rb') as magData:
                AMPERE_dataMag = pkl.load(magData); #load in the data
            #END WITH
            #!note don't have to change lat/long b/c want it on an even grid!
            keyz = list(AMPERE_dataMag.keys()); #get some keys
            for j in range(0,len(keyz)):
                AMPERE_data[keyz[j]] = np.copy(AMPERE_dataMag[keyz[j]]); #place mag-shifted data in
            #END FOR k
            del AMPERE_dataMag; #clear mem
        else:
            import time
            tic_mag = time.time();
            print('Warning in GRITI_import_AMPERE_prepepared: coordType "mag" requested but mag coordinates have not been calculated for this day. Calculating. '+\
                  'Estimated time is '+textNice(np.round(tEst_mag/60,2))+' min on a good comp.');
            parallel_numCores = cpu_count(); #use multiprocess to get # of CPU cores (will strangle comp)
            alt4mag = 120.; #altitude used
            AMPERE_data_timeUnique = np.unique(AMPERE_data['time']); #get the time uniques
            AMPERE_lat_orig = np.copy(AMPERE_data['lat']); #copy for later
            AMPERE_long_orig = np.copy(AMPERE_data['long']); #copy for later
            with joblib.parallel_backend('loky'): #parallel backend uses 
                with joblib.Parallel(n_jobs=parallel_numCores,pre_dispatch=parallel_numCores,batch_size=1) as parallel_arbiter: 
                    #--- Parallel prep datetime object which is apparently hella hard ---
                    parallel_splitterIndexes = np.int64(np.round(np.linspace(0,AMPERE_data_timeUnique.size,parallel_numCores+1,endpoint=True))); #split up the indexes to be parallelized
                    parallel_list = []; #Prep
                    for j in range(0,parallel_splitterIndexes.size-1):
                        parallel_list.append([AMPERE_data['time'], AMPERE_data_timeUnique[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]],  AMPERE_data['year'][0]]);
                    #END FOR j
                    parallel_time4mag = parallel_arbiter(joblib.delayed(calc_datetimes)(j, k, l) for j, k, l in parallel_list); #will this not destroy the world?
                    time4mag = [None for i in range(0,AMPERE_data['time'].size)]; #preallocate
                    for i in range(0,AMPERE_data['time'].size): #recombine them, they're made of pieces of each (and Nones in between) - since they're objects I am in list purgatory still
                        for j in range(0,len(parallel_time4mag)):
                            if( parallel_time4mag[j][i] != None ):
                                time4mag[i] = parallel_time4mag[j][i]; #load it in
                            #END IF
                        #END FOR j
                    #END FOR i 
                    del parallel_time4mag
                    gc.collect(); #make that parallel_time4mag go away
                    
                    #--- Parallel calc mag coords ---
                    parallel_splitterIndexes = np.int64(np.round(np.linspace(0,AMPERE_data['lat'].size,parallel_numCores+1,endpoint=True))); #split up the indexes to be parallelized
                    parallel_list = []; #Prep
                    for j in range(0,parallel_splitterIndexes.size-1):
                        parallel_list.append([AMPERE_data['lat'][parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]], AMPERE_data['long'][parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]], alt4mag, time4mag[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]]]);
                    #END FOR j
                    parallel_mag_list = parallel_arbiter(joblib.delayed(convert_to_mag)(j, k, l, m, method_code='G2A') for j, k, l, m in parallel_list); #will this not destroy the world?
                    for j in range(0,parallel_splitterIndexes.size-1):
                        AMPERE_data['lat'][parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_mag_list[j][0].astype(dataAccuracy().dtype, order='C', casting='same_kind', subok=True, copy=False); #load it in
                        AMPERE_data['long'][parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_mag_list[j][1].astype(dataAccuracy().dtype, order='C', casting='same_kind', subok=True, copy=False); #load it in
                    #END FOR j
                #END WITH
            #END WITH
            del parallel_list, parallel_mag_list #save some mem 
            
            #Rebase to x y z b/c spherical no go good
            AMPERE_xOrig, AMPERE_yOrig, AMPERE_zOrig = geo_to_xyz(AMPERE_data['lat'], AMPERE_data['long'], np.ones(AMPERE_data['lat'].size)*120.); #calc
            AMPERE_x, AMPERE_y, AMPERE_z = geo_to_xyz(AMPERE_lat_orig, AMPERE_long_orig, np.ones(AMPERE_lat_orig.size)*120.); #calc
            
            #----- Prepare a list of dicts ------
            AMPERE_temp_numData = np.int64(120*50); #number of data per time step
            AMPERE_temp_timeData = AMPERE_time_range.size; #number of time steps
            # FLG_newMag = False;
            keyzCalc = list(AMPERE_data.keys()); #get the keys available
            parallel_list = []; #prep a parallel list
            for k in range(0,len(keyzCalc)):
                if( not keyzCalc[k] in ['year', 'dayNum', 'data rate', 'time', 'hour', 'min', 'sec', 'time aligned', 'lat', 'long', 'lat geo', 'long geo', 'mag info'] ):
                    AMPERE_saveTemp = {
                        # 'mag'+AMPERE_desiredResStr:{
                        #     'lat':AMPERE_save['mag'+AMPERE_desiredResStr]['lat'],
                        #     'long':AMPERE_save['mag'+AMPERE_desiredResStr]['long'],
                        #     },
                        'AMPERE_temp_timeData':AMPERE_temp_timeData,
                        # 'AMPERE_temp_numData':AMPERE_temp_numData,
                        'AMPERE_desired_numData':AMPERE_temp_numData,
                        # 'AMPERE_desired_latNum':50,
                        # 'AMPERE_desired_longNum':120,
                        # 'AMPERE_desiredResStr':AMPERE_desiredResStr,
                        'AMPERE_xOrig':AMPERE_xOrig,
                        'AMPERE_yOrig':AMPERE_yOrig,
                        'AMPERE_zOrig':AMPERE_zOrig,
                        'AMPERE_x':AMPERE_x,
                        'AMPERE_y':AMPERE_y,
                        'AMPERE_z':AMPERE_z,
                        'dataKey':keyzCalc[k],
                        'dataAccuracy':dataAccuracy,
                        # 'FLG_newMag':FLG_newMag,
                        keyzCalc[k]:AMPERE_data[keyzCalc[k]],
                        }; #prime this up
                    # if( FLG_newMag == True ):
                    #     AMPERE_saveTemp['mag orig'] = {
                    #         keyzCalc[k]:AMPERE_save['mag orig'][keyzCalc[k]],
                    #         'lat':AMPERE_save['mag orig']['lat'],
                    #         'long':AMPERE_save['mag orig']['long'],
                    #         }; #needed for new mag calcs
                    # #END IF
                    
                    parallel_list.append(AMPERE_saveTemp); #tack on this dict of data
                #END IF
            #END FOR k
            
            #----- calc in parallel -----
            # AMPERE_desired_results = pool.map(regrid_parallel, parallel_list); #parallel calc the results [this is awfully coded on Windows]
            # AMPERE_desired_results = joblib.Parallel(n_jobs=7)(joblib.delayed(regrid_parallel)(i) for i in parallel_list); #will this not destroy the world?
            with joblib.parallel_backend('loky', inner_max_num_threads=parallel_numCores//2): #this //2 assumption is built for hyperthreaded CPUs.
                AMPERE_desired_results = joblib.Parallel(n_jobs=parallel_numCores//2)(joblib.delayed(regrid_parallel)(i) for i in parallel_list); #will this not destroy the world?
            #END WITH
            
            #-----unpack the AMPERE_desired_results -----
            for j in range(0,len(AMPERE_desired_results)):
                AMPERE_data[parallel_list[j]['dataKey']] = AMPERE_desired_results[j]['geo'].copy(); #load in geo results
                # if( FLG_newMag == True ):
                #     AMPERE_data[parallel_list[j]['dataKey']] = AMPERE_desired_results[j]['mag'].copy(); #load in mag results
                # #END IF
            #END FOR j
            del parallel_list, AMPERE_desired_results #clean up mem
            
            #now save the data to use again
            keyz = list(AMPERE_data.keys()); #keys to the dict
            AMPERE_dataMag = {}; #prep new dict
            for j in range(0,len(keyz)):
                if( not keyz[j] in ['year', 'dayNum', 'data rate', 'time', 'hour', 'min', 'sec', 'time aligned', 'lat', 'long', 'lat geo', 'long geo', 'mag info', 'version'] ):
                    if( AMPERE_data[keyz[j]].size > 1 ):
                        AMPERE_dataMag[keyz[j]] = AMPERE_data[keyz[j]]; #record
                    #END IF
                #END IF
            #END FOR j
            with open(settings_paths['cache']+'\\'+magHash, 'wb') as magData:
                pkl.dump(AMPERE_dataMag, magData); #dump that data
            #END WITH
            del AMPERE_dataMag; #clean up
            
            tEst_mag = (time.time() - tic_mag); #record
            print('Conversion to coordType "mag" took '+textNice(np.round(tEst_mag/60,2))+' min.');
        #END IF
    #END IF
    # AMPERE_data['field-aligned current'][np.abs(AMPERE_data['field-aligned current']) < 0.2] = 0; # Remove low-energy field-aligned currents that may not be reliable

    return AMPERE_data #return the success
#END DEF

def convert_to_mag(lat4mag, long4mag, alt4mag, time4mag, method_code='G2A'):
    temp_lat = np.empty(lat4mag.size,dtype=lat4mag.dtype); #preallocate
    temp_long = np.empty(lat4mag.size,dtype=lat4mag.dtype); #preallocate
    for jk in range(0,lat4mag.size):
        temp_lat[jk], temp_long[jk], _ = aacgmv2.convert_latlon(lat4mag[jk], long4mag[jk], alt4mag, time4mag[jk], method_code=method_code); #converts from geographic to geomagnetic (AACGMv2) - vectorized _arr version has a memory leak atm
        tryCntr = 1; #reset cntr
        while( (np.isnan(temp_long[jk]) | np.isnan(temp_lat[jk])) & (tryCntr <= 20) ):
            #recalc at higher altitude to deal with error (idk how bad this is, but it allows the points to exist at least
            temp_lat[jk], temp_long[jk], _ = aacgmv2.convert_latlon(lat4mag[jk], long4mag[jk], alt4mag+200*tryCntr, time4mag[jk], method_code=method_code); #converts from geographic to geomagnetic (AACGMv2)
            tryCntr += 1 #increment
        #END IF
    #END FOR jh
    return temp_lat, temp_long; #this is a return I guess
#END DEF

def calc_datetimes(timez, timezUnique, year):
    import datetime
    time4mag = [None for i in range(0,timez.size)]; #preallocate
    for i in range(0,timezUnique.size):
        k = np.where(timezUnique[i] == timez)[0]; #get where the data pts at
        #no year support yet
        time4mag_dayNum = np.int32(timezUnique[i]//86400); #get days
        time4mag_hr = np.int32(np.mod(timezUnique[i],86400)//3600); #get hours
        time4mag_min = np.int32(np.mod(timezUnique[i],86400)//60-time4mag_hr*60); #get the minutes
        time4mag_sec = np.int32(np.mod(timezUnique[i],86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
        time4mag_dayNmonth = subfun_dayNum_to_date( (year,time4mag_dayNum) ); #badabing badaboom
        time4mag_datetime = datetime.datetime(year,time4mag_dayNmonth[0,1],time4mag_dayNmonth[0,2], \
                                      hour = time4mag_hr, minute = time4mag_min, second = time4mag_sec); #date time object for aacgmv2   
        for j in range(0,k.size):
            time4mag[k[j]] = time4mag_datetime; #load it up in the right spot
        #END FOR j
    #END FOR i
    return time4mag
#END DEF

#Rebase to xyz - spherical has made a powerful enemy
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
    # AMPERE_temp_numData = AMPERE_saveTemp['AMPERE_temp_numData'];
    AMPERE_desired_numData = AMPERE_saveTemp['AMPERE_desired_numData'];
    # AMPERE_desired_latNum = AMPERE_saveTemp['AMPERE_desired_latNum'];
    # AMPERE_desired_longNum = AMPERE_saveTemp['AMPERE_desired_longNum'];
    # AMPERE_desiredResStr = AMPERE_saveTemp['AMPERE_desiredResStr'];
    AMPERE_xOrig = AMPERE_saveTemp['AMPERE_xOrig'];
    AMPERE_yOrig = AMPERE_saveTemp['AMPERE_yOrig'];
    AMPERE_zOrig = AMPERE_saveTemp['AMPERE_zOrig'];
    AMPERE_x = AMPERE_saveTemp['AMPERE_x'];
    AMPERE_y = AMPERE_saveTemp['AMPERE_y'];
    AMPERE_z = AMPERE_saveTemp['AMPERE_z'];
    dataKey = AMPERE_saveTemp['dataKey'];
    dataAccuracy = AMPERE_saveTemp['dataAccuracy'];
    # FLG_newMag = AMPERE_saveTemp['FLG_newMag'];
    
    #precalc what can be precalcd
    # mag_string = 'mag'+AMPERE_desiredResStr;
    
    AMPERE_return_geo = np.empty(AMPERE_saveTemp['AMPERE_x'].size,dtype=dataAccuracy); #preallocate
    for j in range(0,AMPERE_temp_timeData):
        # if( FLG_newMag ):
        #     #RegridderZen interpolates OK but does not smooth enough (it is perfectly - and I mean perfectly - exact tho if that's ever important)
        #     #this interpolates just right, but it needs a lot of help and padding to prevent NaNs - longitude axis is rolled by 2 to align it correctly, idk why it gets unaligned and idk why its 2 to fix it but idk anymore
        #     #!!DO NOT KNOW IF 2 KEEPS WORKING FOR OTHER VALUES (other than 3)!!
        #     AMPERE_interper_dataOrig = AMPERE_saveTemp[dataKey][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
        #     AMPERE_interper_latOrig = AMPERE_saveTemp['mag orig']['lat'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
        #     AMPERE_interper_longOrig = AMPERE_saveTemp['mag orig']['long'][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData].copy();
        #     AMPERE_interper_longOrig[AMPERE_interper_longOrig<0] = AMPERE_interper_longOrig[AMPERE_interper_longOrig<0]+360; #make it roll around, idk if its req'd but doin it
        #     kr = AMPERE_interper_longOrig == 0; #get the edge
        #     AMPERE_interper_longOrig = np.append(AMPERE_interper_longOrig,np.repeat(360,kr.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
        #     AMPERE_interper_latOrig = np.append(AMPERE_interper_latOrig,AMPERE_interper_latOrig[kr]); #append on lat values (making it cyclical)
        #     AMPERE_interper_dataOrig = np.append(AMPERE_interper_dataOrig,AMPERE_interper_dataOrig[kr]); #append on data values (making it cyclical)
        #     kl = AMPERE_interper_longOrig == 345; #get the edge
        #     AMPERE_interper_longOrig = np.insert(AMPERE_interper_longOrig,0,np.repeat(-15,kl.sum())); #append on 360 to make extra data for the interper to work with (to make it cyclical in a way)
        #     AMPERE_interper_latOrig = np.insert(AMPERE_interper_latOrig,0,AMPERE_interper_latOrig[kl]); #append on lat values (making it cyclical)
        #     AMPERE_interper_dataOrig = np.insert(AMPERE_interper_dataOrig,0,AMPERE_interper_dataOrig[kl]); #append on data values (making it cyclical)
        #     AMPERE_interper_long = AMPERE_saveTemp[mag_string]['long'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData].copy();
        #     AMPERE_interper_long[AMPERE_interper_long<0] = AMPERE_interper_long[AMPERE_interper_long<0]+360; #make it roll around too
        #     AMPERE_return_mag[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = \
        #         np.roll(interpolate.griddata(np.vstack((AMPERE_interper_longOrig, AMPERE_interper_latOrig)).T, AMPERE_interper_dataOrig, \
        #             np.vstack((AMPERE_interper_long, AMPERE_saveTemp[mag_string]['lat'][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T,\
        #             method='linear').reshape(AMPERE_desired_longNum,AMPERE_desired_latNum),2,axis=0).ravel(); #one helluva call to interpolate
        # #END IF
        # if( np.all(AMPERE_desired_latLongSteps == None) ): #note the below works b/c for default stepping 'mag'+AMPERE_desiredResStr is mapped to 'mag orig' so its the same deal
        #     AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_saveTemp['mag orig'][dataKey][j*AMPERE_temp_numData:(j+1)*AMPERE_temp_numData], smoothing=0, kernel='linear');
        #     AMPERE_saveTemp[geo_string][dataKey][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
        # else:
        #this special mode is based off of the mag just calc'd and uses the same gridding function as previous
        #Rebase to x y z b/c spherical no go good
        AMPERE_interper_geo = interpolate.RBFInterpolator( np.vstack((AMPERE_xOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_yOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_zOrig[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T, AMPERE_saveTemp[dataKey][j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], smoothing=0, kernel='linear');
        AMPERE_return_geo[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData] = AMPERE_interper_geo(np.vstack((AMPERE_x[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_y[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData], AMPERE_z[j*AMPERE_desired_numData:(j+1)*AMPERE_desired_numData])).T); #interpolate!  
        #******WARNING****** THIS CODE MAY CRASH YOUR COMPUTER LIKE STRAIGHT SHUT OFF AND TURN BACK ON (BC COMP NOT 100% STABLE) BUT ITS DONE IT ON AN i5-4690K SYSTEM AND FX-8530 SYSTEM (i7-7700K was fine)******************************************************
        #END IF
    #END FOR j
    
    return {'geo':AMPERE_return_geo}
#END DEF