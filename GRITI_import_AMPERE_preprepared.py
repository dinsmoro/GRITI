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
from subfun_strfind import strfind
from subfun_daysInAYear import subfun_daysInAYear

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


def GRITI_import_AMPERE_preprepared(dates,folder,dateRange_dayNum_full_adj,FLG_dataMix=0,FLG_dataPreference=0,FLG_float64=0):
    
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
    
    #AMPERE_jouleHeating_plotLimValu = [1,5]; #erg/(cm^2*sec), low and high limits for plotting joule heating - taken from example video (seems Joule Heating low end limit is 1 - nothing less than 1)
    #AMPERE_jouleHeating_pos = 2; #set the position of the joule heating data
    # AMPERE_data_num = 5; #number of data entries per line (+3 for time/lat/long, added later on)
     #0 Pedersen Conductance / 1 Hall Conductance / 2 Joule Heat (ergs/(cm^2*sec)) / 3 Electric Potential / 4 Field-Algined Current / 5 time (hr) / 6 latitude / 7 longitude
    
                                 
    #==============Prep to look for the data==============                                 
    if( os.path.isdir(folder[1] + '\\' + folder_AMPERE) == 0 ): #check if AMPERE folder exists
        #if not, make it
        os.makedirs(folder[1] + '\\' + folder_AMPERE);
        print("NOTA BENE: Importing AMPERE Func - Created AMPERE directory: {}\n".format(folder[1] + '\\' + folder_AMPERE) );
    #END IF
    
    dateRange_uniqueYears = np.unique(dateRange_dayNum_full[:,0]); #get the unique years involved
    AMPERE_dataAmnt = len(dateRange_dayNum_full[:,1]); #get number of days needed to investigate
    AMPERE_dataAvail_perSource = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for finding where data is available
    AMPERE_dataAvail_toUse = np.zeros( (AMPERE_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for recording which piece of data to use on the day
    #0 = no data available on required days, will quit
    #1 = note data is there, ready to use
    
    for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
        if( os.path.isdir(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_uniqueYears[i]) ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_uniqueYears[i]) );
            print("NOTA BENE: Importing AMPERE Func - Created AMPERE subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_uniqueYears[i]) ));
        #END IF
    #END FOR
    
    #==============Look for data in expected naming formats==============
    for i in range(0,len(folder_fileNameFormat)): #loop through the different name formats    
        
        for j in range(0,AMPERE_dataAmnt): #loop through the different years needed
            
            AMPERE_fileName = folder_fileNameFormat[i].replace('#DN', str(dateRange_dayNum_full[j,1]) ); #replace any #DN with the current day number
            AMPERE_fileName = AMPERE_fileName.replace('#YR', str(dateRange_dayNum_full[j,0]) ); #replace any #YR with the current year
            AMPERE_fileName = AMPERE_fileName.replace('#MO', str(dateRange_full[j,1]).zfill(2) ); #replace any #MO with the current month (padded w/ 0's so always 2 long)
            AMPERE_fileName = AMPERE_fileName.replace('#DY', str(dateRange_full[j,2]).zfill(2) ); #replace any #DY with current day (padded w/ 0's so always 2 long)
            
            AMPERE_dataAvail_perSource[j,i] = os.path.isfile(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[j,0]) + '\\' + AMPERE_fileName); #record if data is there or not
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
            with h5py.File(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, 'r') as AMPERE_file:
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
                    print(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+ \
                        ' file version is '+str(AMPERE_file.attrs['version'])+' which is less than the current file version of '+str(version_alg)+'!'+\
                        '\nThat\'s bad! Renaming file to add a \'_old\' at the end, then gonna crash. Re-run and code will remake the file with the latest alg.\nI didn\'t code in a way to recover from this. Srry </3.');
                    os.rename(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, \
                        folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName+'_old'); #rename
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
            with open(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName, 'r') as AMPERE_file: #open file pointer (with for safety cause life is hard apparently)
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
            with h5py.File(folder[1] + '\\' + folder_AMPERE + '\\' + str(dateRange_full[i,0]) + '\\' + AMPERE_fileName_write, 'w') as AMPERE_file:
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

    return AMPERE_data #return the success