#Function to import TEC data from pre-proccessed sources
#RD on 12/10/2018
    #======FORMAT OF "other sources" DATA======
    #1 - time of day (days, UTC I think??)
    #2 - pierce point lat
    #3 - pierce point long
    #4 - TID estimate (dTEC est?)
    #5 - vTEC from GPS site (calc'd)
    #6 - site (guess site #?)

#To properly use, place pre-proccessed sources in their respective year folders in the data folder you've specified
    
#Info on stuff you can send:
#FLG_dataMix = 0; #0 prevents data sources to mixing to fill a time span, 1 allows data sources to mix to fill a time span
#FLG_dataPreference = 0; #preffered data type (by order as appended below at folder_fileNameFormat)
#TEC_dataLimPercent = 0.05; #0.05 = 5%, cut out times with very low data content (less than 5% of the mean data content #)

import numpy as np
import os
import h5py
from subfun_dayNum_to_date import subfun_dayNum_to_date
from subfun_daysInAYear import subfun_daysInAYear

#-----Testing variables-----
##Date range goes Month-Day-Year
##dateRange = np.array([[2013,5,8],[2013,5,10]],dtype="int16"); #dates are in int16 because they can be
#dateRange = np.array([[2013,5,6],[2013,5,8]],dtype="int16"); #dates are in int16 because they can be
##dateRange = np.array([[2014,12,31],[2015,1,1]],dtype="int16"); #for debug, check year success
##dates better go earlier -> later
##print("{}".format(dateRange))
#folder = [os.getcwd()]; #current working directory, leave it as this call usually
#folder.append('E:\Big Data'); #place to save data files to
##folder var structure: 0 = running folder, 1 = data folder
#FLG_dataPreference = 0; #preffered data type (by order as appended below)
#TEC_dataLimPercent = 0.05; #0.05 = 5%, cut out times with very low data content (less than 5% of the mean data content #)
#from subfun_date_to_dayNum import subfun_date_to_dayNum
#from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#FLG_dataMix=0;
#FLG_dataPreference=0;
#TEC_dataLimPercent = 0.05;
#plotLatRange = [35,50]; #latitude limit for plotting
##-90 to 90 is world, 35 to 50 is good for USA East Coast
#plotLongRange = [-85,-60]; #longitude limit for plotting
##-180 to 180 is world, -85 to -60 is good for USA East Coast

def GRITI_import_TEC_otherSources(dateRange_dayNum_full,dateRange_dayNum_zeroHr,folder,plotLatRange,plotLongRange,TEC_dataRate=30,TEC_minimumTimeGap=5,FLG_dataMix=0,FLG_dataPreference=0,TEC_dataLimPercent = 0.05,TEC_timeTolerance=0.1):
        
    #==============File System Locations==============
    #*************************************************
    #==============Filtered File Layout==============
    #Integer Layout
    #0 = Satellite ID [# that corresponds to GPS sat]
    locInt_sat = 0; #index where sat ID is
    #1 = Year timestamp [years]
    locInt_year = 1; #index where year timestamp is
    #2 = Day Number timestamp [days]
    locInt_dayNum = 2; #index where day number timestamp is
    #3 = Hour timestamp [hrs]
    locInt_hour = 3; #index where hour timestamp is
    #4 = Minute timestamp [mins]
    locInt_min = 4; #index where minute timestamp is
    #5 = Second timestamp [secs]
    locInt_sec = 5; #index where second timestamp is
    locInt_size = 6; #size of the int variable
    
    #Float Layout
    #0 = current time in day format [days] - does not support years
    locFloat_time = 0; #index where time in days is
    #1 = geodedic latitude [arcdeg]
    locFloat_lat = 1; #index where geodedic latitude is
    #2 = longitude [arcdeg]
    locFloat_long = 2; #index where longitude is
    #3 = elevation [deg]
    locFloat_elev = 3; #index where elevation is
    #4 = delta-TEC "kinda de-biased TEC" [TECU]
    locFloat_dTEC = 4; #index where delta-TEC is
    #5 = delta-TEC error [TECU]
    locFloat_dTECerr = 5; #index where the delta-TEC error is
    locFloat_size = 6; #size of float variable
    
    #String Layout
    #[] = Receiver Site Name (there's no dim on this one) - use .decode('UTF-8') to make it a string again
    locString_site = 0; #index where site name is
    locString_size = 1; #size of string layout
    
    print("\nTEC Other Sources - Date range requested (yr/day num format): {}/{} to {}/{}.".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
    
    #==============Constants Needed==============
    folder_TEC = 'TEC'; #name for the TEC folder
    folder_fileNameFormat = []; #prep a list
    # #YR for year
    # #DN for day num
    # #MO for month
    # #DY for day
    #.append more if more arise
    folder_fileNameFormat.append('tid_#YR_#DN_Anth.h5'); #0 - expected file name style
    #folder_fileNameFormat.append('otherFormat_#DY_#MO.h5'); #1 - tack on another (example)
    
                                 
    #==============Prep to look for the data==============
    
    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #get the full date range
                                 
    if( os.path.isdir(folder[1] + '\\' + folder_TEC) == 0 ): #check if TEC folder exists
        #if not, make it
        os.makedirs(folder[1] + '\\' + folder_TEC);
        print("NOTA BENE: Importing TEC Func - Created TEC directory: {}\n".format(folder[1] + '\\' + folder_TEC) );
    #END IF
        
    dateRange_uniqueYears = np.unique(dateRange_dayNum_full[:,0]); #get the unique years involved
    TEC_dataAmnt = len(dateRange_dayNum_full[:,1]); #get number of days needed to investigate
    TEC_dataAvail_perSource = np.zeros( (TEC_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for finding where data is available
    TEC_dataAvail_toUse = np.zeros( (TEC_dataAmnt,len(folder_fileNameFormat)) , dtype='int8');  #preallocate logical array, for recording which piece of data to use on the day
    #0 = no data available on required days, will quit
    #1 = note data is there, ready to use
    
    for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
        if( os.path.isdir(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) );
            print("NOTA BENE: Importing TEC Func - Created TEC subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) ));
        #END IF
    #END FOR
        
        
    #==============Look for data in expected naming formats==============
    for i in range(0,len(folder_fileNameFormat)): #loop through the different name formats    
        
        for j in range(0,TEC_dataAmnt): #loop through the different years needed
            
            TEC_fileName = folder_fileNameFormat[i].replace('#DN', str(dateRange_dayNum_full[j,1]).zfill(3) ); #replace any #DN with the current day number
            TEC_fileName = TEC_fileName.replace('#YR', str(dateRange_dayNum_full[j,0]) ); #replace any #YR with the current year
            TEC_fileName = TEC_fileName.replace('#MO', str(dateRange_full[j,1]).zfill(2) ); #replace any #MO with the current month
            TEC_fileName = TEC_fileName.replace('#DY', str(dateRange_full[j,2]).zfill(2) ); #replace any #DY with current day
            
            TEC_dataAvail_perSource[j,i] = os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileName); #record if data is there or not
        #END FOR j
    #END FOR i
        
    FLG_dataAvail_entireSpan = 0; #flag that needs to be set to 1 with this data finding alg
    if( FLG_dataMix == 1): #this allows for the data sources to mix to cover more availability
        TEC_dataAvail_entireSpan = np.any( TEC_dataAvail_perSource , axis=1 ); #see if each day has some availability
        
        if( np.all(TEC_dataAvail_entireSpan) ): #if all are 1 this is true
            #only do the work if we're good to go
            TEC_dataAvail_toUse = np.copy(TEC_dataAvail_perSource); #copy this over
            for i in range(0,TEC_dataAmnt): #check each entry for multiples, etc.
                
                if( np.sum(TEC_dataAvail_toUse[i,:]) > 1 ): #if the sum is greater than 1, choose one
                    
                    if( TEC_dataAvail_toUse[i,FLG_dataPreference] == 1 ): #easy, set it
                        TEC_dataAvail_toUse[i,:] = TEC_dataAvail_toUse[i,:] & 0; #0 it out
                        TEC_dataAvail_toUse[i,FLG_dataPreference] = 1; #set it to 1
                    else:
                        print("\n==============ERROR==============");
                        print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                        #return("No"); #return something that will def crash things
                        import sys #yolo import
                        sys.exit(); #def will crash
                        #don't feel like writing a disp, but def I didn't get around to writing alternate stuff
                    #END IF
                #END IF
            #END FOR i    
            
            FLG_dataAvail_entireSpan = 1; #set the flag to good
        #END IF    
    
    else: #otherwise one data source has to cover the entire span
        
        TEC_dataAvail_entireSpan = np.all( TEC_dataAvail_perSource , axis=0 ); #see if each source has complete availability
        
        if( np.any(TEC_dataAvail_entireSpan) ): #if any source has data for all dates, this will be true
            #don't do the work if it's not true
        
            if( np.sum(TEC_dataAvail_entireSpan) > 1): #make sure to choose one data source
                
                if( TEC_dataAvail_entireSpan[FLG_dataPreference] == 1 ): #easy, data type preferred is there
                    TEC_dataAvail_entireSpan = TEC_dataAvail_entireSpan & 0; #0 it out
                    TEC_dataAvail_entireSpan[FLG_dataPreference]  = 1; #put a 1 where the preffered data is there
                else: 
                    #don't feel like writing a disp, but def I didn't get around to writing alternate stuff
                    print("\n==============ERROR==============");
                    print("I didn't get around to writing a thing if option 0 wasn't the choice so this happened and you need to write cod eto fix it. Sorry. Crashin'");
                    #return("No"); #return something that will def crash things
                    import sys #yolo import
                    sys.exit(); #def will crash
                #END IF
            #END IF
                
            TEC_dataAvail_toUse[:,TEC_dataAvail_entireSpan] = 1; #set 
            
            FLG_dataAvail_entireSpan = 1; #set the flag to good
        #END IF
        
    #END IF
        
    if( FLG_dataAvail_entireSpan == 0 ): #make sure there's data
        print("\n==============ERROR==============");
        print("There is no data available from {}/{}/{} to {}/{}/{} in YR/M/D format.\nFLG_dataMix is set to {} (0 means all data comes from a single source, 1 means data can mix from sources).".format(dateRange_full[0,0],dateRange_full[0,1],dateRange_full[0,2],dateRange_full[-1,0],dateRange_full[-1,1],dateRange_full[-1,2],FLG_dataMix));
        print("Printing file name formats supported:");
        print("{}".format(folder_fileNameFormat)); #print for error
        print("Printing available data matrix (made of dates and file name formats):");
        print("{}".format(TEC_dataAvail_perSource)); #print for error - lets user know available days
        print("Will exit via returning FLG_TECloc = 2 so other alg can run\n");
        # return(0 , 0 , 0 , 2); #return something that will def crash things
        return(0, 2); #return failure
#        import sys #yolo import
#        sys.exit(); #def will crash
    #END IF
        
    
    #==============Read in the data==============
    
    for i in range(0,TEC_dataAmnt ):
        
        TEC_sourceIndex = np.where(TEC_dataAvail_toUse[i,:]==1)[0]; #get the location of the index (corresponds to which source)
        if( np.any(TEC_sourceIndex == FLG_dataPreference) ): #if any of the ones that have data are the chosen source, then use it
            TEC_sourceIndex = np.int64(FLG_dataPreference); #just set it n forget it
        else: #otherwise just use the next one in line (don't have a hierarchy yet)
            TEC_sourceIndex = TEC_sourceIndex[0]; #just get the first one in line
        #END IF
        TEC_fileName = folder_fileNameFormat[ TEC_sourceIndex ].replace('#DN', str(dateRange_dayNum_full[i,1]).zfill(3) ); #replace any #DN with the current day number
        TEC_fileName = TEC_fileName.replace('#YR', str(dateRange_dayNum_full[i,0]) ); #replace any #YR with the current year
        TEC_fileName = TEC_fileName.replace('#MO', str(dateRange_full[i,1]).zfill(2) ); #replace any #MO with the current month
        TEC_fileName = TEC_fileName.replace('#DY', str(dateRange_full[i,2]).zfill(2) ); #replace any #DY with current day
        #prep the name needed
        
        #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
        with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileName, 'r') as TEC_file:
            #Read in the HDF5 file from the name list
            
            if( TEC_sourceIndex == 0 ): #supports file setup for source #0
                TEC_dataTemp = TEC_file['/tecs'][()]; #data in /tecs
            else:
                print("\n==============ERROR==============");
                print("Unsupported source index {} used! That's bad!\nGonna crash now.");
                #return("No"); #return something that will def crash things
                import sys #yolo import
                sys.exit(); #def will crash
            #END IF
        #END WITH
        
        if( TEC_sourceIndex == 0 ): #supports file setup for source #0        
            if i == 0: #if first one initialize the array
                
                #-----Remove lat/long combos not in the range requested-----
                k = np.where( (TEC_dataTemp[:,2] > np.max(plotLongRange)) | (TEC_dataTemp[:,2] < np.min(plotLongRange)) | \
                    (TEC_dataTemp[:,1] > np.max(plotLatRange)) | (TEC_dataTemp[:,1] < np.min(plotLatRange)) )[0];
                TEC_dataTemp = np.delete(TEC_dataTemp,k,axis=0); #delete out of lat/long range stuff
                
                #-----Delete low data time slots per day-----
                #Since numpy arrays copy themselves to be deleted, this saves memory as it gets *bad* with multiple days doing it on TEC_float, etc.
                (timeUnique,k) = np.unique(TEC_dataTemp[:,0],return_counts=True); #days, gathers unique times
                #also get the number of occurance for each set
                
                #Cut off time stamps with very little data in the range selected
                TEC_dataAvgNum = len(TEC_dataTemp[:,0])/len(timeUnique); #average data per time
                TEC_dataLim = np.round(TEC_dataLimPercent*TEC_dataAvgNum); #min number before eject time set
                
                k = k < TEC_dataLim; #get the time uniques to delete for lack of data
                k = np.in1d(TEC_dataTemp[:,0],timeUnique[k]); #find where time matches times to delete
                k = np.nonzero(k)[0]; #get the indicies as numpy wants that
                
                TEC_dataTemp = np.delete(TEC_dataTemp,k,axis=0); #delete low data # stuff
                
                #always on below
                k = np.where( np.int16(TEC_dataTemp[:,0]) != (dateRange_dayNum_full[i,1]-1) )[0]; #get the indexes of the days that don't match
                
                TEC_dataTemp = np.delete(TEC_dataTemp,k,axis=0); #delete indexes of the days taht don't match
                        
                TEC_dataTemp[:,0] = TEC_dataTemp[:,0] + 1; #days, adjust because the data comes starting at 0 days not 1 days
#                
#                #I was wrong - don't have satellite info so can't split the tracks
#                #this code replaces the low time slot deletion code by moving the low time slots to regular rates
#                #~~~COMPACT THE TIME STAMPS IN THE LISN DATA~~~
#                #basically, it needs to meet the lowest common demoninator (30 sec per satellite reading)
#                #alternatively - the 30 sec data could be repeated to "fill in" the data better, but that would req a lot more memory (more than 16GB) that I don't have in this day and age, but in the future this is where you'd implement something like that
#                #10 sec data will be compacted by averaging 1-10, 11-20, and 21-30 (shown as 10, 20, 30 in the data) into just 30 at the
#                #position of the 30 sec data pt
#                #now it could be also an average of the position datas, but I'm going with the last position for now since I don't know how the position is chosen to start with so the easier way is the way to go!
#                
#                TEC_dataTemp_hr = np.int16((TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24); #record hour for data
#                TEC_dataTemp_min = np.int16(( (TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24 - TEC_dataTemp_hr)*60); #record min for data
#                TEC_dataTemp_sec = np.int16(np.round((( (TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24 - TEC_dataTemp_hr)*60 - TEC_dataTemp_min)*60)); #record sec for data
#                
#                TEC_dataRate_allowedStamps = np.arange(0,60,TEC_dataRate); #sec, make an array of allowed second timestamps
#                #calculate a more acturate time array using 64 bits instead of the standard 32 bits
#                TEC_time64 = np.float64(dateRange_dayNum_full[i,1]) + TEC_dataTemp_hr/24 + TEC_dataTemp_min/1440 + TEC_dataTemp_sec/86400; #days, calculate hour/min/sec into days and add to the current day
#                timeDeltas = np.abs(np.diff(TEC_time64)*86400); #sec, get the time deltas
#                newSatTracks = np.where( timeDeltas > TEC_minimumTimeGap*60 )[0] + 1; #sec, identify time gaps that indicate a new satellite pass
#                uniqueSats, uniqueSats_reverseIndex = np.unique(TEC_int_temp[:,locInt_sat],return_inverse=True); #get the unqiues and an array the size of the original array that is made of the indexes that correspond to unique
#                newSats = np.where(np.abs(np.diff(uniqueSats_reverseIndex)) != 0)[0] + 1; #identify when a new sat is seen
#                uniqueSites, uniqueSites_reverseIndex = np.unique(TEC_str_temp,return_inverse=True); #get the unqiues and an array the size of the original array that is made of the indexes that correspond to unique
#                newSites = np.where(np.abs(np.diff(uniqueSites_reverseIndex)) != 0)[0] + 1; #identify when a new site occurs
#                
#                newSatTracks = np.unique(np.concatenate( (np.array( (0,) ),newSatTracks, newSats, newSites, np.array( (uniqueSats_reverseIndex.size,) )), axis=0 )); #jam them all together, get the uniques
#                #don't need to do any heuristics to see what is in where, if there's a new one in newSats or newSites it's a new gap but the time stamps are close enough not to register on the time diff one
#                #also tack on the start (0) and end (size of the data array)
#                
#                for i in range(0 , newSatTracks.size-1):
#                    startIndex = newSatTracks[i];
#                    endIndex = newSatTracks[i+1];
#            
#                    track_timeDelta = np.int16(np.round(np.median(np.diff(TEC_time64[startIndex:endIndex]))*86400)); #sec, get the time delta for the track    
#                    
#                    if( track_timeDelta < TEC_dataRate ): 
#                        track_timeSecs = TEC_int_temp[startIndex:endIndex,locInt_sec]; #sec, get the seconds for the track
#                        
#                        #get where the second time stamps are the allowed values
#                        track_Matches = np.zeros( track_timeSecs.size, dtype=np.bool_ ); #preallocate truth table
#                        for j in range(0,TEC_dataRate_allowedStamps.size):
#                            track_Matches = track_Matches | (TEC_dataRate_allowedStamps[j] == track_timeSecs); #get the seconds that match the allowed times
#                        #END FOR j
#                        track_Matches_where = np.concatenate( (np.array( (0,) ),np.where(track_Matches == True)[0]) ); #get the indexes, put 0 at the front
#                        
#                        track_timeDeltas = np.int16(np.round(np.diff(TEC_time64[startIndex:endIndex])*86400)); #sec, get the time delta for the track
#                        track_Gaps_where = np.where(track_timeDeltas > track_timeDelta)[0]+1; #get the indexes for where there are gaps in the data
#                        if( track_Gaps_where.size == 0 ):
#                            track_Gaps_where = -1776; #signify a non-index to avoid errors
#                        #END IF 
#                        
#                        track_dTEC = TEC_float_temp[startIndex:endIndex,locFloat_dTEC]; #TECU, pull out the dTEC for the track
#                        #now work through it and combine as needed
#                        for j in range(0,track_Matches_where.size-1):
#                            if( np.any( (track_Gaps_where > track_Matches_where[j]) & (track_Gaps_where < (track_Matches_where[j+1]+1)) ) == 1 ):
#                                track_Gaps_whereIndiv = np.where( (track_Gaps_where > track_Matches_where[j]) & (track_Gaps_where < (track_Matches_where[j+1]+1)) )[0][0]; #get the index
#                                track_dTEC[track_Gaps_where[track_Gaps_whereIndiv]-1] = np.mean(track_dTEC[track_Matches_where[j]:track_Gaps_where[track_Gaps_whereIndiv]]); #TECU, average them together
#                                track_dTEC[track_Matches_where[j+1]] = np.mean(track_dTEC[track_Gaps_where[track_Gaps_whereIndiv]:track_Matches_where[j+1]+1]); #TECU, average them together
#                                if( track_dTEC[track_Matches_where[j]:track_Gaps_where[track_Gaps_whereIndiv]].size == 0 ):
#                                    sys.exit();
#                                #END IF
#                                if( track_dTEC[track_Gaps_where[track_Gaps_whereIndiv]:track_Matches_where[j+1]+1].size == 0 ):
#                                    sys.exit();
#                                #END IF
#                                #now make sure that the point at [track_Gaps_where[track_Gaps_whereIndiv]-1] is adjusted to be a correct value (treat it like the ending on a false alg below)
#                                if( track_Matches[track_Gaps_where[track_Gaps_whereIndiv]-1] == False ): #fix it up
#                                    #check if the last one is an off-time, and in that case
#                                    #push it to the next time slot <- decided to do this b/c predicting could end up rough and deleting is less data while the time scales are short compared to the 1 hour periodicity I'm looking at
#                                    TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
#                                    TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
#                                    TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
#                                    track_shiftSec = np.where(TEC_dataRate_allowedStamps60 > track_timeSecs[track_Gaps_where[track_Gaps_whereIndiv]-1])[0][0]; #get where the allowed indexes are greater than the current index, and get the first one always
#                                    TEC_int_temp[startIndex+track_Gaps_where[track_Gaps_whereIndiv]-1,locInt_sec] = TEC_dataRate_allowedStamps60[track_shiftSec]; #choose the right one
#                                    if( TEC_int_temp[startIndex+track_Gaps_where[track_Gaps_whereIndiv]-1,locInt_sec] == 60 ):
#                                        #force it back to 0
#                                        TEC_int_temp[startIndex+track_Gaps_where[track_Gaps_whereIndiv]-1,locInt_sec] = 0; #force it back to 0 from 60
#                                    #END
#                                    track_Matches[track_Gaps_where[track_Gaps_whereIndiv]-1] = True; #set the last value to true so it isn't NaN'd
#                                #END IF
#                            else: #otherwise no gaps good to go
#                                track_dTEC[track_Matches_where[j+1]] = np.mean(track_dTEC[track_Matches_where[j]:track_Matches_where[j+1]+1]); #TECU, average them together
#                                if( track_dTEC[track_Matches_where[j]:track_Matches_where[j+1]+1].size == 0 ):
#                                    sys.exit();
#                                #END IF
#                            #END IF
#                        #END FOR j
#                        if( track_Matches[-1] == False ):
#                            #check if the last one is an off-time, and in that case
#                            #push it to the next time slot <- decided to do this b/c predicting could end up rough and deleting is less data while the time scales are short compared to the 1 hour periodicity I'm looking at
#                            track_dTEC[-1] = np.mean(track_dTEC[track_Matches_where[-1]+1:track_dTEC.size]); #TECU, average the trailing ones
#                            if( track_dTEC[track_Matches_where[-1]+1:track_dTEC.size].size == 0 ):
#                                sys.exit();
#                            #END IF
#                            TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
#                            TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
#                            TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
#                            track_shiftSec = np.where(TEC_dataRate_allowedStamps60 > TEC_int_temp[endIndex-1,locInt_sec])[0][0]; #get where the allowed indexes are greater than the current index, and get the first one always
#                            TEC_int_temp[endIndex-1,locInt_sec] = TEC_dataRate_allowedStamps60[track_shiftSec]; #choose the right one
#                            if( TEC_int_temp[endIndex-1,locInt_sec] == 60 ):
#                                #force it back to 0
#                                TEC_int_temp[endIndex-1,locInt_sec] = 0; #force it back to 0 from 60
#                            #END
#                            track_Matches[-1] = True; #set the last value to true so it isn't NaN'd
#                        #END IF
#                        
#                        #keep the lat/long the same as the original values
#                        
#                        #now NaN the values that didn't match for easy deletion
#                        track_dTEC[np.logical_not(track_Matches)] = np.nan; #nan it up
#                        
#                        #now record, clear out NaNs later
#                        #The memory is tied, so it's automagic
#            #                TEC_float_temp[startIndex:endIndex,locFloat_dTEC] = np.copy(track_dTEC); #TECU, record the new condensed TECU
#                    elif( track_timeDelta > TEC_dataRate ):
#                        #need copying code for this
#                        print("ERROR DON\'T HAVE CODE FOR THIS YET");
#                        sys.exit();
#                    #END IF
#                #END FOR i
#                del TEC_time64, timeDeltas #clear more memory
#                
#                #remove NaNs - which indicate time periods that were averaged to be condensed
#                TEC_nonNans = np.where( np.logical_not(np.isnan(TEC_float_temp[:,locFloat_dTEC])) == True)[0]; #find the not NaNs, get the index (maybe faster?)
                
                #-----Save the data as ints, floats, strings-----
                TEC_int = np.zeros( (TEC_dataTemp.shape[0],locInt_size),dtype='int16'); #preallocate
                TEC_float = np.zeros( (TEC_dataTemp.shape[0],locFloat_size),dtype='float32'); #preallocate
#                TEC_str = np.array( '~~~~' ,dtype='S4'); #prep preallocating (this has to be gibberish as the source doesn't have)
#                TEC_str = np.repeat(TEC_str,TEC_dataTemp.shape[0]); #preallocate fully
                TEC_str = np.zeros( (TEC_dataTemp.shape[0],),dtype='S4'); #preallocate
                
                #write ints
                TEC_int[:,locInt_year] = dateRange_full[i,0]; #record year for data
                TEC_int[:,locInt_dayNum] = dateRange_dayNum_full[i,1]; #record daynumber for data
                TEC_int[:,locInt_hour] = np.int16((TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24); #record hour for data
                TEC_int[:,locInt_min] = np.int16(( (TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24 - TEC_int[:,locInt_hour])*60); #record min for data
                TEC_int[:,locInt_sec] = np.int16(np.round((( (TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24 - TEC_int[:,locInt_hour])*60 - TEC_int[:,locInt_min])*60)); #record sec for data
                                
                #write floats
                TEC_float[:,locFloat_dTEC] = TEC_dataTemp[:,3]; #record sTEC data
                TEC_float[:,locFloat_time] = TEC_dataTemp[:,0]; #record time in day format
                TEC_float[:,locFloat_lat] = TEC_dataTemp[:,1]; #record geodedic latitude
                TEC_float[:,locFloat_long] = TEC_dataTemp[:,2]; #record longitude
                TEC_float[:,locFloat_elev] = 90; #set elev to 90 for all data since none exists (prevents incorrect culling)
                
                #write strings
                TEC_str[:] = TEC_dataTemp[:,5].astype('S4'); #record the string data
                
            else: #if next iterations stick data onto array
                
                #-----Remove lat/long combos not in the range requested-----
                k = np.where( (TEC_dataTemp[:,2] > np.max(plotLongRange)) | (TEC_dataTemp[:,2] < np.min(plotLongRange)) | \
                    (TEC_dataTemp[:,1] > np.max(plotLatRange)) | (TEC_dataTemp[:,1] < np.min(plotLatRange)) )[0];
                TEC_dataTemp = np.delete(TEC_dataTemp,k,axis=0); #delete out of lat/long range stuff
                
                #-----Delete low data time slots per day-----
                #Since numpy arrays copy themselves to be deleted, this saves memory as it gets *bad* with multiple days doing it on TEC_float, etc.
                (timeUnique,k) = np.unique(TEC_dataTemp[:,0],return_counts=True); #days, gathers unique times
                #also get the number of occurance for each set
                
                #Cut off time stamps with very little data in the range selected
                TEC_dataAvgNum = len(TEC_dataTemp[:,0])/len(timeUnique); #average data per time
                TEC_dataLim = np.round(TEC_dataLimPercent*TEC_dataAvgNum); #min number before eject time set
                
                k = k < TEC_dataLim; #get the time uniques to delete for lack of data
                k = np.in1d(TEC_dataTemp[:,0],timeUnique[k]); #find where time matches times to delete
                k = np.nonzero(k)[0]; #get the indicies as numpy wants that
                
                TEC_dataTemp = np.delete(TEC_dataTemp,k,axis=0); #delete low data # stuff
                
                #always on below
                k = np.where( np.int16(TEC_dataTemp[:,0]) != (dateRange_dayNum_full[i,1]-1) )[0]; #get the indexes of the days that don't match
                
                TEC_dataTemp = np.delete(TEC_dataTemp,k,axis=0); #delete indexes of the days that don't match
                
                TEC_dataTemp[:,0] = TEC_dataTemp[:,0] + 1; #days, adjust because the data comes starting at 0 days not 1 days
                                
                #-----Save the data as ints, floats, strings-----
                TEC_int = np.concatenate( (TEC_int , np.zeros( (TEC_dataTemp.shape[0],locInt_size),dtype='int16')) , axis=0 ); #mass add 0's
                TEC_float = np.concatenate( (TEC_float , np.zeros( (TEC_dataTemp.shape[0],locFloat_size),dtype='float32')) , axis=0); #mass add 0's
#                TEC_str = np.concatenate( (TEC_str ,  np.repeat(np.array( '~~~~' ,dtype='S4'),TEC_dataTemp.shape[0])) , axis=0); #mass add ~~~~'s
                TEC_str = np.concatenate( (TEC_str ,  np.zeros( (TEC_dataTemp.shape[0],),dtype='S4')) , axis=0); #mass add 0's
                
                #write ints
                TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_year] = dateRange_full[i,0]; #record year for data
                TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_dayNum] = dateRange_dayNum_full[i,1]; #record year for data
                TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_hour] = np.int16((TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24); #record hour for data
                TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_min] = np.int16(( (TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24 - TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_hour])*60); #record min for data
                TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_sec] = np.int16(np.round((( (TEC_dataTemp[:,0] - dateRange_dayNum_full[i,1])*24 - TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_hour])*60 - TEC_int[ TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0] ,locInt_min])*60)); #record sec for data
                
                #write floats
                TEC_float[TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0],locFloat_dTEC] = TEC_dataTemp[:,3]; #record d-vTEC data
                TEC_float[TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0],locFloat_time] = TEC_dataTemp[:,0]; #record time in day forma
                TEC_float[TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0],locFloat_lat] = TEC_dataTemp[:,1]; #record geodedic latitude
                TEC_float[TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0],locFloat_long] = TEC_dataTemp[:,2]; #record longitude
                TEC_float[TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0],locFloat_elev] = 90; #set it again
                                
                #write strings
                TEC_str[TEC_int.shape[0]-TEC_dataTemp.shape[0]: TEC_int.shape[0]] = TEC_dataTemp[:,5].astype('S4'); #record the string data
            #END IF
                
        else:
            print("\n==============ERROR==============");
            print("Unsupported source index {} used! That's bad!\nGonna crash now.");
            #return("No"); #return something that will def crash things
            import sys #yolo import
            sys.exit(); #def will crash
        #END IF
            
    #END FOR i   
    
    #-----Declare some constants used to force the data to the required data rate-----
    TEC_dataRate_allowedStamps = np.arange(0,60,TEC_dataRate); #sec, make an array of allowed second timestamps
    TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
    TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
    TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
    
    allowedStamps_range = np.zeros( (TEC_dataRate_allowedStamps.size,2)); #preallocate
    for m in range(0,TEC_dataRate_allowedStamps.size):
        allowedStamps_range[m, :] = np.array( (TEC_dataRate_allowedStamps[m]-TEC_dataRate*TEC_timeTolerance,TEC_dataRate_allowedStamps[m]+TEC_dataRate*TEC_timeTolerance) ); #get the current allowed time stamp tolerance range
        if( np.any(allowedStamps_range[m,:] < 0) ):
            allowedStamps_range[m, allowedStamps_range[m,:] < 0] = allowedStamps_range[m, allowedStamps_range[m,:] < 0]+60; #keep in the 0-59 sec range
            current_timeInTolerance = (TEC_int[:,locInt_sec] >= allowedStamps_range[m,0]) | (TEC_int[:,locInt_sec] <= allowedStamps_range[m,1]); #get the times within the time stamp tolerance range (so for 0 sec, 57 to 3 sec is taken as 0 sec)
            if( np.any(TEC_int[:,locInt_sec] >= allowedStamps_range[m,0]) ):
                kj = TEC_int[:,locInt_sec] >= allowedStamps_range[m,0]; #get where the minute rolls over
                TEC_int[kj,locInt_min] += 1; #min, the minutes are incremented as well if the minute rolled over
                if( np.any(TEC_int[:,locInt_min] == 60) ):
                    kj = TEC_int[:,locInt_min] == 60; #get where the hour rolls over
                    TEC_int[kj,locInt_min] = 0; #min, set to 0
                    TEC_int[kj,locInt_hour] += 1; #hour, the hours are incremented as well if the hour rolled over
                    if( np.any(TEC_int[:,locInt_hour] == 24) ):
                        kj = TEC_int[:,locInt_hour] == 24; #get where the day rolls over
                        TEC_int[kj,locInt_hour] =  0; #hour, set to 0
                        TEC_int[kj,locInt_dayNum] += 1; #day, the days are incremented as well if the day rolled over
                        #-----Leap Year Detection-----
                        #deal with day limit is based on leap year or not
                        dayLim = np.ones(TEC_int[:,locInt_dayNum].shape,dtype=np.int16)*365; #day number, get the day number limits as 365
                        #adjust leap years to 366
                        leapYears = (np.mod(TEC_int[:,locInt_year],4) == 0) & (np.mod(TEC_int[:,locInt_year],100) != 0) & (np.mod(TEC_int[:,locInt_year],400) == 0)
                        dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
                        if( np.any(TEC_int[:,locInt_dayNum] >= dayLim) ):
                            kj = TEC_int[:,locInt_dayNum] >= dayLim; #get where the year rolls over
                            TEC_int[kj,locInt_dayNum] =  1; #day, set to 1 [the 1st day]
                            TEC_int[kj,locInt_year] += 1; #year, the years are incremented as well if the year rolled over
                        #END IF
                    #END IF
                #END IF
            #END IF
        else:
            current_timeInTolerance = (TEC_int[:,locInt_sec] >= allowedStamps_range[m,0]) & (TEC_int[:,locInt_sec] <= allowedStamps_range[m,1]); #get the times within the time stamp tolerance range (so for 30 sec, 27 to 33 sec is taken as 30 sec)
        #END IF
        TEC_int[current_timeInTolerance,locInt_sec] = TEC_dataRate_allowedStamps[m]; #set the times within the tolerance range to the allowed time stamp for that range
    #END FOR m
    
    #-----make sure hour 24 or min 60 doesn't exist (they hsould be rolled over) - found issue via day 129 data showing up b/c it was recorded as d128/h24/m60...------       
    #---CHECK SECONDS---
    kj = TEC_int[:,locInt_sec] >= 60; #find incorrect time keeping
    if( np.sum(kj) > 0 ):
        #increment the main minute variable where seconds are 60 or more
        TEC_int[kj,locInt_min] += 1; #min, increment time by 1
        TEC_int[kj,locInt_sec] -= 60; #sec, remove 60 time units from the time keeping
        if( np.any(TEC_int[:,locInt_sec] >= 60) ):
            print('ERROR: TIME KEEPING IS REALLY BAD, 60+ SECOND TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ SECOND TIME STAMPS');
            sys.crash(); #call a crash (crash doesn't exist so it crashes)
        #END IF
    #END IF
    #---CHECK MINUTES---
    kj = TEC_int[:,locInt_min] >= 60; #find incorrect time keeping
    if( np.sum(kj) > 0 ):
        #increment the main hour variable where minutes are 60 or more
        TEC_int[kj,locInt_hour] += 1; #hour, increment time by 1
        TEC_int[kj,locInt_min] -= 60; #min, remove 60 time units from the time keeping
        if( np.any(TEC_int[:,locInt_min] >= 60) ):
            print('ERROR: TIME KEEPING IS REALLY BAD, 60+ MINUTE TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ MIN TIME STAMPS');
            sys.crash(); #call a crash (crash doesn't exist so it crashes)
        #END IF
    #END IF
    #---CHECK HOURS---
    kj = TEC_int[:,locInt_hour] >= 24; #find incorrect time keeping
    if( np.sum(kj) > 0 ):
        #increment the main day variable where hours are 24 or more
        TEC_int[kj,locInt_dayNum] += 1; #day number, increment time by 1
        TEC_int[kj,locInt_hour] -= 24; #hour, remove 24 time units from the time keeping
        if( np.any(TEC_int[:,locInt_hour] >= 24) ):
            print('ERROR: TIME KEEPING IS REALLY BAD, 24+ HOUR TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 48+ HOUR TIME STAMPS');
            sys.crash(); #call a crash (crash doesn't exist so it crashes)
        #END IF
    #END IF
    #---CHECK DAYS---
    #deal with day limit is based on leap year or not
    dayLim = np.ones(TEC_int[:,locInt_dayNum].shape,dtype=np.int16)*365; #day number, get the day number limits as 365
    #adjust leap years to 366
    leapYears = (np.mod(TEC_int[:,locInt_year],4) == 0) & (np.mod(TEC_int[:,locInt_year],100) != 0) & (np.mod(TEC_int[:,locInt_year],400) == 0)
    dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
    kj = TEC_int[:,locInt_dayNum] >= dayLim; #find incorrect time keeping
    if( np.sum(kj) > 0 ):
        #increment the main year variable where day number is equal to the day number limit or higher than it
        TEC_int[kj,locInt_year] += 1; #year, increment time by 1
        TEC_int[kj,locInt_dayNum] -= dayLim; #hour, remove the day number limit time units from the time keeping
        if( np.any(TEC_int[:,locInt_dayNum] >= dayLim) ):
            print('ERROR: TIME KEEPING IS REALLY BAD, 365/366 DAY TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 365*2/366*2 DAY TIME STAMPS');
            sys.crash(); #call a crash (crash doesn't exist so it crashes)
        #END IF
    #END IF
    #---luckily years have no limits (that we know of??)---
    #---DELETE data not on the right day (could have incrememted past the right day due to stuff and things like a 24 hour/60 minute time on the current day---
    TEC_logical_onDay = (TEC_int[:,locInt_dayNum] <= dateRange_dayNum_full[-1,1]) & (TEC_int[:,locInt_year] <= dateRange_dayNum_full[-1,0]) &\
        (TEC_int[:,locInt_dayNum] >= dateRange_dayNum_full[0,1]) & (TEC_int[:,locInt_year] >= dateRange_dayNum_full[0,0]); #find when the day reported is the day we want [and the year we want]
    TEC_int = TEC_int[TEC_logical_onDay,:]; #keep only the good stuff
    TEC_float = TEC_float[TEC_logical_onDay,:]; #keep only the good stuff
    TEC_str = TEC_str[TEC_logical_onDay]; #keep only the good stuff
    #---UPDATE float date variable---
    TEC_float[:,locFloat_time] = np.float32(TEC_int[:,locInt_dayNum]) + np.float32(TEC_int[:,locInt_hour])/24 + np.float32(TEC_int[:,locInt_min])/1440 + np.float32(TEC_int[:,locInt_sec])/86400; #days, calculate hour/min/sec into days and add to the current day
    
    #==============Filter Data for Small Data Amounts==============
    
    #--- Convert from monolithic to dict ---
    import gc
    gc.collect(); #force the mem clear (just in case there's something lurkin cause we're gonna double up real quick)
    TEC_data = {
        'lat':TEC_float[:,locFloat_lat],
        'long':TEC_float[:,locFloat_long],
        'elev':TEC_float[:,locFloat_elev],
        'dTEC':TEC_float[:,locFloat_dTEC],
        'dTECerr':TEC_float[:,locFloat_dTECerr],
        };
    del TEC_float #clear mem
    gc.collect(); #force the mem clear
    TEC_data['sat'] = TEC_int[:,locInt_sat];
    TEC_data['year'] = TEC_int[:,locInt_year];
    TEC_data['dayNum'] = TEC_int[:,locInt_dayNum];
    TEC_data['hour'] = TEC_int[:,locInt_hour];
    TEC_data['min'] = TEC_int[:,locInt_min];
    TEC_data['sec'] = TEC_int[:,locInt_sec];
    del TEC_int #clear mem
    gc.collect(); #force the mem clear
    TEC_data['time'] = np.int32(TEC_data['dayNum'])*86400 + np.int32(TEC_data['hour'])*3600 + np.int32(TEC_data['min'])*60 + np.int32(TEC_data['sec']); #set up time in seconds
    TEC_data['time aligned'] = np.int32(TEC_data['dayNum']-dateRange_dayNum_zeroHr[1]-(dateRange_dayNum_zeroHr[0]-TEC_data['year'])*subfun_daysInAYear(TEC_data['year']))*86400 + np.int32(TEC_data['hour'])*3600 + np.int32(TEC_data['min'])*60 + np.int32(TEC_data['sec']); #set up time in seconds
    TEC_data['site'] = TEC_str; #it's not a vector
    del TEC_str #clear mem
    gc.collect(); #force the mem clear
    #get the attributes
    # TEC_data['pierceAlt'] = TEC_pierceAlt;
    # TEC_data['paddedDayMissing'] = TEC_paddedDayMissing;
    TEC_data['data rate'] = TEC_dataRate;
    # TEC_data['maxAmplAllowed'] = TEC_maxAmpAllowed;
    # TEC_data['TEC_savgolFiltPeriodSec'] = TEC_savgolFiltPeriodHr; #this needs to be updated
    
    # return(TEC_int, TEC_float, TEC_str, 1) #return the stuff for the function! (note that the 1 is the flag, keeps the alg from running that downloads data)
    return(TEC_data, 1) #return the stuff for the function! (note that the 1 is the flag, keeps the alg from running that downloads data)
        
        
        
        
    
    
    
    
    
    
    
        
        
        
