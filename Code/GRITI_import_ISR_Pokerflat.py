#Function to import ISR data in the Pokerflat Madgrial format
#RD on 12/10/2018
#-----ISR File Layout------
    #Integer Layout
    #0 = 
    
    #Float Layout
    #0 =
    
    #String Layout
    #[] = 

#To properly use, place pre-proccessed sources in their respective year folders in the data folder you've specified
    
#Info on stuff you can send:
#FLG_ISRcodeType = 97; #97 for alternating code (default) or 115 for single pulse ISR modes
#bs = 2; hr, high-pass cutoff period (any higher period is removed - 2 is usual)
import numpy as np
import glob
import h5py
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_monthWord_to_num import subfun_monthWord_to_num
from scipy import signal
from Code.subfun_strstr import strstr

#-----Testing variables-----
#import os
#import matplotlib.pyplot as plt
#from Code.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
###Date range goes Month-Day-Year
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
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units
#pointAltitude = 300; #km, altitude for centering upon
#filter_cutoffPeriod = 2; #hr, high-pass cutoff period (any higher period is removed - 2 is usual)
##PFSIR_reqDataVersion = 1; #the file endings 001, 002, 003, 004, 005 are the data versions that have various versions of
##correction factors used or even different pulse lengths/other stuff used
##1 corresponds to 001, etc.

def GRITI_import_ISR_Pokerflat(dateRange_dayNum_full,folder,dateRange_dayNum_zeroHr,pointAltitude,filter_cutoffPeriod = 2*3600):

    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #get the full date range
    
    filenamePossibles = glob.glob(folder[1]+"\\ISR\\"+str(dateRange_dayNum_full[0,0])+"\\pfa*"); #make a list of file names GLOB GLOB
    filenameInRange = np.zeros(len(filenamePossibles),dtype=np.bool_); #get a logical vector for if a filename is what is needed
    filenamePossibles = np.asarray(filenamePossibles); #make it a numpy array so it's not awful to work with
    
    for i in range(0,len(filenamePossibles)):
        try: #gonna try to read the file and get the year - if it fails whoops
            with h5py.File(filenamePossibles[i], 'r') as testFile:
                testFile_keys = list(testFile.keys()); #get the folder names in the hdf5 file
                for j in range(0,len(testFile_keys)): #run through all the keys
                    if( testFile_keys[j].lower() == "metadata" ):
                        testFile_metadataLocale = j; #record the place
                    #END IF
                #END FOR j
                testFile_metadataKeys = list(testFile[testFile_keys[testFile_metadataLocale]]); #choose the data folder
                for j in range(0,len(testFile_metadataKeys)): #run through all the keys
                    if( testFile_metadataKeys[j].lower() == "experiment parameters" ):
                        testFile_expParamsLocale = j; #record the place
                    #END IF
                #END FOR j
                testFile_expParamsKeys = list(testFile[testFile_keys[testFile_metadataLocale]+'/'+testFile_metadataKeys[testFile_expParamsLocale]]); #choose the data folder
                for j in range(0,len(testFile_expParamsKeys)): #run through all the keys
                    if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "start time" ):
                        testFile_expParams_start = testFile_expParamsKeys[j][1].decode('UTF-8'); #record the start time
                    #END IF
                    if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "end time" ):
                        testFile_expParams_end = testFile_expParamsKeys[j][1].decode('UTF-8'); #record the end time
                    #END IF
    #                if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "instrument latitude" ):
    #                    testFile_expParams_lat = np.float64(testFile_expParamsKeys[j][1].decode('UTF-8')); #record the latitude
                    #END IF
    #                if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "instrument longitude" ):
    #                    testFile_expParams_long = np.float64(testFile_expParamsKeys[j][1].decode('UTF-8')); #record the longitude
                    #END IF
                #END FOR j
            #END WITH
        except OSError and NameError:
            print("\n==============~Warning~==============");
            print("The current file: "+filenamePossibles[i]+" |did not follow supported pfa...hdf5 file formating. Skipping it.");
        #END TRYING
        testFile_expParams_start_Date = testFile_expParams_start[0:testFile_expParams_start.find(' ')]; #get the start date
        testFile_expParams_end_Date = testFile_expParams_end[0:testFile_expParams_end.find(' ')]; #get the start date
        
        testFile_expParams_start_Date_Split = np.int16(np.array(testFile_expParams_start_Date.split('-'))); #split by dash
        testFile_expParams_end_Date_Split = np.int16(np.array(testFile_expParams_end_Date.split('-'))); #split by dash
        
        filenameInRange[i] = (dateRange_full[-1,0] >= testFile_expParams_start_Date_Split[0]) & (dateRange_full[-1,1] >= testFile_expParams_start_Date_Split[1]) & \
        (dateRange_full[-1,2] >= testFile_expParams_start_Date_Split[2]) & (dateRange_full[0,0] <= testFile_expParams_end_Date_Split[0]) & \
        (dateRange_full[0,1] <= testFile_expParams_end_Date_Split[1]) & (dateRange_full[0,2] <= testFile_expParams_end_Date_Split[2]);
    #    ( strstr(filenamePossibles[i],str(PFSIR_reqDataVersion).zfill(3)).size > 0 ); #make sure data start is before or during the first day desired and the data end is after or during the last day desired
    #END FOR i
    
    if( np.sum(filenameInRange) == 0 ):
        print("\n==============ERROR==============");
        print("Of the files found in the given directory \'"+folder[1]+"\\ISR\\"+"\':");
        print({}.format( filenamePossibles ));
    #    print("None of them covered the date range needed OR used the data version \'"+str(PFSIR_reqDataVersion).zfill(3)+"\': "+str(dateRange_full[0,:])+" to "+str(dateRange_full[-1,:]));
        print("None of them covered the date range needed: "+str(dateRange_full[0,:])+" to "+str(dateRange_full[-1,:]));
        print("Crashin' on purpose, sorry!");
        import sys; #prep for destruction
        sys.exit(); #rekt
    #END IF
    
    filenameInRange = np.where(filenameInRange == True)[0]; #convert to indexes
    filenamePossibles = filenamePossibles[filenameInRange]; #get only the file name that are really possible
    
    for k in range(0,filenamePossibles.size):
        
        filenamePath = filenamePossibles[filenameInRange[k]].item(); #get the current file path
    
        with h5py.File(filenamePath, 'r') as filename:
            filename_expParamsKeys = list(filename["/Metadata/Experiment Parameters/"]); #choose the data folder
            for j in range(0,len(filename_expParamsKeys)): #run through all the keys
                if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "start time" ):
                    ISR_start = filename_expParamsKeys[j][1].decode('UTF-8'); #record the start time
                #END IF
                if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "instrument latitude" ):
                    ISR_lat = np.float64(filename_expParamsKeys[j][1].decode('UTF-8')); #record the latitude
                #END IF
                if( testFile_expParamsKeys[j][0].decode('UTF-8').lower() == "instrument longitude" ):
                    ISR_long = np.float64(filename_expParamsKeys[j][1].decode('UTF-8')); #record the longitude
                    if( ISR_long > 180 ):
                        ISR_long = 360 - ISR_long; #degc, convert from 0-360 longitude range to -180 to 180 longitude range
                    #END IF
                #END IF
            #END FOR j
            
            ISR_start = np.int16(np.array(ISR_start[0:ISR_start.find(' ')].split('-'))); #get the start date and split by dash
            
            filename_expParamsKeys = list(filename["/Data/Array Layout/"]); #get the data types
            
            #PFISR is like a software defined radio (physically it's like a series of tubes, not a big truck)
            #so the one PFISR square has data for multiple elevation angles at the same time, instead of being 1 big dish that has 1 elevation angle
            PFISR_beamIDNum = len(filename_expParamsKeys); #the number of elevation angles means number of datasets
            
            PFISR_beamID_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_POPL_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_SNR_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_range_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_time_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_el_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_az_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_vel_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_dopplar_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_height_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_POPL_bp_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data
            PFISR_SNR_bp_tmp =  PFISR_beamIDNum*[None]; #preallocate list that holds a data        
            
            for j in range(0, PFISR_beamIDNum):
                PFISR_beamID_tmp[j] = np.int64(filename_expParamsKeys[j].split('=')[1].strip()); #get the beam ID number
                
                PFISR_POPL_tmp[j] = 10**(filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/2D Parameters/popl'][:]); #POPL copied (alt to SNR - is electron density e-/m^3)
                try:
                    PFISR_SNR_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/2D Parameters/snp3'][:]; #SNR copied
                except:
                    PFISR_SNR_tmp[j] = np.zeros(PFISR_POPL_tmp[j].shape); #do zeros instead because the file lacks SNP3 which is common
                #END TRY
                PFISR_range_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/range'][:]; #range, km
                PFISR_time_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/timestamps'][:]; #time, sec from 1970
                PFISR_el_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/1D Parameters/elm'][:]; #elevation angle, deg
                PFISR_az_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/1D Parameters/azm'][:]; #azimuth angle, deg
                try:
                    PFISR_vel_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/2D Parameters/vo'][:]; #velocity, km/s
                except:
                    PFISR_vel_tmp[j] = np.zeros(PFISR_POPL_tmp[j].shape); #do zeros instead because the file lacks velocity
                #END TRY
                try:
                    PFISR_dopplar_tmp[j] = filename['/Data/Array Layout/'+filename_expParamsKeys[j]+'/2D Parameters/vdopp'][:]; #dopplar eh?
                except:
                    PFISR_dopplar_tmp[j] = np.zeros(PFISR_POPL_tmp[j].shape); #do zeros instead because the file lacks dopplar
                #END TRY
                
                if( PFISR_POPL_tmp[j].shape[1] != PFISR_range_tmp[j].size ):
                    #the range is on the wrong axis
                    PFISR_POPL_tmp[j] = np.transpose(PFISR_POPL_tmp[j]); #flip it around
                    PFISR_SNR_tmp[j] = np.transpose(PFISR_SNR_tmp[j]); #flip it around
                    PFISR_vel_tmp[j] = np.transpose(PFISR_vel_tmp[j]); #flip it around
                    PFISR_dopplar_tmp[j] = np.transpose(PFISR_dopplar_tmp[j]); #flip it around
                #END IF
                if( PFISR_range_tmp[j].size == PFISR_time_tmp[j].size ):
                    #this is the worst case scenario because we need to find a better way to know what axis is vs what
                    #it's like is it time vs range like it should be, or is it range vs time? WE CAN'T KNOW BASED ON CURRENT DETECTION METHODS!
                    print('AH INCONSISTENT DATA SHAPES AH FIX AH WHY OH WHY')
                    sys.exit();
                #END IF
                
            #END FOR j
        #END WITH
        
        # fileinfo  = h5info(filename); #use to check what data is up to, or use HDFView
        
        #fns=textscan(filename, '%3s %2d %2d %2d %10s'); #rips apart the filename to get the date (for mlh130506g.001.hdf5 it is mlh 13 05 06 g.001.hdf5)
        #year=double(2000+fns{1,2}); #Year from 2nd piece of file name (with same example, 13 or 2013)
        #month=double(fns{1,3}); #Month from 3rd piece of file name (with same example, 05)
        #day=double(fns{1,4}); #Day from 4th piece of file name (with same example, 06)
        #date=[year, month, day, 00, 00, 00]; #makes the date the year, month, day from the file name and says time is 0 hr 0 min 0 sec
        
        
        
        # There is a offset in time. The starting time is the total number of
        # seconds counted from 01 Jan 1970.
        # This experiment started on day month year. Therefore the difference between
        # day month year at 00:00:00 and 01 Jan 1970 at 00:00:00 will give the offset.
        # Julian day is calculated for both dates. The difference between Julian days
        # gives the total number days (taking care of all leap years) and hence the
        # offset time in seconds...
        
        #    ISR_start_dayNum = subfun_date_to_dayNum(ISR_start)[0]; #get day num version
        
        try: #try to open the leap second file
            from datetime import date #used for getting current date
            
            with open(folder[1]+"\\Supporting\\"+"leap-seconds.txt","r") as leapFile: #prep to read 
                leapRead = leapFile.read().split("\n"); #read file, split by \n
            #END WITH
            leapReadDel = np.zeros(len(leapRead),dtype=np.int64); #prep list, keep it in python's weird systems
            for i in range(0,len(leapRead)):
                if( len(leapRead[i]) == 0 ):
                    leapReadDel[i] = i;  #record strings to delete
                #END IF
            #END FOR i
            leapReadDel = np.delete(leapReadDel,np.where(leapReadDel == 0)[0]); #remove the 0's we didn't need
            leapRead = np.delete(leapRead,leapReadDel); #now it's a numpy array
            
            #now gotta make sure file we have is current enough
            leapCurrentDay = subfun_date_to_dayNum(np.int16(np.array(date.today().strftime("%Y %m %d").split(" "))))[0]; #get current day
            leapExpire = leapRead[0][ (leapRead[0].find(":")+3):(len(leapRead[0])+1) ].split(" "); #get the date of expiry
            leapExpire = subfun_date_to_dayNum(np.int16(np.array( [leapExpire[2], subfun_monthWord_to_num(leapExpire[1]) , leapExpire[0]] )))[0]; #create a number version
            
            if( leapExpire[0] < leapCurrentDay[0] ): #if the expiry year is less than the current one, easy catch
                raise ValueError("Raising an error so the except triggers and it goes and gets new data from the website, as the current data is old! Current date:\n"+str(leapCurrentDay)+"\nExpiring date:\n"+str(leapExpire));
            elif( leapExpire[0] == leapCurrentDay[0] ): #if the year is the same
                if( leapExpire[1] < leapCurrentDay[1] ): #check the days within them
                    raise ValueError("Raising an error so the except triggers and it goes and gets new data from the website, as the current data is old! Current date:\n"+str(leapCurrentDay)+"\nExpiring date:\n"+str(leapExpire));
            #END IF
            
            #if we made it through, time to read that data in and then use it!
            leapDates = np.zeros( [(len(leapRead)-2),3] , dtype=np.int16); #prep the dates (minus 2 because the first line is file expiration and the)
            #THE FIRST LINE IN THE LISTED LEAP SECOND DATES IS THE INITIAL ONE THAT DOESN'T COUNT AS ONE OK wiki shows it kinda
            for i in range(2,len(leapRead)):
                leapDates[i-2,:] = np.int16(leapRead[i].split(" ")); #convert the dates as needed
            #END FOR i
            
        except: #fail and do this to get it from 
            from urllib.request import urlopen #only need it here
            web_leapSecond = "https://www.ietf.org/timezones/data/leap-seconds.list"; #build site to go to
            web_page = urlopen(web_leapSecond); #get raw HTML
            web_htmlContent = web_page.read().decode("UTF-8").split("\n"); #read off the HTML from whatever the page holder is
            leapStart = -1; #prep flag/counter recorders
            leapEnd = -1;
            leapExpire = -1;
            for i in range(0,len(web_htmlContent)):
                if( len(web_htmlContent[i]) > 0 ): #prevent "none lines" from being read
                    if( (web_htmlContent[i][0] != "#" ) & ( leapStart == -1 ) ): #catch the first instance of no #'s
                        leapStart = i; #record the place where the comments end
                    #END IF
                    if( (web_htmlContent[i][0] == "#" ) & ( leapStart != -1 ) & ( leapEnd == -1 ) ): #catch where the #'s start again
                        leapEnd = i-1; #record the place where the comments start again
                    #END IF
                    if( "File expires" in web_htmlContent[i] ): #catch the file expires line
                        leapExpire = i; #record that line too
                #END IF
            #END FOR i
            
            #write a text file with the info
            with open(folder[1]+"\\Supporting\\"+"leap-seconds.txt","w") as leapFile: #prep to write 
                leapFile.write(web_htmlContent[leapExpire]+"\n\r"); #write the line with the expiring date
                leapStrings = web_htmlContent[leapStart:(leapEnd+1)]; #pull out the leap strings
                leapDates = np.zeros( [len(leapStrings),3] , dtype=np.int16); #prep the dates
                for i in range(0,len(leapStrings)):
                    leapStrings_temp = leapStrings[i][(leapStrings[i].find("#")+2):(len(leapStrings[i])+1)].split(" ");
                    leapStrings_temp[1] = str(subfun_monthWord_to_num(leapStrings_temp[1])); #convert to a number
                    leapDates[i,:] = np.int16(np.array( [leapStrings_temp[2] , leapStrings_temp[1] , leapStrings_temp[0]] )); #convert to a numpy array as ints
                    leapFile.write( " ".join([leapStrings_temp[2] , leapStrings_temp[1] ,  leapStrings_temp[0]])+"\n\r" ); #write this string in
                #END FOR i
            #END WITH
            leapDates = np.delete(leapDates,0,axis=0); #delete the first entry because it's the START of UTC and NOT an actual added LEAP SECOND
        #END TRY
        
        leapDates_dayNum = subfun_date_to_dayNum(leapDates); #convert it
        
        #Go hard
        #secondsSince1970 = 
        
        secSince_yearDiff = np.arange(1970,dateRange_dayNum_zeroHr[0],1,dtype="int16") #get the year range from 1970 to dataYear-1 b/c slicing but it's OK here, I guess.
        secSince_yearDiffDays = np.zeros(len(secSince_yearDiff),dtype="int16"); #use to record the days in the years
        for i in range(0, len(secSince_yearDiff) ):
            #Leap Year Detection
            if np.mod(secSince_yearDiff[i],4) == 0: #leap year
                #LEAP YEAR! Possibly.
                if (np.mod(secSince_yearDiff[i],100) == 0) and (np.mod(secSince_yearDiff[i],400) != 0):
                    #NO LEAP YEAR
                    #Leap Year Skipped Detected - next will be 2100
                    secSince_yearDiffDays[i] = 365; #no leap year
                else:
                    #Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
                    secSince_yearDiffDays[i] = 366; #it was a leap year
                #END IF
            else: #no leap year if this
                #no leap year
                secSince_yearDiffDays[i] = 365; #no leap year
            #END IF 
        #END FOR i loop per year
            
        secSince = (np.sum(secSince_yearDiffDays) + (dateRange_dayNum_zeroHr[1]-1))*24*3600; #sec, seconds since 1970 Jan 1st.
        #-1 on the day so we don't count the seconds in the current day - only up to the current day!
        
        #NOW ADD IN SOME LEAP SECONDS IT'S SUPER EASY
        secSince_leapWhere = np.where(leapDates_dayNum[:,0] < dateRange_dayNum_zeroHr[0])[0]; #get leap seconds that are relevant for the date
        secSince_leapWhereCheck = np.where(leapDates_dayNum[:,0] == dateRange_dayNum_zeroHr[0])[0]; #get leap seconds that are relevant for the date
        for i in range(0, len(secSince_leapWhereCheck)):
            if( leapDates_dayNum[secSince_leapWhereCheck[i],1] <= dateRange_dayNum_zeroHr[1] ): #find days that are before or on the start day
                secSince_leapWhereCheck[i] = -1; #use this flag to remove them as those leap seconds happened
            #END IF
        #END FOR i
            
        secSince_leapWhereCheck = np.delete(secSince_leapWhereCheck,np.where(secSince_leapWhereCheck == -1)); #delete ones that are ok
        if( len(secSince_leapWhereCheck) != 0 ): #if it's not empty, there's something to delete
            secSince_leapWhereDel = np.empty(len(secSince_leapWhereCheck),dtype=np.int64); #preallocate
            for i in range(0, len(secSince_leapWhereCheck)):
                secSince_leapWhereDel[i] = np.where(secSince_leapWhereCheck[i] == secSince_leapWhere)[0];
            #END FOR i
            secSince_leapWhere = np.delete(secSince_leapWhere,secSince_leapWhereDel); #delete the ones from the list
        #END IF
            
        secSince = secSince + len(secSince_leapWhere); #add in the number of leap seconds needed
        
        for j in range(0, PFISR_beamIDNum):
            # PFISR_time_tmp[j] = (PFISR_time_tmp[j]-secSince)/(3600*24) + dateRange_dayNum_zeroHr[1]; #removes offset from 1970, converts to days and then adjusts to the zero hour day
            #sec version below
            PFISR_time_tmp[j] = PFISR_time_tmp[j]-secSince + dateRange_dayNum_zeroHr[1]*86400; #removes offset from 1970, converts to days and then adjusts to the zero hour day
        #END FOR j
            
        # ========================Replacing NAN Values with 0======================
        for j in range(0, PFISR_beamIDNum):
            PFISR_SNR_tmp[j][np.isnan(PFISR_SNR_tmp[j])] = 0; #replaces NAN with 0
            
            PFISR_POPL_tmp[j][np.isnan(PFISR_POPL_tmp[j])] = 0; #Replaces NAN with 0
            
            PFISR_range_tmp[j][np.isnan(PFISR_range_tmp[j])] = 0; #Replaces NAN with 0
            
            PFISR_vel_tmp[j][np.isnan(PFISR_vel_tmp[j])] = 0; #Replaces NAN with 0
            
            PFISR_dopplar_tmp[j][np.isnan(PFISR_dopplar_tmp[j])] = 0; #Replaces NAN with 0
        #END FOR j
        
        # ============Convert Time to EDT (UT-4 hrs)==============================
        # PFISR_time = PFISR_time - 4/24;
        # MISA_time = MISA_time - 4/24;
        # Disabled currently
        
        # ===============PFISR and MISA log of SNR Calcs=========================
        for j in range(0, PFISR_beamIDNum):
            PFISR_SNR_tmp[j][PFISR_SNR_tmp[j] <= 0] = 0.1; #Removes negative numbers/0's and replaces with 0.1
            #PFISR_SNR_log = np.log10(PFISR_SNR);
        #END FOR j
        
        # ==============PFISR and MISA Height Calcs (from range + angle)=========
        for j in range(0, PFISR_beamIDNum):
            if( np.abs(np.median(PFISR_el_tmp[j])-np.mean(PFISR_el_tmp[j])) < 0.1 ):
                PFISR_height_tmp[j] = PFISR_range_tmp[j]/((1. + ((np.tan( (90. - np.median(PFISR_el_tmp[j]))*np.pi/180 ))**2))**0.5);
            else:
                print("WARNING: For file "+filenamePath+"\nPFISR elevation mean is "+str(np.mean(PFISR_el_tmp))+" but median is "+str(np.median(PFISR_el_tmp))+" which indicates PFISR is moving but it shouldn't be able to - using median");
                PFISR_height_tmp[j] = PFISR_range_tmp[j]/((1. + ((np.tan( (90. - np.median(PFISR_el_tmp[j]))*np.pi/180 ))**2))**0.5);
            #END IF
        #END FOR j
        
        # ===============High Pass Filtering to Excentuate Hourly Periods=========
        for j in range(0, PFISR_beamIDNum):
            #These are unused, but accentuate the periods desired
            #Original file used highpass_fir.m, it has been merged
            # lp=1; # Lower period (in hrs)
            # hp=2; # Higher period (in hrs)
            # lf=(1/hp);  # Lowpass frequency corner (1/hr)
            # hf=(1/lp);  # Highpass frequency corner (1/hr)
            
            
            bs = (1/filter_cutoffPeriod);   # highpass cuttoff frequency (not sure what to make of it)
            #I think it is related to lf above (which is 1/2 )
            #ls = (1/(1/2)); #low-pass cutoff frequency
            
            PFISR_delt = np.median(np.diff(PFISR_time_tmp[j])); #Calculates a time delta based on avg of all the deltas
            
            # ===============Highpass filtering on original PFISR SNR================
            
            n=42; # order of the Hamming window used in the custom function (uses 43 in reality)
            # c = 3.32*pi; #Hamming constant
            # M = n/2;
            # bandwidth = c/M;
            #The above was included to investigate the bandwidth - 1/2. Related to
            #bs/lf?
            
            fp = bs; # stop band stoppin freq (or pass band passing freq, depending on how you look at it)
            
            f= 1/(PFISR_delt); #1/sec, the sampling frequency, based off of the time delta calc'd
            
            wp = 2*fp/f; # Normalizing the frequencies (Matlab is 0 to 1)
            #wl = 2*ls/f; #norm this
            #Calculation of filter coefficients
            # [b,a]=fir1(n,wp,'high'); #This takes the order and lower cut off
            #Uses the default hamming window
            # frequency ORIG
            #W = np.hanning(n+1);
            #W = signal.get_window('hann',n+1);
            # [b,a] = fir1(n,[wp,wl],W); #option for band-pass (30 min to 2 hr)
            #[b,a] = fir1(n,wp,'high',W); #just high-pass (2 hr and lower OK)
            b = signal.firwin(n+1, wp, window='hann', pass_zero=False); #just high-pass (2 hr and lower OK)
            a = 1; #for FIR a is 1
            #Applys Hanning Window for shiggles and giggles
            #NOTE: b is *very slightly* different from MATLAB output.
            
            #equiv to MATLAB's freqz(b,a,512) for checking
            #plt.figure();
            #w, h = signal.freqz(b, worN=512);
            #plt.plot( (w/np.pi), np.log10(np.abs(h))*20, linewidth=2);
            #plt.xlabel('Normalized Freqeuency (x pi rad/sample)');
            #plt.ylabel('Magnitude (dB)');
            #plt.title('Frequency Response');
            #plt.grid(True);
            #plt.show(); #req to make plot show up
            
            try:
                PFISR_SNR_bp_tmp[j] = signal.filtfilt(b,a,PFISR_SNR_tmp[j],axis=0,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
                #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
            except ValueError:
                PFISR_SNR_bp_tmp[j] = np.ones(PFISR_SNR_tmp[j].shape)*np.nan; #fill with NaN's because the PFISR was turned on/off and it's time frame is actually too short to high-pass correctly
            #END TRY
            try:
                PFISR_POPL_bp_tmp[j] = signal.filtfilt(b,a,PFISR_POPL_tmp[j],axis=0,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
                #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
            except ValueError:
                PFISR_POPL_bp_tmp[j] = np.ones(PFISR_POPL_tmp[j].shape)*np.nan; #fill with NaN's because the PFISR was turned on/off and it's time frame is actually too short to high-pass correctly
            #END TRY
            # PFISR_SNR_bp = PFISR_SNR; %No filter option to see what it really looks
            # like
        #END FOR j
               
        # =============== Spectral Estimation of Filtered PFISR SNR at 300 km ===
        # 47 corresponds to 300 km
        # 25 corresponds to 200 km
        # 3 corresponds to 100 km
        # 70 corresponds to 400 km
        # 92 corresponds to 500 km (all approximate)
        for j in range(0, PFISR_beamIDNum):
            PFISR_filtHeight = np.where( np.min(np.abs(PFISR_height_tmp[j] - pointAltitude)) == np.abs(PFISR_height_tmp[j] - pointAltitude) )[0][0]; #See 'table' above (300 km was original)
            
            # Removing the mean from PFISR SNR time series at 300 km....
            PFISR_SNR_bp_tmp[j][:,PFISR_filtHeight] = PFISR_SNR_bp_tmp[j][:,PFISR_filtHeight] - np.mean(PFISR_SNR_bp_tmp[j][:,PFISR_filtHeight]);
            #    PFISR_POPL_bp[:,PFISR_filtHeight] = PFISR_POPL_bp[:,PFISR_filtHeight] - np.mean(PFISR_POPL_bp[:,PFISR_filtHeight]);
        #END FOR j
        
        
        #remove time ranges that aren't within the day range desired
        for j in range(0, PFISR_beamIDNum):
            tmp_cutOut = (PFISR_time_tmp[j] < dateRange_dayNum_full[0,1]*86400) | (PFISR_time_tmp[j] > dateRange_dayNum_full[-1,1]*86400); #get times out of bounds
            tmp_Keep = np.logical_not(tmp_cutOut); #inverse it
            
            if( np.sum(tmp_cutOut) > 0 ):
                PFISR_POPL_tmp[j] =  PFISR_POPL_tmp[j][tmp_Keep,:]; #preallocate list that holds a data
                PFISR_SNR_tmp[j] =  PFISR_SNR_tmp[j][tmp_Keep,:]; #preallocate list that holds a data
                #range isn't based on time
                PFISR_time_tmp[j] =  PFISR_time_tmp[j][tmp_Keep]; #preallocate list that holds a data
                PFISR_el_tmp[j] =  PFISR_el_tmp[j][tmp_Keep]; #preallocate list that holds a data
                PFISR_az_tmp[j] = PFISR_az_tmp[j][tmp_Keep]; #keep data that is in the range
                PFISR_vel_tmp[j] =  PFISR_vel_tmp[j][tmp_Keep,:]; #preallocate list that holds a data
                PFISR_dopplar_tmp[j] =  PFISR_dopplar_tmp[j][tmp_Keep,:]; #preallocate list that holds a data
                #height isn't based on time
                PFISR_POPL_bp_tmp[j] =  PFISR_POPL_bp_tmp[j][tmp_Keep,:]; #preallocate list that holds a data
                PFISR_SNR_bp_tmp[j] =  PFISR_SNR_bp_tmp[j][tmp_Keep,:]; #preallocate list that holds a data   
            #END IF
        #END FOR j
        
        
        # ============Combine into the main variable==============================
        
        if(k == 0):
            #write directly - hopefully python memory stuff doesn't ruin this (it did, .copy() fixes)
            PFISR_beamID = PFISR_beamID_tmp.copy();
            PFISR_POPL = PFISR_POPL_tmp.copy();
            PFISR_SNR = PFISR_SNR_tmp.copy();
            PFISR_range = PFISR_range_tmp.copy();
            PFISR_time = PFISR_time_tmp.copy();
            PFISR_el = PFISR_el_tmp.copy();
            PFISR_az = PFISR_az_tmp.copy();
            PFISR_vel = PFISR_vel_tmp.copy();
            PFISR_dopplar = PFISR_dopplar_tmp.copy();
            PFISR_height = PFISR_height_tmp.copy();
            PFISR_POPL_bp = PFISR_POPL_bp_tmp.copy();
            PFISR_SNR_bp = PFISR_SNR_bp_tmp.copy();
        else:
            #copy in
            for j in range(0, PFISR_beamIDNum):
    #            sys.exit()
    #            #calculate the median cause LISTS why can't numpy suck it up and make a nice cell? Maybe objects? but I don't feel like figuring that out now
    #            mainElvMedian = np.zeros(len(PFISR_el)); #preallocate
    #            for jj in range(0,len(PFISR_el)):
    #                mainElvMedian[jj] = np.median(PFISR_el[jj]); #get the median
    #            #END FOR jj
    #            mainBeamIDIndex = np.where( np.isclose(mainElvMedian,np.median(PFISR_el_tmp[j])) )[0]; #connect where the data should go
                
                mainBeamIDIndex = np.where( PFISR_beamID == PFISR_beamID_tmp[j] )[0]; #connect where the data should go
                if( mainBeamIDIndex.size == 1 ):
                    mainBeamIDIndex = mainBeamIDIndex.item(); #convert to scalar so it's easy
                    #make sure the ranges are the same
                    if( (np.all( np.isin(PFISR_range[mainBeamIDIndex] , PFISR_range_tmp[j]) ) == False) | (np.all( np.isin( PFISR_range_tmp[j] , PFISR_range[mainBeamIDIndex] ) ) == False) ):
                        #if false gotta figure out how to add some ranges
                        #update Main with required ranges
                        addToMain = np.where( np.isin(PFISR_range_tmp[j] , PFISR_range[mainBeamIDIndex] ) == False )[0];
                        if( addToMain.size > 0 ):
                            #updating main variables with new ranges
                            PFISR_range[mainBeamIDIndex] = np.concatenate(PFISR_range[mainBeamIDIndex],PFISR_range_tmp[j][addToMain], axis=0); #stack it on
                            PFISR_POPL[mainBeamIDIndex] = np.concatenate(PFISR_POPL[mainBeamIDIndex], np.zeros( (PFISR_time[mainBeamIDIndex].size,PFISR_range_tmp[j][addToMain].size) ), axis=1); #stack on 0's
                            PFISR_SNR[mainBeamIDIndex] = np.concatenate(PFISR_SNR[mainBeamIDIndex], np.zeros( (PFISR_time[mainBeamIDIndex].size,PFISR_range_tmp[j][addToMain].size) ), axis=1); #stack on 0's
                            PFISR_vel[mainBeamIDIndex] = np.concatenate(PFISR_vel[mainBeamIDIndex], np.zeros( (PFISR_time[mainBeamIDIndex].size,PFISR_range_tmp[j][addToMain].size) ), axis=1); #stack on 0's
                            PFISR_dopplar[mainBeamIDIndex] = np.concatenate(PFISR_dopplar[mainBeamIDIndex], np.zeros( (PFISR_time[mainBeamIDIndex].size,PFISR_range_tmp[j][addToMain].size) ), axis=1); #stack on 0's
                            PFISR_height[mainBeamIDIndex] = np.concatenate(PFISR_height[mainBeamIDIndex],PFISR_height_tmp[j][addToMain], axis=0); #stack it on
                            PFISR_POPL_bp[mainBeamIDIndex] = np.concatenate(PFISR_POPL_bp[mainBeamIDIndex], np.zeros( (PFISR_time[mainBeamIDIndex].size,PFISR_range_tmp[j][addToMain].size) ), axis=1); #stack on 0's
                            PFISR_SNR_bp[mainBeamIDIndex] = np.concatenate(PFISR_SNR_bp[mainBeamIDIndex], np.zeros( (PFISR_time[mainBeamIDIndex].size,PFISR_range_tmp[j][addToMain].size) ), axis=1); #stack on 0's
                            
                            mainRangeSortIndex = np.argsort(PFISR_range[mainBeamIDIndex],axis=0); #sort it
                            if( np.all( mainRangeSortIndex == np.arange(0,PFISR_range[mainBeamIDIndex].size) ) == False ):
                                #if above is false, then the ranges aren't sorted properly
                                PFISR_range[mainBeamIDIndex] = PFISR_range[mainBeamIDIndex][mainRangeSortIndex]; #sort correctly
                                PFISR_POPL[mainBeamIDIndex] = PFISR_POPL[mainBeamIDIndex][:,mainRangeSortIndex]; #sort correctly
                                PFISR_SNR[mainBeamIDIndex] = PFISR_SNR[mainBeamIDIndex][:,mainRangeSortIndex]; #sort correctly
                                PFISR_vel[mainBeamIDIndex] = PFISR_vel[mainBeamIDIndex][:,mainRangeSortIndex]; #sort correctly
                                PFISR_dopplar[mainBeamIDIndex] = PFISR_dopplar[mainBeamIDIndex][:,mainRangeSortIndex]; #sort correctly
                                PFISR_height[mainBeamIDIndex] = PFISR_height[mainBeamIDIndex][mainRangeSortIndex]; #sort correctly
                                PFISR_POPL_bp[mainBeamIDIndex] = PFISR_POPL_bp[mainBeamIDIndex][:,mainRangeSortIndex]; #sort correctly
                                PFISR_SNR_bp[mainBeamIDIndex] = PFISR_SNR_bp[mainBeamIDIndex][:,mainRangeSortIndex]; #sort correctly
                            #END IF
                        #END IF
                        addToTmp = np.where( np.isin(PFISR_range[mainBeamIDIndex] , PFISR_range_tmp[j]) == False )[0];
                        if( addToTmp.size > 0 ):
                            #update tmp with required ranges to be compatible with combining into main variables
                            PFISR_range_tmp[j] = np.concatenate( (PFISR_range_tmp[j], PFISR_range[mainBeamIDIndex][addToTmp]), axis=0); #stack it on
                            PFISR_POPL_tmp[j] = np.concatenate( (PFISR_POPL_tmp[j], np.zeros( (PFISR_time_tmp[j].size,PFISR_range[mainBeamIDIndex][addToTmp].size) ) ), axis=1); #stack on 0's
                            PFISR_SNR_tmp[j] = np.concatenate( (PFISR_SNR_tmp[j], np.zeros( (PFISR_time_tmp[j].size,PFISR_range[mainBeamIDIndex][addToTmp].size) ) ), axis=1); #stack on 0's
                            PFISR_vel_tmp[j] = np.concatenate( (PFISR_vel_tmp[j], np.zeros( (PFISR_time_tmp[j].size,PFISR_range[mainBeamIDIndex][addToTmp].size) ) ), axis=1); #stack on 0's
                            PFISR_dopplar_tmp[j] = np.concatenate( (PFISR_dopplar_tmp[j], np.zeros( (PFISR_time_tmp[j].size,PFISR_range[mainBeamIDIndex][addToTmp].size) ) ), axis=1); #stack on 0's
                            PFISR_height_tmp[j] = np.concatenate( (PFISR_height_tmp[j],PFISR_height[mainBeamIDIndex][addToTmp]), axis=0); #stack it on
                            PFISR_POPL_bp_tmp[j] = np.concatenate( (PFISR_POPL_bp_tmp[j], np.zeros( (PFISR_time_tmp[j].size,PFISR_range[mainBeamIDIndex][addToTmp].size) ) ), axis=1); #stack on 0's
                            PFISR_SNR_bp_tmp[j] = np.concatenate( (PFISR_SNR_bp_tmp[j], np.zeros( (PFISR_time_tmp[j].size,PFISR_range[mainBeamIDIndex][addToTmp].size) ) ), axis=1); #stack on 0's
                            
                            tmpRangeSortIndex = np.argsort(PFISR_range_tmp[j],axis=0); #sort it
                            if( np.all( tmpRangeSortIndex == np.arange(0,PFISR_range_tmp[j].size) ) == False ):
                                #if above is false, then the ranges aren't sorted properly
                                PFISR_range_tmp[j] = PFISR_range_tmp[j][tmpRangeSortIndex]; #sort correctly
                                PFISR_POPL_tmp[j] = PFISR_POPL_tmp[j][:,tmpRangeSortIndex]; #sort correctly
                                PFISR_SNR_tmp[j] = PFISR_SNR_tmp[j][:,tmpRangeSortIndex]; #sort correctly
                                PFISR_vel_tmp[j] = PFISR_vel_tmp[j][:,tmpRangeSortIndex]; #sort correctly
                                PFISR_dopplar_tmp[j] = PFISR_dopplar_tmp[j][:,tmpRangeSortIndex]; #sort correctly
                                PFISR_height_tmp[j] = PFISR_height_tmp[j][tmpRangeSortIndex]; #sort correctly
                                PFISR_POPL_bp_tmp[j] = PFISR_POPL_bp_tmp[j][:,tmpRangeSortIndex]; #sort correctly
                                PFISR_SNR_bp_tmp[j] = PFISR_SNR_bp_tmp[j][:,tmpRangeSortIndex]; #sort correctly
                            #END IF
                        #END IF
                    #END IF
                    #insert the data - the data will add to the TIME dimension. Above takes care of new ranges that pop up, mostly should be time expansion
                    PFISR_POPL[mainBeamIDIndex] = np.concatenate( (PFISR_POPL[mainBeamIDIndex] , PFISR_POPL_tmp[j]) , axis=0); #stack it on
                    PFISR_SNR[mainBeamIDIndex] = np.concatenate( (PFISR_SNR[mainBeamIDIndex] , PFISR_SNR_tmp[j]) , axis=0); #stack it on
                    #PFISR_range was updated above if it was needed to be updated
                    PFISR_time[mainBeamIDIndex] = np.concatenate( (PFISR_time[mainBeamIDIndex] , PFISR_time_tmp[j]) , axis=0); #stack it on
                    PFISR_el[mainBeamIDIndex] = np.concatenate( (PFISR_el[mainBeamIDIndex] , PFISR_el_tmp[j]) , axis=0); #stack it on
                    PFISR_az[mainBeamIDIndex] = np.concatenate( (PFISR_az[mainBeamIDIndex] , PFISR_az_tmp[j]) , axis=0); #stack it on
                    PFISR_vel[mainBeamIDIndex] = np.concatenate( (PFISR_vel[mainBeamIDIndex] , PFISR_vel_tmp[j]) , axis=0); #stack it on
                    PFISR_dopplar[mainBeamIDIndex] = np.concatenate( (PFISR_dopplar[mainBeamIDIndex] , PFISR_dopplar_tmp[j]) , axis=0); #stack it on
                    #PFISR_height was updated above if it was needed to be updated
                    PFISR_POPL_bp[mainBeamIDIndex] = np.concatenate( (PFISR_POPL_bp[mainBeamIDIndex] , PFISR_POPL_bp_tmp[j]) , axis=0); #stack it on
                    PFISR_SNR_bp[mainBeamIDIndex] = np.concatenate( (PFISR_SNR_bp[mainBeamIDIndex] , PFISR_SNR_bp_tmp[j]) , axis=0); #stack it on
                elif( mainBeamIDIndex.size > 1 ):
                    print("\n==============ERROR==============");
                    print("Data import for ISR Poker Flat has egregious inconsistencies and idk how they happened but they did. Check mainBeamIDIndex to start.");
                    import sys; #prep for destruction
                    sys.exit(); #rekt
                else:
                    #new beam ID section - gotta append!
                    PFISR_beamID.append(PFISR_beamID_tmp[j]);
                    PFISR_POPL.append(PFISR_POPL_tmp[j]);
                    PFISR_SNR.append(PFISR_SNR_tmp[j]);
                    PFISR_range.append(PFISR_range_tmp[j]);
                    PFISR_time.append(PFISR_time_tmp[j]);
                    PFISR_el.append(PFISR_el_tmp[j]);
                    PFISR_az.append(PFISR_az_tmp[j]);
                    PFISR_vel.append(PFISR_vel_tmp[j]);
                    PFISR_dopplar.append(PFISR_dopplar_tmp[j]);
                    PFISR_height.append(PFISR_height_tmp[j]);
                    PFISR_POPL_bp.append(PFISR_POPL_bp_tmp[j]);
                    PFISR_SNR_bp.append(PFISR_SNR_bp_tmp[j]);
                #END IF
            #END FOR j
        #END IF
        
        
    #END FOR k
    
    #reorder the elevations so they're in line
    mainElvMedian = np.zeros(len(PFISR_el)); #preallocate
    for jj in range(0,len(PFISR_el)):
        mainElvMedian[jj] = np.median(PFISR_el[jj]); #get the median
    #END FOR jj
    mainElvMedian_sortIndex = np.flip(np.argsort(mainElvMedian)); #get the indexes to sort it right - high to low
    mainElvMedian_sortIndex = np.ndarray.tolist(mainElvMedian_sortIndex);
    PFISR_beamID = [PFISR_beamID[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_POPL = [PFISR_POPL[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_SNR = [PFISR_SNR[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_range = [PFISR_range[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_time = [PFISR_time[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_el = [PFISR_el[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_az = [PFISR_az[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_vel = [PFISR_vel[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_dopplar = [PFISR_dopplar[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_height = [PFISR_height[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_POPL_bp = [PFISR_POPL_bp[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    PFISR_SNR_bp = [PFISR_SNR_bp[jj] for jj in mainElvMedian_sortIndex]; #sort it by elevation!
    
    ##quality check
    #for j in range(0,len(PFISR_el)):
    #    if( np.all( np.median(PFISR_el[j]) == PFISR_el[j]) == False ):
    #        print("\n==============ERROR==============");
    #        print("The median elevation angle "+str(np.median(PFISR_el[j]))+" is not consistent throughout PFISR_el for j#"+str(j));
    #        import sys; #prep for destruction
    #        sys.exit(); #rekt
    #    #END IF
    ##END FOR j
    
    return ISR_lat, ISR_long, PFISR_SNR, PFISR_SNR_bp, PFISR_POPL, PFISR_POPL_bp, PFISR_height, PFISR_time, PFISR_vel, PFISR_el, PFISR_az, PFISR_dopplar, PFISR_filtHeight