#Function to import ISR data in the Haystack Madgrial format (but already downloaded)
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
import os, glob
import h5py
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_monthWord_to_num import subfun_monthWord_to_num
from scipy import signal
from Code.subfun_strstr import strstr

#-----Testing variables-----
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
#from Code.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#FLG_ISRcodeType = 97; #sets the default ISR mode to alternating code
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units
#pointAltitude = 300; #km, altitude for centering upon
#bs = 2; #hr, high-pass cutoff period (any higher period is removed - 2 is usual)
#import matplotlib.pyplot as plt
#import os


def GRITI_import_ISR_Haystack(dateRange_dayNum_full,folder,dateRange_dayNum_zeroHr,pointAltitude,FLG_ISRcodeType = 97,filter_cutoffPeriod = 2*3600):

    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #get the full date range
    
    filenamePossibles = glob.glob( os.path.join(folder[1], "ISR", str(dateRange_dayNum_full[0,0]), "mlh*") ); #make a list of file names GLOB GLOB
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
            print("The current file: "+filenamePossibles[i]+" |did not follow supported mlh...hdf5 file formating. Skipping it.");
        #END TRYING
        testFile_expParams_start_Date = testFile_expParams_start[0:testFile_expParams_start.find(' ')]; #get the start date
        testFile_expParams_end_Date = testFile_expParams_end[0:testFile_expParams_end.find(' ')]; #get the start date
        
        testFile_expParams_start_Date_Split = np.int16(np.array(testFile_expParams_start_Date.split('-'))); #split by dash
        testFile_expParams_end_Date_Split = np.int16(np.array(testFile_expParams_end_Date.split('-'))); #split by dash
        
        filenameInRange[i] = (dateRange_full[0,0] >= testFile_expParams_start_Date_Split[0]) & (dateRange_full[0,1] >= testFile_expParams_start_Date_Split[1]) & \
        (dateRange_full[0,2] >= testFile_expParams_start_Date_Split[2]) & (dateRange_full[-1,0] >= testFile_expParams_end_Date_Split[0]) & \
        (dateRange_full[-1,1] >= testFile_expParams_end_Date_Split[1]) & (dateRange_full[-1,2] >= testFile_expParams_end_Date_Split[2]); #make sure data start is before or during the first day desired and the data end is after or during the last day desired
    #END FOR i
    
    if( np.sum(filenameInRange) > 1 ):
        print("\n==============~Warning~==============");
        print(str(np.sum(filenameInRange))+" files found that cover the date range required. The file names are:");
        filenameChosenIsBad = 1; #prep flag, assume bad
        while( filenameChosenIsBad == 1 ):
            print("{}".format( filenamePossibles[filenameInRange] ));
            print("Please choose which file to use using the numbers 1..."+str(np.sum(filenameInRange)));
            filenameChosen = input(); #get user input on this
            
            if( len(filenameChosen) > np.int64(np.floor(np.log10(np.abs(np.sum(filenameInRange)))) + 1) ): #check the maximum it can be, if the input is 3 characters long but the filenameInRange is 2 (1 character), then no go
                filenameChosenIsBad = 1; #flag on
            elif( filenameChosen.isdigit() == False ): #if it isn't all digits, no go
                filenameChosenIsBad = 1; #flag on
            elif( np.int64(filenameChosen) > np.sum(filenameInRange) ): #if the number is larger than the number range, no go
                filenameChosenIsBad = 1; #flag on
            else:
                filenameChosenIsBad = 0; #flag off, we're good
            #END IF
            
            if( filenameChosenIsBad == 1 ):
                print("\n==============~Warning~==============");
                print("Whatever you input, specifically '"+filenameChosen+"' isn't in the range needed. Try again :). The file names are:");
            #END IF
        #END WHILE
        filenameChosen = np.int64(filenameChosen)-1; #convert to int finally, minus 1 for indexing stuff
        filenamePath = filenamePossibles[filenameInRange][filenameChosen]; #pick out the correct file name finally
    elif( np.sum(filenameInRange) == 0 ):
        print("\n==============ERROR==============");
        print("Of the files found in the given directory '"+os.path.join(folder[1],"ISR")+"':");
        print({}.format( filenamePossibles ));
        print("None of them covered the date range needed: "+str(dateRange_full[0,:])+" to "+str(dateRange_full[-1,:]));
        print("Crashin' on purpose, sorry!");
        import sys; #prep for destruction
        sys.exit(); #rekt
    else:
        filenamePath = filenamePossibles[filenameInRange][0]; #choose the file that's the only one that fits
    #END IF
    
    #filename ='mlh130506g.002.hdf5'; #got this from madrigal site - claims to be final version
    
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
        #ZENITH FIRST
        keyToUse = np.zeros(len(filename_expParamsKeys),dtype='bool'); #preallocate this
        for i in range(0,len(filename_expParamsKeys)): #search each key for the stuff we want (32 is Zenith, code type, and pulse length (hard coded for now))
            keyToUse[i] = (strstr(filename_expParamsKeys[i],'32').size > 0) & (strstr(filename_expParamsKeys[i],str(FLG_ISRcodeType)).size > 0) & (strstr(filename_expParamsKeys[i],'00048').size > 0);
            #get the key to use, it'll be 1 if all these things are checked off
        #END FOR i
        keyToUse = np.where(keyToUse)[0][0]; #get the index (since 0 and 1 don't need extra stuff - the 1 will be chosen)
        Zenith_POPL = 10**(filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/popl'][:]); #Zenith POPL copied (alt to SNR - is electron density e-/m^3)
        try:
            Zenith_SNR = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/snp3'][:]; #Zenith SNR copied
        except:
            Zenith_SNR = np.zeros(Zenith_POPL.shape); #do zeros instead because the file lacks SNP3 which is common
        #END TRY
#        Zenith_POPL = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/popl'][:]; #Zenith POPL copied (alt to SNR - is electron density log10(e-/m^3))
        Zenith_range = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/range'][:]; #Zenith range, km
        Zenith_time = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/timestamps'][:]; #Zenith time, sec from 1970
        Zenith_el = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/1D Parameters/el1'][:]; #Zenith elevation angle, deg
        Zenith_az = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/1D Parameters/az1'][:]; #Zenith azimuth angle, deg
        Zenith_vel = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/vo'][:]; #Zenith velocity, km/s
        Zenith_dopplar = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/vdopp'][:]; #Zenith dopplar eh?
        

        
        #MISA NOW
        keyToUse = np.zeros(len(filename_expParamsKeys),dtype='bool'); #preallocate this
        for i in range(0,len(filename_expParamsKeys)): #search each key for the stuff we want (32 is Zenith, code type, and pulse length (hard coded for now))
            keyToUse[i] = (strstr(filename_expParamsKeys[i],'31').size > 0) & (strstr(filename_expParamsKeys[i],str(FLG_ISRcodeType)).size > 0) & (strstr(filename_expParamsKeys[i],'00048').size > 0);
            #get the key to use, it'll be 1 if all these things are checked off
        #END FOR i
        keyToUse = np.where(keyToUse)[0][0]; #get the index (since 0 and 1 don't need extra stuff - the 1 will be chosen)
        MISA_POPL = 10**(filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/popl'][:]); #MISA POPL copied (alt to SNR - is electron density e-/m^3)
        try:
            MISA_SNR = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/snp3'][:]; #MISA SNR copied
        except:
            MISA_SNR = np.zeros(MISA_POPL.shape); #do zeros instead because the file lacks SNP3 which is common
        #END TRY
#        MISA_POPL = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/popl'][:]; #MISA POPL copied (alt to SNR - is electron density log10(e-/m^3))
        MISA_range = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/range'][:]; #MISA range, km
        MISA_time = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/timestamps'][:]; # MISA time, sec from 1970
        MISA_el = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/1D Parameters/el1'][:]; #MISA elevation angle, deg
        MISA_az = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/1D Parameters/az1'][:]; #MISA azimuth angle, deg
        MISA_vel = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/vo'][:]; #MISA velocity, km/s
        MISA_dopplar = filename['/Data/Array Layout/'+filename_expParamsKeys[keyToUse]+'/2D Parameters/vdopp'][:]; #MISA dopplar eh?
        
        for i in range(0,Zenith_time.size): #range-squared correction for POPL
            Zenith_POPL[:,i] = Zenith_POPL[:,i]/Zenith_range**2;
        #END IF
        for i in range(0,MISA_time.size): #range-squared correction for POPL
            MISA_POPL[:,i] = MISA_POPL[:,i]/MISA_range**2;
        #END IF
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
        
        with open( os.path.join(folder[1],"Supporting","leap-seconds.txt"),"r") as leapFile: #prep to read 
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
        with open( os.path.join(folder[1],"Supporting","leap-seconds.txt"),"w") as leapFile: #prep to write 
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
    
    # Zenith_time = (Zenith_time-secSince)/(3600*24) + dateRange_dayNum_zeroHr[1]; #days, removes offset from 1970, converts to days and then adjusts to the zero hour day
    # MISA_time = (MISA_time-secSince)/(3600*24) + dateRange_dayNum_zeroHr[1]; #days, removes offset from 1970, converts to days and then adjusts to the zero hour day
    #seconds here
    Zenith_time = Zenith_time-secSince + dateRange_dayNum_zeroHr[1]*86400; #sec, removes offset from 1970, converts to days and then adjusts to the zero hour day
    MISA_time = MISA_time-secSince + dateRange_dayNum_zeroHr[1]*86400; #sec, removes offset from 1970, converts to days and then adjusts to the zero hour day
    
    # ========================Replacing NAN Values with 0======================
    Zenith_SNR[np.isnan(Zenith_SNR)] = 0; #eplaces NAN with 0
    
    MISA_SNR[np.isnan(MISA_SNR)] = 0; #Replaces NAN with 0
    
    Zenith_POPL[np.isnan(Zenith_POPL)] = 0; #Replaces NAN with 0
    
    MISA_POPL[np.isnan(MISA_POPL)] = 0; #Replaces NAN with 0
    
    Zenith_range[np.isnan(Zenith_range)] = 0; #Replaces NAN with 0
    
    MISA_range[np.isnan(MISA_range)] = 0; #Replaces NAN with 0
    
    Zenith_vel[np.isnan(Zenith_vel)] = 0; #Replaces NAN with 0
    
    MISA_vel[np.isnan(MISA_vel)] = 0; #Replaces NAN with 0
    
    Zenith_dopplar[np.isnan(Zenith_dopplar)] = 0; #Replaces NAN with 0
    
    MISA_dopplar[np.isnan(MISA_dopplar)] = 0; #Replaces NAN with 0
    
    
    # ============Convert Time to EDT (UT-4 hrs)==============================
    # Zenith_time = Zenith_time - 4/24;
    # MISA_time = MISA_time - 4/24;
    # Disabled currently
    
    # ===============Zenith and MISA log of SNR Calcs=========================
    Zenith_SNR[Zenith_SNR <= 0] = 0.1; #Removes negative numbers/0's and replaces with 0.1
    #Zenith_SNR_log = np.log10(Zenith_SNR);
    
    MISA_SNR[MISA_SNR <= 0] = 0.1; #Removes negative numbers/0's and replaces with 0.1
    #MISA_SNR_log = np.log10(MISA_SNR);
    
    # ==============Zenith and MISA Height Calcs (from range + angle)=========
    if( np.abs(np.median(Zenith_el)-np.mean(Zenith_el)) < 0.1 ):
        Zenith_height = Zenith_range/((1. + ((np.tan( (90. - np.median(Zenith_el))*np.pi/180 ))**2))**0.5);
    else:
        print("WARNING: For file "+filenamePath+"\nZenith elevation mean is "+str(np.mean(Zenith_el))+" but median is "+str(np.median(Zenith_el))+" which indicates Zenith is moving but it shouldn't be able to - using median");
        Zenith_height = Zenith_range/((1. + ((np.tan( (90. - np.median(Zenith_el))*np.pi/180 ))**2))**0.5);
    #END IF
    if( np.abs(np.median(Zenith_el)-np.mean(Zenith_el)) < 0.1 ):
        MISA_height = MISA_range/((1. + ((np.tan( (90. - np.median(MISA_el))*np.pi/180 ))**2))**0.5);
    else:
        print("WARNING: For file "+filenamePath+"\MISA elevation mean is "+str(np.mean(MISA_el))+" but median is "+str(np.median(MISA_el))+" which indicates MISA is moving with time, may mess with data - using median b/c I don't have code for that");
        MISA_height = MISA_range/((1. + ((np.tan( (90. - np.median(MISA_el))*np.pi/180 ))**2))**0.5);
    #END IF
    
    # ===============High Pass Filtering to Excentuate Hourly Periods=========
    #These are unused, but accentuate the periods desired
    #Original file used highpass_fir.m, it has been merged
    # lp=1; # Lower period (in hrs)
    # hp=2; # Higher period (in hrs)
    # lf=(1/hp);  # Lowpass frequency corner (1/hr)
    # hf=(1/lp);  # Highpass frequency corner (1/hr)
    
    
    bs = (1/filter_cutoffPeriod);   # highpass cuttoff frequency (not sure what to make of it)
    #I think it is related to lf above (which is 1/2 )
    #ls = (1/(1/2)); #low-pass cutoff frequency
    
#    Zenith_delt = np.mean(Zenith_time[1:] - Zenith_time[0:-1]); #Calculates a time delta based on avg of all the deltas
    Zenith_delt = np.median(np.diff(Zenith_time)); #Calculates a time delta based on avg of all the deltas

    # ===============Highpass filtering on original Zenith SNR================
    
    n=42; # order of the Hamming window used in the custom function (uses 43 in reality)
    # c = 3.32*pi; #Hamming constant
    # M = n/2;
    # bandwidth = c/M;
    #The above was included to investigate the bandwidth - 1/2. Related to
    #bs/lf?
    
    fp = bs; # stop band stoppin freq (or pass band passing freq, depending on how you look at it)
    
    f= 1/(Zenith_delt); #1/sec, the sampling frequency, based off of the time delta calc'd
    
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
        Zenith_SNR_bp = signal.filtfilt(b,a,Zenith_SNR,axis=1,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
        #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    except ValueError:
        Zenith_SNR_bp = np.ones(Zenith_SNR.shape)*np.nan; #fill with NaN's because the Zenith was turned on/off and it's time frame is actually too short to high-pass correctly
    #END TRY
    try:
        Zenith_POPL_bp = signal.filtfilt(b,a,Zenith_POPL,axis=1,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
        #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    except ValueError:
        Zenith_POPL_bp = np.ones(Zenith_POPL.shape)*np.nan; #fill with NaN's because the Zenith was turned on/off and it's time frame is actually too short to high-pass correctly
    #END TRY
    # Zenith_SNR_bp = Zenith_SNR; %No filter option to see what it really looks like
    try:
        Zenith_vel_bp = signal.filtfilt(b,a,Zenith_vel,axis=1,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
        #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    except ValueError:
        Zenith_vel_bp = np.ones(Zenith_vel.shape)*np.nan; #fill with NaN's because the Zenith was turned on/off and it's time frame is actually too short to high-pass correctly
    #END TRY
    
    
    MISA_delt = np.mean(MISA_time[1:] - MISA_time[0:-1]);
    
    # ===============Highpass filtering on original MISA SNR==================
    n=42; # order of the Hamming window;
    
    if(MISA_delt != Zenith_delt): #If they're not the same, it re-calcs the window
        fp = bs; # stop band stoppin freq (or pass band passing freq, depending on how you look at it)
    
        f = 1/(MISA_delt); #1/hr, the sampling frequency, based off of the time delta calc'd
        #Should be same as Zenith
    
        wp = 2*fp/f; # Normalizing the frequencies
    
        #Calculation of filter coefficients
        #[b,a]=fir1(n,wp,'high'); #This takes the order and lower cut off frequency
        b = signal.firwin(n+1, wp, window='hann', pass_zero=False); #just high-pass (2 hr and lower OK)
        #Uses the default hamming window
        #Should be the same as Zenith
    #END IF
    
    try:
        MISA_SNR_bp = signal.filtfilt(b,a,MISA_SNR,axis=1,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
        #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    except ValueError:
        MISA_SNR_bp = np.ones(MISA_SNR.shape)*np.nan; #fill with NaN's because the MISA was turned on/off and it's time frame is actually too short to high-pass correctly
    #END TRY
    try:
       MISA_POPL_bp = signal.filtfilt(b,a,MISA_POPL,axis=1,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
       #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    except ValueError:
        MISA_POPL_bp = np.ones(MISA_POPL.shape)*np.nan; #fill with NaN's because the MISA was turned on/off and it's time frame is actually too short to high-pass correctly
    #END TRY
    try:
       MISA_vel_bp = signal.filtfilt(b,a,MISA_vel,axis=1,padtype = 'odd', padlen=3*(b.size-1)); #Appies the filter
       #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    except ValueError:
        MISA_vel_bp = np.ones(MISA_vel.shape)*np.nan; #fill with NaN's because the MISA was turned on/off and it's time frame is actually too short to high-pass correctly
    #END TRY
            
    
    # =============== Spectral Estimation of Filtered Zenith SNR at 300 km ===
    # 47 corresponds to 300 km
    # 25 corresponds to 200 km
    # 3 corresponds to 100 km
    # 70 corresponds to 400 km
    # 92 corresponds to 500 km (all approximate)
    Zenith_filtHeight = np.where( np.min(np.abs(Zenith_height - pointAltitude)) == np.abs(Zenith_height - pointAltitude) )[0][0]; #See 'table' above (300 km was original)
    
    # Removing the mean from Zenith SNR time series at 300 km....
    Zenith_SNR_bp[:,Zenith_filtHeight] = Zenith_SNR_bp[:,Zenith_filtHeight] - np.mean(Zenith_SNR_bp[:,Zenith_filtHeight]);
#    Zenith_POPL_bp[:,Zenith_filtHeight] = Zenith_POPL_bp[:,Zenith_filtHeight] - np.mean(Zenith_POPL_bp[:,Zenith_filtHeight]);
    
    
    # =============== Spectral Estimation of Filtered MISA SNR at 330 km =====
    #3330 km MISA approximates to 300 km Zenith. Looking at height, very rough
    #guesstimation
    #51 corresponds to 330 km
    
    MISA_filtHeight = np.where( np.min(np.abs(MISA_height - pointAltitude)) == np.abs(MISA_height - pointAltitude) )[0][0]; #range at which filtering for(see 'table' above-330 orig)
    
    # Removing the mean from MISA SNR time series at 330 km....
    MISA_SNR_bp[:,MISA_filtHeight] = MISA_SNR_bp[:,MISA_filtHeight] - np.mean(MISA_SNR_bp[:,MISA_filtHeight]);
#    MISA_POPL_bp[:,MISA_filtHeight] = MISA_POPL_bp[:,MISA_filtHeight] - np.mean(MISA_POPL_bp[:,MISA_filtHeight]);


    return ISR_lat, ISR_long, Zenith_SNR, Zenith_SNR_bp, Zenith_POPL, Zenith_POPL_bp, Zenith_height, Zenith_time, Zenith_vel, Zenith_vel_bp, Zenith_el, Zenith_az, Zenith_dopplar, Zenith_filtHeight, MISA_SNR, MISA_SNR_bp, MISA_POPL, MISA_POPL_bp, MISA_height, MISA_time, MISA_vel, MISA_vel_bp, MISA_el, MISA_az, MISA_dopplar, MISA_filtHeight
#END DEF












