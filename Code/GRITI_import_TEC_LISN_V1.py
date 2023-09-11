#Function to import TEC data from the LISN FTP network
#RD on 8/23/2018
#
#

#The date range is heavily limited
#The files are labeled as 4-letterYYDDD.dat
#
#The 1st line includes the name of the city and the country, where the receiver is installed
#The 2nd line should indicate the latitude, longitude, and altitude (km). 
#The 3rd line has the brand of the receiver: Novatel or Ashtech. Many times the brand is mistaken, but this does not have any effect on the processing.
#The 4th line contains the type of the input file.  The program that creates the files only uses RINEX files and it is typed: 1dRINEX that means a 1 day RINEX.
#The 5th line has the year, month and day.  This double checks the information of the file name.
#The following lines list the PRN # and estimated # of TEC data entries for that PRN (but it is almost always incorrect due to farther data processing).
# That line is followed by the TEC information for each period (10, 15, or 30 seconds).  The PRNs (satellites) are listed sequentially from 01 to 32.
#The columns in the TEC info are:
#1 - Consecutive numbers
#2 – day number from January 1, 2000 = 1
#3, 4, 5 – hour, min, sec (UT)
#6 – vTEC, see below for calc
#7 – local time
#8 – elevation of looking direction
#9 – latitude of sub-ionospheric point
#10 – Longitude
#11 – slant TEC
#12 – Azimuth

#6 - vTEC is calculated via:
#sf   = 1.0/cos(asin(0.94092*cos(xelev(nk)*degrad)))
#eqvtec(iprn,m) = tec(nk)/sf
#The factor = 0.94092 includes the altitude of 350.0 km (which is the sub-ionospheric point)



import numpy as np
from scipy.signal import savgol_filter
import os
import time
import sys
from urllib.request import urlopen, urlretrieve
import html2text
from subfun_findall import strstr
from subfun_strfind import strfind
#from subfun_downloadProgress import downloadProgress
import h5py
#from subfun_date_to_dayNum import subfun_date_to_dayNum
from subfun_addADay import subfun_addADay
from subfun_dayNum_to_date import subfun_dayNum_to_date
from scipy import signal
#from numba import jit, prange
##-----Testing variables-----
##http://cedar.openmadrigal.org/ftp/fullname/Ross+Dinsmore/email/rld5204@psu.edu/affiliation/PSU/kinst/8000/year/2015/kindat/3505/format/hdf5/
#import matplotlib.pyplot as plt
#from matplotlib import _pylab_helpers
#from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#from subfun_date_to_dayNum import subfun_date_to_dayNum
#
##Date range goes Month-Day-Year
#dateRange = np.array([[2013,5,6],[2013,5,8]],dtype="int16"); #dates are in int16 because they can be
##dateRange = np.array([[2014,7,30],[2014,8,1]],dtype="int16"); #dates are in int16 because they can be
##dateRange = np.array([[2016,11,28],[2016,11,30]],dtype="int16"); #dates are in int16 because they can be
##dateRange = np.array([[2014,12,31],[2015,1,1]],dtype="int16"); #for debug, check year success
##dates better go earlier -> later
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_full_orig = subfun_dayNum_to_date(dateRange_dayNum_full); #record orginal date ranges
##print("{}".format(dateRange))
#folder = [os.getcwd()]; #current working directory, leave it as this call usually
#folder.append('E:\Big Data'); #place to save data files to
##folder var structure: 0 = running folder, 1 = data folder
#FLG_reqPaddedDays = 0; #0 no padded day required, 1 required padded days (better science since filtering on day edges will have data cutoff)
#FLG_deleteOrig = 0; #0 don't delete original raw TEC data file, 1 delete original raw TEC data file
#FLG_deleteUnfilt = 0; #0 don't delete unfiltered TEC data file, 1 delete unfiltered TEC data file
#FLG_overwrite = 0; #0 don't overwrite data even if it is there, 1 overwrite final data even if it is there, 2 overwrite previously downloaded data as well (site files, converted site files)
#FLG_verboseReturn = 0; #0 return just the data arrays for int/float/string, 1 return extra information in the HDF5 file such as the filtering info and pierce point altitude
#FLG_dataAggregation = 0; #0 don't aggregate nearby TEC values with a site/sat combo, 1 aggregate nearby TEC values with a site/sat combo for MORE DATA PER DATA
#minElevation = 30; #deg, min elevation angle accepted
#minimumTimeGap = 5; #min, how long to accept a gap in a satellite observation before calling it a new, unrelated observation
#TEC_dataRate = 30; #sec/datapt, the data rate. No super easy way to get it out of the data automatically, so it's a variable to control
#filter_savGolPeriod = 1*60; #min, filter period
#order_savGol = 1; #order to use
#filter_cutoffPeriod = 2; #hr, cut off period to high pass filter at (high pass removes anything higher period than # (e.g. 2) hr)
##inspired by https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
#deltaTEC_compareValue = 6; #hard coded, seems like a good number and such
#web_base_name = 'Ross+Dinsmore'; #put your name here
#web_base_email = 'rld5204@psu.edu'; #put your email here
#web_base_affil = 'PSU'; #put your affiliation here (school etc)
#plotLatRange = [35,50]; #latitude limit for plotting
##-90 to 90 is world, 35 to 50 is good for USA East Coast
#plotLongRange = [-85,-60]; #longitude limit for plotting
##-180 to 180 is world, -85 to -60 is good for USA East Coast
#TEC_dataAgg_timeAdditionLimit = 1*60; #minutes, time to add on to a satellite's visible time before calling it quits
#TEC_dataAgg_distToPts = 50; #km, distance to extended points for TEC data to exist
#TEC_dataLimPercent = 0.75; #0.75 = 75%, means 75% of the data must be a good fit (within deltaTEC_compareValue)
#TEC_deltaLim = 0.5; #unrelated to above, value between fit/real data to find where ends of the fit where it starts diverging


def GRITI_import_TEC_LISN(dateRange_dayNum_full, folder,plotLatRange,plotLongRange,web_LISN_creden_user,web_LISN_creden_pass, minElevation = 30, minimumTimeGap = 5, deltaTEC_compareValue = 3.5, TEC_dataRate = 30, filter_savGolPeriod = 1*60, order_savGol=1, filter_cutoffPeriod = 2, FLG_reqPaddedDays = 1, FLG_deleteOrig = 0, FLG_deleteUnfilt = 0, FLG_overwrite = 0, FLG_verboseReturn = 0 ,FLG_dataAggregation = 0 , TEC_dataAgg_timeAdditionLimit = 60, TEC_dataAgg_distToPts = 20, TEC_dataLimPercent = 0.05 , TEC_deltaLim = 0.5, TEC_timeTolerance = 0.1, TEC_maxAmpAllowed = 6 ):
    #==============Constants Needed==============
    version_unfilt = 1.1; #unfiltered algorithm version
    #changelog:
    #1 10/8/2019 - initial algorithm
    #1.1 9/14/2020 - fixed hour 24/minute 60 handling (before ignored, now adjusted to +1 day/+1 hour)
    
    version_filt = 1.3; #filtered algorithm version
    #1 10/8/2019 - initial algorithm
    #1.1 9/11/2020 - fixed non-0/30 time step handling (29/59 were big ones)
    #1.2 9/16/2020 - removed highpass filter, only dampens data
    #1.3 9/17/2020 - outlier control introduced
    
    folder_TEC = 'TEC'; #name for the TEC folder
    folder_fileEnding = '_LISN.h5'; #file extension for hdf5 files
    folder_TEC_rawData = 'LISN Raw'; #sub-folder that holds the raw LISN data
    web_base_site = 'ftp://'+str(web_LISN_creden_user)+':'+str(web_LISN_creden_pass)+'@'+'lisn.igp.gob.pe/'; #website to be used
    earthRadius = 6371; #km, earth radius
    TEC_dataAgg_distToPts_degc = TEC_dataAgg_distToPts/earthRadius*180/np.pi; #degc, distance to extended points for TEC data to exist
    TEC_dataAgg_distToPts_degcSq = TEC_dataAgg_distToPts_degc**2; #degc^2, square it to avoid sqrt later
    
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
    
    #==============Unfiltered File Layout==============
    #Integer Layout
    #0 = Satellite ID [# that corresponds to GPS sat]
    locIntUnfilt_sat = 0; #index where sat ID is
    #1 = Year timestamp [years]
    locIntUnfilt_year = 1; #index where year timestamp is
    #2 = Day Number timestamp [days]
    locIntUnfilt_dayNum = 2; #index where day number is
    #3 = Hour timestamp [hrs]
    locIntUnfilt_hour = 3; #index where hour timestamp is
    #4 = Minute timestamp [mins]
    locIntUnfilt_min = 4; #index where minute timestamp is
    #5 = Second timestamp [secs]
    locIntUnfilt_sec = 5; #index where second timestamp is
    locIntUnfilt_size = 6; #size of the int variable
    
    #Float Layout
    #0 = current time in day format [days] - does not support years
    locFloatUnfilt_time = 0; #index where time in days is
    #1 = geodedic latitude [arcdeg]
    locFloatUnfilt_lat = 1; #index where geodedic latitude is
    #2 = longitude [arcdeg]
    locFloatUnfilt_long = 2; #index where longitude is
    #3 = elevation [deg]
    locFloatUnfilt_elev = 3; #index where elevation is
    #4 = line-of-sight TEC [TECU]
    locFloatUnfilt_sTEC = 4; #index where los TEC is
    #5 = error in line-of-sight TEC [TECU]
    locFloatUnfilt_sTECerr = 5; #index where error in los TEC is
    #6 = vertical TEC [TECU] (which is the sTEC combined with a mapping function)
    locFloatUnfilt_vTEC = 6; #index where vertical TEC is
    locFloatUnfilt_size = 7; #size of float variable
    
    #String Layout
    #[] = Receiver Site Name (there's no dim on this one) - use .decode('UTF-8') to make it a string again
    locStringUnfilt_site = 0; #index where site name is
    locStringUnfilt_size = 1; #size of string layout
    
    #==============Adjust dates to be padded==============
    print("\nDate range requested (yr/day num format): {}/{} to {}/{}.".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
    dateRange_dayNum_full_orig = dateRange_dayNum_full; #record orginal date ranges
    dateRange_dayNum_full = subfun_addADay(dateRange_dayNum_full); #call fun to pad a day onto the ends of the range
    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #convert the padded range to date range
    print("Date range used due to padding requirement for filtering (yr/day num format): {}/{} to {}/{}.\n".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
    #The padding is required for filtering since the two days 
    
    print("==============Importing TEC LISN Func - Starting==============");
    #==============Check if data already there==============
    if( os.path.isdir(folder[1] + '\\' + folder_TEC) == 0 ): #check if TEC folder exists
        #if not, make it
        os.makedirs(folder[1] + '\\' + folder_TEC);
        print("NOTA BENE: Importing TEC Madrigal Func - Created TEC directory: {}\n".format(folder[1] + '\\' + folder_TEC) );
    #END IF
        
    dateRange_uniqueYears = np.unique(dateRange_dayNum_full[:,0]); #get the unique years involved
    TEC_dataAmnt = len(dateRange_dayNum_full[:,0]); #get number of days needed to investigate
    TEC_dataAvail = np.zeros( (TEC_dataAmnt,) , dtype='int8');  #preallocate logical array
    #-3 = data is a padded day and it has been passed in the loop, signifying it is donezo
    #-2 = data is a padded day and the non-padded days (the desired days) have been filtered already
    #-1 = no data available on required days, will quit
    #0 = no data available on padded day and flag says its OK to skip. Will impact data availability on edge of days, but middle of day is OK
    #1 = note data is there, filtered and all
    #2 = data exists and needs to be downloaded - OR wasn't downloaded fully and needs to be redownloaded, end goal is the same
    #3 = note data is already downloaded but needs to be converted to the standardized format & naming scheme
    #4 = note that data is already downloaded and converted, but not filtered
    #5 = orig data needs to be downloaded (2) - BUT filtered data is finished so don't filter this day (for supporting other days to be filtered)
    #6 = orig data is downloaded but not converted (3) - BUT filtered data is finished so don't filter this day (for supporting other days to be filtered)
    #7 = unfiltered data is already downloaded and converted (4) - BUT filtered data is finished so don't filter this day (for supporting other days to be filtered)
    
    TEC_dataPath = ["{}\{}\{}".format(a_, b_, c_) for a_, b_, c_ in zip([folder[1]]*TEC_dataAmnt, [folder_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)) ) ]; #get the base path where data will be in
    TEC_dataFileName = ["{}_{}_{}{}".format(a_, b_, c_, d_) for a_, b_, c_, d_ in zip([folder_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)),np.ndarray.tolist(dateRange_dayNum_full[:,1].astype(str)),[folder_fileEnding]*TEC_dataAmnt ) ]; #get the expected filenames
    TEC_dataFileNameUnfilt = ["{}_{}_{}_unfilt{}".format(a_, b_, c_, d_) for a_, b_, c_, d_ in zip([folder_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)),np.ndarray.tolist(dateRange_dayNum_full[:,1].astype(str)),[folder_fileEnding]*TEC_dataAmnt ) ]; #get the expected filenames for unfiltered data (if stopped mid filtering)
    TEC_dataFilePath = ["{}\{}".format(a_, b_) for a_, b_ in zip(TEC_dataPath,TEC_dataFileName ) ]; #get the full path right to the expected files
    TEC_dataFilePathUnfilt = ["{}\{}".format(a_, b_) for a_, b_ in zip(TEC_dataPath,TEC_dataFileNameUnfilt ) ]; #get the full path right to the expected files that are unfiltered
    
    for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
        if( os.path.isdir(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) );
            print("NOTA BENE: Importing TEC Madrigal Func - Created TEC subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) ));
        #END IF
        if( os.path.isdir(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) + '\\' + folder_TEC_rawData ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) + '\\' + folder_TEC_rawData );
            print("NOTA BENE: Importing TEC Madrigal Func - Created LISN raw data TEC subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) + '\\' + folder_TEC_rawData ));
        #END IF
    #END FOR
        
    #==============Check website and locally for data (in all formats)==============
    #    FLG_newYr = -1; #start the new year flag - used to detect if a new year occured
    for i in range(0,len(dateRange_dayNum_full[:,0])): #loop to check if any data exists - data will not be downloaded yet!
        #Since data download takes so long (5GB+) I decided to check if it is there before trying to download anything
        if( (os.path.isfile(TEC_dataFilePath[i]) == 0) ):
            #If data doesn't exist, time to get it
            
            if( (os.path.isfile(TEC_dataFilePathUnfilt[i]) == 0) ):
                #if the unconverted file is not there, download it
                #first - need to check if TEC data is availiable for the days requested
                
    #                if( FLG_newYr != dateRange_full[i,0] ): #check to go get the year's web links and file names and file dates
                    #basically, if within the same year only need to get this once for the whole year's data availiability
    #                    FLG_newYr = dateRange_full[i,0]; #set it so the flag works good
                #above disabled since currently the ftp server is just a direct zone
                
                try: #only need to download once, so checking if it exists
                    rendered_content; #try to call it
                except NameError:
                    page = urlopen(web_base_site); #get raw HTML
                    html_content = page.read(); #read off the HTML from whatever the page holder is
                    charset = page.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                    if( charset is None ):
                        charset = 'utf-8'; #assume utf-8
                    #END IF
                    html_content = html_content.decode(charset); #"decode" the HTML content so it's legible
                    rendered_content = html2text.html2text(html_content); #render the HTML like webpage would and get the real stuff
                    #print("{}".format(rendered_content)); #print for debug
                #END TRY
                
                web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for
                web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
                
                if( web_fileNamesIndex.size > 0 ):
                    TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
                    for j in range(0,web_fileNamesIndex.size):
                        TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                        if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                            TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                        #END IF
                    #END FOR j
                else:
                    #if online doesn't have the data - maybe local still does
                    web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
                    fileNames = os.listdir(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData ); #get the names of the files in the folder
                    
                    fileNames_GoodMask = np.zeros( len(fileNames) ); #a whitelist of good file names
                    for j in range(0,len(fileNames) ):
                        fileNames_Index = strstr(fileNames[j],web_fileNames_dateString); #find if the string is in the file name (if empty, it's not)
                        
                        if( fileNames_Index.size != 0 ):
                            #it's a useful file as its name matches the expected date string
                            fileNames_GoodMask[j] = 1; #note it's a good file
                        #END IF
                    #END FOR j
                    if( np.sum(fileNames_GoodMask) > 0 ):
                        TEC_fileNameOrig_isThere = 1; #set the flag to 1 if a file is local
                    else:
                        TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local
                    #END IF                    
                #END IF
                                
                #Check if the data for the date requested is locally available and up to date
                if( TEC_fileNameOrig_isThere == 1 ):
                    TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                    #use case if download success but failure on conversion - big enough files to warrant this care :)
                    
                #Check if website has data for the date requested
                elif( web_fileNamesIndex.size > 0 ): #basically checks the date vector with all the dates on the website, if one matches we're good
                    TEC_dataAvail[i] = 2; #note data is there but needs to be downloaded
                else:
                    if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == 1 ): 
                        #Catch where no data available
                        print("\n==============ERROR==============");
                        print("There is no data availiable on {}/{}/{} in YR/M/D format. [NOT SUPPORTED IN LISN] Printing availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
    #                        print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                        TEC_dataAvail[i] = -1; #note the error
                    else:
                        #Catch where no data available - but it is a padded day to help with filtering on day edges
                        if( FLG_reqPaddedDays == 0 ):
                            print("\n==============~Warning~==============");
                            print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are NOT required via passed flag - program will not quit but filtering on day edges will be impacted.\n[NOT SUPPORTED IN LISN] Printing availiable data days for the relevant year:".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
    #                            print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                            TEC_dataAvail[i] = 0; #note the error
                        else:
                            print("\n==============ERROR==============");
                            print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are required via passed flag.\n[NOT SUPPORTED IN LISN] Printing availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
    #                            print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                            TEC_dataAvail[i] = -1; #note the error
                        #END IF
                    #END IF
                #END IF
            elif( os.path.isfile(TEC_dataFilePathUnfilt[i]) == 1 ): #check if unfiltered data file exists            
                try: #gonna try to read the file - if we fail, it failed mid download probably
                    with h5py.File(TEC_dataFilePathUnfilt[i], 'r') as testFile:
                        testFile.keys(); #tries to check out some stuff in the file
                    #END WITH
                    with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'r') as TEC_fileUnfilt:
                        TEC_unfilt_version = TEC_fileUnfilt["float"].attrs['version']; #read the version
                    #END WITH
                    if( TEC_unfilt_version >= version_unfilt ):
                        #if the unfiltered file is there (filtering got stopped mid-filter or something)
                        TEC_dataAvail[i] = 4; #note that data is already downloaded and converted, but not filtered
                        #use case if download and conversion success but failure on filter - big enough files to warrant this care :)
                    else:
                        print("\n==============~Warning~==============");
                        print("The file {} uses version {} while version {} is the most current version. The file will be deleted and remade with the new version of the algorithm.\n".format(TEC_dataFilePathUnfilt[i], TEC_unfilt_version , version_unfilt) );
                        #force an OSError to use the rest of the algorithm - uses ***'s which can't be in file names
                        with h5py.File(TEC_dataFilePathUnfilt[i]+'.***notapossiblefile', 'r') as testFile:
                            testFile.keys(); #tries to check out some stuff in the file
                        #END WITH
                    #END IF
                except OSError:
                    print("\n==============~Warning~==============");
                    print("Data wasn't converted to standardized format & naming scheme for {} and will be reconverted. Deleting partial file as well.\n".format(TEC_dataFilePathUnfilt[i]) );
                    os.remove(TEC_dataFilePathUnfilt[i]); #delete partially downloaded file to avoid conflicts
                    
                    web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for
                    web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
                    
                    TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
                    for j in range(0,web_fileNamesIndex.size):
                        TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                        if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                            TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                        #END IF
                    #END FOR j
                    
                    if( TEC_fileNameOrig_isThere == 1 ): #if data is there, yay
                        TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                        #use case if download success but failure on conversion - big enough files to warrant this care :)
                    else: #if not gotta download
                        print("Data wasn't downloaded fully for {}/{}/{} (YR/M/D) and will be redownloaded.\n".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
                        TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
                    #END IF
                #END TRYING
            else: #otherwise orig data downloaded but not converted to standard (fast) format
                web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for
                web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
            
                TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
                for j in range(0,web_fileNamesIndex.size):
                    TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                    if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                        TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                    #END IF
                #END FOR j
                
                if( TEC_fileNameOrig_isThere == 1 ): #if data is there, yay
                    TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                    #use case if download success but failure on conversion - big enough files to warrant this care :)
                else: #if not gotta download
                    print("Data wasn't downloaded fully for {}/{}/{} (YR/M/D) and will be redownloaded.\n".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
                    TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
                #END IF
            #END IF
            
        else:
            with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'r') as TEC_file:
                TEC_version = TEC_file["float"].attrs['version']; #read the version
            #END WITH
            if( TEC_version >= version_filt ):
                TEC_dataAvail[i] = 1; #note data is there
            else:
                print("\n==============~Warning~==============");
                print("The file {} uses version {} while version {} is the most current version. The file will be deleted and remade with the new version of the algorithm.\n".format(TEC_dataFilePath[i], TEC_version , version_filt) );
                #the filtered version is too old and needs to be redone
                TEC_dataAvail[i] = 4; #note that data is already downloaded and converted, but not filtered
            #END IF
        #END IF  
    #END FOR i
    
    #==============Make decisions based on what is downloaded/filtered and what overwrite options are engaged (this was indeed tacked on later)==============
    #Now that that the data has been adjusted, fix it up for days that the padding isn't needed or days that data is arleady there but other days need the original data to be filtered  
    if( FLG_overwrite == 0 ): #with this 0, use whatever is available whenever possible
        
        i = dateRange_dayNum_full[:,0].size-1; #set this so I don't have to change anything in my codes
        #find the last one (it's a padded day) - do this first since we want it done before the middle days are compared
        if( (TEC_dataAvail[i] == 1) & (TEC_dataAvail[i-1] != 1) ): #if a padded day happens to be downloaded, force the raw data to be downloaded as it's only for filtering (as long as the prev day isn't also finished - then day isn't needed)
            web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
            web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
            
            TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
            for j in range(0,web_fileNamesIndex.size):
                TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                    TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                #END IF
            #END FOR j
            
            if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                TEC_dataAvail[i] = 4; #unfiltered data is already downloaded and converted
            elif( TEC_fileNameOrig_isThere == True ): #if this is true, orig data is downloaded but not converted to the faster format
                TEC_dataAvail[i] = 3; #orig data is downloaded but not converted
            else: #otherwise get data from download, since data is 1 - assume it exists on the site cause how else would the data have been filtered to start
                TEC_dataAvail[i] = 2; #orig data needs to be downloaded
            #END IF
        #END IF
        elif( TEC_dataAvail[i-1] == 1 ): #prev day is done, don't need padded day downloaded
            TEC_dataAvail[i] = -2; #don't need to import a padded day if the day before is finished
        #END IF
        
        for i in range(0,dateRange_dayNum_full[:,0].size-1): #loop to logic check data choices - doesn't do last cause we did that first
            if( i == 0 ): #find the first one (it's a padded day)
                if( (TEC_dataAvail[i] == 1) & (TEC_dataAvail[i+1] != 1) ): #if a padded day happens to be downloaded, force the raw data to be downloaded as it's only for filtering (as long as the next day isn't also finished - then day isn't needed)                  
                    web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
                    web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
                    
                    TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
                    for j in range(0,web_fileNamesIndex.size):
                        TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                        if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                            TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                        #END IF
                    #END FOR j
                    
                    if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                        TEC_dataAvail[i] = 4; #unfiltered data is already downloaded and converted
                    elif( TEC_fileNameOrig_isThere == True ): #if this is true, orig data is downloaded but not converted to the faster format
                        TEC_dataAvail[i] = 3; #orig data is downloaded but not converted
                    else: #otherwise get data from download, since data is 1 - assume it exists on the site cause how else would the data have been filtered to start
                        TEC_dataAvail[i] = 2; #orig data needs to be downloaded
                    #END IF
                #END IF
                
                if( TEC_dataAvail[i+1] == 1 ): #next day is done, don't need padded day downloaded
                    TEC_dataAvail[i] = -2; #don't need to import a padded day if the day after is finished
                #END IF
                
            else: #not first or last day, so look forward and back
                if( ( (TEC_dataAvail[i+1] != 1) | (TEC_dataAvail[i-1] != 1) ) & (TEC_dataAvail[i] == 1) ): #if the day before or the day after need to be filtered, and the current day is finished, change the signage                   
                    web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
                    web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
                    
                    TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
                    for j in range(0,web_fileNamesIndex.size):
                        TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                        if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                            TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                        #END IF
                    #END FOR j
                    
                    if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                        TEC_dataAvail[i] = 7; #unfiltered data is already downloaded and converted - BUT filtered data is finished so don't filter this day
                    elif( TEC_fileNameOrig_isThere == True ): #if this is true, orig data is downloaded but not converted to the faster format
                        TEC_dataAvail[i] = 6; #orig data is downloaded but not converted - BUT filtered data is finished so don't filter this day
                    else: #otherwise get data from download, since data is 1 - assume it exists on the site cause how else would the data have been filtered to start
                        TEC_dataAvail[i] = 5; #orig data needs to be downloaded - BUT filtered data is finished so don't filter this day
                    #END IF
                #END IF
            #END IF
            
        #END FOR i
        
    elif( FLG_overwrite == 1 ): #FLG_overwrite == 1 flag is on to overwrite all finished data, so just adjust 1's to 2/3/4 depending on what is downloaded
        
        for i in range(0,dateRange_dayNum_full[:,0].size): #loop to check if any data exists - data will not be downloaded yet!           
            web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
            web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
            
            TEC_fileNameOrig_isThere = 1; #prep the flag at 1, goes to 0 if a file that's online isn't local
            for j in range(0,web_fileNamesIndex.size):
                TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                if( os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0 ):
                    TEC_fileNameOrig_isThere = 0; #set the flag to 0 if a file isn't local but it is online
                #END IF
            #END FOR j
            
            if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                TEC_dataAvail[i] = 4; #unfiltered data is already downloaded and converted
            elif( TEC_fileNameOrig_isThere == True ): #if this is true, orig data is downloaded but not converted to the faster format
                TEC_dataAvail[i] = 3; #orig data is downloaded but not converted
            else: #otherwise get data from download, since data is 1 - assume it exists on the site cause how else would the data have been filtered to start
                if( TEC_dataAvail[i] != 0 ): #prevent 0's from being escalated to 2's
                    TEC_dataAvail[i] = 2; #orig data needs to be downloaded
                #END IF
            #END IF
        #END FOR i        
        
    else: #FLG_overwrite == 2 (or anything that's not 0 or 1, really) then, so everything is getting redownloaded
        
        TEC_dataAvail[ TEC_dataAvail > 0 ] = 2; #set all data that is there in some form (numbers greater than 1) to download from source (#2 option)
        
    #END IF
    
        
    #Error check:
    if( np.any(TEC_dataAvail == -1) == 1 ):
        print('Exiting due to no available data on key days.');
    #        import sys
        sys.exit(); #pull it off
        #return("No."); #Can I? I shall.
    #END IF
            
    #==============Download needed data and convert it==============
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
        
        if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == 1): #check for original data requested
            print("\n==============Now Working on Date: {}/{}/{} (Y/M/D) or {}/{} (Y/D#)==============".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2],dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
        else: #else it is a padded day - no filtering will occur
            print("\n==============Now Working on Date: {}/{}/{} (Y/M/D) or {}/{} (Y/D#) - Padded Date, Filter will be Skipped==============".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2],dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
        #END IF
        
        #if( TEC_dataAvail[i] == 1 ): #1 means data is there and we are good
            #pass; #gotta have something in an if statement
            
        #-----DOWNLOAD THE DATA-----
        if( (TEC_dataAvail[i] == 2) | (TEC_dataAvail[i] == 5) ): #2 means data needs to be downloaded (5 means download it for other days, but the current day is filtered and finished already)
            
    #            #First apply code that makes sure the web stuff is for the correct year
    #            if( FLG_newYr != dateRange_full[i,0] ): #check to go get the year's web links and file names and file dates
    #                #basically, if within the same year only need to get this once for the whole year's data availiability
    #                FLG_newYr = dateRange_full[i,0]; #set it so the flag works good
            
            web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
            web_fileNamesIndex = strstr(rendered_content,web_fileNames_dateString); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
            
            print("Downloading {} LISN files for the current day to \"{}\".\nFrom site {}".format( web_fileNamesIndex.size, folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData, web_base_site ));
            tic = time.time(); #for time testing
            estimatedUpdates = np.arange(np.int64(web_fileNamesIndex.size*.1),web_fileNamesIndex.size-np.int64(web_fileNamesIndex.size*.1),np.int64(web_fileNamesIndex.size*.1));
            for j in range(0,web_fileNamesIndex.size):
                TEC_fileNameOrig = rendered_content[web_fileNamesIndex[j]-4:web_fileNamesIndex[j]+9]; #pull out the web file name
                
                if( (os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig) == 0) | (FLG_overwrite == 2) ):
                    #check to see if the data is local, if it isn't download - or if FLG_overwrite is set to 2
                    urlretrieve(web_base_site+TEC_fileNameOrig, folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + TEC_fileNameOrig); #download the file in question to the data directory
                #END IF
                
                if( np.any(j == estimatedUpdates) ):
                    #write an update approximately every 10%
                    sys.stdout.write("\rDownloaded "+str(j+1)+"/"+str(web_fileNamesIndex.size)+" ("+str(np.round((j+1)/web_fileNamesIndex.size*100))+"% finished) ETA is "+str(np.round( ((time.time() - tic)/60)/((j+1)/web_fileNamesIndex.size) - ((time.time() - tic)/60),2))+" min          " );
                    sys.stdout.flush();
                #END IF
            #END FOR j
            toc = time.time() - tic; #for time testing
            print("\nTime to download {} files: {:.2f} min\n\n".format(web_fileNamesIndex.size,toc/60,)); #extra space at end
            
            if( TEC_dataAvail[i] == 2 ):
                TEC_dataAvail[i] = 3; #move on to the next stage ;) oh yeah it feels wrong but its so easy
            else:
                TEC_dataAvail[i] = 6; #move on to the next stage ;) oh yeah it feels wrong but its so easy - but remembers day is filtered and doesn't need to be redone
            #END IF
        #END IF
            
        #-----CONVERT THE DATA TO FASTER FORMAT-----
        if( (TEC_dataAvail[i] == 3) | (TEC_dataAvail[i] == 6) ): #3 means data was downloaded but needs to be converted (6 means data was downloaded but needs to be converted, but the current day is filtered and finished already)
            
            print("Converting raw LISN data for the current day to {}.\nAt 18 data vectors and 90 sec per vector, expect 27 minutes for conversion to finish.\n".format(TEC_dataFileNameUnfilt[i]));
            
            tic = time.time(); #for time testing
            
            #READ the raw data
            web_fileNames_dateString = str(dateRange_dayNum_full[i,0])[2:]+str(dateRange_dayNum_full[i,1]); #make the date string to look for 
            fileNames = os.listdir(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData ); #get the names of the files in the folder
            
            fileNames_GoodMask = np.zeros( len(fileNames) ); #a whitelist of good file names
            for j in range(0,len(fileNames) ):
                fileNames_Index = strstr(fileNames[j],web_fileNames_dateString); #find if the string is in the file name (if empty, it's not)
                
                if( fileNames_Index.size != 0 ):
                    #it's a useful file as its name matches the expected date string
                    fileNames_GoodMask[j] = 1; #note it's a good file
                #END IF
            #END FOR j
            fileNames_GoodMask = np.where( fileNames_GoodMask == 1 )[0]; #convert to indexes cause lists are dum
            
            fileNames_Good = []; #prep a list
            for j in range(0,fileNames_GoodMask.size ):
                fileNames_Good.append(fileNames[fileNames_GoodMask[j]]); #append the good names
            #END FOR j
            
            estimatedUpdates = np.arange(np.int64(fileNames_GoodMask.size*.1),fileNames_GoodMask.size-np.int64(fileNames_GoodMask.size*.1),np.int64(fileNames_GoodMask.size*.1));
            TEC_fileData_pierceAlt = 350; #km, import pierce point altitude - defined as a constant based on communications
            for j in range(0,fileNames_GoodMask.size ):
                #Import raw LISN TEC data file
                with open(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + fileNames_Good[j], 'r') as TEC_file: #open file pointer (with for safety cause life is hard apparently)
                    TEC_raws = TEC_file.readlines(); #multiple headers of varying sizes make reading this real annoying
                #END with
                
                num_dataEntries_PRNnum = 0; #prep counter
                num_dataEntries_headerLen = 0; #prep counter
                for k in range(0, len(TEC_raws)):
                    num_dataEntries_PRNindex = strstr(TEC_raws[k],'PRN');
                    if( (num_dataEntries_headerLen == 0) & (num_dataEntries_PRNindex.size > 0) ):
                        num_dataEntries_headerLen = num_dataEntries_headerLen = k; #record the header length
                        #the header length is found by the first instance of PRN, which signals the data is starting
                    #END IF
                    if( num_dataEntries_PRNindex.size > 0 ):
                        num_dataEntries_PRNnum = num_dataEntries_PRNnum + 1; #increment
                    #END IF
                #END FOR k
                num_dataEntries = len(TEC_raws) - num_dataEntries_PRNnum - num_dataEntries_headerLen; #get the number of real data entries
                
                #Now, collect the date line from the header and confirm it's the correct day
                dateLine = np.where(strfind(TEC_raws[0:num_dataEntries_headerLen], str(dateRange_full[i,0])) == 1)[0]; #get the index of the date line
                if( dateLine.size == 0 ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the year "+str(dateRange_full[i,0])+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.exit();
                elif( dateLine.size > 1 ):
                    #in this instance one of the location numbers happened to have the year in it, so choose the last one
                    dateLine = dateLine[-1];
                #END IF
                dateLine = TEC_raws[np.asscalar(dateLine)].split(' '); #split by a space
                dateLine = list(filter(None, dateLine)); #remote empty entries (that were just a space)
                dateLine_slashes = np.where(strfind(dateLine,'\\') == 1)[0]; #get the index of strings with \'s in them from \n or \r
                if( dateLine_slashes.size == 1 ):
                    dateLine_slashesIndex = strstr(dateLine[np.asscalar(dateLine_slashes)],'\\'); #get the index of where the \ occurs in the string
                    dateLine[np.asscalar(dateLine_slashes)] = dateLine[np.asscalar(dateLine_slashes)][0:np.asscalar(dateLine_slashesIndex)]; #delete the \
                elif( dateLine_slashes.size > 1 ):
                    for k in range(0, dateLine_slashes.size):
                        dateLine_slashesIndex = strstr(dateLine[dateLine_slashes[k]],'\\'); #get the index of where the \ occurs in the string
                        dateLine[dateLine_slashes[k]] = dateLine[dateLine_slashes[k]][0:np.asscalar(dateLine_slashesIndex)]; #delete the \
                    #END FOR k
                #END IF
                if( dateLine[0] != str(dateRange_full[i,0]) ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the correct year "+str(dateRange_full[i,0])+". The year found was "+dateLine[0]+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.exit();
                #END IF
                if( dateLine[1] != str(dateRange_full[i,1]) ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the correct month "+str(dateRange_full[i,1])+". The month found was "+dateLine[1]+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.exit();
                #END IF
                if( dateLine[2] != str(dateRange_full[i,2]) ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the correct day "+str(dateRange_full[i,2])+". The day found was "+dateLine[2]+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.exit();
                #END IF
                
                TEC_fileData_int_temp = np.zeros( (num_dataEntries,locIntUnfilt_size),dtype='int16'); #preallocate
                TEC_fileData_float_temp = np.zeros( (num_dataEntries,locFloatUnfilt_size),dtype='float32'); #preallocate
                TEC_fileData_string_temp = np.zeros( (num_dataEntries),dtype='S4'); #preallocate
                
                #record known constants
                TEC_fileData_int_temp[:,locIntUnfilt_year] = dateRange_full[i,0]; #record the year
    #                TEC_fileData_int_temp[:,locIntUnfilt_month] = dateRange_full[i,1]; #record the month
    #                TEC_fileData_int_temp[:,locIntUnfilt_day] = dateRange_full[i,2]; #record the day
                TEC_fileData_int_temp[:,locIntUnfilt_dayNum] = dateRange_dayNum_full[i,1]; #record the day number
                TEC_fileData_string_temp[:] = fileNames_Good[j][0:4]; #record the site name
                
                lineCounter = 0; #counts the data lines, ignores the PRN lines
                for k in range(num_dataEntries_headerLen, len(TEC_raws)):
                    PRN_check = strstr(TEC_raws[k],'PRN'); #a PRN is a satellite ID #, for GPS it's 1-32 I think, may change w/ constellation and time
                    if( PRN_check.size > 0 ):
                        #found a PRN delcaration line, time to record it
                        current_PRN = TEC_raws[k].split(' '); #split by space
                        current_PRN = list(filter(None, current_PRN)); #remote empty entries (that were just a space)
                        current_PRN = np.int16(current_PRN[ np.asscalar(np.where(strfind(current_PRN, 'PRN') == 1)[0]) + 1 ]); #pulls out the current PRN based on finding where the phrase 'PRN' is and then getting the number right after that
                    else:
                        #format the line for being able to use it
                        current_line = TEC_raws[k].split(' '); #split by space
                        current_line = list(filter(None, current_line)); #remote empty entries (that were just a space)
                        #run code to remove \'s
                        current_slashes = np.where(strfind(current_line,'\\') == 1)[0]; #get the index of strings with \'s in them from \n or \r
                        if( current_slashes.size == 1 ):
                            current_slashesIndex = strstr(current_line[np.asscalar(current_slashes)],'\\'); #get the index of where the \ occurs in the string
                            current_line[np.asscalar(current_slashes)] = current_line[np.asscalar(current_slashes)][0:np.asscalar(current_slashesIndex)]; #delete the \
                        elif( current_slashes.size > 1 ):
                            for k in range(0, current_slashes.size):
                                current_slashesIndex = strstr(current_line[current_slashes[k]],'\\'); #get the index of where the \ occurs in the string
                                current_line[current_slashes[k]] = current_line[current_slashes[k]][0:np.asscalar(current_slashesIndex)]; #delete the \
                            #END FOR k
                        #END IF
                        if( len(current_line) != 12 ):
                            print("\n==============ERROR==============");
                            print("Reading line #"+str(k+1)+" in "+fileNames_Good[j]+" has found less than 12 numbers. Develop code to deal with this.");
                            sys.exit();
                        #END IF
                        
                        #pull out data from the line
                        current_hr = np.int16(current_line[2]); #get the current hour (UT)
                        current_min = np.int16(current_line[3]); #get the current min (UT)
                        current_sec = np.int16(current_line[4]); #get the current sec (UT)
                        current_vTEC = np.float32(current_line[5]); #get the current vTEC (TECU)
                        current_elev = np.float32(current_line[7]); #get the current elevation angle (deg)
                        current_lat = np.float32(current_line[8]); #get the current pierce-point latitude (arcdeg)
                        current_long = np.float32(current_line[9]); #get the current pierce-point longitude (arcdeg)
                        current_sTEC = np.float32(current_line[10]); #get the current sTEC (TECU)
                        
                        #record the data found
                        #record ints
                        TEC_fileData_int_temp[lineCounter,locIntUnfilt_sat] = current_PRN; #record the current satellite ID (PRN)
                        TEC_fileData_int_temp[lineCounter,locIntUnfilt_hour] = current_hr; #record the current hour
                        TEC_fileData_int_temp[lineCounter,locIntUnfilt_min] = current_min; #record the current min
                        TEC_fileData_int_temp[lineCounter,locIntUnfilt_sec] = current_sec; #record the current sec
                        #record floats
                        # TEC_fileData_float_temp[lineCounter,locFloatUnfilt_time] = dateRange_dayNum_full[i,1] + current_hr/24 + current_min/1440 + current_sec/86400; #days, calculate hour/min/sec into days and add to the current day [this is now done later]
                        TEC_fileData_float_temp[lineCounter,locFloatUnfilt_lat] = current_lat; #record the current lat
                        TEC_fileData_float_temp[lineCounter,locFloatUnfilt_long] = current_long; #record the current long
                        TEC_fileData_float_temp[lineCounter,locFloatUnfilt_elev] = current_elev; #record the current elevation
                        TEC_fileData_float_temp[lineCounter,locFloatUnfilt_vTEC] = current_vTEC; #record the current vTEC
                        TEC_fileData_float_temp[lineCounter,locFloatUnfilt_vTEC] = current_sTEC; #record the current sTEC
    #                        TEC_fileData_float_temp[lineCounter,locFloatUnfilt_TECerr] = current_TECerr; #record the current TEC error
                        #TEC err recording is currently disabled b/c it's not in this dataset
                        
                        lineCounter = lineCounter + 1; #increment the line counter
                    #END IF
                #END FOR k
                
                if( num_dataEntries != lineCounter ):
                    print("\n==============ERROR==============");
                    print("Less data than expected was found when reading the data file "+fileNames_Good[j]+". Fix this issue.");
                    sys.exit();
                #END IF
                
                if( j != 0 ): 
                    #following runs, stack the data on
                    TEC_fileData_int = np.concatenate( (TEC_fileData_int, TEC_fileData_int_temp), axis=0); #stack on the data
                    TEC_fileData_float = np.concatenate( (TEC_fileData_float, TEC_fileData_float_temp), axis=0); #stack on the data
                    TEC_fileData_string = np.concatenate( (TEC_fileData_string, TEC_fileData_string_temp), axis=0); #stack on the data
                else:
                    #first run, make the data vars
                    TEC_fileData_int = np.copy(TEC_fileData_int_temp); #copy the data over
                    TEC_fileData_float = np.copy(TEC_fileData_float_temp); #copy the data over
                    TEC_fileData_string = np.copy(TEC_fileData_string_temp); #copy the data over
                #END IF
                
                if( np.any(j == estimatedUpdates) ):
                    #write an update approximately every 10%
                    sys.stdout.write("\rParsed "+str(j+1)+"/"+str(fileNames_GoodMask.size)+" ("+str(np.round((j+1)/fileNames_GoodMask.size*100))+"% finished) ETA is "+str(np.round( ((time.time() - tic)/60)/((j+1)/fileNames_GoodMask.size) - ((time.time() - tic)/60),2))+" min          " );
                    sys.stdout.flush();
                #END IF
            #END FOR j
            
            #-----make sure hour 24 or min 60 doesn't exist (they hsould be rolled over) - found issue via day 129 data showing up b/c it was recorded as d128/h24/m60...------       
            #---CHECK SECONDS---
            kj = TEC_fileData_int[:,locIntUnfilt_sec] >= 60; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main minute variable where seconds are 60 or more
                TEC_fileData_int[kj,locIntUnfilt_min] += 1; #min, increment time by 1
                TEC_fileData_int[kj,locIntUnfilt_sec] -= 60; #sec, remove 60 time units from the time keeping
                if( np.any(TEC_fileData_int[:,locIntUnfilt_sec] >= 60) ):
                    print('ERROR: TIME KEEPING IS REALLY BAD, 60+ SECOND TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ SECOND TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---CHECK MINUTES---
            kj = TEC_fileData_int[:,locIntUnfilt_min] >= 60; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main hour variable where minutes are 60 or more
                TEC_fileData_int[kj,locIntUnfilt_hour] += 1; #hour, increment time by 1
                TEC_fileData_int[kj,locIntUnfilt_min] -= 60; #min, remove 60 time units from the time keeping
                if( np.any(TEC_fileData_int[:,locIntUnfilt_min] >= 60) ):
                    print('ERROR: TIME KEEPING IS REALLY BAD, 60+ MINUTE TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ MIN TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---CHECK HOURS---
            kj = TEC_fileData_int[:,locIntUnfilt_hour] >= 24; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main day variable where hours are 24 or more
                TEC_fileData_int[kj,locIntUnfilt_dayNum] += 1; #day number, increment time by 1
                TEC_fileData_int[kj,locIntUnfilt_hour] -= 24; #hour, remove 24 time units from the time keeping
                if( np.any(TEC_fileData_int[:,locIntUnfilt_hour] >= 24) ):
                    print('ERROR: TIME KEEPING IS REALLY BAD, 24+ HOUR TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 48+ HOUR TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---CHECK DAYS---
            #deal with day limit is based on leap year or not
            dayLim = np.ones(TEC_fileData_int[:,locIntUnfilt_dayNum].shape,dtype=np.int16)*365; #day number, get the day number limits as 365
            #adjust leap years to 366
            leapYears = (np.mod(TEC_fileData_int[:,locIntUnfilt_year],4) == 0) & (np.mod(TEC_fileData_int[:,locIntUnfilt_year],100) != 0) & (np.mod(TEC_fileData_int[:,locIntUnfilt_year],400) == 0)
            dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
            kj = TEC_fileData_int[:,locIntUnfilt_dayNum] >= dayLim; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main year variable where day number is equal to the day number limit or higher than it
                TEC_fileData_int[kj,locIntUnfilt_year] += 1; #year, increment time by 1
                TEC_fileData_int[kj,locIntUnfilt_dayNum] -= dayLim; #hour, remove the day number limit time units from the time keeping
                if( np.any(TEC_fileData_int[:,locIntUnfilt_dayNum] >= dayLim) ):
                    print('ERROR: TIME KEEPING IS REALLY BAD, 365/366 DAY TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 365*2/366*2 DAY TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---luckily years have no limits (that we know of??)---
            #---DELETE data not on the right day (could have incrememted past the right day due to stuff and things like a 24 hour/60 minute time on the current day---
            TEC_fileData_logical_onDay = (TEC_fileData_int[:,locIntUnfilt_dayNum] == dateRange_dayNum_full[i,1]) & (TEC_fileData_int[:,locIntUnfilt_year] == dateRange_dayNum_full[i,0]); #find when the day reported is the day we want [and the year we want]
            TEC_fileData_int = TEC_fileData_int[TEC_fileData_logical_onDay,:]; #keep only the good stuff
            TEC_fileData_float = TEC_fileData_float[TEC_fileData_logical_onDay,:]; #keep only the good stuff
            TEC_fileData_string = TEC_fileData_string[TEC_fileData_logical_onDay]; #keep only the good stuff
            #---UPDATE float date variable---
            TEC_fileData_float[:,locFloatUnfilt_time] = np.float32(TEC_fileData_int[:,locIntUnfilt_dayNum]) + np.float32(TEC_fileData_int[:,locIntUnfilt_hour])/24 + np.float32(TEC_fileData_int[:,locIntUnfilt_min])/1440 + np.float32(TEC_fileData_int[:,locIntUnfilt_sec])/86400; #days, calculate hour/min/sec into days and add to the current day
            
            #WRITE the unfiltered file
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'w') as TEC_fileUnfilt:
                #TEC_fileUnfilt = h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'w'); #write that file new file
                TEC_file_dataSet_int = TEC_fileUnfilt.create_dataset("int", (TEC_fileData_int.shape[0],TEC_fileData_int.shape[1]), dtype='int16' ,compression="gzip"); #create dataset for the integers
                TEC_file_dataSet_float = TEC_fileUnfilt.create_dataset("float", (TEC_fileData_float.shape[0],TEC_fileData_float.shape[1]), dtype='float32' ,compression="gzip"); #create dataset for the floats
                TEC_file_dataSet_string = TEC_fileUnfilt.create_dataset("string", (TEC_fileData_string.shape[0],), dtype='S4' ,compression="gzip"); #create dataset for the strings
                sys.stdout.write("\rWriting ints to file & {} min".format(np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                TEC_file_dataSet_int[...] = TEC_fileData_int; #write that data
                sys.stdout.write("\rWriting floats to file & {} min".format(np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                TEC_file_dataSet_float[...] = TEC_fileData_float; #write that data
                sys.stdout.write("\rWriting strings to file & {} min".format(np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                TEC_file_dataSet_string[...] = TEC_fileData_string; #write that data
                if( np.isscalar(TEC_fileData_pierceAlt) == 0 ):
                    if( np.all( TEC_fileData_pierceAlt == TEC_fileData_pierceAlt[0]) ):
                        TEC_file_dataSet_float.attrs['piercealt'] = TEC_fileData_pierceAlt[0]; #record the attribute
                    else:
                        print("\n==============~Warning~==============");
                        print("Not all pierce point altitudes reported are the same - I don't have anything to deal with this. Continuing. They should print below:");
                        print("{}".format(TEC_fileData_pierceAlt));
                        TEC_file_dataSet_float.attrs['piercealt'] = TEC_fileData_pierceAlt[0]; #record the attribute
                    #END IF
                else:
                    TEC_file_dataSet_float.attrs['piercealt'] = TEC_fileData_pierceAlt; #record the attribute
                #END IF
                TEC_file_dataSet_float.attrs['version'] = version_unfilt; #record the unfiltered algorithm version
                sys.stdout.write("\rDone writing unfiltered file!\t\t\t\t\t");
                sys.stdout.flush();
            #END WITH
            if( FLG_deleteOrig == 1 ): #only delete if flag is on
                for j in range(0,fileNames_GoodMask.size ):
                    os.remove(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + folder_TEC_rawData + '\\' + fileNames_Good[j]); #delete the downloaded file now that we're done with it (save some space!)
                #END FOR j
            #END IF
            del TEC_fileData_int #clean the memory
            del TEC_fileData_float
            del TEC_fileData_string
            del TEC_fileData_pierceAlt
            del TEC_fileData_int_temp
            del TEC_fileData_float_temp
            del TEC_fileData_string_temp
            
            toc = time.time() - tic; #for time testing
            print("\nTime to convert: {} min\n".format(np.round(toc/60,2))); #extra space at end   
            
            if( TEC_dataAvail[i] == 3 ):
                TEC_dataAvail[i] = 4; #move on to the next stage
            else:
                TEC_dataAvail[i] = 7; #move on to the next stage - but remembers day is filtered and doesn't need to be redone
            #END IF
        #END IF
            
    #END FOR
    
    #==============Filter Data after download & conversion==============
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
        
        #-----FILTER THE DATA-----
        if( TEC_dataAvail[i] == 4 ): #4 means data was downloaded and converted but needs to be filtered (this is where 7 jumps off)
            
            #done below with a more accurate TEC_dataRate (based on the site's data)
    #        windowLen_savGol = np.int64(np.round(filter_savGolPeriod*60/TEC_dataRate)); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
    #        #from conversations with AC ^
    #        if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
    #            windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
    #        #END IF
            
            # ===============High Pass Filtering to Excentuate Hourly Periods=========
            filter_hp = 1/(filter_cutoffPeriod*3600);   # highpass cuttoff frequency (not sure what to make of it)
            filter_n = 42; # order of the Hamming window used in the custom function (uses 43 in reality)
    #        filter_f = 1/TEC_dataRate; #1/sec, the sampling frequency, based off of the data interval
    #        filter_wp = 2*filter_hp/filter_f; # Normalizing the frequencies (Matlab is 0 to 1)
    #        filter_b = signal.firwin(filter_n+1, filter_wp, window='hann', pass_zero=False); #just high-pass (2 hr and lower OK)
    #            filter_b = np.array([1]); #bypass filter
            filter_a = 1; #for FIR a is 1
            #Applys Hanning Window for shiggles and giggles
            #NOTE: b is *very slightly* different from MATLAB output.
            
            #-----Declare some constants used to force the data to the required data rate-----
            TEC_dataRate_allowedStamps = np.arange(0,60,TEC_dataRate); #sec, make an array of allowed second timestamps
            TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
            TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
            TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
            
            #TIME TO IMPORT DATA TO FILTER!
            if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == True ): #if statement to find original days requested
                print("\nFiltering {} to {}.".format(TEC_dataFileNameUnfilt[i],TEC_dataFileName[i]) );
                
                tic = time.time(); #for time testing
                TEC_fileData_paddingWarning = []; #prime this list
                for j in range(-1,2): #run through the days to get the main day and the days around it
                    if( (TEC_dataAvail[i+j] == 4) | (TEC_dataAvail[i+j] == 7) | (TEC_dataAvail[i+j] == -3) ): #make sure data is there before reading it (4 or 7 mean it is there) and -3 also means it is there, just a padded day
                        with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i+j], 'r') as TEC_fileUnfilt:
                            #-----import floats (big boys)-----
                            #this is done in waves to keep the memory usage under control
                            #TEC_fileData_float_temp = TEC_fileUnfilt.get("float").value; #.value is being deppreciated apparently
                            TEC_fileData_float_temp = TEC_fileUnfilt.get("float")[()];
                            #Remove data less than min elevation before including it into the main dataset
                            if( j == -1 ):
                                TEC_goodVals = (TEC_fileData_float_temp[:,locFloatUnfilt_time]-np.int64(TEC_fileData_float_temp[:,locFloatUnfilt_time]))*24 > 12; #get UTC hours 12 and greater to cut down on the data imported (would have to think more about the max theoretical observable GPS satellite time period - 12 is safe.)
                                TEC_goodVals = np.where( (TEC_fileData_float_temp[:,locFloatUnfilt_elev] >= minElevation) & TEC_goodVals )[0]; #get the locations where elevation is above the min elevation
                            elif( j == 1 ):
                                TEC_goodVals = (TEC_fileData_float_temp[:,locFloatUnfilt_time]-np.int64(TEC_fileData_float_temp[:,locFloatUnfilt_time]))*24 < 12; #get UTC hours 12 and less to cut down on the data imported (would have to think more about the max theoretical observable GPS satellite time period - 12 is safe.)
                                TEC_goodVals = np.where( (TEC_fileData_float_temp[:,locFloatUnfilt_elev] >= minElevation) & TEC_goodVals )[0]; #get the locations where elevation is above the min elevation
                            else:
                                TEC_goodVals = np.where( (TEC_fileData_float_temp[:,locFloatUnfilt_elev] >= minElevation) )[0]; #get the locations where elevation is above the min elevation
                            #END IF
                            TEC_fileData_float_temp = TEC_fileData_float_temp[TEC_goodVals,:]; #keep the good stuff
                            try: #floats
                                TEC_fileData_float = np.append(TEC_fileData_float,TEC_fileData_float_temp,axis=0); #if var exists, append new data on it
                            except(NameError):
                                TEC_fileData_float = TEC_fileData_float_temp; #if var didn't exist, time to make it exist
                            #END TRY
                            del TEC_fileData_float_temp #clean the memory
                                        
                            #-----import integers (less big int16 only!)-----
                            #this is done in waves to keep the memory usage under control
                            TEC_fileData_int_temp = TEC_fileUnfilt.get("int")[()]; #grab that data
                            TEC_fileData_int_temp = TEC_fileData_int_temp[TEC_goodVals,:]; #keep the good stuff
                            try: #integers
                                TEC_fileData_int = np.append(TEC_fileData_int,TEC_fileData_int_temp,axis=0); #if var exists, append new data on it
                            except(NameError):
                                TEC_fileData_int = TEC_fileData_int_temp; #if var didn't exist, time to make it exist
                            #END TRY
                            del TEC_fileData_int_temp #clean the memory
                            
                            #-----import strings (big bigger because it's like 32 bytes x1)-----
                            #this is done in waves to keep the memory usage under control
                            TEC_fileData_string_temp = TEC_fileUnfilt.get("string")[()];
                            TEC_fileData_string_temp = TEC_fileData_string_temp[TEC_goodVals]; #keep the good stuff
                            TEC_fileData_pierceAlt_temp = TEC_fileUnfilt["float"].attrs['piercealt']; #read the pierce point altitude attribute
                            try: #strings
                                TEC_fileData_string = np.append(TEC_fileData_string,TEC_fileData_string_temp); #if var exists, append new data on it
                            except(NameError):
                                TEC_fileData_string = TEC_fileData_string_temp; #if var didn't exist, time to make it exist
                            #END TRY
                            del TEC_fileData_string_temp #clean up that memory
                            del TEC_goodVals #clean up that memory
                        #END WITH
                        
                        try: #pierce altitude
                            TEC_fileData_pierceAlt = np.append(TEC_fileData_pierceAlt,TEC_fileData_pierceAlt_temp); #if var exists, append new data on it
                        except(NameError):
                            TEC_fileData_pierceAlt = TEC_fileData_pierceAlt_temp; #if var didn't exist, time to make it exist
                        #END TRY
                                            
                    elif(FLG_reqPaddedDays == 0): #otherwise data isn't there - because strict padding is off and data wasn't there
                        if( j == -1 ): #before padding not there
                            TEC_fileData_paddingWarning.append("Before"); #append to var
                        elif( j == 1): #after padding not there
                            TEC_fileData_paddingWarning.append("After"); #append to var
                        else: #something bad
                            print("\n==============ERROR==============");
                            print("There's no data on a main day requested {}/{} (Y/#D) - which shouldn't have happened. Help. Exiting.\n".format(dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
                            sys.exit(); #yolo outa here
                        #END TRY
                    else:
                        print("\n==============ERROR==============");
                        print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are required via passed flag.\n[NOT SUPPORTED IN LISN] Printing availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i+j,0],dateRange_full[i+j,1],dateRange_full[i+j,2]) );
                        print("{}".format(TEC_dataAvail)); #print for error - lets user know available days
                        print("{}".format(i+j)); #print for error - lets user know available days
    #                    print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                        TEC_dataAvail[i+j] = -1; #note the error
                        sys.exit(); #yolo out
                    #END IF              
                #END FOR j
                
                if( np.all(TEC_fileData_pierceAlt == TEC_fileData_pierceAlt[0]) == 1 ): #make sure pierce alt is all the same
                    TEC_fileData_pierceAlt = TEC_fileData_pierceAlt[0]; #just 1 number is good
                else:
                    print("\n==============~Warning~==============");
                    print("Not all pierce point altitudes reported are the same - I don't have anything to deal with this. Continuing. They should print below:");
                    print("{}".format(TEC_fileData_pierceAlt));
                    TEC_fileData_pierceAlt = TEC_fileData_pierceAlt[0]; #record the attribute
                #END IF
                
                del TEC_fileData_pierceAlt_temp #clean memory
                
                toc = time.time() - tic; #for time testing
                print("Time to import all needed data into memory: {} min".format(round(toc/60,2))); #extra space at end  
                
                
                #==============THIS WILL GO IN A FUNCTION FOR MAX SPEED==============
                #TIME TO FILTER!                
                tic = time.time(); #for time testing
                
                TEC_fileData_dTEC = np.ones( TEC_fileData_float.shape[0] , dtype=np.float64 )*np.nan; #preallocate at 64 bit - nan is delete by default
                TEC_fileData_dTECerr = np.zeros( TEC_fileData_float.shape[0] , dtype=np.float64 ); #preallocate at 64 bit
                
    #            TEC_fileData_string_unique = np.unique(TEC_fileData_string); #get unique site names
                TEC_fileData_string_unique, TEC_fileData_string_unique_indexes = np.unique(TEC_fileData_string , return_inverse=True); #get unique site names and the indexes to speed up searching for where to look
                #inspired by https://stackoverflow.com/questions/33281957/faster-alternative-to-numpy-where
    #            TEC_fileData_string_unique_sortedIndex = np.argsort(TEC_fileData_string, kind='mergesort'); #returns the indexes IF the array TEC_fileData_string was sorted
                TEC_fileData_string_unique_currentSiteArray = np.split(np.argsort(TEC_fileData_string, kind='mergesort'), np.cumsum(np.bincount(TEC_fileData_string_unique_indexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
                #faster than for j in range... currentSite_loc = np.where(j == TEC_fileData_string_unique_indexes)[0]; #get the data locations for 1 site per
                        
                TEC_timeUnique, TEC_timeUniqueIndexes , TEC_timeUniqueCount = np.unique( TEC_fileData_float[:,locFloatUnfilt_time] , return_inverse=True , return_counts=True);
                #Cut off time stamps with very little data in the range selected
    #            TEC_dataAvgNum = np.sum(TEC_timeUniqueCount)/TEC_timeUnique.size; #average data per time
    #            TEC_dataLim = np.round(TEC_dataLimPercent*TEC_dataAvgNum); #min number before eject time set
    #            TEC_timeUnique = TEC_timeUnique[ TEC_timeUniqueCount > TEC_dataLim ]; #remove the low data # stuff
                #disabled for now
                TEC_timeUnique_currentTimeArray = np.split(np.argsort(TEC_fileData_float[:,locFloatUnfilt_time], kind='mergesort'), np.cumsum(np.bincount(TEC_timeUniqueIndexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
    
                cntr = 0; #for controlling the number of plots during debugging
                if( FLG_dataAggregation == 0 ): #with this off, just use what is directly there
                    #cruise through every site (from the unique of TEC_fileData_string)
                    for j in range(0,TEC_fileData_string_unique.size): 
        #            for j in range(0,1): #len(TEC_fileData_string_unique)
        #            for j in range(0,np.int64(len(TEC_fileData_string_unique)/4)): #run through all of the sites
                        #j = 1;
                        
        #                currentSite_loc = np.where(j == TEC_fileData_string_unique_indexes)[0]; #get the data locations for 1 site
                        currentSite_loc = TEC_fileData_string_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                                        
                        currentSat_unique = np.unique(TEC_fileData_int[currentSite_loc,locIntUnfilt_sat]); #get the unique sats in at that current site
                        #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])
                                            
                        for k in range(0,currentSat_unique.size):
        #                for k in range(22,len(currentSat_unique)):
                        #k = 21;
                            currentSat_loc = np.where( currentSat_unique[k] == TEC_fileData_int[currentSite_loc,locIntUnfilt_sat] )[0]; #get the data locations for 1 sat at that 1 site
                            currentTEC = TEC_fileData_float[currentSite_loc,locFloatUnfilt_vTEC][currentSat_loc]; #TECU, get the vTEC data for 1 sat at that 1 site - uses abs illegal double indexing
        #                        currentTECerror = TEC_fileData_float[currentSite_loc,locFloatUnfilt_sTECerr][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
    #                        currentDayNum = TEC_fileData_float[currentSite_loc,locFloatUnfilt_time][currentSat_loc]; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentDayNum = np.float64(TEC_fileData_int[currentSite_loc,locIntUnfilt_dayNum][currentSat_loc]) + TEC_fileData_int[currentSite_loc,locIntUnfilt_hour][currentSat_loc]/24 + TEC_fileData_int[currentSite_loc,locIntUnfilt_min][currentSat_loc]/1440 + TEC_fileData_int[currentSite_loc,locIntUnfilt_sec][currentSat_loc]/86400; #days, calculate hour/min/sec into days and add to the current day but do it at 64 bit for extra good precision in the seconds
    #                            currentElv = TEC_fileData_float[currentSite_loc,locFloatUnfilt_elev][currentSat_loc];
                            currentSec = TEC_fileData_int[currentSite_loc,locIntUnfilt_sec][currentSat_loc]; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentMin = TEC_fileData_int[currentSite_loc,locIntUnfilt_min][currentSat_loc]; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentHour = TEC_fileData_int[currentSite_loc,locIntUnfilt_hour][currentSat_loc]; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing

                            currentDayNum_timeSplits_loc = np.append( np.insert( np.where(np.diff(currentDayNum)*1440 > minimumTimeGap)[0]+1 , 0, 0), currentDayNum.shape ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                                                    
                            for l in range(0,len(currentDayNum_timeSplits_loc)-1):
                            #l = 0;
                                
                                #-----Identify current TEC_dataRate & prep the filters-----
                                if( currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]].size > 1 ):
                                    current_TEC_dataRate = np.int16(np.min(np.round(np.diff(currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]])*86400))); #sec, get the current data rate for the site - this assumes the site keeps the data rate the same for all satellites, move deeper into the loop if not
                                else:
                                    current_TEC_dataRate = 0; #set current data rate to 0 so code doesn't run on it 
                                #END IF
                                # current_TEC_dataRateDiffs = np.round(np.diff(currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]])*86400); #sec, get the data rates between time stamps (can be like 30,30,150,30,.. due to gaps)
                                # current_TEC_dataRate = np.where(np.abs(TEC_dataRate - current_TEC_dataRateDiffs) == \
                                #     np.min(np.abs(TEC_dataRate - current_TEC_dataRateDiffs)))[0]; #sec, get the current data rate for the site - this assumes the site keeps the data rate the same for all satellites, move deeper into the loop if not
                                # if( current_TEC_dataRate.size != 0):
                                #     current_TEC_dataRate = current_TEC_dataRateDiffs[current_TEC_dataRate[0]]; #sec, get the current data rate (value closest to desired data rate, usually should be the desired data rate)
                                # else:
                                #     current_TEC_dataRate = 0; #set current data rate to 0 so code doesn't run on it
                                # #END IF
                                
                                if( current_TEC_dataRate != 0 ):
                                    #-----!!I replaced current_TEC_dataRate with TEC_dataRate because I made code that makes all the data be on the right time cadence!!-----
                                    #prep the Sav-Gol filter for debiasing
                                    windowLen_savGol = np.int64(np.round(filter_savGolPeriod*60/TEC_dataRate)); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
                                    #from conversations with AC ^
                                    if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
                                        windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
                                    #END IF
                                    
                                    # ===============High Pass Filtering to Excentuate Hourly Periods=========
                                    filter_f = 1/TEC_dataRate; #1/sec, the sampling frequency, based off of the data interval
                                    filter_wp = 2*filter_hp/filter_f; # Normalizing the frequencies (Matlab is 0 to 1)
                                    filter_b = signal.firwin(filter_n+1, filter_wp, window='hann', pass_zero=False); #just high-pass (2 hr and lower OK)
                                    #filter_b = np.array([1]); #bypass filter
                                    #Applys Hanning Window for shiggles and giggles
                                    #NOTE: b is *very slightly* different from MATLAB output.
                                else:
                                    #otherwise it's 0 cause the data amount is useless, so just use inf for the stuff that matters
                                    windowLen_savGol = np.inf; #can't be greater than this
                                    filter_b = np.array([1]); #bypass filter
                                #END IF
                            
                                #make sure do work only on the day we want
                                #only get the data for times ranges thare are two times the filter period or longer
                                if( ( (currentDayNum[currentDayNum_timeSplits_loc[l+1]-1]-currentDayNum[currentDayNum_timeSplits_loc[l]])*1440 > filter_savGolPeriod*2) & #make sure the data time period is 2*filterPeriod
                                   ( (currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l]) >= windowLen_savGol) & #and that the data number is greater than the winodw length (in the event data is skipped and it's not perfectly 30 sec)
                                   # ( (currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l]) > (3*(filter_b.size-1)) ) & #and that the data number is greater than the filter padding length requirement [not used anymore]
                                   (np.any(np.int16(currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                    
                                    #copy day and TEC data into their own vars
                                    currentDayNum_singleSatLine = currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
        #                            currentTECerror_singleSatLine = currentTECerror[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
        #                            currentDayNum_singleSatLine_min = (currentDayNum_singleSatLine-np.round(np.mean(currentDayNum_singleSatLine)))*144; #min/10, convert to minutes and /10 so it's more stable when fitting a curve
                                    #currentTEC_singleSatLine_zerod = currentTEC_singleSatLine - np.min(currentTEC_singleSatLine); #TECU, adjusted TECU line that is zero'd at minimum
                                    currentTEC_singleSatLine = currentTEC[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
    #                                    currentElv_singleSatLine = currentElv[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]];
                                    currentSec_singleSatLine = currentSec[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                    currentMin_singleSatLine = currentMin[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                    currentHour_singleSatLine = currentHour[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                    
                                    if( np.int16(currentDayNum_singleSatLine[0]) != np.int16(currentDayNum_singleSatLine[-1]) ): #if this is true, we are looking at different days' data being used
                                        #the problem with that is that it seems the site can have different TEC fit values for different days
                                        #so the lines are continuous - but have a big jump where the different days are fit to different means(?) I'm not sure
                                        #but there's def a big jump we're gonna fix here
                                        currentDayNum_singleSatLine_firstDayLocations = np.where( np.int16(currentDayNum_singleSatLine[0]) == np.int16(currentDayNum_singleSatLine) )[0]; #get index where 
                                        currentDayNum_singleSatLine_secondDayLocations = np.where( np.int16(currentDayNum_singleSatLine[-1]) == np.int16(currentDayNum_singleSatLine) )[0]; #get index where 
                                        currentTEC_singleSatLine_diff = np.diff(currentTEC_singleSatLine); #get the deltas between the TEC values
                                        currentTEC_singleSatLine_firstDayLast = currentTEC_singleSatLine[currentDayNum_singleSatLine_firstDayLocations[-1]]; #get the first day's last TEC value
                                        currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations]; #get the second day's TEC values
                                        if( currentDayNum_singleSatLine_firstDayLocations.size > 1 ): #make sure currentDayNum_singleSatLine_firstDayLocations[-2] exists
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_firstDayLocations[-2]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        else: #if it doesn't, time to use currentDayNum_singleSatLine_secondDayLocations[0] instead of currentDayNum_singleSatLine_firstDayLocations[-2] (so projecting from reverse instead of forward)
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_secondDayLocations[0]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        #END IF
                                        currentTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations] = currentTEC_singleSatLine_secondDay; #put em in
                                    #END IF
                                    
                                    #-----This code bit exists to catch if some weird thing happens where 2 satellite lines are recorded at the same time-----
                                    #it happens
                                    if( np.unique(currentDayNum_singleSatLine).size != currentDayNum_singleSatLine.size ):
                                        current_deltaTEC = np.ones(currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l],dtype=np.float32)*np.nan; #no filter, just filler (NaN to easy identify)
                                    else: #otherwise good to go!
                                        
                                        #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                        
                                        #get the time steps right is the first step!
                                        #~~~COMPACT THE TIME STAMPS IN THE DATA~~~
                                        #basically, it needs to meet the lowest common demoninator (30 sec per satellite reading)
                                        #alternatively - the 30 sec data could be repeated to "fill in" the data better, but that would req a lot more memory (more than 16GB) that I don't have in this day and age, but in the future this is where you'd implement something like that
                                        #10 sec data will be compacted by averaging 1-10, 11-20, and 21-30 (shown as 10, 20, 30 in the data) into just 30 at the
                                        #position of the 30 sec data pt
                                        #now it could be also an average of the position datas, but I'm going with the last position for now since I don't know how the position is chosen to start with so the easier way is the way to go!
                                        
                                        if( np.all(np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps) == False) ):
                                            print('CHECK THIS OUT!!! NO 0/30 AT ALL')
                                        #END IF
                                        
                                        #run this if the data rate is equal or below and if all of the sec time stamps don't match the allowed ones (sometimes there can be a 30 sec data rate but it's on 29/59 instead of 30/0)
                                        if( (current_TEC_dataRate <= TEC_dataRate) & np.any(np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps) == False) ): 
                                            
                                            currentSec_singleSatLine_orig = np.copy(currentSec_singleSatLine); #make a copy to know if it changed
                                            allowedStamps_range = np.zeros( (TEC_dataRate_allowedStamps.size,2)); #preallocate
                                            for m in range(0,TEC_dataRate_allowedStamps.size):
                                                allowedStamps_range[m, :] = np.array( (TEC_dataRate_allowedStamps[m]-TEC_dataRate*TEC_timeTolerance,TEC_dataRate_allowedStamps[m]+TEC_dataRate*TEC_timeTolerance) ); #get the current allowed time stamp tolerance range
                                                if( np.any(allowedStamps_range[m,:] < 0) ):
                                                    allowedStamps_range[m, allowedStamps_range[m,:] < 0] = allowedStamps_range[m, allowedStamps_range[m,:] < 0]+60; #keep in the 0-59 sec range
                                                    current_timeInTolerance = (currentSec_singleSatLine >= allowedStamps_range[m,0]) | (currentSec_singleSatLine <= allowedStamps_range[m,1]); #get the times within the time stamp tolerance range (so for 0 sec, 57 to 3 sec is taken as 0 sec)
                                                    if( np.any(currentSec_singleSatLine >= allowedStamps_range[m,0]) ):
                                                        kj = currentSec_singleSatLine >= allowedStamps_range[m,0]; #get where the minute rolls over
                                                        currentMin_singleSatLine[kj] = currentMin_singleSatLine[kj]+1; #min, the minutes are incremented as well if the minute rolled over
                                                        if( np.any(currentMin_singleSatLine == 60) ):
                                                            kj = currentMin_singleSatLine == 60; #get where the hour rolls over
                                                            currentMin_singleSatLine[kj] = 0; #min, set to 0
                                                            currentHour_singleSatLine[kj] = currentHour_singleSatLine[kj]+1; #hour, the hours are incremented as well if the hour rolled over
                                                            if( np.any(currentHour_singleSatLine == 24) ):
                                                                kj = currentHour_singleSatLine == 24; #get where the hour rolls over
                                                                currentHour_singleSatLine[kj] =  0; #hour, set to 0
                                                                currentDayNum_singleSatLine[kj] += 1; #day, the days are incremented as well if the day rolled over
                                                                #update the main day variable
                                                                #there's 2 day number zones, one is the integer one and one is the float with the other time bits in it
                                                                TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_dayNum][kj] += 1; #set the time stamps to be their corrected values in the main array as well
                                                                #-----Leap Year Detection-----
                                                                currentYear_singleSatLine = TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_year]; #get the current years
                                                                currentDayNum_singleSatLine_leapYearDays = np.ones( currentYear_singleSatLine.shape ,dtype=np.int16 )*365; #preallocate, no leap year days
                                                                for k in range(0,currentYear_singleSatLine.size):
                                                                    if np.mod(currentYear_singleSatLine[k],4) == 0: #leap year
                                                                        #-----Leap Year Skipped Detected - next will be 2100-----
                                                                        if (np.mod(currentYear_singleSatLine[k],100) == 0) and (np.mod(currentYear_singleSatLine[k],400) != 0):
                                                                            #Same alg as no leap year used here
                                                                            pass;
                                                                        #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                                                        else:            
                                                                            currentDayNum_singleSatLine_leapYearDays[k] = 366; #note it is a leap year
                                                                        #END IF
                                                                    #-----No leap year detected-----
                                                                    else: #no leap year if this
                                                                        pass;
                                                                    #END IF
                                                                #END FOR k
                                                                if( np.any(np.int16(currentDayNum_singleSatLine) > currentDayNum_singleSatLine_leapYearDays) ):
                                                                    kj = np.int16(currentDayNum_singleSatLine) > currentDayNum_singleSatLine_leapYearDays; #get where the year rolls over
                                                                    currentDayNum_singleSatLine[kj] -= currentDayNum_singleSatLine_leapYearDays; #day number, subtract off the days [this is a float so there's decimal stuff that needs to stick around]
                                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_dayNum][kj] = 1; #day number, set to 1 [start of 1st day of a year]
                                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_year][kj] += 1; #set the time stamps to be their corrected values in the main array as well
                                                                #END IF
                                                            #END IF
                                                            #update the main hour variable
                                                            TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_hour] = np.copy(currentHour_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                        #END IF
                                                        #update the main minute variable
                                                        TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_min] = np.copy(currentMin_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                    #END IF
                                                else:
                                                    current_timeInTolerance = (currentSec_singleSatLine >= allowedStamps_range[m,0]) & (currentSec_singleSatLine <= allowedStamps_range[m,1]); #get the times within the time stamp tolerance range (so for 30 sec, 27 to 33 sec is taken as 30 sec)
                                                #END IF
                                                currentSec_singleSatLine[current_timeInTolerance] = TEC_dataRate_allowedStamps[m]; #set the times within the tolerance range to the allowed time stamp for that range
                                            #END FOR m
                                            if( np.all(currentSec_singleSatLine == currentSec_singleSatLine_orig) == False ):
                                                TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_sec] = np.copy(currentSec_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                            #END IF
                                            
                                            #get where the second time stamps are the allowed values
                                            track_Matches = np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps); #get the seconds that match the allowed times
                                            if( np.any(track_Matches == False) ): #only run if there are still time stmaps that aren't at allowed times
                                                track_Matches_where = np.concatenate( (np.array( (0,) ),np.where(track_Matches == True)[0]) ); #get the indexes, put 0 at the front     
                                                #0 at the front is OK b/c if an allowed time, it'll average [0] and [0] together to get [0]
                                                #if not an allowed time, it'll average [0] and [1] together like it should - basically it gets it started good
                                                
                                                current_time = np.int32(currentDayNum_singleSatLine)*86400 + np.int32(currentHour_singleSatLine)*3600 + np.int32(currentMin_singleSatLine)*60 + np.int32(currentSec_singleSatLine); #get a strong integer representation of the current times
                                                current_timeDiff = np.diff(current_time); #get the diff between the times
                                                track_gaps = np.where( current_timeDiff > current_TEC_dataRate )[0]; #catch time gaps that shouldn't be smashed together
                                                for m in range(0,track_gaps.size):
                                                    if( track_Matches[track_gaps[m]] == False ):
                                                        #this means there's a gap that starts on a non-allowed time (so it's like a 1:45 -> 2:15 gap)
                                                        #force that 1:45 -> 2:00 so everything is cool
                                                        track_Matches[track_gaps[m]] = True; #set to true b/c now it is on an allowed time
                                                        TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
                                                        TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
                                                        TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
                                                        track_shiftSec = np.where(TEC_dataRate_allowedStamps60 > currentSec_singleSatLine[track_gaps[m]])[0][0]; #get where the allowed indexes are greater than the current index, and get the first one always
                                                        currentSec_singleSatLine[track_gaps[m]] = TEC_dataRate_allowedStamps60[track_shiftSec]; #choose the right one
                                                        if( currentSec_singleSatLine[track_gaps[m]] == 60 ):
                                                            kj = currentSec_singleSatLine[track_gaps[m]] == 60; #get where the minute rolls over
                                                            currentMin_singleSatLine[kj] = currentMin_singleSatLine[kj]+1; #min, the minutes are incremented as well if the minute rolled over
                                                            if( np.any(currentMin_singleSatLine == 60) ):
                                                                kj = currentMin_singleSatLine == 60; #get where the hour rolls over
                                                                currentMin_singleSatLine[kj] = 0; #min, set to 0
                                                                currentHour_singleSatLine[kj] = currentHour_singleSatLine[kj]+1; #hour, the hours are incremented as well if the hour rolled over
                                                                if( np.any(currentHour_singleSatLine == 24) ):
                                                                    kj = currentHour_singleSatLine == 24; #get where the hour rolls over
                                                                    currentHour_singleSatLine[kj] =  0; #hour, set to 0
                                                                    currentDayNum_singleSatLine[kj] = np.int16(currentDayNum_singleSatLine[kj])+1; #day, the days are incremented as well if the day rolled over
                                                                    #update the main day variable
                                                                    #there's 2 day number zones, one is the integer one and one is the float with the other time bits in it
                                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_dayNum][kj] += 1; #set the time stamps to be their corrected values in the main array as well
                                                                    #-----Leap Year Detection-----
                                                                    currentYear_singleSatLine = TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_year]; #get the current years
                                                                    currentDayNum_singleSatLine_leapYearDays = np.ones( currentYear_singleSatLine.shape ,dtype=np.int16 )*365; #preallocate, no leap year days
                                                                    for k in range(0,currentYear_singleSatLine.size):
                                                                        if np.mod(currentYear_singleSatLine[k],4) == 0: #leap year
                                                                            #-----Leap Year Skipped Detected - next will be 2100-----
                                                                            if (np.mod(currentYear_singleSatLine[k],100) == 0) and (np.mod(currentYear_singleSatLine[k],400) != 0):
                                                                                #Same alg as no leap year used here
                                                                                pass;
                                                                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                                                            else:            
                                                                                currentDayNum_singleSatLine_leapYearDays[k] = 366; #note it is a leap year
                                                                            #END IF
                                                                        #-----No leap year detected-----
                                                                        else: #no leap year if this
                                                                            pass;
                                                                        #END IF
                                                                    #END FOR k
                                                                    if( np.any(np.int16(currentDayNum_singleSatLine) > currentDayNum_singleSatLine_leapYearDays) ):
                                                                        kj = np.int16(currentDayNum_singleSatLine) > currentDayNum_singleSatLine_leapYearDays; #get where the year rolls over
                                                                        currentDayNum_singleSatLine[kj] -= currentDayNum_singleSatLine_leapYearDays; #day number, subtract off the days [this is a float so there's decimal stuff that needs to stick around]
                                                                        TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_dayNum][kj] = 1; #day number, set to 1 [start of 1st day of a year]
                                                                        TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_year][kj] += 1; #set the time stamps to be their corrected values in the main array as well
                                                                    #END IF
                                                                #END IF
                                                                #update the main hour variable
                                                                TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_hour] = np.copy(currentHour_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                            #END IF
                                                            #force it back to 0
                                                            currentSec_singleSatLine[track_gaps[m]] = 0; #force it back to 0 from 60
                                                            #update the main minute variable
                                                            TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_min] = np.copy(currentMin_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                        #END IF
                                                        #now record the updated sec int
                                                        TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_sec] = np.copy(currentSec_singleSatLine); #save that adjusted sec into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                                                        #update track_Matches_where
                                                        track_Matches_where = np.concatenate( (np.array( (0,) ),np.where(track_Matches == True)[0]) ); #get the indexes, put 0 at the front     
                                                    #END IF
                                                    #otherwise if gap starts on an allowed time I don't think it matters right now, the code will work fine. i think
                                                #END FOR m
                                                if( track_gaps.size > 0 ):
                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_sec] = np.copy(currentSec_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                #END IF
                                                
                                                #now work through it and combine as needed
                                                for m in range(1,track_Matches_where.size-1):
                                                    currentTEC_singleSatLine[track_Matches_where[m+1]] = np.nanmean(currentTEC_singleSatLine[track_Matches_where[m]+1:track_Matches_where[m+1]+1]); #TECU, average them together
                                                    if( currentTEC_singleSatLine[track_Matches_where[m]:track_Matches_where[m+1]+1].size == 0 ):
                                                        sys.crash();
                                                    #END IF
                                                #END FOR m
                                                
                                                if( track_Matches[-1] == False ):
                                                    #check if the last one is an off-time, and in that case
                                                    #push it to the next time slot <- decided to do this b/c predicting could end up rough and deleting is less data while the time scales are short compared to the 1 hour periodicity I'm looking at
                                                    # currentTEC_singleSatLine[-1] = np.nanmean(currentTEC_singleSatLine[track_Matches_where[-1]+1:currentTEC_singleSatLine.size]); #TECU, average the trailing ones
                                                    # if( currentTEC_singleSatLine[track_Matches_where[-1]+1:currentTEC_singleSatLine.size].size == 0 ):
                                                        # sys.crash();
                                                    #END IF
                                                    TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
                                                    TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
                                                    TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
                                                    track_shiftSec = np.where(TEC_dataRate_allowedStamps60 > currentSec_singleSatLine[-1])[0][0]; #get where the allowed indexes are greater than the current index, and get the first one always
                                                    currentSec_singleSatLine[-1] = TEC_dataRate_allowedStamps60[track_shiftSec]; #choose the right one
                                                    if( currentSec_singleSatLine[-1] == 60 ):
                                                        kj = currentSec_singleSatLine[-1] == 60; #get where the minute rolls over
                                                        currentMin_singleSatLine[kj] = currentMin_singleSatLine[kj]+1; #min, the minutes are incremented as well if the minute rolled over
                                                        if( np.any(currentMin_singleSatLine == 60) ):
                                                            kj = currentMin_singleSatLine == 60; #get where the hour rolls over
                                                            currentMin_singleSatLine[kj] = 0; #min, set to 0
                                                            currentHour_singleSatLine[kj] = currentHour_singleSatLine[kj]+1; #hour, the hours are incremented as well if the hour rolled over
                                                            if( np.any(currentHour_singleSatLine == 24) ):
                                                                kj = currentHour_singleSatLine == 24; #get where the hour rolls over
                                                                currentHour_singleSatLine[kj] =  0; #hour, set to 0
                                                                currentDayNum_singleSatLine[kj] = np.int16(currentDayNum_singleSatLine[kj])+1; #day, the days are incremented as well if the day rolled over
                                                                #update the main day variable
                                                                #there's 2 day number zones, one is the integer one and one is the float with the other time bits in it
                                                                TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_dayNum][kj] += 1; #set the time stamps to be their corrected values in the main array as well
                                                                #-----Leap Year Detection-----
                                                                currentYear_singleSatLine = TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_year]; #get the current years
                                                                currentDayNum_singleSatLine_leapYearDays = np.ones( currentYear_singleSatLine.shape ,dtype=np.int16 )*365; #preallocate, no leap year days
                                                                for k in range(0,currentYear_singleSatLine.size):
                                                                    if np.mod(currentYear_singleSatLine[k],4) == 0: #leap year
                                                                        #-----Leap Year Skipped Detected - next will be 2100-----
                                                                        if (np.mod(currentYear_singleSatLine[k],100) == 0) and (np.mod(currentYear_singleSatLine[k],400) != 0):
                                                                            #Same alg as no leap year used here
                                                                            pass;
                                                                        #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                                                        else:            
                                                                            currentDayNum_singleSatLine_leapYearDays[k] = 366; #note it is a leap year
                                                                        #END IF
                                                                    #-----No leap year detected-----
                                                                    else: #no leap year if this
                                                                        pass;
                                                                    #END IF
                                                                #END FOR k
                                                                if( np.any(np.int16(currentDayNum_singleSatLine) > currentDayNum_singleSatLine_leapYearDays) ):
                                                                    kj = np.int16(currentDayNum_singleSatLine) > currentDayNum_singleSatLine_leapYearDays; #get where the year rolls over
                                                                    currentDayNum_singleSatLine[kj] -= currentDayNum_singleSatLine_leapYearDays; #day number, subtract off the days [this is a float so there's decimal stuff that needs to stick around]
                                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_dayNum][kj] = 1; #day number, set to 1 [start of 1st day of a year]
                                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_year][kj] += 1; #set the time stamps to be their corrected values in the main array as well
                                                                #END IF
                                                            #END IF
                                                            #update the main hour variable
                                                            TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_hour] = np.copy(currentHour_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                        #END IF
                                                        #update the main minute variable
                                                        TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_min] = np.copy(currentMin_singleSatLine); #set the time stamps to be their corrected values in the main array as well
                                                        #force it back to 0
                                                        currentSec_singleSatLine[-1] = 0; #force it back to 0 from 60
                                                    #END IF
                                                    track_Matches[-1] = True; #set the last value to true so it isn't NaN'd
                                                    
                                                    #now record the updated sec int
                                                    TEC_fileData_int[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locIntUnfilt_sec] = np.copy(currentSec_singleSatLine); #save that adjusted sec into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                                                #END IF
                                                
                                                #keep the lat/long the same as the original values
                                                
                                                #now NaN the values that didn't match the req time stamps for easy deletion
                                                currentTEC_singleSatLine[np.logical_not(track_Matches)] = np.nan; #nan it up
                                            #END IF                                 
                                        elif( current_TEC_dataRate > TEC_dataRate ):
                                            #need copying code for this
                                            print("ERROR DON\'T HAVE CODE FOR THIS YEET");
                                            sys.exit();
                                        else:
                                            #if all times are OK, make them all true
                                            track_Matches = np.ones(currentTEC_singleSatLine.shape,dtype=bool); #make them all true
                                        #END IF 
                                        
                                        #  & (track_Matches.sum() > (3*(filter_b.size-1))) [removed from if statement b/c highpass filter removed]
                                        if( (track_Matches.sum() >= windowLen_savGol) ): #only run this if there's still enough data points to do the filters
                                            #windowLen_savGol = np.int64(np.round(len(currentTEC_singleSatLine)/4)) - (np.mod(np.int64(np.round(len(currentTEC_singleSatLine)/4)),2) - 1); #unused
                                            currentPolyYvals = np.zeros(currentTEC_singleSatLine.shape); #preallocate
                                            current_deltaTEC = np.zeros(currentTEC_singleSatLine.shape); #preallocate
                                            currentPolyYvals[track_Matches] = savgol_filter(currentTEC_singleSatLine[track_Matches],windowLen_savGol, order_savGol ); #filter it up
                                            current_deltaTEC[track_Matches] = currentTEC_singleSatLine[track_Matches] - currentPolyYvals[track_Matches];
                                            #now NaN the values that didn't match the req time stamps for easy deletion
                                            current_deltaTEC[np.logical_not(track_Matches)] = np.nan; #nan it up
                                            
                                            #-----This code bit exists to catch where the fit greatly deviates from the start OR end, which seems to happen often-----
                                            #use median/median distance instead of mean/stdev for stability apparently
                                            #inspired by https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
        #                                        current_deltaTEC_diff = np.diff(current_deltaTEC); #get the difference between each value and its following neighbor value
        #                                        current_deltaTEC_dist = np.abs(current_deltaTEC_diff - np.median(current_deltaTEC_diff)); #get the absolute distance between the difference and the median of the difference
        #                                        current_deltaTEC_comparitor = current_deltaTEC_dist/np.median(current_deltaTEC_dist); #scale the dist from the median by the median of the dist fromt he median
        #                                        current_deltaTEC_diffWhere = np.where(current_deltaTEC_comparitor > deltaTEC_compareValue)[0]; #get where the scaled distance is greater than the hard-coded comparator value (set to 3.5 when I did this - is a variable)
        #                                        
        #                                        if( current_deltaTEC_diffWhere.size != 0 ): #makes sure the culling only occurs if it is needed (if current_deltaTEC_diffWhere is empty then everything is good)                   
        #                                            current_deltaTEC_diffWhere_diff = np.diff(current_deltaTEC_diffWhere); #the deltas of the indexes of the deltas that are too great am I right
        #                                            current_deltaTEC_diffWhere_diffWhere = np.where(current_deltaTEC_diffWhere_diff >= np.int64(deltaTEC_compareValue*2))[0]; #where the index difference is greater than deltaTEC_compareValue*2
        #                                            if( current_deltaTEC_diffWhere[0] <= np.int64(deltaTEC_compareValue*2) ): #if true, at the start, the actual values and the fit line greatly diverge
        #                                                if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the last current_deltaTEC_diffWhere since we're coming from the start
        #                                                    current_deltaTEC_startCut = current_deltaTEC_diffWhere[ -1 ] + 1; #get the last difference where the delta was out of the delta-TEC-compareValue range and use it.
        #                                                    #add the +1 cause of that slicing weirdness in python
        #                                                else:
        #                                                    current_deltaTEC_startCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[0] ] + 1; #get the first big difference (where the start divergence ends - before the first big difference the index deltas are small till we hit the start!)
        #                                                    #add the +1 cause of that slicing weirdness in python
        #                                                #END IF
        #                                            else:
        #                                                current_deltaTEC_startCut = 0; #set no values to be cut off
        #                                            #END IF
        #                                            if( current_deltaTEC_diffWhere[-1] >= (current_deltaTEC_diff.size - np.int64(deltaTEC_compareValue*2)) ): #if true, at the end, the actual values and the fit line greatly diverge
        #                                                if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the first current_deltaTEC_diffWhere since we're coming from the end
        #                                                    current_deltaTEC_endCut = current_deltaTEC_diffWhere[0] + 1; #get the first diff where the delta was out of the deltaTEC_compareValue range and use it.
        #                                                    #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
        #                                                else:
        #                                                    current_deltaTEC_endCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[-1]+1 ] + 1; #get the last big difference (where the end divergence ends - after the last big difference the index deltas are small till we run hit the end!)
        #                                                    #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
        #                                                #END IF
        #                                            else:
        #                                                current_deltaTEC_endCut = current_deltaTEC.size; #set no values to be cut off
        #                                            #END IF  
        #                                            
        #                                            current_deltaTEC[0:current_deltaTEC_startCut] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
        #                                            current_deltaTEC[current_deltaTEC_endCut:current_deltaTEC.size] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
        #                                        #END IF
                                            
                                            
                                            #DEBUG PLOTS
        #                                    plt.figure();
        #                                    plt.scatter( currentDayNum_singleSatLine , currentTEC_singleSatLine , 20 , "r" );
        #                                    plt.scatter( currentDayNum_singleSatLine , currentPolyYvals , 20 );
        #                                    plt.title("Fit Sav-Gol: Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
        #                                    plt.xlabel('Time [days]');
        #                                    plt.ylabel('los TEC [TECU]');
        #                                    plt.show();
        #                                                                        
        #                                    plt.figure();
        #                                    plt.scatter( currentDayNum_singleSatLine , current_deltaTEC );
        #                                    plt.title("Delta TEC Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
        #                                    plt.xlabel('Time [days]');
        #                                    plt.ylabel('delta los TEC [TECU]');
        #                                    plt.show();
                                            
                                            #adapted from https://en.wikipedia.org/wiki/Coefficient_of_determination
        #                                        SStot = np.nansum( (currentTEC_singleSatLine - np.mean(currentTEC_singleSatLine))**2 );
        #                                        SSreg = np.nansum( (currentTEC_singleSatLine - currentPolyYvals)**2  );
        #                                        Rsquared = 1 - SSreg/SStot;
                                            
        #                                        print('Rsq = '+str(Rsquared));
                                            
        #                                        cntr_minElevation = 0; #prep
        #                                        while( Rsquared < 0.9 ):
        #                                            
        #                                            cntr_minElevation = cntr_minElevation + 5; #increment by 5 degrees
        #                                            
        #                                            if( (minElevation+cntr_minElevation) < 90 ):
        #                                                
        #                                                currentTEC_singleSatLine[ currentElv_singleSatLine > (minElevation+cntr_minElevation) ] = np.nan; #NaN more data that's below this new elevation
        #                                            
        #                                                #-----Re-run the code from above again on the modified currentTEC_singleSatLine
        #                                                try:
        #                                                    currentPolyYvals = savgol_filter(currentTEC_singleSatLine,windowLen_savGol, order_savGol ); #filter it up
        #                                                    current_deltaTEC = currentTEC_singleSatLine - currentPolyYvals;
        #                                                
        #                                                    #-----This code bit exists to catch where the fit greatly deviates from the start OR end, which seems to happen often-----
        #                                                    #use median/median distance instead of mean/stdev for stability apparently
        #                                                    #inspired by https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
        #                                                    current_deltaTEC_diff = np.diff(current_deltaTEC); #get the difference between each value and its following neighbor value
        #                                                    current_deltaTEC_dist = np.abs(current_deltaTEC_diff - np.median(current_deltaTEC_diff)); #get the absolute distance between the difference and the median of the difference
        #                                                    current_deltaTEC_comparitor = current_deltaTEC_dist/np.median(current_deltaTEC_dist); #scale the dist from the median by the median of the dist fromt he median
        #                                                    current_deltaTEC_diffWhere = np.where(current_deltaTEC_comparitor > deltaTEC_compareValue)[0]; #get where the scaled distance is greater than the hard-coded comparator value (set to 3.5 when I did this - is a variable)
        #                                                    
        #                                                    if( current_deltaTEC_diffWhere.size != 0 ): #makes sure the culling only occurs if it is needed (if current_deltaTEC_diffWhere is empty then everything is good)                   
        #                                                        current_deltaTEC_diffWhere_diff = np.diff(current_deltaTEC_diffWhere); #the deltas of the indexes of the deltas that are too great am I right
        #                                                        current_deltaTEC_diffWhere_diffWhere = np.where(current_deltaTEC_diffWhere_diff >= np.int64(deltaTEC_compareValue*2))[0]; #where the index difference is greater than deltaTEC_compareValue*2
        #                                                        if( current_deltaTEC_diffWhere[0] <= np.int64(deltaTEC_compareValue*2) ): #if true, at the start, the actual values and the fit line greatly diverge
        #                                                            if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the last current_deltaTEC_diffWhere since we're coming from the start
        #                                                                current_deltaTEC_startCut = current_deltaTEC_diffWhere[ -1 ] + 1; #get the last difference where the delta was out of the delta-TEC-compareValue range and use it.
        #                                                                #add the +1 cause of that slicing weirdness in python
        #                                                            else:
        #                                                                current_deltaTEC_startCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[0] ] + 1; #get the first big difference (where the start divergence ends - before the first big difference the index deltas are small till we hit the start!)
        #                                                                #add the +1 cause of that slicing weirdness in python
        #                                                            #END IF
        #                                                        else:
        #                                                            current_deltaTEC_startCut = 0; #set no values to be cut off
        #                                                        #END IF
        #                                                        if( current_deltaTEC_diffWhere[-1] >= (current_deltaTEC_diff.size - np.int64(deltaTEC_compareValue*2)) ): #if true, at the end, the actual values and the fit line greatly diverge
        #                                                            if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the first current_deltaTEC_diffWhere since we're coming from the end
        #                                                                current_deltaTEC_endCut = current_deltaTEC_diffWhere[0] + 1; #get the first diff where the delta was out of the deltaTEC_compareValue range and use it.
        #                                                                #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
        #                                                            else:
        #                                                                current_deltaTEC_endCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[-1]+1 ] + 1; #get the last big difference (where the end divergence ends - after the last big difference the index deltas are small till we run hit the end!)
        #                                                                #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
        #                                                            #END IF
        #                                                        else:
        #                                                            current_deltaTEC_endCut = current_deltaTEC.size; #set no values to be cut off
        #                                                        #END IF  
        #                                                        
        #                                                        current_deltaTEC[0:current_deltaTEC_startCut] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
        #                                                        current_deltaTEC[current_deltaTEC_endCut:current_deltaTEC.size] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
        #                                                    #END IF
        #                                                    
        #                                                    #Recalc Rsquared
        #                                                    #adapted from https://en.wikipedia.org/wiki/Coefficient_of_determination
        #                                                    SStot = np.nansum( (currentTEC_singleSatLine - np.mean(currentTEC_singleSatLine))**2 );
        #                                                    SSreg = np.nansum( (currentTEC_singleSatLine - currentPolyYvals)**2  );
        #                                                    Rsquared = 1 - SSreg/SStot;
        #                                                except(np.linalg.LinAlgError):
        #                                                    current_deltaTEC = np.ones(currentTEC_singleSatLine.size)*np.nan; #nan it all
        #                                                    Rsquared = 1.1; #exit out of the loop
        #                                                #END TRY
        #                                                
        #                                                
        #                                            else:
        #                                                Rsquared = 1.1; #exit out of the loop
        #                                                
        #                                                current_deltaTEC = current_deltaTEC*np.nan; #NaN everything because it hit 90 degrees elevation and it was all bad
        #                                            #END IF
        #                                            
        #                                        #END WHILE
                                                                        
                                            # ===============Highpass filtering================
                                            # current_deltaTEC[track_Matches] = signal.filtfilt(filter_b,filter_a,current_deltaTEC[track_Matches],padtype='odd',padlen=3*(filter_b.size-1)); #Appies the filter
                                            #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
                                            
        #                                    plt.figure();
        #                                    plt.scatter( currentDayNum_singleSatLine , current_deltaTEC );
        #                                    plt.title("HP'd Delta TEC Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
        #                                    plt.xlabel('Time [days]');
        #                                    plt.ylabel('delta los TEC [TECU]');
        #                                    plt.show();
                                            
                                            #-----This code bit exists to catch where the deltaTEC greatly exceeds the expected values-----
                                            #Basically, fit deviates (sometimes greatly) on the ends - so before and after currentTEC_singleSatLine ==  currentPolyYvals (approximately) I'm just going to cut off that data, as I've tried controlling it (above)
                                            #This is to prevent the "runoff" look on the ends where it quickly ascends/descends
        #                                        current_deltaTEC_edgeLim = np.where(np.abs(currentTEC_singleSatLine - currentPolyYvals) < TEC_deltaLim)[0]; #get indexes to limit the edge with
        #                                        if( current_deltaTEC_edgeLim.size >= 2 ):
        #                                            if( np.int64(np.round(current_deltaTEC.size*.25)) > current_deltaTEC_edgeLim[0] ):
        #                                                current_deltaTEC[0:current_deltaTEC_edgeLim[0]] = np.nan; #nan out the data that's after the lines deviate
        #                                            #END IF
        #                                            if( np.int64(np.round(current_deltaTEC.size*.75)) < current_deltaTEC_edgeLim[-1] ):
        #                                                current_deltaTEC[current_deltaTEC_edgeLim[-1]+1:current_deltaTEC.size] = np.nan; #nan out the data that's after the lines deviate
        #                                            #END IF
        #                                        #END IF
                                            
                                            #-----Heavy-handed test to limit stuff
        #                                        current_deltaTEC[np.abs(current_deltaTEC) > 1] = np.nan; #nan out the data greater/less than 1
                                            
            #                                plt.figure();
            #                                plt.scatter( currentDayNum_singleSatLine , current_deltaTEC );
            #                                plt.title("HP'd Delta TEC Data Post-More Culling for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
            #                                plt.xlabel('Time [days]');
            #                                plt.ylabel('delta los TEC [TECU]');
            #                                plt.show();
                                            
                                            # ===============Force mean to 0================
                                            #I decided this is wrong
                                            # current_deltaTEC = current_deltaTEC - np.nanmean(current_deltaTEC); #subtract the mean to make it a mean of 0
                                        
        #                                    plt.figure();
        #                                    plt.scatter( currentDayNum_singleSatLine , current_deltaTEC );
        #                                    plt.title("HP'd & Zero-Mean'd TEC Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
        #                                    plt.xlabel('Time [days]');
        #                                    plt.ylabel('delta los TEC [TECU]');
        #                                    plt.show();
            #                                
            #                                cntr = cntr + 1; #increment cntr
            #                                if( cntr == 11 ):
            #                                    manager = _pylab_helpers.Gcf.get_active()
            #                                    if manager is not None:
            #                                        canvas = manager.canvas
            #                                        if canvas.figure.stale:
            #                                            canvas.draw_idle()
            #                                        #END IF
            #                                        canvas.start_event_loop(10);
            #                                    #END IF
            ##                                    input("\n\tpaused: press Enter to continue\n")
            #                                    os.system("pause")
            #                                    cntr = 0; #reset cntr
            #                                    plt.close('all'); #close all plots
            #                                #END IF
                                        else:
                                            #fill with NaNs b/c too few data pts
                                            current_deltaTEC = np.ones(currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l],dtype=np.float32)*np.nan; #no filter, just filler (NaN to easy identify)                       
                                        #END IF - filter size catch round 2
                                    #END IF
                                else:
                                    #this is if it failed time length catch or filter size catch round 1
        #                            currentTEC_singleSatLine = currentTEC[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                    current_deltaTEC = np.ones(currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l],dtype=np.float32)*np.nan; #no filter, just filler (NaN to easy identify)
                                #END IF
                                
                                #no matter what, write it in!
        #                        TEC_fileData_float[currentSite_loc,locFloatUnfilt_dTEC][currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]]  = current_deltaTEC; #save that deltaTEC into wherever it should go - assigning with double indexing DOES NOT WORK
    #                                TEC_fileData_float[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locFloatUnfilt_dTEC] = current_deltaTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                                TEC_fileData_dTEC[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]]] = np.copy(current_deltaTEC); #save that deltaTEC into wherever it should go - getting around double indexing
                                #WRONG BELOW
        #                        TEC_fileData_float[currentSite_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]],locFloatUnfilt_dTEC] = current_deltaTEC; #save that deltaTEC into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                            #END FOR l
                        #END FOR k
                    #ENF FOR j
                    
                else: #with dataAggregation on, time to try to include other data to make longer satellite tracks
                    #cruise through every site (from the unique of TEC_fileData_string)
    #                    cntr = 0; #for testing
                    for j in range(0,TEC_fileData_string_unique.size): 
                        
        #                currentSite_loc = np.where(j == TEC_fileData_string_unique_indexes)[0]; #get the data locations for 1 site
                        currentSite_loc = TEC_fileData_string_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                                        
                        currentSat_unique = np.unique(TEC_fileData_int[currentSite_loc,locIntUnfilt_sat]); #get the unique sats in at that current site
                        #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])
                        
    #                    currentSite_latRange = np.array( (np.min(TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat]), np.max(TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat]) ) ); #get the lat range
    #                    currentSite_longRange = np.array( (np.min(TEC_fileData_float[currentSite_loc,locFloatUnfilt_long]), np.max(TEC_fileData_float[currentSite_loc,locFloatUnfilt_long]) ) ); #get the long range
    #                    
    #                    currentSite_nearbyLat = TEC_fileData_float[ (TEC_fileData_float[:,locFloatUnfilt_lat] >= currentSite_latRange[0]) & (TEC_fileData_float[:,locFloatUnfilt_lat] <= currentSite_latRange[0])  ,locFloatUnfilt_lat]
                        
                        for k in range(0,currentSat_unique.size):
                            currentSat_loc = np.where( currentSat_unique[k] == TEC_fileData_int[currentSite_loc,locIntUnfilt_sat] )[0]; #get the data locations for 1 sat at that 1 site
                            currentTEC = TEC_fileData_float[currentSite_loc,locFloatUnfilt_vTEC][currentSat_loc]; #TECU, get the vTEC data for 1 sat at that 1 site - uses abs illegal double indexing
        #                        currentTECerror = TEC_fileData_float[currentSite_loc,locFloatUnfilt_sTECerr][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentDayNum = TEC_fileData_float[currentSite_loc,locFloatUnfilt_time][currentSat_loc]; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
    #                        currentElv = TEC_fileData_float[currentSite_loc,locFloatUnfilt_elev][currentSat_loc];
                            
                            currentDayNum_timeSplits_loc = np.append( np.insert( np.where(np.diff(currentDayNum)*1440 > minimumTimeGap)[0]+1 , 0, 0), currentDayNum.shape ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                            
                            currentLat = TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat][currentSat_loc]; #degc, get the lat data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentLong = TEC_fileData_float[currentSite_loc,locFloatUnfilt_long][currentSat_loc]; #degc, get the TEC long for 1 sat at that 1 site - uses abs illegal double indexing
                            
                            for l in range(0,len(currentDayNum_timeSplits_loc)-1):
                                currentDayNum_singleSatLine = currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                currentDayNum_singleSatLine_min = (currentDayNum_singleSatLine-np.round(np.mean(currentDayNum_singleSatLine)))*144; #min/10, convert to minutes and /10 so it's more stable when fitting a curve
                                currentTEC_singleSatLine = currentTEC[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
    
                                currentLat_singleSatLine = currentLat[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                currentLong_singleSatLine = currentLong[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                
    #                            currentElv_singleSatLine = currentElv[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]];
    
    #                            plt.figure();
    #                            plt.scatter( currentDayNum_singleSatLine , currentLat_singleSatLine );
    #                            plt.title("Latitude Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Lat [degc]');
    #                            plt.show();
    #                            
    #                            plt.figure();
    #                            plt.scatter( currentDayNum_singleSatLine , currentLong_singleSatLine );
    #                            plt.title("Longitude Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Long [degc]');
    #                            plt.show();
                                
    #                            plt.figure();
    #                            plt.scatter( currentLong_singleSatLine , currentLat_singleSatLine );
    #                            plt.title("Lat/Long Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Long [degc]');
    #                            plt.ylabel('Lat [degc]');
    #                            plt.show();
                                
    #                            plt.figure();
    #                            plt.scatter( currentDayNum_singleSatLine , currentElv_singleSatLine );
    #                            plt.title("Elevation Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Elevation [deg]');
    #                            plt.show();
                                                                                    
                                currentDataAgg_stepsFwdBack= np.int64(TEC_dataAgg_timeAdditionLimit*60/TEC_dataRate); #get the number of steps to take forward and back in time
                                currentDataAgg_dayNumSSLBefore = TEC_timeUnique[ np.where( np.min(currentDayNum_singleSatLine) == TEC_timeUnique )[0][0]-currentDataAgg_stepsFwdBack:np.where( np.min(currentDayNum_singleSatLine) == TEC_timeUnique )[0][0]  ]; #days, get the times before the observation period
                                currentDataAgg_dayNumSSLAfter = TEC_timeUnique[ (np.where( np.max(currentDayNum_singleSatLine) == TEC_timeUnique )[0][0]+1):np.where( np.max(currentDayNum_singleSatLine) == TEC_timeUnique )[0][0]+currentDataAgg_stepsFwdBack+1  ]; #days, get the times before the observation period
                                
                                #only do this if any of the days are in the correct day
                                if( (np.any(np.int16(currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                    
                                    currentDataAgg_dayNumSSLBefore_min = (currentDataAgg_dayNumSSLBefore - np.round(np.mean(currentDayNum_singleSatLine)))*144; #min/10, get the times before the observation period
                                    currentDataAgg_dayNumSSLAfter_min = (currentDataAgg_dayNumSSLAfter - np.round(np.mean(currentDayNum_singleSatLine)))*144; #min/10, get the times after the observation period
        
                                    currentDataAgg_coefs = np.polynomial.polynomial.polyfit(currentDayNum_singleSatLine_min, currentLat_singleSatLine, 3); #gets the coefficients that fit an x^2 polynomial - what a function to call am I right
                                    currentDataAgg_latGuessBefore = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLBefore_min, currentDataAgg_coefs); #calc the guessed values 
                                    currentDataAgg_latGuessAfter = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLAfter_min, currentDataAgg_coefs); #calc the guessed values 
    #                                    currentDataAgg_latGuess = np.hstack( (currentDataAgg_latGuessBefore, currentLat_singleSatLine , currentDataAgg_latGuessAfter) ); #degc, combine lat guess into one variable
                                                                                                     
                                    currentDataAgg_coefs = np.polynomial.polynomial.polyfit(currentDayNum_singleSatLine_min, currentLong_singleSatLine, 3); #gets the coefficients that fit an x^2 polynomial - what a function to call am I right
                                    currentDataAgg_longGuessBefore = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLBefore_min, currentDataAgg_coefs); #calc the guessed values 
                                    currentDataAgg_longGuessAfter = np.polynomial.polynomial.polyval(currentDataAgg_dayNumSSLAfter_min, currentDataAgg_coefs); #calc the guessed values 
    #                                    currentDataAgg_longGuess = np.hstack( (currentDataAgg_longGuessBefore, currentLong_singleSatLine , currentDataAgg_longGuessAfter) ); #degc, combine lat guess into one variable
                                    
    #                                plt.figure();
    #                                plt.scatter( currentDataAgg_longGuess , currentDataAgg_latGuess );
    #                                plt.title("Predicted Lat/Long Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                plt.xlabel('Long [degc]');
    #                                plt.ylabel('Lat [degc]');
    #                                plt.show();
                                    
                                    currentDataAgg_TECBefore = np.ones(currentDataAgg_stepsFwdBack , dtype=np.float32)*np.nan; #preallocate
                                    cntrSkip = 0; #start up a counter
                                    cntrTime = 0; #start up another counter
                                    #run backwards
                                    currentDataAgg_dayNumSSLBeforeRev = np.flip(currentDataAgg_dayNumSSLBefore); #flip it to make it easy
                                    currentDataAgg_dayNumSSLBeforeRevStartIndex = np.where(TEC_timeUnique == currentDataAgg_dayNumSSLBeforeRev[0])[0][0]; #where to start it off
                                    while( (cntrSkip < np.int64(minimumTimeGap*60/TEC_dataRate)) & (cntrTime < currentDataAgg_stepsFwdBack) ):
                                        
                                        currentDataAgg_currentTimeIndexes = TEC_timeUnique_currentTimeArray[currentDataAgg_dayNumSSLBeforeRevStartIndex+cntrTime]; #increment it up as we go to save processing
                                        
                                        currentDataAgg_closeIndexes = np.where( TEC_dataAgg_distToPts_degcSq >= ((TEC_fileData_float[currentDataAgg_currentTimeIndexes,locFloatUnfilt_lat] - currentDataAgg_latGuessBefore[cntrTime])**2 + (TEC_fileData_float[currentDataAgg_currentTimeIndexes,locFloatUnfilt_long] - currentDataAgg_longGuessBefore[cntrTime])**2) )[0];
                                        
                                        if( currentDataAgg_closeIndexes.size != 0 ): #record some data
                                            currentDataAgg_TECBefore[cntrTime] = np.mean(TEC_fileData_float[currentDataAgg_currentTimeIndexes[currentDataAgg_closeIndexes],locFloatUnfilt_vTEC]); #get the mean of the vTEC around the point we guessed
                                            cntrSkip = 0; #reset
                                        else: #else, note data skipped
                                            cntrSkip = cntrSkip + 1; #increment skip counter
                                        #END IF
                                        cntrTime = cntrTime + 1; #increment total time counter
                                    #END WHILE
                                    currentDataAgg_TECBefore = np.flip(currentDataAgg_TECBefore); #reverse back to correct time sequence
                                    
                                    currentDataAgg_TECAfter = np.ones(currentDataAgg_stepsFwdBack , dtype=np.float32)*np.nan; #preallocate
                                    cntrSkip = 0; #start up a counter
                                    cntrTime = 0; #start up another counter
                                    #run forwards
                                    currentDataAgg_dayNumSSLAfterStartIndex = np.where(TEC_timeUnique == currentDataAgg_dayNumSSLAfter[0])[0][0]; #where to start it off
                                    while( (cntrSkip < np.int64(minimumTimeGap*60/TEC_dataRate)) & (cntrTime < currentDataAgg_stepsFwdBack) ):
                                        
                                        currentDataAgg_currentTimeIndexes = TEC_timeUnique_currentTimeArray[currentDataAgg_dayNumSSLAfterStartIndex + cntrTime]; #increment it up as we go to save processing
                                        
                                        currentDataAgg_closeIndexes = np.where( TEC_dataAgg_distToPts_degcSq >= ((TEC_fileData_float[currentDataAgg_currentTimeIndexes,locFloatUnfilt_lat] - currentDataAgg_latGuessAfter[cntrTime])**2 + (TEC_fileData_float[currentDataAgg_currentTimeIndexes,locFloatUnfilt_long] - currentDataAgg_longGuessAfter[cntrTime])**2) )[0];
                                        
                                        if( currentDataAgg_closeIndexes.size != 0 ): #record some data
                                            currentDataAgg_TECAfter[cntrTime] = np.mean(TEC_fileData_float[currentDataAgg_currentTimeIndexes[currentDataAgg_closeIndexes],locFloatUnfilt_vTEC]); #get the mean of the vTEC around the point we guessed
                                            cntrSkip = 0; #reset
                                        else: #else, note data skipped
                                            cntrSkip = cntrSkip + 1; #increment skip counter
                                        #END IF
                                        cntrTime = cntrTime + 1; #increment total time counter
                                    #END WHILE
                                    
                                    currentDataAgg_nansBefore = np.logical_not(np.isnan(currentDataAgg_TECBefore)); #find the NaNs
                                    currentDataAgg_dayNumSSLBefore = currentDataAgg_dayNumSSLBefore[currentDataAgg_nansBefore]; #delete the NaNs
                                    currentDataAgg_TECBefore = currentDataAgg_TECBefore[currentDataAgg_nansBefore]; #delete the NaNs
                                                                    
                                    currentDataAgg_nansAfter = np.logical_not(np.isnan(currentDataAgg_TECAfter)); #find the NaNs
                                    currentDataAgg_dayNumSSLAfter = currentDataAgg_dayNumSSLAfter[currentDataAgg_nansAfter]; #delete the NaNs
                                    currentDataAgg_TECAfter = currentDataAgg_TECAfter[currentDataAgg_nansAfter]; #delete the NaNs
                                    
                                    currentDataAgg_numberBefore = currentDataAgg_dayNumSSLBefore.size; #the size
                                    currentDataAgg_numberAfter = currentDataAgg_dayNumSSLAfter.size; #the size
                                    currentDataAgg_number = currentDataAgg_numberBefore + currentDataAgg_numberAfter; #sum em for ez
                                    currentDataAgg_dayNumSSLExtended = np.hstack( (currentDataAgg_dayNumSSLBefore , currentDayNum_singleSatLine , currentDataAgg_dayNumSSLAfter ) ); #days, stack em
                                                                    
                                else:
                                    currentDataAgg_numberBefore = 0; #set to 0
                                    currentDataAgg_numberAfter = 0; #set to 0
                                    currentDataAgg_number = 0; #set to 0
                                #END IF
                            
                                #make sure do work only on the day we want
                                #only get the data for times ranges thare are two times the filter period or longer
                                if( ( ((currentDayNum[currentDayNum_timeSplits_loc[l+1]-1]-currentDayNum[currentDayNum_timeSplits_loc[l]])*1440+currentDataAgg_number*TEC_dataRate/60) > filter_savGolPeriod*2) & #make sure the data time period is 2*filterPeriod
                                   ( (currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l]+currentDataAgg_number) >= windowLen_savGol) & #and that the data number is greater than the winodw length (in the event data is skipped and it's not perfectly 30 sec)
                                   ( (currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l]+currentDataAgg_number) > (3*(filter_b.size-1)) ) & #and that the data number is greater than the filter padding length requirement (cannot be equal)
                                   (np.any(np.int16(currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                   
                                    if( np.int16(currentDayNum_singleSatLine[0]) != np.int16(currentDayNum_singleSatLine[-1]) ): #if this is true, we are looking at different days' data being used
                                        #the problem with that is that it seems the site can have different TEC fit values for different days
                                        #so the lines are continuous - but have a big jump where the different days are fit to different means(?) I'm not sure
                                        #but there's def a big jump we're gonna fix here
                                        currentDayNum_singleSatLine_firstDayLocations = np.where( np.int16(currentDayNum_singleSatLine[0]) == np.int16(currentDayNum_singleSatLine) )[0]; #get index where 
                                        currentDayNum_singleSatLine_secondDayLocations = np.where( np.int16(currentDayNum_singleSatLine[-1]) == np.int16(currentDayNum_singleSatLine) )[0]; #get index where 
                                        currentTEC_singleSatLine_diff = np.diff(currentTEC_singleSatLine); #get the deltas between the TEC values
                                        currentTEC_singleSatLine_firstDayLast = currentTEC_singleSatLine[currentDayNum_singleSatLine_firstDayLocations[-1]]; #get the first day's last TEC value
                                        currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations]; #get the second day's TEC values
                                        if( currentDayNum_singleSatLine_firstDayLocations.size > 1 ): #make sure currentDayNum_singleSatLine_firstDayLocations[-2] exists
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_firstDayLocations[-2]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        else: #if it doesn't, time to use currentDayNum_singleSatLine_secondDayLocations[0] instead of currentDayNum_singleSatLine_firstDayLocations[-2] (so projecting from reverse instead of forward)
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_secondDayLocations[0]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        #END IF
                                        currentTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations] = currentTEC_singleSatLine_secondDay; #put em in
                                    #END IF
                                    
                                    #now make sure the TEC data before and after connects to the current TEC data (cause of means and day differences and whatever)                                                        
                                    if( currentDataAgg_numberBefore > 0 ): #make sure there's some before data
                                        currentDataAgg_TECBefore = currentDataAgg_TECBefore*( 2*currentTEC_singleSatLine[0] - currentTEC_singleSatLine[1] )/currentDataAgg_TECBefore[-1]; ##scale it so that the last Before data pt is the same as the first real data pt plus the delta between the first and second points
                                    #END IF
                                    if( currentDataAgg_numberAfter > 0 ): #make sure there's some after data
                                        currentDataAgg_TECAfter = currentDataAgg_TECAfter*( 2*currentTEC_singleSatLine[-1] - currentTEC_singleSatLine[-2] )/currentDataAgg_TECAfter[0]; ##scale it so that the first After data pt is the same as the last real data pt plus the delta between the last and second-to-last points
                                    #END IF
                                    currentDataAgg_TEC = np.hstack( ( currentDataAgg_TECBefore, currentTEC_singleSatLine, currentDataAgg_TECAfter ) ); #stack em into one
                                    
                                    
                                    #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                    
                                    #windowLen_savGol = np.int64(np.round(len(currentTEC_singleSatLine)/4)) - (np.mod(np.int64(np.round(len(currentTEC_singleSatLine)/4)),2) - 1); #unused
                                    
                                    currentPolyYvals = savgol_filter(currentDataAgg_TEC,windowLen_savGol,order_savGol ); #filter it up
                                    current_deltaTEC = currentDataAgg_TEC - currentPolyYvals;
                                    
                                    #-----This code bit exists to catch where the fit greatly deviates from the start OR end, which seems to happen often-----
                                    #use median/median distance instead of mean/stdev for stability apparently
                                    #inspired by https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
                                    current_deltaTEC_diff = np.diff(current_deltaTEC); #get the difference between each value and its following neighbor value
                                    current_deltaTEC_dist = np.abs(current_deltaTEC_diff - np.median(current_deltaTEC_diff)); #get the absolute distance between the difference and the median of the difference
                                    current_deltaTEC_comparitor = current_deltaTEC_dist/np.median(current_deltaTEC_dist); #scale the dist from the median by the median of the dist fromt he median
                                    current_deltaTEC_diffWhere = np.where(current_deltaTEC_comparitor > deltaTEC_compareValue)[0]; #get where the scaled distance is greater than the hard-coded comparator value (set to 3.5 when I did this - is a variable)
                                    
                                    if( current_deltaTEC_diffWhere.size != 0 ): #makes sure the culling only occurs if it is needed (if current_deltaTEC_diffWhere is empty then everything is good)                   
                                        current_deltaTEC_diffWhere_diff = np.diff(current_deltaTEC_diffWhere); #the deltas of the indexes of the deltas that are too great am I right
                                        current_deltaTEC_diffWhere_diffWhere = np.where(current_deltaTEC_diffWhere_diff >= np.int64(deltaTEC_compareValue*2))[0]; #where the index difference is greater than deltaTEC_compareValue*2
                                        if( current_deltaTEC_diffWhere[0] <= np.int64(deltaTEC_compareValue*2) ): #if true, at the start, the actual values and the fit line greatly diverge
                                            if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the last current_deltaTEC_diffWhere since we're coming from the start
                                                current_deltaTEC_startCut = current_deltaTEC_diffWhere[ -1 ] + 1; #get the last difference where the delta was out of the delta-TEC-compareValue range and use it.
                                                #add the +1 cause of that slicing weirdness in python
                                            else:
                                                current_deltaTEC_startCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[0] ] + 1; #get the first big difference (where the start divergence ends - before the first big difference the index deltas are small till we hit the start!)
                                                #add the +1 cause of that slicing weirdness in python
                                            #END IF
                                        else:
                                            current_deltaTEC_startCut = 0; #set no values to be cut off
                                        #END IF
                                        if( current_deltaTEC_diffWhere[-1] >= (current_deltaTEC_diff.size - np.int64(deltaTEC_compareValue*2)) ): #if true, at the end, the actual values and the fit line greatly diverge
                                            if( current_deltaTEC_diffWhere_diffWhere.size == 0 ): #in this instance, all values had deltas less than deltaTEC_compareValue*2 - so just choose the first current_deltaTEC_diffWhere since we're coming from the end
                                                current_deltaTEC_endCut = current_deltaTEC_diffWhere[0] + 1; #get the first diff where the delta was out of the deltaTEC_compareValue range and use it.
                                                #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
                                            else:
                                                current_deltaTEC_endCut = current_deltaTEC_diffWhere[ current_deltaTEC_diffWhere_diffWhere[-1]+1 ] + 1; #get the last big difference (where the end divergence ends - after the last big difference the index deltas are small till we run hit the end!)
                                                #+ 1 again cause difference is diff[n] = value[n+1] - value[n] so the +1 value dissapears I think
                                            #END IF
                                        else:
                                            current_deltaTEC_endCut = current_deltaTEC.size; #set no values to be cut off
                                        #END IF  
                                        
                                        current_deltaTEC[0:current_deltaTEC_startCut] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
                                        current_deltaTEC[current_deltaTEC_endCut:current_deltaTEC.size] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
                                    #END IF
                                    
                                    #DEBUG PLOTS
                                    #if( np.any( np.abs(current_deltaTEC) > 1.5) ):
    #                                plt.figure();
    #                                plt.scatter( currentDayNum_singleSatLine , currentTEC_singleSatLine , 20 , "r" );
    #                                plt.scatter( currentDayNum_singleSatLine , currentPolyYvals , 20 );
    #                                plt.title("Fit Sav-Gol: Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                plt.xlabel('Time [days]');
    #                                plt.ylabel('los TEC [TECU]');
    #                                plt.show();
        #                            
        #                            plt.figure();
        #                            plt.scatter( currentDayNum_singleSatLine , currentElv_singleSatLine );
        #                            plt.title("Elevation Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
        #                            plt.xlabel('Time [days]');
        #                            plt.ylabel('Elevation [deg]');
        #                            plt.show();
        #                           
    #                                if(  currentDataAgg_number > 0 ):
    #                                    plt.figure();
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentDataAgg_TEC , 20 , "r" );
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentPolyYvals , 20 );
    #                                    plt.title("Fit Sav-Gol: Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('los TEC [TECU]');
    #                                    plt.show();
    #                                    
    #                                    plt.figure();
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentDataAgg_TEC , 20 , "r" );
    #                                    plt.scatter( currentDayNum_singleSatLine , currentTEC_singleSatLine , 20 );
    #                                    plt.title("OG TEC B, Extra TEC R: Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('los TEC [TECU]');
    #                                    plt.show();
    #                                    
    #                                    plt.figure();
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , current_deltaTEC );
    #                                    plt.title("Projected Delta TEC Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('delta los TEC [TECU]');
    #                                    plt.show();
    #                                #END IF
                                    
                                    
        
                                    
                                                                
                                    # ===============Highpass filtering================
                                    current_deltaTEC = signal.filtfilt(filter_b,filter_a,current_deltaTEC,padtype='odd',padlen=3*(filter_b.size-1)); #Appies the filter
                                    #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
                                    
                                    # ===============Cut to the correct size================
                                    if( (currentDataAgg_numberBefore >= 0) & (currentDataAgg_numberAfter > 0) ):
                                        current_deltaTEC = current_deltaTEC[currentDataAgg_numberBefore:-currentDataAgg_numberAfter]; #remove any of the added TEC values - they were just used to increase the possible "visible" time of a satellite
                                    elif( (currentDataAgg_numberBefore > 0) & (currentDataAgg_numberAfter == 0) ): #this needs some lil help
                                        current_deltaTEC = current_deltaTEC[currentDataAgg_numberBefore:currentDataAgg_dayNumSSLExtended.size+1]; #remove any of the added TEC values - they were just used to increase the possible "visible" time of a satellite
                                    #END IF
                                    
    #                                if(  currentDataAgg_number > 0 ):
    #                                    plt.figure();
    #                                    plt.scatter( currentDayNum_singleSatLine , current_deltaTEC );
    #                                    plt.title("Trimmed Delta TEC Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('delta los TEC [TECU]');
    #                                    plt.show();
    #                                    cntr = cntr + 1;
    #                                    if( cntr == 11 ):
    #                                        sys.exit();
    #                                    #END IF
    #                                #END IF
                                    
                                    # ===============Force mean to 0================
                                    #I decided this is wrong
                                    # current_deltaTEC = current_deltaTEC - np.nanmean(current_deltaTEC); #subtract the mean to make it a mean of 0
                                    if( np.sum(np.isnan(current_deltaTEC)) == current_deltaTEC.size):
                                        sys.exit();
                                    #END IF
                                                            
                                else:
        #                            currentTEC_singleSatLine = currentTEC[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                                    current_deltaTEC = np.ones(currentDayNum_timeSplits_loc[l+1]-currentDayNum_timeSplits_loc[l],dtype=np.float32)*np.nan; #no filter, just filler (NaN to easy identify)
                                #END IF
                                
                                #no matter what, write it in!
        #                        TEC_fileData_float[currentSite_loc,locFloatUnfilt_dTEC][currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]]  = current_deltaTEC; #save that deltaTEC into wherever it should go - assigning with double indexing DOES NOT WORK
    #                                TEC_fileData_float[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]],locFloatUnfilt_dTEC] = current_deltaTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                                TEC_fileData_dTEC[currentSite_loc[currentSat_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]]] = current_deltaTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                                #WRONG BELOW
        #                        TEC_fileData_float[currentSite_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]],locFloatUnfilt_dTEC] = current_deltaTEC; #save that deltaTEC into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                            #END FOR l
                        #END FOR k
                    #ENF FOR j
                #END IF dataAggregation
                
                toc = time.time() - tic; #for time testing
                print("Time to filter all data: {} min".format(np.round(toc/60,2))); #extra space at end  
                tic = time.time(); #start a new timer
                
                #-----Get where the good data is-----
                #---NAN CONTROL---
                #Time to remove NaNs that were accumulated for TEC_fileData_float[:,0] (they mean data was bad)
    #            TEC_fileData_logical_TECnotnans = np.where(np.logical_not(np.isnan(TEC_fileData_float[:,locFloatUnfilt_dTEC])) == True)[0]; #find the not NaNs (index seems to be a bit faster maybe)
                TEC_fileData_logical_TECnotnans = np.logical_not(np.isnan(TEC_fileData_dTEC)); #find the not NaNs    
                #---NON-DAY CONTROL---
                TEC_fileData_logical_onDay = (TEC_fileData_int[:,locIntUnfilt_dayNum] == dateRange_dayNum_full[i,1]) & (TEC_fileData_int[:,locIntUnfilt_year] == dateRange_dayNum_full[i,0]); #find when the day reported is the day we want and the year
                #---OUTLIER CONTROL---
                #adapted from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list b/c median and median dist are more stable
                distFromMedian = np.abs(TEC_fileData_dTEC - np.median(TEC_fileData_dTEC)); #get distance from median
                medianDistFromMedian = np.nanmedian(distFromMedian); #get median distance from median
                normDistFromMedian = distFromMedian/medianDistFromMedian if medianDistFromMedian else 0.0; #get normalized distance from median
                refDistFromMedian = np.abs(np.nanstd(TEC_fileData_dTEC)*3.5 - np.nanmedian(TEC_fileData_dTEC)); #distance from median for 3.5*stdev TECU, median/mean assumed 0, keep ~92% of the data [in practice was 99.95% of data was kept -> implies not accurate gaussian distribution?]
                refDistFromMedian = refDistFromMedian/medianDistFromMedian; #norm it so it can be used to compare
                TEC_fileData_logical_withinBounds = normDistFromMedian < refDistFromMedian; #get where the normalized distance from median is under the reference ceiling normalized value
                TEC_fileData_logical_combined = np.where( (TEC_fileData_logical_TECnotnans & TEC_fileData_logical_onDay & TEC_fileData_logical_withinBounds) == True)[0]; #combine them, get the index (seems slightly faster maybe)
                
                #SPECIAL STUFF OCCURS TO KEEP MEMORY USAGE UNDER CONTROL
                #-----Int clearing-----
                #explicit
                TEC_int = np.zeros( (TEC_fileData_logical_combined.size,locInt_size), dtype=np.int16); #preallocate the ints
                TEC_int[:,locInt_sat] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_sat]; #copy over sat number
                TEC_int[:,locInt_year] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_year]; #copy over year number
                TEC_int[:,locInt_dayNum] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_dayNum]; #copy over day number
                TEC_int[:,locInt_hour] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_hour]; #copy over hour number
                TEC_int[:,locInt_min] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_min]; #copy over min number
                TEC_int[:,locInt_sec] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_sec]; #copy over sec number
                #non-explicit that relies on identical order between unfilt and filt
    #                TEC_int = TEC_fileData_int[TEC_fileData_logical_combined,0:locInt_size]; #copy it over in one shot
                del TEC_fileData_int;
                #THIS MEMORY CULLING THING DID NOT WORK AS INTENDED
                #NOTE THAT THIS "HACKY" MEMORY CULLING ONLY WORKS BY CUTTING OFF THE LAST COLUMN! So... order matters
    #            TEC_int = np.zeros( (TEC_fileData_logical_combined.size, locInt_size) , dtype=np.int16 ); #prep a new array that saves just what we need
    #            TEC_fileData_int.resize( (TEC_fileData_int.shape[0], TEC_fileData_int.shape[1]-6) , refcheck=False); #cuts off the last 2 variables, the los-TEC and error are not carried over into final output stuff
    #            TEC_int[:,locInt_year] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_year]; #year copy over
    #            TEC_fileData_int.resize( (TEC_fileData_int.shape[0], TEC_fileData_int.shape[1]-1) , refcheck=False); #cuts off the last column (which was the var above)
    #            TEC_int[:,locInt_sat] = TEC_fileData_int[TEC_fileData_logical_combined,locIntUnfilt_sat]; #sat number copy over
    #            del TEC_fileData_int #finally delete the whole var
                #-----Float clearing-----
                #explicit
                TEC_float = np.zeros( (TEC_fileData_logical_combined.size,locFloat_size), dtype=np.float32); #preallocate the float32s
                #TEC_float[:,locFloat_time] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_time]; #do not do this, this is old
                TEC_float[:,locFloat_time] += np.float32(np.int32(TEC_int[:,locInt_dayNum])); #use += to conserve memory used, += modifies in-place
                TEC_float[:,locFloat_time] += np.float32(np.int32(TEC_int[:,locInt_hour]))/24; #do it piece by piece so I don't make a ton of copies of these at once
                TEC_float[:,locFloat_time] += np.float32(np.int32(TEC_int[:,locInt_min]))/1440; #do it as float32 to conserve memory
                TEC_float[:,locFloat_time] += np.float32(np.int32(TEC_int[:,locInt_sec]))/86400; #F in the chat for 64 bit accuracy
                TEC_float[:,locFloat_lat] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_lat]; #copy over day number
                TEC_float[:,locFloat_long] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_long]; #copy over hour number
                TEC_float[:,locFloat_elev] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_elev]; #copy over min number
                TEC_float[:,locFloat_dTEC] = TEC_fileData_dTEC[TEC_fileData_logical_combined]; #copy over sat number
                TEC_float[:,locFloat_dTECerr] = TEC_fileData_dTECerr[TEC_fileData_logical_combined]; #copy over sec number
                #non-explicit that relies on identical order between unfilt and filt - mostly won't work
    #                TEC_float = TEC_fileData_float[TEC_fileData_logical_combined,0:locFloat_size]; #copy it over in one shot
                del TEC_fileData_float, TEC_fileData_dTEC, TEC_fileData_dTECerr #finally delete the whole var
                #DID NOT WORK AS INTENDED
    #            TEC_float = np.zeros( (TEC_fileData_logical_combined.size, locFloat_size) , dtype=np.float32 ); #prep a new array that saves just what we need 
    #            TEC_fileData_float.resize( (TEC_fileData_float.shape[0], TEC_fileData_float.shape[1]-2) , refcheck=False); #cuts off the last 2 variables, the los-TEC and error are not carried over into final output stuff
    #            #resizing alg inspired by https://stackoverflow.com/questions/20580775/efficient-way-to-drop-a-column-from-a-numpy-array
    #            #NOTE THAT THIS "HACKY" MEMORY CULLING ONLY WORKS BY CUTTING OFF THE LAST COLUMN! So... order matters
    #            TEC_float[:,locFloat_elev] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_elev]; #elevation copy over
    #            TEC_fileData_float.resize( (TEC_fileData_float.shape[0], TEC_fileData_float.shape[1]-1) , refcheck=False); #cuts off the last column (which was the var above)
    #            TEC_float[:,locFloat_long] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_long]; #longitude copy over
    #            TEC_fileData_float.resize( (TEC_fileData_float.shape[0], TEC_fileData_float.shape[1]-1) , refcheck=False); #cuts off the last column (which was the var above)
    #            TEC_float[:,locFloat_lat] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_lat]; #geodedic lat copy over
    #            TEC_fileData_float.resize( (TEC_fileData_float.shape[0], TEC_fileData_float.shape[1]-1) , refcheck=False); #cuts off the last column (which was the var above)
    #            TEC_float[:,locFloat_time] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_time]; #day# copy over
    #            TEC_fileData_float.resize( (TEC_fileData_float.shape[0], TEC_fileData_float.shape[1]-1) , refcheck=False); #cuts off the last column (which was the var above)
    #            TEC_float[:,locFloat_dTEC] = TEC_fileData_float[TEC_fileData_logical_combined,locFloatUnfilt_dTEC]; #delta-TEC copy over
    #            del TEC_fileData_float #finally delete the whole var
                #-----String clearing-----
                TEC_string = TEC_fileData_string[TEC_fileData_logical_combined]; #just copy it over, donezo
                del TEC_fileData_string, TEC_fileData_logical_combined, TEC_fileData_logical_TECnotnans, TEC_fileData_logical_onDay #delete the vars we're done here
                                
                #-----TIME TO SAVE FILTERED DATA!-----
                
                #with h5py.File(TEC_dataFilePath[i], 'w') as TEC_file:
                #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
                with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'w') as TEC_file:
                    #TEC_fileUnfilt = h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'w'); #write that file new file
                    TEC_file_dataSet_int = TEC_file.create_dataset("int", (TEC_int.shape[0],TEC_int.shape[1]), dtype='int16' ,compression="gzip"); #create dataset for the integers
                    TEC_file_dataSet_float = TEC_file.create_dataset("float", (TEC_float.shape[0],TEC_float.shape[1]), dtype='float32' ,compression="gzip"); #create dataset for the floats
                    TEC_file_dataSet_string = TEC_file.create_dataset("string", (TEC_string.shape[0],), dtype='S4' ,compression="gzip"); #create dataset for the strings
                    sys.stdout.write("\rWriting ints to file & {} min".format(round((time.time()-tic)/60,2)));
                    sys.stdout.flush();
                    TEC_file_dataSet_int[...] = TEC_int; #write that data
                    sys.stdout.write("\rWriting floats to file & {} min".format(round((time.time()-tic)/60,2)));
                    sys.stdout.flush();
                    TEC_file_dataSet_float[...] = TEC_float; #write that data
                    sys.stdout.write("\rWriting strings to file & {} min".format(round((time.time()-tic)/60,2)));
                    sys.stdout.flush();
                    TEC_file_dataSet_string[...] = TEC_string; #write that data
                    
                    TEC_file_dataSet_float.attrs['piercealt'] = TEC_fileData_pierceAlt; #record the attribute
                    if( len(TEC_fileData_paddingWarning) > 0 ):
                        TEC_file_dataSet_float.attrs['paddedDayMissing'] = " and ".join(TEC_fileData_paddingWarning) + " missing"; #record the attribute
                    else:
                        TEC_file_dataSet_float.attrs['paddedDayMissing'] = "None, before and after padded"; #record the attribute
                    #END IF
                    TEC_file_dataSet_float.attrs['forcedTECdataRateSec'] = TEC_dataRate; #record the attribute
                    TEC_file_dataSet_float.attrs['medianRejectRatio'] = deltaTEC_compareValue; #record the attribute
                    TEC_file_dataSet_float.attrs['savgolFiltPeriodHr'] = filter_savGolPeriod/60; #record the attribute
                    TEC_file_dataSet_float.attrs['highpassFiltPeriodHr'] = filter_cutoffPeriod; #record the attribute
                    TEC_file_dataSet_float.attrs['highpassFiltOrd'] = filter_n+1; #record the attribute
                    TEC_file_dataSet_float.attrs['highpassFiltWindow'] = "FIR Hanning"; #record the attribute
                    TEC_file_dataSet_float.attrs['highpassFiltType'] = "signal.filtfilt, padtype='odd',padlen=3*(b.size-1)"; #record the attribute
                    TEC_file_dataSet_float.attrs['version'] = version_filt; #record the filtered algorithm version
                    sys.stdout.write("\rDone writing filtered file!\t\t\t\t\t");
                    sys.stdout.flush();
                #END WITH
                
                del TEC_int, TEC_float, TEC_string; #clean memory
                
                TEC_dataAvail[i] = 7; #data filtered and done! setting to 7 is good, setting to 1 is bad because 7 implies 4 while 1 doesn't guarantee that at alllllll
                #this is an artform of index juggling buddi
                
                toc = time.time() - tic; #for time testing
                print("\nTime to save filtered data: {} min\n".format(np.round(toc/60,2))); #extra space at end  
                
            else:
                if( TEC_dataAvail[i] == 4 ): #only do this if there is data there
                    TEC_dataAvail[i] = -3; #data is a padded day and the filtering is done for non-padded days, so these are donezo
                #END IF
            #END IF
            
        #END IF  
        
        if( (FLG_deleteUnfilt == 1) & ((TEC_dataAvail[i] == 4) | (TEC_dataAvail[i] == 7) | (TEC_dataAvail[i] == 1) | (TEC_dataAvail[i] == -2)) ): #if delete unfiltered data is on - delete it
            os.remove(TEC_dataFilePathUnfilt[i]); #delete the unfiltered file instead of keeping it for more filtering stuff (saves hard drive space)
        #END IF
        
    #END FOR i
        
        
    #==============Read filtered data==============
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
        
        #-----READ THE FILTERED DATA-----
        if( (TEC_dataAvail[i] == 1) | (TEC_dataAvail[i] == 7) ): #1 means data is there, filtered and all - 7 also means that
            
            with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'r') as TEC_file:
                #-----import floats (big boys)-----
                #this is done in waves to keep the memory usage under control
                
                #TEC_fileData_float_temp = TEC_file.get("float").value; #.value is being deppreciated apparently
                TEC_float_temp = TEC_file.get("float")[()]; #get that float data
                #Remove lat/long combos not in the range requested
                TEC_goodVals = np.where( (TEC_float_temp[:,locFloatUnfilt_long] <= np.max(plotLongRange)) & (TEC_float_temp[:,locFloatUnfilt_long] >= np.min(plotLongRange)) & \
                    (TEC_float_temp[:,locFloatUnfilt_lat] <= np.max(plotLatRange)) & (TEC_float_temp[:,locFloatUnfilt_lat] >= np.min(plotLatRange)) )[0]; #delete out of lat/long range stuff
                TEC_float_temp = TEC_float_temp[TEC_goodVals,:]; #keep the good stuff
                try: #floats
                    TEC_float = np.append(TEC_float,TEC_float_temp,axis=0); #if var exists, append new data on it
                except(NameError):
                    TEC_float = TEC_float_temp; #if var didn't exist, time to make it exist
                #END TRY
                del TEC_float_temp; #clean the memory
                            
                #-----import integers (less big int16 only!)-----
                #this is done in waves to keep the memory usage under control
                TEC_int_temp = TEC_file.get("int")[()]; #grab that data
                TEC_int_temp = TEC_int_temp[TEC_goodVals,:]; #keep the good stuff
                try: #integers
                    TEC_int = np.append(TEC_int,TEC_int_temp,axis=0); #if var exists, append new data on it
                except(NameError):
                    TEC_int = TEC_int_temp; #if var didn't exist, time to make it exist
                #END TRY
                del TEC_int_temp; #clean the memory
                
                #-----import strings (big bigger because it's like 32 bytes x1)-----
                #this is done in waves to keep the memory usage under control
                TEC_string_temp = TEC_file.get("string")[()];
                TEC_string_temp = TEC_string_temp[TEC_goodVals]; #keep the good stuff
                try: #strings
                    TEC_string = np.append(TEC_string,TEC_string_temp); #if var exists, append new data on it
                except(NameError):
                    TEC_string = TEC_string_temp; #if var didn't exist, time to make it exist
                #END TRY
                del TEC_string_temp; #clean up that memory
                TEC_pierceAlt_temp = TEC_file["float"].attrs['piercealt']; #read the pierce point altitude attribute
                TEC_paddedDayMissing_temp = TEC_file["float"].attrs['paddedDayMissing']; #read the padded day missing stuff
                TEC_dataRate_temp = TEC_file["float"].attrs['forcedTECdataRateSec']; #read the assumed TEC data rate
                TEC_medianRejectRatio_temp = TEC_file["float"].attrs['medianRejectRatio']; #read the median rejection ratio
                TEC_savgolFiltPeriodHr_temp = TEC_file["float"].attrs['savgolFiltPeriodHr']; #read the Sav-Gol filter period (hr)
                TEC_highpassFiltPeriodHr_temp = TEC_file["float"].attrs['highpassFiltPeriodHr']; #read the high-pass filter period (hr)
                TEC_highpassFiltOrd_temp = TEC_file["float"].attrs['highpassFiltOrd']; #read the high-pass filter order
                TEC_highpassFiltWindow_temp = TEC_file["float"].attrs['highpassFiltWindow']; #read the high-pass filter window type
                TEC_highpassFiltType_temp = TEC_file["float"].attrs['highpassFiltType']; #read the high-pass filter type
            #END WITH
            
            try: #pierce altitude
                TEC_pierceAlt = np.append(TEC_pierceAlt,TEC_pierceAlt_temp); #if var exists, append new data on it
            except(NameError):
                TEC_pierceAlt = TEC_pierceAlt_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_pierceAlt_temp; #clean the memory
            
            try: #padded day missing stuff
                TEC_paddedDayMissing = np.append(TEC_paddedDayMissing,TEC_paddedDayMissing_temp); #if var exists, append new data on it
            except(NameError):
                TEC_paddedDayMissing = TEC_paddedDayMissing_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_paddedDayMissing_temp; #clean the memory
            
    #            try: #assumed TEC data rate
    #                TEC_dataRate = np.append(TEC_dataRate,TEC_dataRate_temp); #if var exists, append new data on it
    #            except(NameError):
    #                TEC_dataRate = TEC_dataRate_temp; #if var didn't exist, time to make it exist
    #            #END TRY
    #            del TEC_dataRate_temp; #clean the memory
            
            try: #median rejection ratio
                TEC_medianRejectRatio = np.append(TEC_medianRejectRatio,TEC_medianRejectRatio_temp); #if var exists, append new data on it
            except(NameError):
                TEC_medianRejectRatio = TEC_medianRejectRatio_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_medianRejectRatio_temp; #clean the memory
            
            try: #Sav-Gol filter period (hr)
                TEC_savgolFiltPeriodHr = np.append(TEC_savgolFiltPeriodHr,TEC_savgolFiltPeriodHr_temp); #if var exists, append new data on it
            except(NameError):
                TEC_savgolFiltPeriodHr = TEC_savgolFiltPeriodHr_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_savgolFiltPeriodHr_temp; #clean the memory
            
            try: #high-pass filter period (hr)
                TEC_highpassFiltPeriodHr = np.append(TEC_highpassFiltPeriodHr,TEC_highpassFiltPeriodHr_temp); #if var exists, append new data on it
            except(NameError):
                TEC_highpassFiltPeriodHr = TEC_highpassFiltPeriodHr_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_highpassFiltPeriodHr_temp; #clean the memory
            
            try: #high-pass filter order
                TEC_highpassFiltOrd = np.append(TEC_highpassFiltOrd,TEC_highpassFiltOrd_temp); #if var exists, append new data on it
            except(NameError):
                TEC_highpassFiltOrd = TEC_highpassFiltOrd_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_highpassFiltOrd_temp; #clean the memory
            
            try: #high-pass filter window type
                TEC_highpassFiltWindow = np.append(TEC_highpassFiltWindow,TEC_highpassFiltWindow_temp); #if var exists, append new data on it
            except(NameError):
                TEC_highpassFiltWindow = TEC_highpassFiltWindow_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_highpassFiltWindow_temp; #clean the memory
            
            try: #high-pass filter type
                TEC_highpassFiltType = np.append(TEC_highpassFiltType,TEC_highpassFiltType_temp); #if var exists, append new data on it
            except(NameError):
                TEC_highpassFiltType = TEC_highpassFiltType_temp; #if var didn't exist, time to make it exist
            #END TRY
            del TEC_highpassFiltType_temp; #clean the memory
                            
        #END IF    
        
    #END FOR i
            
    if( FLG_verboseReturn == 1 ): #if verborse return, report everything that the files recorded
        return TEC_int, TEC_float, TEC_string, TEC_pierceAlt, TEC_paddedDayMissing, TEC_dataRate, \
            TEC_medianRejectRatio, TEC_savgolFiltPeriodHr, TEC_highpassFiltPeriodHr, TEC_highpassFiltOrd, TEC_highpassFiltWindow, TEC_highpassFiltType
    else:
        return TEC_int, TEC_float, TEC_string
    #END IF