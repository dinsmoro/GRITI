#Function to import TEC data from Madrigal3
#RD on 8/23/2018
#
#

import numpy as np
import os, glob, sys
from time import time
from urllib.request import urlopen, urlretrieve
import html2text
from Code.subfun_strstr import strstr
from Code.subfun_strstrNB import strstrNB
from Code.subfun_downloadProgress import downloadProgress
import h5py
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_addADay import subfun_addADay
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.subfun_daysInAYear import subfun_daysInAYear
from Code.subfun_strfind import strfind
from Code.subfun_textNice import textNice
from Code.subfun_satConv import satConv_toStr_fromLongStr, satConv_toStr_fromInt
from Code.GRITI_import_TEC_support import GRITI_import_TEC_support_filter, GRITI_import_TEC_support_filterVersion
#from numba import jit, prange
#-----Testing variables-----
##http://cedar.openmadrigal.org/ftp/fullname/Ross+Dinsmore/email/rld5204@psu.edu/affiliation/PSU/kinst/8000/year/2015/kindat/3505/format/hdf5/
#import matplotlib.pyplot as plt
#from matplotlib import _pylab_helpers
#from Code.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
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
#FLG_overwrite = 1; #0 don't overwrite data even if it is there, 1 overwrite final data even if it is there, 2 overwrite previously downloaded data as well (site files, converted site files)
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


def GRITI_import_TEC_Madrigal(dates, settings, FLG_justChecking=False):
    #==============Version Declaration==============
    version_unfilt = 3.0; #unfiltered algorithm version
    #changelog:
    #1 10/8/2019 - initial algorithm
    #1.1 9/14/2020 - fixed hour 24/minute 60 handling (before ignored, now adjusted to +1 day/+1 hour)
    #2.0 4/9/2021 - moved from monolithic data bloc to slices of data
    #3.0 2/8/2022 - added support for GNSS satellite type (uses RINEX formating | G: GPS, R: GLONASS, S: SBAS Payload, E: Galileo, C: BeiDou, J: QZSS, I: NavIC)
    
    version_filt = GRITI_import_TEC_support_filterVersion(); #filtered algorithm version
    #1 10/8/2019 - initial algorithm
    #1.1 9/11/2020 - fixed non-0/30 time step handling (29/59 were big ones)    
    #1.2 9/16/2020 - removed highpass filter, only dampens data
    #1.3 9/17/2020 - outlier control introduced
    #1.4 9/17/2020 - linear interpolation to make savgol filter work correctly for gappy data
    #2.0 9/17/2020 - interpolated data used (large groups of receivers would not have data at the same time, leaving vertical "streaks" of no-data)
    #3.0 4/9/2021 - moved from monolithic data bloc to slices of data, moved to dict royale
    #3.1 11/30/2021 - fixed missing pierceAlt info & added geomagnetic coordinate support
    #4.0 2/8/2022 - moved to interpolation for most time keeping controls during a pass & supports different satellite types
    
    #==============Unpack==============
    dateRange_dayNum_full = dates['date range full dayNum']; #unpack
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum']; #unpack
    settings_paths = settings['paths']; #unpack
    settings_config = settings['config']; #unpack
    coordType = settings['map']['coord type']; #unpack
    plotLatRange = settings['TEC import']['lat range']; #unpack
    plotLongRange = settings['TEC import']['long range']; #unpack
    web_base_name = settings['TEC import']['web_base_name']; #unpack
    web_base_email = settings['TEC import']['web_base_email']; #unpack
    web_base_affil = settings['TEC import']['web_base_affil']; #unpack
    minElevation = settings['TEC import']['TEC_minElevation']; #unpack
    minimumTimeGap = settings['TEC import']['TEC_minimumTimeGap']; #unpack
    deltaTEC_compareValue = settings['TEC import']['TEC_deltaTEC_compareValue']; #unpack
    TEC_dataRate = settings['TEC import']['TEC_dataRate']; #unpack
    filter_savGolPeriod = settings['TEC import']['filter_savGolPeriod']; #unpack
    order_savGol = settings['TEC import']['order_savGol']; #unpack
    FLG_reqPaddedDays = settings['TEC import']['FLG_reqPaddedDays']; #unpack
    FLG_deleteOrig = settings['TEC import']['FLG_deleteOrig']; #unpack
    FLG_deleteUnfilt = settings['TEC import']['FLG_deleteUnfilt']; #unpack
    FLG_overwrite = settings['TEC import']['FLG_overwrite']; #unpack
    FLG_dataAggregation = settings['TEC import']['FLG_dataAggregation']; #unpack
    TEC_dataAgg_timeAdditionLimit = settings['TEC import']['TEC_dataAgg_timeAdditionLimit']; #unpack
    TEC_dataAgg_distToPts = settings['TEC import']['TEC_dataAgg_distToPts']; #unpack
    # TEC_dataLimPercent = settings['TEC import']['TEC_dataLimPercent']; #unpack
    # TEC_deltaLim = settings['TEC import']['TEC_deltaLim']; #unpack
    TEC_timeTolerance = settings['TEC import']['TEC_timeTolerance']; #unpack
    TEC_maxAmpAllowed = settings['TEC import']['TEC_maxAmpAllowed']; #unpack
    web_maxRetry = settings['TEC import']['max retry']; #unpack
    web_retryWaitSet = settings['TEC import']['retry wait set']; #unpack
    file_maxRetry = settings['TEC import']['max retry file']; #unpac
    
    #==============Constants Needed==============
    paths_TEC = 'TEC'; #name for the TEC folder
    paths_fileEnding = '_Madrigal.h5'; #file extension for hdf5 files
    web_base_site = 'http://cedar.openmadrigal.org'; #website to be used
    web_base = '/ftp/fullname/' + web_base_name + '/email/' + web_base_email + '/affiliation/' + web_base_affil + '/kinst/8000/year/'; #year goes after this
    web_baseAfter = '/kindat/3505/format/hdf5/'; #this goes after year
    #    web_baseFileDL = 'fullFilename/%252Fopt%252Fmadrigal3%252Fexperiments3%252F2016%252Fgps%252F'; #goes after above, for DL of file - date string in form DDmonYY%252Flos_YYYYMMDD
    #    web_baseFileDLAfter = '.001.h5/'; #after bit above, final piece of the puzzle
    earthRadius = 6371; #km, earth radius
    TEC_dataAgg_distToPts_degc = TEC_dataAgg_distToPts/earthRadius*180/np.pi; #degc, distance to extended points for TEC data to exist
    TEC_dataAgg_distToPts_degcSq = TEC_dataAgg_distToPts_degc**2; #degc^2, square it to avoid sqrt later
    
    dataType_str = 'float32'; #declare the float precision here so can change it easily if need be
    dataType_meth = np.float32; #declare the float precision here so can change it easily if need be
    
    #==============File System Locations==============
    #*************************************************
    #==============Filtered File Layout==============
    #Integer Layout
    #0 = Satellite ID [# that corresponds to GPS sat]
    # locInt_sat = 0; #index where sat ID is
    # #1 = Year timestamp [years]
    # locInt_year = 1; #index where year timestamp is
    # #2 = Day Number timestamp [days]
    # locInt_dayNum = 2; #index where day number timestamp is
    # #3 = Hour timestamp [hrs]
    # locInt_hour = 3; #index where hour timestamp is
    # #4 = Minute timestamp [mins]
    # locInt_min = 4; #index where minute timestamp is
    # #5 = Second timestamp [secs]
    # locInt_sec = 5; #index where second timestamp is
    # locInt_size = 6; #size of the int variable
    
    #Float Layout
    #0 = current time in day format [days] - does not support years
    # locFloat_time = 0; #index where time in days is
    # #1 = geodedic latitude [arcdeg]
    # locFloat_lat = 1; #index where geodedic latitude is
    # #2 = longitude [arcdeg]
    # locFloat_long = 2; #index where longitude is
    # #3 = elevation [deg]
    # locFloat_elev = 3; #index where elevation is
    # #4 = delta-TEC "kinda de-biased TEC" [TECU]
    # locFloat_dTEC = 4; #index where delta-TEC is
    # #5 = delta-TEC error [TECU]
    # locFloat_dTECerr = 5; #index where the delta-TEC error is
    # locFloat_size = 6; #size of float variable
    
    #String Layout
    #[] = Receiver Site Name (there's no dim on this one) - use .decode('UTF-8') to make it a string again
    # locString_site = 0; #index where site name is
    # locString_size = 1; #size of string layout
    
    #==============Unfiltered File Layout==============
    #Integer Layout
    #0 = Satellite ID [# that corresponds to GPS sat]
    # locIntUnfilt_sat = 0; #index where sat ID is
    # #1 = Year timestamp [years]
    # locIntUnfilt_year = 1; #index where year timestamp is
    # #2 = Day Number timestamp [days]
    # locIntUnfilt_dayNum = 2; #index where day number is
    # #3 = Hour timestamp [hrs]
    # locIntUnfilt_hour= 3; #index where hour timestamp is
    # #4 = Minute timestamp [mins]
    # locIntUnfilt_min = 4; #index where minute timestamp is
    # #5 = Second timestamp [secs]
    # locIntUnfilt_sec = 5; #index where second timestamp is
    locIntUnfilt_size = 6; #size of the int variable
    
    #Float Layout
    #0 = current time in day format [days] - does not support years
    # locFloatUnfilt_time = 0; #index where time in days is
    # #1 = geodedic latitude [arcdeg]
    # locFloatUnfilt_lat = 1; #index where geodedic latitude is
    # #2 = longitude [arcdeg]
    # locFloatUnfilt_long = 2; #index where longitude is
    # #3 = elevation [deg]
    # locFloatUnfilt_elev = 3; #index where elevation is
    # #4 = line-of-sight TEC [TECU]
    # locFloatUnfilt_sTEC = 4; #index where los TEC is
    # #5 = error in line-of-sight TEC [TECU]
    # locFloatUnfilt_sTECerr = 5; #index where error in los TEC is
    # #6 = vertical TEC [TECU] (which is the sTEC combined with a mapping function)
    # locFloatUnfilt_vTEC = 6; #index where vertical TEC is
    locFloatUnfilt_size = 7; #size of float variable
    
    #String Layout
    #[] = Receiver Site Name (there's no dim on this one) - use .decode('UTF-8') to make it a string again
    # locStringUnfilt_site = 0; #index where site name is
    # locStringUnfilt_satType = 0; #index where sat type is (GPS/GNLONASS/etc..)
    locStringUnfilt_size = 2; #size of string layout
    
    print("\n==============Importing TEC Madrigal Func - Starting==============");
    #==============Adjust dates to be padded==============
    print("Date range requested (yr/day num format): {}/{} to {}/{}.".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
    dateRange_dayNum_full_orig = dateRange_dayNum_full; #record orginal date ranges
    dateRange_dayNum_full = subfun_addADay(dateRange_dayNum_full); #call fun to pad a day onto the ends of the range
    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #convert the padded range to date range
    print("Date range used due to padding requirement for filtering (yr/day num format): {}/{} to {}/{}.\n".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
    #The padding is required for filtering since the two days 
    
    #==============Check if data already there==============
    if( os.path.isdir(settings_paths['data'] + '\\' + paths_TEC) == 0 ): #check if TEC folder exists
        #if not, make it
        os.makedirs(settings_paths['data'] + '\\' + paths_TEC);
        print("NOTA BENE: Importing TEC Madrigal Func - Created TEC directory: {}\n".format(settings_paths['data'] + '\\' + paths_TEC) );
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
    
    TEC_dataPath = ["{}\{}\{}".format(a_, b_, c_) for a_, b_, c_ in zip([settings_paths['data']]*TEC_dataAmnt, [paths_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)) ) ]; #get the base path where data will be in
    TEC_dataFileName = ["{}_{}_{}{}".format(a_, b_, c_, d_) for a_, b_, c_, d_ in zip([paths_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)),np.ndarray.tolist(dateRange_dayNum_full[:,1].astype(str)),[paths_fileEnding]*TEC_dataAmnt ) ]; #get the expected filenames
    TEC_dataFileNameUnfilt = ["{}_{}_{}_unfilt{}".format(a_, b_, c_, d_) for a_, b_, c_, d_ in zip([paths_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)),np.ndarray.tolist(dateRange_dayNum_full[:,1].astype(str)),[paths_fileEnding]*TEC_dataAmnt ) ]; #get the expected filenames for unfiltered data (if stopped mid filtering)
    TEC_dataFilePath = ["{}\{}".format(a_, b_) for a_, b_ in zip(TEC_dataPath,TEC_dataFileName ) ]; #get the full path right to the expected files
    TEC_dataFilePathUnfilt = ["{}\{}".format(a_, b_) for a_, b_ in zip(TEC_dataPath,TEC_dataFileNameUnfilt ) ]; #get the full path right to the expected files that are unfiltered
    
    for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
        if( os.path.isdir(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_uniqueYears[i]) ) == 0 ): #check if date folders exist
            #doesn't exist, gotta make it
            os.makedirs(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_uniqueYears[i]) );
            print("NOTA BENE: Importing TEC Madrigal Func - Created TEC subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_uniqueYears[i]) ));
        #END IF
    #END FOR
        
    #==============Check website and locally for data (in all formats)==============
    FLG_newYr = -1; #start the new year flag - used to detect if a new year occured
    for i in range(0,len(dateRange_dayNum_full[:,0])): #loop to check if any data exists - data will not be downloaded yet!
        #Since data download takes so long (5GB+) I decided to check if it is there before trying to download anything
        if( os.path.isfile(TEC_dataFilePath[i]) == 0):
            #If data doesn't exist, time to get it
            
            #Check if the file is downloaded but not converted - need to deal with 0 inclusion
            #los_20140227.001.h5 - example file name [001 means final data, 00p is custom for preliminary data]
            TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types
            
            if( (os.path.isfile(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == 0) and (os.path.isfile(TEC_dataFilePathUnfilt[i]) == 0) ):
                #if the unconverted file is not there, download it
                #first - need to check if TEC data is availiable for the days requested
                
                if( FLG_newYr != dateRange_full[i,0] ): #check to go get the year's web links and file names and file dates
                    #basically, if within the same year only need to get this once for the whole year's data availiability
                    FLG_newYr = dateRange_full[i,0]; #set it so the flag works good
                
                    web_fileSearch = web_base_site + web_base + str(dateRange_dayNum_full[i,0]) + web_baseAfter; #build site to go to
                    web_retryNum = 1; #prep the retry num
                    web_retryWait = web_retryWaitSet[0]; #prep the retry wait period
                    FLG_gotPage = False; #lets while exit
                    while( ((web_maxRetry >= web_retryNum) | (web_maxRetry == -1)) & (FLG_gotPage == False) ):
                        try:
                            page = urlopen(web_fileSearch); #get raw HTML
                            FLG_webWork = True; #flag for web working
                            FLG_gotPage = True; #lets page exit
                        except:
                            FLG_webWork = False; #flag for web no work
                            sys.stdout.write('\rNo web access for Madrigal data avail search. Try #{}/{}, will wait {} sec until next try\t'.format(web_retryNum,web_maxRetry,web_retryWait)); #report
                            sys.stdout.flush();
                            time.sleep(web_retryWait);
                            web_retryNum += 1; #increment try
                            if( web_retryWait < web_retryWaitSet[1] ):
                                web_retryWait += web_retryWaitSet[2]; #increment by amount
                            #END IF
                        #END TRYING
                    #END WHILE
                    if( FLG_webWork == True ):
                        if( web_retryNum > 1 ):
                            print('\n'); #for the no web access info
                        #END IF
                        html_content = page.read(); #read off the HTML from whatever the page holder is
                        charset = page.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                        if( charset is None ):
                            charset = 'utf-8'; #assume utf-8
                        #END IF
                        html_content = html_content.decode(charset); #"decode" the HTML content so it's legible
                        rendered_content = html2text.html2text(html_content); #render the HTML like webpage would and get the real stuff
                        #print("{}".format(rendered_content)); #print for debug
                        web_fileNamesIndex = strstrNB(rendered_content,'los_'); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
                        web_fileNotNamesIndex = strstrNB(rendered_content,'%252Flos_')+5; #pull out indexes related to the end of the file links only ~for removing from above~
                        web_fileLinksIndex = strstrNB(rendered_content,'.h5](/ftp')+5; #pull out indexes related to the ~start~ of the file links only
                        web_fileLinksEndIndex = strstrNB(rendered_content,'.h5/)')+4; #pull out indexes related to the ~end~ of the file links olnly   
                    
                        web_fileNamesIndex = np.delete( web_fileNamesIndex,np.where(np.in1d(web_fileNamesIndex, web_fileNotNamesIndex))[0]); #search for the not wanted file name ending and remove them
                        
                        web_fileNames = []; #prep a list for strings, names of the files
                        web_fileNamesDate = []; #prep a list for strings, dates of the files
                        web_fileLinks = []; #prep a list for strings, links to the files
                        for j in range(0,len(web_fileNamesIndex)):
                            web_fileNames.append(rendered_content[web_fileNamesIndex[j]:web_fileNamesIndex[j]+19]); #get the file name recorded
                            web_fileNamesDate.append(rendered_content[web_fileNamesIndex[j]+4:web_fileNamesIndex[j]+12]); #get the file name's date recorded
                            web_fileLinks.append(rendered_content[web_fileLinksIndex[j]:web_fileLinksEndIndex[j]]); #get the file name recorded
                        #END FOR j
                        web_fileNamesDate = np.swapaxes(np.array([np.asarray([w[0:4] for w in web_fileNamesDate]).astype("int16"),np.asarray([w[4:6] for w in web_fileNamesDate]).astype("int16"),np.asarray([w[6:8] for w in web_fileNamesDate]).astype("int16")]),0,1); #pull out the yr/month/day in the file name in M/D/YR format
                    else:
                        print('\nWARNING: No internet access, can\'t check Madrigal for latest data availability.');
                    #END IF
                #END IF
                
                #Check if website has data for the date requested
                if( FLG_webWork == True ):
                    if( np.any(np.all(( web_fileNamesDate - dateRange_full[i,:] )  == 0, axis=1 ) == 1 )): #basically checks the date vector with all the dates on the website, if one matches we're good
                        TEC_dataAvail[i] = 2; #note data is there but needs to be downloaded
                    else:
                        if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == 1 ): 
                            #Catch where no data available
                            print("\n==============ERROR==============");
                            print("There is no data availiable on {}/{}/{} in YR/M/D format. Printing availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
                            print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                            TEC_dataAvail[i] = -1; #note the error
                        else:
                            #Catch where no data available - but it is a padded day to help with filtering on day edges
                            if( FLG_reqPaddedDays == 0 ):
                                print("\n==============~Warning~==============");
                                print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are NOT required via passed flag - program will not quit but filtering on day edges will be impacted.\nPrinting availiable data days for the relevant year:".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
                                print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                                TEC_dataAvail[i] = 0; #note the error
                            else:
                                print("\n==============ERROR==============");
                                print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are required via passed flag.\nPrinting availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
                                print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                                TEC_dataAvail[i] = -1; #note the error
                            #END IF
                        #END IF
                    #END IF
                else:
                    TEC_dataAvail[i] = -1; #note the error, this is just to punt (seems to work if data already avail and calc'd, may not in some instances (e.g. calc needs to occur but raw data avail))
                #END IF
            elif( os.path.isfile(TEC_dataFilePathUnfilt[i]) == 1 ): #check if unfiltered data file exists            
                try: #gonna try to read the file - if we fail, it failed mid download probably
                    with h5py.File(TEC_dataFilePathUnfilt[i], 'r', rdcc_nbytes =500*1024*1024) as testFile:
                        testFile.keys(); #tries to check out some stuff in the file
                    #END WITH
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'r', rdcc_nbytes =500*1024*1024) as TEC_fileUnfilt:
                        TEC_unfilt_version = TEC_fileUnfilt.attrs['version']; #read the version
                    #END WITH
                    if( TEC_unfilt_version >= version_unfilt ):
                        #if the unfiltered file is there (filtering got stopped mid-filter or something)
                        TEC_dataAvail[i] = 4; #note that data is already downloaded and converted, but not filtered
                        #use case if download and conversion success but failure on filter - big enough files to warrant this care :)
                    else:
                        print("\n==============~Warning~==============");
                        print("The file {} uses version {} while version {} is the most current version. The file will be renamed with _oldV{} and remade with the new version of the algorithm.\n".format(TEC_dataFilePathUnfilt[i], TEC_unfilt_version , version_unfilt, TEC_unfilt_version) );
                        os.rename(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], \
                                  settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i]+'_oldV'+str(TEC_unfilt_version).replace('.','p')); #rename
                        #force an OSError to use the rest of the algorithm - uses ***'s which can't be in file names
                        with h5py.File(TEC_dataFilePathUnfilt[i]+'.***notapossiblefile', 'r', rdcc_nbytes =500*1024*1024) as testFile:
                            testFile.keys(); #tries to check out some stuff in the file
                        #END WITH
                    #END IF
                except (OSError , KeyError):
                    try:
                        os.remove(TEC_dataFilePathUnfilt[i]); #delete partially downloaded file to avoid conflicts
                        print("\n==============~Warning~==============");
                        print("Data wasn't converted to standardized format & naming scheme for {} and will be reconverted. Deleting partial file as well.\n".format(TEC_dataFilePathUnfilt[i]) );
                    except:
                        pass;
                    #END TRY
                    try: #gonna try to read the file - if we fail, it failed mid download
                        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r', rdcc_nbytes =500*1024*1024) as testFile:
                            testFile.keys(); #tries to check out some stuff in the file
                        #END WITH
                        TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                        #use case if download success but failure on conversion - big enough files to warrant this care :)
                    except (OSError , KeyError):
                        print("\n==============~Warning~==============");
                        print("Data wasn't downloaded fully for {} and will be redownloaded. Deleting partial file as well.\n".format(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) );
                        os.remove(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete partially downloaded file to avoid conflicts
                        TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
                    #END TRYING
                #END TRYING
            else: #otherwise orig data downloaded but not converted to standard (fast) format
                try: #gonna try to read the file - if we fail, it failed mid download
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r', rdcc_nbytes =500*1024*1024) as testFile:
                        testFile.keys(); #tries to check out some stuff in the file
                    #END WITH
                    TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                    #use case if download success but failure on conversion - big enough files to warrant this care :)
                except (OSError , KeyError):
                    print("\n==============~Warning~==============");
                    print("Data wasn't downloaded fully for {} and will be redownloaded. Deleting partial file as well.\n".format(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) );
                    os.remove(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete partially downloaded file to avoid conflicts
                    TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
                #END TRYING
            #END IF
            
        else:
            try:
                with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'r', rdcc_nbytes =500*1024*1024) as TEC_file:
                    TEC_version = TEC_file.attrs['version']; #read the version
                #END WITH
            except(OSError, KeyError):
                TEC_version = -1; #set an impossibly low version if reading the key fails
            #END TRY
            if( TEC_version >= version_filt ):
                TEC_dataAvail[i] = 1; #note data is there
            else:
                print("\n==============~Warning~==============");
                print("The file {} uses version {} while version {} is the most current version. The file will be renamed with _oldV{} and remade with the new version of the algorithm.\n".format(TEC_dataFilePath[i], TEC_version , version_filt, TEC_version) );
                #the filtered version is too old and needs to be redone
                os.rename(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], \
                        settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i]+'_oldV'+str(TEC_version).replace('.','p')); #rename
                try:
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'r', rdcc_nbytes =500*1024*1024) as TEC_fileUnfilt:
                        TEC_unfilt_version = TEC_fileUnfilt.attrs['version']; #read the version
                    #END WITH
                    if( TEC_unfilt_version >= version_unfilt ):
                        TEC_dataAvail[i] = 4; #note that data is already downloaded and converted, but not filtered
                    else:
                        print("\n==============~Warning~==============");
                        print("The file {} uses version {} while version {} is the most current version. The file will be renamed with _oldV{} and remade with the new version of the algorithm.\n".format(TEC_dataFilePathUnfilt[i], TEC_unfilt_version , version_unfilt, TEC_unfilt_version) );
                        os.rename(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], \
                                  settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i]+'_oldV'+str(TEC_unfilt_version).replace('.','p')); #rename
                        #force an OSError to use the rest of the algorithm - uses ***'s which can't be in file names
                        with h5py.File(TEC_dataFilePathUnfilt[i]+'.***notapossiblefile', 'r', rdcc_nbytes =500*1024*1024) as testFile:
                            testFile.keys(); #tries to check out some stuff in the file
                        #END WITH
                    #END IF
                except (OSError , KeyError):
                    TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types
                    try: #gonna try to read the file - if we fail, it failed mid download
                        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r', rdcc_nbytes =500*1024*1024) as testFile:
                            testFile.keys(); #tries to check out some stuff in the file
                        #END WITH
                        TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                        #use case if download success but failure on conversion - big enough files to warrant this care :)
                    except (OSError , KeyError):
                        print("\n==============~Warning~==============");
                        print("Data wasn't downloaded fully for {} and will be redownloaded. Deleting partial file as well.\n".format(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) );
                        os.remove(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete partially downloaded file to avoid conflicts
                        TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
                    #END TRYING
                #END TRYING
            #END IF
        #END IF  
    #END FOR i
    
    #==============Make decisions based on what is downloaded/filtered and what overwrite options are engaged (this was indeed tacked on later)==============
    #Now that that the data has been adjusted, fix it up for days that the padding isn't needed or days that data is arleady there but other days need the original data to be filtered  
    if( FLG_overwrite == 0 ): #with this 0, use whatever is available whenever possible
        
        i = dateRange_dayNum_full[:,0].size-1; #set this so I don't have to change anything in my codes
        #find the last one (it's a padded day) - do this first since we want it done before the middle days are compared
        if( (TEC_dataAvail[i] == 1) & (TEC_dataAvail[i-1] != 1) ): #if a padded day happens to be downloaded, force the raw data to be downloaded as it's only for filtering (as long as the prev day isn't also finished - then day isn't needed)
            TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types  
            
            if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                TEC_dataAvail[i] = 4; #unfiltered data is already downloaded and converted
            elif( os.path.isfile(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == True ): #if this is true, orig data is downloaded but not converted to the faster format
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
                    TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types    
                    
                    if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                        TEC_dataAvail[i] = 4; #unfiltered data is already downloaded and converted
                    elif( os.path.isfile(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == True ): #if this is true, orig data is downloaded but not converted to the faster format
                        TEC_dataAvail[i] = 3; #orig data is downloaded but not converted
                    else: #otherwise get data from download, since data is 1 - assume it exists on the site cause how else would the data have been filtered to start
                        TEC_dataAvail[i] = 2; #orig data needs to be downloaded
                    #END IF
                #END IF
                
                if( TEC_dataAvail[i+1] == 1 ): #next day is done, don't need padded day downloaded
                    TEC_dataAvail[i] = -2; #don't need to import a padded day if the day after is finished
                #END IF
                
            else: #not first or last day, so look forward and back
                if( ( ((TEC_dataAvail[i+1] != 1) & (TEC_dataAvail[i+1] != -2)) | ((TEC_dataAvail[i-1] != 1) & (TEC_dataAvail[i+1] != -2)) ) & (TEC_dataAvail[i] == 1) ):  #if the day before or the day after need to be filtered, and the current day is finished, change the signage
                # if( ( (TEC_dataAvail[i+1] != 1) | (TEC_dataAvail[i-1] != 1) ) & (TEC_dataAvail[i] == 1) ): #if the day before or the day after need to be filtered, and the current day is finished, change the signage
                    TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types    
                    
                    if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                        TEC_dataAvail[i] = 7; #unfiltered data is already downloaded and converted - BUT filtered data is finished so don't filter this day
                    elif( os.path.isfile(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == True ): #if this is true, orig data is downloaded but not converted to the faster format
                        TEC_dataAvail[i] = 6; #orig data is downloaded but not converted - BUT filtered data is finished so don't filter this day
                    else: #otherwise get data from download, since data is 1 - assume it exists on the site cause how else would the data have been filtered to start
                        TEC_dataAvail[i] = 5; #orig data needs to be downloaded - BUT filtered data is finished so don't filter this day
                    #END IF
                #END IF
            #END IF
            
        #END FOR i
        
    elif( FLG_overwrite == 1 ): #FLG_overwrite == 1 flag is on to overwrite all finished data, so just adjust 1's to 2/3/4 depending on what is downloaded
        
        for i in range(0,dateRange_dayNum_full[:,0].size): #loop to check if any data exists - data will not be downloaded yet!
            TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types            
            
            if( os.path.isfile(TEC_dataFilePathUnfilt[i]) == True ): #if this is true, unfiltered data is already downloaded and converted
                TEC_dataAvail[i] = 4; #unfiltered data is already downloaded and converted
            elif( os.path.isfile(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == True ): #if this is true, orig data is downloaded but not converted to the faster format
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
        sys.crash(); #yeet
    #END IF
                
    #==============Download needed data==============
    #--- this is NOT parallelized since the limit is either the Madrigal server output or your internet input, parallelizing will only mob Madrigal essentially ---
    current_webDownloadSpeed = 10; #MB/s, start off guess for download speed
    current_timePerVect = 120; #set initial sec per vector estimate
    tEst_mag = 16; #min, initial estimate for converison to mag
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download needed data
        
        if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == 1): #check for original data requested
            print("\n==============Now Working on Date: {}/{}/{} (Y/M/D) or {}/{} (Y/D#)==============".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2],dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
        else: #else it is a padded day - no filtering will occur
            print("\n==============Now Working on Date: {}/{}/{} (Y/M/D) or {}/{} (Y/D#) - Padded Date, Filter will be Skipped==============".format(dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2],dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
        #END IF
        
        #if( TEC_dataAvail[i] == 1 ): #1 means data is there and we are good
            #pass; #gotta have something in an if statement
            
        #-----DOWNLOAD THE DATA-----
        if( (TEC_dataAvail[i] == 2) | (TEC_dataAvail[i] == 5) ): #2 means data needs to be downloaded (5 means download it for other days, but the current day is filtered and finished already)
            
            #First apply code that makes sure the web stuff is for the correct year
            if( FLG_newYr != dateRange_full[i,0] ): #check to go get the year's web links and file names and file dates
                #basically, if within the same year only need to get this once for the whole year's data availiability
                FLG_newYr = dateRange_full[i,0]; #set it so the flag works good
            
                web_fileSearch = web_base_site + web_base + str(dateRange_dayNum_full[i,0]) + web_baseAfter; #build site to go to
                page = urlopen(web_fileSearch); #get raw HTML
                html_content = page.read(); #read off the HTML from whatever the page holder is
                charset = page.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                html_content = html_content.decode(charset); #"decode" the HTML content so it's legible
                rendered_content = html2text.html2text(html_content); #render the HTML like webpage would and get the real stuff
                #print("{}".format(rendered_content)); #print for debug
                web_fileNamesIndex = strstr(rendered_content,'los_'); #pull out indexes related to the file name and file links ~can't get file name only reliably directly~
                web_fileNotNamesIndex = strstr(rendered_content,'%252Flos_')+5; #pull out indexes related to the end of the file links only ~for removing from above~
                web_fileLinksIndex = strstr(rendered_content,'.h5](/ftp')+5; #pull out indexes related to the ~start~ of the file links only
                web_fileLinksEndIndex = strstr(rendered_content,'.h5/)')+4; #pull out indexes related to the ~end~ of the file links olnly   
            
                web_fileNamesIndex = np.delete( web_fileNamesIndex,np.where(np.in1d(web_fileNamesIndex, web_fileNotNamesIndex))[0]); #search for the not wanted file name ending and remove them
                
                web_fileNames = []; #prep a list for strings, names of the files
                web_fileNamesDate = []; #prep a list for strings, dates of the files
                web_fileLinks = []; #prep a list for strings, links to the files
                for j in range(0,len(web_fileNamesIndex)):
                    web_fileNames.append(rendered_content[web_fileNamesIndex[j]:web_fileNamesIndex[j]+19]); #get the file name recorded
                    web_fileNamesDate.append(rendered_content[web_fileNamesIndex[j]+4:web_fileNamesIndex[j]+12]); #get the file name's date recorded
                    web_fileLinks.append(rendered_content[web_fileLinksIndex[j]:web_fileLinksEndIndex[j]]); #get the file name recorded
                #END FOR
                web_fileNamesDate = np.swapaxes(np.array([np.asarray([w[4:6] for w in web_fileNamesDate]).astype("int16"),np.asarray([w[6:8] for w in web_fileNamesDate]).astype("int16"),np.asarray([w[0:4] for w in web_fileNamesDate]).astype("int16")]),0,1); #pull out the yr/month/day in the file name in M/D/YR format
            #END IF
                    
            current_webFileNames = [web_fileNames[k] for k in np.where(np.all(( web_fileNamesDate - dateRange_full[i,:] )  == 0, axis=1 ))[0]][0]; #get the current web file name, because it's a list it is real oof
            current_web_fileLinks = [web_fileLinks[k] for k in np.where(np.all(( web_fileNamesDate - dateRange_full[i,:] )  == 0, axis=1 ))[0]][0]; #get the current web file link, because it's a list it is real oof
            current_webFileSize = np.asarray(urlopen(web_base_site+current_web_fileLinks).info()['Content-Length'],dtype=np.int64);
            print("Downloading {} to \"{}\".\nFile size is {:.2f} GB. At {:.2f} MB/s ({:.2f} Mbps) expect it to take {:.2f} min.\nFrom site {}".format( current_webFileNames,settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]), current_webFileSize/(1024**3),current_webDownloadSpeed,current_webDownloadSpeed*8.0,(current_webFileSize/1024**2*(1/current_webDownloadSpeed))/60, web_base_site+current_web_fileLinks ));
            web_retryNum = 1; #prep the retry num
            web_retryWait = web_retryWaitSet[0]; #prep the retry wait period
            FLG_gotPage = False; #lets while exit
            while( ((web_maxRetry >= web_retryNum) | (web_maxRetry == -1)) & (FLG_gotPage == False) ):
                try:
                    tic = time(); #for time testing
                    urlretrieve(web_base_site+current_web_fileLinks,settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + current_webFileNames,reporthook=downloadProgress); #download the file in question to the data directory
                    toc = time() - tic; #for time testing
                    FLG_gotPage = True; #exit the loop
                except:
                    print('\n'); #print a new line
                    sys.stdout.write('\rDownload failed for some reason ¯\_(ツ)_/¯. Try #{}/{}, will wait {} sec\t'.format(web_retryNum,web_maxRetry,web_retryWait)); #report
                    sys.stdout.flush();
                    time.sleep(web_retryWait);
                    web_retryNum += 1; #increment try
                    if( web_retryWait < web_retryWaitSet[1] ):
                        web_retryWait += web_retryWaitSet[2]; #increment by amount
                    #END IF
                #END TRYING
            #END WHILE
            current_webDownloadSpeed = current_webFileSize/(1024**2)/toc; #MB/s, calc current web download speed
            print("\nTime to download: {:.2f} min, download speed: {:.2f} MB/s\n\n".format(toc/60,current_webDownloadSpeed)); #extra space at end
            
            if( TEC_dataAvail[i] == 2 ):
                TEC_dataAvail[i] = 3; #move on to the next stage ;) oh yeah it feels wrong but its so easy
            else:
                TEC_dataAvail[i] = 6; #move on to the next stage ;) oh yeah it feels wrong but its so easy - but remembers day is filtered and doesn't need to be redone
            #END IF
        #END IF
    #END FOR i
            
    #==============Convert data to unfiltered (faster!!) format==============
    # establish which days need to to be converted
    instances_toDo = np.zeros(TEC_dataAvail.size, dtype=np.bool_); #prep
    for i in range(0,len(dateRange_dayNum_full[:,0])): #just unpack needed data
        if( (TEC_dataAvail[i] == 3) | (TEC_dataAvail[i] == 6) ): #3 means data was downloaded but needs to be converted (6 means data was downloaded but needs to be converted, but the current day is filtered and finished already)
            instances_toDo[i] = True; # gotta do it
        #END IF
    #END FOR i
    instances_toDo_where = np.where(instances_toDo)[0]; #get where we gotta do
    
    if( instances_toDo_where.size > 0 ):
        import psutil, joblib # these are here not at the top b/c some bad mojo with importing joblib was happening (says it was not imported, local var referenced before assignment)
        #--- estimate the memory needed to filter 1 day ---
        i = instances_toDo_where[0]; #start off with i=first one
        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC), 'r', rdcc_nbytes =500*1024*1024) as TEC_fileOrig:
            #--- deal with table format of file ---
            # List all groups
            TEC_fileKeys = list(TEC_fileOrig.keys()); #get the folder names in the hdf5 file
            TEC_fileKeys_dataIndex = -1; #holder
            for j in range(0,len(TEC_fileKeys)): #get the index of the 'Data' folder
                if( TEC_fileKeys[j] == 'Data' ):
                    TEC_fileKeys_dataIndex = j; #record the index
                #END IF
            #END FOR j
            if( TEC_fileKeys_dataIndex == -1 ): #ERROR CATCH
                print('ERROR: NO \'Data\' folder found in the HDF5 file \"'+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig+'\"\nData structure not supported, exiting. Please update to support cause otherwise whomp whomp :(\n');
                sys.crash();
            #END IF
            TEC_fileDatalocale = list(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]); #choose the data folder
            TEC_fileDatalocale_tableIndex = -1; #holder
            for j in range(0,len(TEC_fileDatalocale)): #get the index of the 'Data' folder
                if( TEC_fileDatalocale[j] == 'Table Layout' ):
                    TEC_fileDatalocale_tableIndex = j; #record the index
                #END IF
            #END FOR j
            if( TEC_fileDatalocale_tableIndex == -1 ): #ERROR CATCH
                print('ERROR: NO \'Data\'Table Layout\' folder found in the HDF5 file \"'+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig+'\"\nData structure not supported, exiting. Please update to support cause otherwise whomp whomp :(\n');
                sys.crash();
            #END IF
            
            #----- Read in the raw data, size it up -----
            TEC_dataSize = np.int64(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]].shape[0]); #get the size
        #END WITH
        floatBytes = np.dtype(dataType_meth).itemsize; #get the bytes used by the floats, which can be 64
        TEC_dataMem = TEC_dataSize*2*8 + TEC_dataSize*floatBytes*6 + TEC_dataSize*8*1 + TEC_dataSize*1*1 + TEC_dataSize*4*1 + TEC_dataSize*1*1; #bytes, valid for current raw format and unfilt 3.0 (characters are 1 byte per character in numpy)
        #7,6 GB used, check
        comp_availMem = np.int64(psutil.virtual_memory().total*.9) #bytes, use 90% of available RAM (.available is also an option to use what's available - but I feel we gotta go hard to do the big calcs)
        
        numInstances_memLim = np.uint32(comp_availMem//TEC_dataMem); # number of instances that can be fit
        
        if( numInstances_memLim < 1 ):
            numInstances_toDo = 1; #just one, no parallel possible due to memory limit
        else:
            if( numInstances_memLim > settings_config['parallel num threads'] ):
                numInstances_toDo = settings_config['parallel num threads']; #enough memory to use all threads
            else:
                numInstances_toDo = numInstances_memLim; #just use as much memory as possible, threads won't be filled out
            #END IF
        #END IF
        if( instances_toDo_where.size < numInstances_toDo ):
            numInstances_toDo = instances_toDo_where.size; #limit to number of days to work on if that's under the computation limit
        #END IF
        
        #-----CONVERT THE DATA TO FASTER FORMAT-----
        if( numInstances_toDo > 1 ):
            try:
                #--- pack up for parallel ---
                print('WARNING in GRITI_import_TEC: Parallelization possible for converting raw data to unfiltered data, which will speed it up. Downside is that no updates will be printed. See you on the other side!');
                tic = time();
                parallel_list = []; #Prep
                for i in range(0,instances_toDo_where.size): # Every iteration appends a list of inputs to parallel_list, each index of parallel_list will be independently run
                    parallel_list.append([ instances_toDo_where[i], TEC_dataAvail, dateRange_full, dateRange_dayNum_full,
                                           TEC_dataFileNameUnfilt, dataType_str, dataType_meth,
                                           locIntUnfilt_size, locFloatUnfilt_size, locStringUnfilt_size, current_timePerVect,
                                           settings_paths, paths_TEC, version_unfilt,
                                           FLG_deleteOrig, True ]);
                #END FOR i
                
                #--- apply parallel calc on function ---
                with joblib.parallel_backend('loky'):
                    with joblib.Parallel(n_jobs=numInstances_toDo,pre_dispatch=numInstances_toDo,batch_size=1) as parallel_arbiter:    
                        parallel_results = parallel_arbiter(joblib.delayed(GRITI_import_TEC_Madrigal_raw2unfilt)(*sublist) for sublist in parallel_list); #annoyingly I need to unpack every variable input manually
                    #END WITH
                #END WITH
                del parallel_list #save some mem
                #--- unpack parallel results ---
                for i in range(0,len(parallel_results)):
                    TEC_dataAvail[instances_toDo_where[i]] = parallel_results[i][0]; #1st one is data avail value, 2nd doesn't matter for parallel
                #END FOR i
                del parallel_results #save some mem
            except Exception as err:
                print('WARNING IN GRITI_import_TEC_Madrigal: ERROR detected in parallelization attempt. Trying non-parallelization. Error message here:\n'+str(err)); #report something went bad, hopefully it was joblib
                tic = time();
                for i in range(0,instances_toDo_where.size):
                     raw2unfilt_results = GRITI_import_TEC_Madrigal_raw2unfilt(instances_toDo_where[i], TEC_dataAvail, dateRange_full, dateRange_dayNum_full,
                                           TEC_dataFileNameUnfilt, dataType_str, dataType_meth,
                                           locIntUnfilt_size, locFloatUnfilt_size, locStringUnfilt_size, current_timePerVect,
                                           settings_paths, paths_TEC, version_unfilt,
                                           FLG_deleteOrig, False);
                     TEC_dataAvail[instances_toDo_where[i]] = raw2unfilt_results[0]; #1st one is new data avail value
                     current_timePerVect = raw2unfilt_results[1]; #2nd one is new current_timePerVect
                #END FOR i
            #END TRYING
            print('\nTime to convert '+str(instances_toDo_where.size)+' days in series took: '+str(np.round((time()-tic)/60,2))+' min\n'); #extra space at end   
        else:
            # no need to parallelize if just 1 instance can run
            for i in range(0,instances_toDo_where.size):
                 raw2unfilt_results = GRITI_import_TEC_Madrigal_raw2unfilt(instances_toDo_where[i], TEC_dataAvail, dateRange_full, dateRange_dayNum_full,
                                       TEC_dataFileNameUnfilt, dataType_str, dataType_meth,
                                       locIntUnfilt_size, locFloatUnfilt_size, locStringUnfilt_size, current_timePerVect,
                                       settings_paths, paths_TEC, version_unfilt,
                                       FLG_deleteOrig, False);
                 TEC_dataAvail[instances_toDo_where[i]] = raw2unfilt_results[0]; #1st one is new data avail value
                 current_timePerVect = raw2unfilt_results[1]; #2nd one is new current_timePerVect
            #END FOR i
        #END IF
    #END IF
       
    #==============Filter Data after download & conversion==============
    TEC_dataAvail = GRITI_import_TEC_support_filter(settings_paths,settings_config,paths_TEC,TEC_dataAvail,TEC_dataFileName,TEC_dataFileNameUnfilt,dateRange_full,dateRange_dayNum_full,dateRange_dayNum_full_orig, \
                                                    TEC_dataRate,TEC_timeTolerance,TEC_maxAmpAllowed,filter_savGolPeriod,order_savGol,minElevation,minimumTimeGap,dataType_meth, \
                                                    TEC_dataFilePathUnfilt, TEC_dataAgg_timeAdditionLimit, TEC_dataAgg_distToPts_degcSq, deltaTEC_compareValue, \
                                                    FLG_deleteUnfilt,FLG_reqPaddedDays,FLG_dataAggregation); #function so that other import_TEC functions can be always up to date
    
    if( FLG_justChecking == True ):
        print('WARNING in GRITI_import_TEC_Madrigal: justChecking mode is ON. Verified the dates are ready, and that\'s it! Future code may not be justChecking safe and be expecting something to come out of this.');
        return {'justChecking':True} # Return empty dict
    #END IF
        
    #==============Read filtered data==============
    TEC_data = {}; #prep the dict
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
        #-----READ THE FILTERED DATA-----
        if( (TEC_dataAvail[i] == 1) | (TEC_dataAvail[i] == 7) ): #1 means data is there, filtered and all - 7 also means that
            keyz = list(TEC_data.keys()); #get the current keys in TEC_data
            file_retryNum = 1; #prep the retry num
            FLG_gotRead = False; #lets while exit
            FLG_appender = False; #appender flag to activate appending magnetic coordinates (append avoided b/c causes h5py bug that prevents reading files until restart !even if python restarted!)
            while( ((file_maxRetry >= file_retryNum) | (file_maxRetry == -1)) & (FLG_gotRead == False) ): #deal with intermittant file read issues I had on an old comp, it tries
                try:
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'r') as TEC_file:
                        #----- Read the data types in -----
                        keyzNew = list(TEC_file.keys()); #get the saved keys
                        #----- Load in lat & long (or mlat & mlong if requested) -----
                        if( coordType == 'geo' ):
                            TEC_dataTemp_lat = TEC_file.get('lat')[()]; #get that dataset out
                            TEC_dataTemp_long = TEC_file.get('long')[()]; #get that dataset out
                        else:
                            if( 'mlat' in keyzNew ):
                                TEC_dataTemp_lat = TEC_file.get('mlat')[()]; #get that dataset out
                                TEC_dataTemp_long = TEC_file.get('mlong')[()]; #get that dataset out
                            else:
                                FLG_appender = True; #calls appender code to do the work
                            #END IF
                        #END IF
                        if( FLG_appender == False ):
                            #--- First prep the mask and import needed vars for it ---
                            keyzNew.remove('lat'); #remove from the list, manually got it for data corralling
                            keyzNew.remove('long'); #remove from the list, manually got it for data corralling
                            if( 'mlat' in keyzNew ):
                                keyzNew.remove('mlat'); #remove from the list, manually got it for data corralling
                                keyzNew.remove('mlong'); #remove from the list, manually got it for data corralling
                            #END IF
                            TEC_goodVals = np.where( (TEC_dataTemp_long <= np.max(plotLongRange)) & (TEC_dataTemp_long >= np.min(plotLongRange)) & \
                                (TEC_dataTemp_lat <= np.max(plotLatRange)) & (TEC_dataTemp_lat >= np.min(plotLatRange)) )[0]; #delete out of lat/long range stuff
                            TEC_dataTemp_lat = TEC_dataTemp_lat[TEC_goodVals]; #keep the good stuff
                            TEC_dataTemp_long = TEC_dataTemp_long[TEC_goodVals]; #keep the good stuff coordType
                            
                            if( strfind(keyz,'lat',opt=1) > 0 ):
                                TEC_data['lat'].append(TEC_dataTemp_lat); #tack that dataset on
                            else:
                                #otherwise it's a new data type to add in
                                TEC_data['lat'] = [np.copy(TEC_dataTemp_lat)]; #get that dataset out
                            #END IF
                            del TEC_dataTemp_lat; #clean up the memory
                            if( strfind(keyz,'long',opt=1) > 0 ):
                                TEC_data['long'].append(TEC_dataTemp_long); #tack that dataset on
                            else:
                                #otherwise it's a new data type to add in
                                TEC_data['long'] = [np.copy(TEC_dataTemp_long)]; #get that dataset out
                            #END IF
                            del TEC_dataTemp_long; #clean up the memory
                            #--- Now import rest of vars and use mask on them ---
                            for j in range(0,len(keyzNew)):
                                if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                                    TEC_data[keyzNew[j]].append(TEC_file.get(keyzNew[j])[()][TEC_goodVals]); #tack that dataset on, keep only the godo stuff
                                else:
                                    #otherwise it's a new data type to add in
                                    TEC_data[keyzNew[j]] = [TEC_file.get(keyzNew[j])[()][TEC_goodVals]]; #get that dataset out, keep only the good stuff
                                #END IF
                            #END FOR j
                            #--- Read the attributes in ---
                            keyzNew = list(TEC_file.attrs.keys()); #get the attribute keys
                            # keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                            for j in range(0,len(keyzNew)):
                                if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                                    TEC_data[keyzNew[j]].append(TEC_file.attrs[keyzNew[j]]); #get that attribute out and tack it on
                                else:
                                    TEC_data[keyzNew[j]] = [TEC_file.attrs[keyzNew[j]]]; #get that attribute out, prep as a list
                                #END IF
                            #END FOR j
                        #END IF
                    #END WITH
                    if( FLG_appender == True ): #append sparingly because it apparently kills h5py's ability to read files if used to just read files (restart req'd)                    
                        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'a') as TEC_file:
                            #----- Read the data types in -----
                            keyzNew = list(TEC_file.keys()); #get the saved keys
                            #----- Load in lat & long (or mlat & mlong if requested) -----
                            if( coordType == 'geo' ):
                                TEC_dataTemp_lat = TEC_file.get('lat')[()]; #get that dataset out
                                TEC_dataTemp_long = TEC_file.get('long')[()]; #get that dataset out
                            else:
                                if( 'mlat' in keyzNew ):
                                    TEC_dataTemp_lat = TEC_file.get('mlat')[()]; #get that dataset out
                                    TEC_dataTemp_long = TEC_file.get('mlong')[()]; #get that dataset out
                                else:
                                    tic_mag = time();
                                    print('Warning in GRITI_import_TEC_Madrigal: coordType "mag" requested but mag coordinates have not been calculated for this day. Calculating and appending onto file '+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i]+\
                                          '. Estimated time is '+textNice(np.round(tEst_mag,2))+' min on a good comp.');
                                    # import aacgmv2
                                    # import datetime
                                    import gc
                                    from multiprocessing import cpu_count #only for cpu_count
                                    import joblib #lets multiprocess happen w/o an insane reloading of GRITI_main
                                    parallel_numCores = cpu_count(); #use multiprocess to get # of CPU cores
                                    alt4mag = TEC_file.attrs['pierceAlt']; #get the altitude for TEC data
                                    TEC_dataTemp_lat = np.copy(TEC_file.get('lat')[()]); #get that dataset out
                                    TEC_dataTemp_long = np.copy(TEC_file.get('long')[()]); #get that dataset out
                                    TEC_dataTemp_time = TEC_file.get('time')[()]; #get that dataset out
                                    TEC_dataTemp_year = TEC_file.get('year')[()][0]; #get that dataset out
                                    TEC_dataTemp_timeUnique = np.unique(TEC_dataTemp_time); #get the time uniques
                                    with joblib.parallel_backend('loky'): #parallel backend uses 
                                        with joblib.Parallel(n_jobs=parallel_numCores,pre_dispatch=parallel_numCores,batch_size=1) as parallel_arbiter: 
                                            #--- Parallel prep datetime object which is apparently hella hard ---
                                            parallel_splitterIndexes = np.int64(np.round(np.linspace(0,TEC_dataTemp_timeUnique.size,parallel_numCores+1,endpoint=True))); #split up the indexes to be parallelized
                                            parallel_list = []; #Prep
                                            for j in range(0,parallel_splitterIndexes.size-1):
                                                parallel_list.append([TEC_dataTemp_time, TEC_dataTemp_timeUnique[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]],  TEC_dataTemp_year]);
                                            #END FOR j
                                            parallel_time4mag = parallel_arbiter(joblib.delayed(calc_datetimes)(j, k, l) for j, k, l in parallel_list); #will this not destroy the world?
                                            time4mag = [None for i in range(0,TEC_dataTemp_time.size)]; #preallocate
                                            for i in range(0,TEC_dataTemp_time.size): #recombine them, they're made of pieces of each (and Nones in between) - since they're objects I am in list purgatory still
                                                for j in range(0,len(parallel_time4mag)):
                                                    if( parallel_time4mag[j][i] != None ):
                                                        time4mag[i] = parallel_time4mag[j][i]; #load it in
                                                    #END IF
                                                #END FOR j
                                            #END FOR i 
                                            del parallel_time4mag
                                            gc.collect(); #make that parallel_time4mag go away
                                            
                                            #--- Parallel calc mag coords ---
                                            parallel_splitterIndexes = np.int64(np.round(np.linspace(0,TEC_dataTemp_lat.size,parallel_numCores+1,endpoint=True))); #split up the indexes to be parallelized
                                            parallel_list = []; #Prep
                                            for j in range(0,parallel_splitterIndexes.size-1):
                                                parallel_list.append([TEC_dataTemp_lat[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]], TEC_dataTemp_long[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]], alt4mag, time4mag[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]]]);
                                            #END FOR j
                                            parallel_mag_list = parallel_arbiter(joblib.delayed(convert_to_mag)(j, k, l, m, method_code='G2A') for j, k, l, m in parallel_list); #will this not destroy the world?
                                            for j in range(0,parallel_splitterIndexes.size-1):
                                                TEC_dataTemp_lat[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_mag_list[j][0].astype(dataType_meth().dtype, order='C', casting='same_kind', subok=True, copy=False); #load it in
                                                TEC_dataTemp_long[parallel_splitterIndexes[j]:parallel_splitterIndexes[j+1]] = parallel_mag_list[j][1].astype(dataType_meth().dtype, order='C', casting='same_kind', subok=True, copy=False); #load it in
                                            #END FOR j
                                        #END WITH
                                    #END WITH
                                    del parallel_list, parallel_mag_list #save some mem 
                                    # Save those newly calc'd datas
                                    TEC_file.create_dataset('mlat', data=TEC_dataTemp_lat, chunks=TEC_dataTemp_lat.shape, compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #write that data
                                    TEC_file.create_dataset('mlong', data=TEC_dataTemp_long, chunks=TEC_dataTemp_lat.shape, compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #write that data
                                    tEst_mag = (time() - tic_mag)/60; #record
                                    print('Conversion to coordType "mag" & appending to file took '+textNice(np.round(tEst_mag,2))+' min.');
                                #END IF
                            #END IF
                            #--- First prep the mask and import needed vars for it ---
                            keyzNew.remove('lat'); #remove from the list, manually got it for data corralling
                            keyzNew.remove('long'); #remove from the list, manually got it for data corralling
                            if( 'mlat' in keyzNew ):
                                keyzNew.remove('mlat'); #remove from the list, manually got it for data corralling
                                keyzNew.remove('mlong'); #remove from the list, manually got it for data corralling
                            #END IF
                            TEC_goodVals = np.where( (TEC_dataTemp_long <= np.max(plotLongRange)) & (TEC_dataTemp_long >= np.min(plotLongRange)) & \
                                (TEC_dataTemp_lat <= np.max(plotLatRange)) & (TEC_dataTemp_lat >= np.min(plotLatRange)) )[0]; #delete out of lat/long range stuff
                            TEC_dataTemp_lat = TEC_dataTemp_lat[TEC_goodVals]; #keep the good stuff
                            TEC_dataTemp_long = TEC_dataTemp_long[TEC_goodVals]; #keep the good stuff coordType
                            
                            if( strfind(keyz,'lat',opt=1) > 0 ):
                                TEC_data['lat'].append(TEC_dataTemp_lat); #tack that dataset on
                            else:
                                #otherwise it's a new data type to add in
                                TEC_data['lat'] = [np.copy(TEC_dataTemp_lat)]; #get that dataset out
                            #END IF
                            del TEC_dataTemp_lat; #clean up the memory
                            if( strfind(keyz,'long',opt=1) > 0 ):
                                TEC_data['long'].append(TEC_dataTemp_long); #tack that dataset on
                            else:
                                #otherwise it's a new data type to add in
                                TEC_data['long'] = [np.copy(TEC_dataTemp_long)]; #get that dataset out
                            #END IF
                            del TEC_dataTemp_long; #clean up the memory
                            #--- Now import rest of vars and use mask on them ---
                            for j in range(0,len(keyzNew)):
                                if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                                    TEC_data[keyzNew[j]].append(TEC_file.get(keyzNew[j])[()][TEC_goodVals]); #tack that dataset on, keep only the godo stuff
                                else:
                                    #otherwise it's a new data type to add in
                                    TEC_data[keyzNew[j]] = [TEC_file.get(keyzNew[j])[()][TEC_goodVals]]; #get that dataset out, keep only the good stuff
                                #END IF
                            #END FOR j
                            #--- Read the attributes in ---
                            keyzNew = list(TEC_file.attrs.keys()); #get the attribute keys
                            # keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                            for j in range(0,len(keyzNew)):
                                if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                                    TEC_data[keyzNew[j]].append(TEC_file.attrs[keyzNew[j]]); #get that attribute out and tack it on
                                else:
                                    TEC_data[keyzNew[j]] = [TEC_file.attrs[keyzNew[j]]]; #get that attribute out, prep as a list
                                #END IF
                            #END FOR j
                        #END WITH
                    #END IF
                    FLG_gotRead = True; #lets while loop exit
                except OSError as errorText:
                    print('Warning in GRITI_import_TEC_Madrigal: '+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i] +\
                          ' had an OSError when reading data. Try #'+str(file_retryNum)+'/'+str(file_maxRetry)+'. Error text follows:');
                    print(str(errorText));
                    time.sleep(0.1); #wait a tiny lil bit just in case
                    file_retryNum += 1; #increment try
                #END TRY
            #END WHILE
            if( FLG_gotRead == False ):
                print("\n==============ERROR==============");
                print(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i] +\
                      ' failed to read. Renaming the file to "'+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i]+'_corrupted" and crashing. Rerun code to generate new file.');
                os.rename(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], \
                          settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i]+'_corrupted'); #rename
                print('NOTE this may be an OS read error that can be fixed with a system restart. Not sure why it happens but you\'ll need to rename it back from _corrupted and restart your comp. Restarting python instance doesn\'t seem to cut it.');
                sys.crash(); #donezo
            #END IF
        #END IF    
    #END FOR i
    #----- Combine lists into numpy arrays -----
    keyz = list(TEC_data.keys()); #get the current keys
    for j in range(0,len(keyz)):
        if( np.isscalar(TEC_data[keyz[j]][0]) == False ):
            #if not a scalar, apply the logical mask
            TEC_data[keyz[j]] = np.hstack(TEC_data[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
        # else:
        #     tempHolder = np.empty(len(TEC_data[keyz[j]]),dtype=type(TEC_data[keyz[j]])); #prep
        #     for k in range(0,len(TEC_data[keyz[j]])):
        #         tempHolder[k] = TEC_data[keyz[j]][k]; #get them all
        #     #END FOR k
        #     if( np.all(tempHolder == tempHolder[0]) ):
        #         TEC_data[keyz[j]] = TEC_data[keyz[j]][0]; #get 1st one since they're all the same
        #     else:
        #         print('\n==============~Warning~==============');
        #         print('Not all '+keyz[j]+' reported are the same - I don\'t have anything to deal with this. Continuing and using 1st one as only one. They should print below:');
        #         print('{}'.format(TEC_data[keyz[j]]));
        #         TEC_data[keyz[j]] = TEC_data[keyz[j]][0]; #get 1st one even though they're not all the same whoopsie
        #     #END IF
        #END IF
    #END FOR j
    # Calc total sec for the entire time period aligned to the zero hr [debuting year support! I hope it works]
    TEC_data['time aligned'] = np.int32(TEC_data['dayNum']-dateRange_dayNum_zeroHr[1]-(dateRange_dayNum_zeroHr[0]-TEC_data['year'])*subfun_daysInAYear(TEC_data['year']))*86400 + \
        np.int32(TEC_data['hour'])*3600 + np.int32(TEC_data['min'])*60 + np.int32(TEC_data['sec']); #sec, calc total sec for the day range
    # Record the total data rate
    tempHolder = np.empty(len(TEC_data['forcedTECdataRateSec']),dtype=type(TEC_data['forcedTECdataRateSec'])); #prep
    for k in range(0,len(TEC_data['forcedTECdataRateSec'])):
        tempHolder[k] = TEC_data['forcedTECdataRateSec'][k]; #get them all
    #END FOR k
    if( np.all(tempHolder == tempHolder[0]) ):
        TEC_data['data rate'] = TEC_data['forcedTECdataRateSec'][0]; #get 1st one since they're all the same
    else:
        print('\n==============~Warning~==============');
        print('Not all '+'forcedTECdataRateSec'+' reported are the same - I don\'t have anything to deal with this. Continuing and using 1st one as only one. They should print below:');
        print('{}'.format(TEC_data['forcedTECdataRateSec']));
        TEC_data['data rate'] = TEC_data['forcedTECdataRateSec'][0]; #get 1st one even though they're not all the same whoopsie
    #END IF
    
    return TEC_data
#END DEF

def get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC):
    #los_20140227.001.h5 - example file name [001 means final data, 00p is custom for preliminary data]
    TEC_fileNameOrig_format = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.*.h5'; #create expected file name format
    TEC_fileNameOrig_avail = glob.glob(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig_format); #get available raw data
    if( len(TEC_fileNameOrig_avail) == 1 ):
        TEC_fileNameOrig = TEC_fileNameOrig_avail[0][TEC_fileNameOrig_avail[0].rfind('\\')+1:]; #remove the list
    elif( len(TEC_fileNameOrig_avail) > 1 ):
        if( strfind(TEC_fileNameOrig_avail,'.001.',opt=1) == 1 ): #preferred first
            TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name
        elif( strfind(TEC_fileNameOrig_avail,'.00p.',opt=1) == 1 ): #other types of data in order of preference
            TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.00p.h5'; #create expected file name
        #END IF
    else:
        #if no match send out the default 001
        TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name
    #END IF
    return TEC_fileNameOrig
#END DEF

def convert_to_mag(lat4mag, long4mag, alt4mag, time4mag, method_code='G2A'):
    import aacgmv2
    temp_lat = np.empty(lat4mag.size,dtype=lat4mag.dtype); #preallocate
    temp_long = np.empty(lat4mag.size,dtype=lat4mag.dtype); #preallocate
    for jk in range(0,lat4mag.size):
        temp_lat[jk], temp_long[jk], _ = aacgmv2.convert_latlon(lat4mag[jk], long4mag[jk], alt4mag, time4mag[jk], method_code=method_code); #converts from geographic to geomagnetic (AACGMv2) - vectorized _arr version has a memory leak atm
        tryCntr = 1; #reset cntr
        while( (np.isnan(temp_long[jk]) | np.isnan(temp_lat[jk])) & (tryCntr <= 20) ):
            #recalc at higher altitude to deal with error (idk how bad this is, but it allows the points to exist at least
            temp_lat[jk], temp_long[jk], _ = aacgmv2.convert_latlon(lat4mag[jk], long4mag[jk], alt4mag+200*tryCntr, time4mag[jk], method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
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

def GRITI_import_TEC_Madrigal_raw2unfilt(i, TEC_dataAvail, dateRange_full, dateRange_dayNum_full,
                                         TEC_dataFileNameUnfilt, dataType_str, dataType_meth,
                                         locIntUnfilt_size, locFloatUnfilt_size, locStringUnfilt_size, current_timePerVect,
                                         settings_paths, paths_TEC, version_unfilt,
                                         FLG_deleteOrig, FLG_parallel):
    #this as a function supports parallelizing (which is really needed since the Madrigal format data comes in is hamstrung by poor chunk sizing so nothing is maxed out while reading except one thread or something) ---
    #-----CONVERT THE DATA TO FASTER FORMAT-----
    if( (TEC_dataAvail[i] == 3) | (TEC_dataAvail[i] == 6) ): #3 means data was downloaded but needs to be converted (6 means data was downloaded but needs to be converted, but the current day is filtered and finished already)
        
        TEC_fileNameOrig = get_TEC_fileNameOrig(dateRange_full, i, settings_paths, paths_TEC); #call function that deals with various raw data names and which are preferred if there's multiple types            
        
        dataVectMax = locIntUnfilt_size+locFloatUnfilt_size+locStringUnfilt_size+2; #total number of data vectors to import, order: int / float / string
        if( FLG_parallel == False ):
            print("Converting {} to {}.\nAt {} data vectors and {} sec per vector, expect {} minutes for conversion to finish.\n".format(TEC_fileNameOrig,TEC_dataFileNameUnfilt[i],dataVectMax,current_timePerVect,textNice(np.round(current_timePerVect*dataVectMax/60,2))));
            tic = time(); #for time testing
        #END IF
        
        unfilt_dict = {
            'version':version_unfilt, #attributes here, code automagically saves them too
            };
                    
        #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
        if( '.001.' in TEC_fileNameOrig ):
            # file reader for 001 default Madrigal format
            with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r', rdcc_nbytes =500*1024*1024) as TEC_fileOrig:
                #TEC_fileOrig = h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r'); #read that file
            
                # List all groups
                TEC_fileKeys = list(TEC_fileOrig.keys()); #get the folder names in the hdf5 file
                TEC_fileKeys_dataIndex = -1; #holder
                for j in range(0,len(TEC_fileKeys)): #get the index of the 'Data' folder
                    if( TEC_fileKeys[j] == 'Data' ):
                        TEC_fileKeys_dataIndex = j; #record the index
                    #END IF
                #END FOR j
                if( TEC_fileKeys_dataIndex == -1 ): #ERROR CATCH
                    print('ERROR: NO \'Data\' folder found in the HDF5 file \"'+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig+'\"\nData structure not supported, exiting. Please update to support cause otherwise whomp whomp :(\n');
                    sys.crash();
                #END IF
                TEC_fileDatalocale = list(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]); #choose the data folder
                TEC_fileDatalocale_tableIndex = -1; #holder
                for j in range(0,len(TEC_fileDatalocale)): #get the index of the 'Data' folder
                    if( TEC_fileDatalocale[j] == 'Table Layout' ):
                        TEC_fileDatalocale_tableIndex = j; #record the index
                    #END IF
                #END FOR j
                if( TEC_fileDatalocale_tableIndex == -1 ): #ERROR CATCH
                    print('ERROR: NO \'Data\'Table Layout\' folder found in the HDF5 file \"'+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig+'\"\nData structure not supported, exiting. Please update to support cause otherwise whomp whomp :(\n');
                    sys.crash();
                #END IF
                # TEC_fileData = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]]; #import the data shape
                dataVectCntr = 0; #counter to announce progress
                unfilt_dict['sat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'sat_id'].astype("int16"); #import satellite ID
                # unfilt_dict['sat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]].fields('sat_id')[()].astype("int16"); #import satellite ID
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                if( 'gnss_type' in list(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]].dtype.fields.keys()) ): #only read if it is there
                    # satType_temp = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'gnss_type']; #import satellite type
                    unfilt_dict['satType'] = satConv_toStr_fromLongStr(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]].fields('gnss_type')[:]); #convert long string name to single character name
                    # satType_temp = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]].fields('gnss_type')[:]; #import satellite type
                    # satType_tempLen = np.int64(satType_temp.dtype.str[satType_temp.dtype.str.find('|S')+2:]); #get the string length
                    # unfilt_dict['satType'] = np.empty(unfilt_dict['sat'].size, dtype='S1'); #preallocate a single character
                    # unfilt_dict['satType'][satType_temp == np.string_('GPS'.ljust(satType_tempLen))] = 'G'; #record GPS
                    # unfilt_dict['satType'][satType_temp == np.string_('GLONASS'.ljust(satType_tempLen))] = 'R'; #record GLONASS
                    # unfilt_dict['satType'][satType_temp == np.string_('SBAS'.ljust(satType_tempLen))] = 'S'; #record SBAS Payload
                    # unfilt_dict['satType'][satType_temp == np.string_('GALILEO'.ljust(satType_tempLen))] = 'E'; #record Galileo
                    # unfilt_dict['satType'][satType_temp == np.string_('QZSS'.ljust(satType_tempLen))] = 'J'; #record QZSS
                    # unfilt_dict['satType'][satType_temp == np.string_('BEIDOU'.ljust(satType_tempLen))] = 'C'; #record BeiDou
                    # unfilt_dict['satType'][satType_temp == np.string_('NAVIC'.ljust(satType_tempLen))] = 'I'; #record NavIC
                    # del satType_temp; #clean the mem
                else:
                    unfilt_dict['satType'] = np.tile(np.string_('G'),unfilt_dict['sat'].size); #assume all GPS satellite type
                #END TRY
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector satType {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                unfilt_dict['year'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'year'].astype("int16"); #import year
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector year {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                tempMonth = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'month'].astype("int16"); #import month
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector month {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                tempDay = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'day'].astype("int16"); #import day
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector day {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                unfilt_dict['dayNum'] = subfun_date_to_dayNum( np.array( (unfilt_dict['year'],tempMonth,tempDay) ) , 2 ); #convert to day number, option 2 only returns day number
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector dayNum calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                unfilt_dict['hour'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'hour'].astype("int16"); #import hour
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector hour {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                unfilt_dict['min'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'min'].astype("int16"); #import min
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector min {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                unfilt_dict['sec'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'sec'].astype("int16"); #import sec
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sec {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF

                unfilt_dict['lat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'gdlat'].astype(dataType_str); #arcdeg, import pierce point geodedic latitude
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector lat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                unfilt_dict['long'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'glon'].astype(dataType_str); #arcdeg, import pierce point longitude
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector long {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                unfilt_dict['elev'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'elm'].astype(dataType_str); #deg, import pierce point elevation
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector elev {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                unfilt_dict['sTEC'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'los_tec'].astype(dataType_str); #TECU, import pierce point line-of-sight TEC
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sTEC {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                unfilt_dict['sTECerr'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'dlos_tec'].astype(dataType_str); #TECU, import pierce point line-of-sight TEC error
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sTECerr {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                unfilt_dict['vTEC'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'tec'].astype(dataType_str); #TECU, import adjusted vertical TEC (as recommended by AC)
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector vTEC {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
    #                TEC_fileData_float[:,locFloatUnfilt_TECerr] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'dlos_tec']; #TECU, import pierce point line-of-sight TEC error
    #                dataVectCntr = dataVectCntr + 1;
    #                sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
    #                sys.stdout.flush();
                #TEC_fileData_float[:,8] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'rec_bias']; #TECU, import receiver bias
                #TEC_fileData_float[:,9] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'drec_bias']; #TECU, import receiver bias error
                
                unfilt_dict['pierceAlt'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'pierce_alt'].astype(dataType_str); #km, import pierce point altitude
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector pierceAlt {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                #pierce point altitude is the same number, no need to record it a bunch
                            
                unfilt_dict['site'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'gps_site']; #import GPS site name            
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector site {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
            #END WITH
            del tempMonth, tempDay
        elif( '.00p.' in TEC_fileNameOrig ):
            # file reader for preliminary Madrigal format
            with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r', rdcc_nbytes =500*1024*1024) as TEC_fileOrig:
                #TEC_fileOrig = h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r'); #read that file
            
                # List all groups
                TEC_fileKeys = list(TEC_fileOrig.keys()); #get the folder names in the hdf5 file
                TEC_fileKeys_dataIndex = -1; #holder
                for j in range(0,len(TEC_fileKeys)): #get the index of the 'Data' folder
                    if( TEC_fileKeys[j] == 'data' ):
                        TEC_fileKeys_dataIndex = j; #record the index
                    #END IF
                #END FOR j
                if( TEC_fileKeys_dataIndex == -1 ): #ERROR CATCH
                    print('ERROR: NO \'data\' folder found in the HDF5 file \"'+settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig+'\"\nData structure not supported, exiting. Please update to support cause otherwise whomp whomp :(\n');
                    sys.crash();
                #END IF
                # TEC_fileData = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]]; #import the data shape
                dataVectCntr = 0; #counter to announce progress
                             
                # time_temp = np.mod(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['time'].astype("float64"), dateRange_dayNum_full[i,1]-1); #import time float [day 14 has times from 13 to 14 b/c I assume so that Jan 1st can start at 0]
                time_temp = np.mod(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['time'], dateRange_dayNum_full[i,1]-1); #import time float [day 14 has times from 13 to 14 b/c I assume so that Jan 1st can start at 0]
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector time[temp] {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                time_temp *= 24; #hourify
                unfilt_dict['hour'] = np.int16(time_temp); #get the hour
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector hour calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                time_temp *= 60; #minutify
                unfilt_dict['min'] = np.int16(time_temp - unfilt_dict['hour']*60); #get the minute
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector min calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                time_temp *= 60; #secondify
                unfilt_dict['sec'] = np.int16(time_temp - unfilt_dict['min']*60 - unfilt_dict['hour']*3600); #get the second
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sec calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                del time_temp; #clean mem  
                unfilt_dict['year'] = np.ones(unfilt_dict['sec'].size ,dtype=np.int16)*dateRange_dayNum_full[i,0];
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector year calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                unfilt_dict['dayNum'] = np.ones(unfilt_dict['sec'].size ,dtype=np.int16)*dateRange_dayNum_full[i,1];
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector dayNum calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                
                # unfilt_dict['sat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['prn'].astype("int16"); #import satellite ID
                unfilt_dict['sat'] = np.int16(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['prn']); #import satellite ID
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();  
                #END IF
                unfilt_dict['satType'] = satConv_toStr_fromInt(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['sat_type']); #convert integer sat type to single character sat type
                # satType_temp = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['sat_type']; #import satellite type (0 for GPS, 1 for GLONASS)
                # unfilt_dict['satType'] = np.empty(satType_temp.size, dtype='S1'); #preallocate a single character
                # unfilt_dict['satType'][satType_temp == 0] = 'G'; #record GPS
                # unfilt_dict['satType'][satType_temp == 1] = 'R'; #record GLONASS
                # del satType_temp; #clean mem                    
                dataVectCntr += 1; #increment
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2))); #report
                    sys.stdout.flush();
                #END IF
                
                # unfilt_dict['lat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['pplat'].astype(dataType_str); #arcdeg, import pierce point geodedic latitude
                unfilt_dict['lat'] = dataType_meth(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['pplat']); #arcdeg, import pierce point geodedic latitude
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector lat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                # unfilt_dict['long'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['pplon'].astype(dataType_str); #arcdeg, import pierce point longitude
                unfilt_dict['long'] = dataType_meth(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['pplon']); #arcdeg, import pierce point longitude
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector long {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                # unfilt_dict['elev'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['el'].astype(dataType_str); #deg, import pierce point elevation
                unfilt_dict['elev'] = dataType_meth(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['el']); #deg, import pierce point elevation
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector elev {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                # unfilt_dict['sTEC'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['los_tec'].astype(dataType_str); #TECU, import pierce point line-of-sight TEC
                unfilt_dict['sTEC'] = dataType_meth(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['los_tec']); #TECU, import pierce point line-of-sight TEC
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sTEC {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                # unfilt_dict['sTECerr'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['err_los_tec'].astype(dataType_str); #TECU, import pierce point line-of-sight TEC error
                unfilt_dict['sTECerr'] = dataType_meth(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['err_los_tec']); #TECU, import pierce point line-of-sight TEC error
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector sTECerr {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                # unfilt_dict['vTEC'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['vtec'].astype(dataType_str); #TECU, import adjusted vertical TEC (as recommended by AC)
                unfilt_dict['vTEC'] = dataType_meth(TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['vtec']); #TECU, import adjusted vertical TEC (as recommended by AC)
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector vTEC {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
    #                TEC_fileData_float[:,locFloatUnfilt_TECerr] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]][:,'dlos_tec']; #TECU, import pierce point line-of-sight TEC error
    #                dataVectCntr = dataVectCntr + 1;
    #                sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
    #                sys.stdout.flush();
                #TEC_fileData_float[:,8] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]][:,'rec_bias']; #TECU, import receiver bias
                #TEC_fileData_float[:,9] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]][:,'drec_bias']; #TECU, import receiver bias error
                
                unfilt_dict['pierceAlt'] = np.ones(unfilt_dict['sec'].size ,dtype=dataType_str)*350.; #km, import pierce point altitude [assume 350 since not given]
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector pierceAlt {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
                #pierce point altitude is the same number, no need to record it a bunch
                            
                unfilt_dict['site'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex]]['site']; #import GPS site name            
                dataVectCntr = dataVectCntr + 1;
                if( FLG_parallel == False ):
                    sys.stdout.write("\rData vector site {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time()-tic)/60,2)));
                    sys.stdout.flush();
                #END IF
            #END WITH
        #END IF
        
        #-----make sure hour 24 or min 60 doesn't exist (they hsould be rolled over) - found issue via day 129 data showing up b/c it was recorded as d128/h24/m60...------       
        #---CHECK SECONDS---
        kj = unfilt_dict['sec'] >= 60; #find incorrect time keeping
        if( np.sum(kj) > 0 ):
            #increment the main minute variable where seconds are 60 or more
            unfilt_dict['min'][kj] += 1; #min, increment time by 1
            unfilt_dict['sec'][kj] -= 60; #sec, remove 60 time units from the time keeping
            if( np.any(unfilt_dict['sec'] >= 60) ):
                print('ERROR: TIME KEEPING IS REALLY BAD, 60+ SECOND TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ SECOND TIME STAMPS');
                sys.crash(); #call a crash (crash doesn't exist so it crashes)
            #END IF
        #END IF
        #---CHECK MINUTES---
        kj = unfilt_dict['min'] >= 60; #find incorrect time keeping
        if( np.sum(kj) > 0 ):
            #increment the main hour variable where minutes are 60 or more
            unfilt_dict['hour'][kj] += 1; #hour, increment time by 1
            unfilt_dict['min'][kj] -= 60; #min, remove 60 time units from the time keeping
            if( np.any(unfilt_dict['min'] >= 60) ):
                print('ERROR: TIME KEEPING IS REALLY BAD, 60+ MINUTE TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ MIN TIME STAMPS');
                sys.crash(); #call a crash (crash doesn't exist so it crashes)
            #END IF
        #END IF
        #---CHECK HOURS---
        kj = unfilt_dict['hour'] >= 24; #find incorrect time keeping
        if( np.sum(kj) > 0 ):
            #increment the main day variable where hours are 24 or more
            unfilt_dict['dayNum'][kj] += 1; #day number, increment time by 1
            unfilt_dict['hour'][kj] -= 24; #hour, remove 24 time units from the time keeping
            if( np.any(unfilt_dict['hour'] >= 24) ):
                print('ERROR: TIME KEEPING IS REALLY BAD, 24+ HOUR TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 48+ HOUR TIME STAMPS');
                sys.crash(); #call a crash (crash doesn't exist so it crashes)
            #END IF
        #END IF
        #---CHECK DAYS---
        #deal with day limit is based on leap year or not
        dayLim = np.ones(unfilt_dict['dayNum'].shape,dtype=np.int16)*365; #day number, get the day number limits as 365
        #adjust leap years to 366
        leapYears = (np.mod(unfilt_dict['year'],4) == 0) & (np.mod(unfilt_dict['year'],100) != 0) & (np.mod(unfilt_dict['year'],400) == 0)
        dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
        kj = unfilt_dict['dayNum'] >= dayLim; #find incorrect time keeping
        if( np.sum(kj) > 0 ):
            #increment the main year variable where day number is equal to the day number limit or higher than it
            unfilt_dict['year'][kj] += 1; #year, increment time by 1
            unfilt_dict['dayNum'][kj] -= dayLim; #days, remove the day number limit time units from the time keeping
            if( np.any(unfilt_dict['dayNum'] >= dayLim) ):
                print('ERROR: TIME KEEPING IS REALLY BAD, 365/366 DAY TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 365*2/366*2 DAY TIME STAMPS');
                sys.crash(); #call a crash (crash doesn't exist so it crashes)
            #END IF
        #END IF
        del kj
        #---luckily years have no limits (that we know of??)---
        #---DELETE data not on the right day (could have incrememted past the right day due to stuff and things like a 24 hour/60 minute time on the current day---
        keyz = list(unfilt_dict.keys()); #get the keys
        TEC_fileData_logical_onDay = (unfilt_dict['dayNum'] == dateRange_dayNum_full[i,1]) & (unfilt_dict['year'] == dateRange_dayNum_full[i,0]); #find when the day reported is the day we want [and the year we want]
        for j in range(0,len(keyz)):
            if( np.isscalar(unfilt_dict[keyz[j]]) == False ):
                unfilt_dict[keyz[j]] = unfilt_dict[keyz[j]][TEC_fileData_logical_onDay]; #keep only the good stuff
            #END IF
        #END FOR j
        del TEC_fileData_logical_onDay
        #---UPDATE total time variable---
        unfilt_dict['time'] = np.int64(unfilt_dict['hour'])*3600 + np.int64(unfilt_dict['min'])*60 + np.int64(unfilt_dict['sec']); #sec, calculate hour/min/sec into days and add to the current day np.int64(unfilt_dict['dayNum'])*86400 + 
        
        #make sure pierce altitude is consistent
        if( np.isscalar(unfilt_dict['pierceAlt']) == 0 ):
            if( np.all( unfilt_dict['pierceAlt'] == unfilt_dict['pierceAlt'][0]) ):
                unfilt_dict['pierceAlt'] = unfilt_dict['pierceAlt'][0]; #record the attribute
            else:
                print("\n==============~Warning~==============");
                print("Not all pierce point altitudes reported are the same - I don't have anything to deal with this. Continuing. They should print below:");
                print("{}".format(unfilt_dict['pierceAlt']));
                unfilt_dict['pierceAlt'] = unfilt_dict['pierceAlt'][0]; #record the attribute
            #END IF
        #END IF
        if( FLG_parallel == False ):
            print('\nDone reading data from original file.');
        #END IF
        
        keyz = list(unfilt_dict.keys()); #get the keys again
        #--- Get chunk size (it's the size of one of the vectors, since we read in entire vectors at a time) ---
        for j in range(0,len(keyz)):
            if( np.isscalar(unfilt_dict[keyz[j]]) == False ):
                h5pyChunkShape = unfilt_dict[keyz[j]].shape; #get the chunk size of a vector
                break;
            #END IF
        #END FOR j
        #--- Write in the unfiltered data in a new faster format ---
        #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'w', rdcc_nbytes =500*1024*1024) as TEC_fileUnfilt:               
            for j in range(0,len(keyz)):
                if( np.isscalar(unfilt_dict[keyz[j]]) == False ):
                    TEC_fileUnfilt.create_dataset(keyz[j], data=unfilt_dict[keyz[j]], chunks=h5pyChunkShape, compression='gzip', compression_opts=9, shuffle=True, fletcher32=True); #write that data , compression="gzip"
                else:
                    #if size 1, add it as an attribute
                    TEC_fileUnfilt.attrs[keyz[j]] = unfilt_dict[keyz[j]]; #save the attribute
                #END IF
                if( FLG_parallel == False ):
                    sys.stdout.write("\rWriting {} to file & {} min | {} out of {}\t\t\t\t\t".format(keyz[j],np.round((time()-tic)/60,2),j+1,len(keyz)));
                    sys.stdout.flush();
                #END IF
            #END FOR j
            if( FLG_parallel == False ):
                print('\nDone writing unfiltered file!');
            #END IF
        #END WITH
        if( FLG_deleteOrig == 1 ): #only delete if flag is on
            os.remove(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete the downloaded file now that we're done with it (save some space!)
        #END IF
        del unfilt_dict; #clean the memory
        if( FLG_parallel == False ):
            toc = time() - tic; #for time testing
            current_timePerVect = np.round(toc/dataVectMax,2); #update the time estiamte (includes writing files and such)
            print("\nTime to convert: {} min\n".format(np.round(toc/60,2))); #extra space at end   
        #END IF
        
        if( TEC_dataAvail[i] == 3 ):
            TEC_dataAvail[i] = 4; #move on to the next stage
        else:
            TEC_dataAvail[i] = 7; #move on to the next stage - but remembers day is filtered and doesn't need to be redone
        #END IF
    #END IF
    
    return [TEC_dataAvail[i], current_timePerVect]
#END DEF