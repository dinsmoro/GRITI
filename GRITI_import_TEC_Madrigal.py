#Function to import TEC data from Madrigal3
#RD on 8/23/2018
#
#

import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import os
import time
import sys
from urllib.request import urlopen, urlretrieve
import html2text
from subfun_strstr import strstr
from subfun_strstrNB import strstrNB
from subfun_downloadProgress import downloadProgress
import h5py
from subfun_date_to_dayNum import subfun_date_to_dayNum
from subfun_addADay import subfun_addADay
from subfun_dayNum_to_date import subfun_dayNum_to_date
from subfun_daysInAYear import subfun_daysInAYear
from subfun_strfind import strfind
from scipy import signal
#from numba import jit, prange
#-----Testing variables-----
##http://cedar.openmadrigal.org/ftp/fullname/Ross+Dinsmore/email/rld5204@psu.edu/affiliation/PSU/kinst/8000/year/2015/kindat/3505/format/hdf5/
#import matplotlib.pyplot as plt
#from matplotlib import _pylab_helpers
#from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
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


def GRITI_import_TEC_Madrigal(dates, settings):
    #==============Version Declaration==============
    version_unfilt = 2.0; #unfiltered algorithm version
    #changelog:
    #1 10/8/2019 - initial algorithm
    #1.1 9/14/2020 - fixed hour 24/minute 60 handling (before ignored, now adjusted to +1 day/+1 hour)
    #2.0 4/9/2021 - moved from monolithic data bloc to slices of data
    
    version_filt = 3.0; #filtered algorithm version
    #1 10/8/2019 - initial algorithm
    #1.1 9/11/2020 - fixed non-0/30 time step handling (29/59 were big ones)    
    #1.2 9/16/2020 - removed highpass filter, only dampens data
    #1.3 9/17/2020 - outlier control introduced
    #1.4 9/17/2020 - linear interpolation to make savgol filter work correctly for gappy data
    #2.0 9/17/2020 - interpolated data used (large groups of receivers would not have data at the same time, leaving vertical "streaks" of no-data)
    #3.0 4/9/2021 - moved from monolithic data bloc to slices of data, moved to dict royale
    
    #==============Unpack==============
    dateRange_dayNum_full = dates['date range full dayNum']; #unpack
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum']; #unpack
    settings_paths = settings['paths']; #unpack
    plotLatRange = settings['map']['lat range']; #unpack
    plotLongRange = settings['map']['long range']; #unpack
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
    locStringUnfilt_size = 1; #size of string layout
    
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
            #los_20140227.001.h5 - example file name
            TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name 
            
            if( (os.path.isfile(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == 0) and (os.path.isfile(TEC_dataFilePathUnfilt[i]) == 0) ):
                #if the unconverted file is not there, download it
                #first - need to check if TEC data is availiable for the days requested
                
                if( FLG_newYr != dateRange_full[i,0] ): #check to go get the year's web links and file names and file dates
                    #basically, if within the same year only need to get this once for the whole year's data availiability
                    FLG_newYr = dateRange_full[i,0]; #set it so the flag works good
                
                    web_fileSearch = web_base_site + web_base + str(dateRange_dayNum_full[i,0]) + web_baseAfter; #build site to go to
                    page = urlopen(web_fileSearch); #get raw HTML
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
                #END IF
                
                #Check if website has data for the date requested
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
            elif( os.path.isfile(TEC_dataFilePathUnfilt[i]) == 1 ): #check if unfiltered data file exists            
                try: #gonna try to read the file - if we fail, it failed mid download probably
                    with h5py.File(TEC_dataFilePathUnfilt[i], 'r') as testFile:
                        testFile.keys(); #tries to check out some stuff in the file
                    #END WITH
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'r') as TEC_fileUnfilt:
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
                                  settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i]+'_oldV'+str(TEC_unfilt_version)); #rename
                        #force an OSError to use the rest of the algorithm - uses ***'s which can't be in file names
                        with h5py.File(TEC_dataFilePathUnfilt[i]+'.***notapossiblefile', 'r') as testFile:
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
                        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as testFile:
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
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as testFile:
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
                with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'r') as TEC_file:
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
                        settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i]+'_oldV'+str(TEC_version)); #rename
                try:
                    with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'r') as TEC_fileUnfilt:
                        TEC_unfilt_version = TEC_fileUnfilt.attrs['version']; #read the version
                    #END WITH
                    if( TEC_unfilt_version >= version_unfilt ):
                        TEC_dataAvail[i] = 4; #note that data is already downloaded and converted, but not filtered
                    else:
                        print("\n==============~Warning~==============");
                        print("The file {} uses version {} while version {} is the most current version. The file will be renamed with _oldV{} and remade with the new version of the algorithm.\n".format(TEC_dataFilePathUnfilt[i], TEC_unfilt_version , version_unfilt, TEC_unfilt_version) );
                        os.rename(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], \
                                  settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i]+'_oldV'+str(TEC_unfilt_version)); #rename
                        #force an OSError to use the rest of the algorithm - uses ***'s which can't be in file names
                        with h5py.File(TEC_dataFilePathUnfilt[i]+'.***notapossiblefile', 'r') as testFile:
                            testFile.keys(); #tries to check out some stuff in the file
                        #END WITH
                    #END IF
                except (OSError , KeyError):
                    try: #gonna try to read the file - if we fail, it failed mid download
                        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as testFile:
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
            TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name 
            
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
                    TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name 
                    
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
                if( ( (TEC_dataAvail[i+1] != 1) | (TEC_dataAvail[i-1] != 1) ) & (TEC_dataAvail[i] == 1) ): #if the day before or the day after need to be filtered, and the current day is finished, change the signage
                    TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name 
                    
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
            TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name 
            
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
            
    #==============Download needed data and convert it==============
    current_webDownloadSpeed = 10; #MB/s, start off guess for download speed
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
            tic = time.time(); #for time testing
            urlretrieve(web_base_site+current_web_fileLinks,settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + current_webFileNames,reporthook=downloadProgress); #download the file in question to the data directory
            toc = time.time() - tic; #for time testing
            current_webDownloadSpeed = current_webFileSize/(1024**2)/toc; #MB/s, calc current web download speed
            print("\nTime to download: {:.2f} min, download speed: {:.2f} MB/s\n\n".format(toc/60,current_webDownloadSpeed)); #extra space at end
            
            if( TEC_dataAvail[i] == 2 ):
                TEC_dataAvail[i] = 3; #move on to the next stage ;) oh yeah it feels wrong but its so easy
            else:
                TEC_dataAvail[i] = 6; #move on to the next stage ;) oh yeah it feels wrong but its so easy - but remembers day is filtered and doesn't need to be redone
            #END IF
        #END IF
            
        #-----CONVERT THE DATA TO FASTER FORMAT-----
        if( (TEC_dataAvail[i] == 3) | (TEC_dataAvail[i] == 6) ): #3 means data was downloaded but needs to be converted (6 means data was downloaded but needs to be converted, but the current day is filtered and finished already)
            
            TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name
            
            dataVectMax = locIntUnfilt_size+locFloatUnfilt_size+locStringUnfilt_size+2; #total number of data vectors to import, order: int / float / string
            print("Converting {} to {}.\nAt {} data vectors and 120 sec per vector, expect {} minutes for conversion to finish.\n".format(TEC_fileNameOrig,TEC_dataFileNameUnfilt[i],dataVectMax,str(np.round(120*dataVectMax/60,2)).rstrip('0').rstrip('.')));
            
            tic = time.time(); #for time testing
            
            unfilt_dict = {
                'version':version_unfilt, #attributes here, code automagically saves them too
                };
                        
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as TEC_fileOrig:
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
                    sys.exit();
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
                    sys.exit();
                #END IF
                # TEC_fileData = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]]; #import the data shape
                dataVectCntr = 0; #counter to announce progress
                
                unfilt_dict['sat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'sat_id'].astype("int16"); #import satellite ID
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector sat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                unfilt_dict['year'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'year'].astype("int16"); #import year
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector year {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                tempMonth = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'month'].astype("int16"); #import month
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector month {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                tempDay = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'day'].astype("int16"); #import day
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector day {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                unfilt_dict['dayNum'] = subfun_date_to_dayNum( np.array( (unfilt_dict['year'],tempMonth,tempDay) ) , 2 ); #convert to day number, option 2 only returns day number
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector dayNum calc {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                unfilt_dict['hour'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'hour'].astype("int16"); #import hour
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector hour {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                unfilt_dict['min'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'min'].astype("int16"); #import min
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector min {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();
                unfilt_dict['sec'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'sec'].astype("int16"); #import sec
                dataVectCntr += 1; #increment
                sys.stdout.write("\rData vector sec {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2))); #report
                sys.stdout.flush();

                unfilt_dict['lat'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'gdlat'].astype(dataType_str); #arcdeg, import pierce point geodedic latitude
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector lat {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                unfilt_dict['long'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'glon'].astype(dataType_str); #arcdeg, import pierce point longitude
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector long {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                unfilt_dict['elev'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'elm'].astype(dataType_str); #deg, import pierce point elevation
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector elev {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                unfilt_dict['sTEC'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'los_tec'].astype(dataType_str); #TECU, import pierce point line-of-sight TEC
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector sTEC {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                unfilt_dict['sTECerr'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'dlos_tec'].astype(dataType_str); #TECU, import pierce point line-of-sight TEC error
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector sTECerr {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                unfilt_dict['vTEC'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'tec'].astype(dataType_str); #TECU, import adjusted vertical TEC (as recommended by AC)
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector vTEC {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
    #                TEC_fileData_float[:,locFloatUnfilt_TECerr] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'dlos_tec']; #TECU, import pierce point line-of-sight TEC error
    #                dataVectCntr = dataVectCntr + 1;
    #                sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
    #                sys.stdout.flush();
                #TEC_fileData_float[:,8] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'rec_bias']; #TECU, import receiver bias
                #TEC_fileData_float[:,9] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'drec_bias']; #TECU, import receiver bias error
                
                unfilt_dict['pierceAlt'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'pierce_alt'].astype(dataType_str); #km, import pierce point altitude
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector pierceAlt {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
                #pierce point altitude is the same number, no need to record it a bunch
                            
                unfilt_dict['site'] = TEC_fileOrig[TEC_fileKeys[TEC_fileKeys_dataIndex] + "/" + TEC_fileDatalocale[TEC_fileDatalocale_tableIndex]][:,'gps_site']; #import GPS site name            
                dataVectCntr = dataVectCntr + 1;
                sys.stdout.write("\rData vector site {}/{} & {} min\t\t\t\t".format(dataVectCntr,dataVectMax,np.round((time.time()-tic)/60,2)));
                sys.stdout.flush();
            #END WITH
            
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
            #---luckily years have no limits (that we know of??)---
            #---DELETE data not on the right day (could have incrememted past the right day due to stuff and things like a 24 hour/60 minute time on the current day---
            keyz = list(unfilt_dict.keys()); #get the keys
            TEC_fileData_logical_onDay = (unfilt_dict['dayNum'] == dateRange_dayNum_full[i,1]) & (unfilt_dict['year'] == dateRange_dayNum_full[i,0]); #find when the day reported is the day we want [and the year we want]
            for j in range(0,len(keyz)):
                if( np.isscalar(unfilt_dict[keyz[j]]) == False ):
                    unfilt_dict[keyz[j]] = unfilt_dict[keyz[j]][TEC_fileData_logical_onDay]; #keep only the good stuff
                #END IF
            #END FOR j
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
            print('\nDone reading data from original file.');
            
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
                        TEC_fileUnfilt.create_dataset(keyz[j], data=unfilt_dict[keyz[j]], chunks=h5pyChunkShape, compression='gzip'); #write that data , compression="gzip"
                    else:
                        #if size 1, add it as an attribute
                        TEC_fileUnfilt.attrs[keyz[j]] = unfilt_dict[keyz[j]]; #save the attribute
                    #END IF
                    sys.stdout.write("\rWriting {} to file & {} min | {} out of {}\t\t\t\t\t".format(keyz[j],np.round((time.time()-tic)/60,2),j+1,len(keyz)));
                    sys.stdout.flush();
                #END FOR j
                print('\nDone writing unfiltered file!');
            #END WITH
            if( FLG_deleteOrig == 1 ): #only delete if flag is on
                os.remove(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete the downloaded file now that we're done with it (save some space!)
            #END IF
            del unfilt_dict, tempMonth, tempDay; #clean the memory
            toc = time.time() - tic; #for time testing
            print("\nTime to convert: {} min\n".format(np.round(toc/60,2))); #extra space at end   
            
            if( TEC_dataAvail[i] == 3 ):
                TEC_dataAvail[i] = 4; #move on to the next stage
            else:
                TEC_dataAvail[i] = 7; #move on to the next stage - but remembers day is filtered and doesn't need to be redone
            #END IF
        #END IF
    #END FOR i
    
    #==============Filter Data after download & conversion==============
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
        
        #-----FILTER THE DATA-----
        if( TEC_dataAvail[i] == 4 ): #4 means data was downloaded and converted but needs to be filtered (this is where 7 jumps off)
            
            #prep the Sav-Gol filter for debiasing
            windowLen_savGol = np.int64(np.round(filter_savGolPeriod/TEC_dataRate)); #window length, 60 minutes in seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
            #from conversations with AC ^
            if( np.remainder(windowLen_savGol,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
                windowLen_savGol = windowLen_savGol + 1; #add 1 so it's odd
            #END IF
            
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
                unfilt_dict = {}; #prep this
                for j in range(-1,2): #run through the days to get the main day and the days around it
                    if( (TEC_dataAvail[i+j] == 4) | (TEC_dataAvail[i+j] == 7) | (TEC_dataAvail[i+j] == -3) ): #make sure data is there before reading it (4 or 7 mean it is there) and -3 also means it is there, just a padded day
                        with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i+j], 'r', rdcc_nbytes =500*1024*1024) as TEC_fileUnfilt:
                            #----- Read in the unfiltered data, keep only the good data -----
                            keyz = list(unfilt_dict.keys()); #get the saved keys
                            keyzNew = list(TEC_fileUnfilt.keys()); #get the saved keys
                            if( j == -1 ):
                                #j = -1 is guaranteed to be 1st, so can assume unfilt_dict is unfilled yet
                                unfilt_dict['time'] = TEC_fileUnfilt.get('time')[()]; #get that dataset out
                                unfilt_dict['elev'] = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
                                TEC_goodVals = unfilt_dict['time']/3600 > 18; #get UTC hours 18 and greater (6 hour buffer) to cut down on the data imported (would have to think more about the max theoretical observable GPS satellite time period - 6 is safe.)
                                TEC_goodVals = (unfilt_dict['elev'] >= minElevation) & TEC_goodVals; #get the locations where elevation is above the min elevation
                                unfilt_dict['time'] = [unfilt_dict['time'][TEC_goodVals]] #keep the good stuff
                                unfilt_dict['elev'] = [unfilt_dict['elev'][TEC_goodVals]]; #keep the good stuff
                                keyzNew.remove('time'); #remove from the list, manually got it for data corralling
                                keyzNew.remove('elev'); #remove from the list, manually got it for data corralling
                            elif( j == 1 ):
                                unfilt_dictTemp_time = TEC_fileUnfilt.get('time')[()]; #get that dataset out
                                unfilt_dictTemp_elev = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
                                TEC_goodVals = unfilt_dictTemp_time/3600 < 6; #get UTC hours 6 and less (6 hour buffer) to cut down on the data imported (would have to think more about the max theoretical observable GPS satellite time period - 6 is safe.)
                                TEC_goodVals = (unfilt_dictTemp_elev >= minElevation) & TEC_goodVals; #get the locations where elevation is above the min elevation
                                unfilt_dictTemp_time = unfilt_dictTemp_time[TEC_goodVals]; #keep the good stuff
                                unfilt_dictTemp_elev = unfilt_dictTemp_elev[TEC_goodVals]; #keep the good stuff
                                keyzNew.remove('time'); #remove from the list, manually got it for data corralling
                                keyzNew.remove('elev'); #remove from the list, manually got it for data corralling
                                if( strfind(keyz,'time',opt=1) > 0 ):
                                    unfilt_dict['time'].append(unfilt_dictTemp_time); #tack that dataset on
                                    del unfilt_dictTemp_time; #clean up the memory
                                else:
                                    #otherwise it's a new data type to add in
                                    unfilt_dict['time'] = [np.copy(unfilt_dictTemp_time)]; #get that dataset out
                                    del unfilt_dictTemp_time; #clean up the memory
                                #END IF
                                if( strfind(keyz,'elev',opt=1) > 0 ):
                                    unfilt_dict['elev'].append(unfilt_dictTemp_elev); #tack that dataset on
                                    del unfilt_dictTemp_elev; #clean up the memory
                                else:
                                    #otherwise it's a new data type to add in
                                    unfilt_dict['elev'] = [np.copy(unfilt_dictTemp_elev)]; #get that dataset out
                                    del unfilt_dictTemp_elev; #clean up the memory
                                #END IF
                            else:
                                unfilt_dictTemp_elev = TEC_fileUnfilt.get('elev')[()]; #get that dataset out
                                TEC_goodVals = unfilt_dictTemp_elev >= minElevation; #get the locations where elevation is above the min elevation
                                # TEC_goodVals = np.ones(TEC_fileData_float_temp[:,locFloatUnfilt_elev].shape,dtype=bool); #set all to good [elevation clipped later]
                                unfilt_dictTemp_elev = unfilt_dictTemp_elev[TEC_goodVals]; #keep the good stuff
                                keyzNew.remove('elev'); #remove from the list, manually got it for data corralling
                                if( strfind(keyz,'elev',opt=1) > 0 ):
                                    unfilt_dict['elev'].append(unfilt_dictTemp_elev); #tack that dataset on
                                    del unfilt_dictTemp_elev; #clean up the memory
                                else:
                                    #otherwise it's a new data type to add in
                                    unfilt_dict['elev'] = [np.copy(unfilt_dictTemp_elev)]; #get that dataset out
                                    del unfilt_dictTemp_elev; #clean up the memory
                                #END IF
                            #END IF
                            for k in range(0,len(keyzNew)):
                                if( strfind(keyz,keyzNew[k],opt=1) > 0 ):
                                    unfilt_dict[keyzNew[k]].append(TEC_fileUnfilt.get(keyzNew[k])[()][TEC_goodVals]); #tack that dataset on
                                else:
                                    #otherwise it's a new data type to add in
                                    unfilt_dict[keyzNew[k]] = [TEC_fileUnfilt.get(keyzNew[k])[()][TEC_goodVals]]; #get that dataset out, only keep the good stuff
                                #END IF
                            #END FOR k
                            #--- Read the attributes in ---
                            keyzNew = list(TEC_fileUnfilt.attrs.keys()); #get the attribute keys
                            keyzNew = np.delete(np.asarray(keyzNew),np.where(strfind(keyzNew,'version'))[0]).tolist(); #remove version from the list, only needed here
                            for k in range(0,len(keyzNew)):
                                if( strfind(keyz,keyzNew[k],opt=1) > 0 ):
                                    if( np.isclose(unfilt_dict[keyzNew[k]],TEC_fileUnfilt.attrs[keyzNew[k]]) == False ):
                                        #only worry is if the attribute isn't consistent
                                        print('-----Warning-----');
                                        print('Attribute '+keyzNew[k]+' isn\'t the same as the previously recorded value from another file of '+ \
                                            str(unfilt_dict[keyzNew[k]])+' and this file\'s value of '+str(TEC_fileUnfilt.attrs[keyzNew[k]])+ \
                                            '.\n NaN\'ing it and try to sort that out.');
                                        unfilt_dict[keyzNew[k]] = np.nan; #nan that attribute, figure it out later
                                    #END IF
                                else:
                                    unfilt_dict[keyzNew[k]] = TEC_fileUnfilt.attrs[keyzNew[k]]; #get that attribute out
                                #END IF
                            #END FOR k
                            del TEC_goodVals #clean up that memory
                        #END WITH
                                            
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
                        print("There is no data availiable on ~PADDED DAY~ {}/{}/{} in YR/M/D format. Padded days are required via passed flag.\nPrinting availiable data days for the relevant year - will exit on finishing checking all days:".format(dateRange_full[i+j,0],dateRange_full[i+j,1],dateRange_full[i+j,2]) );
                        print("{}".format(TEC_dataAvail)); #print for error - lets user know available days
                        print("{}".format(i+j)); #print for error - lets user know available days
                        print("{}".format(web_fileNamesDate)); #print for error - lets user know available days
                        TEC_dataAvail[i+j] = -1; #note the error
                        sys.exit(); #yolo out
                    #END IF              
                #END FOR j
                #--- Convert to numpy arrays from the dynamically added lists ---
                keyz = list(unfilt_dict.keys()); #get the current keys
                for j in range(0,len(keyz)):
                    if( np.isscalar(unfilt_dict[keyz[j]]) == False ):
                        #if not a scalar, apply the logical mask
                        unfilt_dict[keyz[j]] = np.hstack(unfilt_dict[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
                    #END IF
                #END FOR j
                
                toc = time.time() - tic; #for time testing
                print("Time to import all needed data into memory: {} min".format(round(toc/60,2))); #extra space at end  
                
                #==============THIS WILL GO IN A FUNCTION FOR MAX SPEED==============
                #TIME TO FILTER!
                tic = time.time(); #for time testing
                
                #-----prep dynamic list to hold data as it is accepted-----
                filt_dict = {
                    'sat':[],
                    'year':[],
                    'dayNum':[],
                    'hour':[],
                    'min':[],
                    'sec':[],
                    'time':[],
                    'lat':[],
                    'long':[],
                    'elev':[],
                    'dTEC':[],
                    'dTECerr':[],
                    'site':[],
                    }; #prep a dict to hold the data as it's made
                
    #            unfilt_site_unique = np.unique(unfilt_dict['site']); #get unique site names
                unfilt_site_unique, unfilt_site_unique_indexes = np.unique(unfilt_dict['site'] , return_inverse=True); #get unique site names and the indexes to speed up searching for where to look
                #inspired by https://stackoverflow.com/questions/33281957/faster-alternative-to-numpy-where
    #            unfilt_site_unique_sortedIndex = np.argsort(unfilt_dict['site'], kind='mergesort'); #returns the indexes IF the array unfilt_dict['site'] was sorted
                unfilt_site_unique_currentSiteArray = np.split(np.argsort(unfilt_dict['site'], kind='mergesort'), np.cumsum(np.bincount(unfilt_site_unique_indexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
                #faster than [for j in range... currentSite_loc = np.where(j == unfilt_site_unique_indexes)[0]; #get the data locations for 1 site per]
                        
                TEC_timeUnique, TEC_timeUniqueIndexes , TEC_timeUniqueCount = np.unique( unfilt_dict['time'] , return_inverse=True , return_counts=True);
                #Cut off time stamps with very little data in the range selected
    #            TEC_dataAvgNum = np.sum(TEC_timeUniqueCount)/TEC_timeUnique.size; #average data per time
    #            TEC_dataLim = np.round(TEC_dataLimPercent*TEC_dataAvgNum); #min number before eject time set
    #            TEC_timeUnique = TEC_timeUnique[ TEC_timeUniqueCount > TEC_dataLim ]; #remove the low data # stuff
                #disabled for now
                TEC_timeUnique_currentTimeArray = np.split(np.argsort(unfilt_dict['time'], kind='mergesort'), np.cumsum(np.bincount(TEC_timeUniqueIndexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
                                
                # cntr = 0; #for controlling the number of plots during debugging
                if( FLG_dataAggregation == 0 ): #with this off, just use what is directly there
                    #cruise through every site (from the unique of unfilt_dict['site'])
                    for j in range(0,unfilt_site_unique.size):
                        # print('unfilt_site_unique.size = '+str(unfilt_site_unique.size))
                        currentSite_loc = unfilt_site_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                        currentSat_unique = np.unique(unfilt_dict['sat'][currentSite_loc]); #get the unique sats in at that current site
                        #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])             
                        for k in range(0,currentSat_unique.size):
                            # print('currentSat_unique.size = '+str(currentSat_unique.size))
                            currentSat_loc = np.where( currentSat_unique[k] == unfilt_dict['sat'][currentSite_loc] )[0]; #get the data locations for 1 sat at that 1 site
                            #-----GET NEEDED LOCAL VARS-----
                            #---INTS---
                            currentSat = unfilt_dict['sat'][currentSite_loc][currentSat_loc]; #sat#, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentYear = unfilt_dict['year'][currentSite_loc][currentSat_loc]; #year, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentDayNum = unfilt_dict['dayNum'][currentSite_loc][currentSat_loc]; #day#, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentHour = unfilt_dict['hour'][currentSite_loc][currentSat_loc]; #hr, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentMin = unfilt_dict['min'][currentSite_loc][currentSat_loc]; #min, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentSec = unfilt_dict['sec'][currentSite_loc][currentSat_loc]; #sec, get the day data for 1 sat at that 1 site - uses abs illegal double indexing                        
                            #---FLOATS---
                            currentLat = unfilt_dict['lat'][currentSite_loc][currentSat_loc];
                            currentLong = unfilt_dict['long'][currentSite_loc][currentSat_loc];
                            currentElv = unfilt_dict['elev'][currentSite_loc][currentSat_loc];
                            #---STRINGS---
                            currentSite = unfilt_dict['site'][currentSite_loc][currentSat_loc];
                            #---UNFILT DATA---
                            currentvTEC = unfilt_dict['vTEC'][currentSite_loc][currentSat_loc]; #TECU, get the vTEC data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentsTECerror = unfilt_dict['sTECerr'][currentSite_loc][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentSecTotal = np.int64(unfilt_dict['dayNum'][currentSite_loc][currentSat_loc])*86400 + unfilt_dict['time'][currentSite_loc][currentSat_loc]; #sec, calculate days/hour/min/sec into total seconds but do it at 64 bit for extra good precision in the seconds [no overflow]
                            # currentsTEC = TEC_fileData_float[currentSite_loc,locFloatUnfilt_sTEC][currentSat_loc]; #TECU, get the sTEC data for 1 sat at that 1 site - uses abs illegal double indexing [only for debug]
                            
                            currentTimeSplits_loc = np.append( np.insert( np.where(np.diff(currentSecTotal) > minimumTimeGap)[0]+1 , 0, 0), currentSecTotal.shape ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                            
                            for l in range(0,len(currentTimeSplits_loc)-1):
                                # print('len(currentTimeSplits_loc)-1 = '+str(len(currentTimeSplits_loc)-1))
                                #-----Identify current TEC_dataRate & prep the filters-----
                                if( currentDayNum[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]].size > 1 ):
                                    current_TEC_dataRate = np.diff(currentSecTotal[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #sec, get the current data rate for the site
                                    uni, cnts = np.unique(current_TEC_dataRate,return_counts=True); #get unique rates and their counts
                                    current_TEC_dataRate = uni[cnts >= np.mean(cnts)]; #choose the counts that are larger than the mean (heavily influenced by outlier uniques with big counts (most likely the correct data rate))
                                    if( current_TEC_dataRate.size == 1 ):
                                        current_TEC_dataRate = np.int64(current_TEC_dataRate.item()); #get the current data rate as an integer
                                    else:
                                        current_TEC_dataRate = np.int64(np.min(current_TEC_dataRate)); #choose the minimum one, assuming 15's and 30's mean that the real time is 15 and there's just a lot of skips
                                    #END IF
                                else:
                                    current_TEC_dataRate = 0; #set current data rate to 0 so code doesn't run on it 
                                #END IF
                                
                                #make sure do work only on the day we want
                                #only get the data for times ranges thare are two times the filter period or longer
                                if( ( (currentSecTotal[currentTimeSplits_loc[l+1]-1]-currentSecTotal[currentTimeSplits_loc[l]]) > filter_savGolPeriod*2) & #make sure the data time period is 2*filterPeriod
                                (np.any(currentDayNum[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]] == dateRange_dayNum_full[i,1]) == True) &
                                (current_TEC_dataRate != 0) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                    #-----GET NEEDED LOCAL-LOCAL VARS-----
                                    #---INTS---
                                    currentSat_singleSatLine = currentSat[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentYear_singleSatLine = currentYear[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentDayNum_singleSatLine = currentDayNum[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentHour_singleSatLine = currentHour[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentMin_singleSatLine = currentMin[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentSec_singleSatLine = currentSec[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    #---FLOATS---
                                    currentLat_singleSatLine = currentLat[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentLong_singleSatLine = currentLong[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentElv_singleSatLine = currentElv[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    #---STRINGS---
                                    currentSite_singleSatLine = currentSite[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    #---UNFILT DATA---
                                    currentvTEC_singleSatLine = np.float64(currentvTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #this will increase readability
                                    currentsTECerror_singleSatLine = currentsTECerror[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    currentSecTotal_singleSatLine = currentSecTotal[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    # currentsTEC_singleSatLine = np.float64(currentsTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]); #this will increase readability [only for debug]
                                    
                                    if( currentDayNum_singleSatLine[0] != currentDayNum_singleSatLine[-1] ): #if this is true, we are looking at different days' data being used
                                    # if( np.any(currentDayNum_singleSatLine != dateRange_dayNum_full[i,1]) ): #if this is true, we are looking at different days' data being used [ponder this implementation]
                                        #the problem with that is that it seems the site can have different TEC fit values for different days
                                        #so the lines are continuous - but have a big jump where the different days are fit to different means(?) I'm not sure
                                        #but there's def a big jump we're gonna fix here
                                        currentDayNum_singleSatLine_firstDayLocations = np.where( currentDayNum_singleSatLine[0] == currentDayNum_singleSatLine )[0]; #get index where 
                                        currentDayNum_singleSatLine_secondDayLocations = np.where( currentDayNum_singleSatLine[-1] == currentDayNum_singleSatLine )[0]; #get index where 
                                        currentTEC_singleSatLine_diff = np.diff(currentvTEC_singleSatLine); #get the deltas between the TEC values
                                        currentTEC_singleSatLine_firstDayLast = currentvTEC_singleSatLine[currentDayNum_singleSatLine_firstDayLocations[-1]]; #get the first day's last TEC value
                                        currentTEC_singleSatLine_secondDay = currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations]; #get the second day's TEC values
                                        if( currentDayNum_singleSatLine_firstDayLocations.size > 1 ): #make sure currentDayNum_singleSatLine_firstDayLocations[-2] exists
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_firstDayLocations[-2]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        else: #if it doesn't, time to use currentDayNum_singleSatLine_secondDayLocations[0] instead of currentDayNum_singleSatLine_firstDayLocations[-2] (so projecting from reverse instead of forward)
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_secondDayLocations[0]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        #END IF
                                        currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations] = np.float64(currentTEC_singleSatLine_secondDay); #put em in
                                    #END IF
                                    
                                    #-----This code bit exists to catch if some weird thing happens where 2 satellite lines are recorded at the same time-----
                                    #it happens
                                    if( np.unique(currentSecTotal_singleSatLine).size == currentSecTotal_singleSatLine.size ):
                                        if( np.all(currentSecTotal_singleSatLine[:-1] <= currentSecTotal_singleSatLine[1:]) == False ):
                                            print('NOT SORTED AS MONOTONICALLY INCREASING TIME CORRECTLY, CHECK IT OUT');
                                            sys.crash(); #not a real call, actually crashes
                                        #END IF
                                        
                                        #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                        
                                        #get the time steps right is the first step!
                                        #~~~COMPACT THE TIME STAMPS IN THE DATA~~~
                                        #basically, it needs to meet the lowest common demoninator (30 sec per satellite reading)
                                        #10 sec data will be compacted by averaging 1-10, 11-20, and 21-30 (shown as 10, 20, 30 in the data) into just 30 at the
                                        #position of the 30 sec data pt
                                        #now it could be also an average of the position datas, but I'm going with the last position for now since I don't know how the position is chosen to start with so the easier way is the way to go!
                                        
                                        if( np.all(np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps) == False) ):
                                            print('CHECK THIS OUT!!! NO 0/30 AT ALL')
                                            sys.crash()
                                        #END IF
                                        
                                        #run this if the data rate is equal or below and if all of the sec time stamps don't match the allowed ones (sometimes there can be a 30 sec data rate but it's on 29/59 instead of 30/0)
                                        if( (current_TEC_dataRate <= TEC_dataRate) & np.any(np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps) == False) ): 
                                            
                                            #-----CONVERT NEARLY 0/30 TIMESTAMPS INTO 0/30 TIMESTAMPS-----
                                            for m in range(0,TEC_dataRate_allowedStamps.size):
                                                allowedStamps_range = np.array( (TEC_dataRate_allowedStamps[m]-TEC_dataRate*TEC_timeTolerance,TEC_dataRate_allowedStamps[m]+TEC_dataRate*TEC_timeTolerance) ); #get the current allowed time stamp tolerance range
                                                if( np.any(allowedStamps_range < 0) ):
                                                    #Deals with a -3 to 3 time frame around 0 or something
                                                    allowedStamps_range[allowedStamps_range < 0] = allowedStamps_range[allowedStamps_range < 0]+60; #keep in the 0-59 sec range
                                                    current_timeInTolerance = (currentSec_singleSatLine >= allowedStamps_range[0]) | (currentSec_singleSatLine <= allowedStamps_range[1]); #get the times within the time stamp tolerance range (so for 0 sec, 57 to 3 sec is taken as 0 sec)
                                                    #---CHECK SECONDS---
                                                    if( np.any(currentSec_singleSatLine >= allowedStamps_range[0]) ): #catches a 58 and rolls it over to 0 (+1 minute as well)
                                                        kj = currentSec_singleSatLine >= allowedStamps_range[0]; #get where the minute rolls over
                                                        currentMin_singleSatLine[kj] += 1; #min, the minutes are incremented as well if the seconds rolled over     
                                                        #---CHECK MINUTES---
                                                        kj = currentMin_singleSatLine >= 60; #find incorrect time keeping
                                                        if( np.sum(kj) > 0 ):
                                                            #increment the main hour variable where minutes are 60 or more
                                                            currentHour_singleSatLine[kj] += 1; #hour, increment time by 1
                                                            currentMin_singleSatLine[kj] = 0; #min, remove 60 time units from the time keeping
                                                            #---CHECK HOURS---
                                                            kj = currentHour_singleSatLine >= 24; #find incorrect time keeping
                                                            if( np.sum(kj) > 0 ):
                                                                #increment the main day variable where hours are 24 or more
                                                                currentDayNum_singleSatLine[kj] += 1; #day number, increment time by 1 [integer version]
                                                                currentHour_singleSatLine[kj] = 0; #hour, remove 24 time units from the time keeping
                                                                #---CHECK DAYS---
                                                                #deal with day limit is based on leap year or not
                                                                dayLim = np.ones(currentDayNum_singleSatLine.shape,dtype=np.int16)*365; #day number, get the day number limits as 365
                                                                #adjust leap years to 366
                                                                leapYears = (np.mod(currentYear_singleSatLine,4) == 0) & (np.mod(currentYear_singleSatLine,100) != 0) & (np.mod(currentYear_singleSatLine,400) == 0)
                                                                dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
                                                                kj = currentDayNum_singleSatLine >= dayLim; #find incorrect time keeping
                                                                if( np.sum(kj) > 0 ):
                                                                    #increment the main year variable where day number is equal to the day number limit or higher than it
                                                                    currentYear_singleSatLine += 1; #year, increment time by 1
                                                                    currentDayNum_singleSatLine[kj] = 1; #day num, remove the day number limit time units from the time keeping
                                                                #END IF
                                                            #END IF
                                                        #END IF
                                                    #END IF
                                                else:
                                                    current_timeInTolerance = (currentSec_singleSatLine >= allowedStamps_range[0]) & (currentSec_singleSatLine <= allowedStamps_range[1]); #get the times within the time stamp tolerance range (so for 30 sec, 27 to 33 sec is taken as 30 sec)
                                                #END IF
                                                currentSec_singleSatLine[current_timeInTolerance] = TEC_dataRate_allowedStamps[m]; #set the times within the tolerance range to the allowed time stamp for that range
                                            #END FOR m
                                            
                                            #-----AVERAGE SUB-30 SECOND TIMESTAMPS INTO 0/30 TIMESTAMPS-----
                                            #get where the second time stamps are the allowed values
                                            track_Matches = np.isin(currentSec_singleSatLine,TEC_dataRate_allowedStamps); #get the seconds that match the allowed times
                                            if( np.any(track_Matches == False) ): #only run if there are still time stmaps that aren't at allowed times                                                
                                                #the total time needs a recalc here (trying not to recalc too often)
                                                currentSecTotal_singleSatLine = np.int64(currentDayNum_singleSatLine)*86400 + np.int64(currentHour_singleSatLine)*3600 + np.int64(currentMin_singleSatLine)*60 + np.int64(currentSec_singleSatLine); #sec, calculate days/hour/min/sec into total seconds but do it at 64 bit for extra good precision in the seconds [no overflow]
                                                current_timeDiff = np.diff(currentSecTotal_singleSatLine); #get the diff between the times
                                                track_gaps = np.where( current_timeDiff > (current_TEC_dataRate*(1+TEC_timeTolerance)) )[0]; #catch time gaps that shouldn't be smashed together
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
                                                        #---CHECK SECONDS---
                                                        kj = currentSec_singleSatLine >= 60; #find incorrect time keeping
                                                        if( np.sum(kj) > 0 ):
                                                            #increment the main minute variable where seconds are 60 or more
                                                            currentMin_singleSatLine[kj] += 1; #min, the minutes are incremented as well if the seconds rolled over     
                                                            currentSec_singleSatLine[kj] = 0; #sec, remove 60 time units from the time keeping
                                                            #---CHECK MINUTES---
                                                            kj = currentMin_singleSatLine >= 60; #find incorrect time keeping
                                                            if( np.sum(kj) > 0 ):
                                                                #increment the main hour variable where minutes are 60 or more
                                                                currentHour_singleSatLine[kj] += 1; #hour, increment time by 1
                                                                currentMin_singleSatLine[kj] = 0; #min, remove 60 time units from the time keeping
                                                                #---CHECK HOURS---
                                                                kj = currentHour_singleSatLine >= 24; #find incorrect time keeping
                                                                if( np.sum(kj) > 0 ):
                                                                    #increment the main day variable where hours are 24 or more
                                                                    currentDayNum_singleSatLine[kj] += 1; #day number, increment time by 1 [integer version]
                                                                    currentHour_singleSatLine[kj] = 0; #hour, remove 24 time units from the time keeping
                                                                    #---CHECK DAYS---
                                                                    #deal with day limit is based on leap year or not
                                                                    dayLim = np.ones(currentDayNum_singleSatLine.shape,dtype=np.int16)*365; #day number, get the day number limits as 365
                                                                    #adjust leap years to 366
                                                                    leapYears = (np.mod(currentYear_singleSatLine,4) == 0) & (np.mod(currentYear_singleSatLine,100) != 0) & (np.mod(currentYear_singleSatLine,400) == 0)
                                                                    dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
                                                                    kj = currentDayNum_singleSatLine >= dayLim; #find incorrect time keeping
                                                                    if( np.sum(kj) > 0 ):
                                                                        #increment the main year variable where day number is equal to the day number limit or higher than it
                                                                        currentYear_singleSatLine += 1; #year, increment time by 1
                                                                        currentDayNum_singleSatLine[kj] = 1; #day num, remove the day number limit time units from the time keeping
                                                                    #END IF
                                                                #END IF
                                                            #END IF
                                                        #END IF
                                                    #END IF
                                                    #otherwise if gap starts on an allowed time I don't think it matters right now, the code will work fine. i think
                                                #END FOR m
                                                track_Matches_where = np.concatenate( (np.array( (0,) ),np.where(track_Matches == True)[0]) ); #get the indexes, put 0 at the front     
                                                #0 at the front is OK b/c if an allowed time, it'll average [0] and [0] together to get [0]
                                                #if not an allowed time, it'll average [0] and [1] together like it should - basically it gets it started good
                                                #---this is where the averaging actually happens, just took a lot of work to get ready for it---
                                                #now work through it and combine as needed
                                                for m in range(1,track_Matches_where.size-1):
                                                    currentvTEC_singleSatLine[track_Matches_where[m+1]] = np.nanmean(currentvTEC_singleSatLine[track_Matches_where[m]+1:track_Matches_where[m+1]+1]); #TECU, average them together
                                                    if( currentvTEC_singleSatLine[track_Matches_where[m]:track_Matches_where[m+1]+1].size == 0 ):
                                                        sys.crash();
                                                    #END IF
                                                #END FOR m
                                                
                                                #-----MAKE SURE THE LAST TIME STAMP VALUE ISN'T A NON-ALLOWED TIME-----
                                                if( track_Matches[-1] == False ):
                                                    #check if the last one is an off-time, and in that case
                                                    #push it to the next time slot <- decided to do this b/c predicting could end up rough and deleting is less data while the time scales are short compared to the 1 hour periodicity I'm looking at
                                                    # currentvTEC_singleSatLine[-1] = np.nanmean(currentvTEC_singleSatLine[track_Matches_where[-1]+1:currentvTEC_singleSatLine.size]); #TECU, average the trailing ones
                                                    # if( currentvTEC_singleSatLine[track_Matches_where[-1]+1:currentvTEC_singleSatLine.size].size == 0 ):
                                                        # sys.crash();
                                                    #END IF
                                                    TEC_dataRate_allowedStamps60 = np.copy(TEC_dataRate_allowedStamps); #copy
                                                    TEC_dataRate_allowedStamps60[TEC_dataRate_allowedStamps == 0] = 60; #force 0 to be 60 for easy math
                                                    TEC_dataRate_allowedStamps60 = np.sort(TEC_dataRate_allowedStamps60); #sort it 
                                                    track_shiftSec = np.where(TEC_dataRate_allowedStamps60 > currentSec_singleSatLine[-1])[0][0]; #get where the allowed indexes are greater than the current index, and get the first one always
                                                    currentSec_singleSatLine[-1] = TEC_dataRate_allowedStamps60[track_shiftSec]; #choose the right one
                                                    #---CHECK SECONDS---
                                                    kj = currentSec_singleSatLine >= 60; #find incorrect time keeping
                                                    if( np.sum(kj) > 0 ):
                                                        #increment the main minute variable where seconds are 60 or more
                                                        currentMin_singleSatLine[kj] += 1; #min, the minutes are incremented as well if the seconds rolled over     
                                                        currentSec_singleSatLine[kj] = 0; #sec, remove 60 time units from the time keeping
                                                        #---CHECK MINUTES---
                                                        kj = currentMin_singleSatLine >= 60; #find incorrect time keeping
                                                        if( np.sum(kj) > 0 ):
                                                            #increment the main hour variable where minutes are 60 or more
                                                            currentHour_singleSatLine[kj] += 1; #hour, increment time by 1
                                                            currentMin_singleSatLine[kj] = 0; #min, remove 60 time units from the time keeping
                                                            #---CHECK HOURS---
                                                            kj = currentHour_singleSatLine >= 24; #find incorrect time keeping
                                                            if( np.sum(kj) > 0 ):
                                                                #increment the main day variable where hours are 24 or more
                                                                currentDayNum_singleSatLine[kj] += 1; #day number, increment time by 1 [integer version]
                                                                currentHour_singleSatLine[kj] = 0; #hour, remove 24 time units from the time keeping
                                                                #---CHECK DAYS---
                                                                #deal with day limit is based on leap year or not
                                                                dayLim = np.ones(currentDayNum_singleSatLine.shape,dtype=np.int16)*365; #day number, get the day number limits as 365
                                                                #adjust leap years to 366
                                                                leapYears = (np.mod(currentYear_singleSatLine,4) == 0) & (np.mod(currentYear_singleSatLine,100) != 0) & (np.mod(currentYear_singleSatLine,400) == 0)
                                                                dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
                                                                kj = currentDayNum_singleSatLine >= dayLim; #find incorrect time keeping
                                                                if( np.sum(kj) > 0 ):
                                                                    #increment the main year variable where day number is equal to the day number limit or higher than it
                                                                    currentYear_singleSatLine += 1; #year, increment time by 1
                                                                    currentDayNum_singleSatLine[kj] = 1; #day num, remove the day number limit time units from the time keeping
                                                                #END IF
                                                            #END IF
                                                        #END IF
                                                    #END IF
                                                    track_Matches[-1] = True; #set the last value to true so it isn't NaN'd
                                                #END IF
                                            #END IF                                 
                                        elif( current_TEC_dataRate > TEC_dataRate ):
                                            #gotta interpolate a slower TEC data rate to a faster one
                                            #---GET INTERP'D TIMES---
                                            ideal_time = np.arange(currentSecTotal_singleSatLine[0],currentSecTotal_singleSatLine[-1]+TEC_dataRate,TEC_dataRate); #sec, get the ideal time
                                            real_times = np.isin(ideal_time,currentSecTotal_singleSatLine); #get a logical array of the real times
                                            #---GET TEC INTERPER---
                                            TEC_interper = interp1d(currentSecTotal_singleSatLine,currentvTEC_singleSatLine,kind='linear'); #make an interpolator
                                            #---INTERP TEC---
                                            ideal_TEC = TEC_interper(ideal_time); #make TEC for all the times
                                            ideal_TEC[real_times] = currentvTEC_singleSatLine; #make sure the real times are the real data     
                                            currentvTEC_singleSatLine = np.copy(ideal_TEC); #overwrite-filter later
                                            #FILTER LATER TO MAKE CODE JIVE
                                            #---INTERP TEC ERR(???)---
                                            #!!need to deal with mapping sTEC->vTEC on error, vTEC->delta-vTEC on error!!
                                            #note I'm not super sure how to handle TEC error, this is just to get the ball rolling (figure out good way later)
                                            TECerr_interper = interp1d(currentSecTotal_singleSatLine,currentsTECerror_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_TECerr = TECerr_interper(ideal_time); #make TEC for all the times
                                            ideal_TECerr[real_times] = currentsTECerror_singleSatLine; #make sure the real times are the real data
                                            currentsTECerror_singleSatLine = np.copy(ideal_TECerr); #for now keep it
                                            #---INTERP LAT---
                                            #!!THIS MAY BE A BAD FIT!!
                                            lat_interper = interp1d(currentSecTotal_singleSatLine,currentLat_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_lat = lat_interper(ideal_time); #make lat for all the times
                                            ideal_lat[real_times] = currentLat_singleSatLine; #make sure the real times are the real data
                                            currentLat_singleSatLine = np.copy(ideal_lat); #keep it now
                                            #---INTERP LONG---
                                            #!!THIS MAY BE A BAD FIT!!
                                            long_interper = interp1d(currentSecTotal_singleSatLine,currentLong_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_long = long_interper(ideal_time); #make long for all the times
                                            ideal_long[real_times] = currentLong_singleSatLine; #make sure the real times are the real data
                                            currentLong_singleSatLine = np.copy(ideal_long); #keep it now
                                            #---INTERP ELV---
                                            #!!THIS MAY BE A BAD FIT!!
                                            elv_interper = interp1d(currentSecTotal_singleSatLine,currentElv_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_elv = elv_interper(ideal_time); #make lat for all the times
                                            ideal_elv[real_times] = currentElv_singleSatLine; #make sure the real times are the real data
                                            currentElv_singleSatLine = np.copy(ideal_elv); #keep it now
                                            #---CALC DATA THAT IS NOT CONSTANT---
                                            currentDayNum_singleSatLine = np.int16(ideal_time//86400); #get day number
                                            midTime = ideal_time - np.int64(currentDayNum_singleSatLine)*86400; #get mid time w/o days
                                            currentHour_singleSatLine = np.int16(midTime//3600); #get hour
                                            midTime = midTime - np.int64(currentHour_singleSatLine)*3600; #get mid time w/o days
                                            currentMin_singleSatLine = np.int16(midTime//60); #get min
                                            midTime = midTime - np.int64(currentMin_singleSatLine)*60; #get mid time w/o days
                                            currentSec_singleSatLine = np.int16(midTime); #get sec
                                            #Calc day num float w/ all the recovered times above
                                            # currentDayNumF_singleSatLine = np.float64(currentDayNum_singleSatLine) + currentHour_singleSatLine/24 + currentMin_singleSatLine/1440 + currentSec_singleSatLine/86400; #days, calculate hour/min/sec into days and add to the current day but do it at 64 bit for extra good precision in the seconds
                                            #Calc sec int w/ all the recovered times above
                                            currentTotalSec_singleSatLine = np.int32(currentDayNum_singleSatLine)*86400 + np.int32(currentHour_singleSatLine)*3600 + np.int32(currentMin_singleSatLine)*60 + np.int32(currentSec_singleSatLine); #sec, calculate hour/min/sec into sec and add to the current day in sec
                                            #---REPEAT DATA THAT IS CONSTANT---
                                            if( np.all(currentSat_singleSatLine == currentSat_singleSatLine[0]) == False ):
                                                print('yo bad news'); #error check
                                                sys.crash()
                                            #END IF
                                            if( np.all(currentYear_singleSatLine == currentYear_singleSatLine[0]) == False ):
                                                print('yo bad news'); #error check
                                                sys.crash()
                                            #END IF
                                            if( np.all(currentSite_singleSatLine == currentSite_singleSatLine[0]) == False ):
                                                print('yo bad news'); #error check
                                                sys.crash()
                                            #END IF
                                            currentSat_singleSatLine = np.tile(currentSat_singleSatLine[0],(ideal_time.shape)); #copy sat ID (doesn't change)
                                            currentYear_singleSatLine  = np.tile(currentYear_singleSatLine[0],(ideal_time.shape)); #copy year (doesn't change FOR NOW)
                                            currentSite_singleSatLine = np.tile(currentSite_singleSatLine[0],(ideal_time.shape)); #copy site name (doesn't change)
                                            track_Matches = np.ones(currentvTEC_singleSatLine.shape,dtype=bool); #make them all true
                                        else:
                                            #if all times are on the required time cadence, make them all true
                                            track_Matches = np.ones(currentvTEC_singleSatLine.shape,dtype=bool); #make them all true
                                        #END IF                                        
                                        #-----REMOVE UNDEEDED DATA-----
                                        #---INTS---
                                        currentSat_singleSatLine = currentSat_singleSatLine[track_Matches]; #keep the good stuff
                                        currentYear_singleSatLine = currentYear_singleSatLine[track_Matches]; #keep the good stuff
                                        currentDayNum_singleSatLine = currentDayNum_singleSatLine[track_Matches]; #keep the good stuff
                                        currentHour_singleSatLine = currentHour_singleSatLine[track_Matches]; #keep the good stuff
                                        currentMin_singleSatLine = currentMin_singleSatLine[track_Matches]; #keep the good stuff
                                        currentSec_singleSatLine = currentSec_singleSatLine[track_Matches]; #keep the good stuff
                                        #---FLOATS---
                                        # currentDayNumF_singleSatLine = currentDayNumF_singleSatLine[track_Matches]; #keep the good stuff
                                        currentLat_singleSatLine = currentLat_singleSatLine[track_Matches]; #keep the good stuff
                                        currentLong_singleSatLine = currentLong_singleSatLine[track_Matches]; #keep the good stuff
                                        currentElv_singleSatLine = currentElv_singleSatLine[track_Matches]; #keep the good stuff
                                        #---STRINGS---
                                        currentSite_singleSatLine = currentSite_singleSatLine[track_Matches]; #keep the good stuff
                                        #---UNFILT DATA---
                                        currentvTEC_singleSatLine = currentvTEC_singleSatLine[track_Matches]; #keep the good stuff
                                        currentsTECerror_singleSatLine = currentsTECerror_singleSatLine[track_Matches]; #keep the good stuff
                                        # currentSecTotal_singleSatLine = currentSecTotal_singleSatLine[track_Matches]; #keep the good stuff
                                        
                                        #-----RECALC TOTAL TIME VARS TO UPDATE THEM TO CHANGED TIMES-----
                                        currentSecTotal_singleSatLine = np.int64(currentDayNum_singleSatLine)*86400 + np.int64(currentHour_singleSatLine)*3600 + np.int64(currentMin_singleSatLine)*60 + np.int64(currentSec_singleSatLine); #sec, calculate days/hour/min/sec into total seconds but do it at 64 bit for extra good precision in the seconds [no overflow]
                                        
                                        current_timeDiff = np.diff(currentSecTotal_singleSatLine); #get the diff between the times
                                        track_gaps = np.where( current_timeDiff > TEC_dataRate )[0]; #catch time gaps
                                        if( track_gaps.size == 0 ):
                                            currentPolyYvals = savgol_filter(currentvTEC_singleSatLine,windowLen_savGol, order_savGol ); #filter it up
                                            current_deltavTEC = currentvTEC_singleSatLine - currentPolyYvals;
                                            #recalc here so don't have to do it twice
                                            # currentDayNumF_singleSatLine = np.float64(currentDayNum_singleSatLine) + currentHour_singleSatLine/24 + currentMin_singleSatLine/1440 + currentSec_singleSatLine/86400; #days, calculate hour/min/sec into days and add to the current day but do it at 64 bit for extra good precision in the seconds
                                            #Calc sec int w/ all the recovered times above
                                            currentTotalSec_singleSatLine = np.int32(currentDayNum_singleSatLine)*86400 + np.int32(currentHour_singleSatLine)*3600 + np.int32(currentMin_singleSatLine)*60 + np.int32(currentSec_singleSatLine); #sec, calculate hour/min/sec into sec and add to the current day in sec
                                            
                                            #This is for debug viewing/methods paper
                                            # from importlib import reload
                                            # import GRITI_import_TEC_Madrigal_viewSatellitePass
                                            # reload(GRITI_import_TEC_Madrigal_viewSatellitePass)
                                            # from GRITI_import_TEC_Madrigal_viewSatellitePass import GRITI_import_TEC_Madrigal_viewSatellitePass
                                            # GRITI_import_TEC_Madrigal_viewSatellitePass(currentsTEC_singleSatLine, currentvTEC_singleSatLine, currentPolyYvals, current_deltavTEC, \
                                            #     currentElv_singleSatLine, currentSat_singleSatLine, currentSite_singleSatLine, currentYear_singleSatLine, currentDayNum_singleSatLine, \
                                            #     currentHour_singleSatLine, currentMin_singleSatLine, currentSec_singleSatLine, currentLat_singleSatLine, currentLong_singleSatLine, \
                                            #     0.5, minElevation, folder, savedVersion = True); #savedVersion also saves a nicer looking version of the plot
                                        else:
                                            #---GET INTERP'D TIMES---
                                            ideal_time = np.arange(currentSecTotal_singleSatLine[0],currentSecTotal_singleSatLine[-1]+TEC_dataRate,TEC_dataRate); #sec, get the ideal time
                                            real_times = np.isin(ideal_time,currentSecTotal_singleSatLine); #get a logical array of the real times
                                            #---GET TEC INTERPER---
                                            TEC_interper = interp1d(currentSecTotal_singleSatLine,currentvTEC_singleSatLine,kind='linear'); #make an interpolator
                                            #---INTERP TEC---
                                            ideal_TEC = TEC_interper(ideal_time); #make TEC for all the times
                                            ideal_TEC[real_times] = currentvTEC_singleSatLine; #make sure the real times are the real data                                          
                                            #---FILTER INTERP'D TEC---
                                            currentPolyYvals = savgol_filter(ideal_TEC,windowLen_savGol, order_savGol ); #filter it up
                                            current_deltavTEC = ideal_TEC - currentPolyYvals;
                                            #---INTERP TEC ERR(???)---
                                            #!!need to deal with mapping sTEC->vTEC on error, vTEC->delta-vTEC on error!!
                                            #note I'm not super sure how to handle TEC error, this is just to get the ball rolling (figure out good way later)
                                            TECerr_interper = interp1d(currentSecTotal_singleSatLine,currentsTECerror_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_TECerr = TECerr_interper(ideal_time); #make TEC for all the times
                                            ideal_TECerr[real_times] = currentsTECerror_singleSatLine; #make sure the real times are the real data
                                            currentsTECerror_singleSatLine = np.copy(ideal_TECerr); #for now keep it
                                            #---INTERP LAT---
                                            #!!THIS MAY BE A BAD FIT!!
                                            lat_interper = interp1d(currentSecTotal_singleSatLine,currentLat_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_lat = lat_interper(ideal_time); #make lat for all the times
                                            ideal_lat[real_times] = currentLat_singleSatLine; #make sure the real times are the real data
                                            currentLat_singleSatLine = np.copy(ideal_lat); #keep it now
                                            #---INTERP LONG---
                                            #!!THIS MAY BE A BAD FIT!!
                                            long_interper = interp1d(currentSecTotal_singleSatLine,currentLong_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_long = long_interper(ideal_time); #make long for all the times
                                            ideal_long[real_times] = currentLong_singleSatLine; #make sure the real times are the real data
                                            currentLong_singleSatLine = np.copy(ideal_long); #keep it now
                                            #---INTERP ELV---
                                            #!!THIS MAY BE A BAD FIT!!
                                            elv_interper = interp1d(currentSecTotal_singleSatLine,currentElv_singleSatLine,kind='linear'); #make an interpolator
                                            ideal_elv = elv_interper(ideal_time); #make lat for all the times
                                            ideal_elv[real_times] = currentElv_singleSatLine; #make sure the real times are the real data
                                            currentElv_singleSatLine = np.copy(ideal_elv); #keep it now
                                            #---CALC DATA THAT IS NOT CONSTANT---
                                            currentDayNum_singleSatLine = np.int16(ideal_time//86400); #get day number
                                            midTime = ideal_time - np.int64(currentDayNum_singleSatLine)*86400; #get mid time w/o days
                                            currentHour_singleSatLine = np.int16(midTime//3600); #get hour
                                            midTime = midTime - np.int64(currentHour_singleSatLine)*3600; #get mid time w/o days
                                            currentMin_singleSatLine = np.int16(midTime//60); #get min
                                            midTime = midTime - np.int64(currentMin_singleSatLine)*60; #get mid time w/o days
                                            currentSec_singleSatLine = np.int16(midTime); #get sec
                                            #Calc day num float w/ all the recovered times above
                                            # currentDayNumF_singleSatLine = np.float64(currentDayNum_singleSatLine) + currentHour_singleSatLine/24 + currentMin_singleSatLine/1440 + currentSec_singleSatLine/86400; #days, calculate hour/min/sec into days and add to the current day but do it at 64 bit for extra good precision in the seconds
                                            #Calc sec int w/ all the recovered times above
                                            currentTotalSec_singleSatLine = np.int32(currentDayNum_singleSatLine)*86400 + np.int32(currentHour_singleSatLine)*3600 + np.int32(currentMin_singleSatLine)*60 + np.int32(currentSec_singleSatLine); #sec, calculate hour/min/sec into sec and add to the current day in sec
                                            #---REPEAT DATA THAT IS CONSTANT---
                                            if( np.all(currentSat_singleSatLine == currentSat_singleSatLine[0]) == False ):
                                                print('yo bad news'); #error check
                                                sys.crash()
                                            #END IF
                                            if( np.all(currentYear_singleSatLine == currentYear_singleSatLine[0]) == False ):
                                                print('yo bad news'); #error check
                                                sys.crash()
                                            #END IF
                                            if( np.all(currentSite_singleSatLine == currentSite_singleSatLine[0]) == False ):
                                                print('yo bad news'); #error check
                                                sys.crash()
                                            #END IF
                                            currentSat_singleSatLine = np.tile(currentSat_singleSatLine[0],(ideal_time.shape)); #copy sat ID (doesn't change)
                                            currentYear_singleSatLine  = np.tile(currentYear_singleSatLine[0],(ideal_time.shape)); #copy year (doesn't change FOR NOW)
                                            currentSite_singleSatLine = np.tile(currentSite_singleSatLine[0],(ideal_time.shape)); #copy site name (doesn't change)
                                        #END IF
                                        
                                        #-----RECORD DATA IN THE BIG LISTS-----
                                        filt_dict['sat'].append(np.int16(currentSat_singleSatLine)); #tack that data on
                                        filt_dict['year'].append(np.int16(currentYear_singleSatLine)); #tack that data on
                                        filt_dict['dayNum'].append(np.int16(currentDayNum_singleSatLine)); #tack that data on
                                        filt_dict['hour'].append(np.int16(currentHour_singleSatLine)); #tack that data on
                                        filt_dict['min'].append(np.int16(currentMin_singleSatLine)); #tack that data on
                                        filt_dict['sec'].append(np.int16(currentSec_singleSatLine)); #tack that data on
                                        filt_dict['time'].append(np.int32(currentTotalSec_singleSatLine)); #tack that data on
                                        filt_dict['lat'].append(dataType_meth(currentLat_singleSatLine)); #tack that data on
                                        filt_dict['long'].append(dataType_meth(currentLong_singleSatLine)); #tack that data on
                                        filt_dict['elev'].append(dataType_meth(currentElv_singleSatLine)); #tack that data on
                                        filt_dict['dTEC'].append(dataType_meth(current_deltavTEC)); #tack that data on
                                        filt_dict['dTECerr'].append(dataType_meth(currentsTECerror_singleSatLine)); #tack that data on
                                        filt_dict['site'].append(currentSite_singleSatLine); #tack that data on
                                    #END IF
                                #END IF
                            #END FOR l
                        #END FOR k
                        toc = time.time() - tic; #for time testing
                        sys.stdout.write('\rSite {} of {} | {} min so far, {} sec per site, ETA {} min\t\t\t\t'.format(j+1,unfilt_site_unique.size,str(round(toc/60,2)).rstrip('0').rstrip('.'),str(round(toc/(j+1),4)).rstrip('0').rstrip('.'),str(round(toc/(j+1)*(unfilt_site_unique.size-j-1)/60,2)).rstrip('0').rstrip('.'))); #report
                        sys.stdout.flush();
                    #END FOR j
                else: #with dataAggregation on, time to try to include other data to make longer satellite tracks
                    #cruise through every site (from the unique of unfilt_dict['site'])
    #                    cntr = 0; #for testing
                    for j in range(0,unfilt_site_unique.size): 
                        
        #                currentSite_loc = np.where(j == unfilt_site_unique_indexes)[0]; #get the data locations for 1 site
                        currentSite_loc = unfilt_site_unique_currentSiteArray[j]; #pull it out of the pre-calc'd list of data locations for 1 site
                                        
                        currentSat_unique = np.unique(unfilt_dict['sat'][currentSite_loc]); #get the unique sats in at that current site
                        #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])
                        
    #                    currentSite_latRange = np.array( (np.min(TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat]), np.max(TEC_fileData_float[currentSite_loc,locFloatUnfilt_lat]) ) ); #get the lat range
    #                    currentSite_longRange = np.array( (np.min(TEC_fileData_float[currentSite_loc,locFloatUnfilt_long]), np.max(TEC_fileData_float[currentSite_loc,locFloatUnfilt_long]) ) ); #get the long range
    #                    
    #                    currentSite_nearbyLat = TEC_fileData_float[ (TEC_fileData_float[:,locFloatUnfilt_lat] >= currentSite_latRange[0]) & (TEC_fileData_float[:,locFloatUnfilt_lat] <= currentSite_latRange[0])  ,locFloatUnfilt_lat]
                        
                        for k in range(0,currentSat_unique.size):
                            currentSat_loc = np.where( currentSat_unique[k] == unfilt_dict['sat'][currentSite_loc] )[0]; #get the data locations for 1 sat at that 1 site
                            currentvTEC = unfilt_dict['vTEC'][currentSite_loc][currentSat_loc]; #TECU, get the vTEC data for 1 sat at that 1 site - uses abs illegal double indexing
        #                        currentsTECerror = TEC_fileData_float[currentSite_loc,locFloatUnfilt_sTECerr][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentDayNumF = np.float64(unfilt_dict['dayNum'][currentSite_loc][currentSat_loc]) + np.float64(unfilt_dict['time'][currentSite_loc][currentSat_loc])/86400; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
    #                        currentElv = TEC_fileData_float[currentSite_loc,locFloatUnfilt_elev][currentSat_loc];
                            
                            currentTimeSplits_loc = np.append( np.insert( np.where(np.diff(currentDayNumF)*86400 > minimumTimeGap)[0]+1 , 0, 0), currentDayNumF.shape ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                            
                            currentLat = unfilt_dict['lat'][currentSite_loc][currentSat_loc]; #degc, get the lat data for 1 sat at that 1 site - uses abs illegal double indexing
                            currentLong = unfilt_dict['long'][currentSite_loc][currentSat_loc]; #degc, get the TEC long for 1 sat at that 1 site - uses abs illegal double indexing
                            
                            for l in range(0,len(currentTimeSplits_loc)-1):
                                currentDayNumF_singleSatLine = currentDayNumF[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentDayNum_singleSatLine_min = (currentDayNumF_singleSatLine-np.round(np.mean(currentDayNumF_singleSatLine)))*144; #min/10, convert to minutes and /10 so it's more stable when fitting a curve
                                currentvTEC_singleSatLine = currentvTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
    
                                currentLat_singleSatLine = currentLat[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                currentLong_singleSatLine = currentLong[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                
    #                            currentElv_singleSatLine = currentElv[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]];
    
    #                            plt.figure();
    #                            plt.scatter( currentDayNumF_singleSatLine , currentLat_singleSatLine );
    #                            plt.title("Latitude Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Lat [degc]');
    #                            plt.show();
    #                            
    #                            plt.figure();
    #                            plt.scatter( currentDayNumF_singleSatLine , currentLong_singleSatLine );
    #                            plt.title("Longitude Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Long [degc]');
    #                            plt.show();
                                
    #                            plt.figure();
    #                            plt.scatter( currentLong_singleSatLine , currentLat_singleSatLine );
    #                            plt.title("Lat/Long Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Long [degc]');
    #                            plt.ylabel('Lat [degc]');
    #                            plt.show();
                                
    #                            plt.figure();
    #                            plt.scatter( currentDayNumF_singleSatLine , currentElv_singleSatLine );
    #                            plt.title("Elevation Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                            plt.xlabel('Time [days]');
    #                            plt.ylabel('Elevation [deg]');
    #                            plt.show();
                                                                                    
                                currentDataAgg_stepsFwdBack= np.int64(TEC_dataAgg_timeAdditionLimit/TEC_dataRate); #get the number of steps to take forward and back in time
                                currentDataAgg_dayNumSSLBefore = TEC_timeUnique[ np.where( np.min(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]-currentDataAgg_stepsFwdBack:np.where( np.min(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]  ]; #days, get the times before the observation period
                                currentDataAgg_dayNumSSLAfter = TEC_timeUnique[ (np.where( np.max(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]+1):np.where( np.max(currentDayNumF_singleSatLine) == TEC_timeUnique )[0][0]+currentDataAgg_stepsFwdBack+1  ]; #days, get the times before the observation period
                                
                                #only do this if any of the days are in the correct day
                                if( (np.any(np.int16(currentDayNumF[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                    
                                    currentDataAgg_dayNumSSLBefore_min = (currentDataAgg_dayNumSSLBefore - np.round(np.mean(currentDayNumF_singleSatLine)))*144; #min/10, get the times before the observation period
                                    currentDataAgg_dayNumSSLAfter_min = (currentDataAgg_dayNumSSLAfter - np.round(np.mean(currentDayNumF_singleSatLine)))*144; #min/10, get the times after the observation period
        
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
    #                                plt.title("Predicted Lat/Long Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                plt.xlabel('Long [degc]');
    #                                plt.ylabel('Lat [degc]');
    #                                plt.show();
                                    
                                    currentDataAgg_TECBefore = np.ones(currentDataAgg_stepsFwdBack , dtype=dataType_meth)*np.nan; #preallocate
                                    cntrSkip = 0; #start up a counter
                                    cntrTime = 0; #start up another counter
                                    #run backwards
                                    currentDataAgg_dayNumSSLBeforeRev = np.flip(currentDataAgg_dayNumSSLBefore); #flip it to make it easy
                                    currentDataAgg_dayNumSSLBeforeRevStartIndex = np.where(TEC_timeUnique == currentDataAgg_dayNumSSLBeforeRev[0])[0][0]; #where to start it off
                                    while( (cntrSkip < np.int64(minimumTimeGap/TEC_dataRate)) & (cntrTime < currentDataAgg_stepsFwdBack) ):
                                        
                                        currentDataAgg_currentTimeIndexes = TEC_timeUnique_currentTimeArray[currentDataAgg_dayNumSSLBeforeRevStartIndex+cntrTime]; #increment it up as we go to save processing
                                        
                                        currentDataAgg_closeIndexes = np.where( TEC_dataAgg_distToPts_degcSq >= ((unfilt_dict['lat'][currentDataAgg_currentTimeIndexes] - currentDataAgg_latGuessBefore[cntrTime])**2 + (unfilt_dict['long'][currentDataAgg_currentTimeIndexes] - currentDataAgg_longGuessBefore[cntrTime])**2) )[0];
                                        
                                        if( currentDataAgg_closeIndexes.size != 0 ): #record some data
                                            currentDataAgg_TECBefore[cntrTime] = np.mean(unfilt_dict['vTEC'][currentDataAgg_currentTimeIndexes[currentDataAgg_closeIndexes]]); #get the mean of the vTEC around the point we guessed
                                            cntrSkip = 0; #reset
                                        else: #else, note data skipped
                                            cntrSkip = cntrSkip + 1; #increment skip counter
                                        #END IF
                                        cntrTime = cntrTime + 1; #increment total time counter
                                    #END WHILE
                                    currentDataAgg_TECBefore = np.flip(currentDataAgg_TECBefore); #reverse back to correct time sequence
                                    
                                    currentDataAgg_TECAfter = np.ones(currentDataAgg_stepsFwdBack , dtype=dataType_meth)*np.nan; #preallocate
                                    cntrSkip = 0; #start up a counter
                                    cntrTime = 0; #start up another counter
                                    #run forwards
                                    currentDataAgg_dayNumSSLAfterStartIndex = np.where(TEC_timeUnique == currentDataAgg_dayNumSSLAfter[0])[0][0]; #where to start it off
                                    while( (cntrSkip < np.int64(minimumTimeGap/TEC_dataRate)) & (cntrTime < currentDataAgg_stepsFwdBack) ):
                                        
                                        currentDataAgg_currentTimeIndexes = TEC_timeUnique_currentTimeArray[currentDataAgg_dayNumSSLAfterStartIndex + cntrTime]; #increment it up as we go to save processing
                                        
                                        currentDataAgg_closeIndexes = np.where( TEC_dataAgg_distToPts_degcSq >= ((unfilt_dict['lat'][currentDataAgg_currentTimeIndexes] - currentDataAgg_latGuessAfter[cntrTime])**2 + (unfilt_dict['long'][currentDataAgg_currentTimeIndexes] - currentDataAgg_longGuessAfter[cntrTime])**2) )[0];
                                        
                                        if( currentDataAgg_closeIndexes.size != 0 ): #record some data
                                            currentDataAgg_TECAfter[cntrTime] = np.mean(unfilt_dict['vTEC'][currentDataAgg_currentTimeIndexes[currentDataAgg_closeIndexes]]); #get the mean of the vTEC around the point we guessed
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
                                    currentDataAgg_dayNumSSLExtended = np.hstack( (currentDataAgg_dayNumSSLBefore , currentDayNumF_singleSatLine , currentDataAgg_dayNumSSLAfter ) ); #days, stack em
                                                                    
                                else:
                                    currentDataAgg_numberBefore = 0; #set to 0
                                    currentDataAgg_numberAfter = 0; #set to 0
                                    currentDataAgg_number = 0; #set to 0
                                #END IF
                            
                                #make sure do work only on the day we want
                                #only get the data for times ranges thare are two times the filter period or longer
                                if( ( ((currentDayNumF[currentTimeSplits_loc[l+1]-1]-currentDayNumF[currentTimeSplits_loc[l]])*86400+currentDataAgg_number*TEC_dataRate) > filter_savGolPeriod*2) & #make sure the data time period is 2*filterPeriod
                                   ( (currentTimeSplits_loc[l+1]-currentTimeSplits_loc[l]+currentDataAgg_number) >= windowLen_savGol) & #and that the data number is greater than the winodw length (in the event data is skipped and it's not perfectly 30 sec)
                                   # ( (currentTimeSplits_loc[l+1]-currentTimeSplits_loc[l]+currentDataAgg_number) > (3*(filter_b.size-1)) ) & #and that the data number is greater than the filter padding length requirement (cannot be equal)
                                   (np.any(np.int16(currentDayNumF[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]) == dateRange_dayNum_full[i,1]) == True) ):  #check to make sure some of the time range includes the day we want (otherwise we don't need to process it!) & == True for explicitness
                                   
                                    if( np.int16(currentDayNumF_singleSatLine[0]) != np.int16(currentDayNumF_singleSatLine[-1]) ): #if this is true, we are looking at different days' data being used
                                        #the problem with that is that it seems the site can have different TEC fit values for different days
                                        #so the lines are continuous - but have a big jump where the different days are fit to different means(?) I'm not sure
                                        #but there's def a big jump we're gonna fix here
                                        currentDayNum_singleSatLine_firstDayLocations = np.where( np.int16(currentDayNumF_singleSatLine[0]) == np.int16(currentDayNumF_singleSatLine) )[0]; #get index where 
                                        currentDayNum_singleSatLine_secondDayLocations = np.where( np.int16(currentDayNumF_singleSatLine[-1]) == np.int16(currentDayNumF_singleSatLine) )[0]; #get index where 
                                        currentTEC_singleSatLine_diff = np.diff(currentvTEC_singleSatLine); #get the deltas between the TEC values
                                        currentTEC_singleSatLine_firstDayLast = currentvTEC_singleSatLine[currentDayNum_singleSatLine_firstDayLocations[-1]]; #get the first day's last TEC value
                                        currentTEC_singleSatLine_secondDay = currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations]; #get the second day's TEC values
                                        if( currentDayNum_singleSatLine_firstDayLocations.size > 1 ): #make sure currentDayNum_singleSatLine_firstDayLocations[-2] exists
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_firstDayLocations[-2]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        else: #if it doesn't, time to use currentDayNum_singleSatLine_secondDayLocations[0] instead of currentDayNum_singleSatLine_firstDayLocations[-2] (so projecting from reverse instead of forward)
                                            currentTEC_singleSatLine_secondDay = currentTEC_singleSatLine_secondDay*( currentTEC_singleSatLine_firstDayLast + currentTEC_singleSatLine_diff[currentDayNum_singleSatLine_secondDayLocations[0]] )/currentTEC_singleSatLine_secondDay[0]; #scale it so the first data pt of the second day is now the last data pt of the first day plus the delta between the first day's last and second to last pt
                                            #sorry about these confusing explanations but it def makes sense to me rightnow
                                        #END IF
                                        currentvTEC_singleSatLine[currentDayNum_singleSatLine_secondDayLocations] = currentTEC_singleSatLine_secondDay; #put em in
                                    #END IF
                                    
                                    #now make sure the TEC data before and after connects to the current TEC data (cause of means and day differences and whatever)                                                        
                                    if( currentDataAgg_numberBefore > 0 ): #make sure there's some before data
                                        currentDataAgg_TECBefore = currentDataAgg_TECBefore*( 2*currentvTEC_singleSatLine[0] - currentvTEC_singleSatLine[1] )/currentDataAgg_TECBefore[-1]; ##scale it so that the last Before data pt is the same as the first real data pt plus the delta between the first and second points
                                    #END IF
                                    if( currentDataAgg_numberAfter > 0 ): #make sure there's some after data
                                        currentDataAgg_TECAfter = currentDataAgg_TECAfter*( 2*currentvTEC_singleSatLine[-1] - currentvTEC_singleSatLine[-2] )/currentDataAgg_TECAfter[0]; ##scale it so that the first After data pt is the same as the last real data pt plus the delta between the last and second-to-last points
                                    #END IF
                                    currentDataAgg_TEC = np.hstack( ( currentDataAgg_TECBefore, currentvTEC_singleSatLine, currentDataAgg_TECAfter ) ); #stack em into one
                                    
                                    
                                    #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                    
                                    #windowLen_savGol = np.int64(np.round(len(currentvTEC_singleSatLine)/4)) - (np.mod(np.int64(np.round(len(currentvTEC_singleSatLine)/4)),2) - 1); #unused
                                    
                                    currentPolyYvals = savgol_filter(currentDataAgg_TEC,windowLen_savGol,order_savGol ); #filter it up
                                    current_deltavTEC = currentDataAgg_TEC - currentPolyYvals;
                                    
                                    #-----This code bit exists to catch where the fit greatly deviates from the start OR end, which seems to happen often-----
                                    #use median/median distance instead of mean/stdev for stability apparently
                                    #inspired by https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
                                    current_deltaTEC_diff = np.diff(current_deltavTEC); #get the difference between each value and its following neighbor value
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
                                            current_deltaTEC_endCut = current_deltavTEC.size; #set no values to be cut off
                                        #END IF  
                                        
                                        current_deltavTEC[0:current_deltaTEC_startCut] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
                                        current_deltavTEC[current_deltaTEC_endCut:current_deltavTEC.size] = np.nan; #set bad start/end fits to NaN for easy identification for removal later (can do it in one shot, etc.)
                                    #END IF
                                    
                                    #DEBUG PLOTS
                                    #if( np.any( np.abs(current_deltavTEC) > 1.5) ):
    #                                plt.figure();
    #                                plt.scatter( currentDayNumF_singleSatLine , currentvTEC_singleSatLine , 20 , "r" );
    #                                plt.scatter( currentDayNumF_singleSatLine , currentPolyYvals , 20 );
    #                                plt.title("Fit Sav-Gol: Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                plt.xlabel('Time [days]');
    #                                plt.ylabel('los TEC [TECU]');
    #                                plt.show();
        #                            
        #                            plt.figure();
        #                            plt.scatter( currentDayNumF_singleSatLine , currentElv_singleSatLine );
        #                            plt.title("Elevation Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
        #                            plt.xlabel('Time [days]');
        #                            plt.ylabel('Elevation [deg]');
        #                            plt.show();
        #                           
    #                                if(  currentDataAgg_number > 0 ):
    #                                    plt.figure();
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentDataAgg_TEC , 20 , "r" );
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentPolyYvals , 20 );
    #                                    plt.title("Fit Sav-Gol: Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('los TEC [TECU]');
    #                                    plt.show();
    #                                    
    #                                    plt.figure();
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , currentDataAgg_TEC , 20 , "r" );
    #                                    plt.scatter( currentDayNumF_singleSatLine , currentvTEC_singleSatLine , 20 );
    #                                    plt.title("OG TEC B, Extra TEC R: Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('los TEC [TECU]');
    #                                    plt.show();
    #                                    
    #                                    plt.figure();
    #                                    plt.scatter( currentDataAgg_dayNumSSLExtended , current_deltavTEC );
    #                                    plt.title("Projected Delta TEC Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
    #                                    plt.xlabel('Time [days]');
    #                                    plt.ylabel('delta los TEC [TECU]');
    #                                    plt.show();
    #                                #END IF
                                    
                                    
        
                                    
                                                                
                                    # ===============Highpass filtering================
                                    # current_deltavTEC = signal.filtfilt(filter_b,filter_a,current_deltavTEC,padtype='odd',padlen=3*(filter_b.size-1)); #Appies the filter
                                    #extra bits from https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
                                    
                                    # ===============Cut to the correct size================
                                    if( (currentDataAgg_numberBefore >= 0) & (currentDataAgg_numberAfter > 0) ):
                                        current_deltavTEC = current_deltavTEC[currentDataAgg_numberBefore:-currentDataAgg_numberAfter]; #remove any of the added TEC values - they were just used to increase the possible "visible" time of a satellite
                                    elif( (currentDataAgg_numberBefore > 0) & (currentDataAgg_numberAfter == 0) ): #this needs some lil help
                                        current_deltavTEC = current_deltavTEC[currentDataAgg_numberBefore:currentDataAgg_dayNumSSLExtended.size+1]; #remove any of the added TEC values - they were just used to increase the possible "visible" time of a satellite
                                    #END IF
                                    
    #                                if(  currentDataAgg_number > 0 ):
    #                                    plt.figure();
    #                                    plt.scatter( currentDayNumF_singleSatLine , current_deltavTEC );
    #                                    plt.title("Trimmed Delta TEC Data for Site "+unfilt_site_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
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
                                    # current_deltavTEC = current_deltavTEC - np.nanmean(current_deltavTEC); #subtract the mean to make it a mean of 0
                                    if( np.sum(np.isnan(current_deltavTEC)) == current_deltavTEC.size):
                                        sys.exit();
                                    #END IF
                                                            
                                else:
        #                            currentvTEC_singleSatLine = currentvTEC[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]; #this will increase readability
                                    current_deltavTEC = np.ones(currentTimeSplits_loc[l+1]-currentTimeSplits_loc[l],dtype=dataType_meth)*np.nan; #no filter, just filler (NaN to easy identify)
                                #END IF
                                
                                #no matter what, write it in!
        #                        TEC_fileData_float[currentSite_loc,locFloatUnfilt_dTEC][currentSat_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]]  = current_deltavTEC; #save that deltaTEC into wherever it should go - assigning with double indexing DOES NOT WORK
#                                TEC_fileData_float[currentSite_loc[currentSat_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]],locFloatUnfilt_dTEC] = current_deltavTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                                # TEC_fileData_dTEC[currentSite_loc[currentSat_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]]]] = current_deltavTEC; #save that deltaTEC into wherever it should go - getting around double indexing
                                #WRONG BELOW
        #                        TEC_fileData_float[currentSite_loc[currentTimeSplits_loc[l]:currentTimeSplits_loc[l+1]],locFloatUnfilt_dTEC] = current_deltavTEC; #save that deltaTEC into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                            #END FOR l
                        #END FOR k
                        toc = time.time() - tic; #for time testing
                        sys.stdout.write('\rSite {} of {} | {} min so far, {} sec per site, ETA {} min\t\t\t\t'.format(j+1,unfilt_site_unique.size,str(round(toc/60,2)).rstrip('0').rstrip('.'),str(round(toc/(j+1),4)).rstrip('0').rstrip('.'),str(round(toc/(j+1)*(unfilt_site_unique.size-j-1)/60,2)).rstrip('0').rstrip('.'))); #report
                        sys.stdout.flush();                    
                        #END FOR j
                #END IF dataAggregation
                
                #----DELETE UNFILTERED DATA----
                #Save memory, done with it
                del unfilt_dict;
                
                #-----COMBINE LISTS OF GOOD DATA INTO MECHA ARRAYS-----
                keyz = list(filt_dict.keys()); #get the current keys
                for j in range(0,len(keyz)):
                    if( np.isscalar(filt_dict[keyz[j]]) == False ):
                        #if not a scalar, apply the logical mask
                        filt_dict[keyz[j]] = np.hstack(filt_dict[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
                    #END IF
                #END FOR j
                
                toc = time.time() - tic; #for time testing
                print("\nTime to filter all data: {} min".format(np.round(toc/60,2))); #extra space at end  
                
                tic = time.time(); #start a new timer
                #-----Get where the good data is-----
                #---NAN CONTROL---
                #Time to remove NaNs that were accumulated for TEC_fileData_float[:,0] (they mean data was bad)
    #            TEC_fileData_logical_TECnotnans = np.where(np.logical_not(np.isnan(TEC_fileData_float[:,locFloatUnfilt_dTEC])) == True)[0]; #find the not NaNs (index seems to be a bit faster maybe)
                TEC_logical_TECnotnans = np.logical_not(np.isnan(filt_dict['dTEC'])); #find the not NaNs    
                #---NON-DAY CONTROL---
                TEC_logical_onDay = (filt_dict['dayNum'] == dateRange_dayNum_full[i,1]) & (filt_dict['year'] == dateRange_dayNum_full[i,0]); #find when the day reported is the day we want and the year
                #---OUTLIER CONTROL---
                #from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list b/c median and median dist are more stable
                distFromMedian = np.abs(filt_dict['dTEC'] - np.nanmedian(filt_dict['dTEC'])); #get distance from median
                medianDistFromMedian = np.nanmedian(distFromMedian); #get median distance from median
                distFromMedian = distFromMedian/medianDistFromMedian if medianDistFromMedian else 0.0; #get normalized distance from median [reuse var to save memory]
                refDistFromMedian = np.abs(TEC_maxAmpAllowed - np.nanmedian(filt_dict['dTEC'])); #distance from median for 6 TECU [for dTEC that's a lot]
                refDistFromMedian = refDistFromMedian/medianDistFromMedian; #norm it so it can be used to compare
                TEC_logical_withinBounds = distFromMedian < refDistFromMedian; #get where the normalized distance from median is under the reference ceiling normalized value
                del distFromMedian
                #---MIN ELEVATION CONTROL---
                TEC_logical_withinElevation = filt_dict['elev'] >= minElevation; #get the locations where elevation is above the min elevation
                #---COMBINE ALL CONTROLS---
                TEC_logical_combined = np.where(TEC_logical_TECnotnans & TEC_logical_onDay & TEC_logical_withinBounds & TEC_logical_withinElevation)[0]; #combine them, get the index (seems slightly faster maybe)
                del TEC_logical_TECnotnans, TEC_logical_onDay, TEC_logical_withinBounds, TEC_logical_withinElevation; #save RAM
                #---CLEAR OUT STUFF---
                keyz = list(filt_dict.keys()); #get the current keys
                for j in range(0,len(keyz)):
                    if( np.isscalar(filt_dict[keyz[j]]) == False ):
                        #if not a scalar, apply the logical mask
                        filt_dict[keyz[j]] = filt_dict[keyz[j]][TEC_logical_combined]; #keep only the good stuff
                    #END IF
                #END FOR j
                del TEC_logical_combined; #save RAM
                                
                #-----TIME TO SAVE FILTERED DATA!-----
                #with h5py.File(TEC_dataFilePath[i], 'w') as TEC_file:
                h5pyChunkShape = filt_dict['dTEC'].shape; #get the shape of one of the vectors and use it as a chunk (read only whole chunks)
                #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
                with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'w', rdcc_nbytes =500*1024*1024) as TEC_file:
                    for j in range(0,len(keyz)):
                        if( np.isscalar(filt_dict[keyz[j]]) == False ):
                            TEC_file.create_dataset(keyz[j], data=filt_dict[keyz[j]], chunks=h5pyChunkShape, compression='gzip'); #write that data
                        else:
                            #if size 1, add it as an attribute
                            TEC_file.attrs[keyz[j]] = filt_dict[keyz[j]]; #save the attribute
                        #END IF
                        del filt_dict[keyz[j]]; #clean up that memory
                        sys.stdout.write("\rWriting {} to file & {} min | {} out of {}\t\t\t\t\t".format(keyz[j],np.round((time.time()-tic)/60,2),j+1,len(keyz)));
                        sys.stdout.flush();
                    #END FOR j
                    #add on non-data-related attributes
                    if( len(TEC_fileData_paddingWarning) > 0 ):
                        TEC_file.attrs['paddedDayMissing'] = " and ".join(TEC_fileData_paddingWarning) + " missing"; #record the attribute
                    else:
                        TEC_file.attrs['paddedDayMissing'] = "None, before and after padded"; #record the attribute
                    #END IF
                    TEC_file.attrs['TECmaxAmplitudeAllowed'] = TEC_maxAmpAllowed; #record the attribute
                    TEC_file.attrs['forcedTECdataRateSec'] = TEC_dataRate; #record the attribute
                    # TEC_file.attrs['medianRejectRatio'] = deltaTEC_compareValue; #record the attribute
                    TEC_file.attrs['savgolFiltPeriodSec'] = filter_savGolPeriod; #record the attribute
                    # TEC_file.attrs['highpassFiltPeriodHr'] = filter_cutoffPeriod; #record the attribute
                    # TEC_file.attrs['highpassFiltOrd'] = filter_n+1; #record the attribute
                    # TEC_file.attrs['highpassFiltWindow'] = "FIR Hanning"; #record the attribute
                    # TEC_file.attrs['highpassFiltType'] = "signal.filtfilt, padtype='odd',padlen=3*(b.size-1)"; #record the attribute
                    TEC_file.attrs['version'] = version_filt; #record the filtered algorithm version
                    print('\nDone writing filtered file!');
                #END WITH
                
                del filt_dict; #clean memory
                
                TEC_dataAvail[i] = 7; #data filtered and done! setting to 7 is good, setting to 1 is bad because 7 implies 4 while 1 doesn't guarantee that at alllllll
                #this is an artform of index juggling buddi
                
                toc = time.time() - tic; #for time testing
                print('\nTime to save filtered data: {} min'.format(np.round(toc/60,2))); #extra space at end  
                
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
    TEC_data = {}; #prep the dict
    for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
        #-----READ THE FILTERED DATA-----
        if( (TEC_dataAvail[i] == 1) | (TEC_dataAvail[i] == 7) ): #1 means data is there, filtered and all - 7 also means that
            keyz = list(TEC_data.keys()); #get the current keys in TEC_data
            with h5py.File(settings_paths['data'] + '\\' + paths_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileName[i], 'r') as TEC_file:
                #----- Read the data types in -----
                #--- First prep the mask and import needed vars for it ---
                keyzNew = list(TEC_file.keys()); #get the saved keys
                TEC_dataTemp_lat = TEC_file.get('lat')[()]; #get that dataset out
                TEC_dataTemp_long = TEC_file.get('long')[()]; #get that dataset out
                TEC_goodVals = np.where( (TEC_dataTemp_long <= np.max(plotLongRange)) & (TEC_dataTemp_long >= np.min(plotLongRange)) & \
                    (TEC_dataTemp_lat <= np.max(plotLatRange)) & (TEC_dataTemp_lat >= np.min(plotLatRange)) )[0]; #delete out of lat/long range stuff
                TEC_dataTemp_lat = TEC_dataTemp_lat[TEC_goodVals]; #keep the good stuff
                TEC_dataTemp_long = TEC_dataTemp_long[TEC_goodVals]; #keep the good stuff
                keyzNew.remove('lat'); #remove from the list, manually got it for data corralling
                keyzNew.remove('long'); #remove from the list, manually got it for data corralling
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