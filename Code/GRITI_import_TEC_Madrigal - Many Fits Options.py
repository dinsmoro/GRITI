#Function to import TEC data from Madrigal3
#RD on 8/23/2018
#

import numpy as np
import warnings
warnings.simplefilter("ignore", np.RankWarning); #ignore warning in np.polynomial.polynomial.polyfit
#import scipy as sp
from scipy import interpolate
from scipy.optimize import least_squares
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy.optimize import fsolve
from scipy.signal import argrelextrema
from scipy.signal import argrelmin 
from scipy.signal import argrelmax
from scipy.signal import savgol_filter
import os
import time
from urllib.request import urlopen, urlretrieve
import html2text
from subfun_findall import strstr
from subfun_downloadProgress import downloadProgress
import h5py
from subfun_date_to_dayNum import subfun_date_to_dayNum
from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
from subfun_addADay import subfun_addADay
from subfun_dayNum_to_date import subfun_dayNum_to_date
#from numba import jit, prange
#-----Testing variables-----
#http://cedar.openmadrigal.org/ftp/fullname/Ross+Dinsmore/email/rld5204@psu.edu/affiliation/PSU/kinst/8000/year/2015/kindat/3505/format/hdf5/
import sys
import matplotlib.pyplot as plt

#Date range goes Month-Day-Year
#dateRange = np.array([[2013,5,7],[2013,5,7]],dtype="int16"); #dates are in int16 because they can be
dateRange = np.array([[2014,7,31],[2014,7,31]],dtype="int16"); #dates are in int16 because they can be
#dateRange = np.array([[2016,11,28],[2016,11,30]],dtype="int16"); #dates are in int16 because they can be
#dateRange = np.array([[2014,12,31],[2015,1,1]],dtype="int16"); #for debug, check year success
#dates better go earlier -> later
#print("{}".format(dateRange))
folder = [os.getcwd()]; #current working directory, leave it as this call usually
folder.append('E:\Big Data'); #place to save data files to
#folder var structure: 0 = running folder, 1 = data folder
FLG_reqPaddedDays = 0; #0 no padded day required, 1 required padded days (better science since filtering on day edges will have data cutoff)
FLG_deleteUnfilt = 0; #0 don't delete unfiltered TEC data file, 1 delete filtered TEC data file
minElevation = 30; #deg, min elevation angle accepted
minimumTimeGap = 5; #min, how long to accept a gap in a satellite observation before calling it a new, unrelated observation
filterPeriod = .005*60; #min, filter period

dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
print("\nDate range requested (yr/day num format): {}/{} to {}/{}.".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
dateRange_dayNum_full_orig = dateRange_dayNum_full; #record orginal date ranges
dateRange_full_orig = subfun_dayNum_to_date(dateRange_dayNum_full); #record orginal date ranges
dateRange_dayNum_full = subfun_addADay(dateRange_dayNum_full); #call fun to pad a day onto the ends of the range
dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #convert the padded range to date range
print("Date range used due to padding requirement for filtering (yr/day num format): {}/{} to {}/{}.\n".format(dateRange_dayNum_full[0,0],dateRange_dayNum_full[0,1],dateRange_dayNum_full[-1,0],dateRange_dayNum_full[-1,1]) )
#The padding is required for filtering since the two days 

#helper for mininimize and fmin functions to calculate the error
#inspired by answer https://stackoverflow.com/questions/52563440/python-scipy-optimize-minimize-fails-with-valueerror-setting-an-array-element
def funErrorHlpr(coef, x, y):
    errorVect = currentFun(coef,x) - y; #calc error (will be a vector when called by minimize/fmin)
    errorSqSum = (errorVect**2).sum(); #calc sum of error squared as a single value that minimize/fmin can work with
    return errorSqSum;
#END DEF
    
def funErrorHlprFirst(x, coef):
    errorVect = currentFunRootzFinderFirst(coef,x); #calc error (will be a vector when called by minimize/fmin)
    errorSqSum = (errorVect**2).sum(); #calc sum of error squared as a single value that minimize/fmin can work with
    return errorSqSum;
#END DEF
    
def funErrorHlprSecond(x, coef):
    errorVect = currentFunRootzFinderSecond(coef,x); #calc error (will be a vector when called by minimize/fmin)
    errorSqSum = (errorVect**2).sum(); #calc sum of error squared as a single value that minimize/fmin can work with
    return errorSqSum;
#END DEF

#def GRITI_import_TEC(dateRange_dayNum_full,folder,minElevation,FLG_reqPaddedDays = 0, FLG_deleteUnfilt = 0):
#==============Constants Needed==============
folder_TEC = 'TEC'; #name for the TEC folder
folder_fileEnding = '.h5'; #file extension for hdf5 files
web_base_name = 'Ross+Dinsmore'; #put your name here
web_base_email = 'rld5204@psu.edu'; #put your email here
web_base_affil = 'PSU'; #put your affiliation here (school etc)
web_base_site = 'http://cedar.openmadrigal.org'; #website to be used
web_base = '/ftp/fullname/' + web_base_name + '/email/' + web_base_email + '/affiliation/' + web_base_affil + '/kinst/8000/year/'; #year goes after this
web_baseAfter = '/kindat/3505/format/hdf5/'; #this goes after year
web_baseFileDL = 'fullFilename/%252Fopt%252Fmadrigal3%252Fexperiments3%252F2016%252Fgps%252F'; #goes after above, for DL of file - date string in form DDmonYY%252Flos_YYYYMMDD
web_baseFileDLAfter = '.001.h5/'; #after bit above, final piece of the puzzle

#==============Check if data already there==============
if( os.path.isdir(folder[1] + '\\' + folder_TEC) == 0 ): #check if TEC folder exists
    #if not, make it
    os.makedirs(folder[1] + '\\' + folder_TEC);
    print("NOTA BENE: Importing TEC Func - Created TEC directory: {}\n".format(folder[1] + '\\' + folder_TEC) );
#END IF
    
dateRange_uniqueYears = np.unique(dateRange_dayNum_full[:,0]); #get the unique years involved
TEC_dataAmnt = len(dateRange_dayNum_full[:,0]); #get number of days needed to investigate
TEC_dataAvail = np.zeros( (TEC_dataAmnt,) , dtype='int8');  #preallocate logical array
#-1 = no data available on required days, will quit
#0 = no data available on padded day and flag says its OK to skip. Will impact data availability on edge of days, but middle of day is OK
#1 = note data is there, filtered and all
#2 = data is there and needs to be downloaded - OR wasn't downloaded fully and needs to be redownloaded, end goal is the same
#3 = note data is already downloaded but needs to be converted to the standardized format & naming scheme
#4 = note that data is already downloaded and converted, but not filtered
TEC_dataPath = ["{}\{}\{}".format(a_, b_, c_) for a_, b_, c_ in zip([folder[1]]*TEC_dataAmnt, [folder_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)) ) ]; #get the base path where data will be in
TEC_dataFileName = ["{}_{}_{}{}".format(a_, b_, c_, d_) for a_, b_, c_, d_ in zip([folder_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)),np.ndarray.tolist(dateRange_dayNum_full[:,1].astype(str)),[folder_fileEnding]*TEC_dataAmnt ) ]; #get the expected filenames
TEC_dataFileNameUnfilt = ["{}_{}_{}_unfilt{}".format(a_, b_, c_, d_) for a_, b_, c_, d_ in zip([folder_TEC]*TEC_dataAmnt,np.ndarray.tolist(dateRange_dayNum_full[:,0].astype(str)),np.ndarray.tolist(dateRange_dayNum_full[:,1].astype(str)),[folder_fileEnding]*TEC_dataAmnt ) ]; #get the expected filenames for unfiltered data (if stopped mid filtering)
TEC_dataFilePath = ["{}\{}".format(a_, b_) for a_, b_ in zip(TEC_dataPath,TEC_dataFileName ) ]; #get the full path right to the expected files
TEC_dataFilePathUnfilt = ["{}\{}".format(a_, b_) for a_, b_ in zip(TEC_dataPath,TEC_dataFileNameUnfilt ) ]; #get the full path right to the expected files that are unfiltered

for i in range(0,len(dateRange_uniqueYears)): #loop to check if data folder for the year exists
    if( os.path.isdir(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) ) == 0 ): #check if date folders exist
        #doesn't exist, gotta make it
        os.makedirs(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) );
        print("NOTA BENE: Importing TEC Func - Created TEC subdirectory for year {}: {}\n".format(dateRange_uniqueYears[i],folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_uniqueYears[i]) ));
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
        
        if( (os.path.isfile(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) == 0) and (os.path.isfile(TEC_dataFilePathUnfilt[i]) == 0) ):
            #if the unconverted file is not there, download it
            #first - need to check if TEC data is availiable for the days requested
            
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
                #if the unfiltered file is there (filtering got stopped mid-filter or something)
                TEC_dataAvail[i] = 4; #note that data is already downloaded and converted, but not filtered
                #use case if download and conversion success but failure on filter - big enough files to warrant this care :)
            except OSError:
                print("\n==============~Warning~==============");
                print("Data wasn't converted to standardized format & naming scheme for {} and will be reconverted. Deleting partial file as well.\n".format(TEC_dataFilePathUnfilt[i]) );
                os.remove(TEC_dataFilePathUnfilt[i]); #delete partially downloaded file to avoid conflicts
                try: #gonna try to read the file - if we fail, it failed mid download
                    with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as testFile:
                        testFile.keys(); #tries to check out some stuff in the file
                    #END WITH
                    TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                    #use case if download success but failure on conversion - big enough files to warrant this care :)
                except OSError:
                    print("\n==============~Warning~==============");
                    print("Data wasn't downloaded fully for {} and will be redownloaded. Deleting partial file as well.\n".format(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) );
                    os.remove(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete partially downloaded file to avoid conflicts
                    TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
                #END TRYING
            #END TRYING
        else: #otherwise orig data downloaded but not converted to standard (fast) format
            try: #gonna try to read the file - if we fail, it failed mid download
                with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as testFile:
                    testFile.keys(); #tries to check out some stuff in the file
                #END WITH
                TEC_dataAvail[i] = 3; #note data is already downloaded but needs to be converted to the standardized format & naming scheme
                #use case if download success but failure on conversion - big enough files to warrant this care :)
            except OSError:
                print("\n==============~Warning~==============");
                print("Data wasn't downloaded fully for {} and will be redownloaded. Deleting partial file as well.\n".format(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig) );
                os.remove(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete partially downloaded file to avoid conflicts
                TEC_dataAvail[i] = 2; #data wasn't downloaded fully and needs to be redownloaded
            #END TRYING
        #END IF
        
    else:
        TEC_dataAvail[i] = 1; #note data is there
    #END IF  
#END FOR
    
#Error check:
if( np.any(TEC_dataAvail == -1) == 1 ):
    print('Exiting due to no available data on key days.');
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
    if( TEC_dataAvail[i] == 2 ): #2 means data needs to be downloaded
        
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
        print("Downloading {} to {}.\nFile size is {:.2f} GB. At 10 MB/s (80 Mbps) expect it to take {:.2f} min.\nFrom site {}".format( current_webFileNames,folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]), current_webFileSize/(1024**3),(current_webFileSize/1024**2*(1/10))/60, web_base_site+current_web_fileLinks ));
        tic = time.time(); #for time testing
        urlretrieve(web_base_site+current_web_fileLinks,folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + current_webFileNames,reporthook=downloadProgress); #download the file in question to the data directory
        toc = time.time() - tic; #for time testing
        print("\nTime to download: {:.2f} min\n\n".format(toc/60)); #extra space at end
        
        TEC_dataAvail[i] = 3; #move on to the next stage ;) oh yeah it feels wrong but its so easy
    #END IF
        
    #-----CONVERT THE DATA TO FASTER FORMAT-----
    if( TEC_dataAvail[i] == 3 ): #3 means data was downloaded but needs to be converted
        
        TEC_fileNameOrig = 'los_' + str(dateRange_full[i,0]) + str(dateRange_full[i,1]).zfill(2) + str(dateRange_full[i,2]).zfill(2) + '.001.h5'; #create expected file name
        print("Converting {} to {}.\nAt 18 data vectors and 90 sec per vector, expect 27 minutes for conversion to finish.\n".format(TEC_fileNameOrig,TEC_dataFileNameUnfilt[i]));
        
        tic = time.time(); #for time testing
        #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
        with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r') as TEC_fileOrig:
            #TEC_fileOrig = h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig, 'r'); #read that file
        
            # List all groups
            TEC_fileKeys = list(TEC_fileOrig.keys()); #get the folder names in the hdf5 file
            TEC_fileDatalocale = list(TEC_fileOrig[TEC_fileKeys[0]]); #choose the data folder
            TEC_fileData = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]]; #import the data shape
            dataVectCntr = 0; #counter to announce progress
            dataVectMax = 8+7+1; #total number of data vectors to import, order: int / float / string
            
            TEC_fileData_int = np.zeros( (TEC_fileData.shape[0],8),dtype='int16'); #preallocate
            TEC_fileData_int[:,0] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'sat_id'].astype("int16"); #import satellite ID
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,1] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'year'].astype("int16"); #import year
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,2] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'month'].astype("int16"); #import month
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,3] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'day'].astype("int16"); #import day
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,4] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'hour'].astype("int16"); #import hour
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,5] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'min'].astype("int16"); #import min
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,6] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'sec'].astype("int16"); #import sec
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_int[:,7] = subfun_date_to_dayNum( np.array([TEC_fileData_int[:,2],TEC_fileData_int[:,3],TEC_fileData_int[:,1]]) , 2 ); #convert to day number, option 2 only returns day number
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
                        
            TEC_fileData_float = np.zeros( (TEC_fileData.shape[0],7),dtype='float32'); #preallocate
            #0 is for de-biased TEC
            TEC_fileData_float[:,1] = TEC_fileData_int[:,7] + TEC_fileData_int[:,4]/24 + TEC_fileData_int[:,5]/1440 + TEC_fileData_int[:,6]/86400; #days, calculate hour/min/sec into days and add to the current day
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_float[:,2] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'elm']; #deg, import pierce point elevation
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_float[:,3] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'gdlat']; #arcdeg, import pierce point geodedic latitude
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_float[:,4] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'glon']; #arcdeg, import pierce point longitude
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_float[:,5] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'los_tec']; #TECU, import pierce point line-of-sight TEC
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_fileData_float[:,6] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'dlos_tec']; #TECU, import pierce point line-of-sight TEC error
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            #TEC_fileData_float[:,8] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'rec_bias']; #TECU, import receiver bias
            #TEC_fileData_float[:,9] = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'drec_bias']; #TECU, import receiver bias error
            
            TEC_fileData_pierceAlt = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'pierce_alt']; #km, import pierce point altitude
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            #pierce point altitude is the same number, no need to record it a bunch
                        
            TEC_fileData_string = TEC_fileOrig[TEC_fileKeys[0] + "/" + TEC_fileDatalocale[0]][:,'gps_site']; #import GPS site name            
            dataVectCntr = dataVectCntr + 1;
            sys.stdout.write("\rData vector {}/{} & {} min".format(dataVectCntr,dataVectMax,round((time.time()-tic)/60,2)));
            sys.stdout.flush();
        #END WITH
        
        #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
        with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'w') as TEC_fileUnfilt:
            #TEC_fileUnfilt = h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i], 'w'); #write that file new file
            TEC_file_dataSet_int = TEC_fileUnfilt.create_dataset("int", (TEC_fileData_int.shape[0],TEC_fileData_int.shape[1]), dtype='int16' ,compression="gzip"); #create dataset for the integers
            TEC_file_dataSet_float = TEC_fileUnfilt.create_dataset("float", (TEC_fileData_float.shape[0],TEC_fileData_float.shape[1]), dtype='float32' ,compression="gzip"); #create dataset for the floats
            TEC_file_dataSet_string = TEC_fileUnfilt.create_dataset("string", (TEC_fileData_string.shape[0],), dtype='S4' ,compression="gzip"); #create dataset for the strings
            sys.stdout.write("\rWriting ints to file & {} min".format(round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_file_dataSet_int[...] = TEC_fileData_int; #write that data
            sys.stdout.write("\rWriting floats to file & {} min".format(round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_file_dataSet_float[...] = TEC_fileData_float; #write that data
            sys.stdout.write("\rWriting strings to file & {} min".format(round((time.time()-tic)/60,2)));
            sys.stdout.flush();
            TEC_file_dataSet_string[...] = TEC_fileData_string; #write that data
            if( np.all( TEC_fileData_pierceAlt == TEC_fileData_pierceAlt[0]) ):
                TEC_file_dataSet_float.attrs['piercealt'] = TEC_fileData_pierceAlt[0]; #record the attribute
            else:
                print("\n==============~Warning~==============");
                print("Not all pierce point altitudes reported are the same - I don't have anything to deal with this. Continuing. They should print below:");
                print("{}".format(TEC_fileData_pierceAlt));
                TEC_file_dataSet_float.attrs['piercealt'] = TEC_fileData_pierceAlt[0]; #record the attribute
            #END IF
            sys.stdout.write("\rDone writing!\t\t\t\t\t");
            sys.stdout.flush();
        #END WITH
        os.remove(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_fileNameOrig); #delete the downloaded file now that we're done with it (save some space!)
        del TEC_fileData_int #clean the memory
        del TEC_fileData_float
        del TEC_fileData_string
        del TEC_fileData_pierceAlt
        
        toc = time.time() - tic; #for time testing
        print("\nTime to convert: {} min\n".format(round(toc/60,2))); #extra space at end   
        
        TEC_dataAvail[i] = 4; #move on to the next stage
    #END IF
    
    if( (TEC_dataAvail[i] > 4) or (TEC_dataAvail[i] < 0) ): #otherwise we got problems!
        print("\n==============ERROR==============");
        print("TEC_dataAvail was set to {} on {}/{}/{} in Y/M/D format. That's not 0, 1, 2, 3, or 4 so something real wrong. Exiting.".format(TEC_dataAvail[i],dateRange_full[i,0],dateRange_full[i,1],dateRange_full[i,2]) );
        #return("No."); #Can I? I shall.
    #END IF  
    
#END FOR
    
#==============Filter Data after download & conversion==============
for i in range(0,len(dateRange_dayNum_full[:,0])): #to download and unpack OR just unpack needed data
    
    #-----FILTER THE DATA-----
    if( TEC_dataAvail[i] == 4 ): #4 means data was downloaded and converted but needs to be filtered
        
        #TIME TO IMPORT DATA TO FILTER!
        if( np.all(np.any(dateRange_dayNum_full_orig == dateRange_dayNum_full[i,:],axis=0)) == 1 ): #if statement to find original days requested
            print("Filtering {} to {}.\n".format(TEC_dataFileNameUnfilt[i],TEC_dataFileName[i]) );
            
            tic = time.time(); #for time testing
            for j in range(-1,2): #run through the days to get the main day and the days around it
                if(TEC_dataAvail[i+j] == 4): #make sure data is there before reading it
                    with h5py.File(folder[1] + '\\' + folder_TEC + '\\' + str(dateRange_full[i,0]) + '\\' + TEC_dataFileNameUnfilt[i+j], 'r') as TEC_fileUnfilt:
                        #-----import floats (big boys)-----
                        #this is done in waves to keep the memory usage under control
                        #TEC_fileData_float_temp = TEC_fileUnfilt.get("float").value; #.value is being deppreciated apparently
                        TEC_fileData_float_temp = TEC_fileUnfilt.get("float")[()];
                        #filter for min elevation before including it into the main dataset
                        TEC_goodVals = np.where(TEC_fileData_float_temp[:,2] >= minElevation)[0]; #get the locations where min elevation is too low
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
                        del TEC_fileData_string_temp
                    #END WITH
                    
                    try: #pierce altitude
                        TEC_fileData_pierceAlt = np.append(TEC_fileData_pierceAlt,TEC_fileData_pierceAlt_temp); #if var exists, append new data on it
                    except(NameError):
                        TEC_fileData_pierceAlt = TEC_fileData_pierceAlt_temp; #if var didn't exist, time to make it exist
                    #END TRY
                                        
                else: #otherwise data isn't there - because strict padding is off and data wasn't there
                    try: #padding warning
                        if( j == -1 ): #before padding not there
                            TEC_fileData_paddingWarning.append("Before"); #append to var
                        elif( j == 1): #after padding not there
                            TEC_fileData_paddingWarning.append("After"); #append to var
                        else: #something bad
                            print("\n==============ERROR==============");
                            print("There's no data on a main day requested {}/{} (Y/#D) - which shouldn't have happened. Help. Exiting.\n".format(dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
                            #return("No");
                        #END IF
                    except(NameError):
                        if( j == -1 ): #before padding not there
                            TEC_fileData_paddingWarning = []; #prep a list
                            TEC_fileData_paddingWarning.append("Before"); #if var didn't exist, time to make it exist
                        elif( j == 1): #after padding not there
                            TEC_fileData_paddingWarning = []; #prep a list
                            TEC_fileData_paddingWarning.append("After"); #if var didn't exist, time to make it exist
                        else: #something bad
                            print("\n==============ERROR==============");
                            print("There's no data on a main day requested {}/{} (Y/#D) - which shouldn't have happened. Help. Exiting.\n".format(dateRange_dayNum_full[i,0],dateRange_dayNum_full[i,1]) );
                            #return("No");
                        #END IF
                    #END TRY
                #END IF              
            #END FOR
            
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
            print("\nTime to import all needed data into memory: {} min\n".format(round(toc/60,2))); #extra space at end  
            
            
            #==============THIS WILL GO IN A FUNCTION FOR MAX SPEED==============
            #TIME TO FILTER!
            #==============Unfiltered File Layout==============
            #Integer Layout
            #0 = Satellite ID [# that corresponds to GPS sat]
            #1 = Year timestamp [years]
            #2 = Month timestamp [months]
            #3 = Day timestamp [days]
            #4 = Hour timestamp [hrs]
            #5 = Minute timestamp [mins]
            #6 = Second timestamp [secs]
            #7 = Day Number timestamp [days]
            #
            #Float Layout
            #0 = delta-TEC "kinda de-biased TEC" [TECU]
            #1 = current time in day format [days] - does not support years
            #2 = elevation [deg]
            #3 = geodedic latitude [arcdeg]
            #4 = longitude [arcdeg]
            #5 = line-of-sight TEC [TECU]
            #6 = error in line-of-sight TEC [TECU]
            #
            #String Layout
            #[] = Receiver Site Name (there's no dim on this one) - use .decode('UTF-8') to make it a string again
            
            TEC_fileData_string_unique = np.unique(TEC_fileData_string); #get unique site names
            #TEC_fileData_float_temp = TEC_fileData_float[100000:100100];
            #print("{}".format(TEC_fileData_string_unique[0:10]));
            
            cntr = 0; #testing cntr
            tic = time.time(); #for time testing
            #cruise through every site (from the unique of TEC_fileData_string)
#            for j in range(0,len(TEC_fileData_string_unique)): 
            for j in range(0,1): #len(TEC_fileData_string_unique)
#            for j in range(0,np.int64(len(TEC_fileData_string_unique)/4)): #run through all of the sites
                #j = 1;
                
                currentSite_loc = np.where( TEC_fileData_string_unique[j] == TEC_fileData_string )[0]; #get the data locations for 1 site
                #sys.exit()
                currentSat_unique = np.unique(TEC_fileData_int[currentSite_loc,0]); #get the unique sats in at that current site
                #cruise through every sat at a site (from unique of TEC_fileData_int[siteIndex,0])
                
#                for k in range(0,len(currentSat_unique)):
                for k in range(14,len(currentSat_unique)):
                #k = 21;
                    currentSat_loc = np.where( currentSat_unique[k] == TEC_fileData_int[currentSite_loc,0] )[0]; #get the data locations for 1 sat at that 1 site
                    currentTEC = TEC_fileData_float[currentSite_loc,5][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                    currentTECerror = TEC_fileData_float[currentSite_loc,6][currentSat_loc]; #TECU, get the TEC data for 1 sat at that 1 site - uses abs illegal double indexing
                    currentDayNum = TEC_fileData_float[currentSite_loc,1][currentSat_loc]; #days, get the day data for 1 sat at that 1 site - uses abs illegal double indexing
                    currentElv = TEC_fileData_float[currentSite_loc,2][currentSat_loc];
                    
                    currentDayNum_timeSplits_loc = np.append( np.insert( np.where(np.diff(currentDayNum)*1440 > minimumTimeGap)[0]+1 , 0, 0), len(currentDayNum) ); #get the locations where new non-contiguous data occurs - also tack on 0 and end#
                    
                    for l in range(0,len(currentDayNum_timeSplits_loc)-1):
                    #l = 0;
                        #only get the data for times ranges thare are two times the filter period or longer
                        if( (currentDayNum[currentDayNum_timeSplits_loc[l+1]-1]-currentDayNum[currentDayNum_timeSplits_loc[l]])*1440 > filterPeriod*2 ):
                            #copy day and TEC data into their own vars
                            currentDayNum_singleSatLine = currentDayNum[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                            currentTECerror_singleSatLine = currentTECerror[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                            currentDayNum_singleSatLine_min = (currentDayNum_singleSatLine-np.round(np.mean(currentDayNum_singleSatLine)))*144; #min/10, convert to minutes and /10 so it's more stable when fitting a curve
                            #currentTEC_singleSatLine_zerod = currentTEC_singleSatLine - np.min(currentTEC_singleSatLine); #TECU, adjusted TECU line that is zero'd at minimum
                            currentTEC_singleSatLine = currentTEC[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]]; #this will increase readability
                            
                            currentElv_singleSatLine = currentElv[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]];
                            
                            #now that we've got a single contiguous satellite data streak for a single site (data will look like a U), I'm gonna "de-bias" it by fitting a 2nd order polynomial (y = a + b*x + c*x^2) and subtracting that polynomial from the data - that's how I'm getting my delta-TEC
                                               
                            #use this for a first guess - 4th order fit
                            #2nd order can't deal with the non-symmetricalness of the TEC data
                            #will use this to identify minimum of curve in real TEC data
                            currentPolyCoef_firstFit = np.polynomial.polynomial.polyfit(currentDayNum_singleSatLine_min, currentTEC_singleSatLine, 4); #gets the coefficients that fit an x^2 polynomial - what a function to call am I right
                            currentPolyYvals_firstFit = np.polynomial.polynomial.polyval(currentDayNum_singleSatLine_min, currentPolyCoef_firstFit); #calc the guessed values 
                            
        #                    weight_estimatedLocalMin = np.int64(np.round(np.mean(argrelmin(currentPolyYvals)[0]))); #calc local minimum of guessed data
                            #extra functions force it to be 1 value - but it always should be 1 value...
                            weight_estimatedLocalMin = np.where(np.min(currentPolyYvals_firstFit) == currentPolyYvals_firstFit)[0][0]; #calc local minimum of guessed data
                            #there should only be 1 value seriously
                            
                            weight_goalLen = np.int64(np.round(len(currentTEC_singleSatLine)*0.025)); #get 2.5% of length - will use it for weighting
                            weight_searchRange = weight_goalLen*4; #search for the min around 10% of the length (quad goalLen)
                            
                            if( (weight_searchRange < weight_estimatedLocalMin) & (weight_searchRange < (len(currentTEC_singleSatLine)-weight_estimatedLocalMin)) & (weight_searchRange != 0) ):
                                weight_estimatedLocalMin = np.where(np.min(currentTEC_singleSatLine[ (weight_estimatedLocalMin-weight_searchRange) : (weight_estimatedLocalMin+weight_searchRange) ]) == currentTEC_singleSatLine)[0][0]; #find minimum in the local range guessed
                                #update that weight_estimatedLocalMin with a better guess based on the real data
                            #else:
                            #    pass; #otherwise, give up. It's prob pretty close I won't bother
                            #END IF
                            
                            #the most powerful function I could muster
                            currentFun = lambda currentPolyCoef, currentDayNum_singleSatLine_min : currentPolyCoef[0] + currentPolyCoef[1]*(currentDayNum_singleSatLine_min-currentPolyCoef[5]) + ( (np.abs(currentPolyCoef[2])*(currentDayNum_singleSatLine_min-currentPolyCoef[5])) - np.abs(currentPolyCoef[6]*(currentDayNum_singleSatLine_min-currentPolyCoef[5])) )**2 + (currentPolyCoef[3]*(currentDayNum_singleSatLine_min-currentPolyCoef[5]))**3 + (np.abs(currentPolyCoef[4])*(currentDayNum_singleSatLine_min-currentPolyCoef[5]) )**4; #create a "lambda" function whatever that is that has a tuple for the polynomial coefficients in i
                            currentPolyCoef_init = (np.min(currentTEC_singleSatLine), -1, 0.2, -0.2, 0.2,currentDayNum_singleSatLine_min[weight_estimatedLocalMin],.2); #for abs model again                                  
        #                    #guess based on 
        #                    #0 coeff for deg 0 = moves the parabola up/down - minimum TEC value for the bottom of the parabola
        #                    #1 coeff for deg 1 = parabola curve linear bit
        #                    #2 coeff for deg 2 = parabola curve curvey pit
        #                    #3 Coeff for deg 2 abs part = parabola curve curvey pit
        #                    #const 3 = moves the parabola left/right *NOT A DEGREE COEFFICIENT*

                            #only powell could fit the data go figure
                            currentPolyCoef = minimize(funErrorHlpr,currentPolyCoef_init, args=(currentDayNum_singleSatLine_min,currentTEC_singleSatLine),method="Powell" ); #calculate the polynomial coefficients to fit the data
                            currentPolyYvals = currentFun(currentPolyCoef.x,currentDayNum_singleSatLine_min); #calculate the polynomial that was fit to the data
                            
                            currentFunRootzFinderSecond = lambda coef, xrr : 6*(coef[3]**3)*(xrr - coef[5]) + 12*((coef[4])**4)*(xrr - coef[5])**2 + 2*(np.abs(coef[2]) - (((coef[6])**2)*(xrr - coef[5]))/np.abs(coef[6]*(xrr - coef[5])) )**2 + 2*( (((coef[6])**4)*(xrr - coef[5])**2)/(np.abs(coef[6]*(xrr - coef[5]))*(((coef[6])**2)*(xrr - coef[5])**2)) - ((coef[6])**2)/np.abs(coef[6]*(xrr - coef[5])) )*(np.abs(coef[2])*(xrr - coef[5]) - np.abs(coef[6]*(xrr - coef[5])) );
                                
                            current_deltaTEC = currentTEC_singleSatLine - currentPolyYvals; #dTECU, subtract the two things, the rest is the "action noise" or whatever
                            current_deltaTEC_firstFit = currentTEC_singleSatLine - currentPolyYvals_firstFit; #dTECU, first fit attempt, it ideally won't outperform the above fit by a lot
                            
                            fitType = 'big eq';

                            #if 2nd deriv has negative values, then concavity is down which is never what we want
                            #this finds that
                            currentSecondDerivLogical = np.any(currentFunRootzFinderSecond(currentPolyCoef.x,currentDayNum_singleSatLine_min) < 0); 
                            
                            if( (np.sum(current_deltaTEC_firstFit**2) < currentPolyCoef.fun) | currentSecondDerivLogical ): #note: currentPolyCoef.fun is the same as np.sum(current_deltaTEC**2)
                                #detect some filtering issues

                                #try again with a simpler function - this one ditches the double slope in the ^2 term
                                currentFun = lambda currentPolyCoef, currentDayNum_singleSatLine_min : currentPolyCoef[0] + currentPolyCoef[1]*(currentDayNum_singleSatLine_min-currentPolyCoef[4]) + ( (np.abs(currentPolyCoef[2])*(currentDayNum_singleSatLine_min-currentPolyCoef[4])) - np.abs(currentPolyCoef[3]*(currentDayNum_singleSatLine_min-currentPolyCoef[4])) )**2; #create a "lambda" function whatever that is that has a tuple for the polynomial coefficients in i
                                currentPolyCoef_init = (np.min(currentTEC_singleSatLine), -1, 0.2, 0.2,currentDayNum_singleSatLine_min[weight_estimatedLocalMin]); #for abs model again
                                currentPolyCoef = minimize(funErrorHlpr,currentPolyCoef_init, args=(currentDayNum_singleSatLine_min,currentTEC_singleSatLine),method="Powell" ); #calculate the polynomial coefficients to fit the data
                                currentPolyYvals_redo = currentFun(currentPolyCoef.x,currentDayNum_singleSatLine_min); #calculate the polynomial that was fit to the data
                                current_deltaTEC_redo = currentTEC_singleSatLine - currentPolyYvals_redo; #dTECU, subtract the two things, the rest is the "action noise" or whatever
                                
                                if( (currentPolyCoef.fun < np.sum(current_deltaTEC**2)) & currentSecondDerivLogical ):
                                    fitType = 'reduc eq - better fit & 4th deriv wrong dir';
                                elif( (currentPolyCoef.fun < np.sum(current_deltaTEC**2)) ):
                                    fitType = 'reduc eq - better fit, 4th deriv OK';
                                else:
                                    fitType = 'reduc eq - worse fit, 4th deriv wrong dir';
                                #END IF
                                #also fit basic x^2 incase it works better s7t2 f9
                                currentPolyCoef_secOrd = np.polynomial.polynomial.polyfit(currentDayNum_singleSatLine_min, currentTEC_singleSatLine, 2); #gets the coefficients that fit an x^2 polynomial - what a function to call am I right
                                currentPolyYvals_secOrd = np.polynomial.polynomial.polyval(currentDayNum_singleSatLine_min, currentPolyCoef_secOrd); #calc the guessed values 
                                current_deltaTEC_secOrd = currentTEC_singleSatLine - currentPolyYvals_secOrd; #dTECU, subtract the two things, the rest is the "action noise" or whatever
                                if( ((np.sum(current_deltaTEC_secOrd**2)) < currentPolyCoef.fun) & currentSecondDerivLogical ): #if a basic c0 + c1*x + c2*x^2 works better, use it
                                    fitType = "reduc to only! x^2 & 4th deriv wrong dir"
                                    currentPolyYvals_redo = currentPolyYvals_secOrd; #overwrite these values
                                    current_deltaTEC_redo = current_deltaTEC_secOrd;  #overwrite these values
                                elif( (np.sum(current_deltaTEC_secOrd**2)) < currentPolyCoef.fun):
                                    fitType = "reduc to only! x^2 & 4th deriv OK"
                                    currentPolyYvals_redo = currentPolyYvals_secOrd; #overwrite these values
                                    current_deltaTEC_redo = current_deltaTEC_secOrd;  #overwrite these values
                                #END IF
                                
                                #use this if error high because it's likely prev. function did something non-parabolic
                                #basically, if the new error square sum is less than the old, use it! Very rarely will both perform poorly
                                if( ((np.sum(current_deltaTEC_redo**2)) < np.sum(current_deltaTEC**2)) | currentSecondDerivLogical ): #note: now currentPolyCoef.fun is the same as np.sum(current_deltaTEC**2) where current_deltaTEC uses the NEW currentPolyYvals                                
                                    
                                    currentPolyYvals = currentPolyYvals_redo; #overwrite for plotting later
                                     
                                    #later to help with plotting
                                    current_deltaTEC = current_deltaTEC_redo; #overwrite these values
                                #END IF
                            #END IF
                            
                            #windowLen = np.int64(np.round(len(currentTEC_singleSatLine)/4)) - (np.mod(np.int64(np.round(len(currentTEC_singleSatLine)/4)),2) - 1);
                            fitType = 'Sav-Gol';
                            windowLen = 21;
                            windowLen = np.int64(60*60/30+1); #window length, 60 minutes, converted to seconds, and divided by sample rate (plus 1 because odd is required) gets 121 for an hour window length (frame length)
                            #from conversations with AC
                            
                            if( len(currentTEC_singleSatLine) >= windowLen ):
                                currentPolyYvals = savgol_filter(currentTEC_singleSatLine,windowLen,1 ); #filter it up
                            else:
                                currentPolyYvals = np.zeros(len(currentTEC_singleSatLine)); #no filter, just filler
                                fitType = fitType + ' DATA LEN SMALL NOT FILTERED'; #report nothing happened
                            #END IF
                            current_deltaTEC = currentTEC_singleSatLine - currentPolyYvals;
                            
                            
                            #if( np.any( np.abs(current_deltaTEC) > 1.5) ):
                            plt.figure();
                            plt.scatter( currentDayNum_singleSatLine , currentTEC_singleSatLine , 20 , "r" );
                            plt.scatter( currentDayNum_singleSatLine , currentPolyYvals , 20 );
                            plt.title("Fit "+fitType+": Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
                            plt.xlabel('Time [days]');
                            plt.ylabel('los TEC [TECU]');
                            plt.show();
                            
                            plt.figure();
                            plt.scatter( currentDayNum_singleSatLine , currentTEC_singleSatLine , 20 , "r" );
                            plt.scatter( currentDayNum_singleSatLine , currentPolyYvals_firstFit , 20 );
                            plt.title("Fit polynomail 1st: Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
                            plt.xlabel('Time [days]');
                            plt.ylabel('los TEC [TECU]');
                            plt.show();
                            
                            plt.figure();
                            plt.scatter( currentDayNum_singleSatLine , currentElv_singleSatLine );
                            plt.title("Elevation Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
                            plt.xlabel('Time [days]');
                            plt.ylabel('Elevation [deg]');
                            plt.show();
                            
                            plt.figure();
                            plt.scatter( currentDayNum_singleSatLine , current_deltaTEC );
                            plt.title("Delta TEC Data for Site "+TEC_fileData_string_unique[j].decode('UTF-8')+" and sat "+str(currentSat_unique[k])+" and "+str(l+1)+" track");
                            plt.xlabel('Time [days]');
                            plt.ylabel('delta los TEC [TECU]');
                            plt.show();
                            cntr = cntr + 1;
                            if( cntr == 11 ):
                                sys.exit();
                            #END IF
                            #END IF
                            
                            TEC_fileData_float[currentSite_loc[currentDayNum_timeSplits_loc[l]:currentDayNum_timeSplits_loc[l+1]],0] = current_deltaTEC; #save that deltaTEC into wherever it should go - double illegal indexing didn't work but moving it in was the same indexing so safe
                        #END IF
                    #END FOR l
                #END FOR k
            #ENF FOR j
            
            toc = time.time() - tic; #for time testing
            print("\nTime to filter all data: {} min\n".format(round(toc/60,2))); #extra space at end  
            
            #==============Filtered File Layout==============
            #Integer Layout
            #0 = Satellite ID [# that corresponds to GPS sat]
            #1 = Year timestamp [years]
            
            #Float Layout
            #0 = delta-TEC "kinda de-biased TEC" [TECU]
            #1 = current time in day format [days] - does not support years
            #2 = geodedic latitude [arcdeg]
            #3 = longitude [arcdeg]
            #4 = elevation [deg]
            
            #String Layout
            #[] = Receiver Site Name (there's no dim on this one) - use .decode('UTF-8') to make it a string again
            
            #-----TIME TO SAVE FILTERED DATA!-----
            #with h5py.File(TEC_dataFilePath[i], 'w') as TEC_file:
            #    pass;
            #END WITH
            
        #END IF
        
        #TEC_dataAvail[i] = 1; #data filtered and done! prob don't need to set this
    #END IF  
    
    
    if( (FLG_deleteUnfilt == 1) and (TEC_dataAvail[i] == 4) ): #if delete unfiltered data is on - delete it
        os.remove(TEC_dataFilePathUnfilt[i]); #delete the unfiltered file instead of keeping it for more filtering stuff (saves hard drive space)
    #END IF
    
#END FOR i
    
    
#if( os.path.isdir(folder[1] += folder_TEC += str()) == 0 ): #check if starting year folder exists
#    #if not, make it
#    os.makedirs(folder[1] + folder_TEC)








#    return;
