#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
get_airs(dates, dataDir, madrigalInst, madrigalExp, madrigalLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False);
-  dates = [[2012,12,12],[2012,12,13],[2012,12,15]]; # yr/mo/day; can also use yr/day# format [[2012,333],[2012,334],[2012,336]]
-  dataDir = './path/where/you/want/data' the code will put year folders in this directory for all years requested, and per-year data in each year directory
-  NASAearthdataLogin = {'fullname': 'First+Last', # You must convert spaces to +
                         'email': 'your@email.here',
                         'affiliation': 'Some+Where'};
**  web_retryMax: Integer number of retries to attempt when accessing the web; setting to -1 will retry forever
**  web_retryWait: Can be an integer# to wait # sec between tries, or list/tuple [startwait, maxwait, incrementwait] where every retry the code increments the startwait by incrementwait until maxwait is reached/exceeded
**  FLG_datesRange: Set to True to indicate that dates is an inclusive range in the form yr/mo/day [[2012,12,12],[2012,12,28]] or yr/day# [[2012,333],[2012,356]] - code will automatically fill in the dates in between
**  FLG_overwrite: Set to True to redownload files even if they are already downloaded (overwrite them)
**  FLG_verbose: 2 = Info and above, 1 = Warning and above, 0 = Errors only
"""

import numpy as np
from TAS.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
from TAS.subfun_date_to_dayNum import subfun_date_to_dayNum
# from TAS.subfun_dayNum_to_date import subfun_dayNum_to_date
from TAS.get_login import get_login
from TAS.get_got import get_got_webpg, get_got_webfile
from os import path as ospath, makedirs as osmakedirs
from shutil import move as shutilmove, get_terminal_size as shutil_get_terminal_size
from tempfile import gettempdir
from sys import stdout as sysstdout
from re import findall as refindall, compile as recompile
import pickle as pkl
import zlib

def get_airs(dates, dataDir, NASAearthdataLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False, FLG_disableHave=False, FLG_verbose=2):
    # Constants
    web_base_site = 'https://airsl2.gesdisc.eosdis.nasa.gov/opendap/Aqua_AIRS_Level2/AIRS2RET.7.0/'; # Website to be used
    
    # Prep retry wait time
    if( not isinstance(web_retryWait, (list, tuple)) ):
        web_retryWait = [web_retryWait, web_retryWait+1, 0]; # Set so the value to retry is static (always less than itself +1 and incremented by 0)
    # END IF
        
    # --- Make sure dates are good ---
    if( isinstance(dates, (list,tuple)) ):
        dates = np.asarray(dates); # Convert to a numpy array
    # END IF
    if( FLG_datesRange ):
        dates = subfun_dateORdayNum_to_fullRange(dates)[0]; # Convert to the full range from end pts (inclusive)
    # END IF
    if( dates.shape[1] == 3 ):
        # dates_yrmd = dates.copy(); # Keep this
        dates = subfun_date_to_dayNum(dates); # Convert to yr/dayNum format if in yr/mo/day
    # if( dates.shape[1] == 2 ):
    #     # dates_yrdn = dates.copy(); # Keep this
    #     dates = subfun_dayNum_to_date(dates); # Convert to yr/mo/day format if in yr/dayNum
    # END IF    
    
    dates_yrs = np.unique(dates[:, 0]); # Get the unique years
    web_fileInfo = {str(dates_yrs[i]):{} for i in range(0, dates_yrs.size)}; # Record info into a dict
    # Get metadata file if available and record the dates already requested
    for i in range(0, dates_yrs.size):
        web_metaFile_path = ospath.join(dataDir, str(dates_yrs[i]), '.metadata.pkl');
        if( ospath.isfile(web_metaFile_path) ):
            with open(web_metaFile_path, 'rb') as pklz:
                pkl_smoosh = pklz.read(); # Read in the compressed pickle file
            # END WITH
            pkl_obj = zlib.decompress(pkl_smoosh); # Decompress the pickle file
            web_metaFile = pkl.loads(pkl_obj);  # Convert the pickle file into a usable variable
            if( 'date requested' in web_metaFile[str(dates_yrs[i])] ):
                web_fileInfo[str(dates_yrs[i])]['date requested'] = web_metaFile[str(dates_yrs[i])]['date requested'].copy(); # Copy over the saved requested dates
            # END IF
        # END IF
    # END FOR i
    # Use metadata file to remove dates that have got-got
    if( FLG_overwrite == False ): # Only remove dates if overwrite is off
        dates_del = np.zeros(dates.shape[0], dtype=np.bool_); # Prep
        for i in range(dates.shape[0]-1, -1, -1):
            if( 'date requested' in web_fileInfo[str(dates[i, 0])] ):
                k = (dates[:, 0] == dates[i, 0]) & (dates[:, 1] == dates[i, 1]); # Get which were requested
                k_w = np.where(k)[0]; # Get which were requested
                dates_alreadyHave = (web_fileInfo[str(dates[i, 0])]['date requested'][:, None] == dates[k, :]).all(axis=2).any(axis=0); # See if we already have requested and got any dates worth of data
                dates = np.delete(dates, k_w[dates_alreadyHave], axis=0); # Delete data that is already got-got
                if( k.sum() == 0 ):
                    dates_del[i] = True; # Set to delete
                    web_fileInfo.pop(str([str(dates[i, 0])])); # Delete from web_fileInfo dict too
                # END IF
            # END IF
        # END FOR i
        if( np.any(dates_del) ):
            dates = dates[~dates_del, :]; # Delete years with no dates left
            dates_yrs = np.unique(dates[:, 0]); # Get the unique years again
        # END IF
    # END IF
    
    # --- Get the available data from Madrigal ---
    if( dates.shape[0] > 0 ): # Only work if we need the data        
        # # Build the specific page variables
        # web_base = 'something specific'; #year goes after this
        # web_baseAfter = '/kindat/'+str(exp)+'/format/hdf5/'; #this goes after year.
        
        # Activate the proxy/login
        get_login(NASAearthdataLogin['user'], NASAearthdataLogin['pass'], cookiez=None);
        
        # First - need to check if data is availiable for the days requested. AIRS goes daynum-by-daynum
        # web_fileMatches = []; # Prep
        for i in range(0, dates.shape[0]):
            # Madrigal lets you browse the FTP data availability by year, so get per-year info that can be applied for later
            web_fileSearch = web_base_site + str(dates[i,0]) + '/' + str(dates[i,1]) + '/'; #build site to go to
            # Get the webpage info
            rendered_content = get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
            # rendered_content = pd.read_html(StringIO( get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0]; # Get the table from the HTML page            
            
            web_fileMatches_now = refindall(r'"name"\: ".*\.hdf"', rendered_content); # Get the matches
            web_fileInfo[str(dates[i, 0])][str(dates[i, 1])] = {'file name': [], 'file link': []}; # Prep a day
            # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'] = []; # Record file name
            # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'] = []; # Record file link
            tmp_idx2del = []; # Prep
            for j in range(0, len(web_fileMatches_now)):
                if( not 'not for science use' in web_fileMatches_now[j] ):
                    # Extract the file name
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( web_fileMatches_now[j][9:-1]+'.dap.nc4');
                    # Build the link
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch + web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][-1] );
                else:
                    tmp_idx2del.append(j); # Delete
                # END IF
            # END FOR j
            # for j in tmp_idx2del[::-1]:
            #     web_fileMatches_now.pop(j); # Delete
            # # END FOR j
            # web_fileMatches.extend(web_fileMatches_now); # Tack the list onto the gigalist
        # END FOR i
        
        # # Second - Compare dates available with the requested range to determine if want to download [per day, not needed]
        # web_fileGet = {keyz:[False for _ in range(0, len(web_fileInfo[keyz]['date start']))] for keyz in web_fileInfo.keys()}; # Record info into a dict
        # for i in range(0, dates.shape[0]):
        #     # Get the date range for the data
        #     dates_start = np.datetime64('-'.join(strang.zfill(2) for strang in resub(r'\s+', '-', str(dates[i, :])[1:-1]).split('-'))); # Get the start date (hooooooboy on those python shennanigans)
        #     dates_end = dates_start + np.timedelta64(1, 'D'); # By definition the end date is just before the beginning of the next day, represent via non-inclusive comparisons (e.g., no 'or equals')
            
        #     # Check if the data has that date range
        #     for j in range(0, len(web_fileGet[str(dates[i, 0])])):
        #         web_dateStart = web_fileInfo[str(dates[i, 0])]['date start obj'][j].astype('datetime64[D]'); # Round down to nearest day
        #         web_dateEnd = (web_fileInfo[str(dates[i, 0])]['date end obj'][j] + np.timedelta64(1, 'D')).astype('datetime64[D]'); # Round up to nearest day (by adding a day and forcing to days) - non-inclusive here since ends on previous day just before the next
                
        #         web_fileGet[str(dates[i, 0])][j] |= (web_dateStart <= dates_start) & (web_dateEnd >= dates_end); # Record to include it with an or statement
        #     # END FOR j
        # # END FOR i
        
        # Third - Download if not already available
        web_dlSpeed = None; # Prep DL speed calcs
        if( FLG_verbose >= 2 ):
            colWidth = shutil_get_terminal_size().columns//2; # Get half the column width
        # END IF
        for i in range(0, dates.shape[0]):
            web_fileNum = len(web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name']); # Get the number of files for these dates
            if( web_fileNum > 0 ):
                # Make sure year/day data folder there
                dataDir_now = ospath.join(dataDir, str(dates[i, 0]), str(dates[i,1]));
                if( ospath.isdir(dataDir_now) == False ):
                    osmakedirs(dataDir_now); # Make the day folder if it doesn't exist
                # END IF
                # Prep tmp dir so that files that fail mid-download are not stored in the dataset directly
                dataDir_tmp = gettempdir();
                for j in range(0, web_fileNum):
                    web_filePath = ospath.join(dataDir_now, web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][j]); # Build the full file path
                    web_fileTmp = ospath.join(dataDir_tmp, web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][j]); # Build the temporary file path
                    if( (ospath.isfile(web_filePath) == False) or (FLG_overwrite == True) ):
                        # If the file isn't already there, get it
                        web_dlSpeed = get_got_webfile(web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'][j], web_fileTmp, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_showDownloadedAmnt=False, FLG_verbose=FLG_verbose);
                        # If that suceeds, we make it to the move
                        shutilmove(web_fileTmp, web_filePath); # Move the file from temp to its real destination
                    # END IF
                    if( FLG_verbose >= 2 ):
                        donezo = round((j+1)/web_fileNum*colWidth); # Calc how close to being done
                        sysstdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
                        sysstdout.flush();
                    # END IF
                # END FOR j
            else:
                if( FLG_verbose >= 1 ):
                    print('WARNING in get_airs.py: Requested yr/dayNum '+str(dates[i, 0])+'/'+str(dates[i, 1])+' data available.');
                # END IF
            # END IF
        # END FOR i
        
        # Record file about what was requested
        for keyz in web_fileInfo.keys():
            if( ospath.isdir( ospath.join(dataDir, keyz) ) ): # Only make a metadata file if the year is there
                # Metadata file path
                web_metaFile_path = ospath.join(dataDir, keyz, '.metadata.pkl');
                k = np.where(dates[:, 0] == np.int16(keyz))[0]; # Get which were requested
                if( 'date requested' in web_fileInfo[keyz] ):
                    web_metaFile_reqMissing = np.logical_not( (dates[k, :] == web_fileInfo[keyz]['date requested'][:, None]).all(2).any(0) ); # Get dates missing from web_metaFile_req
                    web_fileInfo[keyz]['date requested'] = np.append(web_fileInfo[keyz]['date requested'], dates[web_metaFile_reqMissing, :], axis=0); # Put in the missing dates
                    # Sort by year
                    k = np.argsort( web_fileInfo[keyz]['date requested'][:, 0] );
                    web_fileInfo[keyz]['date requested'] = web_fileInfo[keyz]['date requested'][k, :]; # Appy sort
                    # Sort by day in year
                    for kj in range(0, dates_yrs.size):
                        kf = web_fileInfo[keyz]['date requested'][:, 0] == dates_yrs[kj]; # Get where we're working
                        k[kf] = np.argsort( web_fileInfo[keyz]['date requested'][:, 1] );
                    # END FOR kj
                    web_fileInfo[keyz]['date requested'] = web_fileInfo[keyz]['date requested'][k, :]; # Appy sort
                else:
                    web_fileInfo[keyz]['date requested'] = dates[k, :].copy();
                # END IF
                
                # Save the metadata as a file in the year folder
                pkl_obj = pkl.dumps(web_fileInfo, protocol=pkl.HIGHEST_PROTOCOL); # Pickle object
                pkl_smoosh = zlib.compress(pkl_obj); # Compress the object
                with open(web_metaFile_path, 'wb') as pklz:
                    pklz.write(pkl_smoosh); # Write the compressed pickle to disk
                # END WITH
            # END IF
        # END FOR keyz
    else:
        if( (FLG_disableHave == False) and (FLG_verbose >=2) ):
            print('INFO in get_airs.py: Data has been downloaded already for all dates requested in `'+dataDir+'`.');
        # END IF
    # END IF
# END DEF