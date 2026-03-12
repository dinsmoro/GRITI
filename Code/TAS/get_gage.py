#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
!!! gage is canceled for now because it is a very cursed website !!!
get_gage(dates, dataDir, madrigalInst, madrigalExp, madrigalLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False);
-  dates = [[2012,12,12],[2012,12,13],[2012,12,15]]; # yr/mo/day; can also use yr/day# format [[2012,333],[2012,334],[2012,336]]
-  dataDir = './path/where/you/want/data' the code will put year folders in this directory for all years requested, and per-year data in each year directory
-  earthscopeLogin = {'fullname': 'First+Last', # You must convert spaces to +
                      'email': 'your@email.here',
                      'affiliation': 'Some+Where'};
**  statLocNrange: [lat in cdeg, long in cdeg, range in km] list or tuple, disabled if set to None. Does a radius around a point, gets all stations in that radius. Yields priority to statList if statList is also not set to None. 
**  statList: [list, of, stations] list or tuple of station names, disabled if set to None. Directly request stations, has priority over statLocNrange if statLocNrange is also not set to None.
    **if both statLocNrange and statList are set to None, code will get ALL stations**
**  web_retryMax: Integer number of retries to attempt when accessing the web; setting to -1 will retry forever
**  web_retryWait: Can be an integer# to wait # sec between tries, or list/tuple [startwait, maxwait, incrementwait] where every retry the code increments the startwait by incrementwait until maxwait is reached/exceeded
**  FLG_datesRange: Set to True to indicate that dates is an inclusive range in the form yr/mo/day [[2012,12,12],[2012,12,28]] or yr/day# [[2012,333],[2012,356]] - code will automatically fill in the dates in between
**  FLG_overwrite: Set to True to redownload files even if they are already downloaded (overwrite them)
**  FLG_verbose: 2 = Info and above, 1 = Warning and above, 0 = Errors only
"""

import numpy as np
import pandas as pd
from TAS.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
from TAS.subfun_date_to_dayNum import subfun_date_to_dayNum
# from TAS.subfun_dayNum_to_date import subfun_dayNum_to_date
from TAS.get_login import get_login
from TAS.get_got import get_got_webpg, get_got_webfile
from TAS.positizer import positizer_dist
from os import path as ospath, makedirs as osmakedirs
from shutil import move as shutilmove, get_terminal_size as shutil_get_terminal_size
from tempfile import gettempdir
import sys
from re import findall as refindall, compile as recompile, match as rematch
from io import StringIO
import pickle as pkl
import zlib

def get_gage(dates, dataDir, earthscopeLogin, statLocNrange=None, statList=None, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False, FLG_parallel=True, FLG_verbose=2):
    # Constants
    web_base_site = 'https://gage-data.earthscope.org/archive/gnss/rinex/obs/'; # Website to be used
    
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
    for i in range(0, dates.shape[0]):
        web_metaFile_path = ospath.join(dataDir, str(dates[i, 0]), '.metadata.pkl');
        if( ospath.isfile(web_metaFile_path) ):
            with open(web_metaFile_path, 'rb') as pklz:
                pkl_smoosh = pklz.read(); # Read in the compressed pickle file
            # END WITH
            pkl_obj = zlib.decompress(pkl_smoosh); # Decompress the pickle file
            web_metaFile = pkl.loads(pkl_obj);  # Convert the pickle file into a usable variable
            web_fileInfo[str(dates[i, 0])]['date requested'] = web_metaFile['date requested'].copy(); # Copy over the saved requested dates
        # END IF
    # END FOR i
    # Use metadata file to remove dates that have got-got
    if( FLG_overwrite == False ): # Only remove dates if overwrite is off
        dates_del = np.zeros(dates.shape[0], dtype=np.bool_); # Prep
        for i in range(0, dates.shape[0]):
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
    
    # --- Get the available data from gage ---
    if( dates.shape[0] > 0 ): # Only work if we need the data
        # Get the location of the stations
        if( ospath.isfile(ospath.join(dataDir, '.metadatagage.pklz')) ):
            with open(ospath.join(dataDir, '.metadatagage.pklz'), 'rb') as pklz:
                pkl_smoosh = pklz.read(); # Read in the compressed pickle file
            # END WITH
            pkl_obj = zlib.decompress(pkl_smoosh); # Decompress the pickle file
            statLocs = pkl.loads(pkl_obj);  # Convert the pickle file into a usable variable
            statLocs_stats = list(statLocs.keys()); # Get the keyz
        else:
            statLocs = None; # Set to none if not there; only needed for location-based searches
            statLocs_stats = None;
        # END IF
        
        # Get the stations to use
        if( (FLG_verbose >= 1) and (statLocNrange is not None) and (statList is not None) ):
            print('WARNING in get_gage.py: statLocNrange AND statList provided. Using statList as it is more specific, ignoring statLocNrange input.');
        # END IF
        if( statList is not None ):
            if( not isinstance(statList, (list, tuple)) ):
                statList = [statList]; # List it up for easy
            # END IF
            statLocs2use = statList; # Use it
        elif( statLocNrange is not None ):
            # Assume flat radial distance ignoring altitude
            statLocNrange[-1] *= 1000; # Convert to meters
            statLocs2use = []; # Prep up
        else:
            statLocs2use = 'all'; # Use em all
        # END IF
    
        # # Build the specific page variables
        # web_base = 'something specific'; #year goes after this
        # web_baseAfter = '/kindat/'+str(exp)+'/format/hdf5/'; #this goes after year.
        
        # Activate the proxy/login
        get_login(earthscopeLogin['user'], earthscopeLogin['pass'], cookiez=None);
        
        # First - need to check if data is availiable for the days requested. AIRS goes daynum-by-daynum
        if( FLG_verbose >= 2 ):
            print('INFO in get_gage.py: Getting station availability information. This takes many webpage requests, so it\'ll take a while, sorry.');
            colWidth = shutil_get_terminal_size().columns//2; # Get half the column width
        # END IF
        regexr_siteName = recompile(r'^[a-zA-Z0-9]{4,5}/'); # Compile ahead of time
        things2get = 0; # Prep
        breakpoint()
        for i in range(0, dates.shape[0]):
            # gage lets you browse the FTP data availability by year, so get per-year info that can be applied for later
            web_fileSearch = web_base_site + str(dates[i,0]).zfill(4) + '/' + str(dates[i,1]).zfill(3) + '/'; #build site to go to
            # Get the webpage info
            # rendered_content = get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
            rendered_content = pd.read_html(StringIO( get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][2:-1]['Name'].reset_index(drop=True); # Get the table from the HTML page            
            web_fileMatches = []; # Prep
            for j in range(0, rendered_content.size):
                # Get the matches
                if( regexr_siteName.match(rendered_content[j]) is not None ):
                    web_fileMatches.append(rendered_content[j]); # Tack it on if it's a site name
                # END IF
            # END FOR j
            
            web_fileInfo[str(dates[i, 0])][str(dates[i, 1])] = {'file name': [], 'file link': [], 'station': []}; # Prep a day
            # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'] = []; # Record file name
            # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'] = []; # Record file link
            if( isinstance(statLocs2use, str) ): # Get all stations (it's set to the string 'all')                
                # tmp_idx2del = []; # Prep
                for j in range(0, len(web_fileMatches)):
                    # if( web_fileMatches[j][:-1] in statLocs_stats ):
                    # Extract the file name
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( web_fileMatches[j][:-1]+str(dates[i, 1]).zfill(3)+'0.'+str(dates[i, 0])[2:]+'d.gz');
                    # Build the link
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch + web_fileMatches[j] + web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][-1] );
                    # Record the station
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['station'].append( web_fileMatches[j] );
                    # else:
                    #     if( FLG_verbose >= 2 ):
                    #         print('INFO in get_gage.py: Station `'+web_fileMatches[j][:-1]+'` not in list of known stations. Ditching.');
                    #     # END IF
                    #     tmp_idx2del.append(j); # Tack on to delete if station isn't known (must be false positive)
                    # # END IF
                # END FOR j
                # for j in tmp_idx2del[::-1]:
                #     web_fileMatches.pop(j); # Delete
                # # END FOR j
            else: # Otherwise there are specific stations
                tmp_idx2del = []; # Prep
                for j in range(0, len(statLocs2use)):
                    if( (statLocs2use[j] in statLocs_stats) and (statLocs2use[j]+'/' in web_fileMatches) ):
                        # Extract the file name
                        web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( statLocs2use[j]+str(dates[i, 1]).zfill(3)+'0.'+str(dates[i, 0])[2:]+'d.gz');
                        # Build the link
                        web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch + statLocs2use[j] + web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][-1] );
                        # Record the station
                        web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['station'].append( statLocs2use[j] );
                    else:
                        if( FLG_verbose >= 2 ):
                            if( statLocs2use[j] not in statLocs_stats ):
                                print('\nINFO in get_gage.py: Requested station `'+statLocs2use[j]+'` not in list of known stations. Ditching.');
                            else:
                                print('\nINFO in get_gage.py: Requested station `'+statLocs2use[j]+'` not in available stations for yr/dayNum '+str(dates[i, 0])+'/'+str(dates[i, 1])+'. Ditching.');
                            # END IF
                        # END IF
                        tmp_idx2del.append(j); # Tack on to delete if station isn't known (must be false positive)
                    # END IF
                # END FOR j
                web_fileMatches = statLocs2use.copy(); # Copy over the list of stations, and then delete from them as needed (reusing var, kinda)
                for j in tmp_idx2del[::-1]:
                    web_fileMatches.pop(j); # Delete
                # END FOR j
                if( (FLG_verbose >= 1) and (len(web_fileMatches) == 0) ):
                    print('\nWARNING in get_gage.py: No data available on requested yr/dayNum '+str(dates[i, 0])+'/'+str(dates[i, 1])+'.');
                # END 
                things2get += len(web_fileMatches); # Add it on, used for progress calcs later
            # END IF
            if( FLG_verbose >= 2 ):
                donezo = round((i+1)/dates.shape[0]*colWidth); # Calc how close to being done
                sys.stdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
                sys.stdout.flush();
            # END IF
        # END FOR i
        
        
        # # Second - Compare dates available with the requested range to determine if want to download [per day, not needed]
        # web_fileGet = {keyz:[False for _ in range(0, len(web_fileInfo[keyz]['date start']))] for keyz in web_fileInfo.keys()}; # Record info into a dict
        # # Get the date range for the data
        # dates_start = np.datetime64('-'.join(strang.zfill(2) for strang in resub(r'\s+', '-', str(dates[i, :])[1:-1]).split('-'))); # Get the start date (hooooooboy on those python shennanigans)
        # dates_end = dates_start + np.timedelta64(1, 'D'); # By definition the end date is just before the beginning of the next day, represent via non-inclusive comparisons (e.g., no 'or equals')
        
        # # Check if the data has that date range
        # for j in range(0, len(web_fileGet[str(dates[i, 0])])):
        #     web_dateStart = web_fileInfo[str(dates[i, 0])]['date start obj'][j].astype('datetime64[D]'); # Round down to nearest day
        #     web_dateEnd = (web_fileInfo[str(dates[i, 0])]['date end obj'][j] + np.timedelta64(1, 'D')).astype('datetime64[D]'); # Round up to nearest day (by adding a day and forcing to days) - non-inclusive here since ends on previous day just before the next
            
        #     web_fileGet[str(dates[i, 0])][j] |= (web_dateStart <= dates_start) & (web_dateEnd >= dates_end); # Record to include it with an or statement
        # # END FOR j
        
        
        # Third - Download if not already available
        if( FLG_verbose >= 2 ):
            print('INFO in get_gage.py: Downloading requested files for requested dates and stations.');
        # END IF
        web_dlSpeed = None; # Prep DL speed calcs
        cntr = 1; # Prep progress bar calcs
        web_fileInfo_yr = list(web_fileInfo.keys()); # Get the year keys
        for keyz_yr in web_fileInfo_yr:
            # Prep tmp dir so that files that fail mid-download are not stored in the dataset directly
            dataDir_tmp = gettempdir();
            web_fileInfo_day = list(web_fileInfo[keyz_yr].keys()); # Get the day keys
            for keyz_day in web_fileInfo_day:
                # Make sure year/day data folder there
                dataDir_now = ospath.join(dataDir, str(keyz_yr), str(keyz_day) ); # Prep it up
                if( ospath.isdir(ospath.join(dataDir_now)) == False ):
                    osmakedirs(ospath.join(dataDir_now)); # Make the day folder if it doesn't exist
                # END IF
                web_fileInfo_stats = list(web_fileInfo[keyz_yr].keys()); # Get the station keys
                for keyz_stat in web_fileInfo_stats:
                    web_filePath = ospath.join(dataDir_now, web_fileInfo[keyz_yr][keyz_day]['file name'][j]); # Build the full file path
                    web_fileTmp = ospath.join(dataDir_tmp, web_fileInfo[keyz_yr][keyz_day]['file name'][j]); # Build the temporary file path
                    if( (ospath.isfile(web_filePath) == False) or (FLG_overwrite == True) ):
                        # If the file isn't already there, get it
                        try:
                            web_dlSpeed = get_got_webfile(web_fileInfo[keyz_yr][keyz_day]['file link'][j], web_fileTmp, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_verbose=FLG_verbose);
                            # If that suceeds, we make it to the move
                            shutilmove(web_fileTmp, web_filePath); # Move the file from temp to its real destination
                        except:
                            if( FLG_verbose >= 1 ):
                                print('\nWARNING in get_gage.py: Requested station `'+web_fileInfo[keyz_yr][keyz_day]['file name'][j]+'` at link `'+web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'][j]+'` errored, likely not there. Doing a general seach for files matching the format (there\'s a "0" digit that means some level of data, allowing that to vary).');
                            # END IF
                            rendered_content = pd.read_html(StringIO( get_got_webpg(web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j], web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][2:-1]['Name'].reset_index(drop=True); # Get the table from the HTML page            
                            FLG_foundit = False; # Prep a flag
                            for jk in range(0, rendered_content.size):
                                regexr_fallback = rematch(web_fileInfo[keyz_yr][keyz_day]['station'][j]+str(dates[i, 1]).zfill(3)+r'\d\.'+str(dates[i, 0])[2:]+r'd\.gz', rendered_content[jk]); #Match
                                if( regexr_fallback is not None ):
                                    web_dlSpeed = get_got_webfile(web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j] + '/' + regexr_fallback.group(), web_fileTmp, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_verbose=FLG_verbose);
                                    # If that suceeds, we make it to the move
                                    shutilmove(web_fileTmp, web_filePath); # Move the file from temp to its real destination
                                    FLG_foundit = True; # Found it
                                    break; # Ditch at first match
                                # END IF
                            # END FOR jk
                            if( FLG_foundit == True ):
                                if( FLG_verbose >= 2 ):
                                    print('\nINFO in get_gage.py: Requested station `'+web_fileInfo[keyz_yr][keyz_day]['file name'][j]+'` was FOUND at link `'+web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j] + '/' + regexr_fallback.group()+'`.');
                                # END IF
                            else:
                                if( FLG_verbose >= 1 ):
                                    print('\nWARNING in get_gage.py: Requested station `'+web_fileInfo[keyz_yr][keyz_day]['file name'][j]+'` was not found in FTP folder `'+web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j] + '/'+'`. Moving on.');
                                # END IF
                            # END IF
                        # END TRYING
                    # END IF
                # END FOR keyz_stat
                if( FLG_verbose >= 2 ):
                    donezo = round(cntr/things2get*colWidth); # Calc how close to being done
                    sys.stdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
                    sys.stdout.flush();
                    cntr += 1; # Increment
                # END IF
            # END FOR keyz_day
        # END FOR keyz_yr
        
        # Record file about what was requested
        for keyz in web_fileInfo.keys():
            if( ospath.isdir( ospath.join(dataDir, keyz) ) ): # Only make a metadata file if the year is there
                # Metadata file path
                web_metaFile_path = ospath.join(dataDir, keyz, '.metadata.pkl');
                k = np.where(dates[:, 0] == np.int16(keyz))[0]; # Get which were requested
                if( 'date requested' in web_fileInfo ):
                    web_metaFile_reqMissing = np.logical_not( (dates[k, :] == web_fileInfo['date requested'][:, None]).all(2).any(0) ); # Get dates missing from web_metaFile_req
                    web_fileInfo['date requested'] = np.append(web_fileInfo['date requested'], dates[web_metaFile_reqMissing, :], axis=0); # Put in the missing dates
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
                    web_fileInfo['date requested'] = dates[k, :].copy();
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
        if( FLG_verbose >=2 ):
            print('INFO in get_gage.py: Data has been downloaded already for all dates requested in `'+dataDir_now+'`.');
        # END IF
    # END IF
# END DEF

def get_gage_statList(statLocs_nameNow, statLocs_webBase, regexr_start, regexr_end, regexr_lat, regexr_long, regexr_X, regexr_Y, regexr_Z, earthscopeLogin, web_retryMax, web_retryWait, FLG_verbose):                
    # Activate the proxy/login
    get_login(earthscopeLogin['user'], earthscopeLogin['pass'], cookiez=None);
    rendered_content = get_got_webpg(statLocs_webBase+statLocs_nameNow, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
    regexr_startLoc = regexr_start.search(rendered_content).start(); # Get where it starts
    regexr_endRegexr = regexr_end.search(rendered_content[regexr_startLoc:]); # Get where it ends
    rendered_content = rendered_content[regexr_startLoc:regexr_startLoc+regexr_endRegexr.end()]; # Get the loc info we want
    # Buid the regexes
    regexr_latGroup = regexr_lat.search(rendered_content).group();
    regexr_longGroup = regexr_long.search(rendered_content).group();
    regexr_Xgroup = regexr_X.search(rendered_content).group();
    regexr_Ygroup = regexr_Y.search(rendered_content).group();
    regexr_Zgroup = regexr_Z.search(rendered_content).group();
    # Get the data
    tmp_lat = regexr_latGroup[regexr_latGroup.find('=')+1:].strip(' ').split(); # Get the lat
    tmp_long = regexr_longGroup[regexr_longGroup.find('=')+1:].strip(' ').split(); # Get the long
    tmp_alt = regexr_endRegexr.group()[regexr_endRegexr.group().find('=')+1:].strip(' ').split(); # Get the alt
    tmp_X = regexr_Xgroup[regexr_Xgroup.find('=')+1:].strip(' ').split(); # Get the X coord
    tmp_Y = regexr_Ygroup[regexr_Ygroup.find('=')+1:].strip(' ').split(); # Get the Y coord
    tmp_Z = regexr_Zgroup[regexr_Zgroup.find('=')+1:].strip(' ').split(); # Get the Z coord
    # Build a dict
    statLocs_tmp = \
        {
            'lat': int(tmp_lat[0]) + int(tmp_lat[1])/60 + float(tmp_lat[2])/3600, # Get the lat
            'long': int(tmp_long[0]) + int(tmp_long[1])/60 + float(tmp_long[2])/3600, # Get the long
            'alt': float(tmp_alt[0]), # Get the alt
            'X': float(tmp_X[0]), # Get the alt
            'Y': float(tmp_Y[0]), # Get the alt
            'Z': float(tmp_Z[0]), # Get the alt
        }; # Make an entry
    # Adjust as needed
    if( tmp_lat[3] == 'S' ):
        statLocs_tmp['lat'] *= -1; # Flip to the southern hemisphere
    # END IF
    if( tmp_long[3] == 'W' ):
        statLocs_tmp['long'] *= -1; # Flip to the eastern hemisphere
    # END IF
    if( tmp_alt[1] == 'km' ):
        statLocs_tmp['alt'] *= 1000; # Convert to m if in km
    # END IF
    if( tmp_X[1] == 'km' ):
        statLocs_tmp['X'] *= 1000; # Convert to m if in km
    # END IF
    if( tmp_Y[1] == 'km' ):
        statLocs_tmp['Y'] *= 1000; # Convert to m if in km
    # END IF
    if( tmp_Z[1] == 'km' ):
        statLocs_tmp['Z'] *= 1000; # Convert to m if in km
    # END IF
    
    return statLocs_tmp
# END DEF

# --- For multiprocessing ---
from multiprocessing import Pool, cpu_count
from os import nice as osnice
from platform import system as platformsystem

def parallel_starmap_helper(fn, args, kwargs=None):
    #for multiprocess starmap with or without kwargs
    # MUST be outside of the function that calls it or it just hangs
    if( kwargs is None ):
        return fn(*args)
    else:
        return fn(*args, **kwargs)
    # END IF
# END DEF
    
def beNice(niceness=0): # Make nice level niceness if it's not already
    if( niceness != 0 ):
        nice = osnice(0) # Add zero to nice and return it 
        if nice < niceness: 
            osnice(niceness-nice) # Increment nice to niceness
        # END IF
    # END IF
# END DEF

def parallel_orchastrator(parallel_list, parallel_CPUnum=None, parallel_CPUmargin=0, parallel_niceNess=0): # parallel_CPUmargin is in whole percentages like 10 == 10%, NOT 0.1 == 10%
    # Lowers the amount of parallel helper code needed in the main function
    
    # --- Prep Parallel ---
    # Get the CPU num to use
    coresInstalled = cpu_count(); # Get the number of CPU cores here (probably includes hyperthreads on Intel systems and SMT threads on AMD systems)
    if( parallel_CPUnum is None ):
        if( parallel_CPUmargin == 0 ):
            parallel_CPUnum = coresInstalled; # Use multiprocessing to get the cpu count, it includes tiny cores and does not figure out who is who yolo
        else:
            parallel_CPUnum = int(np.ceil(coresInstalled*(100-parallel_CPUmargin)/100));
            if( (parallel_CPUnum == coresInstalled) and (parallel_CPUnum > 1) ):
                parallel_CPUnum -= 1; # Make sure one core is free if the % didn't manage to get a core free - margin implies you don't want ALL
            # END IF
        # END IF
    else:
        if( coresInstalled < parallel_CPUnum ):
            # from warnings import warn
            # warn('WARNING: Number of parallel processes requested `'+str(parallel_CPUnum)+'` exceeds the number of CPU cores found `'+str(coresInstalled)+'`. Limiting to '+str(coresInstalled)+' parallel processes to not oversubscribe the CPU.');
            parallel_CPUnum = coresInstalled; # Directly limit so cores not oversubscribed
        # END IF
        if( parallel_CPUmargin != 0 and (parallel_CPUnum > int(np.ceil(coresInstalled*(100-parallel_CPUmargin)/100))) ):
            parallel_CPUnum = int(np.ceil(coresInstalled*(100-parallel_CPUmargin)/100)); # Limit it to enforce the margin requested
            if( (parallel_CPUnum == coresInstalled) and (parallel_CPUnum > 1) ):
                parallel_CPUnum -= 1; # Make sure one core is free if the % didn't manage to get a core free - margin implies you don't want ALL
            # END IF
        # END IF
    # END IF
    if( parallel_CPUnum > len(parallel_list) ):
        parallel_CPUnum = len(parallel_list); # Limit parallel CPU number to the number of things to work on
    # END IF
        
    # --- Get if Linux to Activate Nice ---
    if( platformsystem() == 'Linux' ): # Only Linux has a working implementation of the scheduling system "nice"
        #--- Execute function on list of inputs ---
        with Pool(processes=parallel_CPUnum, initializer=beNice, initargs=(parallel_niceNess,)) as executor:
             results = executor.starmap(parallel_starmap_helper, parallel_list); #function you want is replaced with; parallel_starmap_helper helps starmap distribute everything right
        #END WITH
    else:
        #--- Execute function on list of inputs ---
        with Pool(processes=parallel_CPUnum) as executor:
             results = executor.starmap(parallel_starmap_helper, parallel_list); #function you want is replaced with; parallel_starmap_helper helps starmap distribute everything right
        #END WITH   
    # END WITH
    
    return results
# END DEF