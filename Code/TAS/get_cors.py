#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
get_cors(dates, dataDir, madrigalInst, madrigalExp, madrigalLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False);
-  dates = [[2012,12,12],[2012,12,13],[2012,12,15]]; # yr/mo/day; can also use yr/day# format [[2012,333],[2012,334],[2012,336]]
-  dataDir = './path/where/you/want/data' the code will put year folders in this directory for all years requested, and per-year data in each year directory
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
from TAS.get_got import get_got_webpg, get_got_webfile
from TAS.positizer import positizer_dist
from os import path as ospath, makedirs as osmakedirs, listdir as oslistdir
import sys, shutil
from re import findall as refindall, compile as recompile, match as rematch
from io import StringIO
import pickle as pkl
import zlib
# --- For Multiprocessing ---
from multiprocessing import Pool, cpu_count, current_process
from os import nice as osnice
from platform import system as platformsystem

def get_cors(dates, dataDir, statLocNrange=None, statList=None, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False, FLG_parallel=False, FLG_verbose=2):
    # Constants
    web_base_site = 'https://geodesy.noaa.gov/corsdata/rinex/'; # Website to be used
    
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
                # k_w = np.where(k)[0]; # Get which were requested
                # dates_alreadyHave = (web_fileInfo[str(dates[i, 0])]['date requested'][:, None] == dates[k, :]).all(axis=2).any(axis=0); # See if we already have requested and got any dates worth of data
                # dates = np.delete(dates, k_w[dates_alreadyHave], axis=0); # Delete data that is already got-got
                if( k.any() ):
                    dates_del[i] = True; # Set to delete
                # END IF
            # END IF
        # END FOR i
        if( np.any(dates_del) ):
            dates = dates[~dates_del, :]; # Delete years with no dates left
            dates_yrs = np.unique(dates[:, 0]); # Get the unique years again
        # END IF
    # END IF
    
    # --- Get the available data from CORS ---
    if( dates.shape[0] > 0 ): # Only work if we need the data
        # Get the location of the stations
        if( ospath.isfile(ospath.join(dataDir, '.metadataCORS.pklz')) ):
            with open(ospath.join(dataDir, '.metadataCORS.pklz'), 'rb') as pklz:
                pkl_smoosh = pklz.read(); # Read in the compressed pickle file
            # END WITH
            pkl_obj = zlib.decompress(pkl_smoosh); # Decompress the pickle file
            statLocs = pkl.loads(pkl_obj);  # Convert the pickle file into a usable variable
            statLocs_stats = list(statLocs.keys()); # Get the keyz
        else:
            if( FLG_verbose >= 2 ):
                print('INFO in get_cors.py: Generating list of station locations, it takes a bit due to the sheer number of web requests needed. Sorry.');
            # END IF
            
            # Get to figuring out where the stations are
            statLocs_links = pd.read_html(StringIO( get_got_webpg('https://geodesy.noaa.gov/corsdata/coord/', web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][2:-1]['Name'].to_list(); # Get the table from the HTML page            
            statLocs_versions = np.asarray([int(statLocs_links[i][statLocs_links[i].find('_')+1:-1]) for i in range(0, len(statLocs_links))]); # Build the years
            statLocs_version2use = np.where( np.max(statLocs_versions) == statLocs_versions )[0].item(); # Get the most recent revision
            statLocs_webBase = 'https://geodesy.noaa.gov/corsdata/coord/'+'coord_'+str(statLocs_versions[statLocs_version2use]).zfill(2)+'/';
            statLocs_names = pd.read_html(StringIO( get_got_webpg(statLocs_webBase, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][2:-1]; # Get the table from the HTML page     
            statLocs_names = statLocs_names['Name'].to_list();
            regexr_names = recompile(r'^[a-zA-Z0-9]+'); # For matching the name later, to support future more than 4 character station names (it could happen)
            regexr = recompile(r'[a-zA-Z0-9]+_'+str(statLocs_versions[statLocs_version2use]).zfill(2)+'.coord.txt'); # Prep a regexr
            for i in range(len(statLocs_names)-1, -1, -1):
                if( regexr.match(statLocs_names[i]) is None ):
                    statLocs_names.pop(i); # Ditch
                # END IF
            # END FOR i
            regexr_start = recompile(r'ITRF\d+ POSITION'); # Prep a regexr
            regexr_end = recompile(r'ellipsoid height\s+=\s+[\-0-9.]+?\s+(?:m|km)'); # Prep a regexr
            regexr_lat = recompile(r'latitude\s+=\s+[\-0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+\s+[NESW]'); # Prep a regexr
            regexr_long = recompile(r'longitude\s+=\s+[\-0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+\s+[NESW]'); # Prep a regexr
            regexr_X = recompile(r'X\s+=\s+[\-0-9.]+\s+(?:m|km)'); # Prep a regexr
            regexr_Y = recompile(r'Y\s+=\s+[\-0-9.]+\s+(?:m|km)'); # Prep a regexr
            regexr_Z = recompile(r'Z\s+=\s+[\-0-9.]+\s+(?:m|km)'); # Prep a regexr
            statLocs = {}; # Prep
            
            parallel_list = []; # Prep
            for i in range(0, len(statLocs_names)):
                parallel_list.append([get_cors_statList, [statLocs_names[i], statLocs_webBase, regexr_start, regexr_end, regexr_lat, regexr_long, regexr_X, regexr_Y, regexr_Z, web_retryMax, web_retryWait, FLG_verbose]])
            # END FOR i
            if( FLG_parallel ):
                results = parallel_orchastrator(parallel_list, parallel_CPUnum=6); # 6 threads to get the files, limit the amount of server-hammering (it throttles at 6 HARD!)
            else:
                results = []; # Prep
                for i in range(0, len(parallel_list)):
                    results.append( parallel_list[i][0](*parallel_list[i][1]) ); # Call the function without parallelization
                # END FOR i
            # END IF
            for i in range(0, len(statLocs_names)):
                statLocs[regexr_names.match(statLocs_names[i]).group()] = results[i]; # Get it out
            # END FOR i
            del results
            
            # for i in range(0, len(statLocs_names)):
            #     rendered_content = get_got_webpg(statLocs_webBase+statLocs_names[i], web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
            #     regexr_startLoc = regexr_start.search(rendered_content).start(); # Get where it starts
            #     regexr_endRegexr = regexr_end.search(rendered_content[regexr_startLoc:]); # Get where it ends
            #     rendered_content = rendered_content[regexr_startLoc:regexr_startLoc+regexr_endRegexr.end()]; # Get the loc info we want
            #     regexr_latGroup = regexr_lat.search(rendered_content).group();
            #     regexr_longGroup = regexr_long.search(rendered_content).group();
            #     regexr_Xgroup = regexr_X.search(rendered_content).group();
            #     regexr_Ygroup = regexr_Y.search(rendered_content).group();
            #     regexr_Zgroup = regexr_Z.search(rendered_content).group();
            #     tmp_lat = regexr_latGroup[regexr_latGroup.find('=')+1:].strip(' ').split(); # Get the lat
            #     tmp_long = regexr_longGroup[regexr_longGroup.find('=')+1:].strip(' ').split(); # Get the long
            #     tmp_alt = regexr_endRegexr.group()[regexr_endRegexr.group().find('=')+1:].strip(' ').split(); # Get the alt
            #     tmp_X = regexr_Xgroup[regexr_Xgroup.find('=')+1:].strip(' ').split(); # Get the X coord
            #     tmp_Y = regexr_Ygroup[regexr_Ygroup.find('=')+1:].strip(' ').split(); # Get the Y coord
            #     tmp_Z = regexr_Zgroup[regexr_Zgroup.find('=')+1:].strip(' ').split(); # Get the Z coord
            #     statLocs[regexr_names.match(statLocs_names[i]).group()] = \
            #         {
            #             'lat': int(tmp_lat[0]) + int(tmp_lat[1])/60 + float(tmp_lat[2])/3600, # Get the lat
            #             'long': int(tmp_long[0]) + int(tmp_long[1])/60 + float(tmp_long[2])/3600, # Get the long
            #             'alt': float(tmp_alt[0]), # Get the alt
            #             'X': float(tmp_X[0]), # Get the alt
            #             'Y': float(tmp_Y[0]), # Get the alt
            #             'Z': float(tmp_Z[0]), # Get the alt
            #         }; # Make an entry
            #     if( tmp_lat[3] == 'S' ):
            #         statLocs[regexr_names.match(statLocs_names[i]).group()]['lat'] *= -1; # Flip to the southern hemisphere
            #     # END IF
            #     if( tmp_long[3] == 'W' ):
            #         statLocs[regexr_names.match(statLocs_names[i]).group()]['long'] *= -1; # Flip to the eastern hemisphere
            #     # END IF
            #     if( tmp_alt[1] == 'km' ):
            #         statLocs[regexr_names.match(statLocs_names[i]).group()]['alt'] *= 1000; # Convert to m if in km
            #     # END IF
            #     if( tmp_X[1] == 'km' ):
            #         statLocs[regexr_names.match(statLocs_names[i]).group()]['X'] *= 1000; # Convert to m if in km
            #     # END IF
            #     if( tmp_Y[1] == 'km' ):
            #         statLocs[regexr_names.match(statLocs_names[i]).group()]['Y'] *= 1000; # Convert to m if in km
            #     # END IF
            #     if( tmp_Z[1] == 'km' ):
            #         statLocs[regexr_names.match(statLocs_names[i]).group()]['Z'] *= 1000; # Convert to m if in km
            #     # END IF
            # # END FOR i      
            
            # Save the metadata as a file in the year folder
            pkl_obj = pkl.dumps(statLocs, protocol=pkl.HIGHEST_PROTOCOL); # Pickle object
            pkl_smoosh = zlib.compress(pkl_obj); # Compress the object
            with open(ospath.join(dataDir, '.metadataCORS.pklz'), 'wb') as pklz:
                pklz.write(pkl_smoosh); # Write the compressed pickle to disk
            # END WITH
            statLocs_stats = list(statLocs.keys()); # Get the keyz
        # END IF
    
        # Get the stations to use
        if( (FLG_verbose >= 1) and (statLocNrange is not None) and (statList is not None) ):
            print('WARNING in get_cors.py: statLocNrange AND statList provided. Using statList as it is more specific, ignoring statLocNrange input.');
        # END IF
        if( statList is not None ):
            if( not isinstance(statList, (list, tuple)) ):
                statList = [statList]; # List it up for easy
            # END IF
            statLocs2use = statList; # Use it
            statLocs_locs = list( statLocs.keys() ); # Get em all
            # Make sure the stations are there
            for i in range(len(statLocs2use)-1, -1, -1):
                statLocs2use[i] = statLocs2use[i].lower(); # Prep it up
                if( statLocs2use[i] not in statLocs_locs ):
                    if( FLG_verbose >= 1 ):
                        print('WARNING in get_cors.py: Requested station `'+statLocs2use[i]+'` is not in the list of available stations.');
                    # END IF
                    statLocs2use.pop(i); # Ditch it
                # END IF
            # END FOR i
        elif( statLocNrange is not None ):
            # Assume flat radial distance ignoring altitude
            statLocNrange[-1] *= 1000; # Convert to meters
            statLocs2use = []; # Prep up
            for keyz in statLocs.keys():
                if( statLocNrange[-1] >= positizer_dist( statLocNrange[0:2], (statLocs[keyz]['lat'], statLocs[keyz]['long']), solverOverride='haversine' ).item() ): # Test the distance, haversine for speed
                    statLocs2use.append(keyz); # Tack on the key to use
                # END IF
            # END FOR keyz
        else:
            statLocs2use = 'all'; # Use em all
        # END IF
        
        # Check if stations requested are not in the list
        for j in range(0, len(statLocs2use)):
            if( (FLG_verbose >= 1) and (statLocs2use[j] not in statLocs_stats) ):
                print('\nWARNING in get_cors.py: Requested station `'+statLocs2use[j]+'` not in list of known stations, check spelling?');
            # END IF
        # END FOR j
        
        # # Build the specific page variables
        # web_base = 'something specific'; #year goes after this
        # web_baseAfter = '/kindat/'+str(exp)+'/format/hdf5/'; #this goes after year.
        
        # --- This structure worked but CORS is so slow at requests that checking each day was insanely slow, better to just yolo and guess at the file name instead of getting it from the day pages ---
        # --- !! you need this method for `statLocs2use == 'all'` IF you want the few stations that aren't listed in the CORS station directory. There are a few rogue ones, they need to be found on each page manually !! ---
        # # First - need to check if data is availiable for the days requested. CORS goes daynum-by-daynum
        # if( FLG_verbose >= 2 ):
        #     print('INFO in get_cors.py: Getting station availability information. This takes many webpage requests, so it\'ll take a while, sorry.');
        #     colWidth = shutil.get_terminal_size().columns//2; # Get half the column width
        # # END IF
        # regexr_siteName = recompile(r'^[a-zA-Z0-9]{4,5}/'); # Compile ahead of time
        # things2get = 0; # Prep
        # for i in range(0, dates.shape[0]):
        #     # CORS lets you browse the FTP data availability by year, so get per-year info that can be applied for later
        #     web_fileSearch = web_base_site + str(dates[i,0]).zfill(4) + '/' + str(dates[i,1]).zfill(3) + '/'; #build site to go to
        #     # Get the webpage info
        #     # rendered_content = get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
        #     rendered_content = pd.read_html(StringIO( get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][2:-1]['Name'].reset_index(drop=True); # Get the table from the HTML page            
        #     web_fileMatches = []; # Prep
        #     for j in range(0, rendered_content.size):
        #         # Get the matches
        #         if( regexr_siteName.match(rendered_content[j]) is not None ):
        #             web_fileMatches.append(rendered_content[j]); # Tack it on if it's a site name
        #         # END IF
        #     # END FOR j
            
        #     web_fileInfo[str(dates[i, 0])][str(dates[i, 1])] = {'file name': [], 'file link': [], 'station': []}; # Prep a day
        #     # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'] = []; # Record file name
        #     # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'] = []; # Record file link
        #     if( isinstance(statLocs2use, str) ): # Get all stations (it's set to the string 'all')                
        #         # tmp_idx2del = []; # Prep
        #         for j in range(0, len(web_fileMatches)):
        #             # if( web_fileMatches[j][:-1] in statLocs_stats ):
        #             # Extract the file name
        #             web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( web_fileMatches[j][:-1]+str(dates[i, 1]).zfill(3)+'0.'+str(dates[i, 0])[2:]+'d.gz');
        #             # Build the link
        #             web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch + web_fileMatches[j] + web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][-1] );
        #             # Record the station
        #             web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['station'].append( web_fileMatches[j] );
        #             # else:
        #             #     if( FLG_verbose >= 2 ):
        #             #         print('INFO in get_cors.py: Station `'+web_fileMatches[j][:-1]+'` not in list of known stations. Ditching.');
        #             #     # END IF
        #             #     tmp_idx2del.append(j); # Tack on to delete if station isn't known (must be false positive)
        #             # # END IF
        #         # END FOR j
        #         # for j in tmp_idx2del[::-1]:
        #         #     web_fileMatches.pop(j); # Delete
        #         # # END FOR j
        #     else: # Otherwise there are specific stations
        #         tmp_idx2del = []; # Prep
        #         for j in range(0, len(statLocs2use)):
        #             if( (statLocs2use[j] in statLocs_stats) and (statLocs2use[j]+'/' in web_fileMatches) ):
        #                 # Extract the file name
        #                 web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( statLocs2use[j]+str(dates[i, 1]).zfill(3)+'0.'+str(dates[i, 0])[2:]+'d.gz');
        #                 # Build the link
        #                 web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch + statLocs2use[j] + web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][-1] );
        #                 # Record the station
        #                 web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['station'].append( statLocs2use[j] );
        #             else:
        #                 if( FLG_verbose >= 2 ):
        #                     if( statLocs2use[j] not in statLocs_stats ):
        #                         print('\nINFO in get_cors.py: Requested station `'+statLocs2use[j]+'` not in list of known stations. Ditching.');
        #                     else:
        #                         print('\nINFO in get_cors.py: Requested station `'+statLocs2use[j]+'` not in available stations for yr/dayNum '+str(dates[i, 0])+'/'+str(dates[i, 1])+'. Ditching.');
        #                     # END IF
        #                 # END IF
        #                 tmp_idx2del.append(j); # Tack on to delete if station isn't known (must be false positive)
        #             # END IF
        #         # END FOR j
        #         web_fileMatches = statLocs2use.copy(); # Copy over the list of stations, and then delete from them as needed (reusing var, kinda)
        #         for j in tmp_idx2del[::-1]:
        #             web_fileMatches.pop(j); # Delete
        #         # END FOR j
        #         if( (FLG_verbose >= 1) and (len(web_fileMatches) == 0) ):
        #             print('\nWARNING in get_cors.py: No data available on requested yr/dayNum '+str(dates[i, 0])+'/'+str(dates[i, 1])+'.');
        #         # END 
        #         things2get += len(web_fileMatches); # Add it on, used for progress calcs later
        #     # END IF
        #     if( FLG_verbose >= 2 ):
        #         donezo = round((i+1)/dates.shape[0]*colWidth); # Calc how close to being done
        #         sys.stdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
        #         sys.stdout.flush();
        #     # END IF
        # # END FOR i
        
        # # Second - Download if not already available
        # if( FLG_verbose >= 2 ):
        #     print('INFO in get_cors.py: Downloading requested files for requested dates and stations.');
        # # END IF
        # web_dlSpeed = None; # Prep DL speed calcs
        # cntr = 1; # Prep progress bar calcs
        # web_fileInfo_yr = list(web_fileInfo.keys()); # Get the year keys
        # for keyz_yr in web_fileInfo_yr:
        #     web_fileInfo_day = list(web_fileInfo[keyz_yr].keys()); # Get the day keys
        #     for keyz_day in web_fileInfo_day:
        #         # Make sure year/day data folder there
        #         dataDir_now = ospath.join(dataDir, str(keyz_yr), str(keyz_day) ); # Prep it up
        #         if( ospath.isdir(ospath.join(dataDir_now)) == False ):
        #             osmakedirs(ospath.join(dataDir_now)); # Make the day folder if it doesn't exist
        #         # END IF
        #         web_fileInfo_stats = list(web_fileInfo[keyz_yr].keys()); # Get the station keys
        #         web_fileSearch = web_base_site + str(keyz_yr).zfill(4) + '/' + str(keyz_day).zfill(3) + '/'; #build site to go to
        #         for keyz_stat in web_fileInfo_stats:
        #             web_filePath = ospath.join(dataDir_now, web_fileInfo[keyz_yr][keyz_day]['file name'][j]); # Build the full file path
        #             if( (ospath.isfile(web_filePath) == False) or (FLG_overwrite == True) ):
        #                 # If the file isn't already there, get it
        #                 try:
        #                     web_dlSpeed = get_got_webfile(web_fileInfo[keyz_yr][keyz_day]['file link'][j], web_filePath, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_verbose=FLG_verbose);
        #                 except:
        #                     if( FLG_verbose >= 1 ):
        #                         print('\nWARNING in get_cors.py: Requested station `'+web_fileInfo[keyz_yr][keyz_day]['file name'][j]+'` at link `'+web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'][j]+'` errored, likely not there. Doing a general seach for files matching the format (there\'s a "0" digit that means some level of data, allowing that to vary).');
        #                     # END IF
        #                     rendered_content = pd.read_html(StringIO( get_got_webpg(web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j], web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][2:-1]['Name'].reset_index(drop=True); # Get the table from the HTML page            
        #                     FLG_foundit = False; # Prep a flag
        #                     for jk in range(0, rendered_content.size):
        #                         regexr_fallback = rematch(web_fileInfo[keyz_yr][keyz_day]['station'][j]+str(dates[i, 1]).zfill(3)+r'\d\.'+str(dates[i, 0])[2:]+r'd\.gz', rendered_content[jk]); #Match
        #                         if( regexr_fallback is not None ):
        #                             web_dlSpeed = get_got_webfile(web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j] + '/' + regexr_fallback.group(), web_filePath, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_verbose=FLG_verbose);
        #                             FLG_foundit = True; # Found it
        #                             break; # Ditch at first match
        #                         # END IF
        #                     # END FOR jk
        #                     if( FLG_foundit == True ):
        #                         if( FLG_verbose >= 2 ):
        #                             print('\nINFO in get_cors.py: Requested station `'+web_fileInfo[keyz_yr][keyz_day]['file name'][j]+'` was FOUND at link `'+web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j] + '/' + regexr_fallback.group()+'`.');
        #                         # END IF
        #                     else:
        #                         if( FLG_verbose >= 1 ):
        #                             print('\nWARNING in get_cors.py: Requested station `'+web_fileInfo[keyz_yr][keyz_day]['file name'][j]+'` was not found in FTP folder `'+web_fileSearch + web_fileInfo[keyz_yr][keyz_day]['station'][j] + '/'+'`. Moving on.');
        #                         # END IF
        #                     # END IF
        #                 # END TRYING
        #             # END IF
        #         # END FOR keyz_stat
        #         if( FLG_verbose >= 2 ):
        #             donezo = round(cntr/things2get*colWidth); # Calc how close to being done
        #             sys.stdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
        #             sys.stdout.flush();
        #             cntr += 1; # Increment
        #         # END IF
        #     # END FOR keyz_day
        # # END FOR keyz_yr
        
        # First - build list of links and filenames to expect (avoiding asking CORS too often)
        if( FLG_verbose >= 2 ):
            print('INFO in get_cors.py: Getting station availability information. This takes many webpage requests, so it\'ll take a while, sorry.');
            colWidth = shutil.get_terminal_size().columns//2; # Get half the column width
        # END IF
        if( isinstance(statLocs2use, str) ): # Get all stations (it's set to the string 'all')   
            things2get = dates.shape[0]*len(statLocs_stats); # Get the total things to check out
            var2loop = statLocs_stats; # Loop over all stations
        else: # Otherwise there are specific stations
            things2get = dates.shape[0]*len(statLocs2use); # Get the total things to check out
            var2loop = statLocs2use; # Loop over requested stations
        # END IF
        for i in range(0, dates.shape[0]):         
            web_fileInfo[str(dates[i, 0])][str(dates[i, 1])] = {'file name': [], 'file link': [], 'station': []}; # Prep a day
            web_fileSearch = web_base_site + str(dates[i,0]).zfill(4) + '/' + str(dates[i,1]).zfill(3) + '/'; #build site to go to
            for stats in var2loop:
                # if( web_fileMatches[j][:-1] in statLocs_stats ):
                # Extract the file name
                web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( stats+str(dates[i, 1]).zfill(3)+'0.'+str(dates[i, 0])[2:]+'d.gz' );
                # Build the link
                web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch +stats + '/' + web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][-1] );
                # Record the station
                web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['station'].append( stats );
            # END FOR j
        # END FOR i
        
        
        # Second - Download if not already available
        if( FLG_verbose >= 2 ):
            print('INFO in get_cors.py: Downloading requested files for requested dates and stations.');
        # END IF
        
        parallel_list = []; # Prep
        for keyz_yr in web_fileInfo.keys():
            parallel_list.append([get_cors_download, [web_fileInfo[keyz_yr], keyz_yr, dataDir, web_base_site, FLG_overwrite, web_retryMax, web_retryWait, FLG_verbose, True]])
        # END FOR i
        if( FLG_parallel ):
            results = parallel_orchastrator(parallel_list, parallel_CPUnum=6); # 6 threads to get the files, limit the amount of server-hammering (it throtltes even at 3 threads!)
        else:
            results = []; # Prep
            for i in range(0, len(parallel_list)):
                results.append( parallel_list[i][0](*parallel_list[i][1]) ); # Call the function without parallelization
            # END FOR i
        # END IF
        del results # Didn't actually need to record the results here

        
        # Record file about what was requested
        for keyz in web_fileInfo.keys():
            if( ospath.isdir( ospath.join(dataDir, keyz) ) ): # Only make a metadata file if the year is there
                # Metadata file path
                web_metaFile_path = ospath.join(dataDir, keyz, '.metadata.pkl');
                k = np.where(dates[:, 0] == np.int16(keyz))[0]; # Get which were requested
                if( 'date requested' in web_fileInfo[keyz] ):
                    web_metaFile_reqMissing = np.logical_not( (dates[k, :] == web_fileInfo[keyz]['date requested'][:, None]).all(2).any(0) ); # Get dates missing from web_metaFile_req
                    if( np.any(web_metaFile_reqMissing) ):
                        web_fileInfo[keyz]['date requested'] = np.append(web_fileInfo[keyz]['date requested'], dates[k[web_metaFile_reqMissing], :], axis=0); # Put in the missing dates
                        k = np.argsort(web_fileInfo[keyz]['date requested']@np.array((1,1/12,1/(12*31),))); # Do lazy time sorting
                        web_fileInfo[keyz]['date requested'] = web_fileInfo[keyz]['date requested'][k, :]; # Appy sort
                    # END IF
                else:
                    web_fileInfo[keyz]['date requested'] = dates[k, :].copy();
                # END IF
                
                # Save the metadata as a file in the year folder
                pkl_obj = pkl.dumps(web_fileInfo[keyz], protocol=pkl.HIGHEST_PROTOCOL); # Pickle object
                pkl_smoosh = zlib.compress(pkl_obj); # Compress the object
                with open(web_metaFile_path, 'wb') as pklz:
                    pklz.write(pkl_smoosh); # Write the compressed pickle to disk
                # END WITH
            # END IF
        # END FOR keyz
    else:
        if( FLG_verbose >=2 ):
            print('INFO in get_cors.py: Data has been downloaded already for all dates requested in `'+dataDir+'`.');
        # END IF
    # END IF
# END DEF

def get_cors_statList(statLocs_nameNow, statLocs_webBase, regexr_start, regexr_end, regexr_lat, regexr_long, regexr_X, regexr_Y, regexr_Z, web_retryMax, web_retryWait, FLG_verbose):                
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

def get_cors_download(web_fileInfo_yr, keyz_yr, dataDir, web_base_site, FLG_overwrite, web_retryMax, web_retryWait, FLG_verbose, FLG_parallel=False):
    # if( FLG_parallel == True ):
    #     crunkName = current_process().name.replace('Worker', 'Process'); # Get the process #
    # else:
    #     crunkName = 'Process 0'; # Default to process 0
    # # END IF
    
    web_fileSearch_base = web_base_site + str(keyz_yr).zfill(4) + '/'; # Prep the base
    
    web_dlSpeed = None; # Prep DL speed calcs
    # cntr = 1; # Prep progress bar calcs
    missing = []; # Prep a missing var
    web_fileInfo_day = list(web_fileInfo_yr.keys()); # Get the day keys
    for keyz_day in web_fileInfo_day:
        # Make sure year/day data folder there
        dataDir_now = ospath.join(dataDir, str(keyz_yr), str(keyz_day) ); # Prep it up
        if( ospath.isdir(ospath.join(dataDir_now)) == False ):
            osmakedirs(ospath.join(dataDir_now)); # Make the day folder if it doesn't exist
            filesThere = []; # Nothing goin on
        else:
            filesThere = [filez for filez in oslistdir(dataDir_now) if ospath.isfile(ospath.join(dataDir_now, filez))]; # Get the files in the directory
            for j in range(0, len(filesThere)):
                filesThere_split = ospath.splitext(filesThere[j]); # Split it
                if( filesThere_split[1] !=  '.'+str(keyz_yr)[2:]+'d' ):
                    filesThere[j] = filesThere_split[0]; # Remove the compression extension
                # END IF
            # END FOR j
        # END IF
        web_fileInfo_stats = web_fileInfo_yr[keyz_day]['station']; # Get the station keys
        web_fileSearch = web_fileSearch_base + str(keyz_day).zfill(3) + '/'; #build site to go to
        for j in range(0, len(web_fileInfo_stats)):
            web_filePath = ospath.join(dataDir_now, web_fileInfo_yr[keyz_day]['file name'][j]); # Build the full file path
            
            if( ((ospath.isfile(web_filePath) == False) and (ospath.splitext(web_fileInfo_yr[keyz_day]['file name'][j])[0] not in filesThere)) or (FLG_overwrite == True) ):
                # If the file isn't already there, get it
                try:
                    web_dlSpeed = get_got_webfile(web_fileInfo_yr[keyz_day]['file link'][j], web_filePath, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_test4hdf5=False, FLG_errorOn404=True, FLG_verbose=FLG_verbose);
                except:
                    if( FLG_verbose >= 1 ):
                        print('\nWARNING in get_cors.py: Requested station `'+web_fileInfo_yr[keyz_day]['file name'][j]+'` at link `'+web_fileInfo_yr[keyz_day]['file link'][j]+'` errored, likely not there. Doing a general seach for files matching the format (there\'s a "0" digit that means some level of data, allowing that to vary).');
                    # END IF
                    webpage = get_got_webpg(web_fileSearch + web_fileInfo_yr[keyz_day]['station'][j]+'/', web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage info
                    if( webpage is not None ):
                        FLG_foundit = False; # Prep a flag
                        rendered_content = pd.read_html(StringIO( webpage ))[0][2:-1]['Name'].reset_index(drop=True); # Get the table from the HTML page            
                        for jk in range(0, rendered_content.size):
                            regexr_fallback = rematch(web_fileInfo_yr[keyz_day]['station'][j]+str(keyz_day).zfill(3)+r'\d\.'+str(keyz_yr)[2:]+r'd\.gz', rendered_content[jk]); #Match
                            if( regexr_fallback is not None ):
                                try:
                                    web_dlSpeed = get_got_webfile(web_fileSearch + web_fileInfo_yr[keyz_day]['station'][j] + '/' + regexr_fallback.group(), web_filePath, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_test4hdf5=False, FLG_errorOn404=True, FLG_verbose=FLG_verbose);
                                    FLG_foundit = True; # Found it
                                    break; # Ditch at first match
                                except:
                                    pass; # Try again
                                # END TRYING
                            # END IF
                        # END FOR jk
                        if( FLG_foundit == True ):
                            if( FLG_verbose >= 2 ):
                                print('\nINFO in get_cors.py: Requested station `'+web_fileInfo_yr[keyz_day]['file name'][j]+'` was FOUND at link `'+web_fileSearch + web_fileInfo_yr[keyz_day]['station'][j] + '/' + regexr_fallback.group()+'`.');
                            # END IF
                        else:
                            missing.append( [keyz_day, j] ); # Tack on the missing day and entry #
                            if( FLG_verbose >= 1 ):
                                print('\nWARNING in get_cors.py: Requested station `'+web_fileInfo_yr[keyz_day]['file name'][j]+'` was not found in FTP folder `'+web_fileSearch + web_fileInfo_yr[keyz_day]['station'][j] + '/'+'`. Moving on.');
                            # END IF
                        # END IF
                    # END IF
                # END TRYING
            # END IF
        # END FOR keyz_stat
        # if( FLG_verbose >= 2 ):
        #     donezo = round(cntr/things2get*colWidth); # Calc how close to being done
        #     sys.stdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
        #     sys.stdout.flush();
        #     cntr += 1; # Increment
        # # END IF
    # END FOR keyz_day
    
    return missing
# END DEF

# --- For multiprocessing ---
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