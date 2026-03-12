#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
get_saber(dates, dataDir, madrigalInst, madrigalExp, madrigalLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False);
-  dates = [[2012,12,12],[2012,12,13],[2012,12,15]]; # yr/mo/day; can also use yr/day# format [[2012,333],[2012,334],[2012,336]]
-  dataDir = './path/where/you/want/data' the code will put year folders in this directory for all years requested, and per-year data in each year directory
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
from TAS.get_got import get_got_webpg, get_got_webfile
from os import path as ospath, makedirs as osmakedirs
from shutil import move as shutilmove, get_terminal_size as shutil_get_terminal_size
from tempfile import gettempdir
from sys import stdout as sysstdout
import pandas as pd
from io import StringIO
from re import compile as recompile
import pickle as pkl
import zlib

# Get saber data from its data source automagically
def get_saber(dates, dataDir, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False, FLG_disableHave=False, FLG_verbose=2):
    # Constants
    web_base_site = 'https://cdaweb.gsfc.nasa.gov/pub/data/timed/saber/level2a_netCDF/'; # Website to be used
    
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
            if( 'date requested' in web_metaFile[str(dates[i, 0])] ):
                web_fileInfo[str(dates[i, 0])]['date requested'] = web_metaFile[str(dates[i, 0])]['date requested'].copy(); # Copy over the saved requested dates
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
        
        # First - need to check if data is availiable for the days requested. AIRS goes daynum-by-daynum
        # web_fileMatches = []; # Prep
        for i in range(0, dates.shape[0]):
            # Madrigal lets you browse the FTP data availability by year, so get per-year info that can be applied for later
            web_fileSearch = web_base_site + str(dates[i,0]).zfill(4) + '/' + str(dates[i,1]).zfill(3) + '/'; #build site to go to
            # Get the webpage info
            # rendered_content = get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
            rendered_content = pd.read_html(StringIO( get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0][1:]; # Get the table from the HTML page            
            
            web_fileMatches_now = rendered_content['Name'].to_list(); # Get the matches
            web_fileInfo[str(dates[i, 0])][str(dates[i, 1])] = {'file name': [], 'file link': []}; # Prep a day
            # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'] = []; # Record file name
            # web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'] = []; # Record file link
            # tmp_idx2del = []; # Prep
            for j in range(0, len(web_fileMatches_now)):
                if( not 'not for science use' in web_fileMatches_now[j] ):
                    # Extract the file name
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].append( web_fileMatches_now[j] );
                    # Build the link
                    web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].append( web_fileSearch + web_fileMatches_now[j] );
                # else:
                #     tmp_idx2del.append(j); # Delete
                # END IF
            # END FOR j
            # for j in tmp_idx2del[::-1]:
            #     web_fileMatches_now.pop(j); # Delete
            # # END FOR j
            # web_fileMatches.extend(web_fileMatches_now); # Tack the list onto the gigalist
        
            # Remove multiple versions, choose latest version
            regexr_start = recompile(r'^[a-zA-Z0-9_]+\.'); # Get the start
            regexr_end = recompile(r'[0-9]+\.nc\.gz$'); # Get similarly-named to look for
            idxr2del = []; # Prep a list to delete after the fact due to list science
            idxr2keep = []; # Prep a list to keep to avoid doublework
            for j in range(0, len(web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'])):
                if( (j not in idxr2del) and (j not in idxr2keep) ): # Only do this if it's not already slated for deletion or been checked
                    regexr_now = recompile(r'^'+regexr_start.match(web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][j]).group()+r'[0-9]+\.nc\.gz$'); # Get similarly-named to look for
                    # Get info on the val now
                    matcht_val_tmp = regexr_end.search(web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][j]).group(); # Get the val now and refine it
                    matcht_val = [int(matcht_val_tmp[:matcht_val_tmp.find('.')])]; # Get the number
                    # Prep checking for matches
                    matcht = [web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][j]]; # Records the match name
                    matcht_idx =[j]; # Records the index of the match
                    # Check before j val now
                    bitOstr = web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][:j]; # Roll through bitOstr before j occurs
                    for jj in range(0, len(bitOstr)):
                        tmpval = regexr_now.match(bitOstr[jj]); # Check for matches
                        if( tmpval is not None ):
                            matcht.append(tmpval.group()); # Get the match
                            matcht_idx.append(jj); # Aligned with j
                            matcht_val_tmp = regexr_end.search(tmpval.group()).group(); # Get the val now and refine it
                            matcht_val.append(int(matcht_val_tmp[:matcht_val_tmp.find('.')])); # Get the number
                        # END IF
                    # END FOR jj
                    # Check after j val now
                    bitOstr = web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'][j+1:]; # Roll through bitOstr after j occurs
                    for jj in range(0, len(bitOstr)):
                        tmpval = regexr_now.match(bitOstr[jj]); # Check for matches
                        if( tmpval is not None ):
                            matcht.append(tmpval.group()); # Get the match
                            matcht_idx.append(jj+j+1); # Offset from j
                            matcht_val_tmp = regexr_end.search(tmpval.group()).group(); # Get the val now and refine it
                            matcht_val.append(int(matcht_val_tmp[:matcht_val_tmp.find('.')])); # Get the number
                        # END IF
                    # END FOR jj
                    if( len(matcht_val) > 1 ):
                        # Divine which to keep
                        matcht_val = np.asarray(matcht_val, dtype=int); # Arrayify it
                        matcht_val_where = np.where(matcht_val == matcht_val.max())[0].item(); # Get where the max is
                        idxr2keep.append(matcht_idx[matcht_val_where]); # Choose the one to keep
                        idxr2del.extend(matcht_idx[:matcht_val_where]+matcht_idx[matcht_val_where+1:]); # Get the ones to ditch
                    else:
                        idxr2keep.append(j); # Keep it by default
                    # END IF
                # END IF
            # END FOR j
            for j in idxr2del[::-1]:
                # web_fileMatches.pop(j); # Delete
                web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file name'].pop(j); # Delete
                web_fileInfo[str(dates[i, 0])][str(dates[i, 1])]['file link'].pop(j); # Delete
            # END FOR j
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
                    print('WARNING in get_saber.py: Requested yr/dayNum '+str(dates[i, 0])+'/'+str(dates[i, 1])+' data available.');
                # END IF
            # END IF
        # EDN FOR i
        
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
            print('INFO in get_saber.py: Data has been downloaded already for all dates requested in `'+dataDir+'`.');
        # END IF
    # END IF
# END DEF

# Remove masked arrays and replace with NaN arrays since mask array support is spotty and if supported proper handling isn't always there
def demask_saber(sbrdct):
    for keyz in sbrdct.keys():
        # Undo masking with NaNs; we are free
        if( np.ma.isMaskedArray(sbrdct[keyz]) ):
            if( np.issubdtype(sbrdct[keyz].dtype, np.inexact) and np.any(sbrdct[keyz].mask) ):
                sbrdct[keyz].fill_value = np.nan; # Set the fill to NaN, more useful
                sbrdct[keyz].data[sbrdct[keyz].mask] = np.nan; # Set the data directly to NaNs as well
            # END IF
        # END IF
    # END FOR keyz
    
    return sbrdct
# END DEF

def dist2median(tmpSplit, dist2meanMultiplier=2.5): # dist2meanMultiplier has equivalent function to # of standard deviations
    tmpAvg_median = np.nanmedian(tmpSplit, axis=0);
    # Assemble dist from median, filter out outliers, do a better mean
    tmpAvg_median_dist2median = np.abs(tmpSplit - tmpAvg_median); # Get the distance to the median
    tmpAvg_median_medainOfDist2median = np.median(tmpAvg_median_dist2median); # Bamzo
    k = tmpAvg_median_dist2median >= tmpAvg_median_medainOfDist2median*dist2meanMultiplier; # Data to ditch, more than 2.5 times median distance to median
    tmpSplit_medMasked = tmpSplit.copy(); # Copy this
    tmpSplit_medMasked[k] = np.nan; # More mask
    tmpAvg_meanBetter = np.nanmean(tmpSplit_medMasked, axis=0); # Ensure it's closer
    
    return tmpAvg_meanBetter
# END DEF

# Fix the saber data up, it has issues - the date range given SHOULD have a date range before AND after so cross-day passes can be padded together
def fix_saber(sbrdct, dateRange_now, timeSplitDiffStdev_minimum=np.timedelta64(35, 'm'), FLG_lastDay=False):
    # Do date calcs
    dateVal = str(dateRange_now)[1:-1].replace(' ',''); # Unique date value for the dict
    keyz4date = sorted([int(keyz) for keyz in list(sbrdct[dateVal].keys())]); # Get that sorted up
    keyz4date = [str(keyz) for keyz in keyz4date]; # Convert back to expected strings
    dateVal_B = str(dateRange_now - np.array( (0,1), dtype=np.int16))[1:-1].replace(' ',''); # Unique date value for the dict
    keyz4date_B = list(sbrdct[dateVal_B].keys()); # Get the before keyz
    if 'aligned' in keyz4date_B: keyz4date_B.remove('aligned'); # Remove 'aligned' from the list of keyz if it is there (likely will be)
    keyz4date_B = sorted([int(keyz) for keyz in keyz4date_B]); # Get that sorted up
    keyz4date_B = [str(keyz) for keyz in keyz4date_B]; # Convert back to expected strings
    dateVal_A = str(dateRange_now + np.array( (0,1), dtype=np.int16))[1:-1].replace(' ',''); # Unique date value for the dict
    keyz4date_A = sorted([int(keyz) for keyz in list(sbrdct[dateVal_A].keys())]); # Get that sorted up
    keyz4date_A = [str(keyz) for keyz in keyz4date_A]; # Convert back to expected strings
    
    # Prep
    cntr = 0; # Start a pass counter
    cntr_clr = 0; # Start a color counter
    sbrdct[dateVal]['aligned'] = {}; # Prep an aligned dictinary
    # Convert each SABER observation from its multiple passes across the equator per visit ID into individual passes
    FLG_ditch_end = False; # Prep a flag to ditch the end pass because it's been added into the future pass
    for j in range(len(keyz4date)-1, -1, -1): # Reverse is important for the end ditching logic
        FLG_tack_B = False; # Reset tack_B flag (before day flag)
        FLG_tack_D = False; # Reset tack_D flag (before pass flag, same day)
        FLG_tack_A = False; # Reset tack_A flag (after day flag)
        keyz = keyz4date[j]; # Get the nowkey
        # Unpack to make it easier, see docs https://saber.gats-inc.com/saber_doc/SABER_level2A_NCDUMP_v2.pdf
        ad = sbrdct[dateVal][keyz]['tpAD']; # 0 = ascending, 1 = descending
        dn = sbrdct[dateVal][keyz]['tpDN']; # 0 = day, 1 = night, 2 = twilight
        md = sbrdct[dateVal][keyz]['mode']; # 0 = down, 1 = up
        dateF = sbrdct[dateVal][keyz]['date']; # get the date, in case it goes over the day [this seems independent of timeM]
        timeM = sbrdct[dateVal][keyz]['time']; # millisec, time since midnight [FOR THE FIRST DAY RECORDED, DOES NOT RESET WHEN DAYS CHANGE! ALSO ONLY RELATED TO dateVal NOT dateF WHICH JUST RECORDS THE ACTUAL DAY OF dateVal+timeM would output]
        # timeF = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateF.tolist()[0]]*dateF.size ]; # Build up
        timeF = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal]*timeM.shape[0] ]; # Build up
        timeF = np.tile( np.array([valz[0]+'-01-01' for valz in timeF], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeF], dtype='timedelta64[D]'), (timeM.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
        time = timeF + np.ma.masked_array(np.array(timeM, 'timedelta64[ms]'), mask=timeM.mask); # Build a masked datetime object based off the time and all such things
        timeLTM = sbrdct[dateVal][keyz]['tpSolarLT']; # millisec, local solar time
        timeLTF = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal]*timeLTM.shape[0] ]; # Build up
        timeLTF = np.tile( np.array([valz[0]+'-01-01' for valz in timeLTF], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeLTF], dtype='timedelta64[D]'), (timeLTM.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
        timeLT = np.ma.masked_array(timeLTF, mask=timeLTM.mask) + np.ma.masked_array(np.array(timeLTM, 'timedelta64[ms]'), mask=timeLTM.mask); # Build a masked datetime object based off the time and all such things
        latSC = sbrdct[dateVal][keyz]['sclatitude']; # degc, space craft latitude
        longSC = sbrdct[dateVal][keyz]['sclongitude']; # degc, space craft longitude
        longSC[longSC > 180] -= 360; # Flip to -180 to 180
        altSC = sbrdct[dateVal][keyz]['scaltitude']; # km, space craft altitude
        latTan = sbrdct[dateVal][keyz]['tplatitude']; # degc, tangent point latitude
        longTan = sbrdct[dateVal][keyz]['tplongitude']; # degc, tangent point longitude
        longTan[longTan > 180] -= 360; # Flip to -180 to 180
        altTan = sbrdct[dateVal][keyz]['tpaltitude']; # km, tangent point altitude
        elv = sbrdct[dateVal][keyz]['elevation']; # milliradians, elevation angle
        temp = sbrdct[dateVal][keyz]['ktemp']; # K, kinentic temperature
        
        # Establish breaks between satellite passes
        ad_splits = np.append(np.insert(np.where(np.diff(ad) != 0 )[0] + 1, 0, 0), len(ad)); # Get where the splits are
        if( FLG_ditch_end ):
            ad_splits = ad_splits[:-1]; # Ditch the end b/c it's already been added into the pass after this
        # END IF
                
        # --- BEFORE set --- <- primary mechanism for adding on, either pull before from day before OR pull before from before pass
        if( j == 0 ):
            # Check for data on Before days
            keyz_B = keyz4date_B[-1]; # Get the b4 keyz
            if( abs(int(keyz_B) - int(keyz)) == 1 ): # Ensure orbit #s are consecutive
                # Unpack to make it easier, see docs https://saber.gats-inc.com/saber_doc/SABER_level2A_NCDUMP_v2.pdf
                ad_B = sbrdct[dateVal_B][keyz_B]['tpAD']; # 0 = ascending, 1 = descending
                dn_B = sbrdct[dateVal_B][keyz_B]['tpDN']; # 0 = day, 1 = night, 2 = twilight
                md_B = sbrdct[dateVal_B][keyz_B]['mode']; # 0 = down, 1 = up
                dateF_B = sbrdct[dateVal_B][keyz_B]['date']; # get the date, in case it goes over the day [this seems independent of timeM]
                timeM_B = sbrdct[dateVal_B][keyz_B]['time']; # millisec, time since midnight [FOR THE FIRST DAY RECORDED, DOES NOT RESET WHEN DAYS CHANGE! ALSO ONLY RELATED TO dateVal NOT dateF WHICH JUST RECORDS THE ACTUAL DAY OF dateVal+timeM would output]
                # timeF_B = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateF_B.tolist()[0]]*dateF_B.size ]; # Build up
                timeF_B = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal_B]*timeM_B.shape[0] ]; # Build up
                timeF_B = np.tile( np.array([valz[0]+'-01-01' for valz in timeF_B], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeF_B], dtype='timedelta64[D]'), (timeM_B.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                time_B = timeF_B + np.ma.masked_array(np.array(timeM_B, 'timedelta64[ms]'), mask=timeM_B.mask); # Build a masked datetime object based off the time and all such things
                timeLTM_B = sbrdct[dateVal_B][keyz_B]['tpSolarLT']; # millisec, local solar time
                timeLTF_B = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal_B]*timeLTM_B.shape[0] ]; # Build up
                timeLTF_B = np.tile( np.array([valz[0]+'-01-01' for valz in timeLTF_B], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeLTF_B], dtype='timedelta64[D]'), (timeLTM_B.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                timeLT_B = np.ma.masked_array(timeLTF_B, mask=timeLTM_B.mask) + np.ma.masked_array(np.array(timeLTM_B, 'timedelta64[ms]'), mask=timeLTM_B.mask); # Build a masked datetime object based off the time and all such things
                latSC_B = sbrdct[dateVal_B][keyz_B]['sclatitude']; # degc, space craft latitude
                longSC_B = sbrdct[dateVal_B][keyz_B]['sclongitude']; # degc, space craft longitude
                longSC_B[longSC_B > 180] -= 360; # Flip to -180 to 180
                altSC_B = sbrdct[dateVal_B][keyz_B]['scaltitude']; # km, space craft altitude
                latTan_B = sbrdct[dateVal_B][keyz_B]['tplatitude']; # degc, tangent point latitude
                longTan_B = sbrdct[dateVal_B][keyz_B]['tplongitude']; # degc, tangent point longitude
                longTan_B[longTan_B > 180] -= 360; # Flip to -180 to 180
                altTan_B = sbrdct[dateVal_B][keyz_B]['tpaltitude']; # km, tangent point altitude
                elv_B = sbrdct[dateVal_B][keyz_B]['elevation']; # milliradians, elevation angle
                temp_B = sbrdct[dateVal_B][keyz_B]['ktemp']; # K, kinentic temperature
                
                # Establish breaks between satellite passes
                ad_splits_B = np.append(np.insert(np.where(np.diff(ad_B) !=0 )[0] + 1, 0, 0), len(ad_B)); # Get where the splits are
                
                timeSplit_B = time_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(); # Get it out
                
                timeSplitDiff_B = np.diff(timeSplit_B[:, 0].data); # Diff the time
                timeSplitDiffMedian_B = np.median(timeSplitDiff_B); # Get the median time diff
                timeSplitDiffStdev_B = np.std(timeSplitDiff_B.astype(int)).astype('timedelta64[ms]'); # Get the standard deviation of the time diffs
                if( timeSplitDiffStdev_B*3 > timeSplitDiffStdev_minimum ):
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_B*3; # Use this one
                else:
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_minimum; # Use this one
                # END IF
                if( np.abs((timeSplit_B[-1, 0] + timeSplitDiffMedian_B) - time[ad_splits[0]:ad_splits[1]][0, 0]) <= timeSplitDiffStdev_2use ):
                    # Tack it on if we're good to go!
                    FLG_tack_B = True; # Tack later
                else:
                    breakpoint()
                # END IF
            else:
                print('WARNING in gw.py: Missing orbits detected: Dateval '+dateVal+' orbit# '+keyz+' has previous orbit# [on previous dateval '+dateVal_B+'] as '+keyz_B+' which is NOT the previous orbit.');
            # END IF
        else:
            # Otherwise, we're going to move a pass from the previous pass to our current pass
            
            # Check for data on After days
            keyz_D = keyz4date[j-1]; # Get the af keyz
            if( abs(int(keyz_D) - int(keyz)) == 1 ): # Ensure orbit #s are consecutive
                # Unpack to make it easier, see docs https://saber.gats-inc.com/saber_doc/SABER_level2A_NCDUMP_v2.pdf
                ad_D = sbrdct[dateVal][keyz_D]['tpAD']; # 0 = ascending, 1 = descending
                dn_D = sbrdct[dateVal][keyz_D]['tpDN']; # 0 = day, 1 = night, 2 = twilight
                md_D = sbrdct[dateVal][keyz_D]['mode']; # 0 = down, 1 = up
                dateF_D = sbrdct[dateVal][keyz_D]['date']; # get the date, in case it goes over the day [this seems independent of timeM]
                timeM_D = sbrdct[dateVal][keyz_D]['time']; # millisec, time since midnight [FOR THE FIRST DAY RECORDED, DOES NOT RESET WHEN DAYS CHANGE! ALSO ONLY RELATED TO dateVal NOT dateF WHICH JUST RECORDS THE ACTUAL DAY OF dateVal+timeM would output]
                # timeF_D = [[str(valz)[:4], int(str(valz)[4:])]  for valz in [dateF_D.tolist()[0]]*dateF_D.size ]; # Build up
                timeF_D = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal]*timeM_D.shape[0] ]; # Build up
                timeF_D = np.tile( np.array([valz[0]+'-01-01' for valz in timeF_D], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeF_D], dtype='timedelta64[D]'), (timeM_D.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                time_D = timeF_D + np.ma.masked_array(np.array(timeM_D, 'timedelta64[ms]'), mask=timeM_D.mask); # Build a masked datetime object based off the time and all such things
                timeLTM_D = sbrdct[dateVal][keyz_D]['tpSolarLT']; # millisec, local solar time
                timeLTF_D = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal]*timeLTM_D.shape[0] ]; # Build up
                timeLTF_D = np.tile( np.array([valz[0]+'-01-01' for valz in timeLTF_D], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeLTF_D], dtype='timedelta64[D]'), (timeLTM_D.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                timeLT_D = np.ma.masked_array(timeLTF_D, mask=timeLTM_D.mask) + np.ma.masked_array(np.array(timeLTM_D, 'timedelta64[ms]'), mask=timeLTM_D.mask); # Build a masked datetime object based off the time and all such things
                latSC_D = sbrdct[dateVal][keyz_D]['sclatitude']; # degc, space craft latitude
                longSC_D = sbrdct[dateVal][keyz_D]['sclongitude']; # degc, space craft longitude
                longSC_D[longSC_D > 180] -= 360; # Flip to -180 to 180
                altSC_D = sbrdct[dateVal][keyz_D]['scaltitude']; # km, space craft altitude
                latTan_D = sbrdct[dateVal][keyz_D]['tplatitude']; # degc, tangent point latitude
                longTan_D = sbrdct[dateVal][keyz_D]['tplongitude']; # degc, tangent point longitude
                longTan_D[longTan_D > 180] -= 360; # Flip to -180 to 180
                altTan_D = sbrdct[dateVal][keyz_D]['tpaltitude']; # km, tangent point altitude
                elv_D = sbrdct[dateVal][keyz_D]['elevation']; # milliradians, elevation angle
                temp_D = sbrdct[dateVal][keyz_D]['ktemp']; # K, kinentic temperature
                
                # Establish breaks between satellite passes
                ad_splits_D = np.append(np.insert(np.where(np.diff(ad_D) !=0 )[0] + 1, 0, 0), len(ad_D)); # Get where the splits are
                
                timeSplit_D = time_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(); # Get it out
                
                timeSplitDiff_D = np.diff(timeSplit_D[:, 0].data); # Diff the time
                timeSplitDiffMedian_D = np.median(timeSplitDiff_D); # Get the median time diff
                timeSplitDiffStdev_D = np.std(timeSplitDiff_D.astype(int)).astype('timedelta64[ms]'); # Get the standard deviation of the time diffs
                if( timeSplitDiffStdev_D*3 > timeSplitDiffStdev_minimum ):
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_D*3; # Use this one
                else:
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_minimum; # Use this one
                # END IF
                if( np.abs((timeSplit_D[-1, 0] + timeSplitDiffMedian_D) - time[ad_splits[0]:ad_splits[1]][0, 0]) <= timeSplitDiffStdev_2use ):
                    # Tack it on if we're good to go!
                    FLG_tack_D = True; # Tack later
                    FLG_ditch_end = True; # Ditch that end pass next time around, we've added it to the pass now
                else:
                    breakpoint() # Occurances seem legit missing data
                # END IF
            else:
                print('WARNING in gw.py: Missing orbits detected: Dateval '+dateVal+' orbit# '+keyz+' has previous orbit# as '+keyz_D+' which is NOT the previous orbit.');
            # END IF
        # END IF
        
        # --- AFTER set --- <--- only for end days OR removing ending passes that are caught by the following day yoinking them into it
        if( (j == (len(keyz4date) - 1)) and FLG_lastDay ): # Only look after on the last date, other days are going to pull before data
            # Check for data on After days
            keyz_A = keyz4date_A[0]; # Get the af keyz
            if( abs(int(keyz_A) - int(keyz)) == 1 ): # Ensure orbit #s are consecutive
                # Unpack to make it easier, see docs https://saber.gats-inc.com/saber_doc/SABER_level2A_NCDUMP_v2.pdf
                ad_A = sbrdct[dateVal_A][keyz_A]['tpAD']; # 0 = ascending, 1 = descending
                dn_A = sbrdct[dateVal_A][keyz_A]['tpDN']; # 0 = day, 1 = night, 2 = twilight
                md_A = sbrdct[dateVal_A][keyz_A]['mode']; # 0 = down, 1 = up
                dateF_A = sbrdct[dateVal_A][keyz_A]['date']; # get the date, in case it goes over the day [this seems independent of timeM]
                timeM_A = sbrdct[dateVal_A][keyz_A]['time']; # millisec, time since midnight [FOR THE FIRST DAY RECORDED, DOES NOT RESET WHEN DAYS CHANGE! ALSO ONLY RELATED TO dateVal NOT dateF WHICH JUST RECORDS THE ACTUAL DAY OF dateVal+timeM would output]
                # timeF_A = [[str(valz)[:4], int(str(valz)[4:])]  for valz in [dateF_A.tolist()[0]]*dateF_A.size ]; # Build up
                timeF_A = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal_A]*timeM_A.shape[0] ]; # Build up
                timeF_A = np.tile( np.array([valz[0]+'-01-01' for valz in timeF_A], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeF_A], dtype='timedelta64[D]'), (timeM_A.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                time_A = timeF_A + np.ma.masked_array(np.array(timeM_A, 'timedelta64[ms]'), mask=timeM_A.mask); # Build a masked datetime object based off the time and all such things
                timeLTM_A = sbrdct[dateVal_A][keyz_A]['tpSolarLT']; # millisec, local solar time
                timeLTF_A = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal_A]*timeLTM_A.shape[0] ]; # Build up
                timeLTF_A = np.tile( np.array([valz[0]+'-01-01' for valz in timeLTF_A], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeLTF_A], dtype='timedelta64[D]'), (timeLTM_A.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                timeLT_A = np.ma.masked_array(timeLTF_A, mask=timeLTM_A.mask) + np.ma.masked_array(np.array(timeLTM_A, 'timedelta64[ms]'), mask=timeLTM_A.mask); # Build a masked datetime object based off the time and all such things
                latSC_A = sbrdct[dateVal_A][keyz_A]['sclatitude']; # degc, space craft latitude
                longSC_A = sbrdct[dateVal_A][keyz_A]['sclongitude']; # degc, space craft longitude
                longSC_A[longSC_A > 180] -= 360; # Flip to -180 to 180
                altSC_A = sbrdct[dateVal_A][keyz_A]['scaltitude']; # km, space craft altitude
                latTan_A = sbrdct[dateVal_A][keyz_A]['tplatitude']; # degc, tangent point latitude
                longTan_A = sbrdct[dateVal_A][keyz_A]['tplongitude']; # degc, tangent point longitude
                longTan_A[longTan_A > 180] -= 360; # Flip to -180 to 180
                altTan_A = sbrdct[dateVal_A][keyz_A]['tpaltitude']; # km, tangent point altitude
                elv_A = sbrdct[dateVal_A][keyz_A]['elevation']; # milliradians, elevation angle
                temp_A = sbrdct[dateVal_A][keyz_A]['ktemp']; # K, kinentic temperature
                
                # Establish breaks between satellite passes
                ad_splits_A = np.append(np.insert(np.where(np.diff(ad_A) !=0 )[0] + 1, 0, 0), len(ad_A)); # Get where the splits are
                
                timeSplit_A = time_A[ad_splits_A[0]:ad_splits_A[1]].copy(); # Get it out
                
                timeSplitDiff_A = np.diff(timeSplit_A[:, 0].data); # Diff the time
                timeSplitDiffMedian_A = np.median(timeSplitDiff_A); # Get the median time diff
                timeSplitDiffStdev_A = np.std(timeSplitDiff_A.astype(int)).astype('timedelta64[ms]'); # Get the standard deviation of the time diffs
                if( timeSplitDiffStdev_A*3 > timeSplitDiffStdev_minimum ):
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_A*3; # Use this one
                else:
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_minimum; # Use this one
                # END IF
                if( np.abs((timeSplit_A[0, 0] - timeSplitDiffMedian_A) - time[ad_splits[-2]:ad_splits[-1]][-1, 0]) <= timeSplitDiffStdev_2use ):
                    # Tack it on if we're good to go!
                    FLG_tack_A = True; # Tack later
                else:
                    breakpoint()
                # END IF
            else:
                print('WARNING in gw.py: Missing orbits detected: Dateval '+dateVal+' orbit# '+keyz+' has previous orbit# [on following dateval '+dateVal_A+'] as '+keyz_A+' which is NOT the following orbit.');
            # END IF
        elif( (j == (len(keyz4date) - 1)) and (not FLG_lastDay) ): # Case for removing the tail end pass if the next day will catch it
            keyz_Z = keyz4date_A[0]; # Get the af keyz
            if( abs(int(keyz_Z) - int(keyz)) == 1 ): # Ensure orbit #s are consecutive
                # Unpack to make it easier, see docs https://saber.gats-inc.com/saber_doc/SABER_level2A_NCDUMP_v2.pdf
                ad_Z = sbrdct[dateVal_A][keyz_Z]['tpAD']; # 0 = ascending, 1 = descending
                timeM_Z = sbrdct[dateVal_A][keyz_Z]['time']; # millisec, time since midnight [FOR THE FIRST DAY RECORDED, DOES NOT RESET WHEN DAYS CHANGE!]
                timeF_Z = [[str(valz)[:4], int(str(valz)[4:])] for valz in [dateVal_A]*timeM_Z.shape[0] ]; # Build up
                timeF_Z = np.tile( np.array([valz[0]+'-01-01' for valz in timeF_Z], dtype='datetime64') + np.array([(valz[1] - 1) for valz in timeF_Z], dtype='timedelta64[D]'), (timeM_Z.shape[1], 1)).T; # Rock n roll, get a datetime object and tile it out
                time_Z = timeF_Z + np.ma.masked_array(np.array(timeM_Z, 'timedelta64[ms]'), mask=timeM_Z.mask); # Build a masked datetime object based off the time and all such things
                
                # Establish breaks between satellite passes
                ad_splits_Z = np.append(np.insert(np.where(np.diff(ad_Z) !=0 )[0] + 1, 0, 0), len(ad_Z)); # Get where the splits are
                
                timeSplit_Z = time_Z[ad_splits_Z[0]:ad_splits_Z[1]].copy(); # Get it out
                
                timeSplitDiff_Z = np.diff(timeSplit_Z[:, 0].data); # Diff the time
                timeSplitDiffMedian_Z = np.median(timeSplitDiff_Z); # Get the median time diff
                timeSplitDiffStdev_Z = np.std(timeSplitDiff_Z.astype(int)).astype('timedelta64[ms]'); # Get the standard deviation of the time diffs
                if( timeSplitDiffStdev_Z*3 > timeSplitDiffStdev_minimum ):
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_Z*3; # Use this one
                else:
                    timeSplitDiffStdev_2use = timeSplitDiffStdev_minimum; # Use this one
                # END IF
                if( np.abs((timeSplit_Z[0, 0] - timeSplitDiffMedian_Z) - time[ad_splits[-2]:ad_splits[-1]][-1, 0]) <= timeSplitDiffStdev_2use):
                    # In this case ditch the ending split because it'll be caught in the before check on the next day
                    ad_splits = ad_splits[:-1]; # Ditch
                else:
                    breakpoint()
                # END IF
            else:
                print('WARNING in gw.py: Missing orbits detected: Dateval '+dateVal+' orbit# '+keyz+' has following orbit# as '+keyz_A+' which is NOT the following orbit.');
            # END IF
        # END IF
        
        for jj in range(0, ad_splits.size-1):
            sbrdct[dateVal]['aligned'][cntr] = {'pass#':[keyz]}; # Prep a dict for this pass and record the pass#
            
            # Get a pass out
            passz = {
                'ad': ad[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'dn': dn[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'md': md[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'dateF': dateF[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'timeM': timeM[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'timeF': timeF[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'time': time[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'timeLT': timeLT[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'latSC': latSC[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'longSC': longSC[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'altSC': altSC[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'latTan': latTan[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'longTan': longTan[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'altTan': altTan[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'elv': elv[ad_splits[jj]:ad_splits[jj+1]].copy(),
                'temp': temp[ad_splits[jj]:ad_splits[jj+1]].copy(),
                }; # Prep a pass to get
            
            # Tack on as needed
            if( (FLG_tack_B == True) and (jj == 0) ):
                passz_B = {
                    'ad': ad_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'dn': dn_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'md': md_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'dateF': dateF_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'timeM': timeM_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'timeF': timeF_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'time': time_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'timeLT': timeLT_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'latSC': latSC_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'longSC': longSC_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'altSC': altSC_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'latTan': latTan_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'longTan': longTan_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'altTan': altTan_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'elv': elv_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    'temp': temp_B[ad_splits_B[-2]:ad_splits_B[-1]].copy(),
                    }; # Prep a pass to get
                for keyz_now in passz.keys(): # Based off of passz so that an update to passz will trigger a missing key error in passz_B if it is also not updated
                    if( np.ma.is_masked(passz[keyz_now]) ):
                        passz[keyz_now] = np.ma.concatenate((passz_B[keyz_now], passz[keyz_now]), axis=0); # Tack it on at the front
                    else:
                        passz[keyz_now] = np.append(passz_B[keyz_now], passz[keyz_now], axis=0); # Tack it on at the front
                    # END IF
                # END FOR keyz_now
                sbrdct[dateVal]['aligned'][cntr]['pass#'] = [keyz_B, keyz]; # Tack on the keyz involved
                FLG_tack_B = False; # Disable after tacking
            elif( (FLG_tack_A == True) and (jj == ad_splits.size-2) ):
                passz_A = {
                    'ad': ad_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'dn': dn_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'md': md_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'dateF': dateF_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'timeM': timeM_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'timeF': timeF_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'time': time_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'timeLT': timeLT_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'latSC': latSC_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'longSC': longSC_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'altSC': altSC_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'latTan': latTan_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'longTan': longTan_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'altTan': altTan_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'elv': elv_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    'temp': temp_A[ad_splits_A[-2]:ad_splits_A[-1]].copy(),
                    }; # Prep a pass to get
                for keyz_now in passz.keys(): # Based off of passz so that an update to passz will trigger a missing key error in passz_A if it is also not updated
                    if( np.ma.is_masked(passz[keyz_now]) ):
                        passz[keyz_now] = np.ma.concatenate((passz[keyz_now], passz_A[keyz_now]), axis=0); # Tack it on at the front
                    else:
                        passz[keyz_now] = np.append(passz[keyz_now], passz_A[keyz_now], axis=0); # Tack it on at the end
                    # END IF
                # END FOR keyz_now
                sbrdct[dateVal]['aligned'][cntr]['pass#'] = [keyz, keyz_A]; # Tack on the keyz involved
                FLG_tack_A = False; # Disable after tacking
            elif( FLG_tack_D == True ):
                # Tack the previous pass in the same day otherwise, if it is on
                passz_D = {
                    'ad': ad_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'dn': dn_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'md': md_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'dateF': dateF_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'timeM': timeM_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'timeF': timeF_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'time': time_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'timeLT': timeLT_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'latSC': latSC_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'longSC': longSC_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'altSC': altSC_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'latTan': latTan_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'longTan': longTan_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'altTan': altTan_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'elv': elv_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    'temp': temp_D[ad_splits_D[-2]:ad_splits_D[-1]].copy(),
                    }; # Prep a pass to get
                for keyz_now in passz.keys(): # Based off of passz so that an update to passz will trigger a missing key error in passz_D if it is also not updated
                    if( np.ma.is_masked(passz[keyz_now]) ):
                        passz[keyz_now] = np.ma.concatenate((passz_D[keyz_now], passz[keyz_now]), axis=0); # Tack it on at the front
                    else:
                        passz[keyz_now] = np.append(passz_D[keyz_now], passz[keyz_now], axis=0); # Tack it on at the front
                    # END IF
                # END FOR keyz_now
                sbrdct[dateVal]['aligned'][cntr]['pass#'] = [keyz_D, keyz]; # Tack on the keyz involved
                FLG_tack_D = False; # Disable after tacking
            # END IF
            
            if( np.all(passz['ad'] == 0) ):
                passz['direction'] = 'ascending'; # Record ascending
            elif( np.all(passz['ad'] == 1) ):
                passz['direction'] = 'descending'; # Record descending
            else:
                print('ERROR in gw: pass is not ONLY ascending or descending. Investigate.');
                breakpoint()
                pass
            # END IF
            
            diffcheck = np.abs(np.diff(passz['longSC'].data, axis=1)); # Get the diff check
            diffcheck_median = np.nanmedian(diffcheck); # Get the median time diff
            diffcheck_diff2median = np.abs(diffcheck - diffcheck_median); # Get the distance to the median
            diffcheck_diff2median_median = np.nanmedian(diffcheck_diff2median); # Get the median distance to the median
            k = diffcheck > diffcheck_median + diffcheck_diff2median_median*1000; # Prep what to get
            k[diffcheck > 180] = False; # Set flips around the globe to false [opens a gap where if the weird thing starts on a lat/long flip, ignoring possibility since SC doesn't matter rn]
            if( np.any( k ) ):
                if( k.sum() > 1 ): # This targets where the longitude gets spread across the globe for some reason but the data seems fine (incl. lat)
                    k = np.where(k); # Get where
                    k = [[k[0][0], k[1][1]], [k[0][-1], k[1][-1]+1]]; # Select the start and end
                    if( np.abs( passz['longSC'][k[0][0], k[0][1]-1] - passz['longSC'][k[0][0], k[1][1]]) > diffcheck_median + diffcheck_diff2median_median*100 ):
                        # Catch if it drops into masked NaNs so diff doesn't catch it ends with a big diff -> but also can't go end-to-end
                        steps2use = np.nanmedian(np.diff(passz['longSC'][k[0][0], k[0][1]-1-(k[1][1]-k[0][1]):k[0][1]])); # Get the steps to use
                        passz['longSC'][k[0][0], k[0][1]:k[1][1]+1] = passz['longSC'][k[0][0], k[0][1]-1] + steps2use*np.arange(1, k[1][1]-k[0][1]+2); # Just smooth it over
                    else:
                        passz['longSC'][k[0][0], k[0][1]:k[1][1]] = np.linspace(passz['longSC'][k[0][0], k[0][1]-1], passz['longSC'][k[0][0], k[1][1]], num=k[1][1]-k[0][1]); # Just smooth it over
                    # END IF
                    # passz['latSC'][k[0][0], k[0][1]:k[1][1]] = np.linspace(passz['latSC'][k[0][0], k[0][1]-1], passz['latSC'][k[0][0], k[1][1]], num=k[1][1]-k[0][1]); # Just smooth it over [seems to be OK actually]
                else:
                    # This targets one off weird values (found one that looks like it ought to be divided by 100 but wasn't)
                    k = np.where(k); # Get where
                    if( k[1][0] < passz['longSC'].shape[1]//2 ): # Lazy way to get the steps2use by just saying if it's in the first 1/2 go to end of array, if in last 1/2 go to start of array
                        steps2use = np.nanmedian(np.diff(passz['longSC'][k[0][0], k[1][0]+1:])); # Get the steps to use
                        passz['longSC'][k[0][0], k[1][0]] = passz['longSC'][k[0][0], k[1][0]+1] - steps2use; # Take the value to the right and adjust to be the problematic one
                    else:
                        steps2use = np.nanmedian(np.diff(passz['longSC'][k[0][0], :k[1][0]])); # Get the steps to use
                        passz['longSC'][k[0][0], k[1][0]] = passz['longSC'][k[0][0], k[1][0] - 1] + steps2use; # Take the value to the left and adjust to be the problematic one
                    # END IF
                # END IF
            # END IF
            
            # Record the pass
            sbrdct[dateVal]['aligned'][cntr].update(passz); # Throw this into the dict
            
            cntr +=1; # Move to the next record
        # END FOR jj
    # END FOR j
    
    return sbrdct
# END DEF

def bin_saber(data_saber, dateRange, avgr_long_width=7.5, avgr_long_latLim=[-90, 90]):
    # --- Assemble Obs into Longitude Bins as Described by avgr_long_width ---
    avgr_long_timeBins = np.append(np.arange(-180, 180, step=avgr_long_width), 180); # Assemble time bins to fit data into
    avgr_long_timeBins_binned = {str(avgr_long_timeBins[i]): {'all':[]} for i in range(0, avgr_long_timeBins.size - 1)}; # Get the bins set up
    
    for dd in range(0, dateRange.shape[0]):
        dateVal = str(dateRange[dd, :])[1:-1].replace(' ',''); # Unique date value for the dict
        
        for cntr in data_saber[dateVal]['aligned'].keys():
            k = np.where(np.min(np.abs(data_saber[dateVal]['aligned'][cntr]['latTan'])) == np.abs(data_saber[dateVal]['aligned'][cntr]['latTan'])); # Get the index at the equator
            ktay = data_saber[dateVal]['aligned'][cntr]['longTan'][k[0], k[1]]; # Get the time at the equator
            # ktay = datetime.fromtimestamp(ktay.data.astype(int).item() / 1e3, UTC); # Convert to datetime object that can actually extract hour easily why is datetime64 unfinished 
            # ktay = ktay.hour + ktay.minute/60 + ktay.second/3600; # Get the hour in decimal
            # ktg = data_saber[dateVal]['aligned'][cntr]['time'][k[0], k[1]].item();
            ktay = data_saber[dateVal]['aligned'][cntr]['longTan'][k[0], k[1]][0]; # Get the longitude
            # if ktay < 0: ktay += 360; # Yeet [just use -180 to 180 directly]
            # ktay = ktay/15; # Get the hour in decimal from longitude directly not the LT [just use longitude directly]
            k = np.searchsorted( avgr_long_timeBins, ktay, side='left' ) - 1; # Get which bin it goes in
            avgr_long_timeBins_binned[str(avgr_long_timeBins[k])]['all'].append([dateVal, cntr]); # Bin it up
        # END FOR cntr
    # END FOR dd
    
    # --- Split by local day/night/twilight ---
    for keyz in avgr_long_timeBins_binned:
        avgr_long_timeBins_binned[keyz]['localtime day'] = []; # Prep
        avgr_long_timeBins_binned[keyz]['localtime night'] = []; # Prep
        avgr_long_timeBins_binned[keyz]['localtime twilight2day'] = []; # Prep
        avgr_long_timeBins_binned[keyz]['localtime twilight2night'] = []; # Prep
        for passz in avgr_long_timeBins_binned[keyz]['all']:
            avgr_long_timeBins_binned[keyz]['localtime day'].append( data_saber[passz[0]]['aligned'][passz[1]]['dn'] == 0 ); # Get the day mansk
            avgr_long_timeBins_binned[keyz]['localtime day'][-1][data_saber[passz[0]]['aligned'][passz[1]]['dn'].mask] = False; # Mask it out
            avgr_long_timeBins_binned[keyz]['localtime day'][-1] = avgr_long_timeBins_binned[keyz]['localtime day'][-1].data; # Remove the mask
            avgr_long_timeBins_binned[keyz]['localtime night'].append( data_saber[passz[0]]['aligned'][passz[1]]['dn'] == 1 ); # Get the night mask
            avgr_long_timeBins_binned[keyz]['localtime night'][-1][data_saber[passz[0]]['aligned'][passz[1]]['dn'].mask] = False; # Mask it out
            avgr_long_timeBins_binned[keyz]['localtime night'][-1] = avgr_long_timeBins_binned[keyz]['localtime night'][-1].data; # Remove the mask
            k = data_saber[passz[0]]['aligned'][passz[1]]['dn'] == 2; # Get if any 2's
            k[k.mask] = False; # Apply the mask if it is there
            k = k.data; # Only data no more masked arrays
            if( k.any() ):
                ks = np.append(np.insert(np.where(np.abs(k[1:] - k[:-1].astype(int)))[0], 0, 0), k.size); # Look for non-consecutive twilight times
                if( (k.sum() > 1) and (ks.size > 3) ):
                    k = [k[ks[jk]:ks[jk+1]] for jk in range(0, len(ks)-1)]; # Get the different lists [this should have days after and before so it can figure out what kind of twilight this is]
                    tmp_dn = [data_saber[passz[0]]['aligned'][passz[1]]['dn'][ks[jk]:ks[jk+1]] for jk in range(0, len(ks)-1)] + [data_saber[passz[0]]['aligned'][passz[1]]['dn'][ks[-1]:]]; # Split up dn so indexes are relevant
                    for jk in range(len(k)-1, -1, -1): # This loop is to fix my lazyness by removing 0-sized arrays
                        if( k[jk].size == 0 ):
                            k.pop(jk); # Ditch
                            tmp_dn.pop(jk); # Ditch
                        # END IF
                    # END FOR jk
                    for jk in range(len(k)-1, 0, -1): # This loop is to fix my lazyness wrt splitting twilight logic
                        if( (k[jk][0] == True) and (k[jk-1][-1] == True) ):
                            k[jk-1] = np.hstack( (k[jk-1], k[jk]) ); # Tack together
                            tmp_dn[jk-1] = np.hstack( (tmp_dn[jk-1], tmp_dn[jk]) ); # Tack together
                            k.pop(jk); # Ditch
                            tmp_dn.pop(jk); # Ditch
                        elif( np.all(k[jk] == False) or np.all(k[jk-1] == False) ):
                            k[jk-1] = np.hstack( (k[jk-1], k[jk]) ); # Tack together
                            tmp_dn[jk-1] = np.hstack( (tmp_dn[jk-1], tmp_dn[jk]) ); # Tack together
                            k.pop(jk); # Ditch
                            tmp_dn.pop(jk); # Ditch
                        # END IF
                    # END FOR jk
                else:
                    k = [k]; # List it so loop is happy
                    tmp_dn = [data_saber[passz[0]]['aligned'][passz[1]]['dn']]; # List it so loop is happy
                # END IF
                
                # I want to make it clear I am just bandaiding ontop of bandaids here
                tmp_twilight = {'localtime twilight2day': [], 
                                'localtime twilight2night': [] }; # Prep a multi-loop holder for into
                for jk in range(0, len(k)):
                    k_end = np.where(k[jk])[0][-1] + 1; # Get the end of the band
                    if( k_end >= k[jk].size ):
                        k_start = np.where(k[jk])[0][0] - 1; # Get the start of the band
                        if( k_start >= 0 ):
                            if( tmp_dn[jk][k_start] == 0 ): # Day to start, goes to night
                                lt_type = 'localtime twilight2night'; # Good to go
                                lt_typeNOT = 'localtime twilight2day'; # Good to go
                            elif( tmp_dn[jk][k_start] == 1 ): # Night to start, goes to day
                                lt_type = 'localtime twilight2day'; # Good to go
                                lt_typeNOT = 'localtime twilight2night'; # Good to go
                            else:
                                breakpoint()
                            # END IF
                        else:
                            kb = np.where(tmp_dn[jk][np.where(k[jk])[0][0]:] != 2)[0][0]; # Get first non-twilight (2) value
                            if( tmp_dn[jk][kb] == 1 ): # Twilight to start, goes to night
                                lt_type = 'localtime twilight2night'; # Good to go
                                lt_typeNOT = 'localtime twilight2day'; # Good to go
                            elif( tmp_dn[jk][kb] == 0 ): # Twilight to start, goes to day
                                lt_type = 'localtime twilight2day'; # Good to go
                                lt_typeNOT = 'localtime twilight2night'; # Good to go
                            else:
                                breakpoint()
                            # END IF
                        # END IF
                    else:
                        if( tmp_dn[jk][k_end] == 0 ): # Day to end, goes to day
                            lt_type = 'localtime twilight2day'; # Good to go
                            lt_typeNOT = 'localtime twilight2night'; # Good to go
                        elif( tmp_dn[jk][k_end] == 1 ): # Night to end, goes to night
                            lt_type = 'localtime twilight2night'; # Good to go
                            lt_typeNOT = 'localtime twilight2day'; # Good to go
                        else:
                            breakpoint()
                        # END IF
                    # END IF
                    tmp_twilight[lt_type].append( k[jk] ); # Get the twilight mask
                    tmp_twilight[lt_typeNOT].append( np.zeros(k[jk].size, dtype=np.bool_) ); # Get an all-false mask for the unused one to keep count right
                # END FOR jk
                avgr_long_timeBins_binned[keyz]['localtime twilight2day'].append( np.concatenate(tmp_twilight['localtime twilight2day']) ); # Set the twilight2day mask
                avgr_long_timeBins_binned[keyz]['localtime twilight2night'].append( np.concatenate(tmp_twilight['localtime twilight2night']) ); # Set the twilight2night mask
                if( avgr_long_timeBins_binned[keyz]['localtime twilight2night'][-1].size != data_saber[passz[0]]['aligned'][passz[1]]['dn'].size ):
                    breakpoint()
            else:
                avgr_long_timeBins_binned[keyz]['localtime twilight2day'].append( k ); # Get an all-false mask for the unused one to keep count right
                avgr_long_timeBins_binned[keyz]['localtime twilight2night'].append( k ); # Get an all-false mask for the unused one to keep count right
            # END IF      
        # END FOR passz
    # END FOR keyz
    
    # --- Build Averages for Bins and Local Times ---
    avgr_avgz = {str(avgr_long_timeBins[i]): {} for i in range(0, avgr_long_timeBins.size - 1)}; # Get the bins set up
    avgr_avgz['time bins'] = avgr_long_timeBins[:-1]; # Record this as well
    avgr_avgz['time bins keys']  = [str(i) for i in avgr_long_timeBins[:-1]]; # Record the string for the keys
    avgr_avgz['time bins ranges'] = {str(avgr_long_timeBins[i]): str(avgr_long_timeBins[i])+' to '+str(avgr_long_timeBins[i+1]) for i in range(0, len(avgr_long_timeBins[:-1]))}; # Make a dict for the ranges
    for keyz in avgr_avgz['time bins keys']:
            l4d_temp = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold day data
            l4d_alt = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold day data
            l4n_temp = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold night data
            l4n_alt = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold night data
            l4t2d_temp = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold twilight2day data
            l4t2d_alt = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold twilight2day data
            l4t2n_temp = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold day twilight2night
            l4t2n_alt = [None for i in range(0, len(avgr_long_timeBins_binned[keyz]['all']))]; # Hold day twilight2night
            for i in range(0, len(avgr_long_timeBins_binned[keyz]['all'])):
                passz = avgr_long_timeBins_binned[keyz]['all'][i]; # Get the pass
                mask_lt_day = avgr_long_timeBins_binned[keyz]['localtime day'][i]; # Get the day mansk
                mask_lt_nite = avgr_long_timeBins_binned[keyz]['localtime night'][i]; # Get the night mask
                mask_lt_t2d = avgr_long_timeBins_binned[keyz]['localtime twilight2day'][i]; # Get the twilight2day mask
                mask_lt_t2n = avgr_long_timeBins_binned[keyz]['localtime twilight2night'][i]; # Get the twilight2night mask
                
                tmpLat = data_saber[passz[0]]['aligned'][passz[1]]['latTan'].data.copy(); # Get the lat data out
                k = (tmpLat < avgr_long_latLim[0]) | (tmpLat > avgr_long_latLim[1]); # Get the lat mask
                
                tmpTemp = data_saber[passz[0]]['aligned'][passz[1]]['temp'].data.copy(); # Get the temp data out
                tmpTemp[k] = np.nan; # Mask the data outside of the lat range
                
                l4d_temp[i] = tmpTemp.copy(); # Get the temp data in
                l4d_temp[i][~mask_lt_day, :] = np.nan; # NaN out the non-day data
                l4n_temp[i] = tmpTemp.copy(); # Get the temp data in
                l4n_temp[i][~mask_lt_nite, :] = np.nan; # NaN out the non-night data
                l4t2d_temp[i] = tmpTemp.copy(); # Get the temp data in
                l4t2d_temp[i][~mask_lt_t2d, :] = np.nan; # NaN out the non-twilight2day data
                l4t2n_temp[i] = tmpTemp.copy(); # Get the temp data in
                l4t2n_temp[i][~mask_lt_t2n, :] = np.nan; # NaN out the non-twilight2night data
                
                tmpAlt = data_saber[passz[0]]['aligned'][passz[1]]['altTan'].data.copy(); # Get the alt data out
                tmpAlt[k] = np.nan; # Mask the data outside of the lat range

                l4d_alt[i] = tmpAlt.copy(); # Get the temp data in
                l4d_alt[i][~mask_lt_day, :] = np.nan; # NaN out the non-day data
                l4n_alt[i] = tmpAlt.copy(); # Get the temp data in
                l4n_alt[i][~mask_lt_nite, :] = np.nan; # NaN out the non-night data
                l4t2d_alt[i] = tmpAlt.copy(); # Get the temp data in
                l4t2d_alt[i][~mask_lt_t2d, :] = np.nan; # NaN out the non-twilight2day data
                l4t2n_alt[i] = tmpAlt.copy(); # Get the temp data in
                l4t2n_alt[i][~mask_lt_t2n, :] = np.nan; # NaN out the non-twilight2night data
            # END FOR i
            
            # --- Day ---
            # Combine the data to make a large bloc of data
            avgr_avgz[keyz]['localtime day alt'] = np.nanmean(np.concatenate(l4d_alt), axis=0); # Combine, immediately get the average altitude involved
            tmpSplit = np.concatenate(l4d_temp); # Combine, stack the obs together, keep 400 measurements per obs dimension

            # Assemble mean and median temps
            # tmpAvg_mean = np.nanmean(tmpSplit, axis=0);
            # Select which tmpAvg to use
            avgr_avgz[keyz]['localtime day temp'] = dist2median(tmpSplit, dist2meanMultiplier=2.5); # Alias
            
            # --- Night ---
            # Combine the data to make a large bloc of data
            avgr_avgz[keyz]['localtime night alt'] = np.nanmean(np.concatenate(l4n_alt), axis=0); # Combine, immediately get the average altitude involved
            tmpSplit = np.concatenate(l4n_temp); # Combine, stack the obs together, keep 400 measurements per obs dimension

            # Assemble mean and median temps
            # tmpAvg_mean = np.nanmean(tmpSplit, axis=0);
            # Select which tmpAvg to use
            avgr_avgz[keyz]['localtime night temp'] = dist2median(tmpSplit, dist2meanMultiplier=2.5); # Alias
            
            # --- Day ---
            # Combine the data to make a large bloc of data
            avgr_avgz[keyz]['localtime twilight2day alt'] = np.nanmean(np.concatenate(l4t2d_alt), axis=0); # Combine, immediately get the average altitude involved
            tmpSplit = np.concatenate(l4t2d_temp); # Combine, stack the obs together, keep 400 measurements per obs dimension

            # Assemble mean and median temps
            # tmpAvg_mean = np.nanmean(tmpSplit, axis=0);
            # Select which tmpAvg to use
            avgr_avgz[keyz]['localtime twilight2day temp'] = dist2median(tmpSplit, dist2meanMultiplier=2.5); # Alias
            
            # --- Day ---
            # Combine the data to make a large bloc of data
            avgr_avgz[keyz]['localtime twilight2night alt'] = np.nanmean(np.concatenate(l4t2n_alt), axis=0); # Combine, immediately get the average altitude involved
            tmpSplit = np.concatenate(l4t2n_temp); # Combine, stack the obs together, keep 400 measurements per obs dimension

            # Assemble mean and median temps
            # tmpAvg_mean = np.nanmean(tmpSplit, axis=0);
            # Select which tmpAvg to use
            avgr_avgz[keyz]['localtime twilight2night temp'] = dist2median(tmpSplit, dist2meanMultiplier=2.5); # Alias
    # END FOR keyz
    
    return avgr_avgz, avgr_long_timeBins_binned
# END DEF