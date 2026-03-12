#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
get_madrigal(dates, dataDir, madrigalInst, madrigalExp, madrigalLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False);
-  dates = [[2012,12,12],[2012,12,13],[2012,12,15]]; # yr/mo/day; can also use yr/day# format [[2012,333],[2012,334],[2012,336]]
-  dataDir = './path/where/you/want/data' the code will put year folders in this directory for all years requested, and per-year data in each year directory
-  madrigalInst = 'An Instrument'; # See https://cedar.openmadrigal.org/instMetadata under "Name" category for the name to put here
                = 21; # Direct instrument ID # also works, see https://cedar.openmadrigal.org/instMetadata under "kinst" category for the ID #s
-  madrigalExp = 'Some Cool experiment'; # Set the experiment name for the above instrument, browse Madrigal to see experiment names
               = 8000; # Direct experiment ID # also works
-  madrigalLogin = {'fullname': 'First+Last', # You must convert spaces to +
                    'email': 'your@email.here',
                    'affiliation': 'Some+Where'};
**  web_retryMax: Integer number of retries to attempt when accessing the web; setting to -1 will retry forever
**  web_retryWait: Can be an integer# to wait # sec between tries, or list/tuple [startwait, maxwait, incrementwait] where every retry the code increments the startwait by incrementwait until maxwait is reached/exceeded
**  FLG_datesRange: Set to True to indicate that dates is an inclusive range in the form yr/mo/day [[2012,12,12],[2012,12,28]] or yr/day# [[2012,333],[2012,356]] - code will automatically fill in the dates in between
**  FLG_overwrite: Set to True to redownload files even if they are already downloaded (overwrite them)
**  FLG_verbose: 2 = Info and above, 1 = Warning and above, 0 = Errors only
"""

import numpy as np
import pandas as pd
from TAS.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
# from TAS.subfun_date_to_dayNum import subfun_date_to_dayNum
from TAS.subfun_dayNum_to_date import subfun_dayNum_to_date
from TAS.get_got import get_got_webpg, get_got_webfile
from os import path as ospath, makedirs as osmakedirs
from shutil import move as shutilmove, get_terminal_size as shutil_get_terminal_size
from tempfile import gettempdir
from sys import stdout as sysstdout
from re import findall as refindall, sub as resub
from io import StringIO
import pickle as pkl
import zlib

def get_madrigal_alias():
    # Put aliases for experiment names here, such as the experiment 'Jicamarca Oblique mode' is also called 'Faraday rotation with alternating code Long Pulse' in Madrigal's FTP access system
    # Order is {'your name': 'Madrigal FTP name'}
    return {
        'Jicamarca Oblique mode': 'Faraday rotation with alternating code Long Pulse',
        'SCINDA scintillation data': 'SCINDA GPS scintillation data',
        }; # Alias setup
# END DEF

def get_madrigal_inject():
    # Put experiments and their inst #s to inject here
    return {
        'Haystack SCINDA Scintillation Receiver': ['8324', 17579],
    }
# END IF

def get_madrigal(dates, dataDir, madrigalInst, madrigalExp, madrigalLogin, web_retryMax=3, web_retryWait=[1,3,1], FLG_datesRange=False, FLG_overwrite=False, FLG_verbose=2):
    # Constants
    web_base_site = 'https://cedar.openmadrigal.org'; # Madrigal website to be used
    web_inst_site = web_base_site+'/instMetadata'; # Webpage that connects instrument to instrument #
    web_exp_site = web_base_site+'/kindatMetadata'; # Webpage that connects instrument to experiment #
    
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
    # if( dates.shape[1] == 3 ):
    #     dates_yrdn = subfun_date_to_dayNum(dates); # Convert to yr/dayNum format if in yr/mo/day
    if( dates.shape[1] == 2 ):
        # dates_yrdn = dates.copy(); # Keep this
        dates = subfun_dayNum_to_date(dates); # Convert to yr/dayNum format if in yr/mo/day
    # END IF    
    
    dates_yrs = np.unique(dates[:, 0]); # Get the unique years
    web_fileInfo = {str(dates_yrs[i]):{'date start':None,'date start obj':None, 'date end':None,'date end obj':None, 'file name':None, 'file link':None} for i in range(0, dates_yrs.size)}; # Record info into a dict
    # Get metadata file if available and record the dates already requested
    for i in range(0, dates_yrs.size):
        web_metaFile_path = ospath.join(dataDir, str(dates_yrs[i]), '.metadata.pkl');
        if( ospath.isfile(web_metaFile_path) ):
            with open(web_metaFile_path, 'rb') as pklz:
                pkl_smoosh = pklz.read(); # Read in the compressed pickle file
            # END WITH
            pkl_obj = zlib.decompress(pkl_smoosh); # Decompress the pickle file
            web_metaFile = pkl.loads(pkl_obj);  # Convert the pickle file into a usable variable
            web_fileInfo[str(dates_yrs[i])]['date requested'] = web_metaFile['date requested'].copy(); # Copy over the saved requested dates
        # END IF
    # END FOR i
    # Use metadata file to remove dates that have got-got
    if( FLG_overwrite == False ): # Only remove dates if overwrite is off
        dates_yrs_del = np.zeros(dates_yrs.size, dtype=np.bool_); # Prep
        for i in range(0, dates_yrs.size):
            if( 'date requested' in web_fileInfo[str(dates_yrs[i])] ):
                k = np.where(dates[:, 0] == dates_yrs[i])[0]; # Get which were requested
                web_fileInfo[str(dates_yrs[i])]['date requested']
                dates_alreadyHave = (web_fileInfo[str(dates_yrs[i])]['date requested'][:, None] == dates[k, :]).all(axis=2).any(axis=0); # See if we already have requested and got any dates worth of data
                dates = np.delete(dates, k[dates_alreadyHave], axis=0); # Delete data that is already got-got
                k = dates[:, 0] == dates_yrs[i]; # Get which were requested
                if( k.sum() == 0 ):
                    dates_yrs_del[i] = True; # Set to delete
                    web_fileInfo.pop(str(dates_yrs[i])); # Delete from web_fileInfo dict too
                # END IF
            # END IF
        # END FOR i
        if( np.any(dates_yrs_del) ):
            dates_yrs = dates_yrs[~dates_yrs_del]; # Delete years with no dates left
        # END IF
    # END IF
    
    # --- Get the available data from Madrigal ---
    if( dates_yrs.size > 0 ): # Only work if we need the data    
        # Initializing - Convert names to instrument/exp IDs
        # --- Load in Instrument List from Madrigal ---
        if( isinstance(madrigalInst, str) and (not madrigalInst.isdigit()) ):
            # It's the instrument name
            madrigalInst_tbl = pd.read_html(StringIO( get_got_webpg(web_inst_site, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0]; # Get the table from the HTML page
            
            madrigalInst_entry = madrigalInst_tbl.loc[madrigalInst_tbl['Name'].str.contains(madrigalInst)]; # Get the row
            if( madrigalInst_entry.size == 0 ):
                raise Exception('\nERROR in get_madrigal.py: Instrument requested `'+madrigalInst+'` was NOT found on the Madrigal instrument table website `'+web_inst_site+'` under the heading `Name`. Check spelling, it must match.');
            # END IF
            madrigalInst_id = madrigalInst_entry['Instrument id (kinst)'].item(); # Get the instrument ID
        else:
            # It's the instrument ID #
            madrigalInst_id = int(madrigalInst); # Assume it's an integer instrument ID
            
            madrigalInst_tbl = pd.read_html(StringIO( get_got_webpg(web_inst_site, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0]; # Get the table from the HTML page
            
            madrigalInst_entry = madrigalInst_tbl.loc[madrigalInst_tbl.index[madrigalInst_tbl['Instrument id (kinst)'] == madrigalInst_id]]; # Get the row
            if( madrigalInst_entry.size == 0 ):
                raise Exception('\nERROR in get_madrigal.py: Instrument ID# requested `'+str(madrigalInst_id)+'` was NOT found on the Madrigal instrument table website `'+web_inst_site+'` under the heading `kinst`. Check ID#, it must match.');
            # END IF
            madrigalInst = madrigalInst_entry['Name'].item(); # Get the instrument name
        # END IF
        
        # --- Load in Experiment List from Madrigal ---
        if( isinstance(madrigalExp, str) and (not madrigalExp.isdigit()) ):
            # It's the instrument name
            
            # Alias list for Madrigal experiment names that are named differently on the FTP for reasons
            madrigalExp_alias = get_madrigal_alias(); # In a function so it's easy to find
            
            # Apply alias
            madrigalExp_orig = madrigalExp; # Remember the original
            if( madrigalExp in madrigalExp_alias ):
                # if( FLG_verbose >= 1 ): # Not needed to even note
                #     print('WARNING in get_madrigal.py: Experiment requested `'+madrigalExp+'` has an alias name in the FTP version of Madrigal\'s site and will be reported past this point as `'+madrigalExp_alias[madrigalExp]+'`.');
                # # END IF
                madrigalExp = madrigalExp_alias[madrigalExp]; # Alias it to the useable name
            # END IF
            
            madrigalExp_tbl = pd.read_html(StringIO( get_got_webpg(web_exp_site, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0]; # Get the table from the HTML page

            injector = get_madrigal_inject(); # Get the things to inject
            frame2place = {'Kind of data file code (kindat)':[], 'Description':[],
                   'Associated instrument code':[], 'Associated instrument':[]}; # Yep
            for keyz in injector.keys():
                frame2place['Kind of data file code (kindat)'].append(injector[keyz][1]);
                frame2place['Description'].append(madrigalExp);
                frame2place['Associated instrument code'].append(injector[keyz][0]);
                frame2place['Associated instrument'].append(keyz);
            # END FOR keyz
            madrigalExp_tbl = pd.concat([madrigalExp_tbl, pd.DataFrame(frame2place)], ignore_index=True); # Tack it on
            
            # Two step to declare if instrument is missing THEN if experiment is missing
            madrigalExp_entry = madrigalExp_tbl.loc[madrigalExp_tbl['Associated instrument code'] == str(madrigalInst_id)]; # Get the row
            if( madrigalExp_entry.size == 0 ):
                raise Exception('\nERROR in get_madrigal.py: Instrument requested `'+madrigalInst+'` was NOT found on the Madrigal experiment table website `'+web_exp_site+'` under the heading `Associated instrument`. This error shouldn\'t have happened.');
            # END IF
            madrigalExp_entry = madrigalExp_entry.loc[madrigalExp_entry['Description'].str.contains(madrigalExp)]; # Get the row
            if( madrigalExp_entry.size == 0 ):
                raise Exception('\nERROR: in get_madrigal.py: Experiment requested `'+madrigalExp+'` was NOT found on the Madrigal experiment table website `'+web_exp_site+'` under the heading `Description`. Check spelling, it must match.');
            # END IF
            madrigalExp_id = madrigalExp_entry['Kind of data file code (kindat)'].unique(); # Get the instrument ID (there can be multiples)
            if( (FLG_verbose >= 1 ) and (madrigalExp_id.size > 1) ):
                print('WARNING in get_madrigal.py: multiple experiment IDs were found for the requested experiment `'+madrigalExp+'`. Listing IDs: '+str(madrigalExp_id));
            # END IF
        else:
            # It's the instrument ID #
            madrigalExp_id = int(madrigalExp); # Assume it's an integer experiment ID
            
            madrigalExp_tbl = pd.read_html(StringIO( get_got_webpg(web_exp_site, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=False, FLG_verbose=FLG_verbose) ))[0]; # Get the table from the HTML page
            
            injector = get_madrigal_inject(); # Get the things to inject
            frame2place = {'Kind of data file code (kindat)':[], 'Description':[],
                   'Associated instrument code':[], 'Associated instrument':[]}; # Yep
            for keyz in injector.keys():
                frame2place['Kind of data file code (kindat)'].append('injected');
                frame2place['Description'].append(madrigalExp);
                frame2place['Associated instrument code'].append(injector[keyz]);
                frame2place['Associated instrument'].append(keyz);
            # END FOR keyz
            madrigalExp_tbl = pd.concat([madrigalExp_tbl, pd.DataFrame(frame2place)], ignore_index=True); # Tack it on
            
            madrigalExp_entry = madrigalExp_tbl.loc[madrigalExp_tbl.index[(madrigalExp_tbl['Kind of data file code (kindat)'] == madrigalExp_id) & (madrigalExp_tbl['Associated instrument'].str.contains(madrigalInst))]].reset_index(drop=True); # Get the row (experiment IDs can be reused, so need to also specify the instrument)
            if( madrigalExp_entry.size == 0 ):
                raise Exception('\nERROR in get_madrigal.py: Experiment ID requested `'+str(madrigalExp_id)+'` was NOT found on the Madrigal experiment table website `'+web_exp_site+'` under the heading `kindat`. Check the ID#, it must match.');
            # END IF
            madrigalExp = madrigalExp_entry['Description'][0]; # Get the experiment name (only one name should be associated with an experiment ID, even if they typed it twice
        
            # Alias list for Madrigal experiment names that are named differently on the FTP for reasons
            madrigalExp_alias = get_madrigal_alias(); # In a function so it's easy to find
            # Reverse alias to go the other way
            madrigalExp_aliasRev = {valz:keyz for keyz, valz in madrigalExp_alias.items()}; # Reverse the lookup dictionary
            
            # Apply alias
            if( madrigalExp in madrigalExp_aliasRev ):
                madrigalExp_orig = madrigalExp_aliasRev[madrigalExp]; # Alias it to the FTP name
            else:
                madrigalExp_orig = madrigalExp; # Alias to keep it the same
            # END IF
            
            # Finally, make madrigalExp_id iterable (other way, naming an experiment, can yield multiple experiment IDs)
            madrigalExp_id = [madrigalExp_id]; # Iterate it up
        # END IF
        
        for exp in madrigalExp_id:
            # Build the specific page variables
            web_base = '/ftp/fullname/' + madrigalLogin['fullname'] + '/email/' + madrigalLogin['email'] + '/affiliation/' + madrigalLogin['affiliation'] + '/kinst/'+str(madrigalInst_id)+'/year/'; #year goes after this
            web_baseAfter = '/kindat/'+str(exp)+'/format/hdf5/'; #this goes after year
            
            # First - need to check if data is availiable for the days requested. Madrigal goes year-by-year
            for i in range(0, dates_yrs.size):
                # Madrigal lets you browse the FTP data availability by year, so get per-year info that can be applied for later
                web_fileSearch = web_base_site + web_base + str(dates_yrs[i]) + web_baseAfter; #build site to go to
                rendered_content = get_got_webpg(web_fileSearch, web_retryMax=web_retryMax, web_retryWait=web_retryWait, FLG_rendered=True, FLG_verbose=FLG_verbose); # Get the webpage and render it to make it easier to parse
                
                web_fileMatches = refindall(r'\[\n[\s\S]*?'+str(madrigalInst_id)+r'[\s\S]*?'+str(exp)+r'[\s\S]*?\n\n', rendered_content); # Get the matches
                # web_fileMatches = refindall(r'\.[a-zA-Z0-9]+\]\(\/ftp\/fullname[\s\S]*?'+str(madrigalInst_id)+r'[\s\S]*?'+str(exp)+r'[\s\S]*?\.[a-zA-Z0-9]+/\)[\s\S]*?\[\\', rendered_content); # Get the matches
                web_fileInfo[str(dates_yrs[i])]['date start'] = []; # Record date start
                web_fileInfo[str(dates_yrs[i])]['date start obj'] = []; # Record date start obj
                web_fileInfo[str(dates_yrs[i])]['date end'] = []; # Record date end
                web_fileInfo[str(dates_yrs[i])]['date end obj'] = []; # Record date end obj
                web_fileInfo[str(dates_yrs[i])]['file name'] = []; # Record file name
                web_fileInfo[str(dates_yrs[i])]['file link'] = []; # Record file link
                web_fileInfo[str(dates_yrs[i])]['exp name'] = []; # Record real exp name (specific to this exp)
                tmp_idx2del = []; # Prep
                for j in range(0, len(web_fileMatches)):
                    if( not 'not for science use' in web_fileMatches[j] ):
                        # Extract the file name
                        web_fileInfo[str(dates_yrs[i])]['file name'].append( web_fileMatches[j][2:web_fileMatches[j].find(']')] );
                        # Build the link
                        web_fileInfo[str(dates_yrs[i])]['file link'].append( web_base_site + web_fileMatches[j][web_fileMatches[j].find(']') + 2:web_fileMatches[j].rfind(web_fileInfo[str(dates_yrs[i])]['file name'][-1]) + len(web_fileInfo[str(dates_yrs[i])]['file name'][-1]) + 1] );
                        # Get the start date
                        web_fileInfo[str(dates_yrs[i])]['date start'].append( web_fileMatches[j][web_fileMatches[j].find(' From ') + 6:web_fileMatches[j].find(' to ')] );
                        web_fileInfo[str(dates_yrs[i])]['date start obj'].append( np.datetime64(web_fileInfo[str(dates_yrs[i])]['date start'][-1]) ); # Convert to a datetime64 object for ease of use
                        # Get the end date
                        web_fileInfo[str(dates_yrs[i])]['date end'].append( web_fileMatches[j][web_fileMatches[j].find(' to ') + 4:web_fileMatches[j].find(': ')] );
                        web_fileInfo[str(dates_yrs[i])]['date end obj'].append(np.datetime64( web_fileInfo[str(dates_yrs[i])]['date end'][-1]) ); # Convert to a datetime64 object for ease of use
                        # Get the real exp name
                        web_fileInfo[str(dates_yrs[i])]['exp name'].append( web_fileMatches[j][web_fileMatches[j].find(web_fileInfo[str(dates_yrs[i])]['date end'][-1]) + len(web_fileInfo[str(dates_yrs[i])]['date end'][-1]) + 1:web_fileMatches[j].rfind(r'[\n')].replace('\n',' ').strip(' ') ); # Remove spaces and replace line breaks to help with working with it later
                    else:
                        tmp_idx2del.append(j); # Delete
                    # END IF
                # END FOR j
                for j in tmp_idx2del[::-1]:
                    web_fileMatches.pop(j); # Delete
                # END FOR j
                if( len(web_fileMatches) > 0 ):
                    # Remove Prelims that have a Final version
                    tmp_prelim = []; # Prep
                    tmp_fin = [];
                    for j in range(0, len(web_fileInfo[str(dates_yrs[i])]['exp name'])):
                        # Detect identically named experiments and remove ones that are preliminary that have a final version
                        if( ': Preliminary' in web_fileInfo[str(dates_yrs[i])]['exp name'][j] ):
                            tmp_prelim.append( j ); # Count em up
                        elif( ': Final' in web_fileInfo[str(dates_yrs[i])]['exp name'][j] ):
                            tmp_fin.append( web_fileInfo[str(dates_yrs[i])]['exp name'][j][:web_fileInfo[str(dates_yrs[i])]['exp name'][j].find(': ')].strip(' ') ); # Count em up
                        # END IF
                    # END FOR j
                    for j in tmp_prelim[::-1]:
                        if( web_fileInfo[str(dates_yrs[i])]['exp name'][j][:web_fileInfo[str(dates_yrs[i])]['exp name'][j].find(': ')].strip(' ') in tmp_fin ):
                            for subkeyz in web_fileInfo[str(dates_yrs[i])].keys():
                                web_fileInfo[str(dates_yrs[i])][subkeyz].pop(j); # Ditch it
                            # END FOR subkeyz
                            web_fileMatches.pop(j); # Ditch it
                        # END IF
                    # END FOR j
                    # # Remove Identical time runs (likely more prelim but they changed the name between prelim and final) [didn't get to this sorry future me]
                    tmp_dtstrt = np.asarray(web_fileInfo[str(dates_yrs[i])]['date start obj']); # Conv
                    tmp_dtstrt_uniq, tmp_dtstrt_uniqIdx, tmp_dtstrt_uniqCnt = np.unique(tmp_dtstrt, return_index=True, return_counts=True);
                    if( tmp_dtstrt.size != tmp_dtstrt_uniq.size ):
                        k = np.where(tmp_dtstrt_uniqCnt > 1 )[0]; # Get where duplicates
                        for jk in k[::-1]:
                            kk = np.where(tmp_dtstrt_uniq[jk] == tmp_dtstrt)[0]; # Get the duplicates
                            # Check for identical name (repeated entries that look identical otherwise so they're prob identical)
                            tmp_ditchr = []; # Prep
                            for j in range(1, kk.size):
                                if( (web_fileInfo[str(dates_yrs[i])]['exp name'][kk[0]] == web_fileInfo[str(dates_yrs[i])]['exp name'][kk[j]]) and (web_fileInfo[str(dates_yrs[i])]['file name'][kk[0]] == web_fileInfo[str(dates_yrs[i])]['file name'][kk[j]]) ):
                                    tmp_ditchr.append( kk[j] ); # Tack on
                                # EN DIF
                            # END FOR j
                            if( len(tmp_ditchr) > 0 ):
                                for j in tmp_ditchr[::-1]:
                                    for subkeyz in web_fileInfo[str(dates_yrs[i])].keys():
                                        web_fileInfo[str(dates_yrs[i])][subkeyz].pop(j); # Ditch it
                                    # END FOR subkeyz
                                    web_fileMatches.pop(j); # Ditch it
                                # END FOR j
                            # END IF
                        # END FOR jk
                    # END IF
                    tmp_dtend = np.asarray(web_fileInfo[str(dates_yrs[i])]['date end obj']); # Conv
                    tmp_dtend_uniq, tmp_dtend_uniqIdx, tmp_dtend_uniqCnt = np.unique(tmp_dtend, return_inverse=True, return_counts=True);
                    if( tmp_dtend.size != tmp_dtend_uniq.size ):
                        k = np.where(tmp_dtend_uniqCnt > 1 )[0]; # Get where duplicates
                        for jk in k[::-1]:
                            kk = np.where(tmp_dtend_uniq[jk] == tmp_dtend)[0]; # Get the duplicates
                            # Check for identical name (repeated entries that look identical otherwise so they're prob identical)
                            tmp_ditchr = []; # Prep
                            for j in range(1, kk.size):
                                if( (web_fileInfo[str(dates_yrs[i])]['exp name'][kk[0]] == web_fileInfo[str(dates_yrs[i])]['exp name'][kk[j]]) and (web_fileInfo[str(dates_yrs[i])]['file name'][kk[0]] == web_fileInfo[str(dates_yrs[i])]['file name'][kk[j]]) ):
                                    tmp_ditchr.append( kk[j] ); # Tack on
                                # EN DIF
                            # END FOR j
                            if( len(tmp_ditchr) > 0 ):
                                for j in tmp_ditchr:
                                    for subkeyz in web_fileInfo[str(dates_yrs[i])].keys():
                                        web_fileInfo[str(dates_yrs[i])][subkeyz].pop(j); # Ditch it
                                    # END FOR subkeyz
                                    web_fileMatches.pop(j); # Ditch it
                                # END FOR j
                            # END IF
                        # END FOR jk
                    # END IF
                # END IF
            # END FOR i
            
            # Second - Compare dates available with the requested range to determine if want to download
            web_fileGet = {keyz:[False for _ in range(0, len(web_fileInfo[keyz]['date start']))] for keyz in web_fileInfo.keys()}; # Record info into a dict
            for i in range(0, dates.shape[0]):
                # Get the date range for the data
                dates_start = np.datetime64('-'.join(strang.zfill(2) for strang in resub(r'\s+', '-', str(dates[i, :])[1:-1]).split('-'))); # Get the start date (hooooooboy on those python shennanigans)
                dates_end = dates_start + np.timedelta64(1, 'D'); # By definition the end date is just before the beginning of the next day, represent via non-inclusive comparisons (e.g., no 'or equals')
                
                # Check if the data has that date range
                for j in range(0, len(web_fileGet[str(dates[i, 0])])):
                    web_dateStart = web_fileInfo[str(dates[i, 0])]['date start obj'][j].astype('datetime64[D]'); # Round down to nearest day
                    web_dateEnd = (web_fileInfo[str(dates[i, 0])]['date end obj'][j] + np.timedelta64(1, 'D')).astype('datetime64[D]'); # Round up to nearest day (by adding a day and forcing to days) - non-inclusive here since ends on previous day just before the next
                    
                    web_fileGet[str(dates[i, 0])][j] |= (web_dateStart <= dates_start) & (web_dateEnd >= dates_end); # Record to include it with an or statement
                # END FOR j
            # END FOR i
            
            # Third - Download if not already available
            web_dlSpeed = None; # Prep DL speed calcs
            for keyz in web_fileGet.keys():
                if( FLG_verbose >= 2 ):
                    colWidth = shutil_get_terminal_size().columns//2; # Get half the column width
                # END IF
                # Prep tmp dir so that files that fail mid-download are not stored in the dataset directly
                dataDir_tmp = gettempdir();
                # Make sure year data folder there
                if( ospath.isdir(ospath.join(dataDir, keyz)) == False ):
                    osmakedirs(ospath.join(dataDir, keyz)); # Make the year folder if it doesn't exist
                # END IF
                if( len(web_fileGet[keyz]) > 0 ):
                    for j in range(0, len(web_fileGet[keyz])):
                        if( web_fileGet[keyz][j] == True ):
                            web_filePath = ospath.join(dataDir, keyz, web_fileInfo[keyz]['file name'][j]); # Build the full file path
                            web_fileTmp = ospath.join(dataDir_tmp, web_fileInfo[keyz]['file name'][j]); # Build the temporary file path
                            if( (ospath.isfile(web_filePath) == False) or (FLG_overwrite == True) ):
                                # If the file isn't already there, get it
                                web_dlSpeed = get_got_webfile(web_fileInfo[keyz]['file link'][j], web_fileTmp, web_dlSpeed=web_dlSpeed, web_retryMax=3, web_retryWait=[1,3,1], FLG_getFileSize=False, FLG_showDownloadedAmnt=False, FLG_verbose=FLG_verbose);
                                # If that suceeds, we make it to the move
                                shutilmove(web_fileTmp, web_filePath); # Move the file from temp to its real destination
                            # END IF
                        # END IF   
                        if( FLG_verbose >= 2 ):
                            donezo = round((j+1)/len(web_fileGet[keyz])*colWidth); # Calc how close to being done
                            sysstdout.write('\rProgress: '+'■'*donezo+'□'*(colWidth - donezo)); #report
                            sysstdout.flush();
                        # END IF
                    # END FOR j
                else:
                    if( FLG_verbose >= 1 ):
                        print('WARNING in get_madrigal.py: Requested year '+keyz+' for Instrument/Exp combo `'+str(madrigalInst)+'`/`'+str(madrigalExp)+'` (instID#'+str(madrigalInst_id)+'/expID#'+str(exp)+') has no data available for the dates requested:\n'+str(dates[dates[:, 0] == int(keyz)]));
                    # END IF
                # END IF
            # END FOR keyz            
        # END FOR exp
        
        # Record file about what was requested
        for keyz in web_fileInfo.keys():
            if( ospath.isdir( ospath.join(dataDir, keyz) ) ): # Only make a metadata file if the year is there
                # Metadata file path
                web_metaFile_path = ospath.join(dataDir, keyz, '.metadata.pkl');
                k = np.where(dates[:, 0] == np.int16(keyz))[0]; # Get which were requested
                if( 'date requested' in web_fileInfo ):
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
                pkl_obj = pkl.dumps(web_fileInfo[keyz], protocol=pkl.HIGHEST_PROTOCOL); # Pickle object
                pkl_smoosh = zlib.compress(pkl_obj); # Compress the object
                with open(web_metaFile_path, 'wb') as pklz:
                    pklz.write(pkl_smoosh); # Write the compressed pickle to disk
                # END WITH
            # END IF
        # END FOR keyz
    else:
        if( FLG_verbose >=2 ):
            print('INFO in get_madrigal.py: Data has been downloaded already for all dates requested in `'+dataDir+'`.');
        # END IF
    # END IF
# END DEF