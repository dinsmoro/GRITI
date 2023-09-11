"""
#===== LISN_h5_write =====
# Discovers & reads LISN data files then writes them to a single compressed HDF5 file on a per-day basis
# Options:
# path_data - supply a path to the folder that has the LISN original data files as a string, if not supplied code will run on the directory the LISN_h5.py file located in
# path_save - folder to save data to, if left empty is set to path_data
# FLG_recuurisve - if True code will search subdirectories for LISN files, otherwise it will stay within the given directory
# fileNameConvention - file name convention to use when saving the HDF5 file.
    #{YR} inserts the year, {Y} inserts the last 2 year digits, {DN} inserts the day number zfilled to 3 digits, 
    #{MO} inserts the month zfilled to 2 digits, {M} inserts the month not zfilled, {DY} inserts the day zfilled to 2 digits, {D} inserts the day not zfilled
# fileNameEnding - the expected file ending for the orignal LISN data files
# floatType - float precision to use, float32 takes 1/2 space in RAM/disk of float64 and is definitely good enough for TEC precision
# FLG_deleteOrig - if True code will delete the original LISN data files after the HDF5 file has been written, use with caution
# FLG_overwrite - if True code will overwrite already created HDF5 files, if False it will skip already created HDF5 files
# FLG_verbose - if True code will output info into wherever text goes
# RD on 2/21/2022, GPL-3.0 licensed via base of https://github.com/dinsmoro/GRITI
# Example:
#     LISN_h5_write(path_data='E:/SomePath/Data', path_save='E:/SomePath/Save', FLG_recursive=False, fileNameConvention='LISN_{YR}_{DN}.h5', fileNameEnding='.dat', floatType='float32', FLG_deleteOrig=False, FLG_overwrite=False, FLG_verbose=True)
#CMDLINE: python -c "from LISN_h5 import LISN_h5_write; LISN_h5_write(path_data='E:/SomePath/Data', path_save='E:/SomePath/Save', FLG_recursive=False, fileNameConvention='LISN_{YR}_{DN}.h5', fileNameEnding='.dat', floatType='float32', FLG_deleteOrig=False, FLG_overwrite=False, FLG_verbose=True)"
"""
#==============================================================================
"""
#===== LISN_h5_read =====
# Reads pre-processed LISN HDF5 files for a specific date range, if provided (otherwise reads em all)
# Options:
# dateRange - date range in format [[YR_start,DN_start],[YR_end,DN_end]] OR [[YR_start,MO_start,DY_start],[YR_end,MO_end,_DY_end]] can be supplied as a numpy array, list, or tuple
# path_read - supply a path to the folder that has the LISN HDF5 files as a string, if not supplied code will run on the directory the LISN_h5.py file located in
# FLG_recursive - if True code will search subdirectories for LISN files, otherwise it will stay within the given directory
# fileNameConvention - file name convention to use when saving the HDF5 file.
    #{YR} inserts the year, {Y} inserts the last 2 year digits, {DN} inserts the day number zfilled to 3 digits, 
    #{MO} inserts the month zfilled to 2 digits, {M} inserts the month not zfilled, {DY} inserts the day zfilled to 2 digits, {D} inserts the day not zfilled
# FLG_combine - if True code will combine all days into a single dictionary, concatenating same arrays together. False will yield a list of dictionaries for each day
# FLG_ignoreSubGroups - if True will ignore sub-groups (that hold individual site info like receiver site name/lat/long/alt/make) but will speed up the file reading by about 2 seconds.
# file_maxRetry - number of times to retry reading an HDF5 file in the event of a weird OS error, set to -1 to try forever
# FLG_verbose - if True code will output info into wherever text goes
# RD on 2/21/2022, GPL-3.0 licensed via base of https://github.com/dinsmoro/GRITI
# Example:
#    LISN_h5_read(dateRange=[[2017,73],[2017,75]], path_read='E:/SomePath/Save', FLG_recursive=False, fileNameConvention='LISN_{YR}_{DN}.h5', FLG_combine=True, FLG_ignoreSubGroups=False, file_maxRetry=3, FLG_verbose=True)
"""
#==============================================================================
"""
#===== LISN_h5_version =====
# Returns the LISN HDF5 file version it's in a function so _write and _read both know the latest verison and so you don't have to update it in two places
# RD on 2/21/2022, GPL-3.0 licensed via base of https://github.com/dinsmoro/GRITI
"""

import numpy as np
import pandas as pd
import h5py
import sys, os
import pathlib
import time
import re


def LISN_h5_version(): #when LISN_h5_write is updated to change the HDF5 file format, increment version here and note changes
    # Increase integer value in version for breaking format changes, otherwise increment decimal for fixes that don't introduce completely new file formatting
    version_LISN = 1.0; # this value is embedded in the HDF5 file and can be used to validate if HDF5 is old and needs to be created again
    #1.0 2/21/2022 - initial LISN_h5_write version
    return version_LISN;
#END DEF


def LISN_h5_write(path_data=None, path_save=None, FLG_recursive=False, fileNameConvention='LISN_{YR}_{DN}.h5', fileNameEnding='.dat', floatType='float32', FLG_deleteOrig=False, FLG_overwrite=False, FLG_verbose=True):
    #--- Define constants ---
    version_LISN = LISN_h5_version(); #get current version (defined by function at top)
    pierceAlt_LISN = 350; #km, import pierce point altitude - defined as a constant based on communications
    if( floatType == 'float64' ): # to change float precision saved (float32 takes 1/2 the size of float64)
        dataType_meth = np.float64;
        dataType_str = 'float64';
    else:
        dataType_meth = np.float32;
        dataType_str = 'float32';
    #END IF
    
    
    #--- Align current working directory ---
    if( path_data != None):
        os.chdir(path_data); # Update CWD to provided dir
    else:
        # If fldr not supplied, force dir to be where LISN_h5.py is
        os.chdir(pathlib.Path('.').resolve()); # Update CWD to where file is
    #END IF
    path_data = os.getcwd(); # Record CWD
    if( path_save == None ):
        path_save = path_data; #update path_save to be path_data
    #END IF
    if( FLG_verbose == True ):
        print('--- LISN_h5_write ---');
        print('Data directory: '+path_data+'\nSave directory: '+path_save); # Report CWD
    #END IF
    
    
    #--- Discover files ---
    # Only have a few static things about the files to work with - they end with .dat, have 4 letter/# site names, and have some #s in there
    if( FLG_recursive == False ):
        filez = [os.path.join(path_data,file) for file in os.listdir(path_data) if file.lower().endswith(fileNameEnding)]; # Get list of .dat files (NOT recursive), inspired by https://stackoverflow.com/a/3964696/2403531
    else:
        filez = [os.path.join(root, file) for root, dirs, files in os.walk(path_data) for file in files if file.lower().endswith(fileNameEnding)]; # Get list of .dat files (recursive!!), inspired by https://stackoverflow.com/a/3964690/2403531
    #END IF
    
    # Remove files that aren't LISN data files but do have a .dat ending, since that is a common ending
    filez_remove = set(); # Create empty set
    for i in range(0,len(filez)):
        if( filez[i][filez[i].rfind(os.path.sep)+len(os.path.sep)+4:-4].isdigit() == False):
            filez_remove.add(i); # Add the index to the set if where the date is isn't only numbers
        #END IF        
    #END FOR i
    filez = [i for j, i in enumerate(filez) if j not in filez_remove]; # Remove filez that don't seem to be LISN data files
    
    # Record the site & date value
    filez_site = np.empty(len(filez), dtype='S4'); # Preallocate, needs to be in S type which causes b'25ma' b/c that's the type that HDF5 can save for strings
    filez_dateNum = np.empty((len(filez),2), dtype=np.int16); # Preallocate [yr, dayNum]
    filez_dateNumComb = np.empty(len(filez), dtype=np.int32); # Preallocate, used as an ez unique identifier
    for i in range(0,len(filez)):
        filez_name = filez[i][filez[i].rfind(os.path.sep)+len(os.path.sep):-4]; # Get the filez name only
        filez_site[i] = np.string_(filez_name[0:4]); # Get site name
        filez_dateNum[i,0] = np.int16(filez_name[4:6]); # Get year
        filez_dateNum[i,1] = np.int16(filez_name[6:]); # Get dayNum
        filez_dateNumComb[i] = np.int32(filez_name[4:]); # Get combined value which is unique to itself and can be easily used
    #END FOR i
    filez_dateNum[:,0] += 2000; # Make year full (assume no 90's data)
    filez_date = subfun_dayNum_to_date(filez_dateNum); # Convert to month/day format
    
    
    #--- Read, amalgamate, & write ---
    current_timePerVect = 0.25; # sec, time per file (total guess to start off)
    filez_unique, filez_uniqueIndx, filez_uniqueRev = np.unique(filez_dateNumComb, return_index=True, return_inverse=True); # Get uniques and the reverse indexes to create the full array
    for i in range(0,filez_unique.size):
        # Get current files working on
        filez_currIndx = np.where(filez_uniqueIndx[i] == filez_uniqueRev)[0]; # Get where current files are
        
        # Create new file name
        filez_outputName = fileNameConvention.replace('{YR}', str(filez_dateNum[filez_uniqueIndx[i],0]) ); #replace any {YR} with the current year
        filez_outputName = filez_outputName.replace('{Y}', str(filez_dateNum[filez_uniqueIndx[i],0])[2:] ); #replace any {Y} with the current year's last 2 digits
        filez_outputName = filez_outputName.replace('{DN}', str(filez_dateNum[filez_uniqueIndx[i],1]).zfill(3) ); #replace any {DN} with the current day number (padded w/ 0's so always 3 long)
        filez_outputName = filez_outputName.replace('{MO}', str(filez_date[filez_uniqueIndx[i],1]).zfill(2) ); #replace any {MO} with the current month (padded w/ 0's so always 2 long)
        filez_outputName = filez_outputName.replace('{M}', str(filez_date[filez_uniqueIndx[i],1]) ); #replace any {M} with the current month
        filez_outputName = filez_outputName.replace('{DY}', str(filez_date[filez_uniqueIndx[i],2]).zfill(2) ); #replace any {DY} with current day (padded w/ 0's so always 2 long)
        filez_outputName = filez_outputName.replace('{D}', str(filez_date[filez_uniqueIndx[i],2]) ); #replace any {D} with current day
        filez_outputName = os.path.join(path_save,filez_outputName); # Make it a full path
        
        if( os.path.isfile(filez_outputName) ):
            try:
                with h5py.File(filez_outputName, 'r') as h5_LISN:
                    version_LISN_fromFile = h5_LISN.attrs['version']; #read the version
                #END WITH
                if( version_LISN_fromFile == version_LISN ):
                    FLG_alreadyExists = True; # File does exist
                elif( version_LISN_fromFile < version_LISN ):
                    FLG_alreadyExists = False; # FIle is old and obsolete, remake it
                    print('WARNING in LISN_h5_write: LISN file version '+str(version_LISN_fromFile)+' is older than this code\'s version '+str(version_LISN)+'. Renaming the old file '+filez_outputName+'_oldV'+str(version_LISN_fromFile).replace('.','p')+' and remaking it.');
                    os.rename(filez_outputName, filez_outputName+'_oldV'+str(version_LISN_fromFile).replace('.','p')); #rename
                else:
                    print('WARNING in LISN_h5_write: LISN file version '+str(version_LISN_fromFile)+' is newer than this code\'s version '+str(version_LISN)+'. Update this code.');
                    FLG_alreadyExists = True; # File does exist
                #END IF
            except:
                FLG_alreadyExists = False; # File is bjorked and can't be read
            #END TRY
        else:
            FLG_alreadyExists = False; # File does not already exist
        #END IF
        
        if( (FLG_alreadyExists == False) | (FLG_overwrite == True) ):
            dataNum = filez_currIndx.size; #total number of steps
            if( FLG_verbose == True ):
                print('Converting raw LISN data for date '+str(filez_dateNum[filez_uniqueIndx[i],0])+'/'+str(filez_dateNum[filez_uniqueIndx[i],1])+' ('+str(filez_date[filez_uniqueIndx[i],0])+'/'+str(filez_date[filez_uniqueIndx[i],1])+'/'+str(filez_date[filez_uniqueIndx[i],2])+') to '+filez_outputName+
                      '.\nAt '+str(dataNum)+' files and '+str(current_timePerVect)+' sec per file, expect '+textNice(np.round(current_timePerVect*dataNum/60,2))+' minutes for conversion to finish.\n');
                tic = time.time(); #for time testing
                estimatedUpdates = np.arange(np.int64(dataNum*.1),dataNum-np.int64(dataNum*.1),np.int64(dataNum*.1));
            #END IF
            
            dict_LISN = {
                'version':version_LISN, # Attributes here, code automagically saves them too
                'site info':{}, # Dict that holds info about each site that's in the site header file
                '#':[], #1 - Consecutive numbers
                'day num J2000':[], #2 – day number from January 1, 2000 = 1
                'sat':[], #PRN mid-header line
                'hour':[], #3, 4, 5 – hour, min, sec (UT)
                'min':[],
                'sec':[],
                'LT':[], #7 – local time
                'lat':[], #9 – latitude of sub-ionospheric point
                'long':[], #10 – Longitude
                'elev':[], #8 – elevation of looking direction
                'az':[], #12 – Azimuth
                'sTEC':[], #11 – slant TEC
                'vTEC':[], #6 – vTEC
                'site':[],
                }; #prep a dict to hold the data as it's read in (lists expand better than numpy arrays)
                        
            for j in range(0,dataNum):
                #Import raw LISN TEC data file
                with open(filez[filez_currIndx[j]], 'r') as LISN_file: #open file pointer (with for safety cause life is hard apparently)
                    TEC_raws = LISN_file.readlines(); #multiple headers of varying sizes make reading this real annoying
                #END with
                
                num_dataEntries_PRNnum = 0; #prep counter
                num_dataEntries_headerLen = 0; #prep counter
                for k in range(0, len(TEC_raws)):
                    num_dataEntries_PRNindex = strstr(TEC_raws[k],'PRN');
                    if( (num_dataEntries_headerLen == 0) & (num_dataEntries_PRNindex.size > 0) ):
                        num_dataEntries_headerLen  = k; #record the header length
                        #the header length is found by the first instance of PRN, which signals the data is starting
                        break;
                    #END IF
                    # if( num_dataEntries_PRNindex.size > 0 ):
                    #     num_dataEntries_PRNnum = num_dataEntries_PRNnum + 1; #increment
                    # #END IF
                #END FOR k
                # num_dataEntries = len(TEC_raws) - num_dataEntries_PRNnum - num_dataEntries_headerLen; #get the number of real data entries
                
                #Get site location
                header_siteLocation = TEC_raws[0].lstrip(' ').rstrip('\n').rstrip(' '); # Get site location
                
                #Get lat/long/alt/4th and 5th params that idk
                header_two = list(filter(None, TEC_raws[1].split(' ')));
                header_lat = dataType_meth(header_two[0].lstrip(' ').rstrip('\n'));
                header_long = dataType_meth(header_two[1].lstrip(' ').rstrip('\n'));
                header_alt = dataType_meth(header_two[2].lstrip(' ').rstrip('\n'));
                header_4thParam = dataType_meth(header_two[3].lstrip(' ').rstrip('\n'));
                header_5thParam = dataType_meth(header_two[4].lstrip(' ').rstrip('\n'));
                
                #Get receiver brand
                header_recBrand = TEC_raws[2].lstrip(' ').rstrip('\n').rstrip(' '); # Get receiver brand
                
                #Get RINEX file type
                header_fileType = TEC_raws[3].lstrip(' ').rstrip('\n').rstrip(' '); # Get RINEX file type
                
                #Now, collect the date line from the header and confirm it's the correct day
                dateLine = np.where(strfind(TEC_raws[0:num_dataEntries_headerLen], str(filez_date[filez_uniqueIndx[i],0])) == 1)[0]; #get the index of the date line
                if( dateLine.size == 0 ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the year "+str(filez_date[filez_uniqueIndx[i],0])+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.crash();
                elif( dateLine.size > 1 ):
                    #in this instance one of the location numbers happened to have the year in it, so choose the last one
                    dateLine = dateLine[-1];
                #END IF
                dateLine = TEC_raws[dateLine.item()].split(' '); #split by a space
                dateLine = list(filter(None, dateLine)); #remote empty entries (that were just a space)
                dateLine_slashes = np.where(strfind(dateLine,'\\') == 1)[0]; #get the index of strings with \'s in them from \n or \r
                if( dateLine_slashes.size == 1 ):
                    dateLine_slashesIndex = strstr(dateLine[dateLine_slashes.item()],'\\'); #get the index of where the \ occurs in the string
                    dateLine[dateLine_slashes.item()] = dateLine[dateLine_slashes.item()][0:dateLine_slashesIndex.item()]; #delete the \
                elif( dateLine_slashes.size > 1 ):
                    for k in range(0, dateLine_slashes.size):
                        dateLine_slashesIndex = strstr(dateLine[dateLine_slashes[k]],'\\'); #get the index of where the \ occurs in the string
                        dateLine[dateLine_slashes[k]] = dateLine[dateLine_slashes[k]][0:dateLine_slashesIndex.item()]; #delete the \
                    #END FOR k
                #END IF
                if( dateLine[0] != str(filez_date[filez_uniqueIndx[i],0]) ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the correct year "+str(filez_date[filez_uniqueIndx[i],0])+". The year found was "+dateLine[0]+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.crash();
                #END IF
                if( dateLine[1] != str(filez_date[filez_uniqueIndx[i],1]) ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the correct month "+str(filez_date[filez_uniqueIndx[i],1])+". The month found was "+dateLine[1]+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.crash();
                #END IF
                if( dateLine[2] != str(filez_date[filez_uniqueIndx[i],2]) ):
                    #catch an error
                    print("\n==============ERROR==============");
                    print("Header does not contain the correct day "+str(filez_date[filez_uniqueIndx[i],2])+". The day found was "+dateLine[2]+". If the file is mis-filed, convert this to a warning and add an ability to skip. Otherwise exiting.");
                    sys.crash();
                #END IF
                
                #read the file with pandas which is much faster than line by line then do some shennanigans to deal with mid-file PRN header lines
                rawData = pd.read_csv(filez[filez_currIndx[j]],
                                delim_whitespace=True, skipinitialspace=True, header=None, 
                                names=['#','dayNum_J2000','hour','min','sec','vTEC','LT','elev','lat','long','sTEC','az'],
                                dtype={'#':'object','dayNum_J2000':np.int64,'hour':np.int64,'min':np.float64,'sec':np.float64,
                                       'vTEC':np.float64,'LT':np.float64,'elev':np.float64,'lat':np.float64,
                                       'long':np.float64,'sTEC':np.float64,'az':np.float64},
                                skiprows=np.arange(0,num_dataEntries_headerLen)); #read w/ pandas
                midHeaders = np.where(np.isnan(rawData['az']))[0]; #get mid header locs by looking for NaN at last column 
                midHeaders_extra = np.append(midHeaders,rawData['az'].size); #add on last index
                num_dataEntries = rawData['az'].size-midHeaders.size; #get the number of data entries
                temp_sat = np.empty(rawData['az'].size, dtype='int16'); #preallocate new column (full size to match indexing despite knowing exact # of data entries)
                for k in range(0,midHeaders.size):
                    temp_sat[midHeaders_extra[k]:midHeaders_extra[k+1]] = rawData['dayNum_J2000'][midHeaders[k]]; #broadcast in the PRN #
                #END FOR k
                temp_sat = np.delete(temp_sat, midHeaders); #ditch the header lines (don't need to incorporate this vector into pandas)
                rawData = rawData.dropna(thresh=4, axis=0); #ditch header lines via dropping nans
                # temp_sat = np.empty( num_dataEntries, dtype='int16'); #preallocate
                temp_numIndx = np.int32(rawData['#']); #load in
                temp_dayNumJ2000 = np.int32(rawData['dayNum_J2000']); #load in
                temp_hour = np.int16(rawData['hour']); #load in
                temp_min = np.int16(rawData['min']); #preallocate
                temp_sec = np.int16(rawData['sec']); #preallocate
                temp_LT = np.int16(rawData['LT']); #preallocate
                temp_lat = dataType_meth(rawData['lat']); #preallocate
                temp_long = dataType_meth(rawData['long']); #preallocate
                temp_elev = dataType_meth(rawData['elev']); #preallocate
                temp_az = dataType_meth(rawData['az']); #preallocate
                temp_sTEC = dataType_meth(rawData['sTEC']); #preallocate
                temp_vTEC = dataType_meth(rawData['vTEC']); #preallocate
                # temp_sTECerr = np.zeros( num_dataEntries, dtype='float32'); #preallocate
                temp_site = np.empty( num_dataEntries, dtype='S4'); #preallocate
                temp_site[:] = filez_site[filez_currIndx[j]]; #record the site name
                        
                # #old line-by-line which is v slow                
                # temp_sat = np.empty( num_dataEntries, dtype='int16'); #preallocate
                # temp_hour = np.empty( num_dataEntries, dtype='int16'); #preallocate
                # temp_min = np.empty( num_dataEntries, dtype='int16'); #preallocate
                # temp_sec = np.empty( num_dataEntries, dtype='int16'); #preallocate
                # temp_lat = np.empty( num_dataEntries, dtype=dataType_str); #preallocate
                # temp_long = np.empty( num_dataEntries, dtype=dataType_str); #preallocate
                # temp_elev = np.empty( num_dataEntries, dtype=dataType_str); #preallocate
                # temp_sTEC = np.empty( num_dataEntries, dtype=dataType_str); #preallocate
                # temp_vTEC = np.empty( num_dataEntries, dtype=dataType_str)*np.nan; #preallocate
                # # temp_sTECerr = np.zeros( num_dataEntries, dtype='float32'); #preallocate
                # temp_site = np.empty( num_dataEntries, dtype='S4'); #preallocate
                
                # #record known constants
                # temp_site[:] = fileNames_Good[j][0:4]; #record the site name
                
                # lineCounter = 0; #counts the data lines, ignores the PRN lines
                # for k in range(num_dataEntries_headerLen, len(TEC_raws)):
                #     PRN_check = strstr(TEC_raws[k],'PRN'); #a PRN is a satellite ID #, for GPS it's 1-32 I think, may change w/ constellation and time
                #     if( PRN_check.size > 0 ):
                #         #found a PRN delcaration line, time to record it
                #         current_PRN = TEC_raws[k].split(' '); #split by space
                #         current_PRN = list(filter(None, current_PRN)); #remote empty entries (that were just a space)
                #         current_PRN = np.int16(current_PRN[ np.asscalar(np.where(strfind(current_PRN, 'PRN') == 1)[0]) + 1 ]); #pulls out the current PRN based on finding where the phrase 'PRN' is and then getting the number right after that
                #     else:
                #         #format the line for being able to use it
                #         current_line = TEC_raws[k].split(' '); #split by space
                #         current_line = list(filter(None, current_line)); #remote empty entries (that were just a space)
                #         #run code to remove \'s
                #         current_slashes = np.where(strfind(current_line,'\\') == 1)[0]; #get the index of strings with \'s in them from \n or \r
                #         if( current_slashes.size == 1 ):
                #             current_slashesIndex = strstr(current_line[np.asscalar(current_slashes)],'\\'); #get the index of where the \ occurs in the string
                #             current_line[np.asscalar(current_slashes)] = current_line[np.asscalar(current_slashes)][0:np.asscalar(current_slashesIndex)]; #delete the \
                #         elif( current_slashes.size > 1 ):
                #             for k in range(0, current_slashes.size):
                #                 current_slashesIndex = strstr(current_line[current_slashes[k]],'\\'); #get the index of where the \ occurs in the string
                #                 current_line[current_slashes[k]] = current_line[current_slashes[k]][0:np.asscalar(current_slashesIndex)]; #delete the \
                #             #END FOR k
                #         #END IF
                #         if( len(current_line) != 12 ):
                #             print("\n==============ERROR==============");
                #             print("Reading line #"+str(k+1)+" in "+fileNames_Good[j]+" has found less than 12 numbers. Develop code to deal with this.");
                #             sys.crash();
                #         #END IF
                        
                #         #pull out data from the line
                #         current_hr = np.int16(current_line[2]); #get the current hour (UT)
                #         current_min = np.int16(current_line[3]); #get the current min (UT)
                #         current_sec = np.int16(current_line[4]); #get the current sec (UT)
                #         current_vTEC = np.float32(current_line[5]); #get the current vTEC (TECU)
                #         current_elev = np.float32(current_line[7]); #get the current elevation angle (deg)
                #         current_lat = np.float32(current_line[8]); #get the current pierce-point latitude (arcdeg)
                #         current_long = np.float32(current_line[9]); #get the current pierce-point longitude (arcdeg)
                #         current_sTEC = np.float32(current_line[10]); #get the current sTEC (TECU)
                        
                #         #record the data found
                #         #record ints
                #         temp_sat[lineCounter] = current_PRN; #record the current satellite ID (PRN)
                #         temp_hour[lineCounter] = current_hr; #record the current hour
                #         temp_min[lineCounter] = current_min; #record the current min
                #         temp_sec[lineCounter] = current_sec; #record the current sec
                #         #record floats
                #         temp_lat[lineCounter] = current_lat; #record the current lat
                #         temp_long[lineCounter] = current_long; #record the current long
                #         temp_elev[lineCounter] = current_elev; #record the current elevation
                #         temp_vTEC[lineCounter] = current_vTEC; #record the current vTEC
                #         temp_sTEC[lineCounter] = current_sTEC; #record the current sTEC
                #         #TEC err recording is currently disabled b/c it's not in this dataset
                        
                #         lineCounter += 1; #increment the line counter
                #     #END IF
                # #END FOR k
                
                # if( num_dataEntries != lineCounter ):
                #     print("\n==============ERROR==============");
                #     print("Less data than expected was found when reading the data file "+fileNames_Good[j]+". Fix this issue.");
                #     sys.crash();
                # #END IF
                
                # Copy in data into the data columns
                dict_LISN['sat'].append(temp_sat); #copy the data over
                dict_LISN['#'].append(temp_numIndx);
                dict_LISN['day num J2000'].append(temp_dayNumJ2000);
                dict_LISN['hour'].append(temp_hour);
                dict_LISN['min'].append(temp_min);
                dict_LISN['sec'].append(temp_sec);
                dict_LISN['LT'].append(temp_LT);
                dict_LISN['lat'].append(temp_lat);
                dict_LISN['long'].append(temp_long);
                dict_LISN['elev'].append(temp_elev);
                dict_LISN['az'].append(temp_az);
                dict_LISN['sTEC'].append(temp_sTEC);
                # dict_LISN['sTECerr'].append(np.zeros( num_dataEntries,dtype=dataType_str)*np.nan); #preallocate as nan
                dict_LISN['vTEC'].append(temp_vTEC);
                dict_LISN['site'].append(temp_site);
                
                # Copy data into the site info dict
                dict_LISN['site info'][filez_site[filez_currIndx[j]].decode("utf-8")] = {
                    'site location':header_siteLocation,
                    'lat':header_lat,
                    'long':header_long,
                    'alt':header_alt, #km
                    'unknown 4th param':header_4thParam, #not sure what these are!
                    'unknown 5th param':header_5thParam, #not sure what these are!
                    'receiver brand':header_recBrand,
                    'file type':header_fileType,
                    };
                            
                if( (FLG_verbose == True) & (np.any(j == estimatedUpdates)) ):
                    #write an update approximately every 10%
                    sys.stdout.write("\rParsed "+str(j+1)+"/"+str(dataNum)+" ("+str(np.round((j+1)/dataNum*100))+"% finished) ETA is "+str(np.round( ((time.time() - tic)/60)/((j+1)/dataNum) - ((time.time() - tic)/60),2))+" min          " );
                    sys.stdout.flush();
                #END IF
            #END FOR j
            #Collapse lists into a vector
            keyz = list(dict_LISN.keys()); #get the current keys
            for j in range(0,len(keyz)):
                if( ((np.isscalar(dict_LISN[keyz[j]]) == False) & (type(dict_LISN[keyz[j]]) is not dict) ) & (type(dict_LISN[keyz[j]]) is not dict) ):
                    dict_LISN[keyz[j]] = np.hstack(dict_LISN[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
                #END IF
            #END FOR j
            #Make identical data
            dict_LISN['year'] = np.tile(np.int16(filez_date[filez_uniqueIndx[i],0]), dict_LISN['vTEC'].size); #same year
            dict_LISN['day num'] = np.tile(np.int16(filez_dateNum[filez_uniqueIndx[i],1]), dict_LISN['vTEC'].size); #same day 
            dict_LISN['sat type'] = np.tile(np.string_('G'), dict_LISN['vTEC'].size); #assume all GPS satellite type
            dict_LISN['pierce alt'] = pierceAlt_LISN; #save as a scalar
            
            del temp_sat, temp_numIndx, temp_dayNumJ2000, temp_hour, temp_min, temp_sec, temp_LT, \
                temp_lat, temp_long, temp_elev, temp_az, temp_sTEC, temp_vTEC, temp_site; #clean mem
            
            #-----make sure hour 24 or min 60 doesn't exist (they hsould be rolled over) - found issue via day 129 data showing up b/c it was recorded as d128/h24/m60...------       
            #---CHECK SECONDS---
            kj = dict_LISN['sec'] >= 60; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main minute variable where seconds are 60 or more
                dict_LISN['min'][kj] += 1; #min, increment time by 1
                dict_LISN['sec'][kj] -= 60; #sec, remove 60 time units from the time keeping
                if( np.any(dict_LISN['sec'] >= 60) ):
                    print('ERROR in LISN_h5_write: TIME KEEPING IS REALLY BAD, 60+ SECOND TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ SECOND TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---CHECK MINUTES---
            kj = dict_LISN['min'] >= 60; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main hour variable where minutes are 60 or more
                dict_LISN['hour'][kj] += 1; #hour, increment time by 1
                dict_LISN['min'][kj] -= 60; #min, remove 60 time units from the time keeping
                if( np.any(dict_LISN['min'] >= 60) ):
                    print('ERROR in LISN_h5_write: TIME KEEPING IS REALLY BAD, 60+ MINUTE TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 120+ MIN TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---CHECK HOURS---
            kj = dict_LISN['hour'] >= 24; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main day variable where hours are 24 or more
                dict_LISN['day num'][kj] += 1; #day number, increment time by 1
                dict_LISN['hour'][kj] -= 24; #hour, remove 24 time units from the time keeping
                if( np.any(dict_LISN['hour'] >= 24) ):
                    print('ERROR in LISN_h5_write: TIME KEEPING IS REALLY BAD, 24+ HOUR TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 48+ HOUR TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---CHECK DAYS---
            #Deal with day limit is based on leap year or not
            dayLim = np.ones(dict_LISN['day num'].shape,dtype=np.int16)*365; #day number, get the day number limits as 365
            #adjust leap years to 366
            leapYears = (np.mod(dict_LISN['year'],4) == 0) & (np.mod(dict_LISN['year'],100) != 0) & (np.mod(dict_LISN['year'],400) == 0)
            dayLim[leapYears] += 1; #day number, increment time by 1 for leap years
            kj = dict_LISN['day num'] >= dayLim; #find incorrect time keeping
            if( np.sum(kj) > 0 ):
                #increment the main year variable where day number is equal to the day number limit or higher than it
                dict_LISN['year'][kj] += 1; #year, increment time by 1
                dict_LISN['day num'][kj] -= dayLim; #hour, remove the day number limit time units from the time keeping
                if( np.any(dict_LISN['day num'] >= dayLim) ):
                    print('ERROR in LISN_h5_write: TIME KEEPING IS REALLY BAD, 365/366 DAY TIMES WERE REMOVED, YET THEY\'RE STILL THERE.\nTHAT MEANS THAT THERE ARE 365*2/366*2 DAY TIME STAMPS');
                    sys.crash(); #call a crash (crash doesn't exist so it crashes)
                #END IF
            #END IF
            #---luckily years have no limits (that we know of??)---
            #---DELETE data not on the right day (could have incrememted past the right day due to stuff and things like a 24 hour/60 minute time on the current day---
            keyz = list(dict_LISN.keys()); #get the keys
            LISN_fileData_logical_onDay = (dict_LISN['day num'] == filez_dateNum[filez_uniqueIndx[i],1]) & (dict_LISN['year'] == filez_dateNum[filez_uniqueIndx[i],0]); #find when the day reported is the day we want [and the year we want]
            if( LISN_fileData_logical_onDay.size != LISN_fileData_logical_onDay.sum() ):
                print("\n==============~Warning~==============");
                print('Not all data on date '+str(filez_dateNum[filez_uniqueIndx[i],0])+'/'+str(filez_dateNum[filez_uniqueIndx[i],1])+' ('+str(filez_date[filez_uniqueIndx[i],0])+'/'+str(filez_date[filez_uniqueIndx[i],1])+'/'+str(filez_date[filez_uniqueIndx[i],2])+') is on that day - timestamps must be mislabeled within the file? Deleting rogue timestamps.');
                for j in range(0,len(keyz)):
                    if( (np.isscalar(dict_LISN[keyz[j]]) == False) & (type(dict_LISN[keyz[j]]) is not dict)  ):
                        dict_LISN[keyz[j]] = dict_LISN[keyz[j]][LISN_fileData_logical_onDay]; #keep only the good stuff
                    #END IF
                #END FOR j
            #END IF
            #---UPDATE int date variable---
            # dict_LISN['time'] = np.int64(dict_LISN['hour'])*3600 + np.int64(dict_LISN['min'])*60 + np.int64(dict_LISN['sec']); #sec, calculate hour/min/sec into days and add to the current day np.int64(dict_LISN['day num'])*86400 + 
            
            #make sure pierce altitude is consistent
            if( np.isscalar(dict_LISN['pierce alt']) == 0 ):
                if( np.all( dict_LISN['pierce alt'] == dict_LISN['pierce alt'][0]) ):
                    dict_LISN['pierce alt'] = dict_LISN['pierce alt'][0]; #record the attribute
                else:
                    print("\n==============~Warning~==============");
                    print("Not all pierce point altitudes reported are the same - I don't have anything to deal with this. Continuing. They should print below:");
                    print("{}".format(dict_LISN['pierce alt']));
                    dict_LISN['pierce alt'] = dict_LISN['pierce alt'][0]; #record the attribute
                #END IF
            #END IF
            if( FLG_verbose == True ):
                print('\nDone reading data from original files.');
            #END IF
            
            keyz = list(dict_LISN.keys()); #get the keys again
            #--- Get chunk size (it's the size of one of the vectors, since we read in entire vectors at a time) ---
            for j in range(0,len(keyz)):
                if( (np.isscalar(dict_LISN[keyz[j]]) == False) & (type(dict_LISN[keyz[j]]) is not dict)  ):
                    h5pyChunkShape = tuple(np.asarray(dict_LISN[keyz[j]].shape)//10); #get the chunk size of a vector (10 chunks to load from memory, good middle ground)
                    # h5pyChunkShape = dict_LISN[keyz[j]].shape; #get the chunk size of a vector (1 chunk to load into memory, does not seem to lead to any major speed up)
                    break;
                #END IF
            #END FOR j
            
            #--- Write in the unfiltered data in a new faster format ---
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(filez_outputName, 'w', rdcc_nbytes =500*1024*1024) as file_LISN:               
                for j in range(0,len(keyz)):
                    if( np.isscalar(dict_LISN[keyz[j]]) == False ):
                        if( (type(dict_LISN[keyz[j]]) == list) | (type(dict_LISN[keyz[j]]) == np.ndarray)  ):
                            file_LISN.create_dataset(keyz[j], data=dict_LISN[keyz[j]], chunks=h5pyChunkShape, compression='gzip'); #write that data , compression="gzip"
                        else:
                            file_LISN_sub = file_LISN.create_group(keyz[j]); #create a handle for a subgroup
                            sub_keyz = list(dict_LISN[keyz[j]].keys()); #keys to the dict
                            for k in range(0,len(sub_keyz)):
                                if( type(dict_LISN[keyz[j]][sub_keyz[k]]) is dict ):
                                   file_LISN_subSub = file_LISN_sub.create_group(sub_keyz[k]); #create a handle for a subgroup
                                   subSub_keyz = list(dict_LISN[keyz[j]][sub_keyz[k]].keys()); #keys to the dict
                                   for op in range(0,len(subSub_keyz)):
                                       if( np.isscalar(dict_LISN[keyz[j]][sub_keyz[k]][subSub_keyz[op]]) == False ):
                                            file_LISN_subSub.create_dataset(subSub_keyz[op], data=dict_LISN[keyz[j]][sub_keyz[k]][subSub_keyz[op]], chunks=tuple(np.asarray(dict_LISN[keyz[j]][sub_keyz[k]][subSub_keyz[op]].shape)//10), compression="gzip"); #write that data
                                           # file_LISN_subSub.create_dataset(subSub_keyz[op], data=dict_LISN[keyz[j]][sub_keyz[k]][subSub_keyz[op]], chunks=dict_LISN[keyz[j]][sub_keyz[k]][subSub_keyz[op]].shape, compression="gzip"); #write that data
                                       else:
                                           file_LISN_subSub.attrs[subSub_keyz[op]] = dict_LISN[keyz[j]][sub_keyz[k]][subSub_keyz[op]]; #save the attribute
                                       #END IF
                                   #END FOR op
                                else:
                                    if( np.isscalar(dict_LISN[keyz[j]][sub_keyz[k]]) == False ):
                                        file_LISN_sub.create_dataset(sub_keyz[k], data=dict_LISN[keyz[j]][sub_keyz[k]], chunks=dict_LISN[keyz[j]][sub_keyz[k]].shape//10, compression="gzip"); #write that data
                                    else:
                                        file_LISN_sub.attrs[sub_keyz[k]] = dict_LISN[keyz[j]][sub_keyz[k]]; #save the attribute
                                    #END IF
                                #END IF
                            #END FOR k
                        #END IF
                    else:
                        #if size 1, add it as an attribute
                        file_LISN.attrs[keyz[j]] = dict_LISN[keyz[j]]; #save the attribute
                    #END IF
                    if( FLG_verbose == True ):
                        sys.stdout.write("\rWriting {} to file & {} min | {} out of {}\t\t\t\t\t".format(keyz[j],np.round((time.time()-tic)/60,2),j+1,len(keyz)));
                        sys.stdout.flush();
                    #END IF
                #END FOR j
                if( FLG_verbose == True ):
                    print('\nDone writing LISN HDF5 file!');
                #END IF
            #END WITH
            
            del dict_LISN; #clean the memory
            
            if( FLG_verbose == True ):
                toc = time.time() - tic; #for time testing
                current_timePerVect = np.round(toc/dataNum,2); #update the time estiamte (includes writing files and such)
                print("\nTime to convert: {} min\n".format(np.round(toc/60,2))); #extra space at end  
            #END IF
        #END IF
    #END FOR i
    

    #--- Delete original files if req'd ---
    if( FLG_deleteOrig == True ): # Only delete if flag is on
        if( FLG_verbose == True ):
            print('Deleting '+str(filez.size)+' original LISN .dat files.');
        #END IF
        for i in range(0,filez.size ):
            os.remove(filez[i]); # Delete the files
        #END FOR i
    #END IF  
#END DEF


def LISN_h5_read(dateRange=None, path_read=None, FLG_recursive=False, fileNameConvention='LISN_{YR}_{DN}.h5', FLG_combine=True, FLG_ignoreSubGroups=False, file_maxRetry=3, FLG_verbose=True):
    #--- Define constants ---
    version_LISN = LISN_h5_version(); #get current version (defined by function at top)
    fileNameEnding = fileNameConvention[fileNameConvention.rfind('.'):]; #get file name ending from fileNameConvention
    
    
    #--- Align current working directory ---
    if( path_read != None):
        os.chdir(path_read); # Update CWD to provided dir
    else:
        # If fldr not supplied, force dir to be where LISN_h5.py is
        os.chdir(pathlib.Path('.').resolve()); # Update CWD to where file is
    #END IF
    path_read = os.getcwd(); # Record CWD
    if( FLG_verbose == True ):
        print('--- LISN_h5_read ---');
        print('Current Working Directory: '+path_read); # Report CWD
    #END IF
    
    
    #--- Analyze dateRange provided ---
    if( dateRange != None ):
        if( (type(dateRange) == list) | (type(dateRange) == tuple) ):
            dateRange = np.array(dateRange); #convert to numpy array
        #END IF
        if( (dateRange.shape[0] == 2) & (dateRange.shape[1] == 2) & (dateRange.ndim == 2) ):
            if( FLG_verbose == True ):
                print('Start date (YR/DN, inclusive):\t'+str(dateRange[0,:])+'\nEnd date (YR/DN, inclusive):\t'+str(dateRange[-1,:])+'\n'); # Report date range
            #END IF
            _, dateRange_full = subfun_dateORdayNum_to_fullRange(dateRange); #get full date range
        elif( (dateRange.shape[0] == 2) & (dateRange.shape[1] == 3) & (dateRange.ndim == 2) ):
            if( FLG_verbose == True ):
                print('Start date (YR/MO/DY, inclusive):\t'+str(dateRange[0,:])+'\nEnd date (YR/MO/DY, inclusive):\t'+str(dateRange[-1,:])+'\n'); # Report date range
            #END IF
            _, dateRange_full = subfun_dateORdayNum_to_fullRange(dateRange); #convert to dayNum from date & get full date range
        elif( (dateRange.shape[0] == 3) & (dateRange.shape[1] == 2) & (dateRange.ndim == 2) ):
            if( FLG_verbose == True ):
                print('Start date (YR/MO/DY, inclusive):\t'+str(dateRange[:,0])+'\nEnd date (YR/MO/DY, inclusive):\t'+str(dateRange[:,-1])+'\n'); # Report date range
            #END IF
            _, dateRange_full = subfun_dateORdayNum_to_fullRange(dateRange.T); #convert to dayNum from date & get full date range
        else:
            print('\nERROR in LISN_h5_read: dateRange provided is not a date range that is recognized. Printing provided dateRange & crashing:\n'+str(dateRange));
            sys.crash(); #crash
        #END IF
    else:
        if( FLG_verbose == True ):
            print('No date range provided, reading & importing in all LISN files that can be found in the current working directory.');
        #END IF
    #END IF
    
    
    #--- Discover files ---
    # Only have a few static things about the files to work with - they end with .dat, have 4 letter/# site names, and have some #s in there
    if( FLG_recursive == False ):
        filez = [os.path.join(path_read,file) for file in os.listdir(path_read) if file.lower().endswith(fileNameEnding)]; # Get list of .h5 files (NOT recursive), inspired by https://stackoverflow.com/a/3964696/2403531
    else:
        filez = [os.path.join(root, file) for root, dirs, files in os.walk(path_read) for file in files if file.lower().endswith(fileNameEnding)]; # Get list of .h5 files (recursive!!), inspired by https://stackoverflow.com/a/3964690/2403531
    #END IF
    
    # Remove files that aren't LISN data files but do have a .h5 ending, since that is a common ending
    filez_remove = set(); # Create empty set
    fileNameConvention_regex = fileNameConvention.replace('{YR}','(\d+)'); #replace with regex for digits
    fileNameConvention_regex = fileNameConvention_regex.replace('{Y}','(\d+)'); #replace with regex for digits
    fileNameConvention_regex = fileNameConvention_regex.replace('{DN}','(\d+)'); #replace with regex for digits
    fileNameConvention_regex = fileNameConvention_regex.replace('{MO}','(\d+)'); #replace with regex for digits
    fileNameConvention_regex = fileNameConvention_regex.replace('{M}','(\d+)'); #replace with regex for digits
    fileNameConvention_regex = fileNameConvention_regex.replace('{DY}','(\d+)'); #replace with regex for digits
    fileNameConvention_regex = fileNameConvention_regex.replace('{D}','(\d+)'); #replace with regex for digits
    fileNameConvention_regexMatchOrder = []; #prep
    fileNameConvention_regexMatchOrderName = []; #prep
    fileNameConvention_Y = False;
    fileNameConvention_DNmode = False;
    fileNameConvention_MDmode = False;
    fileNameConvention_MDmode_M = False; #these help activate MDmode
    fileNameConvention_MDmode_D = False; #these help activate MDmode
    if( fileNameConvention.find('{YR}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{YR}'));
        fileNameConvention_regexMatchOrderName.append('{YR}');
        fileNameConvention_Y = True; # make sure year is here
    #END IF
    if( fileNameConvention.find('{Y}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{Y}'));
        fileNameConvention_regexMatchOrderName.append('{Y}');
        fileNameConvention_Y = True; # make sure year is here
    #END IF
    if( fileNameConvention.find('{DN}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{DN}'));
        fileNameConvention_regexMatchOrderName.append('{DN}');
        fileNameConvention_DNmode = True; # turn on DN mode
    #END IF
    if( fileNameConvention.find('{MO}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{MO}'));
        fileNameConvention_regexMatchOrderName.append('{MO}');
        fileNameConvention_MDmode_M = True;
    #END IF
    if( fileNameConvention.find('{M}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{M}'));
        fileNameConvention_regexMatchOrderName.append('{M}');
        fileNameConvention_MDmode_M = True;
    #END IF
    if( fileNameConvention.find('{DY}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{DY}'));
        fileNameConvention_regexMatchOrderName.append('{DY}');
        fileNameConvention_MDmode_D = True;
    #END IF
    if( fileNameConvention.find('{D}') >= 0 ):
        fileNameConvention_regexMatchOrder.append(fileNameConvention.find('{D}'));
        fileNameConvention_regexMatchOrderName.append('{D}');
        fileNameConvention_MDmode_D = True;
    #END IF
    fileNameConvention_regexMatchIndx = np.argsort(fileNameConvention_regexMatchOrder); #get indexes for the order of the {}'s
    
    if( (fileNameConvention_MDmode_M == True) & (fileNameConvention_MDmode_D == True) ):
        fileNameConvention_MDmode = True;
        filez_date = np.empty((len(filez),3), dtype=np.int16); # Preallocate [yr, mo, dy]
    #END IF
    if( fileNameConvention_DNmode == True ):
        filez_dateNum = np.empty((len(filez),2), dtype=np.int16); # Preallocate [yr, dayNum]
    #END IF
    
    for i in range(0,len(filez)):
        filez_match = re.search(fileNameConvention_regex,filez[i][filez[i].rfind(os.path.sep)+len(os.path.sep):]); #regex search for the values we're looking for
        if( filez_match is None ):
            filez_remove.add(i); # Add the index to the set if where the date is isn't only numbers
        else:
            for j in range(0,fileNameConvention_regexMatchIndx.size):
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{YR}' ):
                    if( fileNameConvention_DNmode == True ):
                        filez_dateNum[i,0] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                    if( fileNameConvention_MDmode == True ):
                        filez_date[i,0] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{Y}' ):
                    if( fileNameConvention_DNmode == True ):
                        filez_dateNum[i,0] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                    if( fileNameConvention_MDmode == True ):
                        filez_date[i,0] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{DN}' ):
                    if( fileNameConvention_DNmode == True ):
                        filez_dateNum[i,1] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{MO}' ):
                    if( fileNameConvention_MDmode == True ):
                        filez_date[i,1] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{M}' ):
                    if( fileNameConvention_MDmode == True ):
                        filez_date[i,1] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{DY}' ):
                    if( fileNameConvention_MDmode == True ):
                        filez_date[i,2] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
                if( fileNameConvention_regexMatchOrderName[fileNameConvention_regexMatchIndx[j]] == '{D}' ):
                    if( fileNameConvention_MDmode == True ):
                        filez_date[i,2] = np.int16(filez_match.groups()[j]); #insert
                    #END IF
                #END IF
            #END FOR j
        #END IF        
    #END FOR i
    filez = [i for j, i in enumerate(filez) if j not in filez_remove]; # Remove filez that don't seem to be LISN data files
    if( ((fileNameConvention_DNmode == True) | (fileNameConvention_MDmode == True)) & (fileNameConvention_Y == True) ):
        if( (fileNameConvention_MDmode == True) & (len(filez_remove) != 0) ):
            filez_date = np.delete(filez_date,np.array(list(filez_remove)),axis=0); # Also remove these since they were calc'd before file removal
        #END IF
        if( fileNameConvention_DNmode == True ):
            if( len(filez_remove) != 0 ):
                filez_dateNum = np.delete(filez_dateNum,np.array(list(filez_remove)),axis=0); # Also remove these since they were calc'd before file removal
            #END IF
        else:
            filez_dateNum = subfun_dayNum_to_date(filez_date); # Convert to month/day format
        #END IF        
        if( fileNameConvention.find('{Y}') >= 0 ):
            filez_dateNum[:,0] += 2000; # Make year full (assume no 90's data)
        #END IF
    else:
        for i in range(0,len(filez)):
            file_retryNum = 1; #prep the retry num
            FLG_gotRead = False; #lets while exit
            while( ((file_maxRetry >= file_retryNum) | (file_maxRetry == -1)) & (FLG_gotRead == False) ): #deal with intermittant file read issues I had on an old comp, it tries
                try:
                    with h5py.File(filez[i], 'r') as LISN_file:
                        temp_year = LISN_file.get('year')[()]; #get that dataset out
                        temp_dayNum = LISN_file.get('day num')[()]; #get that dataset out
                    #END WITH
                    if( np.all(np.unique(temp_year) == temp_year[0]) ):
                        filez_dateNum[i,0] = temp_year[0]; #record year
                    else:
                        print('ERROR in LISN_h5_read: Year is not consistent within "'+filez[i]+'". Look at this and fix the file, crashing.');
                        sys.crash();
                    #END IF
                    if( np.all(np.unique(temp_dayNum) == temp_dayNum[0]) ):
                        filez_dateNum[i,1] = temp_dayNum[0]; #record dayNum
                    else:
                        print('ERROR in LISN_h5_read: day num is not consistent within "'+filez[i]+'". Look at this and fix the file, crashing.');
                        sys.crash();
                    #END IF
                    del temp_year, temp_dayNum; #clean mem
                except OSError as errorText:
                    print('Warning in LISN_h5_read: '+filez[i] +\
                          ' had an OSError when reading data. Try #'+str(file_retryNum)+'/'+str(file_maxRetry)+'. Error text follows:');
                    print(str(errorText));
                    time.sleep(0.1); #wait a tiny lil bit just in case
                    file_retryNum += 1; #increment try
                #END TRY
            #END WHILE
            if( FLG_gotRead == False ):
                print("\n==============ERROR==============");
                print(filez[i] +\
                      ' failed to read. Renaming the file to "'+filez[i]+'_corrupted" and crashing. Rerun code to generate new file.');
                os.rename(filez[i], \
                          filez[i]+'_corrupted'); #rename
                print('NOTE this may be an OS read error that can be fixed with a system restart. Not sure why it happens but you\'ll need to rename it back from _corrupted and restart your comp. Restarting python instance doesn\'t seem to cut it.');
                sys.crash(); #donezo
            #END IF
        #END FOR i
    #END IF
        
    
    #--- Catch errors ---
    #Catch no date range provided
    if( dateRange == None ):
        dateRange_full = filez_dateNum; #copy over file list as date range
    #END IF
    
    # Confirm data availability for the date range requested
    dataAvail = isin_row(dateRange_full, filez_dateNum); #check row-wise if data is avail
    if( np.any(dataAvail == False) ):
        print('\nWARNING in LISN_h5_read: No .h5 file for dates:\n'+str(dateRange_full[~dataAvail,:])+'\nCalling LISN_h5_write to hopefully fill in missing data.\n');
        LISN_h5_write(path_data=path_read, FLG_recursive=True, fileNameConvention=fileNameConvention, FLG_verbose=FLG_verbose);
        dateRange_full_dateFormat = subfun_dayNum_to_date(dateRange_full); #convert just in case it is needed
        dataAvail_missing = np.where(~dataAvail)[0];
        for i in range(0,dataAvail_missing.size):
            dataAvail_missing_fileName = fileNameConvention.replace('{YR}', str(dateRange_full[dataAvail_missing[i],0]) ); #replace any {YR} with the current year
            dataAvail_missing_fileName = dataAvail_missing_fileName.replace('{Y}', str(dateRange_full[dataAvail_missing[i],0])[2:] ); #replace any {Y} with the current year's last 2 digits
            dataAvail_missing_fileName = dataAvail_missing_fileName.replace('{DN}', str(dateRange_full[dataAvail_missing[i],1]).zfill(3) ); #replace any {DN} with the current day number (padded w/ 0's so always 3 long)
            dataAvail_missing_fileName = dataAvail_missing_fileName.replace('{MO}', str(dateRange_full_dateFormat[dataAvail_missing[i],1]).zfill(2) ); #replace any {MO} with the current month (padded w/ 0's so always 2 long)
            dataAvail_missing_fileName = dataAvail_missing_fileName.replace('{M}', str(dateRange_full_dateFormat[dataAvail_missing[i],1]) ); #replace any {M} with the current month
            dataAvail_missing_fileName = dataAvail_missing_fileName.replace('{DY}', str(dateRange_full_dateFormat[dataAvail_missing[i],2]).zfill(2) ); #replace any {DY} with current day (padded w/ 0's so always 2 long)
            dataAvail_missing_fileName = dataAvail_missing_fileName.replace('{D}', str(dateRange_full_dateFormat[dataAvail_missing[i],2]) ); #replace any {D} with current day
            dataAvail_missing_fileName = os.path.join(path_read,dataAvail_missing_fileName); #tack on CWD
            if( os.path.isfile(dataAvail_missing_fileName) ):
                dataAvail[dataAvail_missing[i]] = True; #set data to available
                filez.append(dataAvail_missing_fileName); #append the data on as available
                filez_dateNum = np.vstack((filez_dateNum,np.int16(dateRange_full[dataAvail_missing[i],:]))); #append the date on as available
            else:
                print('\nERROR in LISN_h5_read: "'+dataAvail_missing_fileName+'" does not exist, even after calling LISN_h5_write. That means there is no original data for that day. Subsequent missing dates may also be lacking data but this code will crash before checking that.');
                sys.crash(); #crash
            #END IF
        #END FOR i
        print('All data filled in and ready for use.');
    #END IF
    # Confirm data version
    FLG_anyBad = False;
    for i in range(0,len(filez)):
        file_retryNum = 1; #prep the retry num
        FLG_gotRead = False; #lets while exit
        while( ((file_maxRetry >= file_retryNum) | (file_maxRetry == -1)) & (FLG_gotRead == False) ): #deal with intermittant file read issues I had on an old comp, it tries
            try:
                with h5py.File(filez[i], 'r') as h5_LISN:
                    version_LISN_fromFile = h5_LISN.attrs['version']; #read the version
                #END WITH
                if( version_LISN_fromFile == version_LISN ):
                    pass;
                elif( version_LISN_fromFile < version_LISN ):
                    FLG_anyBad = True; # version bad
                else:
                    print('WARNING in LISN_h5_read: LISN file version is newer than this code\'s version. Update this code.');
                    pass;
                #END IF
                FLG_gotRead = True;
            except OSError as errorText:
                print('Warning in LISN_h5_read: '+filez[i] +\
                      ' had an OSError when reading data. Try #'+str(file_retryNum)+'/'+str(file_maxRetry)+'. Error text follows:');
                print(str(errorText));
                time.sleep(0.1); #wait a tiny lil bit just in case
                file_retryNum += 1; #increment try
            #END TRY
        #END WHILE
        if( FLG_gotRead == False ):
            print("\n==============ERROR==============");
            print(filez[i] +\
                  ' failed to read. Renaming the file to "'+filez[i]+'_corrupted" and crashing. Rerun code to generate new file.');
            os.rename(filez[i], \
                      filez[i]+'_corrupted'); #rename
            print('NOTE this may be an OS read error that can be fixed with a system restart. Not sure why it happens but you\'ll need to rename it back from _corrupted and restart your comp. Restarting python instance doesn\'t seem to cut it.');
            # sys.crash(); #donezo
            FLG_anyBad = True; # Note it bad
        #END IF
    #END FOR i
    if( FLG_anyBad == True ):
        print('\nWARNING in LISN_h5_read: Old versions of .h5 files detected. Calling LISN_h5_write to hopefully fill in old versions of the data.');
        LISN_h5_write(path_data=path_read, FLG_recursive=True, fileNameConvention=fileNameConvention, FLG_verbose=FLG_verbose);
        print('LISN_h5_write complete, hopefully all .h5 files are at current version.');
    #END IF
    FLG_anyBad = False;
    for i in range(0,len(filez)):
        file_retryNum = 1; #prep the retry num
        FLG_gotRead = False; #lets while exit
        while( ((file_maxRetry >= file_retryNum) | (file_maxRetry == -1)) & (FLG_gotRead == False) ): #deal with intermittant file read issues I had on an old comp, it tries
            try:
                with h5py.File(filez[i], 'r') as h5_LISN:
                    version_LISN_fromFile = h5_LISN.attrs['version']; #read the version
                #END WITH
                if( version_LISN_fromFile == version_LISN ):
                    pass;
                elif( version_LISN_fromFile < version_LISN ):
                    FLG_anyBad = True; # version bad
                    print('File "'+filez[i]+'" has bad version '+str(version_LISN_fromFile)+' that is lower than current version '+str(version_LISN)+'. This will cause a controlled crash shortly.');
                else:
                    print('WARNING in LISN_h5_read: LISN file version is newer than this code\'s version. Update this code.');
                    pass;
                #END IF
                FLG_gotRead = True;
            except OSError as errorText:
                print('Warning in LISN_h5_read: '+filez[i] +\
                      ' had an OSError when reading data. Try #'+str(file_retryNum)+'/'+str(file_maxRetry)+'. Error text follows:');
                print(str(errorText));
                time.sleep(0.1); #wait a tiny lil bit just in case
                file_retryNum += 1; #increment try
            #END TRY
        #END WHILE
        if( FLG_gotRead == False ):
            print("\n==============ERROR==============");
            print(filez[i] +\
                  ' failed to read. Renaming the file to "'+filez[i]+'_corrupted" and crashing. Rerun code to generate new file.');
            os.rename(filez[i], \
                      filez[i]+'_corrupted'); #rename
            print('NOTE this may be an OS read error that can be fixed with a system restart. Not sure why it happens but you\'ll need to rename it back from _corrupted and restart your comp. Restarting python instance doesn\'t seem to cut it.');
            # sys.crash(); #donezo
            FLG_anyBad = True; # Note it bad
        #END IF
    #END FOR i
    if( FLG_anyBad == True ):
        print('\nERROR in LISN_h5_read: Old versions of .h5 files still detected. Crashing b/c original data is lacking for files listed above.');
        sys.crash();
    #END IF
    
    
    def HDF5_recursiveRead(HDF5_handle):
        dictHolder = {}; #prep a dict to hold things
        keyz = list(HDF5_handle.keys()); #get the saved keys
        #--- Now import rest of vars and use mask on them ---
        for j in range(0,len(keyz)):
            if( isinstance(HDF5_handle[keyz[j]], h5py.Dataset) ):
                dictHolder[keyz[j]] = HDF5_handle.get(keyz[j])[()]; #get that dataset out
            elif( isinstance(HDF5_handle[keyz[j]], h5py.Group) ):
                dictHolder[keyz[j]] = {}; #prep a dict
                dictHolder[keyz[j]] = HDF5_recursiveRead(HDF5_handle.get(keyz[j])); #recursively call this function
            else:
                print('\nERROR in HDF5_recursiveRead: "'+str(HDF5_handle)+'" has unsupported HDF5 data type: "'+str(HDF5_handle.get(keyz[j]))+'".');
                sys.crash();
            #END IF
        #END FOR j
        #--- Read the attributes in ---
        keyz_attrs = list(HDF5_handle.attrs.keys()); #get the attribute keys
        for j in range(0,len(keyz_attrs)):
            dictHolder[keyz_attrs[j]] = HDF5_handle.attrs[keyz_attrs[j]]; #get that attribute out
        #END FOR j
        return dictHolder
    #END DEF
    
    def HDF5_readNoSubGroups(HDF5_handle):
        dictHolder = {}; #prep a dict to hold things
        keyz = list(HDF5_handle.keys()); #get the saved keys
        #--- Now import rest of vars and use mask on them ---
        for j in range(0,len(keyz)):
            if( isinstance(HDF5_handle[keyz[j]], h5py.Dataset) ):
                dictHolder[keyz[j]] = HDF5_handle.get(keyz[j])[()]; #get that dataset out
            elif( isinstance(HDF5_handle[keyz[j]], h5py.Group) ):
                pass; #no sub groups speeds up reading
            else:
                print('\nERROR in HDF5_readNoSubGroups: "'+str(HDF5_handle)+'" has unsupported HDF5 data type: "'+str(HDF5_handle.get(keyz[j]))+'".');
                sys.crash();
            #END IF
        #END FOR j
        #--- Read the attributes in ---
        keyz_attrs = list(HDF5_handle.attrs.keys()); #get the attribute keys
        for j in range(0,len(keyz_attrs)):
            dictHolder[keyz_attrs[j]] = HDF5_handle.attrs[keyz_attrs[j]]; #get that attribute out
        #END FOR j
        return dictHolder
    #END DEF
    
    #--- Read files ---
    LISN_read_temp = [None for i in range(0,dateRange_full.shape[0])]; #preallocate dicts
    for i in range(0,dateRange_full.shape[0]):
        filez_currIndx = np.where((dateRange_full[i,0] == filez_dateNum[:,0]) & (dateRange_full[i,1] == filez_dateNum[:,1]))[0].item(); #get index where file is at
        file_retryNum = 1; #prep the retry num
        FLG_gotRead = False; #lets while exit
        while( ((file_maxRetry >= file_retryNum) | (file_maxRetry == -1)) & (FLG_gotRead == False) ): #deal with intermittant file read issues I had on an old comp, it tries
            try:
                if( FLG_ignoreSubGroups == False ):
                    with h5py.File(filez[filez_currIndx], 'r') as LISN_file:
                        LISN_read_temp[i] = HDF5_recursiveRead(LISN_file); #call the HDF5 recursive reader to read the file and save it as the same dict that went in to make it
                    #END WITH
                else:
                    with h5py.File(filez[filez_currIndx], 'r') as LISN_file:
                        LISN_read_temp[i] = HDF5_readNoSubGroups(LISN_file); #call the HDF5 recursive reader to read the file and save it as the same dict that went in to make it
                    #END WITH
                #END IF
                FLG_gotRead = True;
            except OSError as errorText:
                print('Warning in LISN_h5_read: '+filez[i] +\
                      ' had an OSError when reading data. Try #'+str(file_retryNum)+'/'+str(file_maxRetry)+'. Error text follows:');
                print(str(errorText));
                time.sleep(0.1); #wait a tiny lil bit just in case
                file_retryNum += 1; #increment try
            #END TRY
        #END WHILE
        if( FLG_gotRead == False ):
            print("\n==============ERROR==============");
            print(filez[i] +\
                  ' failed to read. Renaming the file to "'+filez[i]+'_corrupted" and crashing. Rerun code to generate new file.');
            os.rename(filez[i], \
                      filez[i]+'_corrupted'); #rename
            print('NOTE this may be an OS read error that can be fixed with a system restart. Not sure why it happens but you\'ll need to rename it back from _corrupted and restart your comp. Restarting python instance doesn\'t seem to cut it.');
            sys.crash(); #donezo
        #END IF
    #END FOR i
    
    if( FLG_combine == True ):
        LISN_read = {}; #prep
        for i in range(0,len(LISN_read_temp)):
            keyz =list(LISN_read.keys()); #get the saved keys
            keyzNew = list(LISN_read_temp[i].keys()); #get the saved keys
            for j in range(0,len(keyzNew)):
                if( strfind(keyz,keyzNew[j],opt=1) > 0 ):
                    LISN_read[keyzNew[j]].append(LISN_read_temp[i][keyzNew[j]]); #tack that dataset on
                    del LISN_read_temp[i][keyzNew[j]]; #save memory
                else:
                    #otherwise it's a new data type to add in
                    LISN_read[keyzNew[j]] = [LISN_read_temp[i][keyzNew[j]]]; #get that dataset out
                    del LISN_read_temp[i][keyzNew[j]]; #save memory
                #END IF
            #END FOR j
        #END FOR i
        del LISN_read_temp; #clear it out completely
        keyz = list(LISN_read.keys()); #get the current keys
        for j in range(0,len(keyz)):
            if( (type(LISN_read[keyz[j]][0]) == np.ndarray) | (type(LISN_read[keyz[j]][0]) == list) ):
                LISN_read[keyz[j]] = np.hstack(LISN_read[keyz[j]]); #convert from list to array (cause list can be stacked fast but is awful for using)
            elif( type(LISN_read[keyz[j]][0]) == dict ):
                LISN_read[keyz[j]] = {dictKeyz:dictValues for dictz in LISN_read[keyz[j]] for dictKeyz,dictValues in dictz.items()}; #keep only unique site entries, combine together, inspired by https://stackoverflow.com/a/31745620/2403531
            elif( np.isscalar(LISN_read[keyz[j]][0]) == True ):
                if( np.all(LISN_read[keyz[j]] == LISN_read[keyz[j]][0]) ):
                    LISN_read[keyz[j]] = LISN_read[keyz[j]][0]; #keep 1st one
                else:
                    print('WARNING in LISN_h5_read: All scalars are not the same for key '+keyz[j]+' but they should be.\nPrinting list of scalars and keeping 1st: '+str(LISN_read[keyz[j]]));
                    LISN_read[keyz[j]] = LISN_read[keyz[j]][0]; #keep 1st one
                #END IF
            #END IF
        #END FOR j
    else:
        LISN_read = LISN_read_temp; #no need to do anything
    #END IF
    
    return LISN_read
#END DEF


#=================== SUPPORTING FUNCTIONS =====================================
#GOAL: Day Number to Date
#RD on 8/23/2018, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
#
#INPUT: dateRange_dayNum as [YR,dayNum] format for as many as you like - must be [#dates , 2] sized
#OUTPUT: date in [YRstart/M/D , YRend/M/D] format, for each date input

def subfun_dayNum_to_date(dateRange):
    #==============Catch input format issues==============
    if( isinstance(dateRange,tuple) | isinstance(dateRange,list) ):
        dateRange = np.array(dateRange); #convert to numpy array
    #END IF
    if( dateRange.ndim == 1 ):
        dateRange = dateRange[..., np.newaxis]; #add on a new axis b/c code demands it
    #END IF
    
    if (len(dateRange[0,:]) != 2) and (len(dateRange[:,0]) == 2):
        #Catch where a [3,arbitrary] was sent but this needs to operate on a [arbitrary,3] sized array
        dateRange = np.swapaxes(dateRange,0,1); #flip axes
    elif (len(dateRange[0,:]) != 2) and (len(dateRange[:,0]) != 2):
        #Catch where something formatted completely wrong was sent
        print("\n==============ERROR==============");
        print("Input size was [{},{}] and it needs to be [arbitrary,2] - exiting dayNum to date fun!\n".format( len(dateRange[0,:]),len(dateRange[:,0]) ) );
        # return "No"; #can I? I will
        import sys
        sys.crash(); #even better
    #END IF
    
    #==============Prep constants and preallocate==============
    monthDays = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]); #preps number of days in a month
    monthDays_Leap = np.array([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]); #preps number of days in a month
    
    
    dateRange_converted = np.zeros( (len(dateRange[:,0]),3),dtype="int16"); #preallocate (named so to make sharing code easier)

    #==============Convert dates given in M/D/YR format to dayNum/YR format==============
    for i in range(0, len(dateRange[:,0]) ):
        #-----Leap Year Detection-----
        if np.mod(dateRange[i,0],4) == 0: #leap year
            #-----Leap Year Skipped Detected - next will be 2100-----
            if (np.mod(dateRange[i,0],100) == 0) and (np.mod(dateRange[i,0],400) != 0):
                #Same alg as no leap year used here
                
                tempDayCntr = dateRange[i,1]; #get the current day number
                tempMonCntr = 0; #counts the months
                while (tempDayCntr > 0) and (tempMonCntr <= 12): #makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
                    tempMonCntr = tempMonCntr + 1; #increment month used
                    tempDayCntr = tempDayCntr - monthDays[tempMonCntr-1]; #days, remove days in month from year
                #END WHILE
                if tempDayCntr > 0: #check for errors
                    print("Day count greater than 0 and days left in year: {} days left. Please check that input of {} year and {} day number.".format(tempDayCntr,dateRange[i,1],dateRange[i,0]));
                    return "No."; #can I? I will
                #END IF
                tempDayCntr = tempDayCntr + monthDays[tempMonCntr-1]; #days, calculate the day in the month we found
            
            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
            else:            
                tempDayCntr = dateRange[i,1]; #get the current day number
                tempMonCntr = 0; #counts the months
                while tempDayCntr > 0 and tempMonCntr <= 12: #makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
                    tempMonCntr = tempMonCntr + 1; #increment month used
                    tempDayCntr = tempDayCntr - monthDays_Leap[tempMonCntr-1]; #days, remove days in month from year
                #END WHILE
                if tempDayCntr > 0: #check for errors
                    print("Day count greater than 0 and days left in year: {} days left. Please check that input of {} year and {} day number.".format(tempDayCntr,dateRange[i,1],dateRange[i,0]));
                    return "No."; #can I? I will
                #END IF
                tempDayCntr = tempDayCntr + monthDays_Leap[tempMonCntr-1]; #days, calculate the day in the month we found
            
            #END IF
        #-----No leap year detected-----
        else: #no leap year if this
            
            tempDayCntr = dateRange[i,1]; #get the current day number
            tempMonCntr = 0; #counts the months
            while (tempDayCntr > 0) and (tempMonCntr <= 12): #makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
                tempMonCntr = tempMonCntr + 1; #increment month used
                tempDayCntr = tempDayCntr - monthDays[tempMonCntr-1]; #days, remove days in month from year
            #END WHILE
            if tempDayCntr > 0: #check for errors
                print("Day count greater than 0 and days left in year: {} days left. Please check that input of {} year and {} day number.".format(tempDayCntr,dateRange[i,1],dateRange[i,0]));
                return "No."; #can I? I will
            #END IF
            tempDayCntr = tempDayCntr + monthDays[tempMonCntr-1]; #days, calculate the day in the month we found
        
        #END IF
        
        dateRange_converted[i,:] = [dateRange[i,0],tempMonCntr,tempDayCntr]; #pack up the date and ship out
    #END FOR  
    
    return(dateRange_converted); #return success
#END DEF

#GOAL: Date to Day Number
#RD on 8/23/2018, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
#
#INPUT: date in [YRstart/M/D , YRend/M/D] format, as many as you like - must be [#dates , 3] sized
#OUTPUT: dateRange_dayNum as [YR/dayNum] format for each date input
#options!: 
#0 = Output as [dayNum/YR] 
#1 = Output as [YR/dayNum] [DEFAULT!]
#2 = Output as just [dayNum]

def subfun_date_to_dayNum(dateRange,options = 1):  
    #==============Catch input format issues==============
    if( isinstance(dateRange,tuple) | isinstance(dateRange,list) ):
        dateRange = np.array(dateRange); #convert to numpy array
    #END IF
    if( dateRange.ndim == 1 ):
        dateRange = dateRange[..., np.newaxis]; #add on a new axis b/c code demands it
    #END IF

    if (len(dateRange[0,:]) != 3) and (len(dateRange[:,0]) == 3):
        #Catch where a [3,arbitrary] was sent but this needs to operate on a [arbitrary,3] sized array
        # print("\n==============~Warning~==============");
        # #print("Input size was [{},{}] and it is assumed that the dimensions are flipped -> adjusting to [{},{}] as [arbitrary,3] format is needed.\n".format( len(dateRange[0,:]),len(dateRange[:,0]),len(dateRange[:,0]),len(dateRange[0,:]) ) ); #.format version
        # print("Input size was ["+str(len(dateRange[0,:]))+","+str(len(dateRange[:,0]))+"] and it is assumed that the dimensions are flipped -> adjusting to ["+str(len(dateRange[:,0]))+","+str(len(dateRange[0,:]))+"] as [arbitrary,3] format is needed.\n"); #non-.format version for nopython calls
        dateRange = np.swapaxes(dateRange,0,1); #flip axes
        
    elif (len(dateRange[0,:]) != 3) and (len(dateRange[:,0]) != 3):
        #Catch where something formatted completely wrong was sent
        print("\n==============ERROR==============");
        # print("Input size was [{},{}] and it needs to be [arbitrary,3] - exiting date to dayNum fun!\n".format( len(dateRange[0,:]),len(dateRange[:,0]) ) ); #.format version
        print("Input size was ["+str(len(dateRange[0,:]))+","+str(len(dateRange[:,0]))+"}] and it needs to be [arbitrary,3] - exiting date to dayNum fun!\n"); #non-.format version for nopython calls
        # return("No"); #can I? I will
        import sys
        sys.crash(); #even better
    #END IF

    #Split code base into the seperate if's for a possible speed-up
    #==============Convert dates given in M/D/YR format to dayNum/YR format==============    
    if( options == 0 ):
    #-----Option 0 returns the default output of [dayNum/Yr] format-----
        
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) , 2 ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) , 2 ) ) ); #preallocate , numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1), : ] = np.array([tempMonth + dateRange_unique[i,2],dateRange_unique[i,0]]); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success
    
    elif( options == 1 ):
    #-----Option 1 returns the output in the [Yr/dayNum] format-----
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) , 2 ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) , 2 ) ) ); #preallocate, numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1), : ] = np.array([dateRange_unique[i,0],tempMonth + dateRange_unique[i,2]]); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success
    
    elif( options == 2 ):
    #-----Option 2 returns the output in the [dayNum] only format-----
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) ) ) ); #preallocate , numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1) ] = np.array( [tempMonth + dateRange_unique[i,2]] ); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success
    
    else:
        #-----Anything else returns the default output [dayNum/Yr] format-----
        print("\n==============~Warning~==============");
        #print("Incorrect option number of {} input. Returning in default form of [dayNum, YR]\nAllowed options: 0 == [dayNum/YR], 1 == [YR/dayNum], 2 == [dayNum]\n".format(options)); #.format version
        print("Incorrect option number of "+str(options)+" input. Returning in default form of [YR, dayNum]\nAllowed options: 0 == [dayNum/YR], 1 == [YR/dayNum], 2 == [dayNum]\n"); #non-.format version for nopython calls
        
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) , 2 ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) , 2 ) ) ); #preallocate , numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1), : ] = np.array([dateRange_unique[i,0],tempMonth + dateRange_unique[i,2]]); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success?
    #END IF
#END DEF

#GOAL: Date OR Day Number to the full range between the range
#RD on 8/23/2018, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
#
#INPUT: date range in [YRstart/M/D , YRend/M/D] format, [2 , 3] sized ~OR~ date range in [YRstart/dayNum , YRend/dayNum] format, [2,2] sized
#OUTPUT: dateRange_dayNum_full and dateRange_dayNum_full as [#datesInRange,3] and [#datesInRange,2], respectively, format

def subfun_dateORdayNum_to_fullRange(dateRange):

    #==============Catch input format issues==============
    if( isinstance(dateRange,tuple) | isinstance(dateRange,list) ):
        dateRange = np.array(dateRange); #convert to numpy array
    #END IF
    if( dateRange.ndim == 1 ):
        dateRange = dateRange[..., np.newaxis]; #add on a new axis b/c code demands it
    #END IF

    if (len(dateRange[0,:]) == 2) and (len(dateRange[:,0]) == 3):
        #Catch where a [3,2] was sent and this needs a [2,3] to work on. Assuming everything else is good in it!
        # print("\n==============~Warning~==============");
        # print("Input size was [{},{}] and it is assumed that the dimensions are flipped -> adjusting to [{},{}] as [2,3] (m/d/yr) format is required (other than [2,2] daynum/Yr format)\n".format( len(dateRange[:,0]),len(dateRange[0,:]),len(dateRange[0,:]),len(dateRange[:,0]) ) );
        # dateRange = np.reshape(dateRange,(len(dateRange[0,:]),len(dateRange[:,0])));
        dateRange = np.swapaxes(dateRange,0,1); #flip axes
        
    elif (len(dateRange[0,:]) != 3) and (len(dateRange[:,0]) != 2 and (len(dateRange[0,:]) != 2) ):
        #Catch where a completely wrong size was sent as a [2,3] or [2,2] is required
        print("\n==============ERROR==============");
        print("Input size was [{},{}] and it needs to be [2,3] or [2,2] - exiting date to dayNum fun!\n".format( len(dateRange[:,0]),len(dateRange[0,:]) ) );
        # return "No"; #good luck with this return
        import sys
        sys.crash(); #even better
    #END IF
    
    if (len(dateRange[:,0]) == 2) and (len(dateRange[0,:]) == 3):
        #Catch if a Yr/M/D was sent but this will use Yr/dayNum format because it's easier to math with them
        dateRange = subfun_date_to_dayNum(dateRange); #call it dateRange to make programming easy, but now it is really dateRange_dayNum
    #END IF

    #==============Convert date range to the full range in between, covering years==============
    dateRange_dayNum_full = dateRange[1,0] - dateRange[0,0]; #yr difference
    if dateRange_dayNum_full == 0:
        dateRange_dayNum_full = dateRange[1,1] - dateRange[0,1]; #day difference
        dateRange_dayNum_full = np.hstack( [np.tile(dateRange[0,0],(dateRange_dayNum_full+1,1)) , np.reshape(np.arange(dateRange[0,1],dateRange[1,1]+1,1,dtype="int16"), (-1,1))] ); #create full date range from min/max days since within same year
        #way cooler in matlab, but works - gets the days in between as an array - note that reshapeing a (-1,1) means that -1 automatically calcs size
    #     dateRange_yearRange = dateRange(1,0); #set the year range as one year
    else:
        dateRange_yearRange = np.arange(dateRange[0,0],dateRange[1,0]+1,1,dtype="int16") #get the full year range from min to max
        for i in range(0, len(dateRange_yearRange) ):
            #Leap Year Detection
            if np.mod(dateRange_yearRange[i],4) == 0: #leap year
                #LEAP YEAR! Possibly.
                if (np.mod(dateRange_yearRange[i],100) == 0) and (np.mod(dateRange_yearRange[i],400) != 0):
                    #NO LEAP YEAR
                    #Leap Year Skipped Detected - next will be 2100
                    if np.isscalar(dateRange_dayNum_full) == 1: #see if date range is being used a temp var or not
                        dateRange_dayNum_full = np.hstack( [np.tile(dateRange_yearRange[i],(365,1)) , np.reshape(np.arange(1,365+1,1,dtype="int16"),(-1,1))] ); #create the date range if not it yet
                    else:
                        dateRange_dayNum_full = np.vstack( [dateRange_dayNum_full, np.hstack( [np.tile(dateRange_yearRange[i],(365,1)) , np.reshape(np.arange(1,365+1,1,dtype="int16"),(-1,1))] ) ] ); #if exists, tack on
                    #END IF
                else:
                    #Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
                    if np.isscalar(dateRange_dayNum_full) == 1: #see if date range is being used a temp var or not
                        dateRange_dayNum_full = np.hstack( [np.tile(dateRange_yearRange[i],(366,1)), np.reshape(np.arange(1,366+1,1,dtype="int16"),(-1,1))] ); #create the date range if not it yet
                    else:
                        dateRange_dayNum_full = np.vstack( [dateRange_dayNum_full, np.hstack( [np.tile(dateRange_yearRange[i],(366,1)) , np.reshape(np.arange(1,366+1,1,dtype="int16"),(-1,1))] ) ] ); #if exists, tack on
                    #END IF
                #END IF
            else: #no leap year if this
                #no leap year
                if np.isscalar(dateRange_dayNum_full) == 1: #see if date range is being used a temp var or not
                    dateRange_dayNum_full = np.hstack( [np.tile(dateRange_yearRange[i],(365,1)) , np.reshape(np.arange(1,365+1,1,dtype="int16"),(-1,1))] ); #create the date range if not it yet
                else:
                    
                    dateRange_dayNum_full = np.vstack( [dateRange_dayNum_full, np.hstack( [np.tile(dateRange_yearRange[i],(365,1)) , np.reshape(np.arange(1,365+1,1,dtype="int16"),(-1,1))] ) ] ); #if exists, tack on
                #END IF
            #END IF 
       #END FOR loop per year
       
        #dateRange_dayNum_full = np.vstack( (dateRange_dayNum_full,np.hstack( (np.reshape(np.arange(1,dateRange[1,1]+1,1,dtype="int16"), (-1,1)),np.tile(dateRange[1,0],(dateRange[1,1] - 1 +1,1))) )) );
        #I cut this later, didn't need to make this
        
        dateRange_dayNum_full_min = np.asscalar( np.where( (dateRange_dayNum_full[:,1] == dateRange[0,1]) & (dateRange_dayNum_full[:,0] == dateRange[0,0]) )[0] ); #get the min day to start on
        dateRange_dayNum_full_max = np.asscalar( np.where( (dateRange_dayNum_full[:,1] == dateRange[1,1]) & (dateRange_dayNum_full[:,0] == dateRange[1,0]) )[0] ); #get the max day to start on
        dateRange_dayNum_full = dateRange_dayNum_full[dateRange_dayNum_full_min:dateRange_dayNum_full_max+1, : ]; #cut out the extra
    #END IF
        
    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #convert all the day number dates to Yr/M/D dates
    
    return(dateRange_full,dateRange_dayNum_full); #return both
#END DEF

"""
Makes numbers look nice
RD, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
"""
def textNice(number):
    if( not isinstance(number,str) ):
        text = str(number); #convert to a string
    else:
        text = number; #already a string
    #END IF
    
    if( (text[-1] != '0') | (text.find('.') > -1) ):
        text = text.rstrip('0'); #only strip 0 if there's a decimal (prevents removing like 10 -> 1)
    #END IF
    text = text.rstrip('.'); #always want to strip a . on the end of a number b/c useless info if 10. for plotting
    
    return text #success
#END DEF

#GOAL: Find all occurances of a string B within a string A
#RD on 8/27/2018, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
#Greatly inspired by https://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
#Syntax looks like & acts like C
#
#INPUT: string to find (B) and string to search (A)
#OUTPUT: Indexes of found string (B) in string (A)

def strstr(strSearch, strFind):
    strSearch = "%r"%strSearch; #this is a workaround to allow for finding of \ in \n [this is very important don't disable it]
    #this requires strSearch[1:-1] to be used
    
    cntr = 0; #counter
    strFound = []; #init list that holds found indexes per string
    while cntr != -1:
        #use [1:-1] for "%r" mode above
        cntr = strSearch[1:-1].find(strFind, cntr); #search string for strFind (B) starting at cntr, which essentially moves the search up the strSearch (A) by 1 each loop
        # cntr = strSearch.find(strFind, cntr); #search string for strFind (B) starting at cntr, which essentially moves the search up the strSearch (A) by 1 each loop
        if cntr != -1: #watch for when strFind (B) is not found
            strFound.append(cntr); #record cntr
            cntr = cntr + 1; #Increment to next char
            #To make it go faster, use +len(strFind) instead of +1 - but it will not find overlapping strings
        #END IF
    #END WHILE
    
    strFound = np.asarray(strFound); #I hate lists 
    return strFound; #time to exit when -1 hits as that means no string was found
#END DEF

#GOAL: Search an array or list of strings for a single string
#RD on 8/27/2018, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
#Greatly inspired by https://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
#Syntax looks like & acts like C
#
#INPUT: array or list of strings to search (A), string to find in said array/list (B)
#OUTPUT: Boolean array of true/false, opt==1 causes the number of times it was found
#option opt, 0 == returns array size of strSearch with # occurances found, 1 == returns the total number of occurances found as 1 integer value

def strfind(strSearch, strFind, opt=0):
    if( type(strSearch).__module__ != np.__name__ ):
        strSearch = np.array(strSearch); #convert to array
    #END IF
    
    strFound = np.zeros(strSearch.size,dtype=np.int64); #preallocate array to find
    for i in range(0,strSearch.size):
        strFound[i] = strstr(strSearch[i],strFind).size; #record how many times it was found
    #END FOR
    
    if(opt == 1):
        strFound = np.sum(strFound); #sum em up
    #END IF
    
    return strFound #returns how many were found where (or just how many were found)
#END DEF

#RD, licensed GPL-3.0 @ https://github.com/dinsmoro/GRITI
def isin_row(a, b): #inspired by https://stackoverflow.com/a/67213552/2403531
    #calcs to see if each row of a is in b, yields numpy array that is the same length as a
    return np.asarray([row in set(map(tuple, b)) for row in map(tuple, a)])
#END DEF
