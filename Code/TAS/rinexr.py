#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# RINEX reader, only supports 3.04 .##0 files right now
import numpy as np
from scipy.integrate import solve_ivp
from re import search as research, match as rematch, compile as recompile, findall as refindall
from TAS.subfun_h5_writer import subfun_h5_writer_handler
from TAS.positizer import conv_ecef2geo, conv_aer2ecef
import TAS.dezee as dezee
from multiprocessing import Pool, cpu_count, set_start_method
from platform import system as platform_system
from copy import deepcopy
from datetime import datetime
from warnings import warn, simplefilter as warnings_simplefilter
warnings_simplefilter('always', UserWarning); # Global activation

#!!! OBSERVATION FUNCTION
def rinexr_obs(filePath, FLG_makeHDF5=False, FLG_makeHDF5_path=None, FLG_enableLeapSec=True, FLG_silent=False, FLG_parallel=True, parallel_CPUnum=None, parallel_threshold=1000, FLG_parallel_silent=False):
    # filePath = '/Users/RD36116/Documents/Projects/NovAtel/RINEX/BMGX20150005F_2025-04-21_11-45-41.25O'
    fileVersion = '1.0.0'; # Version 1
    # Check version with:
    # def verT(ver): # per https://stackoverflow.com/a/11887825, downsides ignored because external packages just for this capability are not worth it here
    #     return tuple(map(int, (ver.split("."))))
    # # END DEF
    
    # ----- Read file into memory -----
    if( filePath[filePath.rfind('.'):] == '.Z' ): # "Compact RINEX"
        with dezee.open(filePath) as filez:
            fileContent = filez.read().decode('utf-8').splitlines(); # Unzip and decode, then split
        # END WITH
    else:
        with open(filePath, 'r') as filez:
            fileContent = filez.read().splitlines();
            # with open(filePath, 'r') as filez:
            #     for linez in filez:
            #         # Process each line here
            #         print(linez) # Example: Print each line
            #     # Alternate for reading in batches
            #     # linez = [next(filez) for _ in range(batchRead)]
        # END WITH
    # END IF
    
    # ----- Identify header end location -----
    cntr = 0;
    while( cntr != -1 ):
        if( research(r'END OF HEADER\s*$', fileContent[cntr]) is not None ):
            headerEnd = cntr; # Rinex v2.10 at least also ends in the same way, so this should be "universal"
            cntr = -1; # Ditch
        else:
            cntr +=1 ; # Move on
        # END IF
    # END WHILE
    dataStart = headerEnd + 1; # Know where the data begins!~
    
    # ----- Identify RINEX version -----
    rinx = {'data':{}}; # A dict to hold the data
    # RINEX v2.10 at least also uses this header entry, so this should be "universal"
    for i in range(0, headerEnd):
        if( research(r'RINEX VERSION', fileContent[i]) is not None ):
             rinx['rinex version'] = fileContent[i][:20].strip(' '); # Go at it in a fixed-width way
             isObs = fileContent[i][20:40].strip(' '); # Go at it in a fixed-width way
             
             line2look = rematch(r'^\s*\d\.\d{1,2}', rinx['rinex version']); # Look for RINEX version
             if( line2look is None ):
                 raise Exception('ERROR in rinexr: RINEX VERSION line did NOT have a detectable RINEX version in file '+str(filePath)+'. Printing line:\n******\n'+fileContent[i]+'\n******');
             # END IF
             if( not ((isObs[0].lower() == 'o') or ('observation' in isObs.lower())) ):
                 raise Exception('ERROR in rinexr: the header does NOT declare that this is an observation file. Ditching! Here is the header line found:\n******\n'+fileContent[i]+'\n******');
             # END IF
             break; # Escape
        # END IF
    # END FOR i
    
    # ----- Extract all needed header info -----
    if( rinx['rinex version'] == '2.11' ):
        breakpoint()
        
    elif( rinx['rinex version'] == '3.04' ):
        # !! This does NOT care about FIXED WIDTH stuff and relies on SPACES for QUICK IMPLEMENTATION
        # !! REDO it if needed, see https://web.archive.org/web/20240721234005/https://gage.upc.edu/en/learning-materials/library/gnss-format-descriptions for help (amazing stuff from gAGE @ Technical University of Catalonia (UPC)
        # !! the recomps below will still be relevant for finding lines, but the per-line handling will need fixed width refurbishment
        recomp_pos = recompile(r'APPROX POSITION XYZ\s*$');
        recomp_obsTypes = recompile(r'SYS / # / OBS TYPES\s*$');
        recomp_sigStrUnit = recompile(r'SIGNAL STRENGTH UNIT\s*$');
        recomp_interval = recompile(r'INTERVAL\s*$');
        recomp_time_firstObs = recompile(r'TIME OF FIRST OBS\s*$');
        recomp_time_lastObs = recompile(r'TIME OF LAST OBS\s*$');
        recomp_phaseShift = recompile(r'SYS / PHASE SHIFT\s*$');
        recomp_Rslot2freq = recompile(r'GLONASS SLOT / FRQ #\s*$');
        recomp_RcodeBias = recompile(r'GLONASS COD/PHS/BIS\s*$');
        recomp_leapSec = recompile(r'LEAP SECONDS \s*$');
        # Roll through the whole thing and get everything relevant out
        for i in range(0, headerEnd):
            # Work through the header
            # --- XYZ Position ---
            line2look = research(recomp_pos, fileContent[i]);
            if( line2look is not None ):
                val2work = research(r'^\s*[0-9\-\.]+\s+[0-9\-\.]+\s+[0-9\-\.]+', fileContent[i]).group().strip(' ').split();
                rinx['pos'] = np.empty( 3, dtype=np.float64);
                for j in range(0, len(val2work)):
                    rinx['pos'][j] = float(val2work[j]); # Convert and record
                # END FOR j
            # END IF
            # --- Obs Types ---
            line2look = research(recomp_obsTypes, fileContent[i]);
            if( line2look is not None ):
                if( 'obs types' not in rinx ): rinx['obs types'] = {}; # Add in dynamically
                val2work = fileContent[i][:line2look.start()].strip(' ');
                val2work_newLineRex = rematch(r'^[a-zA-Z]\s+[0-9]+\s+', val2work); # Exists if this is a new line for a sat system
                if( val2work_newLineRex is not None ):
                    rinx['obs types'][val2work[0]] = val2work[val2work_newLineRex.end():].split(); # Put in as a list in a sub-dict of G/R/E/C/J/I/S
                else:
                    val2work_newLineRex = rematch(r'^[a-zA-Z]\s+[0-9]+\s+', fileContent[i-1]); # Go back to the previious line
                    rinx['obs types'][val2work_newLineRex.group().strip()[0]].extend(val2work.strip().split()); # Extend the extended line
                # END IF
            # END IF
            # --- Signal Strength Unit ---
            line2look = research(recomp_sigStrUnit, fileContent[i]);
            if( line2look is not None ):
                val2work = fileContent[i][:line2look.start()].strip(' ');
                rinx['sig str unit'] = val2work; # Just record as-is
            # END IF
            # --- Interval ---
            line2look = research(recomp_interval, fileContent[i]);
            if( line2look is not None ):
                val2work = fileContent[i][:line2look.start()].strip(' ');
                rinx['rate'] = float(val2work); # Record it just as-is
            # END IF
            # --- Time of First Obs ---
            line2look = research(recomp_time_firstObs, fileContent[i]);
            if( line2look is not None ):
                val2work = fileContent[i][:line2look.start()].strip(' ').split();
                if( val2work[5].find('.') == 1 ): val2work[5] = '0'+val2work[5]; # Tack on an extra 0 if it's 1 digit like 0.0000 -> 00.0000
                rinx['time first obs'] = np.array((val2work[0]+'-'+val2work[1]+'-'+val2work[2]+\
                    'T'+val2work[3]+':'+val2work[4]+':'+val2work[5]), dtype='datetime64[ns]'); # Create a numpy datetime64[ns] object
                rinx['time first obs source'] =  val2work[6]; # Record the source
            # END IF
            # --- Time of Last Obs ---
            line2look = research(recomp_time_lastObs, fileContent[i]);
            if( line2look is not None ):
                val2work = fileContent[i][:line2look.start()].strip(' ').split();
                if( val2work[5].find('.') == 1 ): val2work[5] = '0'+val2work[5]; # Tack on an extra 0 if it's 1 digit like 0.0000 -> 00.0000
                rinx['time last obs'] = np.array((val2work[0]+'-'+val2work[1]+'-'+val2work[2]+\
                    'T'+val2work[3]+':'+val2work[4]+':'+val2work[5]), dtype='datetime64[ns]'); # Create a numpy datetime64[ns] object
                rinx['time last obs source'] =  val2work[6]; # Record the source
            # END IF
            # --- Phase Shift (multi-line possible!) ---
            line2look = research(recomp_phaseShift, fileContent[i]);
            if( line2look is not None ):
                if( 'phase shift' not in rinx ): rinx['phase shift'] = {}; # Add in dynamically
                val2work = fileContent[i][:line2look.start()].strip(' ');
                if( rematch(r'[a-zA-Z]\s+', val2work) is not None ):
                    # This is an opener line
                    satType = val2work[0]; # Get the sat type
                    if( satType not in rinx['phase shift'].keys() ):
                        rinx['phase shift'][satType] = {}; # Make it
                    # END IF
                    code2work = val2work[1:rematch(r'^[a-zA-Z]\s+[a-zA-Z0-9]+', val2work).end()].strip(' '); # Get the code for the sat
                    if( code2work not in rinx['phase shift'][satType].keys() ):
                        rinx['phase shift'][satType][code2work] = {}; # Make it
                    # END IF
                    rinx['phase shift'][satType][code2work]['phase shift'] = float(rematch(r'^\s*[\-0-9\.]+', val2work[val2work.find(code2work)+len(code2work):]).group().strip(' '));
                    # rinx['phase shift'][satType][code2work]['#sats'] = rematch(r'^\s*[0-9]+', val2work[val2work.find(rinx['phase shift'][satType][code2work]['phase shift'])+len(rinx['phase shift'][satType][code2work]['phase shift']):]).group().strip(' ');
                    rinx['phase shift'][satType][code2work]['sats'] = val2work[rematch(r'^\s*[a-zA-Z]\s+[a-zA-Z0-9]+\s+[\-0-9\.]+\s+[0-9]+', val2work).end():].strip(' ').split();
                else:
                    # Otherwise this is a continuation of the line
                    # The satType and code2work vars are still relevant, use em
                    rinx['phase shift'][satType][code2work]['sats'].extend(val2work.strip(' ').split()); # Tack em on, easy!
                # END IF
            # END IF
            # --- GLONASS (R) slot<->frequency mapping (multi-line possible!) ---
            line2look = research(recomp_Rslot2freq, fileContent[i]);
            if( line2look is not None ):
                if( 'GLONASS slot/freq' not in rinx ): rinx['GLONASS slot/freq'] = {}; # Add in dynamically
                val2work = fileContent[i][:line2look.start()].strip(' ');
                if( rematch(r'[0-9]+\s+', val2work) is not None ):
                    # This is an opener line
                    satNum_loc = rematch(r'[0-9]+\s+', val2work); # Get it, will reuse
                    satNum = satNum_loc.group().strip(' '); # Get the sat type
                    rinx['GLONASS slot/freq']['#sats'] = int(satNum); # Record
                    if( 'mapping' not in rinx['GLONASS slot/freq'] ): rinx['GLONASS slot/freq']['mapping'] = {}; # Add in dynamically
                    satMap = val2work[satNum_loc.end():].strip(' ').split(); # Work it out
                    for j in range(0, len(satMap), 2):
                        rinx['GLONASS slot/freq']['mapping'][satMap[j]] = None; # Get the sat name set
                    # END FOR j
                    for j in range(1, len(satMap), 2):
                        rinx['GLONASS slot/freq']['mapping'][satMap[j-1]] = int(satMap[j]); # Get the freq offset (-7 to 6)
                    # END FOR j
                else:
                    # Otherwise this is a continuation of the line
                    satMap = val2work.strip(' ').split(); # Work it out
                    for j in range(0, len(satMap), 2):
                        rinx['GLONASS slot/freq']['mapping'][satMap[j]] = None; # Get the sat name set
                    # END FOR j
                    for j in range(1, len(satMap), 2):
                        rinx['GLONASS slot/freq']['mapping'][satMap[j-1]] = int(satMap[j]); # Get the freq offset (-7 to 6)
                    # END FOR j
                # END IF
            # END IF
            # --- GLONASS phase code bias ---
            line2look = research(recomp_RcodeBias, fileContent[i]);
            if( line2look is not None ):
                if( 'GLONASS code bias' not in rinx ): rinx['GLONASS code bias'] = {}; # Add in dynamically
                val2work = fileContent[i][:line2look.start()].strip(' ').split();
                for j in range(0, len(val2work), 2):
                    rinx['GLONASS code bias'][val2work[j]] = None; # Get the sat name set
                # END FOR j
                for j in range(1, len(val2work), 2):
                    rinx['GLONASS code bias'][val2work[j-1]] = float(val2work[j]); # Get the freq offset (-7 to 6)
                # END FOR j
            # END IF
            # --- Leap Seconds ---
            # line2look = research(recomp_leapSec, fileContent[i]);
            # if( line2look is not None ):
            #     val2work = fileContent[i][:line2look.start()].strip(' ').split();
            #     rinx['leap sec'] = int(val2work[0]); # Record it just as-is
            #     rinx['leap sec past&future'] = int(val2work[1]); # Record it just as-is
            #     rinx['GPS week'] = int(val2work[2]); # Record it just as-is
            #     if( 'GPS' in val2work[3] ):
            #         rinx['GPS dayNum'] = int(val2work[3][:val2work[3].find('GPS')]); # Record it just as-is
            #     else:
            #         raise Exception('ERROR in rinexr: unsupported day num in the leap second section of the RINEX header. Prob. BeiDou, never seen it so never coded it.');
            #     # END IF
            # # END IF
            line2look = research(recomp_leapSec, fileContent[i]);
            if( line2look is not None ):
                rinx['leap sec'] = int(fileContent[i][:6].strip(' ')); # Fixed width,
                rinx['leap sec past&future'] = int(fileContent[i][6:12].strip(' ')); # Fixed width,
                rinx['GNSS week'] = int(fileContent[i][12:18].strip(' ')); # Fixed width,
                rinx['GNSS dayNum'] = int(fileContent[i][18:24].strip(' ')); # Fixed width,
                rinx['GNSS identifier'] = fileContent[i][24:27].strip(' '); # Fixed width,
                if( rinx['GNSS identifier'] == '' ):
                    rinx['GNSS identifier'] = 'GPS'; # Blank defaults to GPS
                # END IF
            # END IF
        # END FOR i
        if( FLG_enableLeapSec ):
            # Convert GPS time to UTC (add in leap seconds)
            if( (rinx['time first obs source'] == 'GPS') and ('leap sec' in rinx) ):
                rinx['time first obs'] += np.timedelta64(rinx['leap sec'], 'ns'); # Convert to UTC from GPS by adding in leap seconds
            # END IF
            if( (rinx['time last obs source'] == 'GPS') and ('leap sec' in rinx) ):
                rinx['time last obs'] += np.timedelta64(rinx['leap sec'], 'ns'); # Convert to UTC from GPS by adding in leap seconds
            # END IF
        # END IF
    else:
        raise Exception('ERROR in rinexr: RINEX version of '+rinx['rinex version']+' found in file '+str(filePath)+' but that version is not currently supported.');
    # END IF
    
    # ----- Figure out new data lines -----
    rinexr_dataLines = []; # Prep a list
    if( rinx['rinex version'] != '3.04' ):
        breakpoint() # Add support for not 3.04 RINEX files at a future date
    elif( rinx['rinex version'] == '3.04' ):
        for i in range(dataStart, len(fileContent)):
            if( fileContent[i][0] == '>' ):
                rinexr_dataLines.append(i); # Record the data line
            # END IF
        # END FOR i
        rinexr_dataLines.append(len(fileContent)); # Append final line to help with loops later
    else:
        raise Exception('ERROR in rinexr: RINEX version of '+rinx['rinex version']+' found in file '+str(filePath)+' but that version is not currently supported.');
    # END IF
    
    # ----- Figure out if want parallel -----
    if( FLG_parallel > 0 ):
        if( parallel_CPUnum is None ):
            parallel_CPUnum = cpu_count(); # Use multiprocessing to get the cpu count, it includes tiny cores and does not figure out who is who yolo
        else:
            if( (FLG_silent == False) and (FLG_parallel_silent == False) and (cpu_count() < parallel_CPUnum) ):
                warn('WARNING in rinexr: Number of parallel processes requested `'+str(parallel_CPUnum)+'` exceeds the number of CPU cores found `'+str(cpu_count())+'`. This is probably less-than-ideal for performance.');
            # END IF
        # END IF
        if( (len(rinexr_dataLines)-1)/parallel_CPUnum < parallel_threshold ):
            FLG_parallel = False; # Disable, not enough data to make it worth it
            if( (FLG_silent == False) and (FLG_parallel_silent == False) ):
                warn('WARNING in rinexr: Parallel requested but with CPU num of '+str(parallel_CPUnum)+' and detected data line groups of '+str((len(rinexr_dataLines)-1))+' the estimated parallel ratio ('+str(np.round((len(rinexr_dataLines)-1)/parallel_CPUnum,2))+') is less than the parallel threshhold of '+str(parallel_threshold)+'.\nTurning off parallelization as the overhead is not worth it.');
            # END IF
        else:
            parallel_sets = np.arange(0, len(rinexr_dataLines)-1, int((len(rinexr_dataLines)-1)/parallel_CPUnum));
            parallel_sets[-1] = len(rinexr_dataLines)-1; # Go to the end
            # Get the OS
            osname = platform_system();
            if( osname == 'Darwin' ):
                set_start_method('fork', force=True); # Fork is faster than spawn and what we're doing won't run into fork issues on Mac
            # END IF
        # END IF
    # END IF
    
    # ----- Build the recording entries and record that data! -----
    if( rinx['rinex version'] == '3.04' ):
        # ----- Call the rinexReader_obs function as needed! -----
        if( FLG_parallel == 0 ):
            rinx = rinexReader_obs(rinx, fileContent, rinexr_dataLines); # Do it without parallelization
        else:
            # Prepare for parallel
            parallel_list = []; #Prep
            for j in range(len(parallel_sets)-1, 0, -1): # Do it backwards so we can delete
                 parallel_list.append([ rinexReader_obs, [deepcopy(rinx), fileContent[rinexr_dataLines[parallel_sets[j-1]]:rinexr_dataLines[parallel_sets[j]]], np.asarray(rinexr_dataLines[parallel_sets[j-1]:parallel_sets[j]+1])-rinexr_dataLines[parallel_sets[j-1]] ] ]); # pack up all needed function calls
                 # del fileContent[rinexr_dataLines[parallel_sets[j-1]]:rinexr_dataLines[parallel_sets[j]]]; # Save RAM
            #END FOR j
            
            #--- execute function on list of inputs ---
            with Pool(processes=cpu_count()) as executor:
                 results = executor.starmap(parallel_starmap_helper, parallel_list); #function you want is replaced with; parallel_starmap_helper helps starmap distribute everything right
            #END WITH
            
            #--- get the results out ---
            # for j in range(0, len(results)):
            for j in range(len(results)-1, -1, -1): # Do it backwards cause did it backwards up there
                getResults4Dict(results[j]['data'], rinx['data']); # In place
            # END FOR j
            # Collapse all of the data that's list-of-stuff
            collapse4dict(rinx['data']); # In place
        # END IF
    else:
        raise Exception('ERROR in rinexr: RINEX version of '+rinx['rinex version']+' found in file '+str(filePath)+' but that version is not currently supported.');
    # END IF
    if( FLG_enableLeapSec ):
        rinx['data']['time'] += np.timedelta64(rinx['leap sec'], 'ns'); # UTC, Adjust the entire time array by the leap seconds to bring it to UTC
    # END IF
    
    if( FLG_makeHDF5 == True ):
        if( FLG_makeHDF5_path is None ):
            FLG_makeHDF5_path = filePath+'.h5'; # Reuse the var
        # END IF
        if( FLG_silent == False ):
            from time import time
            tic = time(); # Start timing
            print('INFO in rinexr_obs: RINEX obs file read, saving as HDF5 file now. This can take a bit for big ones.');
        # END IF
        subfun_h5_writer_handler(FLG_makeHDF5_path, rinx, h5_path_prefix=None, keyz_ignore=None, topLevel_attributes={'version':fileVersion}, FLG_overWrite=False);
        if( FLG_silent == False ):
            toc = time();
            print('Saving as an HDF5 file took '+str(np.round(toc-tic, 2))+' sec / '+str(np.round((toc-tic)/60, 2))+' min / '+str(np.round((toc-tic)/3600, 2))+' hr.');
        # END IF
    # END IF
    
    return rinx
# END DEF

def rinexReader_obs(rinx, fileContent, rinexr_dataLines):
    rinx['data']['time'] = np.zeros( (len(rinexr_dataLines)-1) , dtype='datetime64[ns]');
    rinx['data']['epoch flag'] = np.zeros( (len(rinexr_dataLines)-1) , dtype=np.uint8); # small scale integer
    rinx['data']['receiver clock offset'] = np.nan*np.zeros( (len(rinexr_dataLines)-1) , dtype=np.float64); # receiver clock offset
    for gnssType in rinx['obs types'].keys():
        # Rev up the GNSS stuff
        rinx['data'][gnssType] = {}; # New dict
        # Figure out how many sat obs there are TOTAL
        rinx['data'][gnssType]['sats']  = set([]); # Prep a set which will help with only capturing unique sats
        for keyz in rinx['phase shift'][gnssType].keys():
            rinx['data'][gnssType]['sats'].update(rinx['phase shift'][gnssType][keyz]['sats']); # Update the set
        # END FOR keyz
        rinx['data'][gnssType]['sats'] = sorted(list(rinx['data'][gnssType]['sats'] )); # Convert back to list, sort it
        rinx['data'][gnssType]['sats2index'] = {rinx['data'][gnssType]['sats'][jk]: jk for jk in range(0, len(rinx['data'][gnssType]['sats']))}; # Make a sat name to index converter
        # Preallocate the needed arrays
        for i in range(0, len(rinx['obs types'][gnssType])):
            rinx['data'][gnssType][rinx['obs types'][gnssType][i]] = np.nan*np.zeros( (len(rinexr_dataLines)-1, len(rinx['data'][gnssType]['sats'])), dtype=np.float64); # it req float64
            if( rinx['obs types'][gnssType][i][0] == 'C' ):
                # Has an associated signal strength indicator (SSI)
                rinx['data'][gnssType][rinx['obs types'][gnssType][i]+'_SSI'] = np.zeros( (len(rinexr_dataLines)-1, len(rinx['data'][gnssType]['sats'])), dtype=np.uint8); # small scale integer
            elif( rinx['obs types'][gnssType][i][0] == 'L' ):
                # Has an associated signal strenght indicator (SSI) and loss of lock indicator (LLI)
                rinx['data'][gnssType][rinx['obs types'][gnssType][i]+'_SSI'] = np.zeros( (len(rinexr_dataLines)-1, len(rinx['data'][gnssType]['sats'])), dtype=np.uint8); # small scale integer
                rinx['data'][gnssType][rinx['obs types'][gnssType][i]+'_LLI'] = np.zeros( (len(rinexr_dataLines)-1, len(rinx['data'][gnssType]['sats'])), dtype=np.uint8); # small scale integer
            # END IF
        # END FOR i
    # END FOR gnssType

    # ----- Record that data! -----
    # Declare fixed length things
    headerFormat = (0, 5, 3, 3, 3, 3, 11, 3, 3, 6, 15); # Lengths, skips 1st '>' so starts at 0
    headerFormat_cumsum = np.cumsum(headerFormat); # Cumulative sum of the lengths
    dataFormat = 3; # Satellite name size, first part of line before the repeating data format size
    dataFormat_repeat = (14, 1, 1); # This repeats after the initial data size
    dataFormat_repeat_cumsum = np.cumsum(dataFormat_repeat); # Cumulative sum of the lengths
    dataFormat_repeat_sum = np.sum(dataFormat_repeat).item(); # Sum of the data sizes
    dataFormat_repeat_recomp = recompile(r'.'*dataFormat_repeat_sum);            
    
    for i in range(0, len(rinexr_dataLines)-1):
        # ----- HEADER -----
        headerLine = fileContent[rinexr_dataLines[i]][1:]; # Get the header line, lop off the start
        # --- First step, get the time ---
        sec2sec = headerLine[headerFormat_cumsum[5]:headerFormat_cumsum[6]].strip(); # Get the sec, need to standardize it for numpy
        if( sec2sec.find('.') == 1 ):
            sec2sec = '0'+sec2sec; # Tack on an extra 0 if it's 1 digit like 0.0000 -> 00.0000
        # END IF
        rinx['data']['time'][i] = \
            np.array((headerLine[headerFormat_cumsum[0]:headerFormat_cumsum[1]].strip()+'-'+headerLine[headerFormat_cumsum[1]:headerFormat_cumsum[2]].strip()+'-'+headerLine[headerFormat_cumsum[2]:headerFormat_cumsum[3]].strip()+\
            'T'+headerLine[headerFormat_cumsum[3]:headerFormat_cumsum[4]].strip()+':'+headerLine[headerFormat_cumsum[4]:headerFormat_cumsum[5]].strip()+':'+sec2sec), dtype='datetime64[ns]'); # GPS time, Create a numpy datetime64[ns] object
        # --- Epoch flag ---
        rinx['data']['epoch flag'][i] = int(headerLine[headerFormat_cumsum[6]:headerFormat_cumsum[7]].strip());
        # --- Epoch sat # ---
        dataLine_satNum = int(headerLine[headerFormat_cumsum[7]:headerFormat_cumsum[8]].strip()); # Don't record, just use for error checking
        # --- Receiver Clock Offset ---
        rinx['data']['receiver clock offset'][i] = np.float64(headerLine[headerFormat_cumsum[9]:headerFormat_cumsum[10]].strip()); 
        
        # ----- DATA -----
        dataLines = fileContent[rinexr_dataLines[i]+1:rinexr_dataLines[i+1]];
        if( len(dataLines) != dataLine_satNum ): raise Exception('ERROR in rinexr: Length of detected data lines does NOT match declared number of data lines in the data header line. Data header line and data lines follows:\n***\n'+fileContent[rinexr_dataLines[i]]+'\n'+'\n'.join(dataLines)+'\n***');
        for j in range(0, len(dataLines)):
            satName = dataLines[j][0:dataFormat];
            if( len(dataLines[j][dataFormat:]) == dataFormat_repeat_sum*len(rinx['obs types'][satName[0]]) ):
                dataTypes = refindall(dataFormat_repeat_recomp, dataLines[j][dataFormat:]); # Got em
                for kj in range(0, len(rinx['obs types'][satName[0]])):
                    if( dataTypes[kj][0:dataFormat_repeat[0]].strip(' ') != '' ):
                        rinx['data'][satName[0]][rinx['obs types'][satName[0]][kj]][i, rinx['data'][satName[0]]['sats2index'][satName]] = np.float64(dataTypes[kj][0:dataFormat_repeat[0]]); # Load in that data
                        if( (dataTypes[kj][dataFormat_repeat_cumsum[1]:dataFormat_repeat_cumsum[2]] != ' ') and (rinx['obs types'][satName[0]][kj][0] == 'C') ):
                            # Has an associated signal strength indicator (SSI)
                            rinx['data'][satName[0]][rinx['obs types'][satName[0]][kj]+'_SSI'][i, rinx['data'][satName[0]]['sats2index'][satName]] = np.uint8(dataTypes[kj][dataFormat_repeat_cumsum[1]:dataFormat_repeat_cumsum[2]]); # Get it
                        elif( rinx['obs types'][satName[0]][kj][0] == 'L' ):
                            # Has an associated signal strenght indicator (SSI) and loss of lock indicator (LLI)
                            if( dataTypes[kj][dataFormat_repeat_cumsum[1]:dataFormat_repeat_cumsum[2]] != ' ' ):
                                rinx['data'][satName[0]][rinx['obs types'][satName[0]][kj]+'_SSI'][i, rinx['data'][satName[0]]['sats2index'][satName]] = np.uint8(dataTypes[kj][dataFormat_repeat_cumsum[1]:dataFormat_repeat_cumsum[2]]); # Get it
                            # END IF
                            if( dataTypes[kj][dataFormat_repeat_cumsum[0]:dataFormat_repeat_cumsum[1]] != ' ' ):
                                rinx['data'][satName[0]][rinx['obs types'][satName[0]][kj]+'_LLI'][i, rinx['data'][satName[0]]['sats2index'][satName]] = np.uint8(dataTypes[kj][dataFormat_repeat_cumsum[0]:dataFormat_repeat_cumsum[1]]); # Get it
                            # END IF
                        # END IF
                    # END IF
                # END FOR kj
            else:
                raise Exception('ERROR in rinexr: Data line length does not match expected length. '+\
                                'Sat is '+satName+', obs type # is '+str(len(rinx['obs types'][satName[0]]))+\
                                ', fixed data length is '+str(dataFormat_repeat_sum)+\
                                '.\nTotal data length expected is '+str(dataFormat_repeat_sum*len(rinx['obs types'][satName[0]]))+\
                                ', total data length found is '+str(len(dataLines[j][dataFormat:]))+\
                                '. Data line on line '+str(rinexr_dataLines[i]+2+j)+' follows:\n***\n'+\
                                dataLines[j]+'\n***');
            # END IF
        # END FOR j
    # END FOR i
    
    return rinx
# END DEF


#!!! NAVIGATION FUNCTION
def rinexr_nav(filePath, FLG_makeHDF5=False, FLG_makeHDF5_path=None, FLG_enableLeapSec=True, FLG_silent=False, FLG_parallel=True, parallel_CPUnum=None, parallel_threshold=1000, FLG_parallel_silent=False):
    # filePath = '/Users/RD36116/Documents/Projects/NovAtel/RINEX/BMGX20150005F_2025-04-21_11-45-41.25O'
    fileVersion = '1.0.0'; # Version 1
    # Check version with:
    # def verT(ver): # per https://stackoverflow.com/a/11887825, downsides ignored because external packages just for this capability are not worth it here
    #     return tuple(map(int, (ver.split("."))))
    # # END DEF
    
    # ----- Read file into memory -----
    with open(filePath, 'r') as filez:
        fileContent = filez.read().splitlines();
        # with open(filePath, 'r') as filez:
        #     for linez in filez:
        #         # Process each line here
        #         print(linez) # Example: Print each line
        #     # Alternate for reading in batches
        #     # linez = [next(filez) for _ in range(batchRead)]
    # END WITH
    
    # ----- Identify header end location -----
    cntr = 0;
    while( cntr != -1 ):
        if( research(r'END OF HEADER\s*$', fileContent[cntr]) is not None ):
            headerEnd = cntr; # Rinex v2.10 at least also ends in the same way, so this should be "universal"
            cntr = -1; # Ditch
        else:
            cntr +=1 ; # Move on
        # END IF
    # END WHILE
    dataStart = headerEnd + 1; # Know where the data begins!~
    
    # ----- Identify RINEX version -----
    rinx = {'nav':{}}; # A dict to hold the data
    # RINEX v2.10 at least also uses this header entry, so this should be "universal"
    for i in range(0, headerEnd):
        if( research(r'RINEX VERSION', fileContent[i]) is not None ):
             rinx['rinex version'] = fileContent[i][:20].strip(' '); # Go at it in a fixed-width way
             isObs = fileContent[i][20:40].strip(' ').lower(); # Go at it in a fixed-width way
             navType = fileContent[i][40:60].strip(' ').lower(); # Go at it in a fixed-width way
             
             line2look = rematch(r'^\s*\d\.\d{1,2}', rinx['rinex version']); # Look for RINEX version
             if( line2look is None ):
                 raise Exception('ERROR in rinexr: RINEX VERSION line did NOT have a detectable RINEX version in file '+str(filePath)+'. Printing line:\n******\n'+fileContent[i]+'\n******');
             # END IF
             if( not ((isObs[0] == 'n') or ('navigation' in isObs)) ):
                 raise Exception('ERROR in rinexr: the header does NOT declare that this is an navigation file. Ditching! Here is the header line found:\n******\n'+fileContent[i]+'\n******');
             # END IF
             break; # Escape
        # END IF
    # END FOR i
    
    # ----- Extract all needed header info -----
    if( rinx['rinex version'] == '3.04' ):
        # For ez file reference see https://web.archive.org/web/20240721234005/https://gage.upc.edu/en/learning-materials/library/gnss-format-descriptions for help (amazing stuff from gAGE @ Technical University of Catalonia (UPC)
        recomp_iono = recompile(r'IONOSPHERIC CORR\s*$');
        recomp_time = recompile(r'TIME SYSTEM CORR\s*$');
        recomp_leapSec = recompile(r'LEAP SECONDS \s*$');
        
        # Roll through the whole thing and get everything relevant out
        for i in range(0, headerEnd):
            # Work through the header
            # --- Ionospheric Correction ---
            line2look = research(recomp_iono, fileContent[i]);
            if( line2look is not None ):
                val2work = fileContent[i][:5].strip(' '); # Fixed width, name of the correction
                rinx[val2work] = {}; # Dictionary time!
                if( val2work.lower() == 'gpsa' ):
                    for j in range(0, 4):
                        rinx[val2work]['alpha'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width,
                    # END FOR j
                elif( val2work.lower() == 'gpsb' ):
                    for j in range(0, 4):
                        rinx[val2work]['beta'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width, 
                    # END FOR j
                elif( val2work.lower() == 'gal' ):
                    for j in range(0, 3):
                        rinx[val2work]['ai'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width,
                    # END FOR j
                elif( val2work.lower() == 'qzsa' ):
                    for j in range(0, 4):
                        rinx[val2work]['alpha'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width,
                    # END FOR j
                elif( val2work.lower() == 'qzsb' ):
                    for j in range(0, 4):
                        rinx[val2work]['beta'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width, 
                    # END FOR j
                elif( val2work.lower() == 'bdsa' ):
                    val2work_prn = int(fileContent[i][55:58].strip(' ')); # Get the satellite number (PRN)
                    if val2work_prn not in rinx[val2work]: rinx[val2work][val2work_prn] = {}; # Another sub-dict if not already present
                    for j in range(0, 4):
                        rinx[val2work][val2work_prn]['alpha'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width,
                    # END FOR j
                    rinx[val2work][val2work_prn]['time mark code alpha'] = fileContent[i][53:55].strip(' '); # Time code, A=00h-01h, etc.
                    rinx[val2work][val2work_prn]['time mark alpha'] = str(ord(rinx[val2work][val2work_prn]['time mark code alpha'].lower()) - 96-1).zfill(2) + ' to ' + str(ord(rinx[val2work][val2work_prn]['time mark code alpha'].lower()) - 96).zfill(2); # Convert to usable time range in hrs, in BDS time
                elif( val2work.lower() == 'bdsb' ):
                    val2work_prn = int(fileContent[i][55:58].strip(' ')); # Get the satellite number (PRN)
                    if val2work_prn not in rinx[val2work]: rinx[val2work][val2work_prn] = {}; # Another sub-dict!
                    for j in range(0, 4):
                        rinx[val2work][val2work_prn]['beta'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width,
                    # END FOR j
                    rinx[val2work][val2work_prn]['time mark code beta'] = fileContent[i][53:55].strip(' '); # Time code, A=00h-01h, etc.
                    rinx[val2work][val2work_prn]['time mark beta'] = str(ord(rinx[val2work][val2work_prn]['time mark code beta'].lower()) - 96-1).zfill(2) + ' to ' + str(ord(rinx[val2work][val2work_prn]['time mark code beta'].lower()) - 96).zfill(2); # Convert to usable time range in hrs, in BDS time
                elif( val2work.lower() == 'irna' ):
                    for j in range(0, 4):
                        rinx[val2work]['alpha'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width,
                    # END FOR j
                elif( val2work.lower() == 'irnb' ):
                    for j in range(0, 4):
                        rinx[val2work]['beta'+str(j)] = float(fileContent[i][5+12*j:5+12*(j+1)].strip(' ').replace('D','E')); # Fixed width, 
                    # END FOR j
                else:
                    raise Exception('ERROR in rinexr: In file '+str(filePath)+' unsupported ionospheric correction type of `'+val2work+'` found. Printing line:\n*****\n'+fileContent[i]+'\n*****');
                # END IF
            # END IF
            # --- Time System Correction ---
            line2look = research(recomp_time, fileContent[i]);
            if( line2look is not None ):
                val2work = fileContent[i][:5].strip(' '); # Fixed width, name of the correction
                rinx[val2work] = {}; # Dictionary time!
                rinx[val2work]['a0'] = float(fileContent[i][5:22].strip(' ').replace('D','E')); # Fixed width,
                rinx[val2work]['a1'] = float(fileContent[i][22:38].strip(' ').replace('D','E')); # Fixed width,
                rinx[val2work]['ref time'] = fileContent[i][38:45].strip(' '); # Fixed width,
                if( rinx[val2work]['ref time'] == '' ):
                    rinx[val2work]['ref time'] = None; # Set to none
                else:
                    rinx[val2work]['ref time'] = int(rinx[val2work]['ref time']); # Integer it if it exists
                # END IF
                rinx[val2work]['continuous week num'] = fileContent[i][45:50].strip(' '); # Fixed width,
                if( rinx[val2work]['continuous week num'] == '' ):
                    rinx[val2work]['continuous week num'] = None; # Set to none
                else:
                    rinx[val2work]['continuous week num'] = int(rinx[val2work]['continuous week num']); # Integer it if it exists
                # END IF
                rinx[val2work]['augmentation system'] = fileContent[i][50:57].strip(' '); # Fixed width,
                if( rinx[val2work]['augmentation system'] == '' ):
                    rinx[val2work]['augmentation system'] = None; # Set to none
                # END IF
                rinx[val2work]['UTC identifier'] = fileContent[i][57:60].strip(' '); # Fixed width,
                if( rinx[val2work]['UTC identifier'] == '' ):
                    rinx[val2work]['UTC identifier'] = None; # Set to none
                else:
                    rinx[val2work]['UTC identifier'] = int(rinx[val2work]['UTC identifier']); # Integer it if it exists
                # END IF
            # END IF
            # --- Leap Seconds ---
            line2look = research(recomp_leapSec, fileContent[i]);
            if( line2look is not None ):
                rinx['leap sec'] = int(fileContent[i][:6].strip(' ')); # Fixed width,
                rinx['leap sec past&future'] = int(fileContent[i][6:12].strip(' ')); # Fixed width,
                rinx['GNSS week'] = int(fileContent[i][12:18].strip(' ')); # Fixed width,
                rinx['GNSS dayNum'] = int(fileContent[i][18:24].strip(' ')); # Fixed width,
                rinx['GNSS identifier'] = fileContent[i][24:27].strip(' '); # Fixed width,
                if( rinx['GNSS identifier'] == '' ):
                    rinx['GNSS identifier'] = 'GPS'; # Blank defaults to GPS
                # END IF
            # END IF
        # END FOR i
    else:
        raise Exception('ERROR in rinexr: RINEX version of '+rinx['rinex version']+' found in file '+str(filePath)+' but that version is not currently supported.');
    # END IF
    
    # ----- Figure out new data lines -----
    rinexr_dataLines = []; # Prep a list
    # recomp_newline = recompile(r'^\s*['+navType+navType.upper()+r']\d{1,3}'); # Accept the heavier regex for more robustness since nav files are way smaller (for just the navType)
    recomp_newline = recompile(r'^\s*[a-zA-Z]\d{1,3}'); # Accept the heavier regex for more robustness since nav files are way smaller (for all types, this hedges against a "mixed" one if it could exist)
    if( rinx['rinex version'] == '3.04' ):
        for i in range(dataStart, len(fileContent)):
            if( rematch(recomp_newline, fileContent[i]) is not None ):
                rinexr_dataLines.append(i); # Record the data line
            # END IF
        # END FOR i
        rinexr_dataLines.append(len(fileContent)); # Append final line to help with loops later
    else:
        raise Exception('ERROR in rinexr: RINEX version of '+rinx['rinex version']+' found in file '+str(filePath)+' but that version is not currently supported.');
    # END IF
    
    # ----- Figure out if want parallel -----
    if( FLG_parallel > 0 ):
        if( parallel_CPUnum is None ):
            parallel_CPUnum = cpu_count(); # Use multiprocessing to get the cpu count, it includes tiny cores and does not figure out who is who yolo
        else:
            if( (FLG_silent == False) and (FLG_parallel_silent == False) and (cpu_count() < parallel_CPUnum) ):
                warn('WARNING in rinexr: Number of parallel processes requested `'+str(parallel_CPUnum)+'` exceeds the number of CPU cores found `'+str(cpu_count())+'`. This is probably less-than-ideal for performance.');
            # END IF
        # END IF
        if( (len(rinexr_dataLines)-1)/parallel_CPUnum < parallel_threshold ):
            FLG_parallel = False; # Disable, not enough data to make it worth it
            if( (FLG_silent == False) and (FLG_parallel_silent == False) ):
                warn('WARNING in rinexr: Parallel requested but with CPU num of '+str(parallel_CPUnum)+' and detected data line groups of '+str((len(rinexr_dataLines)-1))+' the estimated parallel ratio ('+str(np.round((len(rinexr_dataLines)-1)/parallel_CPUnum,2))+') is less than the parallel threshhold of '+str(parallel_threshold)+'.\nTurning off parallelization as the overhead is not worth it.');
            # END IF
        else:
            parallel_sets = np.arange(0, len(rinexr_dataLines)-1, int((len(rinexr_dataLines)-1)/parallel_CPUnum));
            parallel_sets[-1] = len(rinexr_dataLines)-1; # Go to the end
            # Get the OS
            osname = platform_system();
            if( osname == 'Darwin' ):
                set_start_method('fork', force=True); # Fork is faster than spawn and what we're doing won't run into fork issues on Mac
            # END IF
        # END IF
    # END IF
    
    # ----- Build the recording entries and record that data! -----
    if( rinx['rinex version'] == '3.04' ):
        # ----- Call the rinexReader_obs function as needed! -----
        if( FLG_parallel == 0 ):
            rinx = rinexReader_nav(rinx, fileContent, rinexr_dataLines, navType); # Do it without parallelization
        else:
            # Prepare for parallel
            parallel_list = []; #Prep
            for j in range(len(parallel_sets)-1, 0, -1): # Do it backwards so we can delete
                 parallel_list.append([ rinexReader_nav, [deepcopy(rinx), fileContent[rinexr_dataLines[parallel_sets[j-1]]:rinexr_dataLines[parallel_sets[j]]], np.asarray(rinexr_dataLines[parallel_sets[j-1]:parallel_sets[j]+1])-rinexr_dataLines[parallel_sets[j-1]], navType ] ]); # pack up all needed function calls
                 # del fileContent[rinexr_dataLines[parallel_sets[j-1]]:rinexr_dataLines[parallel_sets[j]]]; # Save RAM
            #END FOR j
        
            #--- execute function on list of inputs ---
            with Pool(processes=cpu_count()) as executor:
                 results = executor.starmap(parallel_starmap_helper, parallel_list); #function you want is replaced with; parallel_starmap_helper helps starmap distribute everything right
            #END WITH
            
            #--- get the results out ---
            # for j in range(0, len(results)):
            for j in range(len(results)-1, -1, -1): # Do it backwards cause did it backwards up there
                getResults4Dict(results[j]['nav'], rinx['nav']); # In place
            # END FOR j
            # Collapse all of the data that's list-of-stuff
            collapse4dict(rinx['nav']); # In place
        # END IF
    else:
        raise Exception('ERROR in rinexr: RINEX version of '+rinx['rinex version']+' found in file '+str(filePath)+' but that version is not currently supported.');
    # END IF
    if( FLG_enableLeapSec and ('NO DATA' not in rinx) ):
        rinx['nav']['time'] += np.timedelta64(rinx['leap sec'], 'ns'); # UTC, Adjust the entire time array by the leap seconds to bring it to UTC
    # END IF
    
    if( (FLG_makeHDF5 == True) and ('NO DATA' not in rinx) ):
        if( FLG_makeHDF5_path is None ):
            FLG_makeHDF5_path = filePath+'.h5'; # Reuse the var
        # END IF
        if( FLG_silent == False ):
            from time import time
            tic = time(); # Start timing
            print('INFO in rinexr_nav: RINEX nav file read, saving as HDF5 file now.');
        # END IF
        subfun_h5_writer_handler(FLG_makeHDF5_path, rinx, h5_path_prefix=None, keyz_ignore=None, topLevel_attributes={'version':fileVersion}, FLG_overWrite=False);
        if( FLG_silent == False ):
            toc = time();
            print('Saving as an HDF5 file took '+str(np.round(toc-tic, 2))+' sec / '+str(np.round((toc-tic)/60, 2))+' min / '+str(np.round((toc-tic)/3600, 2))+' hr.');
        # END IF
    # END IF
    
    return rinx
# END DEF

def rinexReader_nav(rinx, fileContent, rinexr_dataLines, navType):
    # ----- Ensure there is data to be read -----
    if( (len(rinexr_dataLines)-1) != 0 ):
        # ----- Declare fixed length things -----
        headerFormat = (3, 5, 3, 3, 3, 3, 3, 19, 19, 19);# Header line lengths
        headerFormat_cumsum = np.cumsum(headerFormat); # Cumulative sum of the lengths
        dataFormat = (4, 19, 19, 19, 19);# Data line lengths
        dataFormat_cumsum = np.cumsum(dataFormat); # Cumulative sum of the lengths
            
        # ----- Figure out which sats are here, and how many times -----
        satz_cntr = {};
        for i in range(0, len(rinexr_dataLines)-1):
            # --- First step, sat ---
            satz = fileContent[rinexr_dataLines[i]][0:headerFormat[0]]; # Get the sat name
            rinx['nav'][satz] = {}; # Make a dict
            if( satz not in satz_cntr):
                satz_cntr[satz] = 1; # Set to 1
            else:
                satz_cntr[satz] += 1; # Increment
            # END IF
        # END FOR i
        satz_filler = {satz: 0 for satz in satz_cntr.keys()}; # Prep
        
        # ----- Preallocate known quantities -----
        for satz in satz_cntr.keys():
            rinx['nav'][satz]['time'] = np.zeros( (satz_cntr[satz]), dtype='datetime64[ns]'); # Every one has the time at least
            if( satz[0].lower() == 'g' ): # GPS
                quantz = ['clock bias', # sec, SV clock bias
                          'clock drift', # sec/sec, SV clock drift
                          'clock drift rate', # sec/sec^2, SV clock drift rate
                          'IODE', # Issue of data ephemeris
                          'CRS', # m, radius correction sine componenet
                          'deltaN', # rad/s
                          'Mo', # rad
                          'CUC', # rad, argument of latitude correction cosine component
                          'ecc', # eccentricity
                          'CUS', # rad, argument of latitude correction sine component
                          'sqrtA', # m^(1/2), sqrt of orbit semi major axis
                          'TOE', # sec (of GPS week), time of ephemeris
                          'CIC', # rad, inclination correction cosine component
                          'RAAN', # rad, Right Ascenion of the Ascending Node
                          'CIS', # rad, inclination correction sine component
                          'Io', # rad, initial inclination
                          'CRC', # m, radius correction of cosine component
                          'AOP', # rad, Argument of Pericenter
                          'RAAN rate', # rad/s, Nodal precession rate
                          'Io rate', # rad/s, inclination rate
                          'L2 codes', # Codes on L2 channel
                          'GPS week', # GPS week number, continuous not mod(1024)
                          'L2P data flag', # L2P data flag
                          'SV accuracy', # m, SV accuracy
                          'SV health', # 0 = SV health OK, not 0 not OK
                          'TGD', # sec, Total group delay
                          'IODC', # issue of data clock
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message
                          'fit interval', # hrs, curve fit interval used by GPS control segment in determining ephemeris parameters (0 = 4 hr, 1 = greater than 4 hr)
                          ]; # Names for easy looping
            elif( satz[0].lower() == 'r' ): # GLONASS
                quantz = ['clock bias', # sec, SV clock bias
                          'rel freq bias', # sec, SV relative frequency bias
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message
                          'Xpos', # km, position of satellite X component, likely in PZ-90.11 since 2014
                          'Xvel', # km/s, velocity of satellite X component, likely in PZ-90.11 since 2014
                          'Xacc', # km/s^2, acceleration of satellite X component, likely in PZ-90.11 since 2014
                          'SV health', # 0 = SV health OK, not 0 not OK
                          'Ypos', # km, position of satellite Y component, likely in PZ-90.11 since 2014
                          'Yvel', # km/s, velocity of satellite Y component, likely in PZ-90.11 since 2014
                          'Yacc', # km/s^2, acceleration of satellite Y component, likely in PZ-90.11 since 2014
                          'FCN', # Satellite's frequency channel number (FCN) (-7 to +6) since 2005
                          'Zpos', # km, position of satellite Z component, likely in PZ-90.11 since 2014
                          'Zvel', # km/s, velocity of satellite Z component, likely in PZ-90.11 since 2014
                          'Zacc', # km/s^2, acceleration of satellite Z component, likely in PZ-90.11 since 2014
                          'infoAge', # days, "Age of operative information"
                          ]; # Names for easy looping
            elif( satz[0].lower() == 'e' ): # Galileo
                quantz = ['clock bias', # sec, SV clock bias
                          'clock drift', # sec/sec, SV clock drift
                          'clock drift rate', # sec/sec^2, SV clock drift rate
                          'IODE', # Issue of data ephemeris
                          'CRS', # m, radius correction sine componenet
                          'deltaN', # rad/s
                          'Mo', # rad
                          'CUC', # rad, argument of latitude correction cosine component
                          'ecc', # eccentricity
                          'CUS', # rad, argument of latitude correction sine component
                          'sqrtA', # m^(1/2), sqrt of orbit semi major axis
                          'TOE', # sec (of GAL week), time of ephemeris
                          'CIC', # rad, inclination correction cosine component
                          'RAAN', # rad, Right Ascenion of the Ascending Node (OMEGA)
                          'CIS', # rad, inclination correction sine component
                          'Io', # rad, initial inclination
                          'CRC', # m, radius correction of cosine component
                          'AOP', # rad, Argument of Pericenter (Omega)
                          'RAAN rate', # rad/s, Nodal precession rate (OMEGA DOT)
                          'Io rate', # rad/s, inclination rate (IDOT)
                          'data sources', # 10 bit code that means things about the data sources (signals like E1, E5a)
                          'GAL week', # Galileo week number, continuous not mod(1024)
                          'BLANK', # Blank
                          'SISA', # m, Signal in space accuracy (SISA)
                          'SV health', # 9 bits for the different signals
                          'BGD1', # sec, broadcast group delay (E5a/E1)
                          'BGD2', # sec, broadcast group delay (E5b/E1)
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message
                          ]; # Names for easy looping
            elif( satz[0].lower() == 's' ): # SBAS
                quantz = ['clock bias', # sec, SV clock bias
                          'rel freq bias', # sec, SV relative frequency bias
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message
                          'Xpos', # km, position of satellite X component, coordinate system depends on system supported by SBAS
                          'Xvel/s', # km, velocity of satellite X component, coordinate system depends on system supported by SBAS
                          'Xacc/s^2', # km, acceleration of satellite X component, coordinate system depends on system supported by SBAS
                          'SV health', # 0 = SV health OK, not 0 not OK
                          'Ypos', # km, position of satellite Y component, coordinate system depends on system supported by SBAS
                          'Yvel/s', # km, velocity of satellite Y component, coordinate system depends on system supported by SBAS
                          'Yacc/s^2', # km, acceleration of satellite Y component, coordinate system depends on system supported by SBAS
                          'URA', # User range accuracy index, seems to be an index that converts to meters
                          'Zpos', # km, position of satellite Z component, coordinate system depends on system supported by SBAS
                          'Zvel/s', # km, velocity of satellite Z component, coordinate system depends on system supported by SBAS
                          'Zacc/s^2', # km, acceleration of satellite Z component, coordinate system depends on system supported by SBAS
                          'issue of data navigation', # days, D0229 8 first bits after Message Type if MT9 (???)
                          ]; # Names for easy looping
            elif( satz[0].lower() == 'j' ): # QZSS
                quantz = ['clock bias', # sec, SV clock bias
                          'clock drift', # sec/sec, SV clock drift
                          'clock drift rate', # sec/sec^2, SV clock drift rate
                          'IODE', # Issue of data ephemeris
                          'CRS', # m, radius correction sine componenet
                          'deltaN', # rad/s
                          'Mo', # rad
                          'CUC', # rad, argument of latitude correction cosine component
                          'ecc', # eccentricity
                          'CUS', # rad, argument of latitude correction sine component
                          'sqrtA', # m^(1/2), sqrt of orbit semi major axis
                          'TOE', # sec (of GPS week), time of ephemeris
                          'CIC', # rad, inclination correction cosine component
                          'RAAN', # rad, Right Ascenion of the Ascending Node
                          'CIS', # rad, inclination correction sine component
                          'Io', # rad, initial inclination
                          'CRC', # m, radius correction of cosine component
                          'AOP', # rad, Argument of Pericenter
                          'RAAN rate', # rad/s, Nodal precession rate
                          'Io rate', # rad/s, inclination rate
                          'L2 codes', # Codes on L2 channel
                          'GPS week', # GPS week number, continuous not mod(1024)
                          'L2P data flag', # L2P data flag; set to 1 since QZSS doesn't track it
                          'SV accuracy', # m, SV accuracy
                          'SV health', # 0 = SV health OK, not 0 not OK
                          'TGD', # sec, Total group delay
                          'IODC', # issue of data clock
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message
                          'fit interval', # hrs, curve fit interval used by GPS control segment in determining ephemeris parameters (0 = 4 hr, 1 = greater than 4 hr)
                          ]; # Names for easy looping
            elif( satz[0].lower() == 'c' ): # BeiDou
                quantz = ['clock bias', # sec, SV clock bias
                          'clock drift', # sec/sec, SV clock drift
                          'clock drift rate', # sec/sec^2, SV clock drift rate
                          'AODE', # Age of data ephemeris
                          'CRS', # m, radius correction sine componenet
                          'deltaN', # rad/s
                          'Mo', # rad
                          'CUC', # rad, argument of latitude correction cosine component
                          'ecc', # eccentricity
                          'CUS', # rad, argument of latitude correction sine component
                          'sqrtA', # m^(1/2), sqrt of orbit semi major axis
                          'TOE', # sec (of BDS week), time of ephemeris
                          'CIC', # rad, inclination correction cosine component
                          'RAAN', # rad, Right Ascenion of the Ascending Node (OMEGA)
                          'CIS', # rad, inclination correction sine component
                          'Io', # rad, initial inclination
                          'CRC', # m, radius correction of cosine component
                          'AOP', # rad, Argument of Pericenter (Omega)
                          'RAAN rate', # rad/s, Nodal precession rate (OMEGA DOT)
                          'Io rate', # rad/s, inclination rate (IDOT)
                          'BLANK', # Blank
                          'BDS week', # BeiDou week number, to go with TOE
                          'BLANK', # Blank
                          'SV accuracy', # m, SV accuracy; see documentation for formula
                          'SV health', # 9 bits for the different signals
                          'TGD1', # sec, total group delay (B1/B3)
                          'TGD2', # sec, broadcast group delay (B2/B3)
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message
                          'AODC', # sec, age of data clock
                          ]; # Names for easy looping
            elif( satz[0].lower() == 'i' ): # NavIC
                quantz = ['clock bias', # sec, SV clock bias
                          'clock drift', # sec/sec, SV clock drift
                          'clock drift rate', # sec/sec^2, SV clock drift rate
                          'IODEC', # Issue of data ephemeris and clock
                          'CRS', # m, radius correction sine componenet
                          'deltaN', # rad/s
                          'Mo', # rad
                          'CUC', # rad, argument of latitude correction cosine component
                          'ecc', # eccentricity
                          'CUS', # rad, argument of latitude correction sine component
                          'sqrtA', # m^(1/2), sqrt of orbit semi major axis
                          'TOE', # sec (of IRNSS week), time of ephemeris
                          'CIC', # rad, inclination correction cosine component
                          'RAAN', # rad, Right Ascenion of the Ascending Node
                          'CIS', # rad, inclination correction sine component
                          'Io', # rad, initial inclination
                          'CRC', # m, radius correction of cosine component
                          'AOP', # rad, Argument of Pericenter
                          'RAAN rate', # rad/s, Nodal precession rate
                          'Io rate', # rad/s, inclination rate
                          'BLANK', # Blank
                          'IRNSS week', # IRNSS week number (to go with TOE), continuous not mod(1024); same as GPS
                          'BLANK', # Blank
                          'user accuracy', # m, user range accuracy; see documentation for formula
                          'SV health', # 0 = L5 and S OK; 1 = L5 and S not OK; 2 = L5 not OK, S OK; 3 = L5 OK, S not OK
                          'TGD', # sec, Total group delay
                          'BLANK', # Blank
                          'trans time', # sec of GPS/GAL/BDS week, transmission time of message; seconds of IRNSS week are derived from Z-count in Hand-Over-World (HOW) (???)
                          ]; # Names for easy looping
            else:
                raise Exception('ERROR in rinexr: Unsupported sat type of '+satz+' found.');
            # END IF
            for quantities in quantz:
                rinx['nav'][satz][quantities] = np.zeros( (satz_cntr[satz]), dtype=np.float64);
            # END FOR quantities
        # END FOR satz
        
        # ----- Record that data! -----
        for i in range(0, len(rinexr_dataLines)-1):
            # ----- HEADER -----
            headerLine = fileContent[rinexr_dataLines[i]]; # Get the header line, lop off the start
            # --- First step, sat ---
            satz = headerLine[0:headerFormat[0]]; # Get the sat name
            
            # --- Second step, get the time ---
            sec2sec = headerLine[headerFormat_cumsum[5]:headerFormat_cumsum[6]].strip(); # Get the sec, need to standardize it for numpy
            if( sec2sec.find('.') == 1 ):
                sec2sec = '0'+sec2sec; # Tack on an extra 0 if it's 1 digit like 0.0000 -> 00.0000
            # END IF
            rinx['nav'][satz]['time'][satz_filler[satz]] = \
                np.array((headerLine[headerFormat_cumsum[0]:headerFormat_cumsum[1]].strip()+'-'+headerLine[headerFormat_cumsum[1]:headerFormat_cumsum[2]].strip()+'-'+headerLine[headerFormat_cumsum[2]:headerFormat_cumsum[3]].strip()+\
                'T'+headerLine[headerFormat_cumsum[3]:headerFormat_cumsum[4]].strip()+':'+headerLine[headerFormat_cumsum[4]:headerFormat_cumsum[5]].strip()+':'+sec2sec), dtype='datetime64[ns]'); # GPS time, Create a numpy datetime64[ns] object
            
            # --- Third step, get the data included in the header line and the rest of the data ---
            dataLines = fileContent[rinexr_dataLines[i]+1:rinexr_dataLines[i+1]];
            if( satz[0].lower() in ['g', 'e', 'j', 'c', 'i'] ): # GPS, Galileo, BeiDou (all share the same header; everything else is automated)
                # --- HEADER DATA ---
                rinx['nav'][satz]['clock bias'][satz_filler[satz]] = np.float64(headerLine[headerFormat_cumsum[6]:headerFormat_cumsum[7]].strip().replace('D','E')); # sec, SV clock bias
                rinx['nav'][satz]['clock drift'][satz_filler[satz]] = np.float64(headerLine[headerFormat_cumsum[7]:headerFormat_cumsum[8]].strip().replace('D','E')); # sec/sec, SV clock drift
                rinx['nav'][satz]['clock drift rate'][satz_filler[satz]] = np.float64(headerLine[headerFormat_cumsum[8]:headerFormat_cumsum[9]].strip().replace('D','E')); # sec/sec^2, SV clock drift rate
               
                # --- Per-system Offset ---
                quantz_cntr = 3; # Start after the header
                
            elif( satz[0].lower() in ['r', 's'] ): # GLONASS
                # --- HEADER DATA ---
                rinx['nav'][satz]['clock bias'][satz_filler[satz]] = np.float64(headerLine[headerFormat_cumsum[6]:headerFormat_cumsum[7]].strip().replace('D','E')); # sec, SV clock bias
                rinx['nav'][satz]['rel freq bias'][satz_filler[satz]] = np.float64(headerLine[headerFormat_cumsum[7]:headerFormat_cumsum[8]].strip().replace('D','E')); # sec/sec, SV clock drift
                rinx['nav'][satz]['trans time'][satz_filler[satz]] = np.float64(headerLine[headerFormat_cumsum[8]:headerFormat_cumsum[9]].strip().replace('D','E')); # sec/sec^2, SV clock drift rate
                
                # --- Per-system Offset ---
                quantz_cntr = 3; # Start after the header
                
            else:
                raise Exception('ERROR in rinexr: Unsupported sat type of '+satz+' found.');
            # END IF
            
            # --- DATA DATA ---
            for j in range(0, len(dataLines)):
                for jk in range(0, len(dataFormat_cumsum)-1):
                    if( quantz[quantz_cntr] != 'BLANK' ): # Ignores intentional blank spaces
                        val2heck = dataLines[j][dataFormat_cumsum[jk]:dataFormat_cumsum[jk+1]].strip().replace('D','E'); # do it in an easy way 
                        if( val2heck != '' ):
                            rinx['nav'][satz][quantz[quantz_cntr]][satz_filler[satz]] = np.float64(val2heck); # Convert if not empty
                        # END IF
                    # END IF
                    quantz_cntr += 1; # Increment the quantity counter
                    if( quantz_cntr >= len(quantz) ):
                        break; # Ditch it
                    # END IF
                # END FOR jk
            # END FOR j
            satz_filler[satz] += 1; # Record this sat record happened
        # END FOR i
        
        # --- PREVENT WEIRD TIME-STAMPS ---
        for satz in satz_cntr.keys():
            if( satz_cntr[satz] > 1 ): # Work only needed if more than 1 record
                # - Catch Duplicate Times -
                _, uqIdx, uqCnt = np.unique( rinx['nav'][satz]['time'][::-1], return_index=True, return_counts=True ); # Indexing is flipped to keep last unique element, not first since last is most recent if duplicates
                if( np.any(uqCnt > 1) ):
                    for keyz in rinx['nav'][satz]:
                        rinx['nav'][satz][keyz] = rinx['nav'][satz][keyz][::-1][uqIdx]; # Unique as needed [keep last instance b/c more relevant probably]
                    # END FOR keyz
                # END IF
                
                # - Catch Unsorted Times -
                k = np.argsort(rinx['nav'][satz]['time']); # Sort it up
                if( not np.array_equal(k,  np.arange(0, rinx['nav'][satz]['time'].size)) ):
                    for keyz in rinx['nav'][satz]:
                        rinx['nav'][satz][keyz] = rinx['nav'][satz][keyz][k]; # Sort as needed
                    # END FOR keyz
                # END IF
                
                # - Catch Colliding Times (within 5 mintues) -
                timeDiff_matrix = (np.tile(rinx['nav'][satz]['time'], (rinx['nav'][satz]['time'].size, 1)).T - np.tile(rinx['nav'][satz]['time'], (rinx['nav'][satz]['time'].size, 1)))/np.timedelta64(5, 'm'); # The "within 5 minutes" is enforced here
                ditchr = np.ones( rinx['nav'][satz]['time'].size, dtype=np.bool_ ); # prep
                for i in range(0, rinx['nav'][satz]['time'].size):
                    k = (np.abs(timeDiff_matrix[i, i:]) < 1) & (timeDiff_matrix[i, i:] != 0); # Compare
                    if( np.any(k) ):
                        k = np.where(k)[0]; # Get em
                        for j in range(0, k.size):
                            if( rinx['nav'][satz]['trans time'][i] < rinx['nav'][satz]['trans time'][i+k[j]] ): # Choose the one that was transmitted earlier to ditch
                                ditchr[i+k[j]] = False; # Ditch it
                            else:
                                ditchr[i] = False; # Ditch it
                            # END IF
                        # END FOR j
                    # END IF
                # END FOR i
                if( np.any(ditchr == False) ):
                    for keyz in rinx['nav'][satz]:
                        rinx['nav'][satz][keyz] = rinx['nav'][satz][keyz][ditchr]; # Ditch as needed
                    # END FOR keyz
                # END IF
            # END IF
        # END FOR satz
    else:
        rinx['NO DATA'] = 'NO DATA'; # No data, no record
    # END IF
        
    return rinx
# END DEF

#!!! NAVIGATION FUNCTION READERS
def rinexr_navigator_celestial2cartesian( timestamps, rinx_nav_gnss, sat2sat ):    
    # Implements Table 20-IV of https://www.gps.gov/technical/icwg/IS-GPS-200M.pdf
    ephemeri = rinx_nav_gnss['nav'][sat2sat]['TOE'].size; # Get how many ephemeri have been provided
    rinx_navd = {
        'pos':np.empty((ephemeri, 3), dtype=np.float64),
        'vel':np.empty((ephemeri, 3), dtype=np.float64),
        'acc':np.empty((ephemeri, 3), dtype=np.float64),
        }; # Prep a dict
    for i in range(0, ephemeri):
        a = 6378137.; # m, WGS84 equatorial radius
        mew = 3.986005*1E14; # meters^3/sec^2, WGS 84 value of the earth's gravitational constant for GPS user
        rotRate = 7.2921151467*1E-5; # rad/sec WGS 84 value of the earth's rotation rate
        J2 = 0.0010826262; # oblate Earth gravity coefficient
        semiA = rinx_nav_gnss['nav'][sat2sat]['sqrtA'][i]**2; # m, semi major axis
        meanMotion =  np.sqrt(mew/semiA**3) + rinx_nav_gnss['nav'][sat2sat]['deltaN'][i]; # rad/s, corrected mean motion
        tk = rinx_nav_gnss['nav'][sat2sat]['trans time'][i] - rinx_nav_gnss['nav'][sat2sat]['TOE'][i]; # sec of GPS week, time between transmission time and ephemeris reference epoch
        tk = 0; # Above translates from when the message is relevant for to when the message was sent; wellllll I'm using the time I scrape which is when the message is for. So IGNORE that, set this to 0
        if tk > 302400: tk -= 604800; # Enforce GPS week
        if tk < -302400: tk += 604800; # Enforce GPS week
        meanAnomaly = rinx_nav_gnss['nav'][sat2sat]['Mo'][i] + meanMotion*tk; # rad, mean anomaly
        def eccentricAnomalyFinder(eccentricAnomaly, meanAnomaly, eccentricity): # Solve function to solve            
            return eccentricAnomaly + (meanAnomaly - eccentricAnomaly + eccentricity*np.sin(eccentricAnomaly))/(1 - eccentricity*np.cos(eccentricAnomaly))
        # END DEF
        eccentricAnomaly = meanAnomaly
        for jj in range(0, 5): # iterative solve
            eccentricAnomaly = eccentricAnomaly + (meanAnomaly - eccentricAnomaly + rinx_nav_gnss['nav'][sat2sat]['ecc'][i]*np.sin(eccentricAnomaly))/(1 - rinx_nav_gnss['nav'][sat2sat]['ecc'][i]*np.cos(eccentricAnomaly)); # Just loop it 5 times or whatever - remember, divisor is the derivative b/c newton-raphson method
        # END FOR jj        
        trueAnomaly = 2*np.arctan( np.sqrt((1+rinx_nav_gnss['nav'][sat2sat]['ecc'][i])/(1-rinx_nav_gnss['nav'][sat2sat]['ecc'][i]))*np.tan(eccentricAnomaly/2)); # rad, true anomaly
        argOfLat_uncorr = trueAnomaly + rinx_nav_gnss['nav'][sat2sat]['AOP'][i]; # rad, uncorrected argument of latitude; AOP == Omega (prob)
        corr_lat = rinx_nav_gnss['nav'][sat2sat]['CUS'][i]*np.sin(2*argOfLat_uncorr) + rinx_nav_gnss['nav'][sat2sat]['CUC'][i]*np.cos(2*argOfLat_uncorr); # rad, arg of latitude correction
        corr_inc = rinx_nav_gnss['nav'][sat2sat]['CIS'][i]*np.sin(2*argOfLat_uncorr) + rinx_nav_gnss['nav'][sat2sat]['CIC'][i]*np.cos(2*argOfLat_uncorr); # rad, inclination correction
        corr_rad = rinx_nav_gnss['nav'][sat2sat]['CRS'][i]*np.sin(2*argOfLat_uncorr) + rinx_nav_gnss['nav'][sat2sat]['CRC'][i]*np.cos(2*argOfLat_uncorr); # m, radius correction
        argOfLat = argOfLat_uncorr + corr_lat; # rad, corrected argument of latitude
        incl = rinx_nav_gnss['nav'][sat2sat]['Io'][i] + corr_inc + rinx_nav_gnss['nav'][sat2sat]['Io rate'][i]*tk; # rad, corrected inclination
        radius = semiA*(1 - rinx_nav_gnss['nav'][sat2sat]['ecc'][i]*np.cos(eccentricAnomaly)) + corr_rad; # m, corrected radius
        xPos_orbit = radius*np.cos(argOfLat); # m, x position in orbital plane
        yPos_orbit = radius*np.sin(argOfLat); # m, y position in orbital plane
        RAAN = rinx_nav_gnss['nav'][sat2sat]['RAAN'][i] + (rinx_nav_gnss['nav'][sat2sat]['RAAN rate'][i] - rotRate)*tk - rotRate*rinx_nav_gnss['nav'][sat2sat]['TOE'][i]; # rad, corrected RAAN right ascension of the ascending node !! this may not be RAAN !!
        rinx_navd['pos'][i, 0] = xPos_orbit*np.cos(RAAN) - yPos_orbit*np.cos(incl)*np.sin(RAAN); # m, ECEF x pos
        rinx_navd['pos'][i, 1] = xPos_orbit*np.sin(RAAN) + yPos_orbit*np.cos(incl)*np.cos(RAAN); # m, ECEF y pos
        rinx_navd['pos'][i, 2] = yPos_orbit*np.sin(incl); # m, ECEF z pos
        # now for velocity
        eccentricAnomalyRate = meanMotion/(1 - rinx_nav_gnss['nav'][sat2sat]['ecc'][i]*np.cos(eccentricAnomaly)); # rad/s, eccentric anomaly rate
        trueAnomalyRate = eccentricAnomalyRate*np.sqrt(1 - rinx_nav_gnss['nav'][sat2sat]['ecc'][i]**2)/(1 - rinx_nav_gnss['nav'][sat2sat]['ecc'][i]*np.cos(eccentricAnomaly)); # rad/s, true anomaly rate
        inclRate = rinx_nav_gnss['nav'][sat2sat]['Io rate'][i] + 2*trueAnomalyRate*(rinx_nav_gnss['nav'][sat2sat]['CIS'][i]*np.cos(2*argOfLat_uncorr) - rinx_nav_gnss['nav'][sat2sat]['CIC'][i]*np.sin(2*argOfLat_uncorr)); # rad/s corrected inclination rate
        argOfLatRate = trueAnomalyRate + 2*trueAnomalyRate*(rinx_nav_gnss['nav'][sat2sat]['CUS'][i]*np.cos(2*argOfLat_uncorr) - rinx_nav_gnss['nav'][sat2sat]['CUC'][i]*np.sin(2*argOfLat_uncorr)); # rad/s, corrected argument of latitude rate
        radiusRate = rinx_nav_gnss['nav'][sat2sat]['ecc'][i]*semiA*eccentricAnomalyRate*np.sin(eccentricAnomaly) + 2*trueAnomalyRate*(rinx_nav_gnss['nav'][sat2sat]['CRS'][i]*np.cos(2*argOfLat_uncorr) - rinx_nav_gnss['nav'][sat2sat]['CRC'][i]*np.sin(2*argOfLat_uncorr)); # m/s, corrected radius rate
        RAANrate = rinx_nav_gnss['nav'][sat2sat]['RAAN rate'][i] - rotRate; # rad/s, RAAN rate
        xVel_orbit = radiusRate*np.cos(argOfLat) - radius*argOfLatRate*np.sin(argOfLat); # m/s, in-plane x velocity
        yVel_orbit = radiusRate*np.sin(argOfLat) + radius*argOfLatRate*np.cos(argOfLat); # m/s, in-plane x velocity
        rinx_navd['vel'][i, 0] = -xPos_orbit*RAANrate*np.sin(RAAN) + xVel_orbit*np.cos(RAAN) - yVel_orbit*np.sin(RAAN)*np.cos(incl) - yPos_orbit*(RAANrate*np.cos(RAAN)*np.cos(incl) - inclRate*np.sin(RAAN)*np.sin(incl)); # m/s, ECEF x component of velocity
        rinx_navd['vel'][i, 1] = xPos_orbit*RAANrate*np.cos(RAAN) + xVel_orbit*np.sin(RAAN) + yVel_orbit*np.cos(RAAN)*np.cos(incl) - yPos_orbit*(RAANrate*np.sin(RAAN)*np.cos(incl) + inclRate*np.cos(RAAN)*np.sin(incl)); # m/s, ECEF y component of velocity
        rinx_navd['vel'][i, 2] = yVel_orbit*np.sin(incl) + yPos_orbit*inclRate*np.cos(incl); # m/s, ECEF z component of velocity
        np.linalg.norm(rinx_navd['vel'][i, :]/1000)
        np.sqrt(mew/semiA)
        # now for acceleration
        F = -(3*J2*mew*a**2)/(2*radius**4); # oblate Earth acceleration factor
        rinx_navd['acc'][i, 0] = -mew*(rinx_navd['pos'][i, 0]/radius**3) + F*((1 - 5*(rinx_navd['pos'][i, 2]/radius)**2)*(rinx_navd['pos'][i, 0]/radius)) + 2*rinx_navd['vel'][i, 1]*rotRate + rinx_navd['pos'][i, 0]*rotRate**2; # m/s^2, ECEF x component of acceleration
        rinx_navd['acc'][i, 1] = -mew*(rinx_navd['pos'][i, 1]/radius**3) + F*((1 - 5*(rinx_navd['pos'][i, 2]/radius)**2)*(rinx_navd['pos'][i, 1]/radius)) - 2*rinx_navd['vel'][i, 0]*rotRate + rinx_navd['pos'][i, 1]*rotRate**2; # m/s^2, ECEF y component of acceleration
        rinx_navd['acc'][i, 2] = -mew*(rinx_navd['pos'][i, 2]/radius**3) + F*((3 - 5*(rinx_navd['pos'][i, 2]/radius)**2)*(rinx_navd['pos'][i, 2]/radius)); # m/s^2, ECEF z component of acceleration
    # END FOR i
    # finally record the time this is relevant for
    rinx_navd['time'] = rinx_nav_gnss['nav'][sat2sat]['time']; # Record
    
    # if( sat2sat == 'G05' ):
    #     breakpoint()
    #     pass
    
    return rinx_navd
# END DEF

def rinexr_navigator_cartesian( timestamps, rinx_nav_gnss, sat2sat ):
    # Easy since cartesian already
    ephemeri = rinx_nav_gnss['nav'][sat2sat]['trans time'].size; # Get how many ephemeri have been provided
    rinx_navd = {
        'pos':np.empty((ephemeri, 3), dtype=np.float64),
        'vel':np.empty((ephemeri, 3), dtype=np.float64),
        'acc':np.empty((ephemeri, 3), dtype=np.float64),
    }; # Prep a dict
    # Drop it in
    for i in range(0, ephemeri):
        rinx_navd['pos'][i, 0] = rinx_nav_gnss['nav'][sat2sat]['Xpos'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['pos'][i, 1] = rinx_nav_gnss['nav'][sat2sat]['Ypos'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['pos'][i, 2] = rinx_nav_gnss['nav'][sat2sat]['Zpos'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['vel'][i, 0] = rinx_nav_gnss['nav'][sat2sat]['Xvel'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['vel'][i, 1] = rinx_nav_gnss['nav'][sat2sat]['Yvel'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['vel'][i, 2] = rinx_nav_gnss['nav'][sat2sat]['Zvel'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['acc'][i, 0] = rinx_nav_gnss['nav'][sat2sat]['Xacc'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['acc'][i, 1] = rinx_nav_gnss['nav'][sat2sat]['Yacc'][i]*1000; # m, load it directly, convert to meters
        rinx_navd['acc'][i, 2] = rinx_nav_gnss['nav'][sat2sat]['Zacc'][i]*1000; # m, load it directly, convert to meters
    # END FOR i
    # finally record the time this is relevant for
    rinx_navd['time'] = rinx_nav_gnss['nav'][sat2sat]['time']; # Record
    
    return rinx_navd
# END DEF

def orbitPropagator_ODE(t, X, Re, mew, we, J02, xDD_xtra, yDD_xtra, zDD_xtra):
    # From grad school snooping around EOMs, idr much anymore and don't ahve the hrs to brush up - just rock banging rn
    # I have J3 info somewhere, if needed, find it me
    dX = np.empty( 6, dtype=np.float64); # Preallocate
    
    #X[0] = x
    #X[1] = y
    #X[2] = z
    #X[3] = xDot
    #X[4] = yDot
    #X[5] = zDot
    
    r = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2); # Radial component
    
    dX[0] = X[3]; # xDot
    dX[1] = X[4]; # yDot
    dX[2] = X[5]; # zDot
    # dX[3]= -mew*X[0]/r**3 - 3/2*J02*mew*Re**2*X[0]/r**5*(1 - 5*X[2]**2/r**2) + we**2*X[0] + 2*we*dX[1] + xDD_xtra; # xDD
    # # there is an extra term Jxam + Jxas which are perturbations due to the moon/sun but holy they were happening
    # # http://gauss.gge.unb.ca/GLONASS.ICD-98.pdf pg 39
    # # decided to get around the Jxam + Jxas by solving for them (know everything else at t0, pos, vel, acc, solve for them and assume constant for whole time)
    # dX[4]= -mew*X[1]/r**3 - 3/2*J02*mew*Re**2*X[1]/r**5*(1 - 5*X[2]**2/r**2) + we**2*X[1] + 2*we*dX[0] + yDD_xtra; # yDD
    # dX[5]= -mew*X[2]/r**3 - 3/2*J02*mew*Re**2*X[2]/r**5*(1 - 5*X[2]**2/r**2) + zDD_xtra; # zDD
    # # https://web.archive.org/web/20161020203029/http://russianspacesystems.ru/wp-content/uploads/2016/08/ICD_GLONASS_eng_v5.1.pdf pg 55 is what this is based on w/ J02, minor changes
    
    dX[3] = -mew*X[0]/r**3 - 3/2*J02*mew*Re**2*X[0]/r**5*(1 - 5*X[2]**2/r**2) + xDD_xtra; # xDD
    # there is an extra term Jxam + Jxas which are perturbations due to the moon/sun but holy they were happening
    #http://gauss.gge.unb.ca/GLONASS.ICD-98.pdf pg 39
    # decided to get around the Jxam + Jxas by solving for them (know everything else at t0, pos, vel, acc, solve for them and assume constant for whole time)
    dX[4] = -mew*X[1]/r**3 - 3/2*J02*mew*Re**2*X[1]/r**5*(1 - 5*X[2]**2/r**2) + yDD_xtra; # yDD
    dX[5] = -mew*X[2]/r**3 - 3/2*J02*mew*Re**2*X[2]/r**5*(3 - 5*X[2]**2/r**2) + zDD_xtra; # zDD

    return dX
# END DEF


def rinexr_orbitPropagator( rinx, X0_all, dataRate, FLG_ECEF=True, FLG_testing=False ):
    # dataRate is in Hz
    
    Re = 6378136.3; # Earth's mean equatorial radius (m)
    mew = 3.986004415*10**14; # mewEarth (m^3/s^2)
    # mew = 3.986005*10**14; # m^3/s^2, Earth's universal gravitational parameter https://www.navcen.uscg.gov/pubs/gps/icd200/ICD200Cw1234.pdf page 116 WSG 84 given
    # mew_km = mew/1000**3; # km^3/s^2, Earth's universal gravitational parameter
    we = 7.2921151467*10**-5; # rad/s, Earth rotation rate https://www.navcen.uscg.gov/pubs/gps/icd200/ICD200Cw1234.pdf page 116 WSG 84 given
    J02 = 1.082626925638815*10**-3; # J2
    
    # --- Remove super old ephemeri that will cause issues ---
    k_s = np.abs(X0_all['time'] - rinx[1][0])/np.timedelta64(1, 's'); # Get time between start and ephemeri
    k_e = np.abs(X0_all['time'] - rinx[1][-1])/np.timedelta64(1, 's'); # Get time between end and ephemeri
    k = (k_s > 86400) & (k_e > 86400); # Remove based on if time delta is greater than a day
    if( not np.all(k) and np.any(k) ):
        from copy import deepcopy
        X0_all = deepcopy(X0_all);
        X0_all['time'] = X0_all['time'][~k]; # Only keep nearby data
        X0_all['pos'] = X0_all['pos'][~k, :]; # Only keep nearby data
        X0_all['vel'] = X0_all['vel'][~k, :]; # Only keep nearby data
        X0_all['acc'] = X0_all['acc'][~k, :]; # Only keep nearby data
    elif( np.all(k)):
        print('All are old, figure this out')
        breakpoint()
        pass
    # END IF
    
    if( FLG_ECEF == True ):
        from copy import deepcopy
        X0_all = deepcopy(X0_all);
        # Right now major dependency on Astropy (which does a very good job)
        # Look at https://web.archive.org/web/20250812192655/https://github.com/eribean/Geneci/blob/main/geneci/conversions.py for a paired down implementation to drop the Astropy req (maybe Astropy method == FLG_highPrecision)
        from astropy import coordinates as coord
        from astropy import units as u
        from astropy.time import Time
        astro_pos = X0_all['pos']*u.m;
        astro_vel = X0_all['vel']*u.m/u.s;
        astro_time = Time(X0_all['time']);
        
        astro_ecef = coord.ITRS(x=astro_pos[:, 0], y=astro_pos[:, 1], z=astro_pos[:, 2], v_x=astro_vel[:, 0], v_y=astro_vel[:, 1], v_z=astro_vel[:, 2], representation_type='cartesian', differential_type='cartesian', obstime=astro_time); # Generate ITRS (ECEF) astropy coordinates
        astro_eci = astro_ecef.transform_to(coord.GCRS(obstime=astro_time)); # Convert from ITRS (ECEF) to GCRS (ECI)
        
        X0_all['pos'] = astro_eci.cartesian.xyz.T.value; # m, convert to ECI
        X0_all['vel'] = astro_eci.cartesian.differentials['s'].d_xyz.T.to(u.m/u.s).value; # m, convert to ECI (for some reason these become km/s by default)
    # END IF 
    
    X0 = np.empty( 6, dtype=np.float64); # PRep
    # rinx2tinx_min = np.where( np.min(np.abs(rinx[1][0] - X0_all['time'])) == np.abs(rinx[1][0] - X0_all['time']) )[0].item(); # Get the time
    rinx2tinx_min = np.where( rinx[1][0] < X0_all['time'] )[0]; # Get the time min (seems important to propagate from a 0 time, so backpropagating preferable to propagating from N time)
    if( rinx2tinx_min.size != 0 ):
        rinx2tinx_min = rinx2tinx_min[0]; # Get the first one
    else:
        rinx2tinx_min = 0; # First one
    # END IF
    
    rinx2tinx_max = np.where( rinx[1][-1] < X0_all['time'] )[0]; # Get the time max
    if( rinx2tinx_max.size != 0 ):
        rinx2tinx_max = rinx2tinx_max[0]; # Get the first one
    else:
        rinx2tinx_max = X0_all['time'].size; # Max it out
    # END IF
    
    # --- For Testing, Extend to last relevant ephemeri ---
    if( (FLG_testing == True) and np.all(rinx[1][-1] < X0_all['time']) ):
        stepy = np.median(np.diff(rinx[1]));
        rinx[1] = np.append(rinx[1], np.arange(rinx[1][-1] + stepy, X0_all['time'][0]+stepy, stepy)); # Extend it
    # END IF

    tinx_ref = X0_all['time'][rinx2tinx_min:rinx2tinx_max];
    tinx_splits = np.append( np.searchsorted(rinx[1], tinx_ref), rinx[1].size );
    if( tinx_splits[0] != 0 ):
        tinx_splits = np.insert(tinx_splits, 0, 0);
    elif( np.sum(tinx_splits == 0) > 1 ):
        tinx_splits = tinx_splits[np.where(tinx_splits == 0)[0][-1]:]; # Ditch extra 0's, due to a bunch of nav long before this obs
    # END IF
    tinx_diffCheck = np.where(np.diff(tinx_splits) == 1)[0]; # Catches a 1-size-step which doesn't work b/c it's 1. Soln: fold that 1-size-step into the previous step
    if( tinx_diffCheck.size > 0 ):
        for i in range(0, tinx_diffCheck.size):
            tinx_splits = np.hstack( (tinx_splits[:tinx_diffCheck[i]], tinx_splits[tinx_diffCheck[i]+1:]) ); # Remove the 1-size-step
    # END IF
    
    X = np.empty( (rinx[1].size, 6), dtype=np.float64); # Preallocate
    # Roll through most relevant ephemeri
    for i in range(0, tinx_splits.size-1):
        rinx2tinx_where = rinx[1][tinx_splits[i]] < tinx_ref; # Get the emphemeri to propagate from
        if( not np.all(rinx2tinx_where) ):
            rinx2tinx = np.where(~rinx2tinx_where)[0][-1] + rinx2tinx_min; # Get the first index and translate it to the full X0_all range with rinx2tinx_min
            FLG_backProp = False; # Backwards propagation off
        else:
            # If all are true, start time is less than all reference times, which means back propagation is required
            rinx2tinx = rinx2tinx_min; # If a negative indexer was requested, use the 1st one
            FLG_backProp = True; # Backwards propagation on
        # END IF
        
        # rinx2tinx = i + rinx2tinx_min; # Get indexer for the emphemeri
        # if( rinx2tinx < 0 ):
        #     rinx2tinx = 0; # If a negative indexer was requested, use the 1st one (0)
        #     FLG_backProp = True; # Backwards propagation on
        # else:
        #     FLG_backProp = False; # Backwards propagation off
        # # END IF
        
        # Get it out
        X0[0:3] = X0_all['pos'][rinx2tinx, :]; # Get it
        X0[3:6] = X0_all['vel'][rinx2tinx, :]; # Get it
        # if( FLG_zeroMode ):
        #     X0[5] *= -1 # Flip
        # X0[6:9] = X0_all['acc'][rinx2tinx, :]; # Get it
        
        # r = np.sqrt(X0[0]**2 + X0[1]**2 + X0[2]**2); # Radial component
        
        # # there is an extra term Jxam + Jxas which are perturbations due to the moon/sun but holy they were happening calcs
        # # http://gauss.gge.unb.ca/GLONASS.ICD-98.pdf pg 39
        # # decided to get around the Jxam + Jxas by solving for them (know everything else at t0, pos, vel, acc, solve for them and assume constant for whole time)
        # # xDD_xtra = navCorralled(j,2,3) + R_mew*X0[0]/R_r**3 - 3/2*R_C20*R_mew*R_ae**2*X0[0]/R_r**5*(1 - 5*X0[2]**2/R_r**2) - we**2*X0[0] - 2*we*X0[4]; #calc excess xDD accels
        # # yDD_xtra = navCorralled(j,3,3) + R_mew*X0[1]/R_r**3 - 3/2*R_C20*R_mew*R_ae**2*X0[1]/R_r**5*(1 - 5*X0[2]**2/R_r**2) - we**2*X0[1] - 2*we*X0[3];
        # # zDD_xtra = navCorralled(j,4,3) + R_mew*X0[2]/R_r**3 - 3/2*R_C20*R_mew*R_ae**2*X0[2]/R_r**5*(1 - 5*X0[2]**2/R_r**2);
        # # Reworked using J02 (same as C20 afaik actually you know) https://web.archive.org/web/20161020203029/http://russianspacesystems.ru/wp-content/uploads/2016/08/ICD_GLONASS_eng_v5.1.pdf pg 55
        # xDD_xtra = X0_all['acc'][rinx2tinx, 0] + mew*X0[0]/r**3 + 3/2*J02*mew*Re**2*X0[0]/r**5*(1 - 5*X0[2]**2/r**2) - we**2*X0[0] - 2*we*X0[4]; #calc excess xDD accels
        # yDD_xtra = X0_all['acc'][rinx2tinx, 1] + mew*X0[1]/r**3 + 3/2*J02*mew*Re**2*X0[1]/r**5*(1 - 5*X0[2]**2/r**2) - we**2*X0[1] - 2*we*X0[3];
        # zDD_xtra = X0_all['acc'][rinx2tinx, 2] + mew*X0[2]/r**3 + 3/2*J02*mew*Re**2*X0[2]/r**5*(1 - 5*X0[2]**2/r**2);
        # # reworked eqs to solve for what the effective extra perturbations due to Sun/Moon/Anything Else would be based on what is known at t0
        # # decided this was a horrible idea because the accels don't change dirs with the sat I guess. documents say these ests. are only good for +/-15 min anyways
            
        # [t,X] = ode45(@(t,X) rinex_sFUN_CALCodeR_Razzle(t,X,R_ae,R_mew,R_J02,0,0,0) ,tspan , X0, ODE_options); %solve ODE and such things
        xDD_xtra = 0;
        yDD_xtra = 0;
        zDD_xtra = 0;
        # xDD_xtra = X0_all['acc'][rinx2tinx, 0];
        # yDD_xtra = X0_all['acc'][rinx2tinx, 1];
        # zDD_xtra = X0_all['acc'][rinx2tinx, 2];
        
        tinx = (rinx[1][tinx_splits[i]:tinx_splits[i+1]] - X0_all['time'][rinx2tinx])/np.timedelta64(1, 's'); # Convert directly to seconds from whenever the data was taken
        if( not np.any(np.isclose(tinx, 0.)) ):
            tinx_ephemeralWhere = np.where(np.min(np.abs(tinx)) == np.abs(tinx))[0].item(); # Get the closest to 0
            if( len(tinx) > 1 ):
                tinx_ephemeralInsert = np.arange(tinx[tinx_ephemeralWhere] + np.median(np.abs(np.diff(tinx))), 0. + np.median(np.abs(np.diff(tinx))), np.median(np.abs(np.diff(tinx))) ); # What to insert
            else:
                tinx_ephemeralInsert = np.array( (0. ) ); # Just 0 if length is 1 - can't determine time spacing
            # END IF
            # tinx_ephemeralWhere = np.searchsorted( tinx, tinx_ephemeralInsert); # Search sorted to find where to be
            tinx = np.insert(tinx, tinx_ephemeralWhere + 1, tinx_ephemeralInsert ); # Insert as needed
            if( FLG_backProp ):
                tinx = np.flip(tinx); # Flip it, clip the front
                X[tinx_splits[i]:tinx_splits[i+1], :] = solve_ivp(orbitPropagator_ODE, (tinx[0], tinx[-1]), X0, args=(Re, mew, we, J02, xDD_xtra, yDD_xtra, zDD_xtra), t_eval = tinx, first_step=1/dataRate, rtol=1.0e-12, atol=1.0e-12).y.T[tinx_ephemeralInsert.size: :];
            else: # No flip, clip the back
                X[tinx_splits[i]:tinx_splits[i+1], :] = solve_ivp(orbitPropagator_ODE, (tinx[0], tinx[-1]), X0, args=(Re, mew, we, J02, xDD_xtra, yDD_xtra, zDD_xtra), t_eval = tinx, first_step=1/dataRate, rtol=1.0e-12, atol=1.0e-12).y.T[:tinx_ephemeralWhere+1, :];
            # END IF
        else:
            if( FLG_backProp ):
                tinx = np.flip(tinx); # Flip it
            # END IF
            X[tinx_splits[i]:tinx_splits[i+1], :] = solve_ivp(orbitPropagator_ODE, (tinx[0], tinx[-1]), X0, args=(Re, mew, we, J02, xDD_xtra, yDD_xtra, zDD_xtra), t_eval = tinx, first_step=1/dataRate, rtol=1.0e-12, atol=1.0e-12).y.T;
        # END IF
        
        if( FLG_backProp ):
            X[tinx_splits[i]:tinx_splits[i+1], :] = np.flipud(X[tinx_splits[i]:tinx_splits[i+1], :]); # Flip the positions about
        # END IF
    # END FOR i
    
    if( FLG_ECEF == True ):
        # Convert back
        astro_pos = X[:, 0:3]*u.m;
        astro_vel = X[:, 3:6]*u.m/u.s;
        astro_time = Time(rinx[1]);
        
        astro_eci = coord.GCRS(x=astro_pos[:, 0], y=astro_pos[:, 1], z=astro_pos[:, 2], v_x=astro_vel[:, 0], v_y=astro_vel[:, 1], v_z=astro_vel[:, 2], representation_type='cartesian', differential_type='cartesian', obstime=astro_time); # Generate GCRS (ECI) astropy coordinates
        astro_ecef = astro_eci.transform_to(coord.ITRS(obstime=astro_time)); # Convert from GCRS (ECI) to ITRS (ECEF)
        
        X[:, 0:3] = astro_ecef.cartesian.xyz.T.value; # m, convert to ECEF
        X[:, 3:6] = astro_ecef.cartesian.differentials['s'].d_xyz.T.to(u.m/u.s).value; # m, convert to ECEF (for some reason these become km/s by default)
    # END IF

    return X
# END DEF

#!!! NAVIGATION-FROM-OBS FUNCTION READER
def EOMs_wJ2_drag( time, X ): # X is an initial state vector
    NUMSTATE = int(np.round(1/2*(np.sqrt(4*X.size + 1) - 1))); # Solve for the number of state vectors
    dX = np.zeros( X.size);
    #X[1] = x
    #X[2] = y
    #X[3] = z
    #X[4] = xDot
    #X[5] = yDot
    #X[6] = zDot
    #X[7] = mewEarth
    #X[8] = J2
    #X[9] = CD drag coef.
    #X[10] = xSite1
    #X[11] = ySite1
    #X[12] = zSite1

    r = np.sqrt(X[0]**2+X[1]**2+X[2]**2);
    Re = 6378136.3; # Earth's mean equatorial radius (m)
    w = 2*np.pi/(86400); # Earth's rotation rate (rad/s)
    va = np.sqrt((X[3]+w*X[1])**2+(X[4]-w*X[0])**2+X[5]**2);
    rho0 = 3.614e-13; # reference density
    h = 88667.0; # [m] reference altitude
    r0 = 7e5+Re;# reference radius

    a=3.0; #m^2 s/c cross-sectional area
    m=970; #s/c mass [kg]
    # a = 8; # m^2, s/c cross-sectional area | Block IIA=3.5, IIR=5, IIF=6.75, III=8
    # m = 2269.323; # kg, s/c mass | Block II=1660, IIA=1816, IIR=1080, IIF=1465, III=2269 ! most are on-orbit weights, II/IIA prob aren't
    # Galileo mass https://www.gsc-europa.eu/support-to-developers/galileo-satellite-metadata

    rhoa = rho0*np.exp(-(r-r0)/h);
    k = 0.5*rhoa*(a/m);

    A = np.zeros( (NUMSTATE,NUMSTATE) );
    A[0,3] = 1;
    A[1,4] = 1;
    A[2,5] = 1;
    A[3,0] = -(X[6]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-1))+(3*X[6]*X[0]**2/r**5)*(1-2.5*X[7]*(Re/r)**2*\
             (7*(X[2]/r)**2-1))+k*va*X[8]*X[0]*(X[3]+w*X[1])/(r*h)-k*rhoa*X[8]*(-w*X[4]+w**2*X[0])*(X[3]+w*X[1])/va;
    A[3,1] = (3*X[6]*X[0]*X[1]/r**5)*(1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-1))+k*X[8]*va*X[1]*(X[3]+w*X[1])/\
        (r*h)-k*X[8]*(w*X[3]+w**2*X[1])*(X[3]+w*X[1])/va-k*X[8]*va*w;
    A[3,2] = (3*X[6]*X[0]*X[2]/r**5)*(1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-3))+k*X[8]*va*X[2]*(X[3]+w*X[1])/(r*h);
    A[3,3] = -k*X[8]*((X[3]+w*X[1])**2)/va-k*X[8]*va;
    A[3,4] = -k*X[8]*(X[4]-w*X[0])*(X[3]+w*X[1])/va;
    A[3,5] = -k*X[8]*X[5]*(X[3]+w*X[1])/va;
    A[3,6] = -(X[0]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-1));
    A[3,7] = (1.5*X[6]*X[0]/r**3)*(Re/r)**2*(5*(X[2]/r)**2-1);
    A[3,8] = -k*va*(X[3]+w*X[1]);
    A[4,0] = (3*X[6]*X[0]*X[1]/r**5)*(1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-1))+k*X[8]*va*X[0]*\
        (X[4]-w*X[0])/(r*h)-k*X[8]*(w**2*X[0]-w*X[4])*(X[4]-w*X[0])/va+k*X[8]*va*w;
    A[4,1] = -(X[6]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-1))+(3*X[6]*X[1]**2/r**5)*(1-2.5*X[7]*\
        (Re/r)**2*(7*(X[2]/r)**2-1))+k*X[8]*va*X[1]*(X[4]-w*X[0])/(r*h)-k*X[8]*(w*X[3]+w**2*X[1])*\
        (X[4]-X[0]*w)/va;
    A[4,2] = (3*X[6]*X[1]*X[2]/r**5)*(1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-3))+k*X[8]*va*X[2]*(X[4]-w*X[0])/(r*h);
    A[4,3] = -k*X[8]*(X[4]-w*X[0])*(X[3]+w*X[1])/va;
    A[4,4] = -k*X[8]*((X[4]-w*X[0])**2)/va-k*va*X[8];
    A[4,5] = -k*X[8]*X[5]*(X[4]-w*X[0])/va;
    A[4,6] = -(X[1]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-1));
    A[4,7] = (1.5*X[6]*X[1]/r**3)*(Re/r)**2*(5*(X[2]/r)**2-1);
    A[4,8] = -k*va*(X[4]-w*X[0]);

    A[5,0] = (3*X[6]*X[0]*X[2]/r**5)*(1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-3))+k*va*X[8]*X[5]*X[0]/\
        (r*h)-k*X[8]*X[5]*(w**2*X[0]-w*X[4])/va;
    A[5,1] = (3*X[6]*X[1]*X[2]/r**5)*(1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-3))+k*va*X[8]*X[5]*X[1]/\
        (r*h)-k*X[8]*X[5]*(w*X[3]-w**2*X[1])/va;
    A[5,2] = -(X[6]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-3))+(3*X[6]*X[2]**2/r**5)*\
        (1-2.5*X[7]*(Re/r)**2*(7*(X[2]/r)**2-5))+k*va*X[8]*X[2]*X[5]/(r*h);
    A[5,3] = -k*X[8]*X[5]*(X[3]+w*X[1])/va;
    A[5,4] = -k*X[8]*X[5]*(X[4]-w*X[0])/va;
    A[5,5] = -k*X[8]*X[5]**2/va-k*X[8]*va;
    A[5,6] = -(X[2]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-3));
    A[5,7] = (1.5*X[6]*X[2]/r**3)*(Re/r)**2*(5*(X[2]/r)**2-3);
    A[5,8] = -k*va*X[5];


    dX[0] = X[3]; #xDot
    dX[1] = X[4]; #yDot
    dX[2] = X[5]; #zDot
    dX[3] = -(X[6]*X[0]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-1)) \
        -0.5*X[8]*(a/m)*rho0*np.exp(-(r-r0)/h)*va*(X[3]+w*X[1]); #xDD
    dX[4] = -(X[6]*X[1]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-1)) \
        -0.5*X[8]*(a/m)*rho0*np.exp(-(r-r0)/h)*va*(X[4]-w*X[0]); #yDD
    dX[5] = -(X[6]*X[2]/r**3)*(1-1.5*X[7]*(Re/r)**2*(5*(X[2]/r)**2-3)) \
        -0.5*X[8]*(a/m)*rho0*np.exp(-(r-r0)/h)*va*(X[5]); #zDD

    # Converts the ICs into a matrix for multiplication (phidot=A*phi)
    # X_temp = np.empty( (NUMSTATE, NUMSTATE) );
    # for kk in range(0,NUMSTATE):
    #     for kk2 in range(0,NUMSTATE):
    #         X_temp[kk,kk2] = X[NUMSTATE*kk+kk2+NUMSTATE]; # This +NUMSTATE was added since this was wrong for sure
    #     # END FOR kk2
    # # END FOR kk
    X_temp = X[NUMSTATE:].reshape(NUMSTATE, NUMSTATE); # Vectorized

    # Performs the multiplication of phidot=A*phi
    with np.errstate(divide='ignore', invalid='ignore'):
        dX_temp = A@X_temp;
    # END WITH

    # # Puts phidot into the last 324 elements in the dX array
    # for kk in range(0,NUMSTATE):
    #     for kk2 in range(0,NUMSTATE):
    #         dX[kk*NUMSTATE+kk2+NUMSTATE] = dX_temp[kk,kk2]; # This +NUMSTATE was added since this was wrong for sure
    #     # END FOR kk2
    # # END FOR kk
    dX[NUMSTATE:] = dX_temp.ravel(); # Vectorized

    return dX
# END DEF


def funBatchHtilde(X,Xs,XsStation,rhoStar,rhoDotStar,theta,w,NUMSTATE):
    Htilde = np.zeros( (2,NUMSTATE) );

    #X[0] = x
    #X[1] = y
    #X[2] = z
    #X[3] = xDot
    #X[4] = yDot
    #X[5] = zDot

    #Xs[0] = x for station denoted by XsStation
    #Xs[1] = y for station denoted by XsStation
    #Xs[2] = z for station denoted by XsStation

    #XsStation = obs_posID[i] (0 to n)

    #rhoStar = range estimation that is nominal
    #rhoDotStar = range rate estimation that is nominal

    #theta = angle station is at to bring to Inertial Earth frame

    #Htilde[0,0] = drho/dx
    #Htilde[0,1] = drho/dy
    #Htilde[0,2] = drho/dz
    #Htilde[0,3] = 0
    #Htilde[0,4] = 0
    #Htilde[0,5] = 0
    #Htilde[0,6] = 0
    #Htilde[0,7] = 0
    #Htilde[0,8] = 0
    #Htilde[0,9] = drho/dxSite1
    #Htilde[0,10] = drho/dySite1
    #Htilde[0,11] = drho/dzSite1

    #Htilde[1,0] = drhoDot/dx
    #Htilde[1,1] = drhoDot/dy
    #Htilde[1,2] = drhoDot/dz
    #Htilde[1,3] = drhoDot/dxDot
    #Htilde[1,4] = drhoDot/dyDot
    #Htilde[1,5] = drhoDot/dzDot
    #Htilde[1,6] = 0
    #Htilde[1,7] = 0
    #Htilde[1,8] = 0
    #Htilde[1,9] = drhoDot/dxSite1
    #Htilde[1,10] = drhoDot/dySite1
    #Htilde[1,11] = drhoDot/dzSite1

    Htilde[0,0] = (X[0] - Xs[0]*np.cos(theta) + Xs[1]*np.sin(theta))/rhoStar;
    Htilde[0,1] = (X[1] - Xs[1]*np.cos(theta) - Xs[0]*np.sin(theta))/rhoStar;
    Htilde[0,2] = (X[2] - Xs[2])/rhoStar;
    #Only first 3 are non 0, until the station portion

    Htilde[0,3*XsStation+9] = (Xs[0] - X[0]*np.cos(theta) - X[1]*np.sin(theta))/rhoStar;
    Htilde[0,3*XsStation+10] = (Xs[1] - X[1]*np.cos(theta) + X[0]*np.sin(theta))/rhoStar;
    Htilde[0,3*XsStation+11] = (Xs[2] - X[2])/rhoStar;
    #This dynamically calculates the index for Htilde
    #So more stations can be added easily

    rhoStarSq = rhoStar**2; #So mult calcs are not needed
    Htilde[1,0] = (X[3] + Xs[1]*w*np.cos(theta) + Xs[0]*w*np.sin(theta))/rhoStar \
        - rhoDotStar*(X[0] - Xs[0]*np.cos(theta) + Xs[1]*np.sin(theta))/(rhoStarSq);
    Htilde[1,1] = (X[4] - Xs[0]*w*np.cos(theta) + Xs[1]*w*np.sin(theta))/rhoStar \
        - rhoDotStar*(X[1] - Xs[1]*np.cos(theta) - Xs[0]*np.sin(theta))/(rhoStarSq);
    Htilde[1,2] = X[5]/rhoStar - rhoDotStar*(X[2] - Xs[2])/rhoStarSq;

    Htilde[1,3] = (X[0] - Xs[0]*np.cos(theta) + Xs[1]*np.sin(theta))/rhoStar;
    Htilde[1,4] = (X[1] - Xs[1]*np.cos(theta) - Xs[0]*np.sin(theta))/rhoStar;
    Htilde[1,5] = (X[2] - Xs[2])/rhoStar;

    Htilde[1,3*XsStation+9] = (-X[3]*np.cos(theta) - X[4]*np.sin(theta) - \
        w*(X[1]*np.cos(theta) - X[0]*np.sin(theta)))/rhoStar - \
        rhoDotStar*(Xs[0] - X[0]*np.cos(theta) - X[1]*np.sin(theta))/rhoStarSq;

    Htilde[1,3*XsStation+10] = (X[3]*np.sin(theta) - X[4]*np.cos(theta) + \
        w*(X[0]*np.cos(theta) + X[1]*np.sin(theta)))/rhoStar - \
        rhoDotStar*(Xs[1] - X[1]*np.cos(theta) + X[0]*np.sin(theta))/rhoStarSq;

    Htilde[1,3*XsStation+11] = -(X[5]/rhoStar) - rhoDotStar*(Xs[2] - X[2])/rhoStarSq;

    #This dynamically calculates the index for Htilde

    #So more stations can be added easily

    return Htilde
# END DEF



def rinexr_navigator4obs( rinx, obs_pos, X0=None, obs_posID=None ):
    # This does navigation based on obs, no need for the navigation files
    # Only uses the range to determine satellite XYZ
    #AE 597C FINAL PROJECT BATCH PROCESSOR
    #THIS IS A BATCH PROCESSOR FOR BULK OBSERVATION DATA
    #'star' means nominal trajectory
    #'bar' is like something somewhere between the error and the state itself

    # Let's do this
    pinx = rinx[0]; # Get out the prange
    tinx = (rinx[1] - rinx[1][0])/np.timedelta64(1, 's'); # Convert directly to seconds
    # usefulDatetime = datetime.utcfromtimestamp(rinx[1][0].tolist() / 1e9);
    # tinx += usefulDatetime.hour*3600 + usefulDatetime.minute*60 + usefulDatetime.second + 86400/2; # Seconds it
    # Convert obs_pos to an initial guess (20,000 km overhead)
    if( obs_posID is None ):
        obs_posID = np.zeros( rinx[0].size, dtype=np.int32); # Make some obs station IDs
    # END IF
    X0_posGuess = np.empty( (3,3) )*np.nan; # Preallocate
    if( X0 is None ):
        X0_posGuess[0, :] = conv_aer2ecef(0, 20*np.pi/180, rinx[0][0], obs_pos[0][0], obs_pos[0][1], obs_pos[0][2], \
                                    FLG_radIn_azel=True, FLG_rangekm=False, FLG_radIn_latlong=True, \
                                    FLG_altkm=False, FLG_locECEF=True, FLG_kmOut=False, \
                                    FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=2)[0, :];
        X0_posGuess[1, :] = np.array( (1000., 1000., -1000.) ); # Random stuff
    else:
        # Find emphermi closest to rinx[1][0]
        rinx2tinx = np.where( np.min(np.abs(rinx[1][0] - X0['time'])) == np.abs(rinx[1][0] - X0['time']) )[0].item(); # Get the time
        # Apply that as X0_posGuess
        X0_posGuess[0, :] = X0['pos'][rinx2tinx, :]; # Get it
        X0_posGuess[1, :] = X0['vel'][rinx2tinx, :]; # Get it
        X0_posGuess[2, :] = X0['acc'][rinx2tinx, :]; # Get it
    # END IF

    # # --- REFERENCE DATA ---
    # # Read obs data
    # import pandas as pd
    # obsDat = pd.read_csv('obsdat.txt', sep=r'\s+', names=['time','obsID','prange','prangeRate'], header=None); # read data file
    # obsNum = obsDat.shape[0]; # find number of observations
    # obsType = obsDat.shape[1]-2; # find number of observations

    # # Try to hook up
    # pinx = obsDat['prange'].to_numpy();
    # tinx = obsDat['time'].to_numpy();
    # obs_pos = np.empty( (obsDat['obsID'].unique().size, 3), dtype=np.float64);
    # obs_pos[0, 0] = -5127510.0; # X(10) = xSite1 (101 - Pacific ship)
    # obs_pos[0, 1] = -3794160.0; # X(11) = ySite1 (m)
    # obs_pos[0, 2] = 0.0; # X(12) = zSite1
    # obs_pos[1, 0] = 3860910.0; # X(13) = xSite2 (337 - Turkey)
    # obs_pos[1, 1] = 3238490.0; # X(14) = ySite2 (m)
    # obs_pos[1, 2] = 3898094.0; # X(15) = zSite2
    # obs_pos[2, 0] = 549505.0; # X(16) = xSite3 (394 - Greenland)
    # obs_pos[2, 1] = -1380872.0; # X(17) = ySite3 (m)
    # obs_pos[2, 2] = 6182197.0; # X(18) = zSite3
    # _, obs_posID = np.unique(obsDat['obsID'].to_numpy(), return_inverse=True);
    
    # X0_posGuess = np.empty( (2, 3), dtype=np.float64); # PRep
    # X0_posGuess[0, 0] = 757700.0; # X(1) = x (m)
    # X0_posGuess[0, 1] = 5222607.0; # X(2) = y
    # X0_posGuess[0, 2] = 4851500.0; # X(3) = z
    # X0_posGuess[1, 0] = 2213.21; # X(4) = xDot (m/s)
    # X0_posGuess[1, 1] = 4678.34; # X(5) = yDot
    # X0_posGuess[1, 2] = -5371.30; # X(6) = zDot
    # # --- REFERENCE DATA ---

    # Prep Numbers
    NUMSTATE = 9+len(obs_pos)*3; #Number of states being solved for

    #***************Constants***************
    # mE = 5.97219*10**24; #kg, mass of Earth
    # aE = 6378.13649; #km, equatorial radius of Earth
    aE = 6378136.3; #m, equatorial radius of Earth
    # G = 6.67384*10**-11; #m^3/(kg*s^2), Gravitational Constant
    # G = G/1000**3; #kg^3/(kg*s^2), Grav Constant, diff units
    # w = 2*pi/(86400); #Earth's rotation rate (rad/s)
    w = 2*np.pi/(23.9344696*3600); #Earth's rotation rate, sidereal day (rad/s)

    # mew = G*mE; #m^3/s^2, gravitational constant of Earth
    # solarYr = 365.242; #days, solar year time
    # solarYrSec = solarYr*24*3600; #s, solar year time
    # J2 = 1.08264*10**-6; #unitless oblateness term, "Earth Grav Pot Derived from Sat Motion," Kozai 1966

    #Read obs data
    obsNum = pinx.size; # Number of observations
    obsType = 2; # Number of observation types are 1 (prange) and 2 (prange rate) -> we only provide prange, but this also estimates prange rate

    ## Initial Conditions Here
    #Building initial condition vector
    X0 = np.zeros( NUMSTATE, dtype=np.float64); #Preallocate for seriousness
    X0[0] = X0_posGuess[0, 0]; #X[1] = x (m)
    X0[1] = X0_posGuess[0, 1]; #X[2] = y
    X0[2] = X0_posGuess[0, 2]; #X[3] = z
    X0[3] = X0_posGuess[1, 0]; #X[4] = xDot (m/s)
    X0[4] = X0_posGuess[1, 1]; #X[5] = yDot
    X0[5] = X0_posGuess[1, 2]; #X[6] = zDot
    # X0[0] = 757700.0; # X(1) = x (m)
    # X0[1] = 5222607.0; # X(2) = y
    # X0[2] = 4851500.0; # X(3) = z
    # X0[3] = 2213.21; # X(4) = xDot (m/s)
    # X0[4] = 4678.34; # X(5) = yDot
    # X0[5] = -5371.30; # X(6) = zDot
    X0[6] = 3.986004415*10**14; #X[7] = mewEarth (m^3/s^2)
    X0[7] = 1.082626925638815*10**-3; #X[8] = J2
    X0[8] = 2.0; #X[9] = CD drag coef.
    for i in range(0, len(obs_pos)):
        X0[9+3*i] = obs_pos[i][0]; #X[10] = xSite1 (m)
        X0[10+3*i] = obs_pos[i][1]; #X[11] = ySite1 (m)
        X0[11+3*i] = obs_pos[i][2]; #X[12] = zSite1
    # END FOR i
    phi0 = np.eye(NUMSTATE).ravel(); #state transition matrix IC & flatten it out

    X0 = np.hstack( (X0, phi0) ); #put them together for final mega I.C.
    #X0 state matrix I.C. built

    X0bar = X0[0:NUMSTATE];
    P0bar = np.eye(NUMSTATE); # a priori covariance matrix
    P0bar[0,0] = 1*10**6;
    P0bar[1,1] = 1*10**6;
    P0bar[2,2] = 1*10**6;
    P0bar[3,3] = 1*10**6;
    P0bar[4,4] = 1*10**6;
    P0bar[5,5] = 1*10**6; #given in hints of problem on site
    P0bar[6,6] = 1*10**20;
    P0bar[7,7] = 1*10**6;
    P0bar[8,8] = 1*10**6;
    for i in range(0, len(obs_pos)):
        P0bar[9+3*i,9+3*i] = 1*10**6;
        P0bar[10+3*i,10+3*i] = 1*10**6;
        P0bar[11+3*i,11+3*i] = 1*10**6;
        # P0bar[9+3*i,9+3*i] = 1*10**-10;
        # P0bar[10+3*i,10+3*i] = 1*10**-10;
        # P0bar[11+3*i,11+3*i] = 1*10**-10;
    # END FOR i
    P0bar[9+3*0,9+3*0] = 1*10**-10;
    P0bar[10+3*0,10+3*0] = 1*10**-10;
    P0bar[11+3*0,11+3*0] = 1*10**-10;
    # P is NUMSTATE x NUMSTATE sized

    x0bar = np.zeros( NUMSTATE ); #a priori state deviation
    #All 0's, given in hints on site

    sigmaRho = 0.01; #m, error in range obs
    sigmaRhoDot = 0.001; #m/s, error in range rate obs
    R = np.array( ((sigmaRho**2,),) );
    #R is inverse of Weighting Matrix W, given in hints on site

    # sec = tinx[-1]; #sec, end time of the obs
    # tD = 20; #sec, time per step
    # tspan = np.arange(0, sec+tD, tD);
    tspan = tinx; # I want to know when the obs happen
    # options = odeset('RelTol', 1.0e-9,'InitialStep', 5.0e-2);
    # options=odeset('RelTol',1.0e-12,'AbsTol',1.0e-12, 'InitialStep', 1.0e0);
    #ODE solver prepared

    ## Prep for stuff to happen
    #PREALLOCATION
    tmp = np.zeros( NUMSTATE**2 ); #mid step for forming phi into square
    phi = np.zeros( (obsNum, NUMSTATE, NUMSTATE) ); #for state trans matrix
    rhoStar = np.zeros( obsNum ); #for nominal range est
    rhoDotStar = np.zeros( obsNum ); #for nominal range rate est
    # H = np.zeros( (NUMSTATE,NUMSTATE,obsNum) ); #for H matrix
    errOb = np.zeros( (NUMSTATE, obsNum) ); #for error of obs to expected obs
    y = np.zeros( (obsNum,obsType) ); #for obs error, NOT related to y direction

    # Do some batch
    loopLim = 10; #Number of iterations
    X0Holder = np.empty( (loopLim+1, NUMSTATE), dtype=X0.dtype)*np.nan; # Prep
    X0Holder[0, :] = X0[0:NUMSTATE]; #Holds original I.C.'s
    cntr = 0; rhoRMS = 1; # PRep 2 in a row
    while( (cntr < loopLim) and (rhoRMS > 0.01) ): #This loop iterates entire data set over and over
        print('Iteration '+str(cntr+1)+' \n');

        if( cntr == 0 ): #Initialization run
            #IF P0bar = 0, lambda + N = 0; IF x0bar = 0, N = 0 alone
            #IF no obs at t0, lambda = P0bar^-1, N = (P0bar^-1)*x0bar
            #IF obs at t0, lambda = P0^-1+H0'*W0*H0, N = H0'*W0*y0 +
            #(P0bar^-1)*x0bar (not sure about last bit)

            # lambda3 = eye(size(P0bar,1))/P0bar;
            # lambdaz = P0bar**-1;
            cholV = np.linalg.cholesky(P0bar).T; #Cholesky decomp for stability (match matlab needs a transpose)
            # ! Think about steps to improve np.linalg.inv(cholV) stability
            invChol = np.dot(np.eye(NUMSTATE), np.linalg.inv(cholV)); # inverse the Cholesky decomp of P0bar, equiv to np.eye(NUMSTATE)/cholV in MATLAB
            # Note np.dot(np.eye(NUMSTATE), np.linalg.pinv(cholV)) and np.linalg.lstsq(cholV, np.eye(NUMSTATE))[0] were NOT good enough
            #Matlab stability improvement, faster than ^-1
            lambdaz = invChol*invChol.T; #This way is faster than C'*C for some reason
            #Cholesky is faster than direct inverse too (entire thing)

            N = lambdaz@x0bar;
            #lambda is information matrix, defined as lambda = P^-1

        else: #After initialization run
            # lambdaz = P0bar**-1;
            cholV = np.linalg.cholesky(P0bar).T; #Cholesky decomp for stability (match matlab needs a transpose)
            # ! Think about steps to improve np.linalg.inv(cholV) stability
            invChol = np.dot(np.eye(NUMSTATE), np.linalg.inv(cholV)); # inverse the Cholesky decomp of P0bar, equiv to np.eye(NUMSTATE)/cholV in MATLAB
            # Note np.dot(np.eye(NUMSTATE), np.linalg.pinv(cholV)) and np.linalg.lstsq(cholV, np.eye(NUMSTATE))[0] were NOT good enough
            #Matlab stability improvement, faster than ^-1
            lambdaz = invChol@invChol.T; #This way is faster than C'*C for some reason

            N = lambdaz@x0bar;
        # END IF

        #Using ODE45 to solve function!
        X = solve_ivp(EOMs_wJ2_drag, (tspan[0], tspan[-1]), X0, t_eval = tspan, first_step=1.0e0, rtol=1.0e-12, atol=1.0e-12).y.T;
        # options=odeset('RelTol',1.0e-12,'AbsTol',1.0e-12, 'InitialStep', 1.0e0);
        # [t,X] = ode45(@(t,X) J2DragEOMs(t,X) ,tspan , X0, options);
        #A matrix is used in that function

        # for i in range(0, tspan.size): #THIS forms Phi into a square matrix
        #     for j in range( NUMSTATE, (NUMSTATE**2+NUMSTATE) ):
        #         tmp[j-NUMSTATE] = X[i, j]; #gets the phi STM (state trans matrix) bits
        #     # END FOR j
        #     # tmp = X[i, NUMSTATE:].ravel(); # Vectorized
        #     phi[i,:,:] = tmp.reshape(NUMSTATE,NUMSTATE).T; #reshapes in right way & puts in state trans matrix
        # # END FOR i
        phi = X[:, NUMSTATE:].reshape(tspan.size, NUMSTATE, NUMSTATE); # Gigavector

        for i in range(0, obsNum): #This loop runs through datas
            tL = tinx[i]; #Reads time for this time

            iState = np.where( tL == tspan)[0].item(); #This is the index for the state
            #There are more estimations than obs, so cheaper to run through obs
            #and just match estimation to the obs (less iterations than running
            #through estimations and matching to obs)

            Xs = X[iState, 9+3*obs_posID[i]:12+3*obs_posID[i]]; #This chooses the right station locale

            theta = w*tL; #Calculates theta (angle Earth is at from start)

            rhoStar[i] = np.sqrt( X[iState, 0]**2 + Xs[0]**2 + X[iState, 1]**2 + Xs[1]**2 \
                + X[iState, 2]**2 + Xs[2]**2 - 2*(X[iState,0]*Xs[0] + \
                X[iState, 1]*Xs[1])*np.cos(theta) + 2*(X[iState, 0]*Xs[1] - \
                X[iState, 1]*Xs[0])*np.sin(theta) - 2*X[iState, 2]*Xs[2] ); #est. range
            rhoDotStar[i] = (1/rhoStar[i])*( (X[iState, 0]*X[iState, 3]) + \
                (X[iState, 1]*X[iState, 4]) + (X[iState, 2]*X[iState, 5]) + \
                np.cos(theta)*(-Xs[0]*X[iState, 3] - Xs[1]*X[iState, 4] + \
                Xs[1]*X[iState, 0]*w - Xs[0]*X[iState, 1]*w) + \
                np.sin(theta)*(Xs[0]*X[iState, 0]*w + Xs[1]*X[iState, 1]*w + \
                Xs[1]*X[iState, 3] - Xs[0]*X[iState, 4]) - \
                X[iState, 5]*Xs[2] ); #est. range rate

            # Gstar = np.array( ((rhoStar[i], ), ) ); #Gstar holds estimated nominal obs

            # YL = pinx[i]; #Reads observations for this time (CAP Y)
            #YL OBS is 1 = roh, 2 = rohDot

            cholV = np.linalg.cholesky(R).T; #Cholesky decomp for stability (match matlab needs a transpose)
            # ! Think about steps to improve np.linalg.inv(cholV) stability
            invChol = np.dot(np.eye(R.shape[0]), np.linalg.inv(cholV)); # inverse the Cholesky decomp of P0bar, equiv to np.eye(NUMSTATE)/cholV in MATLAB
            # Note np.dot(np.eye(NUMSTATE), np.linalg.pinv(cholV)) and np.linalg.lstsq(cholV, np.eye(NUMSTATE))[0] were NOT good enough
            WL = invChol@invChol.T; #local weighting matrix

            # y[i,:] = YL - Gstar; #obs delta
            y[i,:] = pinx[i] - rhoStar[i]; #obs delta
            #y OBS delta is 1 = roh, 2 = rohdot

            Htilde = funBatchHtilde(X[iState,0:6],Xs,obs_posID[i],rhoStar[i],rhoDotStar[i],theta,w,NUMSTATE);
            #Htilde matrix that connects (through H) obs to predicted

            H = Htilde@phi[iState,:,:]; #calcs H matrix for this time period

            # lambdaz = lambdaz + H.T@WL@H; #calc'n that lamda
            lambdaz +=  H.T@H*WL.item(); #calc'n that lamda
            
            N += H.T@y[i,:]*WL.item(); #calc'n that N
            
            # try:
            cholV = np.linalg.cholesky(lambdaz).T; #Cholesky decomp for stability (match matlab needs a transpose)
            
            # ! Think about steps to improve np.linalg.inv(cholV) stability
            invChol = np.dot(np.eye(NUMSTATE), np.linalg.inv(cholV)); # inverse the Cholesky decomp of P0bar, equiv to np.eye(NUMSTATE)/cholV in MATLAB
            # Note np.dot(np.eye(NUMSTATE), np.linalg.pinv(cholV)) and 
            #     np.linalg.lstsq(cholV, np.eye(NUMSTATE))[0] were NOT good enough
            # cholV**-1;
            P = invChol@invChol.T; #Covariance matrix

            # P = lambdaz**-1; #By definition, inverse of knowledge matrix
            # except:
            #     lambdaz -=  H.T@H*WL.item(); # rollb ack
            #     N -= H.T@y[i,:]*WL.item(); # rollb ack
            # # END TRYING
        # END FOR i

        #x0delt error adjustment to make best guess
        x0delt = P@N; #x0delt
        #Original: = (lambda**-1)*N (already calc'd it)

        # errOb = y[i,:] - H@x0delt; #Error calc (error of the error)
        # RMS = (errOb.T@errOb*WL.item()/(y.shape[1]*y.shape[0]))**(1/2); #RMS

        # sigmax0 = sqrt(P(1,1)); #m, standard deviation x
        # sigmay0 = sqrt(P(2,2)); #m/s, standard deviation y
        # sigmaz0 = sqrt(P(3,3)); #m/s, standard deviation y
        # corrx0v0 = P(1,2)/(sigmax0*sigmav0); #correlation matrix

        rhoRMS = np.sqrt(np.mean(y[:,0]**2)) #m, this is of the adjusted y obs, not the original obs
        # rhoDotRMS = np.sqrt(np.mean(y[:,1]**2)) #m/s, this is of the adjusted y obs, not the orig obs

        X0[0:NUMSTATE] = X0[0:NUMSTATE] + x0delt; #nominal is now updated to best guess
        x0bar = x0bar + x0delt; #nominal error adjusted
        X0Holder[cntr+1, :] = X0[0:NUMSTATE]; #Records X0 adjustment

        from matplotlib import pyplot as plt
        print(rhoRMS)
        plt.figure()
        plt.plot(y[:, 0])
        plt.title(r'Residual (error) of Range ($\rho$) versus Observation Number')
        plt.xlabel('Obs. Number');
        plt.ylabel('Range Residual (m)');
        plt.show(block=True)

        cntr += 1; # Increment
    # END FOR cntr

    return X0Holder[cntr, :]
# END DEF


#!!! Support Functions
def parallel_starmap_helper(fn, args, kwargs=None):
    #for multiprocess starmap with or without kwargs
    # MUST be outside of the function that calls it or it just hangs
    if( kwargs is None ):
        return fn(*args)
    else:
        return fn(*args, **kwargs)
    # END IF
#END DEF

def getResults4Dict(dict_in, dict_out): # This one is in-place
    for keyz in dict_in:
        if( isinstance(dict_in[keyz], dict) ):
            if( keyz not in dict_out ):
                dict_out[keyz] = {}; # Prime it up
            # END IF
            getResults4Dict(dict_in[keyz], dict_out[keyz]); # Recurse
        elif( isinstance(dict_in[keyz], (np.ndarray, list, tuple)) ):
            if( keyz in dict_out ):
                dict_out[keyz].append(dict_in[keyz]); # Tack it on
            else:
                dict_out[keyz] = [dict_in[keyz]]; # Start buildin
            # END IF
        elif( isinstance(dict_in[keyz], (int, float, str, bytes)) ):
            dict_out[keyz] = dict_in[keyz]; # Copy over directly if flat, designed to overwrite if there's a conflict
        # END IF
    # END FOR keyz
# END DEF

def collapse4dict(dict_in): # This one is in-place
    for keyz in dict_in:
        if( isinstance(dict_in[keyz], dict) ):
            collapse4dict(dict_in[keyz]); # Recurse
        elif( isinstance(dict_in[keyz], list) ):
            if( isinstance(dict_in[keyz][0], np.ndarray) ):
                if( dict_in[keyz][0].ndim == 1 ):
                    dict_in[keyz] = np.hstack(dict_in[keyz]); # SMASH EM
                else:
                    dict_in[keyz] = np.vstack(dict_in[keyz]); # SMASH EM DIFFERENTLY
                # END IF
            else:
                # List should be a list of sats, so don't need to record differences
                FLG_samsies = True; # Assume true
                firstList = dict_in[keyz][0]; # Get first list
                for listy in dict_in[keyz][1:]: # Roll through the rest
                    if listy != firstList:
                        FLG_samsies = False;
                    # END IF
                # END FOR listy
                if( FLG_samsies == True ):
                    tmpy  = set([]); # Prep a set which will help with only capturing unique sats
                    for bits in dict_in[keyz]:
                        tmpy.update(bits); # Update the set
                    # END FOR keyz
                    dict_in[keyz] = sorted(list(tmpy)); # Convert back to list, sort it
                else:
                    dict_in[keyz] = sum(dict_in[keyz], []); # Put em all together in one big list
                # END IF
            # END IF
        # END IF
    # END FOR keyz
# END DEF