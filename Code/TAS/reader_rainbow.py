#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os import path as ospath, remove as osremove
from shutil import copyfileobj as shutilcopyfileobj, unpack_archive as shutilunpack_archive
from glob import glob
import io, gzip, tarfile
import h5py, netCDF4
from re import search as research
from TAS.reader_h5 import reader_h5
from TAS.reader_nc4 import reader_nc4

# Unified HDF5 and NetCDF4 reader with built-in de-zipping in memory. (TODO: add unzip, read, delete unzipped as option)
def reader_rainbow(file2read, path_prefix=None, diction=None, FLG_oneRead=False, \
                   FLG_singleListRemove=True, FLG_table2pandas=False, FLG_dezipInMem=True):
    file2read_name = file2read; # Copy that name
    
    # Check if zipped, dezip in memory
    FLG_unzippt = False; # If true, deletes file2read which should be renamed to the unzipped file name
    typez = ospath.splitext(file2read_name)[1].lower().replace('.',''); # Lower it, remove periods
    file2read_rider = []; # Prep a list to hold the tar'd names IF it turns out to be tar'd
    if( (typez == 'tgz') or (research(r'.tar.gz$', file2read_name) is not None) ):
        if( FLG_dezipInMem ):
            file2read = []; # Prep a list, we don't know how big it will be
            with tarfile.open(file2read_name, 'r:gz') as finput:
                for tarz in finput: # Roll through the files in the tar file
                    file2read.append(finput.extractfile(tarz.name).read()); # Extract and read in-place, it's now BYTES
                    file2read_rider.append(tarz.name); # Get the tar'd name
                # END FOR tarz
            # END WITH
        else:
            tmp_fldr = file2read_name[:-4] if (typez == 'tgz') else file2read_name[:-7]; # Get the folder using a naughty one line if statement
            shutilunpack_archive(file2read_name, file2read_name[:-4] if (typez == 'tgz') else file2read_name[:-7]); # Unpack the files
            file2read = glob(tmp_fldr+'/*'); # Get the files dropped out
            file2read_rider = [ospath.basename(filez) for filez in file2read]; # Get the file basenames
            FLG_unzippt = True; # Delete the unzipped file at the end
        # END IF
        # We update to the latest typez later
    elif( typez in ['gz', 'gzip'] ):
        typez = ospath.splitext(ospath.splitext(file2read_name)[0])[1].lower().replace('.',''); # Get the real typez now
        if( FLG_dezipInMem ):
            # Replace the file2read path with a byte object that's the file in memory
            with gzip.open(file2read, 'rb') as finput:
                file2read = finput.read(); # Read in-place, it's now BYTES
            # END WITH
        else:
            with gzip.open(file2read, 'rb') as finput:
                with open(ospath.splitext(file2read)[0], 'wb') as foutput:
                    shutilcopyfileobj(finput, foutput); # Unzip and save it - not in memory but slower b/c disk and adds disk writes
                # END WITH
            # END WITH
            file2read = ospath.splitext(file2read)[0]; # Read the unzipped file directly now
            FLG_unzippt = True; # Delete the unzipped file at the end
        # END IF
    # END IF
    
    # Check if HDF5 or NetCDF4
    if( len(file2read_rider) == 0 ): # If file2read_rider has nothing, then there's no tarball multi-file action
        if( typez in ['h5', 'hdf5'] ):
            if( 'bytes' in str(type(file2read)) ):
                file2read = io.BytesIO(file2read); # Convert to bytesIO for h5py
            # END IF
            with h5py.File(file2read, 'r', rdcc_nbytes=500*1024*1024) as readingRainbow:
                diction = reader_h5(readingRainbow, h5_path_prefix=path_prefix, diction=diction, \
                                           FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, \
                                           FLG_table2pandas=FLG_table2pandas);
            # END WITH
        elif( (typez in ['nc', 'nc4', 'netcdf4']) or (research(r'[0-9]+_nc$', typez) is not None) ):
            if( "'str'" in str(type(file2read)) ):
                with netCDF4.Dataset(file2read, 'r') as readingRainbow:
                    diction = reader_nc4(readingRainbow, nc4_path_prefix=path_prefix, diction=diction, \
                                         FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, \
                                         FLG_table2pandas=FLG_table2pandas);
                # END WITH
            elif( 'bytes' in str(type(file2read)) ):
                with netCDF4.Dataset('inmem.nc', 'r', memory=file2read) as readingRainbow:
                    diction = reader_nc4(readingRainbow, nc4_path_prefix=path_prefix, diction=diction, \
                                         FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, \
                                         FLG_table2pandas=FLG_table2pandas);
                # END WITH
            else:
                raise Exception('ERROR in reader_rainbow.py: File `'+file2read_name+'` is NOT known file format (`'+str(type(file2read))+'` provided).');
            # END IF
        else:
            raise Exception('ERROR in reader_rainbow.py: File `'+file2read_name+'` is NOT known file type (`'+typez+'` provided).');
        # END IF
    else:
        if( diction is None ):
            diction = {}; # Prep for this
        # END IF
        
        for i in range(0, len(file2read)):
            typez = ospath.splitext(file2read_rider[i])[1].lower().replace('.',''); # Get the real typez now
            
            if file2read_rider[i] not in diction:
                diction[file2read_rider[i]] = {}; # Prep this sub-key
            # END IF
            
            if( typez in ['h5', 'hdf5'] ):
                if( 'bytes' in str(type(file2read[i])) ):
                    file2read = io.BytesIO(file2read[i]); # Convert to bytesIO for h5py
                # END IF
                with h5py.File(file2read[i], 'r', rdcc_nbytes=500*1024*1024) as readingRainbow:
                    diction[file2read_rider[i]] = reader_h5(readingRainbow, h5_path_prefix=path_prefix, diction=diction[file2read_rider[i]], \
                                               FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, \
                                               FLG_table2pandas=FLG_table2pandas);
                # END WITH
            elif( (typez in ['nc', 'nc4', 'netcdf4']) or (research(r'[0-9]+_nc$', typez) is not None) ):
                if( "'str'" in str(type(file2read[i])) ):
                    with netCDF4.Dataset(file2read[i], 'r') as readingRainbow:
                        diction[file2read_rider[i]] = reader_nc4(readingRainbow, nc4_path_prefix=path_prefix, diction=diction[file2read_rider[i]], \
                                             FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, \
                                             FLG_table2pandas=FLG_table2pandas);
                    # END WITH
                elif( 'bytes' in str(type(file2read[i])) ):
                    with netCDF4.Dataset('inmem.nc', 'r', memory=file2read[i]) as readingRainbow:
                        diction[file2read_rider[i]] = reader_nc4(readingRainbow, nc4_path_prefix=path_prefix, diction=diction[file2read_rider[i]], \
                                             FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, \
                                             FLG_table2pandas=FLG_table2pandas);
                    # END WITH
                else:
                    raise Exception('ERROR in reader_rainbow.py: File `'+file2read_rider[i]+'` is NOT known file format (`'+str(type(file2read[i]))+'` provided).');
                # END IF
            else:
                raise Exception('ERROR in reader_rainbow.py: File `'+file2read_rider[i]+'` is NOT known file type (`'+typez+'` provided).');
            # END IF
        # END FOR i
    # END IF
    
    if( FLG_unzippt ):
        if( len(file2read_rider) == 0 ): # If file2read_rider has nothing, then there's no tarball multi-file action
            osremove(file2read); # Delete the unzippt file
        else:
            for i in range(0, len(file2read)):
                osremove(file2read[i]); # Delete the unzippt file
            # END FOR i
        # END IF
    # END IF
    
    return diction
# END DEF