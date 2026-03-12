#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Takes a dict and writes it to an h5py file, recursive powers are here as well
# h5_path_prefix is only on the top level, ignored for recursive dict checks
# keyz_ignore is only on the top level, ignored for recusive dict checks

from numpy import isscalar, ndarray
from h5py import File as h5pyFile
from os import path as ospath
# try:
#     from pandas import DataFrame
#     FLG_pd = True; # Use Pandas, keep behind a flag for lazy if statements
# except:
#     FLG_pd = False; # No Pandas, keep behind a flag for lazy if statements
# # END TRYING

def subfun_h5_writer_handler(h5_filePath, diction, h5_path_prefix=None, h5_path_prefix_create=False, keyz_ignore=None, topLevel_attributes=None, FLG_overWrite=False, FLG_append=False):
    # For backwards compatibility
    subfun_h5_writer(h5_filePath, diction, h5_path_prefix=h5_path_prefix, h5_path_prefix_create=h5_path_prefix_create, keyz_ignore=keyz_ignore, topLevel_attributes=topLevel_attributes, FLG_overWrite=FLG_overWrite, FLG_append=FLG_overWrite);
# END DEF

def subfun_h5_writer(h5_filePath, diction, h5_path_prefix=None, h5_path_prefix_create=False, keyz_ignore=None, topLevel_attributes=None, FLG_overWrite=False, FLG_append=False):
    if( (ospath.isfile(h5_filePath) == True) and (FLG_append == True) ):
        with h5pyFile( h5_filePath, 'a') as h5_file:
            # Write the dictionary to the file
            subfun_h5_writer_recusor(h5_file, diction, h5_path_prefix=h5_path_prefix, h5_path_prefix_create=h5_path_prefix_create, keyz_ignore=keyz_ignore);
            # Sprinkle on top level attributes if they're provided
            if( topLevel_attributes is not None ):
                for keyz in topLevel_attributes.keys():
                    h5_file.attrs[keyz] = topLevel_attributes[keyz]; #save the attributes in the topLevel_attributes var if it's here
                # END FOR keyz
            # END IF
        # END WITH
    # END IF
    if( (ospath.isfile(h5_filePath) == False) or (FLG_overWrite == True) ):
        with h5pyFile( h5_filePath, 'w') as h5_file:
            # Write the dictionary to the file
            subfun_h5_writer_recusor(h5_file, diction, h5_path_prefix=h5_path_prefix, h5_path_prefix_create=h5_path_prefix_create, keyz_ignore=keyz_ignore);
            # Sprinkle on top level attributes if they're provided
            if( topLevel_attributes is not None ):
                for keyz in topLevel_attributes.keys():
                    h5_file.attrs[keyz] = topLevel_attributes[keyz]; #save the attributes in the topLevel_attributes var if it's here
                # END FOR keyz
            # END IF
        # END WITH
    # END IF
# END DEF
        
def subfun_h5_writer_recusor(h5_file, diction, h5_path_prefix=None, h5_path_prefix_create=False, keyz_ignore=None):
    #--- prep work ---
    keyz = list(diction.keys()); #keys to the dict
    if( keyz_ignore != None ):
        if( (isinstance(keyz_ignore,list) == False) & (isinstance(keyz_ignore,tuple) == False) ):
            keyz_ignore = [keyz_ignore]; #wrap it
        elif( (isinstance(keyz_ignore,list) == False) & (isinstance(keyz_ignore,tuple) == True) ):
            keyz_ignore = list(keyz_ignore); #convert it in case we need to mess with it
        #END IF
        for keyz_per in keyz_ignore:
            keyz.remove(keyz_per); #ditch keyz requested to be ditched
        #END FOR keyz_per
    #END IF
    
    if( h5_path_prefix != None ):
        if( h5_path_prefix_create ):
            h5_file = h5_file.create_group(str(h5_path_prefix)); #create a handle for a subgroup, deal with non-string handles
        else:
            h5_file = h5_file[h5_path_prefix]; #go to the prefix req'd
        # END IF
    #END IF
    
    #--- load everything into the file ---
    for keyz_per in keyz:
        if( isscalar(diction[keyz_per]) == False ):
            if( isinstance(diction[keyz_per], dict) ):
                if( isinstance(keyz_per, str) ):
                    h5_file_subgrp = h5_file.create_group(keyz_per.replace('/','!!internal!!_FWDSLASH')); #create a handle for a subgroup
                else:
                    h5_file_subgrp = h5_file.create_group(str(keyz_per)); #create a handle for a subgroup, deal with non-string handles
                # END IF
                subfun_h5_writer_recusor(h5_file_subgrp, diction[keyz_per], h5_path_prefix=None, keyz_ignore=None); #recurse!
            elif( isinstance(diction[keyz_per], ndarray) ):
                if( diction[keyz_per].ndim == 0 ):
                    if( '<M8' in diction[keyz_per].__array_interface__['typestr']): # This is to support any sub-time unit [ns] or [ms] or [D] or whatever
                        # Special handling for numpy datetime64 objects
                        #if size 1, add it as an attribute
                        h5_file.attrs[keyz_per] = diction[keyz_per].item(); #save the attribute                        
                        h5_file.attrs[keyz_per+'!!internal!!_dtype'] = diction[keyz_per].__array_interface__['typestr']; #record the true string datatype
                    else:
                        # Regular
                        #if size 1, add it as an attribute
                        h5_file.attrs[keyz_per] = diction[keyz_per].item(); #save the attribute
                    # END IF
                else:
                    if( '<M8' in diction[keyz_per].__array_interface__['typestr']): # This is to support any sub-time unit [ns] or [ms] or [D] or whatever
                        # Special handling for numpy datetime64 objects
                        dataset = h5_file.create_dataset(keyz_per, data=diction[keyz_per].view('<i8'), chunks=diction[keyz_per].shape, compression="gzip", compression_opts=9, shuffle=True, fletcher32=True); #write that data
                        dataset.attrs[keyz_per+'!!internal!!_dtype'] = diction[keyz_per].__array_interface__['typestr']; #record the true string datatype
                    else:
                        # Regular
                        h5_file.create_dataset(keyz_per, data=diction[keyz_per], chunks=diction[keyz_per].shape, compression="gzip", compression_opts=9, shuffle=True, fletcher32=True); #write that data
                    # END IF
                # END IF
            elif( ('pandas' in str(type(diction[keyz_per]))) and ('DataFrame' in str(type(diction[keyz_per]))) ):
                # Pandas array, need to save as if it was a dictionary with metadata to recover as a pandas table later
                tmpDict = {}; # Prep
                for subkeyz in diction[keyz_per].keys():
                    tmpDict[subkeyz] = diction[keyz_per][subkeyz].to_numpy(); # Put the column into
                # END FOR subkeyz
                tmpPd = diction[keyz_per].copy(); # Keep a copy for later
                diction[keyz_per] = tmpDict; # Replace with the dict version
                h5_file_subgrp = h5_file.create_group(str(keyz_per)+'_!!pandas!internal!!_dict2df'); #create a handle for a subgroup
                subfun_h5_writer_recusor(h5_file_subgrp, diction[keyz_per], h5_path_prefix=None, keyz_ignore=None); #recurse!
                diction[keyz_per] = tmpPd.copy(); # Bring it back
                del tmpPd, tmpDict
            elif( diction[keyz_per] is None ):
                # None needs special handling
                h5_file.attrs[keyz_per+'!!internal!!_None'] = 'None'; # Record the true string datatype
            else: # Imagining this captures lists and tuples and maybe pandas will work too who knows not I
                # List did not accept into its heart compression and other things
                dataset = h5_file.create_dataset(keyz_per, data=diction[keyz_per], chunks=len(diction[keyz_per])); #write that data
                dataset.attrs[keyz_per+'!!internal!!_type'] = str(type(diction[keyz_per])); # Record the true string datatype
            #END IF
        else:
            #if size 1, add it as an attribute
            h5_file.attrs[keyz_per] = diction[keyz_per]; #save the attribute
        #END IF
    #END FOR keyz_per
#END DEF