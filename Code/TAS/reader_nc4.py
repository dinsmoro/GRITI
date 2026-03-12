#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Reads NC4 files and shoves them into a dict, recursive powers are here as well

import netCDF4
from numpy import isclose, nan, ndarray, void, empty



def reader_nc4(nc4_file, nc4_path_prefix=None, diction=None, FLG_oneRead=False, FLG_singleListRemove=True, FLG_table2pandas=False):
    if( (nc4_path_prefix != None) and ('netCDF4.' in type(nc4_file).__module__) ):
        nc4_file = nc4_file[nc4_path_prefix]; #go to the prefix req'd
    #END IF
    if( diction == None ):
        diction = {}; #prepare dictionary
    #END IF
    if( '!!pandas!internal!!_df' in diction ):
        from pandas import read_hdf as pd_read_hdf
        diction = {}; # clear out the dictionary, this is a pandas datafile not one made by hand
        diction['df'] = pd_read_hdf(nc4_file); # Read in the hdf5 file using pandas which will turn it back into a dataframe (df)
    # END IF
    if( FLG_table2pandas ):
        from pandas import DataFrame
    # END IF
    
    if( 'netCDF4.' in type(nc4_file).__module__ ): # Only dig into it if it's an HDF5 module component
        #--- check for attributes 1st since they're not in with groups/datasets to preempt recursion ---
        attrs_dict = {}; # Build a dict for the whole thing
        for key in nc4_file.ncattrs(): # Yeehaw
            attrs_dict[key] = nc4_file.getncattr(key); #Get out the key, build the dict
        # END FOR key
        for key, item in attrs_dict.items():
            if( '!!internal!!_FWDSLASH' in key ):
                key = key.replace('!!internal!!_FWDSLASH','/'); # Restore the forward slash lost
            # END IF
            if( '!!internal!!_' not in key ):
                if( key in diction.keys() ):
                    if( isclose(diction[key], item) == False ):
                        #only worry is if the attribute isn't consistent
                        print('-----Warning in reader_nc4-----');
                        print('Attribute '+key+' isn\'t the same as the previously recorded value from another file of '+ \
                            str(diction[key])+' and this ('+str(nc4_file)+') file\'s value of '+str(item)+ \
                            '.\n NaN\'ing it and try to sort that out.');
                        diction[key] = nan; #nan that attribute, figure it out later
                    #END IF
                else:
                    if( isinstance(item, str) ):
                       diction[key] = item; #get that attribute out
                    elif( isinstance(item, bytes) ):
                        diction[key] = item.decode('UTF-8'); #get that attribute out, it's a byte string
                    elif( isinstance(item, ndarray) ):
                        if( key+'!!internal!!_dtype' in attrs_dict.keys() ): # An internal key on the dataset for special Numpy stuff (time)
                            diction[key] = item.astype(attrs_dict[key+'!!internal!!_dtype']); # Convert as requested
                        else:
                            if( item.ndim == 0 ):
                                diction[key] = item.item(); #get that attribute out, remove numpy whatever
                            else:
                                diction[key] = item; # It's an array as an attribute, cursed shoulda been a dataset
                            # END IF
                        # END IF
                    else:
                        if( key+'!!internal!!_dtype' in attrs_dict.keys() ): # An internal key on the dataset for special Numpy stuff (time)
                            diction[key] = item.astype(attrs_dict[key+'!!internal!!_dtype']); # Convert as requested
                        else:
                            diction[key] = item; # If it's something weird, just get it
                        # END IF
                    # END IF
                # END IF
            elif( '!!internal!!_None' in key ): # An internal key on the dataset for None entries
                diction[key[:key.find('!!internal!!_None')]] = None; # Rebuild the None as needed
            #END IF
            #END IF
        #END FOR key, item
        
        #--- check for everything else ---
        for key, item in nc4_file.variables.items():  
            if( '!!internal!!_FWDSLASH' in key ):
                key = key.replace('!!internal!!_FWDSLASH','/'); # Restore the forward slash lost
            # END IF
            if( key in diction.keys() ):
                if( isinstance(item[()], ndarray) and isinstance(item[()][0], void) ): # Rely on lazy python if statement comprehension
                    # It's a table!
                    arr2arr = item[()].copy(); # Copying it out greatly speeds up work - apparently it reads off disk otherwise (woof)
                    tmpNamez = arr2arr.dtype.names; # Get the names of the table headers
                    tmpDtypes = arr2arr.dtype; # Get the names of the table dtypes
                    if( len(tmpNamez) == 1 ):
                        # It's just text most likely
                        if( FLG_table2pandas == False ):
                            tmpStrang = [None for _ in range(0, len(arr2arr))]; # Preallocate
                            for i in range(0, len(arr2arr)):
                                if( ('bytes' == tmpDtypes[0].name[:5]) or ('object' == tmpDtypes[0].name) ): # Bytes means string to decode, hold it in a list
                                    tmpStrang[i] = arr2arr[i][0].decode('utf-8').strip(' '); # Decode the string
                                else:
                                    tmpStrang[i] = arr2arr[i][0]; # Get the thing
                                # END IF
                            # END FOR i
                            # diction[key].append('\n'.join(tmpStrang)); #get that string out after some work [combines, decided safer not to even if it's most likely text]
                            diction[key].append(tmpStrang); #get that string out after some work
                        else:
                            # Prep a Pandas dataframe to hold them
                            tmpDF = DataFrame(arr2arr);
                            for keyz in tmpDF.keys():
                                if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
                                    tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
                                # END IF
                            # END FOR keyz
                            diction[key].append(tmpDF); #get that dataset out after some work
                        # END IF
                    else:
                        # It's a table array
                        if( FLG_table2pandas == False ):
                            tmpDict = {}; # Prep a dict to hold them
                            tmpDict_decode = [False for _ in range(0, len(tmpNamez))]; # Prep
                            for j in range(0, len(tmpNamez)):
                                if( ('bytes' == tmpDtypes[j].name[:5]) or ('object' == tmpDtypes[j].name) ): # Bytes means string to decode, hold it in a list
                                    tmpDict[tmpNamez[j]] = [None for _ in range(0, len(arr2arr))]; # Hold it in a list
                                    tmpDict_decode[j] = True; # Set to decode
                                else:
                                    tmpDict[tmpNamez[j]] = empty( len(arr2arr), dtype=tmpDtypes[j]); # Otherwise it's a datatype to directly save
                                # END IF
                            # END FOR j
                            
                            for j in range(0, len(tmpNamez)):
                                if( tmpDict_decode[j] ):
                                    for i in range(0, len(arr2arr)):
                                        tmpDict[tmpNamez[j]][i] = arr2arr[i][j].decode('utf-8'); # Unpack into the dict after decoding the bytes to string
                                    # END FOR i
                                else:
                                    for i in range(0, len(arr2arr)):
                                        tmpDict[tmpNamez[j]][i] = arr2arr[i][j]; # Unpack into the dict directly
                                    # END FOR i
                                # END IF
                            # END FOR j
                            diction[key].append(tmpDict); #get that dataset out after some work
                        else:
                            # Prep a Pandas dataframe to hold them
                            tmpDF = DataFrame(arr2arr);
                            for keyz in tmpDF.keys():
                                if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
                                    tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
                                # END IF
                            # END FOR keyz
                            diction[key].append(tmpDF); #get that dataset out after some work
                        # END IF
                    # END IF
                else:
                    diction[key].append(item[()]); #tack that dataset on
                # END IF
            else:
                # item.attrs.item()
                #otherwise it's a new data type to add in
                if( isinstance(item[()], ndarray) and isinstance(item[()][0], void) ): # Rely on lazy python if statement comprehension
                    # It's likely a table!
                    arr2arr = item[()].copy(); # Copying it out greatly speeds up work - apparently it reads off disk otherwise (woof)
                    tmpNamez = arr2arr.dtype.names; # Get the names of the table headers
                    tmpDtypes = arr2arr.dtype; # Get the names of the table dtypes
                    if( len(tmpNamez) == 1 ):
                        # It's just text most likely
                        if( FLG_table2pandas == False ):
                            tmpStrang = [None for _ in range(0, len(arr2arr))]; # Preallocate
                            for i in range(0, len(arr2arr)):
                                if( ('bytes' == tmpDtypes[0].name[:5]) or ('object' == tmpDtypes[0].name) ): # Bytes means string to decode, hold it in a list
                                    tmpStrang[i] = arr2arr[i][0].decode('utf-8').strip(' '); # Decode the string
                                else:
                                    tmpStrang[i] = arr2arr[i][0]; # Get the thing
                                # END IF
                            # END FOR i
                            # diction[key] = ['\n'.join(tmpStrang)]; #get that string out after some work [combines, decided safer not to even if it's most likely text]
                            diction[key] = [tmpStrang]; #get that string out after some work
                        else:
                            # Prep a Pandas dataframe to hold them
                            tmpDF = DataFrame(arr2arr);
                            for keyz in tmpDF.keys():
                                if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
                                    tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
                                # END IF
                            # END FOR keyz
                            diction[key] = [tmpDF]; #get that dataset out after some work
                        # END IF
                    else:
                        # It's a table array
                        if( FLG_table2pandas == False ):
                            tmpDict = {}; # Prep a dict to hold them
                            tmpDict_decode = [False for _ in range(0, len(tmpNamez))]; # Prep
                            for j in range(0, len(tmpNamez)):
                                if( ('bytes' == tmpDtypes[j].name[:5]) or ('object' == tmpDtypes[j].name) ): # Bytes means string to decode, hold it in a list
                                    tmpDict[tmpNamez[j]] = [None for _ in range(0, len(arr2arr))]; # Hold it in a list
                                    tmpDict_decode[j] = True; # Set to decode
                                else:
                                    tmpDict[tmpNamez[j]] = empty( len(arr2arr), dtype=tmpDtypes[j]); # Otherwise it's a datatype to directly save
                                # END IF
                            # END FOR j
                            
                            for j in range(0, len(tmpNamez)):
                                if( tmpDict_decode[j] ):
                                    for i in range(0, len(arr2arr)):
                                        tmpDict[tmpNamez[j]][i] = arr2arr[i][j].decode('utf-8'); # Unpack into the dict after decoding the bytes to string
                                    # END FOR i
                                else:
                                    for i in range(0, len(arr2arr)):
                                        tmpDict[tmpNamez[j]][i] = arr2arr[i][j]; # Unpack into the dict directly
                                    # END FOR i
                                # END IF
                            # END FOR j
                            diction[key] = [tmpDict]; #get that dataset out after some work
                        else:
                            # Prep a Pandas dataframe to hold them
                            tmpDF = DataFrame(arr2arr);
                            for keyz in tmpDF.keys():
                                if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
                                    tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
                                # END IF
                            # END FOR keyz
                            diction[key] = [tmpDF]; #get that dataset out after some work
                        # END IF
                    # END IF
                else:
                    diction[key] = [item[()]]; #get that dataset out
                # END IF
            #END IF
            attrs_dict = {}; # Build a dict for the whole thing
            for key in item.ncattrs(): # Yeehaw
                attrs_dict[key] = item.getncattr(key); #Get out the key, build the dict
            # END FOR key
            for keyA, itemA in attrs_dict.items(): # Check attributes on the datset itself
                if( keyA == key+'!!internal!!_type' ): # An internal key on the dataset for Python built-ins
                    if( itemA == "<class 'list'>" ): # List
                        if( all(isinstance(val, bytes) for val in diction[key][-1]) ): # If it's all bytes, it's a list of strings
                            diction[key][-1] = [val.decode('UTF-8') for val in diction[key][-1] ]; # Convert to list while making each byte object into a string
                        else:
                            diction[key][-1] = list(diction[key][-1]); # Convert to list
                        # END IF
                    elif( itemA == "<class 'tuple'>" ):
                        if( all(isinstance(val, bytes) for val in diction[key][-1]) ): # If it's all bytes, it's a list of strings
                            diction[key][-1] = (val.decode('UTF-8') for val in diction[key][-1] ); # Convert to list while making each byte object into a string
                        else:
                            diction[key][-1] = tuple(diction[key][-1]); # Convert to list
                        # END IF
                    # END IF
                elif( keyA == key+'!!internal!!_dtype' ): # An internal key on the dataset for special Numpy stuff (time)
                    diction[key][-1] = diction[key][-1].astype(itemA); # Convert as requested
                # END IF
            # END FOR keyA, itemA
            if( FLG_oneRead ):
                diction[key] = diction[key][-1]; # Unlock the lock
            # END IF
        #END FOR key, item
        for key, item in nc4_file.groups.items():  
            if( '!!internal!!_FWDSLASH' in key ):
                key = key.replace('!!internal!!_FWDSLASH','/'); # Restore the forward slash lost
            # END IF
            
            if( key not in diction.keys() ):
                diction[key] = {}; #create a sub-dict if not already there
            #END IF
            if( any(['pandas_type' in tup[0] for tup in list(item.ncattrs())]) and any(['pandas_version' in tup[0] for tup in list(item.ncattrs())]) ):
                # This is actually a pandas, and we must eject to read it
                diction['!!pandas!internal!!_df'] = None;
                break
            else:
                diction[key] = reader_nc4(item, nc4_path_prefix=None, diction=diction[key], FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, FLG_table2pandas=FLG_table2pandas); #recurse
            # END IF
        #END IF
        
    # END IF
    
    if( FLG_singleListRemove ):
        # Roll through the diction keys and look for lists of length 1 and remove the list
        for keyz in diction.keys():
            if( isinstance(diction[keyz], (list, tuple)) ):
                diction[keyz] = diction[keyz][0]; # Get it out
            # END IF
        # END FOR keyz
    # END IF
    
    return diction
#END DEF
