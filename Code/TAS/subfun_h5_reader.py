# Reads h5py files and shoves them into a dict, recursive powers are here as well

from TAS.reader_h5 import reader_h5

def subfun_h5_reader(h5_file, h5_path_prefix=None, diction=None, FLG_oneRead=False, FLG_singleListRemove=True, FLG_table2pandas=False):
    return reader_h5(h5_file, h5_path_prefix=h5_path_prefix, diction=diction, FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, FLG_table2pandas=FLG_table2pandas); # Alias it over
# END DEF

# import h5py
# from numpy import isclose, nan, ndarray, void, empty

# def subfun_h5_reader(h5_file, h5_path_prefix=None, diction=None, FLG_oneRead=False, FLG_singleListRemove=True, FLG_table2pandas=False):
#     #--- prep work ---
#     if( 'h5py.' not in type(h5_file).__module__ ):
#     # if( not isinstance( h5_file, h5py._hl.files.File) ):
#     #if 'h5py' not in getattr(h5_file, '__module__', None): # Lazy but I hope it works
#         # Recurse with an opened h5_file if h5_file is a (hopefully) path
#         with h5py.File( h5_file, 'r') as h5_file_reader:
#             diction = subfun_h5_reader(h5_file_reader, h5_path_prefix=h5_path_prefix, diction=diction, FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, FLG_table2pandas=FLG_table2pandas)
#         # END WITH
#     # END IF
#     if( (h5_path_prefix != None) and ('h5py.' in type(h5_file).__module__) ):
#         h5_file = h5_file[h5_path_prefix]; #go to the prefix req'd
#     #END IF
#     if( diction == None ):
#         diction = {}; #prepare dictionary
#     #END IF
#     if( '!!pandas!internal!!_df' in diction ):
#         from pandas import read_hdf as pd_read_hdf
#         diction = {}; # clear out the dictionary, this is a pandas datafile not one made by hand
#         diction['df'] = pd_read_hdf(h5_file); # Read in the hdf5 file using pandas which will turn it back into a dataframe (df)
#     # END IF
#     if( FLG_table2pandas ):
#         from pandas import DataFrame
#     # END IF
    
#     if( 'h5py.' in type(h5_file).__module__ ): # Only dig into it if it's an HDF5 module component
#         #--- check for attributes 1st since they're not in with groups/datasets to preempt recursion ---
#         attrs_dict = dict(h5_file.attrs.items()); # Yeehaw
#         for key, item in h5_file.attrs.items():
#             if( '!!internal!!_FWDSLASH' in key ):
#                 key = key.replace('!!internal!!_FWDSLASH','/'); # Restore the forward slash lost
#             # END IF
#             if( '!!internal!!_' not in key ):
#                 if( key in diction.keys() ):
#                     if( isclose(diction[key], item) == False ):
#                         #only worry is if the attribute isn't consistent
#                         print('-----Warning in subfun_h5_reader-----');
#                         print('Attribute '+key+' isn\'t the same as the previously recorded value from another file of '+ \
#                             str(diction[key])+' and this ('+str(h5_file)+') file\'s value of '+str(item)+ \
#                             '.\n NaN\'ing it and try to sort that out.');
#                         diction[key] = nan; #nan that attribute, figure it out later
#                     #END IF
#                 else:
#                     if( isinstance(item, str) ):
#                        diction[key] = item; #get that attribute out
#                     elif( isinstance(item, bytes) ):
#                         diction[key] = item.decode('UTF-8'); #get that attribute out, it's a byte string
#                     elif( isinstance(item, ndarray) ):
#                         if( key+'!!internal!!_dtype' in attrs_dict.keys() ): # An internal key on the dataset for special Numpy stuff (time)
#                             diction[key] = item.astype(attrs_dict[key+'!!internal!!_dtype']); # Convert as requested
#                         else:
#                             if( item.ndim == 0 ):
#                                 diction[key] = item.item(); #get that attribute out, remove numpy whatever
#                             else:
#                                 diction[key] = item; # It's an array as an attribute, cursed shoulda been a dataset
#                             # END IF
#                         # END IF
#                     else:
#                         if( key+'!!internal!!_dtype' in attrs_dict.keys() ): # An internal key on the dataset for special Numpy stuff (time)
#                             diction[key] = item.astype(attrs_dict[key+'!!internal!!_dtype']); # Convert as requested
#                         else:
#                             diction[key] = item; # If it's something weird, just get it
#                         # END IF
#                     # END IF
#                 # END IF
#             elif( '!!internal!!_None' in key ): # An internal key on the dataset for None entries
#                 diction[key[:key.find('!!internal!!_None')]] = None; # Rebuild the None as needed
#             #END IF
#             #END IF
#         #END FOR key, item
        
#         #--- check for everything else ---
#         for key, item in h5_file.items():  
#             if( '!!internal!!_FWDSLASH' in key ):
#                 key = key.replace('!!internal!!_FWDSLASH','/'); # Restore the forward slash lost
#             # END IF
#             if( isinstance(item, h5py.Dataset) ): #dataset catch
#                 if( key in diction.keys() ):
#                     if( isinstance(item[()], ndarray) and isinstance(item[()][0], void) ): # Rely on lazy python if statement comprehension
#                         # It's a table!
#                         arr2arr = item[()].copy(); # Copying it out greatly speeds up work - apparently it reads off disk otherwise (woof)
#                         tmpNamez = arr2arr.dtype.names; # Get the names of the table headers
#                         tmpDtypes = arr2arr.dtype; # Get the names of the table dtypes
#                         if( len(tmpNamez) == 1 ):
#                             # It's just text most likely
#                             if( FLG_table2pandas == False ):
#                                 tmpStrang = [None for _ in range(0, len(arr2arr))]; # Preallocate
#                                 for i in range(0, len(arr2arr)):
#                                     if( ('bytes' == tmpDtypes[0].name[:5]) or ('object' == tmpDtypes[0].name) ): # Bytes means string to decode, hold it in a list
#                                         tmpStrang[i] = arr2arr[i][0].decode('utf-8').strip(' '); # Decode the string
#                                     else:
#                                         tmpStrang[i] = arr2arr[i][0]; # Get the thing
#                                     # END IF
#                                 # END FOR i
#                                 # diction[key].append('\n'.join(tmpStrang)); #get that string out after some work [combines, decided safer not to even if it's most likely text]
#                                 diction[key].append(tmpStrang); #get that string out after some work
#                             else:
#                                 # Prep a Pandas dataframe to hold them
#                                 tmpDF = DataFrame(arr2arr);
#                                 for keyz in tmpDF.keys():
#                                     if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
#                                         tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
#                                     # END IF
#                                 # END FOR keyz
#                                 diction[key].append(tmpDF); #get that dataset out after some work
#                             # END IF
#                         else:
#                             # It's a table array
#                             if( FLG_table2pandas == False ):
#                                 tmpDict = {}; # Prep a dict to hold them
#                                 tmpDict_decode = [False for _ in range(0, len(tmpNamez))]; # Prep
#                                 for j in range(0, len(tmpNamez)):
#                                     if( ('bytes' == tmpDtypes[j].name[:5]) or ('object' == tmpDtypes[j].name) ): # Bytes means string to decode, hold it in a list
#                                         tmpDict[tmpNamez[j]] = [None for _ in range(0, len(arr2arr))]; # Hold it in a list
#                                         tmpDict_decode[j] = True; # Set to decode
#                                     else:
#                                         tmpDict[tmpNamez[j]] = empty( len(arr2arr), dtype=tmpDtypes[j]); # Otherwise it's a datatype to directly save
#                                     # END IF
#                                 # END FOR j
                                
#                                 for j in range(0, len(tmpNamez)):
#                                     if( tmpDict_decode[j] ):
#                                         for i in range(0, len(arr2arr)):
#                                             tmpDict[tmpNamez[j]][i] = arr2arr[i][j].decode('utf-8'); # Unpack into the dict after decoding the bytes to string
#                                         # END FOR i
#                                     else:
#                                         for i in range(0, len(arr2arr)):
#                                             tmpDict[tmpNamez[j]][i] = arr2arr[i][j]; # Unpack into the dict directly
#                                         # END FOR i
#                                     # END IF
#                                 # END FOR j
#                                 diction[key].append(tmpDict); #get that dataset out after some work
#                             else:
#                                 # Prep a Pandas dataframe to hold them
#                                 tmpDF = DataFrame(arr2arr);
#                                 for keyz in tmpDF.keys():
#                                     if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
#                                         tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
#                                     # END IF
#                                 # END FOR keyz
#                                 diction[key].append(tmpDF); #get that dataset out after some work
#                             # END IF
#                         # END IF
#                     else:
#                         diction[key].append(item[()]); #tack that dataset on
#                     # END IF
#                 else:
#                     # item.attrs.item()
#                     #otherwise it's a new data type to add in
#                     if( isinstance(item[()], ndarray) and isinstance(item[()][0], void) ): # Rely on lazy python if statement comprehension
#                         # It's likely a table!
#                         arr2arr = item[()].copy(); # Copying it out greatly speeds up work - apparently it reads off disk otherwise (woof)
#                         tmpNamez = arr2arr.dtype.names; # Get the names of the table headers
#                         tmpDtypes = arr2arr.dtype; # Get the names of the table dtypes
#                         if( len(tmpNamez) == 1 ):
#                             # It's just text most likely
#                             if( FLG_table2pandas == False ):
#                                 tmpStrang = [None for _ in range(0, len(arr2arr))]; # Preallocate
#                                 for i in range(0, len(arr2arr)):
#                                     if( ('bytes' == tmpDtypes[0].name[:5]) or ('object' == tmpDtypes[0].name) ): # Bytes means string to decode, hold it in a list
#                                         tmpStrang[i] = arr2arr[i][0].decode('utf-8').strip(' '); # Decode the string
#                                     else:
#                                         tmpStrang[i] = arr2arr[i][0]; # Get the thing
#                                     # END IF
#                                 # END FOR i
#                                 # diction[key] = ['\n'.join(tmpStrang)]; #get that string out after some work [combines, decided safer not to even if it's most likely text]
#                                 diction[key] = [tmpStrang]; #get that string out after some work
#                             else:
#                                 # Prep a Pandas dataframe to hold them
#                                 tmpDF = DataFrame(arr2arr);
#                                 for keyz in tmpDF.keys():
#                                     if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
#                                         tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
#                                     # END IF
#                                 # END FOR keyz
#                                 diction[key] = [tmpDF]; #get that dataset out after some work
#                             # END IF
#                         else:
#                             # It's a table array
#                             if( FLG_table2pandas == False ):
#                                 tmpDict = {}; # Prep a dict to hold them
#                                 tmpDict_decode = [False for _ in range(0, len(tmpNamez))]; # Prep
#                                 for j in range(0, len(tmpNamez)):
#                                     if( ('bytes' == tmpDtypes[j].name[:5]) or ('object' == tmpDtypes[j].name) ): # Bytes means string to decode, hold it in a list
#                                         tmpDict[tmpNamez[j]] = [None for _ in range(0, len(arr2arr))]; # Hold it in a list
#                                         tmpDict_decode[j] = True; # Set to decode
#                                     else:
#                                         tmpDict[tmpNamez[j]] = empty( len(arr2arr), dtype=tmpDtypes[j]); # Otherwise it's a datatype to directly save
#                                     # END IF
#                                 # END FOR j
                                
#                                 for j in range(0, len(tmpNamez)):
#                                     if( tmpDict_decode[j] ):
#                                         for i in range(0, len(arr2arr)):
#                                             tmpDict[tmpNamez[j]][i] = arr2arr[i][j].decode('utf-8'); # Unpack into the dict after decoding the bytes to string
#                                         # END FOR i
#                                     else:
#                                         for i in range(0, len(arr2arr)):
#                                             tmpDict[tmpNamez[j]][i] = arr2arr[i][j]; # Unpack into the dict directly
#                                         # END FOR i
#                                     # END IF
#                                 # END FOR j
#                                 diction[key] = [tmpDict]; #get that dataset out after some work
#                             else:
#                                 # Prep a Pandas dataframe to hold them
#                                 tmpDF = DataFrame(arr2arr);
#                                 for keyz in tmpDF.keys():
#                                     if( ('bytes' == tmpDF[keyz].dtype.name[:5]) or ('object' == tmpDF[keyz].dtype.name) ): # Bytes means string to decode, hold it in a list
#                                         tmpDF[keyz] = tmpDF[keyz].str.decode('utf-8'); # Convert byte to string
#                                     # END IF
#                                 # END FOR keyz
#                                 diction[key] = [tmpDF]; #get that dataset out after some work
#                             # END IF
#                         # END IF
#                     else:
#                         diction[key] = [item[()]]; #get that dataset out
#                     # END IF
#                 #END IF
#                 for keyA, itemA in item.attrs.items(): # Check attributes on the datset itself
#                     if( keyA == key+'!!internal!!_type' ): # An internal key on the dataset for Python built-ins
#                         if( itemA == "<class 'list'>" ): # List
#                             if( all(isinstance(val, bytes) for val in diction[key][-1]) ): # If it's all bytes, it's a list of strings
#                                 diction[key][-1] = [val.decode('UTF-8') for val in diction[key][-1] ]; # Convert to list while making each byte object into a string
#                             else:
#                                 diction[key][-1] = list(diction[key][-1]); # Convert to list
#                             # END IF
#                         elif( itemA == "<class 'tuple'>" ):
#                             if( all(isinstance(val, bytes) for val in diction[key][-1]) ): # If it's all bytes, it's a list of strings
#                                 diction[key][-1] = (val.decode('UTF-8') for val in diction[key][-1] ); # Convert to list while making each byte object into a string
#                             else:
#                                 diction[key][-1] = tuple(diction[key][-1]); # Convert to list
#                             # END IF
#                         # END IF
#                     elif( keyA == key+'!!internal!!_dtype' ): # An internal key on the dataset for special Numpy stuff (time)
#                         diction[key][-1] = diction[key][-1].astype(itemA); # Convert as requested
#                     # END IF
#                 # END FOR keyA, itemA
#                 if( FLG_oneRead ):
#                     diction[key] = diction[key][-1]; # Unlock the lock
#                 # END IF
#             elif( isinstance(item, h5py.Group) ): #group catch (recurisve time)
#                 if( key not in diction.keys() ):
#                     diction[key] = {}; #create a sub-dict if not already there
#                 #END IF
#                 if( any(['pandas_type' in tup[0] for tup in list(item.attrs.items())]) and any(['pandas_version' in tup[0] for tup in list(item.attrs.items())]) ):
#                     # This is actually a pandas, and we must eject to read it
#                     diction['!!pandas!internal!!_df'] = None;
#                     break
#                 else:
#                     diction[key] = subfun_h5_reader(item, h5_path_prefix=None, diction=diction[key], FLG_oneRead=FLG_oneRead, FLG_singleListRemove=FLG_singleListRemove, FLG_table2pandas=FLG_table2pandas); #recurse
#                 # END IF
#             #END IF
#         #END FOR key, item
#     # END IF
    
#     if( FLG_singleListRemove ):
#         # Roll through the diction keys and look for lists of length 1 and remove the list
#         for keyz in diction.keys():
#             if( isinstance(diction[keyz], (list, tuple)) ):
#                 diction[keyz] = diction[keyz][0]; # Get it out
#             # END IF
#         # END FOR keyz
#     # END IF
    
#     return diction
# #END DEF