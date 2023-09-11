# Takes a dict and writes it to an h5py file, recursive powers are here as well
# h5_path_prefix is only on the top level, ignored for recursive dict checks
# keyz_ignore is only on the top level, ignored for recusive dict checks

from numpy import isscalar

def subfun_h5_writer(h5_file, diction, h5_path_prefix=None, keyz_ignore=None):
    
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
        h5_file = h5_file[h5_path_prefix]; #go to the prefix req'd
    #END IF
    
    
    #--- load everything into the file ---
    for keyz_per in keyz:
        if( isscalar(diction[keyz_per]) == False ):
            if( type(diction[keyz_per]) is dict ):
                h5_file_subgrp = h5_file.create_group(keyz_per); #create a handle for a subgroup
                subfun_h5_writer(h5_file_subgrp, diction[keyz_per], h5_path_prefix=None, keyz_ignore=None); #recurse!
            else:
                h5_file.create_dataset(keyz_per, data=diction[keyz_per], chunks=diction[keyz_per].shape, compression="gzip", compression_opts=9, shuffle=True, fletcher32=True); #write that data
            #END IF
        else:
            #if size 1, add it as an attribute
            h5_file.attrs[keyz_per] = diction[keyz_per]; #save the attribute
        #END IF
    #END FOR keyz_per
#END DEF