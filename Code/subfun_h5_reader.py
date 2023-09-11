# Reads h5py files and shoves them into a dict, recursive powers are here as well

import h5py
from numpy import isclose, nan

def subfun_h5_reader(h5_file, h5_path_prefix=None, diction=None):
    #--- prep work ---
    if( h5_path_prefix != None ):
        h5_file = h5_file[h5_path_prefix]; #go to the prefix req'd
    #END IF
    if( diction == None ):
        diction = {}; #prepare dictionary
    #END IF
    
    #--- check for attributes 1st since they're not in with groups/datasets to preempt recursion ---
    for key, item in h5_file.attrs.items():
        if( key in diction.keys() ):
            if( isclose(diction[key], item) == False ):
                #only worry is if the attribute isn't consistent
                print('-----Warning in subfun_h5_reader-----');
                print('Attribute '+key+' isn\'t the same as the previously recorded value from another file of '+ \
                    str(diction[key])+' and this ('+str(h5_file)+') file\'s value of '+str(item)+ \
                    '.\n NaN\'ing it and try to sort that out.');
                diction[key] = nan; #nan that attribute, figure it out later
            #END IF
        else:
            diction[key] = item; #get that attribute out
        #END IF
    #END FOR key, item
    
    #--- check for everything else ---
    for key, item in h5_file.items():             
        if( isinstance(item, h5py.Dataset) ): #dataset catch
            if( key in diction.keys() ):
                diction[key].append(item[()]); #tack that dataset on
            else:
                #otherwise it's a new data type to add in
                diction[key] = [item[()]]; #get that dataset out
            #END IF
        elif( isinstance(item, h5py.Group) ): #group catch (recurisve time)
            if( key not in diction.keys() ):
                diction[key] = {}; #create a sub-dict if not already there
            #END IF
            diction[key] = subfun_h5_reader(item, h5_path_prefix=None, diction=diction[key]); #recurse
        #END IF
    #END FOR key, item
    
    return diction
#END DEF