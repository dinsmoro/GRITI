#returns the stuff at whatever subdict you want
def dictDelver_get(dictDelver, dictDelver_subdicts, FLG_crashOnFail=False):
    #inspired by https://stackoverflow.com/a/49290758
    dictDelver_alias = dictDelver; #make an alias in here
    for dic in range(0,len(dictDelver_subdicts)): #delve through all
        if( dictDelver_subdicts[dic] in dictDelver_alias ):
            dictDelver_alias = dictDelver_alias[dictDelver_subdicts[dic]]
        else:
            if( FLG_crashOnFail == False ):
                print('WARNING in dictDelver: Dict Key '+str(dictDelver_subdicts[dic])+' not in the supplied dict.');
                return dictDelver
            else:
                print('ERROR in dictDelver: Dict Key '+str(dictDelver_subdicts[dic])+' not in the supplied dict.');
                import sys
                sys.crash(); #byee
            #END IF
        #END IF
    #END FOR dic
    
    return dictDelver_alias #return subdict you want
#END DEF

#inserts something into a dict at whatever subdict you want (and returns the new dict)
def dictDelver_insert(dictDelver, dictDelver_subdicts, dictDelver_insert, createMissing=False, FLG_crashOnFail=False):
    #inspired by https://stackoverflow.com/a/49290758
    dictDelver_alias = dictDelver; #make an alias in here
    for dic in range(0,len(dictDelver_subdicts)-1): #delve through all but the last key (which is set)
        if( dictDelver_subdicts[dic] in dictDelver_alias ):
            dictDelver_alias = dictDelver_alias[dictDelver_subdicts[dic]]
        elif( createMissing ):
            dictDelver_alias = dictDelver_alias.setdefault(dictDelver_subdicts[dic], {}); #create the key as a dictionary
        else:
            if( FLG_crashOnFail == False ):
                print('WARNING in dictDelver_insert: Dict Key '+str(dictDelver_subdicts[dic])+' not in the supplied dict and createMissing=False. Value was NOT inserted.');
                return dictDelver
            else:
                print('ERROR in dictDelver_insert: Dict Key '+str(dictDelver_subdicts[dic])+' not in the supplied dict and createMissing=False. Value was NOT inserted.');
                import sys
                sys.crash(); #byee
            #END IF
        #END IF
    #END FOR dic
    if( (dictDelver_subdicts[-1] in dictDelver_alias) or createMissing ):
        dictDelver_alias[dictDelver_subdicts[-1]] = dictDelver_insert;
    elif( (dictDelver_subdicts[-1] not in dictDelver_alias) and (not createMissing) ):
        if( FLG_crashOnFail == False ):
            print('WARNING in dictDelver_insert: Final Dict Key '+str(dictDelver_subdicts[-1])+' not in the supplied dict and createMissing=False. Value was NOT inserted.');
        else:
            print('ERROR in dictDelver_insert: Final Dict Key '+str(dictDelver_subdicts[-1])+' not in the supplied dict and createMissing=False. Value was NOT inserted.');
            import sys
            sys.crash(); #byee
    #END IF
    
    return dictDelver #return the original that's been updated
#END DEF


from numpy import isscalar, isclose, nan
def dictDelver_saver(dict2save, dict2save_subDict=None, diction=None, FLG_saveTopLevelNonDicts=True):
    #--- prep work ---
    if( diction == None ):
        diction = {}; #prepare dictionary
    #END IF
    
    #--- save top-level stuff that aren't dicts ---
    if( FLG_saveTopLevelNonDicts ):
        for key in dict2save.keys():
            if( not isinstance(dict2save[key], dict) ): #ignore dicts
                if( key in diction.keys() ):
                    if( isscalar(dict2save[key]) == False ):
                        diction[key].append(dict2save[key]); #tack that dataset on
                    else:
                        if( isclose(diction[key], dict2save[key]) == False ):
                            #only worry is if the attribute isn't consistent
                            print('-----Warning in dictDelver_saver-----');
                            print('Attribute '+key+' isn\'t the same as the previously recorded value from another dict of '+ \
                                str(diction[key])+' and this dict2save\'s value of '+str(dict2save[key])+ \
                                '.\n NaN\'ing it and try to sort that out.');
                            diction[key] = nan; #nan that attribute, figure it out later
                        #END IF
                    #END IF
                else:
                    #otherwise it's a new data type to add in
                    if( isscalar(dict2save[key]) == False ):
                        diction[key] = [dict2save[key]]; #get that dataset out
                    else:
                        diction[key] = dict2save[key]; #get that scalar out
                    #END IF
                #END IF
            #END IF
        #END FOR key
    #END IF
    
    #--- delve to the subdict to save and then save it to the diction as if it was top-level
    if( dict2save_subDict != None ):
        dict2save = dict2save[dict2save_subDict]; #initial delve if needed
    #END IF
    for key in dict2save.keys():
        if( not isinstance(dict2save[key], dict) ): #ignore dicts
            if( key in diction.keys() ):
                if( isscalar(dict2save[key]) == False ):
                    diction[key].append(dict2save[key]); #tack that dataset on
                else:
                    if( isclose(diction[key], dict2save[key]) == False ):
                        #only worry is if the attribute isn't consistent
                        print('-----Warning in dictDelver_saver-----');
                        print('Attribute '+key+' isn\'t the same as the previously recorded value from another dict of '+ \
                            str(diction[key])+' and this dict2save\'s value of '+str(dict2save[key])+ \
                            '.\n NaN\'ing it and try to sort that out.');
                        diction[key] = nan; #nan that attribute, figure it out later
                    #END IF
                #END IF
            else:
                #otherwise it's a new data type to add in
                if( isscalar(dict2save[key]) == False ):
                    diction[key] = [dict2save[key]]; #get that dataset out
                else:
                    diction[key] = dict2save[key]; #get that scalar out
                #END IF
            #END IF
        else:
            #delve!
            if( key not in diction.keys() ):
                diction[key] = {}; #create a sub-dict if not already there
            #END IF
            diction[key] = dictDelver_saver(dict2save[key], dict2save_subDict=None, diction=diction[key], FLG_saveTopLevelNonDicts=False); #recurse time
        #END IF
    #END FOR key
    
    return diction
#END DEF