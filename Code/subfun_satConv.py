'''
Converts a satillite type array to its single character format that follows RINEX formatting
RINEX formating | G: GPS, R: GLONASS, S: SBAS Payload, E: Galileo, C: BeiDou, J: QZSS, I: NavIC
satConv_toStr_fromLongStr converts from the name to single character
satConv_toStr_fromInt converts from a number designation to single character
satConv_toInt_fromStr converts from single character to the number designation (for cases where calcs need to be fast)
'''

import numpy as np

def satConv_toStr_fromLongStr(satType_longStr,FLG_unicode=False):
    # You probably need satType_staticLen if FLG_unicode is False b/c then it's using numpy's static string method which has a set length (e.g. S8 == 8 long always)
    if( type(satType_longStr) == list ):
        import sys
        print('ERROR in satConv_toStr_fromLongStr: list satType_longStr not supported go write some support');
        sys.crash() 
    else:
        if( FLG_unicode == False ):
            satType_staticLen = satType_longStr.dtype.itemsize; #get the static length
            satType_str = np.empty(satType_longStr.size,dtype='S1'); #preallocate
            satType_str[satType_longStr == np.string_('GPS'.ljust(satType_staticLen))] = np.string_('G'); #record GPS
            satType_str[satType_longStr == np.string_('GLONASS'.ljust(satType_staticLen))] = np.string_('R'); #record GLONASS
            satType_str[satType_longStr == np.string_('SBAS'.ljust(satType_staticLen))] = np.string_('S'); #record SBAS Payload
            satType_str[satType_longStr == np.string_('GALILEO'.ljust(satType_staticLen))] = np.string_('E'); #record Galileo
            satType_str[satType_longStr == np.string_('QZSS'.ljust(satType_staticLen))] = np.string_('J'); #record QZSS
            satType_str[satType_longStr == np.string_('BEIDOU'.ljust(satType_staticLen))] = np.string_('C'); #record BeiDou
            satType_str[satType_longStr == np.string_('NAVIC'.ljust(satType_staticLen))] = np.string_('I'); #record NavIC
        else:
            import sys
            print('ERROR in satConv_toStr_fromLongStr: Unicode satType_longStr not supported go write some support');
            sys.crash() 
        #END IF
    #END IF
    return satType_str
#END DEF

def satConv_toStr_fromInt(satType_int):
    satType_str = np.empty(satType_int.size,dtype='S1'); #preallocate
    satType_str[satType_int == 0] = np.string_('G'); #record GPS
    satType_str[satType_int == 1] = np.string_('R'); #record GLONASS
    satType_str[satType_int == 2] = np.string_('S'); #record SBAS Payload
    satType_str[satType_int == 3] = np.string_('E'); #record Galileo
    satType_str[satType_int == 4] = np.string_('J'); #record QZSS
    satType_str[satType_int == 5] = np.string_('C'); #record BeiDou
    satType_str[satType_int == 6] = np.string_('I'); #record NavIC
    return satType_str
#END DEF

def satConv_toInt_fromStr(satType_str):
    satType_int = np.empty(satType_str.size,dtype=np.int16); #preallocate
    satType_int[satType_str == np.string_('G')] = 0; #record GPS
    satType_int[satType_str == np.string_('R')] = 1; #record GLONASS
    satType_int[satType_str == np.string_('S')] = 2; #record SBAS Payload
    satType_int[satType_str == np.string_('E')] = 3; #record Galileo
    satType_int[satType_str == np.string_('J')] = 4; #record QZSS
    satType_int[satType_str == np.string_('C')] = 5; #record BeiDou
    satType_int[satType_str == np.string_('I')] = 6; #record NavIC
    return satType_int
#END DEF