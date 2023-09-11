#convert to mag
import numpy as np
import aacgmv2 #install with: pip install aacgmv2 [need buildtools what a pain]
import sys, os
def convert_to_mag(lat, long, alt4mag, time4mag, FLG_fixFails=True):
    stdout_old = sys.stdout; #copy current stdout, inspired by https://stackoverflow.com/a/8447352/2403531 [DOES NOT WORK UGH]
    sys.stdout = open(os.devnull, 'w'); #new stdout pit for the "runtime errors"
    #check lat and long (only lat, assume long matches)
    if( np.ndim(lat) == 0 ): #only check lat, all gas no brakes
        lat = np.atleast_1d(lat); #superior method for ensuring a 1d array
        long = np.atleast_1d(long);
    elif( isinstance(lat,list) | isinstance(lat,tuple) ): #only check lat, all gas no brakes
        lat = np.asarray(lat);
        long = np.asarray(long);
    #END IF
    #check alt4mag
    if( np.ndim(alt4mag) == 0 ):
        alt4mag = np.atleast_1d(alt4mag);
    elif( isinstance(alt4mag,list) | isinstance(alt4mag,tuple) ): #only check lat, all gas no brakes
        alt4mag = np.asarray(alt4mag);
    #END IF
    if( (alt4mag.size == 1) & (lat.size > 1) ):
        alt4mag = np.ones(lat.size)*alt4mag[0]; #make same size
    #END IF
    if( (not isinstance(time4mag,list)) & (not isinstance(time4mag,tuple)) ): #needs to be list since weird object
        if( isinstance(time4mag,np.ndarray) ):
            if( (time4mag.size == 1) & (lat.size > 1) ):
                time4mag = [time4mag[0] for _ in range(lat.size)]; #make same size as lat
            else:
                time4mag = time4mag.tolist(); #somehow?
            #END IF
        else:
            if( lat.size == 1 ):
                time4mag = [time4mag];
            else:
                time4mag = [time4mag for _ in range(lat.size)]; #make same size as lat
            #END IF
        #END IF
    #END IF
    lat_mag = np.empty( (lat.size) )*np.nan; #preallocate
    long_mag = np.empty( (lat.size) )*np.nan; #preallocate
    alt_mag = np.empty( (lat.size) )*np.nan; #preallocate
    if( FLG_fixFails == False ):
        for i in range(0,lat.size):
            [lat_mag[i], long_mag[i], alt_mag[i]] = aacgmv2.convert_latlon(lat[i], long[i], alt4mag[i], time4mag[i], method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
        #END FOR i
    else:
        for i in range(0,lat.size):
            inceasor = 0; #km, altitude increasor
            while( np.isnan(alt_mag[i]) ):
                [lat_mag[i], long_mag[i], alt_mag[i]] = aacgmv2.convert_latlon(lat[i], long[i], alt4mag[i]+inceasor, time4mag[i], method_code='G2A'); #converts from geographic to geomagnetic (AACGMv2)
                inceasor += 100; #increase every failure
            #END WHILE
        #END FOR i
    #END IF

    sys.stdout = stdout_old; #restore prev stdout
    return lat_mag, long_mag
#END DEF

def convert_to_mlt(long_mag, time4mag):
    #check long_mag
    if( np.ndim(long_mag) == 0 ): #only check lat, all gas no brakes
        long_mag = np.atleast_1d(long_mag);
    elif( isinstance(long_mag,list) | isinstance(long_mag,tuple) ): #only check lat, all gas no brakes
        long_mag = np.asarray(long_mag);
    #END IF
    if( (not isinstance(time4mag,list)) & (not isinstance(time4mag,tuple)) ): #needs to be list since weird object
        if( isinstance(time4mag,np.ndarray) ):
            if( (time4mag.size == 1) & (long_mag.size > 1) ):
                time4mag = [time4mag[0] for _ in range(long_mag.size)]; #make same size as lat
            else:
                time4mag = time4mag.tolist(); #somehow?
            #END IF
        else:
            if( long_mag.size == 1 ):
                time4mag = [time4mag];
            else:
                time4mag = [time4mag for _ in range(long_mag.size)]; #make same size as lat
            #END IF
        #END IF
    #END IF
    mlt = aacgmv2.convert_mlt(long_mag+180, time4mag, m2a=False); #convert from maglong to MLT in hrs
    
    return mlt
#END DEF