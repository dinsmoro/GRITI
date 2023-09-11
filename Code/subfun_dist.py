"""
Geodesic distance calcs
"""

#========== WGS-84 using Vincenty alg (prefer geopy's karney alg, this has edge issues - gained usage b/c smol and simple for old comps| note is accurate for most cases) ==========
'''
Inspired by https://stackoverflow.com/a/20615459/2403531 & https://en.wikipedia.org/wiki/Vincenty's_formulae
THIS DOES NOT SUPPORT VECTORS
Calculates geodetic distance between two points specified by latitude/longitude using Vincenty inverse formula for ellipsoids
For stuff that Vincenty inverse formula fails at, falls back to haverside formula which assumes a sphere but will work
Note that references to an auxiliary sphere are about using an adjusted sphere to approximate Earth while making math easier (general gist, don't fite me)
INPUT: 
    lat1, long1: first point in decimal degrees
    lat2, long2: second point in decimal degrees
OUTPUT: 
    dist: distance in m between input pts
'''
import numpy as np
# from numba import jit, prange#, float32, float64, int64
#float32(float32,float32,float64,float64,float32,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64)

# @jit(nopython=True,nogil=True,parallel=True,fastmath=False)
def distWGS84_notBad(lat1, long1, lat2, long2):
    #--- This stuff is to make sure the math is ez and works good (numba hates it) ---
    if( np.isscalar(lat1) ):
        lat1 = np.array((lat1)); #this is to make indexing work
    #END IF
    if( np.isscalar(long1) ):
        long1 = np.array((long1));
    #END IF
    if( np.isscalar(lat2) ):
        lat2 = np.array((lat2));
    #END IF
    if( np.isscalar(long2) ):
        long2 = np.array((long2));
    #END IF
    if( (lat1.size == 1) & (lat2.size != 1) ):
        lat1 = np.repeat(lat1,lat2.size); #makes arrays for indexing
        long1 = np.repeat(long1,lat2.size);
    elif( (lat2.size == 1) & (lat1.size != 1) ):
        lat2 = np.repeat(lat2,lat1.size);
        long2 = np.repeat(long2,lat1.size);
    #END IF
    
    #From https:#en.wikipedia.org/wiki/World_Geodetic_System
    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    f = 1/298.257223563; #unitless, flattening factor
    f_minOne = 1-f; #precalc, reused
    b = a*f_minOne; #m, derived unit from a*(1-f) WGS-84 radius at poles (semi-minor axis of elipsoid)
    
    L = (long2 - long1)*np.pi/180; #rad, diff between longitudes
    U1 = np.arctan(f_minOne*np.tan(lat1*np.pi/180)); #rad, reduced latitude (latitude on the auxiliary sphere)
    U2 = np.arctan(f_minOne*np.tan(lat2*np.pi/180)); #rad, reduced latitude
    U1Sin = np.sin(U1); #precalc, reused
    U1Cos = np.cos(U1);
    U2Sin = np.sin(U2);
    U2Cos = np.cos(U2);
  
    lambchop = L.copy(); #rad, initial definition, diff in longitude on auxiliary sphere
    lambchop_prev = L+100; #keep it away
    iterLimit = 100;
    s = np.empty(L.size); #preallocate
    for i in range(0,L.size):
        while( (np.abs(lambchop[i] - lambchop_prev[i]) > 1e-12) & (iterLimit > 0) ):
            lambchopSin = np.sin(lambchop[i]);
            lambchopCos = np.cos(lambchop[i]);
            sigmaSin = np.sqrt((U2Cos[i] * lambchopSin) * (U2Cos[i] * lambchopSin) + (U1Cos[i] * U2Sin[i] - U1Sin[i] * U2Cos[i] * lambchopCos) * (U1Cos[i] * U2Sin[i] - U1Sin[i] * U2Cos[i] * lambchopCos));
            if( sigmaSin == 0 ):
                s[i] = 0; # co-incident points
                break; #get out of while loop
            #END IF
            sigmaCos = U1Sin[i] * U2Sin[i] + U1Cos[i] * U2Cos[i] * lambchopCos;
            sigma = np.arctan2(sigmaSin, sigmaCos);
            alphaSin = U1Cos[i] * U2Cos[i] * lambchopSin / sigmaSin;
            alphaSqCos = 1 - alphaSin * alphaSin;
            cos2SigmaM = sigmaCos - 2 * U1Sin[i] * U2Sin[i] / alphaSqCos;
            if( np.isnan(cos2SigmaM) ):
                cos2SigmaM = 0; # equatorial line: alphaSqCos=0 (§6)
            #END IF
            C = f / 16 * alphaSqCos * (4 + f * (4 - 3 * alphaSqCos));
            lambchop_prev[i] = lambchop[i];
            lambchop[i] = L[i] + (1 - C) * f * alphaSin * (sigma + C * sigmaSin * (cos2SigmaM + C * sigmaCos * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
            iterLimit += -1; #increment negatively
        #END WHILE
        if( iterLimit == 0 ):
            #if failed to converge, fall back to modified haversine formula yolo
            # return np.nan # formula failed to converge
            a = np.sin((lat2[i]-lat1[i])*np.pi/180 / 2)**2 + np.cos(lat1[i]*np.pi/180) * np.cos(lat2[i]*np.pi/180) * np.sin((long2[i]-long1[i])*np.pi/180 / 2)**2; #modified haversine formula
            s[i] = 6371 * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)); #modified haversine formula
        else:
            uSq = alphaSqCos * (a * a - b * b) / (b * b);
            A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
            B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
            sigmaDelta = B * sigmaSin * (cos2SigmaM + B / 4 * (sigmaCos * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sigmaSin * sigmaSin) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
            s[i] = b * A * (sigma - sigmaDelta);
        #END IF
    #END FOR i
  
    return s
#END DEF

def distWGS84_Vincenty_single(lat1, long1, lat2, long2):
    #From https:#en.wikipedia.org/wiki/World_Geodetic_System
    a = 6378137; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    f = 1/298.257223563; #unitless, flattening factor
    f_minOne = 1-f; #precalc, reused
    b = a*f_minOne; #m, derived unit from a*(1-f) WGS-84 radius at poles (semi-minor axis of elipsoid)
    
    L = (long2 - long1)*np.pi/180; #rad, diff between longitudes
    U1 = np.arctan(f_minOne*np.tan(lat1*np.pi/180)); #rad, reduced latitude (latitude on the auxiliary sphere)
    U2 = np.arctan(f_minOne*np.tan(lat2*np.pi/180)); #rad, reduced latitude
    U1Sin = np.sin(U1); #precalc, reused
    U1Cos = np.cos(U1);
    U2Sin = np.sin(U2);
    U2Cos = np.cos(U2);
  
    lambchop = L; #rad, initial definition, diff in longitude on auxiliary sphere
    lambchop_prev = 100;
    iterLimit = 100;
    while( (np.abs(lambchop - lambchop_prev) > 1e-12) & (iterLimit > 0) ):
        lambchopSin = np.sin(lambchop);
        lambchopCos = np.cos(lambchop);
        sigmaSin = np.sqrt((U2Cos * lambchopSin) * (U2Cos * lambchopSin) + (U1Cos * U2Sin - U1Sin * U2Cos * lambchopCos) * (U1Cos * U2Sin - U1Sin * U2Cos * lambchopCos));
        if( sigmaSin == 0 ):
            return 0; # co-incident points
        #END IF
        sigmaCos = U1Sin * U2Sin + U1Cos * U2Cos * lambchopCos;
        sigma = np.arctan2(sigmaSin, sigmaCos);
        alphaSin = U1Cos * U2Cos * lambchopSin / sigmaSin;
        alphaSqCos = 1 - alphaSin * alphaSin;
        cos2SigmaM = sigmaCos - 2 * U1Sin * U2Sin / alphaSqCos;
        if( np.isnan(cos2SigmaM) ):
            cos2SigmaM = 0; # equatorial line: alphaSqCos=0 (§6)
        #END IF
        C = f / 16 * alphaSqCos * (4 + f * (4 - 3 * alphaSqCos));
        lambchop_prev = lambchop;
        lambchop = L + (1 - C) * f * alphaSin * (sigma + C * sigmaSin * (cos2SigmaM + C * sigmaCos * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
        iterLimit += -1; #increment negatively
    #END WHILE
  
    if( iterLimit == 0 ):
        #if failed to converge, fall back to modified haversine formula
        # return np.nan # formula failed to converge
        a = np.sin((lat2-lat1)*np.pi/180 / 2)**2 + np.cos(lat1*np.pi/180) * np.cos(lat2*np.pi/180) * np.sin((long2-long1)*np.pi/180 / 2)**2; #modified haversine formula
        s = 6371 * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)); #modified haversine formula
    else:
        uSq = alphaSqCos * (a * a - b * b) / (b * b);
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
        sigmaDelta = B * sigmaSin * (cos2SigmaM + B / 4 * (sigmaCos * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sigmaSin * sigmaSin) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        s = b * A * (sigma - sigmaDelta);
    #END IF
  
    return s
#END DEF

def dist_haversine(lat1, long1, lat2, long2):
    dist_sphere = np.sin((lat2-lat1)*np.pi/180 / 2)**2 + np.cos(lat1*np.pi/180) * np.cos(lat2*np.pi/180) * np.sin((long2-long1)*np.pi/180 / 2)**2; #modified haversine formula
    dist_sphere = 6371 * 2 * np.arctan2(np.sqrt(dist_sphere), np.sqrt(1 - dist_sphere)); #modified haversine formula
    return dist_sphere
#END DEF