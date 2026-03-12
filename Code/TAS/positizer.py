#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geodesic distance calcs

dist = positizer_dist( loc_set1, loc_set2 )
INPUT: 
    (lat1, long1): first location in cdegrees
    (lat2, long2): second location in cdegrees
    ++ must be a numpy array of shape (n x 2) or a numpy array/tuple/list of length 2, list-of-list is OK ( lat1, lat2, lat3), (long1, long2, long3))
OUTPUT: 
    dist: distance in m between input pts
    
loc_set_dest, az_dest = positizer_dest( loc_set_origin, dist_from_orgin, az_origin )
INPUT: 
    (lat_origin, long_origin): originating location in cdegrees
    ++ must be a numpy array of shape (n x 2) or a numpy array/tuple/list of length 2, list-of-list is OK ( lat1, lat2, lat3), (long1, long2, long3))
    dist_from_orgin: distance from origin in m
    az_origin: originating azimuth in degrees
    ++ one or two inputs can be single values and are expanded to the other vector input's size 
    ++    (e.g., one loc_set_origin, many dist_from_origin, one az_origin)
OUTPUT: 
    loc_set_dest: (lat_dest, long_dest, destination location in decimal cdegrees
    az_dest: destination azimuth in degrees
"""
import numpy as np
# from scipy.optimize import minimize # Not a hard req
# from scipy.optimize import least_squares # Not a hard req
from scipy.optimize import leastsq # Not a hard req
from numba import jit # Not a hard req

# Get distance between two sets of lat/long points
def positizer_dist( loc_set1, loc_set2, solverOverride=None, FLG_inputRad=False ):
    if( solverOverride == None ):
        # First stage, try for best alg (Karney https://link.springer.com/article/10.1007/s00190-012-0578-z)
        disty = karney();
        if( disty == None ):
            # Second stage, try for second alg (Vincenty https://link.springer.com/article/10.1007/s00190-002-0263-8)
            disty = dist_director( loc_set1, loc_set2, dist_vincenty, FLG_inputRad=FLG_inputRad );
        # END IF
    else:
        solverOverride = solverOverride.replace(' ','').replace('-','').replace('_','').lower(); # Normalize
        if( solverOverride == 'karney' ):
            disty = karney();
        elif( solverOverride == 'vincenty' ):
            disty = dist_director( loc_set1, loc_set2, dist_vincenty, FLG_inputRad=FLG_inputRad );
        elif( solverOverride == 'haversine' ):
            disty = dist_director( loc_set1, loc_set2, dist_haversine, FLG_inputRad=FLG_inputRad );
        elif( solverOverride == 'greatcircle' ): # You should use haversine unless you really need the extra cycles
            disty = dist_director( loc_set1, loc_set2, dist_greatcircle, FLG_inputRad=FLG_inputRad );
        else:
            print('WARNING in positizer.positizer_dist: `solverOverride` provided that does NOT match available `karney`, `vincenty`, `haversine. Printing `solverOverride`: `'+str(solverOverride)+'`.');
            # First stage, try for best alg (Karney https://link.springer.com/article/10.1007/s00190-012-0578-z)
            disty = karney();
            if( disty == None ):
                # Second stage, try for second alg (Vincenty https://link.springer.com/article/10.1007/s00190-002-0263-8)
                disty = dist_director( loc_set1, loc_set2, dist_vincenty, FLG_inputRad=FLG_inputRad );
            # END IF
        # END IF
    # END IF
    
    return disty
# END DEF

# Get destination (loc_set_dest)
def positizer_dest( loc_set_origin, dist, az1, FLG_inputRad=False, FLG_returnRad=False ):
    # First stage, try for best alg (Karney https://link.springer.com/article/10.1007/s00190-012-0578-z)
    desty = karney(); # https://github.com/pbrod/karney/tree/main
    if( karney != None ):
        # Second stage, try for second alg (Vincenty https://link.springer.com/article/10.1007/s00190-002-0263-8)
        desty, az2 = dest_director( loc_set_origin, dist, az1, dest_vincenty, FLG_inputRad=FLG_inputRad, FLG_returnRad=FLG_returnRad );
    # END IF
    
    return desty, az2
# END DEF

def normalizor( var2norm, dimGoal, FLG_protectInputs=True):
    if( dimGoal == 2 ):
        # Common function for normalizing location arrays to a (n x 2) array with column 0 lat and column 1 long
        if( isinstance(var2norm, tuple) or isinstance(var2norm, list) ):
            if( isinstance(var2norm[0], tuple) or isinstance(var2norm[0], list) or isinstance(var2norm[0], np.ndarray) ):
                # Multi-dimension where I assume it'll be ( (lat list), (long list) )
                var2norm = np.asarray(var2norm).T; # Transpose fixes lats and longs being mixed
            else:
                # Just a simple set of (lat, long)
                # var2norm = np.asarray(var2norm).reshape(1,2); # Reformat so matrix math applies
                var2norm = np.expand_dims( np.asarray(var2norm), axis=0); # Reformat so matrix math applies
            # END IF
        elif( isinstance(var2norm, np.ndarray) ):
            if( var2norm.ndim == 2 ):
                if( (var2norm.shape[0] == 2) and (var2norm.shape[1] != 2) ): # Can't heuristic its way out of a 2x2 matrix
                    var2norm = var2norm.T; # Flip it so it's the correct way for what the matrix math is expecting
                elif( (var2norm.shape[0] == 3) and (var2norm.shape[1] != 3) ): # Can't heuristic its way out of a 3x3 matrix
                    var2norm = var2norm.T; # Flip it so it's the correct way for what the matrix math is expecting
                # END IF
            elif( var2norm.ndim == 1 ):
                # var2norm = var2norm.reshape(1,2); # Reformat so matrix math applies
                var2norm = np.expand_dims(var2norm, axis=0); # Reformat so matrix math applies
            else:
                raise Exception('ERROR in positizer.normalizor: requested 2d array is a numpy array with too many dimensions ('+str(var2norm.ndim)+').')
            # END IF
        else:
            raise Exception('ERROR in positizer.normalizor: requested 2d array is unknown type of `'+str(type(var2norm))+'`.')
        # END IF
    elif( dimGoal == 1 ):
        if( isinstance(var2norm, tuple) or isinstance(var2norm, list) ):
            var2norm = np.atleast_1d(var2norm); # Make a 1d array
        elif( isinstance(var2norm, np.ndarray) ):
            if( var2norm.ndim == 2 ):
                if( var2norm.shape[0] == 1 ):
                    var2norm = np.ravel(var2norm); # Flatten it
                elif( var2norm.shape[1] == 1 ):
                    var2norm = np.ravel(var2norm); # Flatten it
                else:
                     raise Exception('ERROR in positizer.normalizor: requested 1d array will not be flattened to 1d with these dimensions (one of them must be 1): `'+str(var2norm.shape)+'`.')
                # END IF
            elif( var2norm.ndim != 1 ):
                raise Exception('ERROR in positizer.normalizor: requested 1d array will not be flattened to 1d due to too many dimensions ('+str(var2norm.ndim)+') and this shape: `'+str(var2norm.shape)+'`.')
            # END IF
        elif( np.isscalar(var2norm) ):
            var2norm = np.atleast_1d(var2norm); # Make a 1d array
        # END IF
    else:
        raise Exception('ERROR in positizer.normalizor: unsupported dimGoal of `'+str(dimGoal)+'` set.')
    # END IF
    
    if( FLG_protectInputs ):
        return var2norm.copy()
    else:
        return var2norm
    # END IF
# END DEF


def karney():
    # Karney https://link.springer.com/article/10.1007/s00190-012-0578-z
    try:
        import karney
        print('ERROR: I haven\'t implemented this sorry.')
        return None
    except ModuleNotFoundError:
        # Try another alg
        return None
    # END TRYING
# END DEF

@jit(nopython=True,nogil=True,parallel=True,fastmath=False)
def dist_vincenty( loc_set1, loc_set2 ):
    # Inspired by https://stackoverflow.com/a/20615459/2403531 & https://en.wikipedia.org/wiki/Vincenty's_formulae
    # Calculates geodetic distance between two points specified by latitude/longitude using Vincenty inverse formula for ellipsoids
    # For stuff that Vincenty inverse formula fails at, falls back to haverside formula which assumes a sphere but will work
    # Note that references to an auxiliary sphere are about using an adjusted sphere to approximate Earth while making math easier
    
    # From https:#en.wikipedia.org/wiki/World_Geodetic_System
    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    f = 1/298.257223563; #unitless, flattening factor
    f_minOne = 1-f; #precalc, reused (f_minOne = b/a)
    b = a*f_minOne; #m, derived unit from a*(1-f) WGS-84 radius at poles (semi-minor axis of elipsoid)
    
    L = loc_set2[:,1] - loc_set1[:,1]; #rad, diff between longitudes
    U1 = np.arctan(f_minOne*np.tan(loc_set1[:,0])); #rad, reduced latitude (latitude on the auxiliary sphere)
    U2 = np.arctan(f_minOne*np.tan(loc_set2[:,0])); #rad, reduced latitude
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
            ahv = np.sin((loc_set2[:,0][i]-loc_set1[:,0][i])/ 2)**2 + np.cos(loc_set1[:,0][i]) * np.cos(loc_set2[:,0][i]) * np.sin((loc_set2[:,1][i]-loc_set1[:,1][i])/ 2)**2; #modified haversine formula
            s[i] = a * 2 * np.arctan2(np.sqrt(ahv), np.sqrt(1 - ahv)); #modified haversine formula
        else:
            uSq = alphaSqCos * (a * a - b * b) / (b * b);
            A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
            B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
            sigmaDelta = B * sigmaSin * (cos2SigmaM + B / 4 * (sigmaCos * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sigmaSin * sigmaSin) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
            s[i] = b * A * (sigma - sigmaDelta);
        #END IF
    #END FOR i
  
    return s
# END DEF

def dist_haversine( loc_set1, loc_set2 ): # More numerically stable than the greatcircle version
    dist_sphere = np.sin((loc_set2[:,0]-loc_set1[:,0])/ 2)**2 + np.cos(loc_set1[:,0]) * np.cos(loc_set2[:,0]) * np.sin((loc_set2[:,1]-loc_set1[:,1])/ 2)**2; #modified haversine formula
    dist_sphere = 6378137. * 2 * np.arctan2(np.sqrt(dist_sphere), np.sqrt(1 - dist_sphere)); #m, modified haversine formula
    return dist_sphere
# END DEF

def dist_greatcircle( loc_set1, loc_set2 ): # Requires 64-bit precision for small distances (~1 km)
    return 6378137. * np.arccos( np.sin(loc_set1[:,0])*np.sin(loc_set2[:,0]) + np.cos(loc_set1[:,0])*np.cos(loc_set2[:,0])*np.cos(loc_set1[:,1]-loc_set2[:,1]) ); #m, great circle
# END DEF

def dist_director( loc_set1, loc_set2, fun, FLG_inputRad=False, FLG_skipInputFiltering=False ):
    # Makes it so functions can directly use numpy vector math without worry
    loc_set1 = normalizor(loc_set1, 2); # Normalize the input for the functions
    loc_set2 = normalizor(loc_set2, 2); # Normalize the input for the functions
    
    loc_set1_size = loc_set1.shape[0]; # Size of this
    loc_set2_size = loc_set2.shape[0]; # Size of this
    
    # Replicate size 1 ones to fit the rest
    if( (loc_set1_size != 1) and (loc_set2_size == 1) ):
        loc_set2 = np.ones( (loc_set1_size, 2), dtype=loc_set2.dtype)*loc_set2[0, :]; # Resize to match
    elif( (loc_set2_size != 1) and (loc_set1_size == 1) ):
        loc_set1 = np.ones( (loc_set2_size, 2), dtype=loc_set1.dtype)*loc_set1[0, :]; # Resize to match
    # END IF
    
    # Call relevant function, convert loc_set1 and loc_set2 to radians in one go for efficiency later
    if( FLG_inputRad ):
        return fun( loc_set1, loc_set2 )
    else:
        return fun( loc_set1*np.pi/180., loc_set2*np.pi/180. )
    # END IF
# END DEF

@jit(nopython=True,nogil=True,parallel=True,fastmath=False)
def dest_vincenty( loc_set_origin, dist, az1, FLG_returnRad=False ):
    # inspired by https://en.wikipedia.org/wiki/Vincenty's_formulae#Direct_problem
    # Calculates destination lat/long given a starting lat/long, azimuth, and distance from said starting lat/long
    
    #From https:#en.wikipedia.org/wiki/World_Geodetic_System
    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    f = 1/298.257223563; #unitless, flattening factor
    f_minOne = 1-f; #precalc, reused (f_minOne = b/a)
    b = a*f_minOne; #m, derived unit from a*(1-f) WGS-84 radius at poles (semi-minor axis of elipsoid)
    
    U1 = np.arctan(f_minOne*np.tan(loc_set_origin[:, 0])); #rad, reduced latitude (latitude on the auxiliary sphere)
    U1Sin = np.sin(U1); #precalc, reused
    U1Cos = np.cos(U1);
    az1Sin = np.sin(az1);
    az1Cos = np.cos(az1);
    sigma1 = np.arctan2(np.tan(U1), az1Cos); #rad
    alphaSin = U1Cos*az1Sin;
    alphaCosSq = np.cos(np.arcsin(alphaSin))**2; #precalc, reused
    uSq = (1. - alphaSin**2)*(a**2 - b**2)/b**2;
    A = 1. + uSq / 16384. * (4096. + uSq * (-768. + uSq * (320. - 175. * uSq)));
    B = uSq / 1024. * (256. + uSq * (-128. + uSq * (74. - 47. * uSq)));
    
    # Prep for solving
    siggy_c = dist/(b*A); # Consistent sigma value to adjust
    siggy = siggy_c.copy(); # initial value for intermediate sigma
    siggy_prev = siggy+100.; #keep it away
    iterLimit = 100;
    loc_set_dest = np.empty(loc_set_origin.shape); #preallocate
    az2 = np.empty(az1.size, dtype=az1.dtype); #preallocate
    siggy2m = np.empty(siggy.size, dtype=siggy.dtype); #preallocate
    for i in range(0,dist.size):
        while( (np.abs(siggy[i] - siggy_prev[i]) > 1e-12) & (iterLimit > 0) ):
            siggy2m[i] = 2.*sigma1[i] + siggy[i];
            siggyDelta = B[i]*np.sin(siggy[i]*(np.cos(siggy2m[i])+1./4.*B[i]*(np.cos(siggy[i])*(-1. + 2.*np.cos(siggy2m[i])**2) - 1./6.*B[i]*np.cos(siggy2m[i])*(-3. + 4.*np.sin(siggy[i])**2)*(-3. + 4.*np.cos(siggy2m[i])**2)))); # Big 'un
            siggy_prev[i] = siggy[i]; # Remember
            siggy[i] = siggy_c[i] + siggyDelta; # Adjust
        # END WHILE
    # END FOR i
    siggySin = np.sin(siggy); # precalc, reused
    siggyCos = np.cos(siggy);
    U1Sin_siggySin = U1Sin*siggySin;
    U1Cos_siggyCos = U1Cos*siggyCos;
    siggy2mCos = np.cos(siggy2m);
    
    loc_set_dest[:, 0] = np.arctan2(U1Sin*siggyCos + U1Cos*siggySin*az1Cos, f_minOne*np.sqrt(alphaSin**2 + (U1Sin_siggySin - U1Cos_siggyCos*az1Cos)**2)); # rad, latitude of destination
    lambchop = np.arctan2(siggySin*az1Sin, U1Cos_siggyCos - U1Sin_siggySin*az1Cos); # rad, diff in longitude on auxiliary sphere
    C = f/16.*alphaCosSq*(4. + f*(4. - 3.*alphaCosSq));
    L = lambchop - (1 - C)*f*alphaSin*(siggy + C*siggySin*(siggy2mCos + C*siggyCos*(-1 + 2.*siggy2mCos**2))); # rad
    loc_set_dest[:, 1] = L + loc_set_origin[i, 1]; # rad, longitude
    az2 = np.arctan2(alphaSin, -U1Sin_siggySin + U1Cos_siggyCos*az1Cos); # rad, destination azimuth 
    
    if( FLG_returnRad ):
        return loc_set_dest, az2
    else:
        return loc_set_dest*180./np.pi, az2*180./np.pi
    # END IF
# END DEF

def dest_director( loc_set_origin, dist, az1, fun, FLG_inputRad=False, FLG_returnRad=False, FLG_skipInputFiltering=False ):
    # Makes it so functions can directly use numpy vector math without worry
    loc_set_origin = normalizor(loc_set_origin, 2); # Normalize the input for the functions
    dist = normalizor(dist, 1); # Normalize the input for the functions
    az1 = normalizor(az1, 1); # Normalize the input for the functions
    
    loc_set_origin_size = loc_set_origin.shape[0]; # Size of this
    dist_size = dist.size; # Size of that
    az1_size = az1.size; # Size of this too
    
    # Replicate size 1 ones to fit the rest
    if( (loc_set_origin_size != 1) and ((dist_size == 1) or (az1_size == 1)) ):
        if( dist_size == 1 ):
            dist = np.ones(loc_set_origin_size, dtype=dist.dtype)*dist[0]; # Resize to match
        # END IF
        if( az1_size == 1 ):
            az1 = np.ones(loc_set_origin_size, dtype=az1.dtype)*az1[0]; # Resize to match
        # END IF
    elif( (dist_size != 1) and ((loc_set_origin_size == 1) or (az1_size == 1)) ):
        if( loc_set_origin_size == 1 ):
            loc_set_origin = np.ones( (dist_size, 2), dtype=loc_set_origin.dtype)*loc_set_origin[0, :]; # Resize to match
        # END IF
        if( az1_size == 1 ):
            az1 = np.ones(dist_size, dtype=az1.dtype)*az1[0]; # Resize to match
        # END IF
    elif( (az1_size != 1) and ((loc_set_origin_size == 1) or (dist_size == 1)) ):
        if( loc_set_origin_size == 1 ):
            loc_set_origin = np.ones( (az1_size, 2), dtype=loc_set_origin.dtype)*loc_set_origin[0, :]; # Resize to match
        # END IF
        if( dist_size == 1 ):
            dist = np.ones(az1_size, dtype=dist.dtype)*dist[0]; # Resize to match
        # END IF
    # END IF
    
    # Call relevant function, convert loc_set1 and loc_set2 to radians in one go for efficiency later
    if( FLG_inputRad ):
        return fun( loc_set_origin, dist, az1, FLG_returnRad=FLG_returnRad )
    else:
        return fun( loc_set_origin*np.pi/180, dist, az1*np.pi/180, FLG_returnRad=FLG_returnRad )
    # END IF
# END DEF

def conv_geo2ecef(lat_rad, long_rad, alt, FLG_radIn=False, FLG_kmIn=False, FLG_kmOut=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded = None):
    # WGS-84 datum, lat/long -> ECEF XYZ
    if( FLG_skipInputFiltering == False ):
        # Normalize inputs just in case
        if( (long_rad is not None) and (alt is not None) ):
            lat_rad = normalizor(lat_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            long_rad = normalizor(long_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            alt = normalizor(alt, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (long_rad is None) and (alt is not None) ):
            # Secret overloaded lat_rad == (lat_rad, long_rad) while long_rad is None
            lat_rad = normalizor(lat_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            alt = normalizor(alt, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded lat_rad == (lat_rad, long_rad, alt) while long_rad is None and alt is None
            lat_rad = normalizor(lat_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            alt = np.ravel(lat_rad[:, 2].copy()); # Get out the alt
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        lat_rad_size = lat_rad.size; # Size of this
        long_rad_size = long_rad.size; # Size of that
        alt_size = alt.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (lat_rad_size != 1) and ((long_rad_size == 1) or (alt_size == 1)) ):
            if( long_rad_size == 1 ):
                long_rad = np.ones(lat_rad_size, dtype=long_rad.dtype)*long_rad[0]; # Resize to match
            # END IF
            if( alt_size == 1 ):
                alt = np.ones(lat_rad_size, dtype=alt.dtype)*alt[0]; # Resize to match
            # END IF
        elif( (long_rad_size != 1) and ((lat_rad_size == 1) or (alt_size == 1)) ):
            if( lat_rad_size == 1 ):
                lat_rad = np.ones(long_rad_size, dtype=lat_rad.dtype)*lat_rad[0]; # Resize to match
            # END IF
            if( alt_size == 1 ):
                alt = np.ones(long_rad_size, dtype=alt.dtype)*alt[0]; # Resize to match
            # END IF
        elif( (alt_size != 1) and ((lat_rad_size == 1) or (long_rad_size == 1)) ):
            if( lat_rad_size == 1 ):
                lat_rad = np.ones(alt_size, dtype=lat_rad.dtype)*lat_rad[0]; # Resize to match
            # END IF
            if( long_rad_size == 1 ):
                long_rad = np.ones(alt_size, dtype=long_rad.dtype)*long_rad[0]; # Resize to match
            # END IF
        # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0; # Default
    # END IF
    
    # Convert as needed        
    if( FLG_radIn == False ):
        lat_rad = lat_rad * np.pi/180; # rad, convert from deg to rad (not *= b/c if int var it fails)
        long_rad = long_rad * np.pi/180; # rad, convert from deg to rad (not *= b/c if int var it fails)
    # END IF
    if( FLG_kmIn ):
        alt *= 1000; # m, convert from km to m
    # END IF
    
    # WGS-84 datum, lat/long -> ECEF XYZ
    # Based on https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
    # From https:#en.wikipedia.org/wiki/World_Geodetic_System
    # lat/long MUST be in radians, alt MUST be in meters
    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    f = 1/298.257223563; #unitless, flattening factor (f = 1 - b/a)
    eSq = (2 - f)*f; # unitless, Eccentricity of elispod (1 - b^2/a^2)
    lat_rad_sin = np.sin(lat_rad); # Reused
    lat_rad_cos = np.cos(lat_rad); # Reused
    N = a/np.sqrt(1 - eSq*lat_rad_sin**2); # m, reused
    N_N_alt = N + alt; # m, Reused
    X = N_N_alt*lat_rad_cos*np.cos(long_rad); # m, X component for ECEF with WGS-84 datum
    Y = N_N_alt*lat_rad_cos*np.sin(long_rad); # m, Y component for ECEF with WGS-84 datum
    Z = ((1 - eSq)*N + alt)*lat_rad_sin; # m, Z component for ECEF with WGS-84 datum
    
    # Convert as needed
    if( FLG_kmOut ):
        X /= 1000;
        Y /= 1000;
        Z /= 1000;
    # END IF
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return X, Y, Z
    elif( FLG_overloaded == 1 ):
        return np.vstack( (X, Y) ).T, Z
    elif( FLG_overloaded == 2 ):
        return np.vstack( (X, Y, Z) ).T
    # END IF
    
    # For check, this is how to call pyproj
    # from pyproj import Transformer as pyproj_Transformer
    # pyprojTransformer = pyproj_Transformer.from_crs(
    #     {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
    #     {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
    #     );
    # X_check, Y_check, Z_check = pyprojTransformer.transform(long_rad, lat_rad, alt, radians = True);
# END DEF

def conv_ecef2geo(X, Y, Z, FLG_kmIn=False, FLG_kmOut=False, FLG_radOut=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=None):
    # WGS-84 datum, ECEF XYZ -> lat/long
    if( FLG_skipInputFiltering == False ):
        # Normalize inputs just in case
        if( (Y is not None) and (Z is not None) ):
            X = normalizor(X, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Y = normalizor(Y, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Z = normalizor(Z, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (Y is None) and (Z is not None) ):
            # Secret overloaded X == (X, Y) while Y is None
            X = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Y = np.ravel(X[:, 1].copy()); # Get out the Y
            X = np.ravel(X[:, 0]); # Convert the X
            Z = normalizor(Z, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded X == (X, Y, Z) while Y is None and Z is None
            X = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Z = np.ravel(X[:, 2].copy()); # Get out the Z
            Y = np.ravel(X[:, 1].copy()); # Get out the Y
            X = np.ravel(X[:, 0]); # Convert the X
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        X_size = X.size; # Size of this
        Y_size = Y.size; # Size of that
        Z_size = Z.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (X_size != 1) and ((Y_size == 1) or (Z_size == 1)) ):
            if( Y_size == 1 ):
                Y = np.ones(X_size, dtype=Y.dtype)*Y[0]; # Resize to match
            # END IF
            if( Z_size == 1 ):
                Z = np.ones(X_size, dtype=Z.dtype)*Z[0]; # Resize to match
            # END IF
        elif( (Y_size != 1) and ((X_size == 1) or (Z_size == 1)) ):
            if( X_size == 1 ):
                X = np.ones(Y_size, dtype=X.dtype)*X[0]; # Resize to match
            # END IF
            if( Z_size == 1 ):
                Z = np.ones(Y_size, dtype=Z.dtype)*Z[0]; # Resize to match
            # END IF
        elif( (Z_size != 1) and ((X_size == 1) or (Y_size == 1)) ):
            if( X_size == 1 ):
                X = np.ones(Z_size, dtype=X.dtype)*X[0]; # Resize to match
            # END IF
            if( Y_size == 1 ):
                Y = np.ones(Z_size, dtype=Y.dtype)*Y[0]; # Resize to match
            # END IF
        # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0; # Default
    # END IF
    
    # Convert as needed
    if( FLG_kmIn ):
        X *= 1000; # m, convert from km to m
        Y *= 1000; # m, convert from km to m
        Z *= 1000; # m, convert from km to m
    # END IF
    
    # WGS-84 datum, ECEF XYZ -> lat/long
    # Based on https://link.springer.com/article/10.1007/s00190-010-0419-x
    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    # X, Y, Z MUST be in meters
    f = 1/298.257223563; #unitless, flattening factor (f = 1 - b/a)
    eSq = (2 - f)*f; # unitless, Eccentricity of elispoid (1 - b^2/a^2) squared
    eQuad = eSq*eSq; # unitless,  Eccentricity of elispoid (1 - b^2/a^2) quadrophenia
    rho = np.sqrt(X**2 + Y**2);
    p = (rho/a)**2;
    q = (1-eSq)/a**2*Z**2;
    r = (p + q - eQuad)/6;
    s = eQuad*p*q;
    sigmaSqrt= np.sqrt(8*r**3 + s);
    sSqrt = np.sqrt(s);
    u = r + 1/2*((sigmaSqrt + sSqrt)**(2/3) + (sigmaSqrt - sSqrt)**(2/3));
    v = np.sqrt(u**2 + eQuad*q);
    uNv = u + v;
    w = eSq*(uNv - q)/(2*v);
    k = (uNv)/(w + np.sqrt(w**2 + uNv));
    d = k*rho/(k + eSq);
    g = np.sqrt(d**2 + Z**2);
    alt = (k + eSq - 1)/k*g; # m, altitude result
    # lat_rad = 2*np.arctan2(Z, d + g);# rad
    lat_rad = 2*np.arctan(Z/(d + g)); # rad, Looks like arctan2 is NOT needed - speed up!
    # Longitude gets special handling
    long_rad = np.empty(Z.size, dtype=Z.dtype)*np.nan;# rad, set to nan (extra work) to detect if the below statements failed
    kj = (np.sqrt(2) - 1)*np.abs(Y) < (rho + X);
    # long_rad[kj] = 2*np.arctan2(Y[kj], rho[kj] + X[kj]);# rad
    long_rad[kj] = 2*np.arctan(Y[kj]/(rho[kj] + X[kj])); # rad, Looks like arctan2 is NOT needed - speed up!
    kjr = kj.copy();
    kj = ((rho + Y) < ((np.sqrt(2) + 1)*np.abs(X))) & ~kjr;
    # long_rad[kj] = 2*np.arctan2(X[kj], rho[kj] - Y[kj]) - np.pi/2;# rad
    long_rad[kj] = 2*np.arctan(X[kj]/(rho[kj] - Y[kj])) - np.pi/2; # rad, Looks like arctan2 is NOT needed - speed up!
    kjr = kj | kjr; # Combine
    kj = ((rho - Y) < (np.sqrt(2) + 1)*np.abs(X)) & ~kjr;
    # long_rad[kj] = np.pi/2 - 2*np.arctan2(X[kj], rho[kj] + Y[kj]);# rad
    long_rad[kj] = np.pi/2 - 2*np.arctan(X[kj]/(rho[kj] + Y[kj])); # rad, Looks like arctan2 is NOT needed - speed up!
    # if np.any(np.isnan(long_rad)): raise Exception('ERROR in positizer.conv_ecef2geo: variable long_rad has NaNs in it which should NOT have happened but did, so this alg needs a looksee. Sorry.');
    kjp = long_rad > np.pi; # Catch things to flip back into -pi to pi range
    kjn = long_rad < -np.pi; # Catch things to flip back into -pi to pi range
    long_rad[kjp] -= 2*np.pi; # flip it around
    long_rad[kjn] += 2*np.pi; # flip it around
    
    # Convert as needed
    if( FLG_kmOut ):
        alt /= 1000; # km, convert from m to km
    # END IF
    if( FLG_radOut == False ):
        lat_rad *= 180/np.pi; # deg, convert from rad to deg
        long_rad *= 180/np.pi; # deg, convert from rad to deg
    # END IF
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return lat_rad, long_rad, alt
    elif( FLG_overloaded == 1 ):
        return np.vstack( (lat_rad, long_rad) ).T, alt
    elif( FLG_overloaded == 2 ):
        return np.vstack( (lat_rad, long_rad, alt) ).T
    # END IF
    return 

    # For check, this is how to call pyproj
    # from pyproj import Transformer as pyproj_Transformer
    # pyprojTransformer = pyproj_Transformer.from_crs(
    #     {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
    #     {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
    #     );
    # long_rad_check, lat_rad_check, alt_check = pyprojTransformer.transform(X, Y, Z, radians=True)
# END DEF


# ENU is a vector, so that's why you need to define the starting XYZ and the ending XYZenu
def conv_ecef2enu(X, Y, Z, Xenu, Yenu, Zenu, FLG_kmIn=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=None):
    # WGS-84 datum, ENU (east-north-up) -> ECEF XYZ
    if( FLG_skipInputFiltering == False ):
        # Normalize inputs just in case
        if( (Y is not None) and (Z is not None) ):
            X = normalizor(X, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Y = normalizor(Y, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Z = normalizor(Z, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (Y is None) and (Z is not None) ):
            # Secret overloaded X == (X, Y) while Y is None
            X = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Y = np.ravel(X[:, 1].copy()); # Get out the Y
            X = np.ravel(X[:, 0]); # Convert the X
            Z = normalizor(Z, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded X == (X, Y, Z) while Y is None and Z is None
            # raise Exception('ERROR in positizer.conv_geo2ecef: I didn\'t code positizer.normalizor for 3 dims yet, sorry.');
            X = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Z = np.ravel(X[:, 2].copy()); # Get out the Z
            Y = np.ravel(X[:, 1].copy()); # Get out the Y
            X = np.ravel(X[:, 0]); # Convert the X
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Normalize inputs just in case
        if( (Yenu is not None) and (Zenu is not None) ):
            Xenu = normalizor(Xenu, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Yenu = normalizor(Yenu, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Zenu = normalizor(Zenu, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (Yenu is None) and (Zenu is not None) ):
            # Secret overloaded Xenu == (Xenu, Yenu) while Yenu is None
            Xenu = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Yenu = np.ravel(Xenu[:, 1].copy()); # Get out the Yenu
            Xenu = np.ravel(Xenu[:, 0]); # Convert the Xenu
            Zenu = normalizor(Zenu, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded Xenu == (Xenu, Yenu, Zenu) while Yenu is None and Zenu is None
            Xenu = normalizor(Xenu, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Zenu = np.ravel(Xenu[:, 2].copy()); # Get out the Zenu
            Yenu = np.ravel(Xenu[:, 1].copy()); # Get out the Yenu
            Xenu = np.ravel(Xenu[:, 0]); # Convert the Xenu
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0; # Default
    # END IF
        
    # Convert as needed        
    if( FLG_kmIn == True ):
        X *= 1000; # m, km -> m conversion
        Y *= 1000;
        Z *= 1000;
        Xenu *= 1000; # m, km -> m conversion
        Yenu *= 1000;
        Zenu *= 1000;
    # END IF
    
    lat_rad, long_rad, _ = conv_ecef2geo(X, Y, Z, FLG_kmIn=False, FLG_kmOut=False, FLG_radOut=True, FLG_skipInputFiltering=False, FLG_overloaded=0); # Need lat_rad/long_rad so calc from X Y Z
        
    latSin = np.sin(lat_rad); # Pre-calc
    latCos = np.cos(lat_rad); # Pre-calc
    longSin = np.sin(long_rad); # Pre-calc
    longCos = np.cos(long_rad); # Pre-calc
    
    # E = -longSin*(Xenu - X) + longCos*(Yenu - Y); # Individual calculations
    # N = -latSin*longCos*(Xenu - X) + -latSin*longSin*(Yenu - Y) + latCos*(Zenu - Z);
    # U = latCos*longCos*(Xenu - X) + latCos*longSin*(Yenu - Y) + latSin*(Zenu - Z);
    (E, N, U) = np.einsum('ijk,jk->ik', np.array( ((-longSin, longCos, np.zeros(longCos.size, dtype=longCos.dtype)), (-latSin*longCos, -latSin*longSin, latCos), (latCos*longCos, latCos*longSin, latSin)) ), (np.array( (Xenu, Yenu, Zenu) ) - np.array( (X, Y, Z) )), optimize = True); # Matrix mult via einsum to deal with that extra dim
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return E, N, U
    elif( FLG_overloaded == 1 ):
        return np.vstack( (E, N) ).T, U
    elif( FLG_overloaded == 2 ):
        return np.vstack( (E, N, U) ).T
    # END IF
    return 
# END DEF

def conv_enu2ecef(lat_rad, long_rad, alt, E, N, U, FLG_radIn=False, FLG_kmIn=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=None):
    # WGS-84 datum, ENU (east-north-up) -> ECEF XYZ
    if( FLG_skipInputFiltering == False ):        
        # Normalize inputs just in case
        if( (long_rad is not None) and (alt is not None) ):
            lat_rad = normalizor(lat_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            long_rad = normalizor(long_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            alt = normalizor(alt, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (long_rad is None) and (alt is not None) ):
            # Secret overloaded lat_rad == (lat_rad, long_rad) while long_rad is None
            lat_rad = normalizor(lat_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            alt = normalizor(alt, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded lat_rad == (lat_rad, long_rad, alt) while long_rad is None and alt is None
            lat_rad = normalizor(lat_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            alt = np.ravel(alt[:, 2].copy()); # Get out the alt
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Normalize inputs just in case
        if( (N is not None) and (U is not None) ):
            E = normalizor(E, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            N = normalizor(N, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            U = normalizor(U, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (N is None) and (U is not None) ):
            # Secret overloaded E == (E, N) while N is None
            E = normalizor(E, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            N = np.ravel(E[:, 1].copy()); # Get out the N
            E = np.ravel(E[:, 0]); # Convert the E
            U = normalizor(U, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded E == (E, N, Y) while N is None and U is None
            E = normalizor(E, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            U = np.ravel(U[:, 2].copy()); # Get out the U
            N = np.ravel(E[:, 1].copy()); # Get out the N
            E = np.ravel(E[:, 0]); # Convert the E
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0; # Default
    # END IF
    
    # Convert as needed        
    if( FLG_radIn == False ):
        lat_rad = lat_rad * np.pi/180; # rad, convert from deg to rad (not *= b/c if int var it fails)
        long_rad = long_rad * np.pi/180; # rad, convert from deg to rad (not *= b/c if int var it fails)
    # END IF
    if( FLG_kmIn == True ):
        alt *= 1000; # m, km -> m conversion
        E *= 1000; # m, km -> m conversion
        N *= 1000;
        U *= 1000;
    # END IF
    
    X, Y, Z = conv_geo2ecef(lat_rad, long_rad, alt, FLG_radIn=True, FLG_kmIn=False, FLG_kmOut=False, FLG_skipInputFiltering=False, FLG_overloaded=0); # Need X Y Z so calc from lat_rad/long_rad/alt
        
    latSin = np.sin(lat_rad); # Pre-calc
    latCos = np.cos(lat_rad); # Pre-calc
    longSin = np.sin(long_rad); # Pre-calc
    longCos = np.cos(long_rad); # Pre-calc
    
    # Xenu2 = -longSin*E - latSin*longCos*N + latCos*longCos*U + X; # Individual calculations
    # Yenu2 = longCos*E - latSin*longSin*N + latCos*longSin*U + Y;
    # Zenu2 = latCos*N + latSin*U + Z;
    (Xenu, Yenu, Zenu) = np.einsum('ijk,jk->ik', np.array( ((-longSin, -latSin*longCos, latCos*longCos), (longCos, -latSin*longSin, latCos*longSin), (np.zeros(longCos.size, dtype=longCos.dtype), latCos, latSin)) ), np.array( (E, N, U) ), optimize = True) + np.array( (X, Y, Z) ); # Matrix mult via einsum to deal with that extra dim    
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return Xenu, Yenu, Zenu
    elif( FLG_overloaded == 1 ):
        return np.vstack( (Xenu, Yenu) ).T, Zenu
    elif( FLG_overloaded == 2 ):
        return np.vstack( (Xenu, Yenu, Zenu) ).T
    # END IF
    return 
# END DEF

def conv_aer2ecef(az_rad, el_rad, rangez, lat_rad, long_rad, alt, FLG_radIn_azel=False, FLG_rangekm=False, FLG_radIn_latlong=False, FLG_altkm=False, FLG_locECEF=False, FLG_kmOut=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=None):
    # Azimuth/Elevation/Range from location (lat/long/alt OR X/Y/Z in ECEF) -> ECEF X/Y/Z at point
    if( FLG_skipInputFiltering == False ):
        # Normalize inputs just in case
        if( (el_rad is not None) and (rangez is not None) ):
            az_rad = normalizor(az_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            el_rad = normalizor(el_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            rangez = normalizor(rangez, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (el_rad is None) and (rangez is not None) ):
            # Secret overloaded az_rad == (az_rad, el_rad) while el_rad is None
            az_rad = normalizor(az_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            el_rad = np.ravel(az_rad[:, 1].copy()); # Get out the el_rad
            az_rad = np.ravel(az_rad[:, 0]); # Convert the az_rad
            rangez = normalizor(rangez, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded az_rad == (az_rad, el_rad, rangez) while el_rad is None and rangez is None
            az_rad = normalizor(az_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            rangez = np.ravel(rangez[:, 2].copy()); # Get out the rangez
            el_rad = np.ravel(az_rad[:, 1].copy()); # Get out the el_rad
            az_rad = np.ravel(az_rad[:, 0]); # Convert the az_rad
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        az_rad_size = az_rad.size; # Size of this
        el_rad_size = el_rad.size; # Size of that
        rangez_size = rangez.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (az_rad_size != 1) and ((el_rad_size == 1) or (rangez_size == 1)) ):
            if( el_rad_size == 1 ):
                el_rad = np.ones(az_rad_size, dtype=el_rad.dtype)*el_rad[0]; # Resize to match
            # END IF
            if( rangez_size == 1 ):
                rangez = np.ones(az_rad_size, dtype=rangez.dtype)*rangez[0]; # Resize to match
            # END IF
        elif( (el_rad_size != 1) and ((az_rad_size == 1) or (rangez_size == 1)) ):
            if( az_rad_size == 1 ):
                az_rad = np.ones(el_rad_size, dtype=az_rad.dtype)*az_rad[0]; # Resize to match
            # END IF
            if( rangez_size == 1 ):
                rangez = np.ones(el_rad_size, dtype=rangez.dtype)*rangez[0]; # Resize to match
            # END IF
        elif( (rangez_size != 1) and ((az_rad_size == 1) or (el_rad_size == 1)) ):
            if( az_rad_size == 1 ):
                az_rad = np.ones(rangez_size, dtype=az_rad.dtype)*az_rad[0]; # Resize to match
            # END IF
            if( el_rad_size == 1 ):
                el_rad = np.ones(rangez_size, dtype=el_rad.dtype)*el_rad[0]; # Resize to match
            # END IF
        # END IF
        
        # Normalize inputs just in case
        if( (long_rad is not None) and (alt is not None) ):
            lat_rad = normalizor(lat_rad, 1); # Normalize the input for the functions
            long_rad = normalizor(long_rad, 1); # Normalize the input for the functions
            alt = normalizor(alt, 1); # Normalize the input for the functions
        elif( (long_rad is None) and (alt is not None) ):
            # Secret overloaded lat_rad == (lat_rad, long_rad) while long_rad is None
            lat_rad = normalizor(lat_rad, 2); # Normalize the input for the functions
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            alt = normalizor(alt, 1); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded lat_rad == (lat_rad, long_rad, alt) while long_rad is None and alt is None
            lat_rad = normalizor(lat_rad, 2); # Normalize the input for the functions
            alt = np.ravel(lat_rad[:, 2].copy()); # Get out the alt
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        lat_rad_size = lat_rad.size; # Size of this
        long_rad_size = long_rad.size; # Size of that
        alt_size = alt.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (lat_rad_size != 1) and ((long_rad_size == 1) or (alt_size == 1)) ):
            if( long_rad_size == 1 ):
                long_rad = np.ones(lat_rad_size, dtype=long_rad.dtype)*long_rad[0]; # Resize to match
            # END IF
            if( alt_size == 1 ):
                alt = np.ones(lat_rad_size, dtype=alt.dtype)*alt[0]; # Resize to match
            # END IF
        elif( (long_rad_size != 1) and ((lat_rad_size == 1) or (alt_size == 1)) ):
            if( lat_rad_size == 1 ):
                lat_rad = np.ones(long_rad_size, dtype=lat_rad.dtype)*lat_rad[0]; # Resize to match
            # END IF
            if( alt_size == 1 ):
                alt = np.ones(long_rad_size, dtype=alt.dtype)*alt[0]; # Resize to match
            # END IF
        elif( (alt_size != 1) and ((lat_rad_size == 1) or (long_rad_size == 1)) ):
            if( lat_rad_size == 1 ):
                lat_rad = np.ones(alt_size, dtype=lat_rad.dtype)*lat_rad[0]; # Resize to match
            # END IF
            if( long_rad_size == 1 ):
                long_rad = np.ones(alt_size, dtype=long_rad.dtype)*long_rad[0]; # Resize to match
            # END IF
        # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0; # Default
    # END IF
    
    # Convert as needed
    if( FLG_radIn_azel == False ):
        az_rad *= np.pi/180; # rad, convert from deg to rad
        el_rad *= np.pi/180; # rad, convert from deg to rad
    # END IF
    
    if( FLG_locECEF == False ):
        # Convert as needed
        if( FLG_radIn_latlong == False ):
            lat_rad *= np.pi/180; # rad, convert from deg to rad
            long_rad *= np.pi/180; # rad, convert from deg to rad
        # END IF
        X, Y, Z = conv_geo2ecef(lat_rad, long_rad, alt, FLG_radIn=True, FLG_kmIn=FLG_altkm, FLG_kmOut=FLG_rangekm, FLG_skipInputFiltering=True, FLG_overloaded = 0)
    else:
        X = np.copy(lat_rad); # Masquerade
        Y = np.copy(long_rad); # Masquerade
        Z = np.copy(alt); # Masquerade
        lat_rad, long_rad, alt = conv_ecef2geo(X, Y, Z, FLG_kmIn=FLG_altkm, FLG_kmOut=FLG_altkm, FLG_radOut=True, FLG_skipInputFiltering=True, FLG_overloaded=0); # FLG_altkm == true means that X,Y,Z masquerading are in km
        if( (FLG_rangekm == True) and (FLG_altkm == False) ): # Overload FLG_altkm to mean X/Y/Z is in km or not and that for a conversion
            X /= 1000; # km, convert from m to km, protect original variable
            Y /= 1000; # km, convert from m to km
            Z /= 1000; # km, convert from m to km
        # END IF
    # END IF

    # Precalc things
    latSin = np.sin(lat_rad); # Pre-calc
    latCos = np.cos(lat_rad); # Pre-calc
    longSin = np.sin(long_rad); # Pre-calc
    longCos = np.cos(long_rad); # Pre-calc
    
    # Circucous, AER -> ENU
    erangez = rangez * np.cos(el_rad); # Reused, so pre-calc
    E = erangez * np.sin(az_rad)
    N = erangez * np.cos(az_rad)
    U = rangez * np.sin(el_rad)
    # ENU -> ECEF
    (Xend, Yend, Zend) = np.einsum('ijk,jk->ik', np.array( ((-longSin, -latSin*longCos, latCos*longCos), (longCos, -latSin*longSin, latCos*longSin), (np.zeros(longCos.size, dtype=longCos.dtype), latCos, latSin)) ), np.array( (E, N, U) ), optimize = False) + np.array( (X, Y, Z) ); # Matrix mult via einsum to deal with that extra dim  (einsum optimize was slower, so this means this is optimal)   
        
    # Convert as needed
    if( FLG_kmOut and (FLG_rangekm == False) ):
        Xend /= 1000; # FLG_rangekm==False forces calcs into m, but output requested to be km
        Yend /= 1000;
        Zend /= 1000;
    elif( (FLG_kmOut == False) and FLG_rangekm ):
        Xend *= 1000; # FLG_rangekm==True forces calcs into km, but output requested to be m 
        Yend *= 1000;
        Zend *= 1000;
    # END IF
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return Xend, Yend, Zend
    elif( FLG_overloaded == 1 ):
        return np.vstack( (Xend, Yend) ).T, Zend
    elif( FLG_overloaded == 2 ):
        return np.vstack( (Xend, Yend, Zend) ).T
    # END IF
# END DEF

def conv_ecef2geo_looper_altOnly(X, Y, Z):
    # WGS-84 datum, ECEF XYZ -> lat/long
    # Very fast for scalars only, has no options, must be meters
    
    # WGS-84 datum, ECEF XYZ -> lat/long
    # Based on https://link.springer.com/article/10.1007/s00190-010-0419-x
    a = 6378137.; #m, WGS-84 equatorial radius (semi-major axis of elipsoid)
    # X, Y, Z MUST be in meters
    f = 1/298.257223563; #unitless, flattening factor (f = 1 - b/a)
    eSq = (2 - f)*f; # unitless, Eccentricity of elispoid (1 - b^2/a^2) squared
    eQuad = eSq*eSq; # unitless,  Eccentricity of elispoid (1 - b^2/a^2) quadrophenia
    rho = np.sqrt(X**2 + Y**2);
    p = (rho/a)**2;
    q = (1-eSq)/a**2*Z**2;
    r = (p + q - eQuad)/6;
    s = eQuad*p*q;
    sigmaSqrt= np.sqrt(8*r**3 + s);
    sSqrt = np.sqrt(s);
    u = r + 1/2*((sigmaSqrt + sSqrt)**(2/3) + (sigmaSqrt - sSqrt)**(2/3));
    v = np.sqrt(u**2 + eQuad*q);
    uNv = u + v;
    w = eSq*(uNv - q)/(2*v);
    k = (uNv)/(w + np.sqrt(w**2 + uNv));
    d = k*rho/(k + eSq);
    g = np.sqrt(d**2 + Z**2);
    alt = (k + eSq - 1)/k*g; # m, altitude result
    
    return alt
# END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def conv_aer2ecef_looper(azSin_elCos, azCos_elCos, elSin, rangez, latSin, latCos, longSin, longCos, X, Y, Z, alt):
    # Azimuth/Elevation/Range from location (lat/long/alt AND X/Y/Z in ECEF) -> ECEF X/Y/Z at point
    # For solving at altitude, has no options, must be rad/meters, uses pre-calc'd XYZ
    # Very fast for scalars only
    
    # Circucous, AER -> ENU
    E = rangez * azSin_elCos
    N = rangez * azCos_elCos
    U = rangez * elSin
    # ENU -> ECEF
    return -longSin*E - latSin*longCos*N + latCos*longCos*U + X, longCos*E - latSin*longSin*N + latCos*longSin*U + Y, latCos*N + latSin*U + Z
# END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def conv_aer2ecef_atAlt(az_rad, el_rad, alt_not_range, lat_rad, long_rad, alt, FLG_radIn_azel=False, FLG_rangekm=False, FLG_radIn_latlong=False, FLG_altkm=False, FLG_locECEF=False, FLG_kmOut=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=None):
    # Azimuth/Elevation/Range from location (lat/long/alt OR X/Y/Z in ECEF) -> ECEF X/Y/Z at point
    if( FLG_skipInputFiltering == False ):
        # Normalize inputs just in case
        if( (el_rad is not None) and (alt_not_range is not None) ):
            az_rad = normalizor(az_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            el_rad = normalizor(el_rad, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            alt_not_range = normalizor(alt_not_range, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (el_rad is None) and (alt_not_range is not None) ):
            # Secret overloaded az_rad == (az_rad, el_rad) while el_rad is None
            az_rad = normalizor(az_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            el_rad = np.ravel(az_rad[:, 1].copy()); # Get out the el_rad
            az_rad = np.ravel(az_rad[:, 0]); # Convert the az_rad
            alt_not_range = normalizor(alt_not_range, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded az_rad == (az_rad, el_rad, alt_not_range) while el_rad is None and alt_not_range is None
            az_rad = normalizor(az_rad, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            alt_not_range = np.ravel(alt_not_range[:, 2].copy()); # Get out the alt_not_range
            el_rad = np.ravel(az_rad[:, 1].copy()); # Get out the el_rad
            az_rad = np.ravel(az_rad[:, 0]); # Convert the az_rad
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        az_rad_size = az_rad.size; # Size of this
        el_rad_size = el_rad.size; # Size of that
        alt_not_range_size = alt_not_range.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (az_rad_size != 1) and ((el_rad_size == 1) or (alt_not_range_size == 1)) ):
            if( el_rad_size == 1 ):
                el_rad = np.ones(az_rad_size, dtype=el_rad.dtype)*el_rad[0]; # Resize to match
            # END IF
            if( alt_not_range_size == 1 ):
                alt_not_range = np.ones(az_rad_size, dtype=alt_not_range.dtype)*alt_not_range[0]; # Resize to match
            # END IF
        elif( (el_rad_size != 1) and ((az_rad_size == 1) or (alt_not_range_size == 1)) ):
            if( az_rad_size == 1 ):
                az_rad = np.ones(el_rad_size, dtype=az_rad.dtype)*az_rad[0]; # Resize to match
            # END IF
            if( alt_not_range_size == 1 ):
                alt_not_range = np.ones(el_rad_size, dtype=alt_not_range.dtype)*alt_not_range[0]; # Resize to match
            # END IF
        elif( (alt_not_range_size != 1) and ((az_rad_size == 1) or (el_rad_size == 1)) ):
            if( az_rad_size == 1 ):
                az_rad = np.ones(alt_not_range_size, dtype=az_rad.dtype)*az_rad[0]; # Resize to match
            # END IF
            if( el_rad_size == 1 ):
                el_rad = np.ones(alt_not_range_size, dtype=el_rad.dtype)*el_rad[0]; # Resize to match
            # END IF
        # END IF
        
        # Normalize inputs just in case
        if( (long_rad is not None) and (alt is not None) ):
            lat_rad = normalizor(lat_rad, 1); # Normalize the input for the functions
            long_rad = normalizor(long_rad, 1); # Normalize the input for the functions
            alt = normalizor(alt, 1); # Normalize the input for the functions
        elif( (long_rad is None) and (alt is not None) ):
            # Secret overloaded lat_rad == (lat_rad, long_rad) while long_rad is None
            lat_rad = normalizor(lat_rad, 2); # Normalize the input for the functions
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            alt = normalizor(alt, 1); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded lat_rad == (lat_rad, long_rad, alt) while long_rad is None and alt is None
            lat_rad = normalizor(lat_rad, 2); # Normalize the input for the functions
            alt = np.ravel(lat_rad[:, 2].copy()); # Get out the alt
            long_rad = np.ravel(lat_rad[:, 1].copy()); # Get out the long_rad
            lat_rad = np.ravel(lat_rad[:, 0]); # Convert the lat_rad
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        lat_rad_size = lat_rad.size; # Size of this
        long_rad_size = long_rad.size; # Size of that
        alt_size = alt.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (lat_rad_size != 1) and ((long_rad_size == 1) or (alt_size == 1)) ):
            if( long_rad_size == 1 ):
                long_rad = np.ones(lat_rad_size, dtype=long_rad.dtype)*long_rad[0]; # Resize to match
            # END IF
            if( alt_size == 1 ):
                alt = np.ones(lat_rad_size, dtype=alt.dtype)*alt[0]; # Resize to match
            # END IF
        elif( (long_rad_size != 1) and ((lat_rad_size == 1) or (alt_size == 1)) ):
            if( lat_rad_size == 1 ):
                lat_rad = np.ones(long_rad_size, dtype=lat_rad.dtype)*lat_rad[0]; # Resize to match
            # END IF
            if( alt_size == 1 ):
                alt = np.ones(long_rad_size, dtype=alt.dtype)*alt[0]; # Resize to match
            # END IF
        elif( (alt_size != 1) and ((lat_rad_size == 1) or (long_rad_size == 1)) ):
            if( lat_rad_size == 1 ):
                lat_rad = np.ones(alt_size, dtype=lat_rad.dtype)*lat_rad[0]; # Resize to match
            # END IF
            if( long_rad_size == 1 ):
                long_rad = np.ones(alt_size, dtype=long_rad.dtype)*long_rad[0]; # Resize to match
            # END IF
        # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0; # Default
    # END IF
    
    # Convert as needed
    if( FLG_radIn_azel == False ):
        az_rad *= np.pi/180; # rad, convert from deg to rad
        el_rad *= np.pi/180; # rad, convert from deg to rad
    # END IF
    
    if( FLG_locECEF == False ):
        # Convert as needed
        if( FLG_radIn_latlong == False ):
            lat_rad *= np.pi/180; # rad, convert from deg to rad
            long_rad *= np.pi/180; # rad, convert from deg to rad
        # END IF
        X, Y, Z = conv_geo2ecef(lat_rad, long_rad, alt, FLG_radIn=True, FLG_kmIn=FLG_altkm, FLG_kmOut=FLG_rangekm, FLG_skipInputFiltering=True, FLG_overloaded = 0)
    else:
        X = np.copy(lat_rad); # Masquerade
        Y = np.copy(long_rad); # Masquerade
        Z = np.copy(alt); # Masquerade
        lat_rad, long_rad, alt = conv_ecef2geo(X, Y, Z, FLG_kmIn=FLG_altkm, FLG_kmOut=FLG_altkm, FLG_radOut=True, FLG_skipInputFiltering=True, FLG_overloaded=0); # FLG_altkm == true means that X,Y,Z masquerading are in km
        if( (FLG_rangekm == True) and (FLG_altkm == False) ): # Overload FLG_altkm to mean X/Y/Z is in km or not and that for a conversion
            X /= 1000; # km, convert from m to km, protect original variable
            Y /= 1000; # km, convert from m to km
            Z /= 1000; # km, convert from m to km
        # END IF
    # END IF
    
    # Solve for the correct range that yields the altitude requested
    # There's definitely a way to do this directly !IN ELLIPSOID LAND!, but it looks pretty hard. This was ez pz
    # def altFinder(rangeEst, az_rad, el_rad, lat_rad, long_rad, alt, alt_not_range): # Solve function to solve
    #     Xest, Yest, Zest = conv_aer2ecef(az_rad, el_rad, rangeEst, lat_rad, long_rad, alt, FLG_radIn_azel=True, FLG_radIn_latlong=True, FLG_skipInputFiltering=True, FLG_overloaded=0);
    #     _, _, altEst = conv_ecef2geo(Xest, Yest, Zest, FLG_skipInputFiltering=True, FLG_overloaded=0);
        
    #     return np.abs(altEst - alt_not_range)
    # # END DEF
    def altFinder(rangeEst, azSin_elCos, azCos_elCos, elSin, latSin, latCos, longSin, longCos, X, Y, Z, alt, alt_not_range): # Solve function to solve
        Xest, Yest, Zest = conv_aer2ecef_looper(azSin_elCos, azCos_elCos, elSin, rangeEst, latSin, latCos, longSin, longCos, X, Y, Z, alt);
        
        return np.abs(conv_ecef2geo_looper_altOnly(Xest, Yest, Zest) - alt_not_range)
    # END DEF
    
    # Precalc
    latSin = np.sin(lat_rad).item();
    latCos = np.cos(lat_rad).item();
    longSin = np.sin(long_rad).item(); 
    longCos = np.cos(long_rad).item();
    azSin_elCos = np.sin(az_rad)*np.cos(el_rad);
    azCos_elCos = np.cos(az_rad)*np.cos(el_rad);
    elSin = np.sin(el_rad);
    X = X.item();
    Y = Y.item();
    Z = Z.item();
    alt = alt.item();
    rangeEst = alt_not_range/elSin;
    
    for i in range(0, len(rangeEst)):
        # rangeEst_solved = minimize(altFinder, np.array((rangeEst[i],)), args=(np.array((az_rad[i],)), np.array((el_rad[i],)), lat_rad, long_rad, alt, np.array((alt_not_range[i],))), method='nelder-mead', options={'xatol': 1e-8, 'disp': False}); # Use a scipy minimizer
        # rangeEst[i] = minimize(altFinder, rangeEst[i], args=(azSin_elCos[i], azCos_elCos[i], elSin[i], latSin, latCos, longSin, longCos, X, Y, Z, alt, alt_not_range[i]), method='COBYLA', options={'disp': False}).x; # Use a scipy minimizer & get the solution out
        # rangeEst[i] = least_squares(altFinder, rangeEst[i], args=(azSin_elCos[i], azCos_elCos[i], elSin[i], latSin, latCos, longSin, longCos, X, Y, Z, alt, alt_not_range[i]), jac='2-point', method='lm' ).x; # Use scipy least squares b/c this is simple and it is faster
        rangeEst[i] = leastsq(altFinder, rangeEst[i], args=(azSin_elCos[i], azCos_elCos[i], elSin[i], latSin, latCos, longSin, longCos, X, Y, Z, alt, alt_not_range[i]))[0].item(); # Use scipy legacy least squares b/c it's very bare-bones and ideal for this
    # END FOR i
    Xend, Yend, Zend = conv_aer2ecef(az_rad, el_rad, rangeEst, lat_rad, long_rad, alt, FLG_radIn_azel=True, FLG_radIn_latlong=True, FLG_skipInputFiltering=True, FLG_overloaded=0)

    # Convert as needed
    if( FLG_kmOut and (FLG_rangekm == False) ):
        Xend /= 1000; # FLG_rangekm==False forces calcs into m, but output requested to be km
        Yend /= 1000;
        Zend /= 1000;
    elif( (FLG_kmOut == False) and FLG_rangekm ):
        Xend *= 1000; # FLG_rangekm==True forces calcs into km, but output requested to be m 
        Yend *= 1000;
        Zend *= 1000;
    # END IF
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return Xend, Yend, Zend
    elif( FLG_overloaded == 1 ):
        return np.vstack( (Xend, Yend) ).T, Zend
    elif( FLG_overloaded == 2 ):
        return np.vstack( (Xend, Yend, Zend) ).T
    # END IF
# END DEF

def conv_ecef2aer(X, Y, Z, Xend, Yend, Zend, FLG_kmIn=False, FLG_kmOut=False, FLG_radOut=False, FLG_skipInputFiltering=False, FLG_protectInputs=True, FLG_overloaded=None):
    # X/Y/Z start and X/Y/Z end in ECEF -> Azimuth/Elevation/Range from location (X/Y/Z start in ECEF)
    if( FLG_skipInputFiltering == False ):
        # Normalize inputs just in case
        if( (Y is not None) and (Z is not None) ):
            X = normalizor(X, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Y = normalizor(Y, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Z = normalizor(Z, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (Y is None) and (Z is not None) ):
            # Secret overloaded X == (X, Y) while Y is None
            X = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Y = np.ravel(X[:, 1].copy()); # Get out the Y
            X = np.ravel(X[:, 0]); # Convert the X
            Z = normalizor(Z, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded X == (X, Y, Z) while Y is None and Z is None
            # raise Exception('ERROR in positizer.conv_geo2ecef: I didn\'t code positizer.normalizor for 3 dims yet, sorry.');
            X = normalizor(X, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Z = np.ravel(X[:, 2].copy()); # Get out the Z
            Y = np.ravel(X[:, 1].copy()); # Get out the Y
            X = np.ravel(X[:, 0]); # Convert the X
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        X_size = X.size; # Size of this
        Y_size = Y.size; # Size of that
        Z_size = Z.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (X_size != 1) and ((Y_size == 1) or (Z_size == 1)) ):
            if( Y_size == 1 ):
                Y = np.ones(X_size, dtype=Y.dtype)*Y[0]; # Resize to match
            # END IF
            if( Z_size == 1 ):
                Z = np.ones(X_size, dtype=Z.dtype)*Z[0]; # Resize to match
            # END IF
        elif( (Y_size != 1) and ((X_size == 1) or (Z_size == 1)) ):
            if( X_size == 1 ):
                X = np.ones(Y_size, dtype=X.dtype)*X[0]; # Resize to match
            # END IF
            if( Z_size == 1 ):
                Z = np.ones(Y_size, dtype=Z.dtype)*Z[0]; # Resize to match
            # END IF
        elif( (Z_size != 1) and ((X_size == 1) or (Y_size == 1)) ):
            if( X_size == 1 ):
                X = np.ones(Z_size, dtype=X.dtype)*X[0]; # Resize to match
            # END IF
            if( Y_size == 1 ):
                Y = np.ones(Z_size, dtype=Y.dtype)*Y[0]; # Resize to match
            # END IF
        # END IF
        
        # Normalize inputs just in case
        if( (Yend is not None) and (Zend is not None) ):
            Xend = normalizor(Xend, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Yend = normalizor(Yend, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Zend = normalizor(Zend, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
        elif( (Yend is None) and (Zend is not None) ):
            # Secret overloaded X == (X, Y) while Y is None
            Xend = normalizor(Xend, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Yend = np.ravel(Xend[:, 1].copy()); # Get out the Y
            Xend = np.ravel(Xend[:, 0]); # Convert the X
            Zend = normalizor(Zend, 1, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 1; # Flag for secret overloading so output is the same as input
            # END IF
        else:
            # Secret overloaded X == (X, Y, Z) while Y is None and Z is None
            # raise Exception('ERROR in positizer.conv_geo2ecef: I didn\'t code positizer.normalizor for 3 dims yet, sorry.');
            Xend = normalizor(Xend, 2, FLG_protectInputs=FLG_protectInputs); # Normalize the input for the functions
            Zend = np.ravel(Xend[:, 2].copy()); # Get out the Z
            Yend = np.ravel(Xend[:, 1].copy()); # Get out the Y
            Xend = np.ravel(Xend[:, 0]); # Convert the X
            if( FLG_overloaded is None ): # Only change the overload if it wasn't specified
                FLG_overloaded = 2; # Flag for secret overloading so output is the same as input
            # END IF
        # END IF
        
        # Get the sizes
        Xend_size = Xend.size; # Size of this
        Yend_size = Yend.size; # Size of that
        Zend_size = Zend.size; # Size of this too
        
        # Replicate size 1 ones to fit the rest
        if( (Xend_size != 1) and ((Yend_size == 1) or (Zend_size == 1)) ):
            if( Yend_size == 1 ):
                Yend = np.ones(Xend_size, dtype=Yend.dtype)*Yend[0]; # Resize to match
            # END IF
            if( Zend_size == 1 ):
                Zend = np.ones(Xend_size, dtype=Zend.dtype)*Zend[0]; # Resize to match
            # END IF
        elif( (Yend_size != 1) and ((Xend_size == 1) or (Zend_size == 1)) ):
            if( Xend_size == 1 ):
                Xend = np.ones(Yend_size, dtype=Xend.dtype)*Xend[0]; # Resize to match
            # END IF
            if( Zend_size == 1 ):
                Zend = np.ones(Yend_size, dtype=Zend.dtype)*Zend[0]; # Resize to match
            # END IF
        elif( (Zend_size != 1) and ((Xend_size == 1) or (Yend_size == 1)) ):
            if( Xend_size == 1 ):
                Xend = np.ones(Zend_size, dtype=Xend.dtype)*Xend[0]; # Resize to match
            # END IF
            if( Yend_size == 1 ):
                Yend = np.ones(Zend_size, dtype=Yend.dtype)*Yend[0]; # Resize to match
            # END IF
        # END IF
        
        # if( (X_size == 1) and (Xend_size != 1) ):
        #     X = np.ones(Xend_size, dtype=X.dtype)*X[0]; # Resize to match
        #     Y = np.ones(Yend_size, dtype=Y.dtype)*Y[0]; # Resize to match
        #     Z = np.ones(Zend_size, dtype=Z.dtype)*Z[0]; # Resize to match
        # elif( (X_size != 1) and (Xend_size == 1) ):
        #     Xend = np.ones(X_size, dtype=Xend.dtype)*Xend[0]; # Resize to match
        #     Yend = np.ones(Y_size, dtype=Yend.dtype)*Yend[0]; # Resize to match
        #     Zend = np.ones(Z_size, dtype=Zend.dtype)*Zend[0]; # Resize to match
        # # END IF
    # END IF
    if( FLG_overloaded is None ):
        FLG_overloaded = 0;
    # END IF
    
    # Convert as needed
    if( FLG_kmIn == True ):
        X *= 1000; # m, km -> m conversion
        Y *= 1000;
        Z *= 1000;
        Xend *= 1000; # m, km -> m conversion
        Yend *= 1000;
        Zend *= 1000;
    # END IF
    
    lat_rad, long_rad, alt = conv_ecef2geo(X, Y, Z, FLG_kmIn=False, FLG_kmOut=False, FLG_radOut=True, FLG_skipInputFiltering=True, FLG_overloaded=0); # FLG_altkm == true means that X,Y,Z masquerading are in km

    # Precalc things
    latSin = np.sin(lat_rad); # Pre-calc
    latCos = np.cos(lat_rad); # Pre-calc
    longSin = np.sin(long_rad); # Pre-calc
    longCos = np.cos(long_rad); # Pre-calc
    
    # Circucous, ECEF -> ENU
    (E, N, U) = np.einsum('ijk,jk->ik', np.array( ((-longSin, longCos, np.zeros(longCos.size, dtype=longCos.dtype)), (-latSin*longCos, -latSin*longSin, latCos), (latCos*longCos, latCos*longSin, latSin)) ), (np.array( (Xend, Yend, Zend) ) - np.array( (X, Y, Z) )), optimize = False); # Matrix mult via einsum to deal with that extra dim, optimize false b/c was optimal already and optimizing slowed it
    
    # ENU -> AER
    rangez = np.sqrt( E**2 + N**2 + U**2 );
    el_rad = np.arcsin(U/rangez);
    az_rad = np.arctan2(E, N); # NEEDS arctan2
    az_rad[az_rad < 0] += 2*np.pi; # Make positive
    
    # Convert as needed
    if( FLG_kmOut ):
        rangez /= 1000;
    # END IF
    if( FLG_radOut == False ):
        az_rad *= 180/np.pi;
        el_rad *= 180/np.pi;
    # END IF
    
    # Output using secret FLG_overloaded flag so output is the same as input
    if( FLG_overloaded == 0 ):
        return az_rad, el_rad, rangez
    elif( FLG_overloaded == 1 ):
        return np.vstack( (az_rad, el_rad) ).T, rangez
    elif( FLG_overloaded == 2 ):
        return np.vstack( (az_rad, el_rad, rangez) ).T
    # END IF
# END DEF