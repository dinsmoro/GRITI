#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Spherical, linear (it's that - in places) RBF interpolation
#distance is spherical calc'd
#RBF kernel is linear (-eps*r) with eps = 1
#rbf_weights 1st to get weights, then feed it into rbf_interp
#rbf_interp's data_lat/data_long should be the same as weight's data_lat/data_long
#data_lat/data_long/data_orig should all be vectors of lat/long/data combos

#spherical units must be radians for lat and long
#spherical distance is approximated via the great circle function

import numpy as np
from numba import jit
from scipy.linalg import lu_factor, lu_solve

# ========= Spherical Distance (lat/long at same altitude) =========
# Orchastrator for spherical RBF, includes a batching alg to approximate closeness without me having to solve closeness ktree stuff
def rbf_interpSpherical(pt_lat, pt_long, data_lat, data_long, data_orig, FLG_weights=False, nearPts=None, FLG_nearIndv=True, FLG_batchum=False, FLG_inputRad=False):
    # Alignment
    if( pt_long is None ):
        if( pt_lat.shape[0] == 2 ):
            pt_long = pt_lat[1, :].copy(); # Get out the long from the overloaded lat
            pt_lat = pt_lat[0, :]; # Clip to just the lat
        else:
            pt_long = pt_lat[:, 1].copy(); # Get out the long from the overloaded lat
            pt_lat = pt_lat[:, 0]; # Clip to just the lat
        # END IF
    else:
        if( pt_lat.ndim == 2 ):
            pt_lat = pt_lat.ravel(); # Reduce an unneeded dimension
        # END IF
        if( pt_long.ndim == 2 ):
            pt_long = pt_long.ravel(); # Reduce an unneeded dimension
        # END IF
    # END IF
    if( data_long is None ):
        if( data_lat.shape[0] == 2 ):
            data_long = data_lat[1, :].copy(); # Get out the long from the overloaded lat
            data_lat = data_lat[0, :]; # Clip to just the lat
        else:
            data_long = data_lat[:, 1].copy(); # Get out the long from the overloaded lat
            data_lat = data_lat[:, 0]; # Clip to just the lat
        # END IF
    else:
        if( data_lat.ndim == 2 ):
            data_lat = data_lat.ravel(); # Reduce an unneeded dimension
        # END IF
        if( data_long.ndim == 2 ):
            data_long = data_long.ravel(); # Reduce an unneeded dimension
        # END IF
    # END IF
    if( FLG_inputRad == False ):
        pt_lat = pt_lat*np.pi/180; # rad, convert from deg
        pt_long = pt_long*np.pi/180; # rad, convert from deg
        data_lat = data_lat*np.pi/180; # rad, convert from deg
        data_long = data_long*np.pi/180; # rad, convert from deg
    # END IF
    FLG_ndimRevert = False; # Prime a flag
    if( data_orig.ndim == 1 ):
        FLG_ndimRevert = True; # Set a flag
        pt_interp = np.empty( (pt_lat.size, 1), dtype=data_orig.dtype); # Preallocate
        data_orig = np.expand_dims(data_orig, axis=-1); # Add on an axis to work with it better
    else:
        pt_interp = np.empty( (pt_lat.size, data_orig.shape[1]), dtype=data_orig.dtype); # Preallocate
    # END IF
    
    # Pre-calc the distances because they're needed for later anyway
    distances = distance_calc_Spherical_noNumba(pt_lat, pt_long, data_lat, data_long)
    
    if( (nearPts is None) or (nearPts is False) ):
        if( FLG_weights == False ):
            lu_factor_precalc = lu_factor(distance_calc_identicalSpherical(data_lat, data_long), overwrite_a=True, check_finite=False); # Precalculate
            for j in range(0, pt_interp.shape[1]):
                pt_interp[:, j] = rbf_interp_applySpherical( distances, None, None, None, \
                    lu_solve(lu_factor_precalc, -data_orig[:, j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
            # END FOR j
        else:
            # In this case data_orig is masquerading as weights
            for j in range(0, pt_interp.shape[1]):
                pt_interp[:, j] = rbf_interp_applySpherical( distances, None, None, None, data_orig[:, :, j]); #less powerful function call
            # END FOR j
        # END IF
    else:
        # from positizer import positizer_dist
        if( FLG_weights == True ):
            raise Exception('ERROR in subfun_rbf: FLG_weights is TRUE which indicates that data_orig is actually a weight matrix, but that is NOT supported by nearPts != None. Crashing.');
        # END IF
        
        if( FLG_batchum == False ):
            if( nearPts >= data_lat.size ):
                # If data_lat.size is smaller than nearPts, using all data_lat/data_long pts of data_orig
                lu_factor_precalc = lu_factor(distance_calc_identicalSpherical(data_lat, data_long), overwrite_a=True, check_finite=False); # Precalculate
                for j in range(0, pt_interp.shape[1]):
                    pt_interp[:, j] = rbf_interp_applySpherical( distances, None, None, None, \
                        lu_solve(lu_factor_precalc, -data_orig[:,j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                # END FOR j
            else:
                # Goal here is to limit the rbf search to the relevant indexes while still being able to matrix math the whole thing, more memory but ideally less than the whole thing
                if( FLG_nearIndv == False ): # This combines all nearby points together to do it once [slower]
                    relevant_indxz = np.zeros( data_lat.size, dtype=np.bool_); # Logical array for the relevant indexz to work with
                    for k in range(0, pt_lat.size):
                        # relevant_dist = positizer_dist( (pt_lat[k], pt_long[k]), (data_lat, data_long), solverOverride='haversine', FLG_inputRad=True ); # Get the distance between this location and the pts in the grid/pile of pts
                        relevant_dist = distances[k, :]; # Get the distance between this location and the pts in the grid/pile of pts
                        relevant_indxz[np.sort(np.argpartition(relevant_dist, nearPts)[:nearPts])] = True; # Update the nearPts# of relevant indexes
                    # END FOR k
                    # rbf only on relevant indexes based on nearPts# specifier
                    lu_factor_precalc = lu_factor(distance_calc_identicalSpherical(data_lat[relevant_indxz], data_long[relevant_indxz]), overwrite_a=True, check_finite=False); # Precalculate
                    for j in range(0, pt_interp.shape[1]):
                        pt_interp[:, j] = rbf_interp_applySpherical( pt_lat, pt_long, data_lat[relevant_indxz], data_long[relevant_indxz], \
                            lu_solve(lu_factor_precalc, -data_orig[relevant_indxz, j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                    # END FOR j
                else:
                    relevant_indxz = np.zeros( data_lat.size, dtype=np.bool_); # Logical array for the relevant indexz to work with
                    for k in range(0, pt_lat.size):
                        # relevant_dist = positizer_dist( (pt_lat[k], pt_long[k]), (data_lat, data_long), solverOverride='haversine', FLG_inputRad=True ); # Get the distance between this location and the pts in the grid/pile of pts
                        relevant_dist = distances[k, :]; # Get the distance between this location and the pts in the grid/pile of pts
                        relevant_indxz[np.sort(np.argpartition(relevant_dist, nearPts)[:nearPts])] = True; # Update the nearPts# of relevant indexes
                        # rbf only on relevant indexes based on nearPts# specifier
                        lu_factor_precalc = lu_factor(distance_calc_identicalSpherical(data_lat[relevant_indxz], data_long[relevant_indxz]), overwrite_a=True, check_finite=False); # Precalculate
                        for j in range(0, pt_interp.shape[1]):
                            pt_interp[k, j] = rbf_interp_applySpherical( pt_lat[k], pt_long[k], data_lat[relevant_indxz], data_long[relevant_indxz], \
                                lu_solve(lu_factor_precalc, -data_orig[relevant_indxz, j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                        # END FOR j
                        relevant_indxz &= False; # Reset
                    # END FOR k
                # END IF
            # END IF
        else:
            # Goal here is to batch shared closest pts together for speed, this saves the most memory but is slower because batches are small so matrix-math-parallelism is small
            batchy_indxz = np.sort(np.argpartition(distances, nearPts, axis=1)[:,:nearPts],axis=1); # Get the nearest distances
            batchy_indxz_unique, batchy_indxz_uniqueIndxz = np.unique(batchy_indxz, axis=0, return_inverse=True); # Get the uniques
            
            for k in range(0, batchy_indxz_unique.shape[0]):
                pt_batchr = np.where(batchy_indxz_uniqueIndxz == k)[0]; # Get where pt indexes line up to the current batch
                
                lu_factor_precalc = lu_factor(distance_calc_identicalSpherical(data_lat[batchy_indxz_unique[k,:]], data_long[batchy_indxz_unique[k,:]]), overwrite_a=True, check_finite=False); # Precalculate
                for j in range(0, pt_interp.shape[1]):
                    pt_interp[pt_batchr, j] = rbf_interp_applySpherical( pt_lat[pt_batchr], pt_long[pt_batchr], data_lat[batchy_indxz_unique[k,:]], data_long[batchy_indxz_unique[k,:]], \
                        lu_solve(lu_factor_precalc, -data_orig[batchy_indxz_unique[k,:], j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                # END FOR j
            # END FOR k
        # END IF
    # END IF
    
    if( FLG_ndimRevert == True ):
        pt_interp = pt_interp.ravel(); # Drop to a vector as the input was a vector but the loop supports 2D so it was 2D's preemptively
    # END IF
    
    return pt_interp
# END DEF

def rbf_interp_applySpherical(pt_lat, pt_long, data_lat, data_long, weights):
    if( pt_long is None ):
        # pt_lat is actually overloaded and is the distnace matrix already
        return pt_lat @ -weights
    else:
        #-- great circle --
        return np.arccos( np.sin(pt_lat.reshape(-1,1))*np.sin(data_lat.reshape(1,-1)) + np.cos(pt_lat.reshape(-1,1))*np.cos(data_lat.reshape(1,-1))*np.cos(data_long.reshape(1,-1)-pt_long.reshape(-1,1)) ) @ -weights
    # END IF
# END DEF

#essential for weighting calcs
@jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def distance_calc_identicalSpherical(lat1, long1):
    dist_identical = np.zeros((lat1.size,lat1.size),dtype=lat1.dtype);
    for jj in range(0,lat1.size-1):
        #this reduces the effort for identical vector inputs by 1/2 or something like that
        # dist_identical[jj+1:,jj] = np.sqrt(np.sin((lat1[jj]-lat1[jj+1:])/2)**2 + np.cos(lat1[jj]) * np.cos(lat1[jj+1:]) * np.sin((long1[jj]-long1[jj+1:])/2)**2); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        # dist_identical[jj+1:,jj] = np.sqrt(np.square(np.sin((lat1[jj]-lat1[jj+1:])/2)) + np.cos(lat1[jj]) * np.cos(lat1[jj+1:]) * np.square(np.sin((long1[jj]-long1[jj+1:])/2))); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        dist_identical[jj+1:,jj] = np.arccos( np.sin(lat1[jj])*np.sin(lat1[jj+1:]) + np.cos(lat1[jj])*np.cos(lat1[jj+1:])*np.cos(long1[jj]-long1[jj+1:]) ); #great circle
        dist_identical[jj,jj+1:] = dist_identical[jj+1:,jj]; #mirror it, faster than adding the transpose
    #END FOR jj
    return dist_identical
#END DEF

def distance_calc_identicalSpherical_noNumba(lat1, long1):
    # Numba-free option, slightly slower than numba:
    sinz = np.sin(lat1);
    cosz = np.cos(lat1);
    dist_identical = np.arccos(np.multiply.outer(sinz, sinz) + np.multiply.outer(cosz, cosz)*np.cos(np.subtract.outer(long1, long1)));
    np.fill_diagonal(dist_identical, 0); # THis is in-place
    # # Numba-free option that doesn't redo work, but is somehow slowest (prob slicing stuff?):
    # ar1, ar2 = np.triu_indices(lat1.size, k=1)
    # dist_identical = np.zeros( (lat1.size, lat1.size), dtype=np.float64); # Prep
    # dist_identical[ar1, ar2] = np.arccos( sinz[ar1]*sinz[ar2] + cosz[ar1]*cosz[ar2]*np.cos(long1[ar1] - long1[ar2]) )
    # al1, al2 = np.tril_indices(lat1.size, k=-1)
    # dist_identical.T[ar1, ar2] = dist_identical[ar1, ar2]
    return dist_identical
# END DEF

#essential for weighting calcs
@jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def distance_calc_Spherical(lat1, long1, lat2, long2):
    dist_identical = np.empty((lat1.size,lat2.size),dtype=lat1.dtype);
    for jj in range(0,lat1.size):            #this reduces the effort for identical vector inputs by 1/2 or something like that
        # dist_identical[jj+1:,jj] = np.sqrt(np.sin((lat1[jj]-lat2/2)**2 + np.cos(lat1[jj]) * np.cos(lat2[) * np.sin((long1[jj]-long2)/2)**2); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        # dist_identical[jj+1:,jj] = np.sqrt(np.square(np.sin((lat1[jj]-lat2)/2)) + np.cos(lat1[jj]) * np.cos(lat2) * np.square(np.sin((long1[jj]-long2)/2))); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        dist_identical[jj,:] = np.arccos( np.sin(lat1[jj])*np.sin(lat2) + np.cos(lat1[jj])*np.cos(lat2)*np.cos(long1[jj]-long2) ); #great circle
    # END FOR jj
    return dist_identical
#END DEF

def distance_calc_Spherical_noNumba(lat1, long1, lat2, long2):
    # Numba-free option, slightly faster than numba:
    sinz1 = np.sin(lat1);
    cosz1 = np.cos(lat1);
    sinz2 = np.sin(lat2);
    cosz2 = np.cos(lat2);
    dist_identical = np.arccos(np.multiply.outer(sinz1, sinz2) + np.multiply.outer(cosz1, cosz2)*np.cos(np.subtract.outer(long1, long2)));
    
    return dist_identical
# END DEF


# def rbf_interp_weights(data_lat, data_long, data_orig):
#     # Remember the -data_orig represents the kernel which is -K, where K is 1 here
#     return lu_solve(lu_factor(distance_calc_identicalSpherical(data_lat, data_long), overwrite_a=True, check_finite=False), -data_orig, trans=0, overwrite_b=True, check_finite=False); #most powerful function call
# #END DEF


# ========= Euclidean Distance (X/Y or X/Y/Z) =========
# Orchastrator for spherical RBF, includes a batching alg to approximate closeness without me having to solve closeness ktree stuff
def rbf_interp(pt_XYZ, data_XYZ, data_orig, FLG_weights=False, nearPts=None, FLG_nearIndv=True, FLG_batchum=True, FLG_inputRad=False):
    # Alignment
    FLG_ndimRevert = False; # Prime a flag
    if( data_orig.ndim == 1 ):
        FLG_ndimRevert = True; # Set a flag
        pt_interp = np.empty( (pt_XYZ.shape[0], 1), dtype=data_orig.dtype); # Preallocate
        data_orig = np.expand_dims(data_orig, axis=-1); # Add on an axis to work with it better
    else:
        pt_interp = np.empty( (pt_XYZ.shape[0], data_orig.shape[1]), dtype=data_orig.dtype); # Preallocate
    # END IF
    
    # Pre-calc the distances because they're needed for later anyway
    distances = np.zeros( (pt_XYZ.shape[0], data_XYZ.shape[0]), dtype=data_XYZ.dtype) # Preallocate to hold euclidean distance
    for i in range(0, data_XYZ.shape[1]): # Roll through each dim
        distances += np.subtract.outer(pt_XYZ[:, i], data_XYZ[:, i])**2 # Build euclidean distance dim-by-dim
    distances = np.sqrt(distances); # Sqrt to get the distance
    
    if( (nearPts is None) or (nearPts is False) ):
        if( FLG_weights == False ):
            lu_factor_precalc = lu_factor(distance_calc_identical_noNumba(data_XYZ), overwrite_a=True, check_finite=False); # Precalculate
            for j in range(0, pt_interp.shape[1]):
                pt_interp[:, j] = rbf_interp_apply( distances, None, \
                    lu_solve(lu_factor_precalc, -data_orig[:, j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
            # END FOR j
        else:
            # In this case data_orig is masquerading as weights
            for j in range(0, pt_interp.shape[1]):
                pt_interp[:, j] = rbf_interp_apply( distances, None, data_orig[:, :, j]); #less powerful function call
            # END FOR j
        # END IF
    else:
        if( FLG_weights == True ):
            raise Exception('ERROR in subfun_rbf: FLG_weights is TRUE which indicates that data_orig is actually a weight matrix, but that is NOT supported by nearPts != None. Crashing.');
        # END IF
        
        if( nearPts >= data_XYZ.shape[0] ):
            # If data_lat.size is smaller than nearPts, using all data_lat/data_long pts of data_orig
            lu_factor_precalc = lu_factor(distance_calc_identical_noNumba(data_XYZ), overwrite_a=True, check_finite=False); # Precalculate
            for j in range(0, pt_interp.shape[1]):
                pt_interp[:, j] = rbf_interp_apply( distances, None, \
                    lu_solve(lu_factor_precalc, -data_orig[:,j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
            # END FOR j
        else:
            if( FLG_batchum == False ):
                # Goal here is to limit the rbf search to the relevant indexes while still being able to matrix math the whole thing, more memory but ideally less than the whole thing
                if( FLG_nearIndv == False ): # This combines all nearby points together to do it once [slower]
                    relevant_indxz = np.zeros( data_XYZ.shape[0], dtype=np.bool_); # Logical array for the relevant indexz to work with
                    for k in range(0, pt_XYZ.shape[0]):
                        relevant_dist = distances[k, :]; # Get the distance between this location and the pts in the grid/pile of pts
                        relevant_indxz[np.sort(np.argpartition(relevant_dist, nearPts)[:nearPts])] = True; # Update the nearPts# of relevant indexes
                    # END FOR k
                    # rbf only on relevant indexes based on nearPts# specifier
                    lu_factor_precalc = lu_factor(distance_calc_identical_noNumba(data_XYZ[relevant_indxz, :]), overwrite_a=True, check_finite=False); # Precalculate
                    for j in range(0, pt_interp.shape[1]):
                        pt_interp[:, j] = rbf_interp_apply( pt_XYZ, data_XYZ[relevant_indxz, :], \
                            lu_solve(lu_factor_precalc, -data_orig[relevant_indxz, j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                    # END FOR j
                else:
                    # This does each point individually with its N nearest pts instead of combining them together
                    relevant_indxz = np.zeros( data_XYZ.shape[0], dtype=np.bool_); # Logical array for the relevant indexz to work with
                    for k in range(0, pt_XYZ.shape[0]):
                        relevant_dist = distances[k, :]; # Get the distance between this location and the pts in the grid/pile of pts
                        relevant_indxz[np.sort(np.argpartition(relevant_dist, nearPts)[:nearPts])] = True; # Update the nearPts# of relevant indexes
                    
                        # rbf only on relevant indexes based on nearPts# specifier
                        lu_factor_precalc = lu_factor(distance_calc_identical_noNumba(data_XYZ[relevant_indxz, :]), overwrite_a=True, check_finite=False); # Precalculate
                        for j in range(0, pt_interp.shape[1]):
                            pt_interp[k, j] = rbf_interp_apply( np.expand_dims(pt_XYZ[k, :], axis=0), data_XYZ[relevant_indxz, :], \
                                lu_solve(lu_factor_precalc, -data_orig[relevant_indxz, j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                        # END FOR j
                        
                        relevant_indxz &= False; # Reset
                    # END FOR k
                # END IF
            else:
                # Goal here is to batch shared closest pts together for speed, this saves the most memory but is slower because batches are small so matrix-math-parallelism is small
                batchy_indxz = np.sort(np.argpartition(distances, nearPts, axis=1)[:,:nearPts],axis=1); # Get the nearest distances
                batchy_indxz_unique, batchy_indxz_uniqueIndxz = np.unique(batchy_indxz, axis=0, return_inverse=True); # Get the uniques
                
                for k in range(0, batchy_indxz_unique.shape[0]):
                    pt_batchr = np.where(batchy_indxz_uniqueIndxz == k)[0]; # Get where pt indexes line up to the current batch
                    
                    lu_factor_precalc = lu_factor(distance_calc_identical_noNumba(data_XYZ[batchy_indxz_unique[k,:], :]), overwrite_a=True, check_finite=False); # Precalculate
                    for j in range(0, pt_interp.shape[1]):
                        pt_interp[pt_batchr, j] = rbf_interp_apply( pt_XYZ[pt_batchr, :], data_XYZ[batchy_indxz_unique[k,:], :], \
                            lu_solve(lu_factor_precalc, -data_orig[batchy_indxz_unique[k,:], j], trans=0, overwrite_b=True, check_finite=False)); #most powerful function call
                    # END FOR j
                # END FOR k
            # END IF
        # END IF
    # END IF
    
    if( FLG_ndimRevert == True ):
        pt_interp = pt_interp.ravel(); # Drop to a vector as the input was a vector but the loop supports 2D so it was 2D's preemptively
    # END IF
    
    return pt_interp
# END DEF

def rbf_interp_apply(pt_XYZ, data_XYZ, weights):
    if( data_XYZ is None ):
        # pt_XYZ is overloaded and is actually the required distance matrix already
        return pt_XYZ @ -weights
    else:
        #-- Euclidean distance --
        distances = np.zeros( (pt_XYZ.shape[0], data_XYZ.shape[0]), dtype=data_XYZ.dtype) # Preallocate to hold euclidean distance
        for i in range(0, data_XYZ.shape[1]): # Roll through each dim
            distances += np.subtract.outer(pt_XYZ[:, i], data_XYZ[:, i])**2 # Build euclidean distance dim-by-dim
        # END FOR i
        return np.sqrt(distances) @ -weights
    # END IF
#END DEF

#essential for weighting calcs
@jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def distance_calc_identical(data_XYZ):
    dist_identical = np.zeros((data_XYZ.shape[0],data_XYZ.shape[0]),dtype=data_XYZ.dtype);
    for jj in range(0,data_XYZ.shape[0]-1):
        #this reduces the effort for identical vector inputs by 1/2 or something like that
        dist_identical[jj+1:,jj] = np.sqrt( np.sum( (data_XYZ[jj, :] - data_XYZ[jj+1:, :])**2, axis=-1) ); # Euclidean distance calc that supports N dimensions
        dist_identical[jj,jj+1:] = dist_identical[jj+1:,jj]; #mirror it, faster than adding the transpose
    #END FOR jj
    return dist_identical
#END DEF

def distance_calc_identical_noNumba(data_XYZ):    
    # This is faster than broadcasting, faster than numba by a hair
    dist_identical = np.zeros( (data_XYZ.shape[0], data_XYZ.shape[0]), dtype=data_XYZ.dtype) # Preallocate to hold euclidean distance
    for i in range(0, data_XYZ.shape[1]): # Roll through each dim
        dist_identical += np.subtract.outer(data_XYZ[:, i], data_XYZ[:, i])**2 # Build euclidean distance dim-by-dim
    dist_identical = np.sqrt(dist_identical); # Sqrt to get the distance
    
    # # Numba-free option that doesn't redo work, but is somehow slowest (prob slicing stuff?):
    # ar1, ar2 = np.triu_indices(data_XYZ.shape[0], k=1);
    # dist_identical = np.zeros( (data_XYZ.shape[0], data_XYZ.shape[0]), dtype=data_XYZ.dtype); # Prep
    # dist_identical[ar1, ar2] = np.sqrt(np.sum((data_XYZ[ar1] - data_XYZ[ar2])**2, axis=-1));
    # al1, al2 = np.tril_indices(data_XYZ.shape[0], k=-1);
    # dist_identical.T[ar1, ar2] = dist_identical[ar1, ar2];
    
    return dist_identical
# END DEF

# def rbf_interp_weights(data_lat, data_long, data_orig):
#     # Remember the -data_orig represents the kernel which is -K, where K is 1 here
#     return lu_solve(lu_factor(distance_calc_identicalSpherical(data_lat, data_long), overwrite_a=True, check_finite=False), -data_orig, trans=0, overwrite_b=True, check_finite=False); #most powerful function call
# #END DEF


# ========= Inverse Distance Weighting =========
# This is faster than RBF but worse in quality
def idw_interp(pt_XYZ, data_XYZ, data_orig, power=2, nearPts=None):
    """
    Interpolates a value at a target point using Inverse Distance Weighting (Shepard's method).

    Args:
        data_XYZ: A numpy array of known point coordinates (dim x n).
        data_orig: A numpy vector of values at known points (n).
        pt_XYZ: A numpy array of the target point coordinates (dim x d).
        power: The power parameter controlling the distance weighting (default: 2).
        nearPts: The number of nearby points to use for interpolation for each target point coordinate set (default: None). Usually slower, can save memory.

    Returns:
        The interpolated value at the target points.
    """
    distances = np.zeros( (data_XYZ.shape[1], pt_XYZ.shape[1]), dtype=pt_XYZ.dtype) # Preallocate to hold euclidean distance
    for i in range(0, data_XYZ.shape[0]): # Roll through each dim
        distances += np.subtract.outer(data_XYZ[i,:], pt_XYZ[i,:])**2 # Build euclidean distance dim-by-dim
    # END FOR i
        
    if( (nearPts is None) or (nearPts is False) or (data_XYZ.shape[1] <= nearPts) ):
        
        if( power == 2 ):
            weights = 1.0 / ( distances + 1e-9); # Adding a small constant to avoid division by zero
        else:
            weights = 1.0 / ( distances ** (power/2) + 1e-9); # Adding a small constant to avoid division by zero
        # END IF
        
        # Normalize weights
        weights /= np.sum(weights, axis=0);
        
        pt_interp = weights.T @ data_orig; # np.dot also works here, but @ for matricies feels better
    else:
        pt_interp = np.empty( pt_XYZ.shape[1], dtype=data_orig.dtype); # Preallocate
        relevant_indxz = np.zeros( data_XYZ.shape[1], dtype=np.bool_); # Logical array for the relevant indexz to work with
        for k in range(0, pt_XYZ.shape[1]):
            # relevant_dist = distances[:, k]; # Get the distance between this location and the pts in the grid/pile of pts
            relevant_indxz[np.sort(np.argpartition(distances[:, k], nearPts)[:nearPts])] = True; # Update the nearPts# of relevant indexes
            
            if( power == 2 ):
                weights = 1.0 / ( distances[relevant_indxz, k] + 1e-9); # Adding a small constant to avoid division by zero
            else:
                weights = 1.0 / ( distances[relevant_indxz, k] ** (power/2) + 1e-9); # Adding a small constant to avoid division by zero
            # END IF
            
            # Normalize weights
            weights /= np.sum(weights, axis=0);
            
            pt_interp[k] = np.dot(weights.T, data_orig[relevant_indxz]).item();
            
            relevant_indxz &= False; # Reset
        # END FOR k
    # END IF
    
    return pt_interp
# END DEF