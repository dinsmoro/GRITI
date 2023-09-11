# Spherical, linear (it's that - in places) RBF interpolation
#distance is spherical calc'd
#RBF kernel is linear (-eps*r) with eps = 1
#rbf_weights 1st to get weights, then feed it into rbf_interp
#rbf_interp's data_lat/data_long should be the same as weight's data_lat/data_long
#data_lat/data_long/data_orig should all be vectors of lat/long/data combos

#spherical units must be radians for lat and long
#spherical distance is approximated via the GR8-haversine function
#the GR8-haversine function is the haversine function with everything but the sqrt thrown out
#it gets the job done relatively - which is all that matters

import numpy as np
from numba import jit
# import scipy.linalg as scilin
# from scipy.linalg import lstsq as sp_lstsq
# from scipy.linalg import inv as sci_inv
# from scipy.linalg import cho_solve
# from scipy.linalg import solve_triangular
# from scipy.sparse.linalg import spsolve
# from scipy.linalg.lapack import dgesv
# from scipy.linalg import lu_factor, lu_solve

#faster in function that main code, prob b/c of slicing crap
# @jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def rbf_interp(pt_lat, pt_long, data_lat, data_long, weights):
    # data_out = np.empty(pt_lat.size,dtype=pt_lat.dtype);
    # for j in range(0,pt_lat.size):
    #     # data_out[j] = rbf_interp_calc(pt_lat[j], pt_long[j], data_lat, data_long, weights);
    #     data_out[j] = -(distance_calc(pt_lat[j], pt_long[j], data_lat, data_long)@weights).item();
    # #END FOR j
    # return data_out
    
    # return distance_calc(pt_lat, pt_long, data_lat, data_long)@-weights;
    
    # return distance_calc(pt_lat.reshape(-1,1), pt_long.reshape(-1,1), data_lat.reshape(1,-1), data_long.reshape(1,-1))@-weights;
    # return -np.matmul(distance_calc(pt_lat.reshape(-1,1), pt_long.reshape(-1,1), data_lat.reshape(1,-1), data_long.reshape(1,-1)),weights);
    
    #baseline - no matrix math (to confirm making the distance matrix is the thing that takes forever)
    # return np.sqrt(np.square(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)) + np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1)) * np.square(np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)))
    # return np.square(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)) + np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1)) * np.square(np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2))

    #one call return very reduced haversine formula matrix multiplied by the -weights (- is really the RBF kernel, which should wrap the np.sqrt() bit if it was more complex than -)
    # return np.sqrt(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)**2 + np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1)) * np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)**2)@-weights
    # return np.dot(np.sqrt(np.square(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)) + np.square(np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)) * np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1))), -weights)
    #------this is the fastest one there is------
    # return np.sqrt(np.square(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)) + np.square(np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)) * np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1))) @ -weights

    # return np.matmul(np.sqrt(np.square(np.sin((data_lat[None,:]-pt_lat[:,None])/2)) + np.square(np.sin((data_long[None,:]-pt_long[:,None])/2)) * np.cos(pt_lat[:,None]) * np.cos(data_lat[None,:])) , -weights)

    #-- alpha max beta min alg -- (slower)
    # lats_sq = np.abs(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)); #for alpha max beta min alg
    # longs_sq = np.abs(np.sqrt(np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1)))*np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)); #for the alphamax beta min
    # return (0.96043387010342*np.maximum.reduce((lats_sq,longs_sq)) + 0.397824734759316*np.minimum.reduce((lats_sq,longs_sq))) @ -weights; #alpha max beta min alg in action

    #-- great circle -- (fastest and most accurate)
    return np.arccos( np.sin(pt_lat.reshape(-1,1))*np.sin(data_lat.reshape(1,-1)) + np.cos(pt_lat.reshape(-1,1))*np.cos(data_lat.reshape(1,-1))*np.cos(data_long.reshape(1,-1)-pt_long.reshape(-1,1)) ) @ -weights
    
    #-gr8 circle approx via cos approx, didn't work out-
    # gr8circ = np.sin(pt_lat.reshape(-1,1))*np.sin(data_lat.reshape(1,-1)) + np.cos(pt_lat.reshape(-1,1))*np.cos(data_lat.reshape(1,-1))*np.cos(data_long.reshape(1,-1)-pt_long.reshape(-1,1)); #mid step
    # return ((-0.69813170079773212 * gr8circ * gr8circ - 0.87266462599716477) * gr8circ + 1.5707963267948966) @ -weights #approximation of arccos

    #-attempt to avoid powers, about same speed-
    # sinCalc_lat = np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2);
    # sinCalc_long = np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2);
    # return np.sqrt(sinCalc_lat*sinCalc_lat + np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1)) * sinCalc_long*sinCalc_long)@-weights
    
    #-attempt to avoid sine on 6000x6000, the needed math was not worth whatever dark magicks np.sin does-
    # return np.sqrt(np.square(sin_ap((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)) + np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1)) * np.square(sin_ap((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)))@-weights
        
    #---attempt to avoid sqrt which wasn't slower than a shoddy polynomial approximation---
    # midStep = np.square(np.sin((data_lat.reshape(1,-1)-pt_lat.reshape(-1,1))/2)) + np.square(np.sin((data_long.reshape(1,-1)-pt_long.reshape(-1,1))/2)) * np.cos(pt_lat.reshape(-1,1)) * np.cos(data_lat.reshape(1,-1))
    # return (3.4531713*midStep*midStep*midStep - 4.895681*midStep*midStep + 2.8887049*midStep + 0.0454234) @ -weights #not accurate enough at all it seems

    #---alternate form only using cosines and no squares (also slower)
    # pt_lat = pt_lat.reshape(-1,1);
    # pt_long = pt_long.reshape(-1,1);
    # data_lat = data_lat.reshape(1,-1);
    # data_long = data_long.reshape(1,-1);
    # x_p_y = data_lat + pt_lat;
    # x_m_y = data_lat - pt_lat;
    # a_m_b = data_long - pt_long;
    # return np.sqrt((-np.cos(a_m_b - x_p_y) - np.cos(a_m_b + x_m_y) - np.cos(a_m_b - x_m_y) - np.cos(a_m_b + x_p_y) - 2*np.cos(x_m_y) + 2*np.cos(x_p_y) + 4)/8) @ -weights
#END DEF


# def sin_ap(deg): #only good 180 to -180 (don't even think about 270)
#     return np.sign(deg)*(720*np.abs(deg)-4*deg*deg)/(40500-180*np.abs(deg)+deg*deg) #https://en.wikipedia.org/wiki/Bhaskara_I%27s_sine_approximation_formula
# #END DEF

# def cos_ap(deg): #only good 180 to -180
#     reg = np.mod(90-deg,180); #for cosine gotta roll - approx is only good for 180 to -180 and this shift makes less than 90 not work w/o mod
#     cosish = (720*np.abs(reg)-4*reg*reg)/(40500-180*np.abs(reg)+reg*reg);
#     cosish[np.abs(deg)>90] *= -1; #negate needed stuff
#     return cosish #https://en.wikipedia.org/wiki/Bhaskara_I%27s_sine_approximation_formula
# #END DEF

# #used directly by rbf_interp which is the only one that uses it so commented out (also in a numba implementation that was slower)
# @jit(nopython=True,nogil=True,parallel=False,fastmath=True)
# def distance_calc(lat1, long1, lat2, long2):
#     #haversine
#     # dist_sphere = np.sin((lat2-lat1)/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin((long2-long1)/2)**2; #modified haversine formula
#     # return np.arctan2(np.sqrt(dist_sphere), np.sqrt(1 - dist_sphere)) #modified haversine formula
#     # return np.arctan(np.sqrt(dist_sphere)/np.sqrt(1 - dist_sphere)) #modified haversine formula

#     # return z(45∘−(z−1)(14∘+3.83∘z))
#     # np.pi/4*x - x*(fabs(x) - 1)*(0.2447 + 0.0663*fabs(x)
    
#     # return np.sqrt(np.sin((lat2-lat1)/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin((long2-long1)/2)**2); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
    
#     # lat1_g, lat2_g = np.meshgrid(lat1,lat2);
#     # long1_g, long2_g = np.meshgrid(long1,long2);
#     # return np.sqrt(np.sin((lat2_g-lat1_g)/2)**2 + np.cos(lat1_g) * np.cos(lat2_g) * np.sin((long2_g-long1_g)/2)**2); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
    
#     dist_sphere = np.empty((lat1.size,lat2.size),dtype=lat1.dtype);
#     for i in prange(0,lat1.size):
#         for j in range(0,lat2.size):
#             dist_sphere[i,j] = np.sqrt(np.sin((lat2[j]-lat1[i])/2)**2 + np.sin((long2[j]-long1[i])/2)**2 * np.cos(lat1[i]) * np.cos(lat2[j])); #faster to do this in parallel manually it seems
#         #END FOR j
#     #END FOR i
#     return dist_sphere
# #END DEF

#essential for weighting calcs
@jit(nopython=True,nogil=True,parallel=False,fastmath=True)
def distance_calc_identical(lat1, long1):
    dist_identical = np.zeros((lat1.size,lat1.size),dtype=lat1.dtype);
    for jj in range(0,lat1.size-1):
        #this reduces the effort for identical vector inputs by 1/2 or something like that
        # dist_sphere = np.sin((lat1[jj]-lat1[jj+1:])/2)**2 + np.cos(lat1[jj]) * np.cos(lat1[jj+1:]) * np.sin((long1[jj]-long1[jj+1:])/2)**2; #modified haversine formula
        # dist_identical[jj+1:,jj] = np.arctan2(np.sqrt(dist_sphere), np.sqrt(1 - dist_sphere)) #modified haversine formula
        # dist_identical[jj+1:,jj] = np.arctan(np.sqrt(dist_sphere)/np.sqrt(1 - dist_sphere)) #modified haversine formula
        
        # lats = np.sin((lat1[jj]-lat1[jj+1:])/2);
        # longs = np.sin((long1[jj]-long1[jj+1:])/2);  
        # dist_identical[jj+1:,jj] = np.sqrt(lats*lats + np.cos(lat1[jj]) * np.cos(lat1[jj+1:]) * longs * longs); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        # lats_sq = np.abs(np.sin((lat1[jj]-lat1[jj+1:])/2)); #for alpha max beta min alg
        # longs_sq = np.abs(np.sqrt(np.cos(lat1[jj]) * np.cos(lat1[jj+1:]))*np.sin((long1[jj]-long1[jj+1:])/2)); #for the alphamax beta min
        # dist_identical[jj+1:,jj] = 0.96043387010342*np.max((lats_sq,longs_sq)) + 0.397824734759316*np.min((lats_sq,longs_sq)); #alpha max beta min alg in action
        # dist_identical[jj+1:,jj] = np.sqrt(np.sin((lat1[jj]-lat1[jj+1:])/2)**2 + np.cos(lat1[jj]) * np.cos(lat1[jj+1:]) * np.sin((long1[jj]-long1[jj+1:])/2)**2); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        # dist_identical[jj+1:,jj] = np.sqrt(np.square(np.sin((lat1[jj]-lat1[jj+1:])/2)) + np.cos(lat1[jj]) * np.cos(lat1[jj+1:]) * np.square(np.sin((long1[jj]-long1[jj+1:])/2))); #reduced haversine formula (seems sqrt is the most important part, divined via squinting)
        dist_identical[jj+1:,jj] = np.arccos( np.sin(lat1[jj])*np.sin(lat1[jj+1:]) + np.cos(lat1[jj])*np.cos(lat1[jj+1:])*np.cos(long1[jj]-long1[jj+1:]) ); #great circle
        dist_identical[jj,jj+1:] = dist_identical[jj+1:,jj]; #mirror it, faster than adding the transpose
    #END FOR jj
    # dist_identical += dist_identical.T; #add on the other half (very broken in jit mode)
    # return dist_identical + dist_identical.T
    return dist_identical
#END DEF

#since just -radius much faster to just - it I guess (yolo SPEED, replace the -'s around with rbf_kernel to use another kernel future me)
# @jit(nopython=True,nogil=True,parallel=False,fastmath=True)
# def rbf_kernel(radius):
#     #if want to use epislon "shape param", it's really -epislon*radius - but I use 1 and neglect for speed
#     return -radius
# #END DEF

#called in code directly now since it is one line, can be uncommented and used ez pz tho (very minor speed up being in code)
# # @jit(nopython=True,nogil=True,parallel=False,fastmath=True)
# def rbf_interp_weights(data_lat, data_long, data_orig):
#     return -np.linalg.solve(distance_calc_identical(data_lat, data_long), data_orig) #best so far
#     # return np.linalg.solve(distance_calc(data_lat.reshape(-1,1), data_long.reshape(-1,1), data_lat.reshape(1,-1), data_long.reshape(1,-1)), -data_orig) #best so far
#     # return np.matmul(sci_inv(rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T))),data_orig)
#     # return scilin.solve(rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T)), data_orig)
#     # return np.linalg.inv(rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T)))@data_orig #big slow
#     # return sp_lstsq(rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T)), data_orig , lapack_driver='gelsy', check_finite=False)
#     # return np.linalg.lstsq(rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T)), data_orig)
    
#     #baseline
#     # return distance_calc(data_lat.reshape(-1,1), data_long.reshape(-1,1), data_lat.reshape(1,-1), data_long.reshape(1,-1)) #just call distance calc
#     # return distance_calc_identical(data_lat, data_long) #just call distance calc
    
#     # #cholesky decomp solve
#     # rbf_mid = rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T));
#     # L = np.linalg.cholesky(rbf_mid@rbf_mid.T + np.eye(data_orig.size))
#     # return cho_solve((L,True), data_orig)
#     # return solve_triangular(L, data_orig)
    
#     # return spsolve(rbf_kernel(distance_calc(data_lat, data_long, data_lat.T, data_long.T)), data_orig)
    
#     # _, _, coeffs, _ = dgesv(-distance_calc_identical(data_lat, data_long), data_orig.copy(), overwrite_a=True, overwrite_b=True)
#     # return coeffs #about as fast as np.linalg.solve
    

#     # return lu_solve(lu_factor(-distance_calc_identical(data_lat, data_long),overwrite_a=True, check_finite=False), data_orig) #about as fast as np.linalg.solve as well
# #END DEF

#merged into rbf_interp so actually useless
# @jit(nopython=True,nogil=True,parallel=False,fastmath=True)
# def rbf_interp_calc(pt_lat, pt_long, data_lat, data_long, weights):
#     # return np.dot(rbf_kernel(distance_calc(pt_lat, pt_long, data_lat, data_long)), weights).item()
#     # return np.matmul(rbf_kernel(distance_calc(pt_lat, pt_long, data_lat, data_long)), weights).item() #not supported by numba
#     # return rbf_kernel(distance_calc(pt_lat, pt_long, data_lat, data_long)).dot(weights).item()
#     return (-distance_calc(pt_lat, pt_long, data_lat, data_long)@weights).item() #all about the same, this seems a tad faster
#     #otions that may be faster: a.dot(b)
#     # np.matmul(a,b)
#     # a@b
# #END DEF


#see cupy.matmul

# #--- from https://stackoverflow.com/a/64198479/2403531 ---
# from numba import cuda, float32

# # Controls threads per block and shared memory usage.
# # The computation will be done on blocks of TPBxTPB elements.
# TPB = 16

# @cuda.jit
# def fast_matmul(A, B, C):
#     # Define an array in the shared memory
#     # The size and type of the arrays must be known at compile time
#     sA = cuda.shared.array(shape=(TPB, TPB), dtype=float32)
#     sB = cuda.shared.array(shape=(TPB, TPB), dtype=float32)

#     x, y = cuda.grid(2)

#     tx = cuda.threadIdx.x
#     ty = cuda.threadIdx.y
#     bpg = cuda.gridDim.x    # blocks per grid

#     # Each thread computes one element in the result matrix.
#     # The dot product is chunked into dot products of TPB-long vectors.
#     tmp = float32(0.)
#     for i in range(bpg):
#         # Preload data into shared memory
#         sA[ty, tx] = 0
#         sB[ty, tx] = 0
#         if y < A.shape[0] and (tx+i*TPB) < A.shape[1]:
#           sA[ty, tx] = A[y, tx + i * TPB]
#         if x < B.shape[1] and (ty+i*TPB) < B.shape[0]:
#           sB[ty, tx] = B[ty + i * TPB, x]

#         # Wait until all threads finish preloading
#         cuda.syncthreads()

#         # Computes partial product on the shared memory
#         for j in range(TPB):
#             tmp += sA[ty, j] * sB[j, tx]

#         # Wait until all threads finish computing
#         cuda.syncthreads()
#     if y < C.shape[0] and x < C.shape[1]:
#         C[y, x] = tmp
# #END DEF