"""
Accelerates the point averaging with NUMBA
"""
import numpy as np
# from numba import jit, prange

# Numba crashes the python kernel when this is run so it's ditched for now !!make range->prange if bringing Numba back!!
# @jit(nopython=True,nogil=True,parallel=True,fastmath=False) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_TEC_avgPtACCELERATOR(TEC_timeUnique,TEC_timeUnique_size,tempTime,tempTEC,
    dataReject,dataRejectOrig,dataRejectLimit,dataRejectLimitOrig,dataRejectMax):
    
    avgPt_vTEC = np.zeros( (TEC_timeUnique_size) ); #preallocate
    for t in range(0,TEC_timeUnique_size): #prange if JIT time!
        
        j = tempTime == TEC_timeUnique[t]; #find the time that connect to the time we are at
        
        if( np.any(j) ):
            tempTempTEC = tempTEC[j]; #get the delta-vTEC involved, do some filtering
            
            # if( tempTempTEC.size > 1 ): #jit parallel/prange get mad at if statement AND silently NaN divide-by-zero errors (non-parallel jit needs if statements to avoid divide-by-zero, but it is slower than parallel)
            
            if(tempTempTEC.size-1 > 0): #only do variance stuff if it's more than 1 number b/c how do you calc the variance for 1
                tempMean = np.mean( tempTempTEC ); #get the mean
                # tempVar = np.var( tempTempTEC, ddof=1 ); #get variance, ddof=1 to reduce bias (matlab does automatically)
                tempVar = np.sum((tempTempTEC - np.mean(tempTempTEC))**2)/(tempTempTEC.size-1); #get de-biased variance (NUMBA not built-in)
            
                #loop to prevent too much data being nixed
                while( (tempTempTEC.size - np.sum(((tempMean+dataReject*tempVar) > tempTempTEC) & ((tempMean-dataReject*tempVar) < tempTempTEC) ) )*100/tempTempTEC.size > dataRejectLimit ):
                    dataReject = dataReject*1.05; #increase the data rejection ratio to include more data
                    if( dataReject > dataRejectMax ):
                        dataRejectLimit = 100; #if this goes on for a bit then it is halted with this manuever
                    #END IF
        #             (length(tempTempTEC) - sum(((tempMean+dataReject*tempVar) > tempTempTEC) & ((tempMean-dataReject*tempVar) < tempTempTEC) ) )*100/length(tempTempTEC)
        #             dataReject
        #             dataRejectLimit
                #END WHILE
                
                if( dataReject <= dataRejectMax ): #if this DID NOT occur leave the data, it is too sparse or too varied to deal with effectively
                    tempTempTEC = tempTempTEC[ ((tempMean+dataReject*tempVar) > tempTempTEC) & ((tempMean-dataReject*tempVar) < tempTempTEC) ]; #delete some extraneous data
                    #this is normal operation (data isn't too sparse)
                #END IF
                
                dataReject = dataRejectOrig; #reset dataReject
                dataRejectLimit = dataRejectLimitOrig; #reset dataRejectLimit
            #END IF
            
            # if( tempTempTEC.size != 0 ):
            avgPt_vTEC[t] = np.mean(tempTempTEC); #TECU, average of the delta-vTEC points in that radius at the same time into one value
            # elif(tempTempTEC.size == 1):
            #     avgPt_vTEC[t] = np.mean(tempTempTEC); #TECU, average of the delta-vTEC points in that radius at the same time into one value
            # else:
            #     avgPt_vTEC[t] = np.nan; #set to NaN b/c nothing is left
            # #END IF
        else:
            avgPt_vTEC[t] = np.nan; #return NaN if no data
        #END
        
    #END FOR t
    
    return avgPt_vTEC