#GOAL: Run any angle averaging strips on TEC data in requested area
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: gif_Grid

import numpy as np #import in here I dunno
from numba import jit, prange

# @jit(nopython=True,nogil=True,parallel=True,fastmath=True) #nopython=True, nogil=True, parallel=True, cache=True , nogil=True, parallel=True ,fastmath=True
def GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,gif_Grid_Lat,gif_Grid_Long,gif_Grid_Lat_Spaces,gif_Grid_Long_Spaces, \
                                        gif_Grid_Lat_Delta,gif_Grid_Long_Delta,dataRejectOrig,dataRejectLimitOrig,dataRejectMax):

    gif_Grid = np.nan*np.ones( (gif_Grid_Long_Spaces+1,gif_Grid_Lat_Spaces+1) ); #preallocate vTEC grid
    
    for j in prange(0,gif_Grid_Lat_Spaces+1): #gets the vTEC in the right area
        for l in range(0,gif_Grid_Long_Spaces+1):
#             jk = (pplat_portion <= (gif_Grid_Lat[j] + gif_Grid_Lat_Delta/2)) & \
#                 (pplat_portion >= (gif_Grid_Lat[j] - gif_Grid_Lat_Delta/2)) & \
#                 (pplong_portion <= (gif_Grid_Long[l] + gif_Grid_Long_Delta/2)) & \
#                 (pplong_portion >= (gif_Grid_Long[l]) - gif_Grid_Long_Delta/2));
            #combines to get just the correct indexes needed within
            #half a resolution up and down from the point
            jk = (pplat_portion < (gif_Grid_Lat[j] + gif_Grid_Lat_Delta)) & \
                (pplat_portion >= (gif_Grid_Lat[j])) & \
                (pplong_portion < (gif_Grid_Long[l] + gif_Grid_Long_Delta)) & \
                (pplong_portion >= (gif_Grid_Long[l]));
            #tedius checks have shown pcolor's boxes follow above math (verified in python too)
            if( np.sum(jk) != 0 ):
                temp = vTEC_portion[jk]; #choose just the data desired
                tempMean = np.mean(temp); #get the mean
                tempVar = np.var(temp); #get variance

                #loop to prevent too much data being nixed
                dataRejectLimit = dataRejectOrig;
                dataReject = dataRejectLimitOrig;
                while( (temp.size - np.sum(((tempMean+dataReject*tempVar) > temp) & ((tempMean-dataReject*tempVar) < temp) ) )*100/temp.size > dataRejectLimit ):
                    dataReject = dataReject*1.05; #increase the data rejection ratio to include more data
                    if( dataReject > dataRejectMax ):
                        dataRejectLimit = 100; #if this goes on for a bit then it is halted with this manuever
                    #END IF
                #END WHILE

                if( dataReject > dataRejectMax ): #if this occured leave the data, it is too sparse or too varied to deal with effectively
                    gif_Grid[l,j] = tempMean; #TECU, average data at the long/lat desired at the desired time
                else:
                    temp = temp[ ((tempMean+dataReject*tempVar) > temp) & ((tempMean-dataReject*tempVar) < temp) ]; #delete some extraneous data
                    #this is normal operation
                    gif_Grid[l,j] = np.mean(temp); #TECU, average data at the long/lat desired at the desired time
                #END IF
            #END IF
        #END FOR l
    #END FOR j
    
    return gif_Grid