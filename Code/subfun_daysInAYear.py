#GOAL: Number of days in a Year
#RD on 8/23/2018
#
#INPUT: year
#OUTPUT: number days in the year
import numpy as np #import in here I dunno
#from numba import jit, prange

#@jit(parallel=True, fastmath=True) #try to activate some parallel speed up stuff (doesn't work really) #,nogil=True, parallel=True, fastmath=True
def subfun_daysInAYear(year):
    if( np.isscalar(year) ):
        #-----Leap Year Detection-----
        if np.mod(year,4) == 0: #leap year
            if (np.mod(year,100) == 0) and (np.mod(year,400) != 0):
                #-----Leap Year Skipped Detected - next will be 2100-----
                numDays = 365; #number of days
            else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                numDays = 366; #number of days
            #END IF
        else: #no leap year if this
            #-----No leap year detected-----
            numDays = 365; #number of days
        #END IF
    else:
        yearShape = year.shape; #get the shape
        year = year.flatten(); #floop it
        numDays = 365*np.ones(year.size,dtype=year.dtype); #prep, default to non-leap years
        numDays[(np.mod(year,4) == 0) & ( (np.mod(year,100) != 0) | (np.mod(year,400) == 0) )] = 366; #upgrade these to leap years
        if( yearShape != year.size ):
            year = year.reshape(yearShape); #reshape lazy
        #END IF
    #END IF

    return numDays #return success
#END DEF