"""
GOAL: Time Indices to Date Range (Dates Involved)
#RD on 6/23/2023
#
#INPUT: Times in seconds aligned to a zero hour date
#OUTPUT: dateRange_dayNum as [YR/dayNum] format for each date input
#options!: 
0 = Output as [YR - DayNum] [DEFAULT!]
1 = Output as [YR - M - D] 
2 = Output as just [dayNum]
"""
import numpy as np #import in here I dunno
# from numba import jit
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date


# @jit(nopython=True,nogil=True,parallel=False,cache=True,fastmath=True)
def subfun_time_to_dateRange(timeIndices, dates_zeroHr = None, FLG_timesWRTzeroHr = False, options = 0):
    
    #==============Catch input format issues==============
    # Stop here if it's a show stopper
    if( (FLG_timesWRTzeroHr == True) & np.all(dates_zeroHr == None) ):
        print('ERROR in subfun_time_to_dateRange: Without a zero hour reference but provided timeIndices being provided with respect to (wrt) a zero hour means this can\'t work. Sad. Returning None.');
        return None
    #END IF
    
    if( options not in [0, 1, 2] ):
        print('WARNING in subfun_time_to_dateRange: Options '+str(options)+' provided which is invalid. Using 0 as default.');
        options = 0; #revert
    #END IF
    
    # Align dates_zeroHr as needed
    if( np.all(dates_zeroHr == None) & (options != 2) ):
        print('WARNING in subfun_time_to_dateRange: Without a zero hour reference only returning just dayNum is possible, so options changed from '+str(options)+' to 2.');
        options = 2; #Without year provided can only output mode 2 (dayNum only)
    elif( np.all(dates_zeroHr != None) & (options != 2) ):
        if( np.isscalar(dates_zeroHr) == True ):
            if( FLG_timesWRTzeroHr == True ):
                print('ERROR in subfun_time_to_dateRange: dates_zeroHr provided is a scalar, which is not enough to rebuild the timeIndices. Sad. Returning None.');
                return None
            #END IF
            dates_zeroHr_year = dates_zeroHr; # Got the year
        else:
            if( dates_zeroHr.ndim != 1):
                print('WARNING in subfun_time_to_dateRange: dates_zeroHr has 2 dimensions which is bad. Guesstimating zero hour but may not succeed.')
                dates_zeroHr = dates_zeroHr[np.int16( np.floor((len(dates_zeroHr[:,0]) - 1)/2) ), :];
                print('Zero hour guesstimated to be '+str(dates_zeroHr)+' and pressing on.');
            #END IF
            dates_zeroHr_year = dates_zeroHr[0]; #1st is year usually
        #END IF    
    #END IF
    
    # Make sure zero hr is in yr/daynum format
    if( dates_zeroHr.size == 3 ):
        dates_zeroHr = subfun_date_to_dayNum(dates_zeroHr)[0]; #convert to yr/daynum format
    #END IF
    
    # Only operate on unique time indices
    timeIndices = np.unique(timeIndices); #only operate on unique indices
    
    # Fix timeIndices if they arrived in zero hr reference form
    if( (FLG_timesWRTzeroHr == True) & np.all(dates_zeroHr != None) ):
        timeIndices_zeroLoc = np.where( np.abs(timeIndices) == np.min(timeIndices) )[0].item(); #get index closest to 0 (ideally 0)
        timeIndices += dates_zeroHr[1]*86400; #conversion to time w/o year in it
    elif( (FLG_timesWRTzeroHr == False) & np.all(dates_zeroHr != None) ):
        timeIndices_zeroLoc = np.where( np.abs(timeIndices - dates_zeroHr[1]*86400) == np.min(np.abs(timeIndices - dates_zeroHr[1]*86400)) )[0]; #get index closest to dayNum time
        if( timeIndices_zeroLoc.size == 1 ):
            timeIndices_zeroLoc = timeIndices_zeroLoc.item(); #gottem
        else:
            print('ERROR in subfun_time_to_dateRange: Time range too ambiguous (there are two matches for the zero hr day num) so the real zero hr loc cannot be determined. Time to add year to the time units. Sad. Returning None.');
            return None
        #END IF
    #END IF
    
    #==============Convert times to date range==============
    dates_dayNum = np.unique(timeIndices//86400).astype(np.int16); #Unique day nums
    dates_dayNum_diff = np.diff(dates_dayNum); #get the diff
    if( (options != 2) & np.any(dates_dayNum_diff > 1) ):
        # Deal with the year changing
        dates_yearArray = np.empty(dates_dayNum.size,dtype=np.int16); #prep it
        splitz = np.append(np.insert(np.where(dates_dayNum_diff > 1)[0]+1,0,0),dates_dayNum.size); #get where splitz occur
        dates_yearLoc = np.where(timeIndices[timeIndices_zeroLoc]//86400 == dates_dayNum)[0].item(); #get location of reference year
        dates_yearOffset = np.where(dates_yearLoc > splitz)[0][0]; #get the offset needed to align the splitz to the reference year 
        for j in range(0,splitz.size-1):
            dates_yearArray[splitz[j]:splitz[j+1]] = dates_zeroHr_year+j-dates_yearOffset; #split it up
        #END FOR j
        dates_dayNum = np.vstack( (dates_yearArray, dates_dayNum) ).T; #stack it up just right
    elif( (options != 2) & np.all(dates_dayNum_diff == 1) ):
        #year constant
        dates_dayNum = np.vstack( (np.ones(dates_dayNum.size,dtype=np.int16)*dates_zeroHr_year, dates_dayNum) ).T; #stack it up just right
    #END IF
    
    
    #Return as needed
    if( options == 0 ):
            return dates_dayNum
    elif( options == 1 ):
            return subfun_dayNum_to_date(dates_dayNum)
    elif( options == 2 ):
            return dates_dayNum
    #END IF
#END DEF