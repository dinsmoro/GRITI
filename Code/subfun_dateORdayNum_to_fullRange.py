#GOAL: Date OR Day Number to the full range between the range
#RD on 8/23/2018
#
#INPUT: date range in [YRstart/M/D , YRend/M/D] format, [2 , 3] sized ~OR~ date range in [YRstart/dayNum , YRend/dayNum] format, [2,2] sized
#OUTPUT: dateRange_dayNum_full and dateRange_dayNum_full as [#datesInRange,3] and [#datesInRange,2], respectively, format
import numpy as np
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date

def subfun_dateORdayNum_to_fullRange(dateRange):

    #==============Catch input format issues==============
    if( isinstance(dateRange,tuple) | isinstance(dateRange,list) ):
        dateRange = np.array(dateRange); #convert to numpy array
    #END IF
    if( dateRange.ndim == 1 ):
        dateRange = dateRange[..., np.newaxis]; #add on a new axis b/c code demands it
    #END IF

    if (len(dateRange[0,:]) == 2) and (len(dateRange[:,0]) == 3):
        #Catch where a [3,2] was sent and this needs a [2,3] to work on. Assuming everything else is good in it!
        # print("\n==============~Warning~==============");
        # print("Input size was [{},{}] and it is assumed that the dimensions are flipped -> adjusting to [{},{}] as [2,3] (m/d/yr) format is required (other than [2,2] daynum/Yr format)\n".format( len(dateRange[:,0]),len(dateRange[0,:]),len(dateRange[0,:]),len(dateRange[:,0]) ) );
        # dateRange = np.reshape(dateRange,(len(dateRange[0,:]),len(dateRange[:,0])));
        dateRange = np.swapaxes(dateRange,0,1); #flip axes
        
    elif (len(dateRange[0,:]) != 3) and (len(dateRange[:,0]) != 2 and (len(dateRange[0,:]) != 2) ):
        #Catch where a completely wrong size was sent as a [2,3] or [2,2] is required
        print("\n==============ERROR==============");
        print("Input size was [{},{}] and it needs to be [2,3] or [2,2] - exiting date to dayNum fun!\n".format( len(dateRange[:,0]),len(dateRange[0,:]) ) );
        # return "No"; #good luck with this return
        import sys
        sys.crash(); #even better
    #END IF
    
    if (len(dateRange[:,0]) == 2) and (len(dateRange[0,:]) == 3):
        #Catch if a Yr/M/D was sent but this will use Yr/dayNum format because it's easier to math with them
        dateRange = subfun_date_to_dayNum(dateRange); #call it dateRange to make programming easy, but now it is really dateRange_dayNum
    #END IF

    #==============Convert date range to the full range in between, covering years==============
    dateRange_dayNum_full = dateRange[1,0] - dateRange[0,0]; #yr difference
    if dateRange_dayNum_full == 0:
        dateRange_dayNum_full = dateRange[1,1] - dateRange[0,1]; #day difference
        dateRange_dayNum_full = np.hstack( [np.tile(dateRange[0,0],(dateRange_dayNum_full+1,1)) , np.reshape(np.arange(dateRange[0,1],dateRange[1,1]+1,1,dtype="int16"), (-1,1))] ); #create full date range from min/max days since within same year
        #way cooler in matlab, but works - gets the days in between as an array - note that reshapeing a (-1,1) means that -1 automatically calcs size
    #     dateRange_yearRange = dateRange(1,0); #set the year range as one year
    else:
        dateRange_yearRange = np.arange(dateRange[0,0],dateRange[1,0]+1,1,dtype="int16") #get the full year range from min to max
        dateRange_dayNum_full = [None for i in range(0, len(dateRange_yearRange) )]; #preallocate
        for i in range(0, len(dateRange_yearRange) ):
            #Leap Year Detection
            if np.mod(dateRange_yearRange[i],4) == 0: #leap year
                #LEAP YEAR! Possibly.
                if (np.mod(dateRange_yearRange[i],100) == 0) and (np.mod(dateRange_yearRange[i],400) != 0):
                    #NO LEAP YEAR
                    #Leap Year Skipped Detected - next will be 2100
                    daysInYear = 365;
                else:
                    #Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
                    daysInYear = 364;
                #END IF
            else: #no leap year if this
                #no leap year
                daysInYear = 365;
            #END IF 
            
            if(i == 0 ): # start year - go to end
                dateRange_dayNum_full[i] = np.hstack( [np.tile(dateRange_yearRange[i],(daysInYear-dateRange[0,1]+1,1)) , np.reshape(np.arange(dateRange[0,1],daysInYear+1,1,dtype="int16"),(-1,1))] ); #start of date to the end of the year
            elif( i == (len(dateRange_yearRange)-1) ): # end year - go to date
                dateRange_dayNum_full[i] = np.hstack( [np.tile(dateRange_yearRange[i],(dateRange[1,1],1)) , np.reshape(np.arange(1,dateRange[1,1]+1,1,dtype="int16"),(-1,1))] ); #start of year to the end date
            else: # a full year covered, go start to end of that year
                dateRange_dayNum_full[i] = np.hstack( [np.tile(dateRange_yearRange[i],(daysInYear,1)) , np.reshape(np.arange(1,daysInYear+1,1,dtype="int16"),(-1,1))] ); #create full date range from min/max days since within same year
            #END IF
        # END FOR i loop per year
        dateRange_dayNum_full = np.vstack(dateRange_dayNum_full); #stack em together
    #END IF
        
    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #convert all the day number dates to Yr/M/D dates
    
    return(dateRange_full,dateRange_dayNum_full); #return both
#END DEF