#GOAL: Day Number to Date
#RD on 8/23/2018
#
#INPUT: dateRange_dayNum as [YR,dayNum] format for as many as you like - must be [#dates , 2] sized
#OUTPUT: date in [YRstart/M/D , YRend/M/D] format, for each date input

def subfun_dayNum_to_date(dateRange):
    import numpy as np #import in here I dunno
#dateRange = np.array([[2014,97],[2014,98]],dtype="int16"); #for debugging

#==============Catch input format issues==============
    if( isinstance(dateRange,tuple) | isinstance(dateRange,list) ):
        dateRange = np.array(dateRange); #convert to numpy array
    #END IF
    if( dateRange.ndim == 1 ):
        dateRange = dateRange[..., np.newaxis]; #add on a new axis b/c code demands it
    #END IF
    
    if (len(dateRange[0,:]) != 2) and (len(dateRange[:,0]) == 2):
        #Catch where a [3,arbitrary] was sent but this needs to operate on a [arbitrary,3] sized array
        # print("\n==============~Warning~==============");
        # print("Input size was [{},{}] and it is assumed that the dimensions are flipped -> adjusting to [{},{}] as [arbitrary,2] format is needed\n".format( len(dateRange[0,:]),len(dateRange[:,0]),len(dateRange[:,0]),len(dateRange[0,:]) ) );
        # dateRange = np.reshape(dateRange,(len(dateRange[0,:]),len(dateRange[:,0])));
        dateRange = np.swapaxes(dateRange,0,1); #flip axes
        
    elif (len(dateRange[0,:]) != 2) and (len(dateRange[:,0]) != 2):
        #Catch where something formatted completely wrong was sent
        print("\n==============ERROR==============");
        print("Input size was [{},{}] and it needs to be [arbitrary,2] - exiting dayNum to date fun!\n".format( len(dateRange[0,:]),len(dateRange[:,0]) ) );
        # return "No"; #can I? I will
        import sys
        sys.crash(); #even better
    #END IF
    
#==============Prep constants and preallocate==============
    monthDays = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]); #preps number of days in a month
    monthDays_Leap = np.array([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]); #preps number of days in a month
    
    
    dateRange_converted = np.zeros( (len(dateRange[:,0]),3),dtype="int16"); #preallocate (named so to make sharing code easier)

#==============Convert dates given in M/D/YR format to dayNum/YR format==============
    for i in range(0, len(dateRange[:,0]) ):
        #-----Leap Year Detection-----
        if np.mod(dateRange[i,0],4) == 0: #leap year
            #-----Leap Year Skipped Detected - next will be 2100-----
            if (np.mod(dateRange[i,0],100) == 0) and (np.mod(dateRange[i,0],400) != 0):
                #Same alg as no leap year used here
                
                tempDayCntr = dateRange[i,1]; #get the current day number
                tempMonCntr = 0; #counts the months
                while (tempDayCntr > 0) and (tempMonCntr <= 12): #makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
                    tempMonCntr = tempMonCntr + 1; #increment month used
                    tempDayCntr = tempDayCntr - monthDays[tempMonCntr-1]; #days, remove days in month from year
                #END WHILE
                if tempDayCntr > 0: #check for errors
                    print("Day count greater than 0 and days left in year: {} days left. Please check that input of {} year and {} day number.".format(tempDayCntr,dateRange[i,1],dateRange[i,0]));
                    return "No."; #can I? I will
                #END IF
                tempDayCntr = tempDayCntr + monthDays[tempMonCntr-1]; #days, calculate the day in the month we found
            
            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
            else:            
                tempDayCntr = dateRange[i,1]; #get the current day number
                tempMonCntr = 0; #counts the months
                while tempDayCntr > 0 and tempMonCntr <= 12: #makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
                    tempMonCntr = tempMonCntr + 1; #increment month used
                    tempDayCntr = tempDayCntr - monthDays_Leap[tempMonCntr-1]; #days, remove days in month from year
                #END WHILE
                if tempDayCntr > 0: #check for errors
                    print("Day count greater than 0 and days left in year: {} days left. Please check that input of {} year and {} day number.".format(tempDayCntr,dateRange[i,1],dateRange[i,0]));
                    return "No."; #can I? I will
                #END IF
                tempDayCntr = tempDayCntr + monthDays_Leap[tempMonCntr-1]; #days, calculate the day in the month we found
            
            #END IF
        #-----No leap year detected-----
        else: #no leap year if this
            
            tempDayCntr = dateRange[i,1]; #get the current day number
            tempMonCntr = 0; #counts the months
            while (tempDayCntr > 0) and (tempMonCntr <= 12): #makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
                tempMonCntr = tempMonCntr + 1; #increment month used
                tempDayCntr = tempDayCntr - monthDays[tempMonCntr-1]; #days, remove days in month from year
            #END WHILE
            if tempDayCntr > 0: #check for errors
                print("Day count greater than 0 and days left in year: {} days left. Please check that input of {} year and {} day number.".format(tempDayCntr,dateRange[i,1],dateRange[i,0]));
                return "No."; #can I? I will
            #END IF
            tempDayCntr = tempDayCntr + monthDays[tempMonCntr-1]; #days, calculate the day in the month we found
        
        #END IF
        
        dateRange_converted[i,:] = [dateRange[i,0],tempMonCntr,tempDayCntr]; #pack up the date and ship out
    #END FOR  
    
    return(dateRange_converted); #return success
#END DEF