#GOAL: Pad a date range with a day before and a day in front (or other things like only before or multiples, depending on the options)
#RD on 9/18/2018
#
#INPUT: dateRange in [YR/dayNum] format [2,arbitrary] shape
#OUTPUT: dateRange in [YR/dayNum] format [2,arbitrary+2] shape (may be +1 or +3 or whatever depending on options)
#options!: 
#padDirection
#0 = Both sides get padded [DEFAULT!] e.g. [2014,212] -> [2014, 211 ; 2014, 212 ; 2014, 213]
#1 = Forward gets padded e.g. [2014, 212] -> [2014, 212 ; 2014, 213]
#2 = Backwards gets padded e.g. [2014, 212] -> [2014, 211 ; 2014, 212]
#padNumber
#1 = e.g. [2014, 212] -> [2014, 211 ; 2014, 212 ; 2014, 213] [DEFAULT!]
#It'll tell you you're wronk if <1 occurs
#It'll tell you you're wronk if non-integers are passed

import numpy as np #import in here I dunno

def subfun_addADay(dateRange,padDirection = 0, padNumber = 1):
#dateRange = np.array([[2014,212]],dtype="int16"); #for debugging
#dateRange = np.array([[2014,212],[2014,213],[2014,214]],dtype="int16"); #for debugging
#padDirection = 3; #for debugging
#padNumber = 1; #for debugging
    
#==============Catch input format issues==============
    if(len(dateRange[0,:]) != 2) and (len(dateRange[:,0]) == 2):
        #Catch where a [2,arbitrary] was sent but this needs to operate on a [arbitrary,2] sized array
        print("\n==============~Warning~==============");
        #print("Input size was [{},{}] and it is assumed that the dimensions are flipped -> adjusting to [{},{}] as [arbitrary,2] format is needed.\n".format( len(dateRange[0,:]),len(dateRange[:,0]),len(dateRange[:,0]),len(dateRange[0,:]) ) ); #.format version
        print("Input size was ["+str(len(dateRange[0,:]))+","+str(len(dateRange[:,0]))+"] and it is assumed that the dimensions are flipped -> adjusting to ["+str(len(dateRange[:,0]))+","+str(len(dateRange[0,:]))+"] as [arbitrary,2] format is needed.\n"); #non-.format version for nopython calls
        dateRange = np.reshape(dateRange,(len(dateRange[0,:]),len(dateRange[:,0])));
        
    elif(len(dateRange[0,:]) != 2) and (len(dateRange[:,0]) != 2):
        #Catch where something formatted completely wrong was sent
        print("\n==============ERROR==============");
        #print("Input size was [{},{}] and it needs to be [arbitrary,2] - exiting date to dayNum fun!\n".format( len(dateRange[0,:]),len(dateRange[:,0]) ) ); #.format version
        print("Input size was ["+str(len(dateRange[0,:]))+","+str(len(dateRange[:,0]))+"}] and it needs to be [arbitrary,2] - exiting date to dayNum fun!\n"); #non-.format version for nopython calls
        #return "No"; #can I? I will
    #END IF
        
    if( isinstance(padNumber, int) == 0 ):
        print("\n==============~Warning~==============");
        print("padNumber was sent as a non-integer '"+str(padNumber)+"' so gonna try to convert it to an integer. Fix this ok thx.");
        padNumber = int(padNumber); #convert to int and hope for the best!
    #END IF
    if( padNumber < 1 ):
        print("\n==============~Warning~==============");
        print("padNumber is "+str(padNumber)+", which is less than 1 so just make it 1. Fix this ok thx.");
        padNumber = 1; #make pad number 1 because its gotta be
    #END IF 

#==============Pad Dates onto the date array==============    
    if( padDirection == 0 ):
    #-----Option 0 returns padding on both sides of the date range by the padNumber requested-----
        
        dateRange_padded = np.zeros([len(dateRange[:,0])+2*padNumber,2],dtype="int16"); #preallocate    
        dateRange_padded[ (padNumber) : (len(dateRange[:,0])+padNumber), :] = dateRange; #insert original data
        
        cntrBefore = padNumber; #run seperate counters for before and after
        FLG_after = 1; #flag for after that facilitates support for padding multiple years forward sucessfully (1 == off shut up)
        for i in range(0,padNumber):
            #==============Tackle before padding==============
            dateRange_padded[i,:] = np.array([dateRange[0,0] , dateRange[0,1] - cntrBefore]); #fill in before date
            cntrBefore = cntrBefore - 1; #subtract off of the counter
            
            while( dateRange_padded[i,1] < 1 ): #while supports padding multiple years if it ever came to that
                #NOTE: the year checks are run on the previous year since we're going back one!
                #-----Leap Year Detection-----
                if np.mod(dateRange_padded[i,0]-1,4) == 0: #leap year
                    #-----Leap Year Skipped Detected - next will be 2100-----
                    if (np.mod(dateRange_padded[i,0]-1,100) == 0) and (np.mod(dateRange_padded[i,0]-1,400) != 0):
                        #Same alg as no leap year used here
                        dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+365]); #adjust to the year before
                    else:
                    #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                        dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+366]); #adjust to the year before
                    #END IF
                        
                #-----No leap year detected-----
                else: #no leap year if this
                    dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+365]); #adjust to the year before
                #END IF
            #END WHILE
            
            #==============Tackle after padding==============
            cntrAfter = i+len(dateRange[:,1])+padNumber; #calc forward counter
            dateRange_padded[cntrAfter,:] = np.array([dateRange[-1,0], dateRange[-1,1] + i+1]); #fill in after date
            #NOW running a check on the current year to see if an increment is needed
            #-----Leap Year Detection-----
            if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                    #Same alg as no leap year used here
                    if( dateRange_padded[cntrAfter,1] > 365 ):
                        FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if( dateRange_padded[cntrAfter,1] > 366 ):
                        FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                    #END IF
                #END IF 
            #-----No leap year detected-----
            else: #no leap year if this
                if( dateRange_padded[cntrAfter,1] > 365 ):
                    FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                #END IF
            #END IF
            
            while( FLG_after == 0 ): #while supports padding multiple years if it ever came to that
                #-----Leap Year Detection-----
                if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                    #-----Leap Year Skipped Detected - next will be 2100-----
                    if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                        #Same alg as no leap year used here
                        dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1, dateRange_padded[cntrAfter,1]-365]); #adjust to the year after
                        #-----Rerun the year check on the next to see if we need to go through this again-----
                        if np.mod(dateRange_padded[i,0],4) == 0: #leap year
                            #-----Leap Year Skipped Detected - next will be 2100-----
                            if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                                #Same alg as no leap year used here
                                if( dateRange_padded[cntrAfter,1] < 366 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            else:
                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                if( dateRange_padded[cntrAfter,1] < 367 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            #END IF
                            
                        #-----No leap year detected-----
                        else: #no leap year if this
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    else:
                    #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                        dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1, dateRange_padded[cntrAfter,1]-366]); #adjust to the year after
                        #-----Rerun the year check on the next to see if we need to go through this again-----
                        if np.mod(dateRange_padded[i,0],4) == 0: #leap year
                            #-----Leap Year Skipped Detected - next will be 2100-----
                            if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                                #Same alg as no leap year used here
                                if( dateRange_padded[cntrAfter,1] < 366 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            else:
                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                if( dateRange_padded[cntrAfter,1] < 367 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            #END IF
                        #-----No leap year detected-----
                        else: #no leap year if this
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    #END IF
                #-----No leap year detected-----
                else: #no leap year if this
                    dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1], dateRange_padded[cntrAfter,1]-365); #adjust to the year after
                    #-----Rerun the year check on the next to see if we need to go through this again-----
                    if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                        #-----Leap Year Skipped Detected - next will be 2100-----
                        if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                            #Same alg as no leap year used here
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        else:
                        #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                            if( dateRange_padded[cntrAfter,1] < 367 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    #-----No leap year detected-----
                    else: #no leap year if this
                        if( dateRange_padded[cntrAfter,1] < 366 ):
                            FLG_after = 1; #we're free
                        #END IF
                    #END IF
                #END IF
            #END WHILE
        #END FOR
        
        return(dateRange_padded); #return success
    
    elif( padDirection == 1 ):
    #-----Option 1 returns forward padding onto the date range by the padNumber requested-----
        
        dateRange_padded = np.zeros([len(dateRange[:,0])+padNumber,2],dtype="int16"); #preallocate    
        dateRange_padded[ 0 : (len(dateRange[:,0])) ,:] = dateRange; #insert original data
        
        FLG_after = 1; #flag for after that facilitates support for padding multiple years forward sucessfully (1 == off shut up)
        for i in range(0,padNumber):        
            #==============Tackle after padding==============
            cntrAfter = i+len(dateRange[:,0]); #calc forward counter
            dateRange_padded[cntrAfter,:] = np.array([dateRange[-1,0], dateRange[-1,1] + i+1]); #fill in after date
            #NOW running a check on the current year to see if an increment is needed
            #-----Leap Year Detection-----
            if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                    #Same alg as no leap year used here
                    if( dateRange_padded[cntrAfter,1] > 365 ):
                        FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if( dateRange_padded[cntrAfter,1] > 366 ):
                        FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                    #END IF
                #END IF 
            #-----No leap year detected-----
            else: #no leap year if this
                if( dateRange_padded[cntrAfter,1] > 365 ):
                    FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                #END IF
            #END IF
            
            while( FLG_after == 0 ): #while supports padding multiple years if it ever came to that
                #-----Leap Year Detection-----
                if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                    #-----Leap Year Skipped Detected - next will be 2100-----
                    if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                        #Same alg as no leap year used here
                        dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1, dateRange_padded[cntrAfter,1]-365]); #adjust to the year after
                        #-----Rerun the year check on the next to see if we need to go through this again-----
                        if np.mod(dateRange_padded[i,0],4) == 0: #leap year
                            #-----Leap Year Skipped Detected - next will be 2100-----
                            if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                                #Same alg as no leap year used here
                                if( dateRange_padded[cntrAfter,1] < 366 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            else:
                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                if( dateRange_padded[cntrAfter,1] < 367 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            #END IF
                            
                        #-----No leap year detected-----
                        else: #no leap year if this
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    else:
                    #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                        dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1], dateRange_padded[cntrAfter,1]-366); #adjust to the year after
                        #-----Rerun the year check on the next to see if we need to go through this again-----
                        if np.mod(dateRange_padded[i,0],4) == 0: #leap year
                            #-----Leap Year Skipped Detected - next will be 2100-----
                            if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                                #Same alg as no leap year used here
                                if( dateRange_padded[cntrAfter,1] < 366 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            else:
                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                if( dateRange_padded[cntrAfter,1] < 367 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            #END IF
                        #-----No leap year detected-----
                        else: #no leap year if this
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    #END IF
                #-----No leap year detected-----
                else: #no leap year if this
                    dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1], dateRange_padded[cntrAfter,1]-365); #adjust to the year after
                    #-----Rerun the year check on the next to see if we need to go through this again-----
                    if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                        #-----Leap Year Skipped Detected - next will be 2100-----
                        if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                            #Same alg as no leap year used here
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        else:
                        #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                            if( dateRange_padded[cntrAfter,1] < 367 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    #-----No leap year detected-----
                    else: #no leap year if this
                        if( dateRange_padded[cntrAfter,1] < 366 ):
                            FLG_after = 1; #we're free
                        #END IF
                    #END IF
                #END IF
            #END WHILE
        #END FOR
        
        return(dateRange_padded); #return success
    
    elif( padDirection == 2 ):
    #-----Option 2 returns padding before the date range by the padNumber requested-----
        
        dateRange_padded = np.zeros([len(dateRange[:,0])+padNumber,2],dtype="int16"); #preallocate    
        dateRange_padded[ (padNumber) : (len(dateRange[:,0])+padNumber) ,:] = dateRange; #insert original data
        
        cntrBefore = padNumber; #run seperate counters for before and after
        for i in range(0,padNumber):
            #==============Tackle before padding==============
            dateRange_padded[i,:] = np.array([dateRange[0,0], dateRange[0,1] - cntrBefore]); #fill in before date
            cntrBefore = cntrBefore - 1; #subtract off of the counter
            
            while( dateRange_padded[i,1] < 1 ): #while supports padding multiple years if it ever came to that
                #NOTE: the year checks are run on the previous year since we're going back one!
                #-----Leap Year Detection-----
                if np.mod(dateRange_padded[i,0]-1,4) == 0: #leap year
                    #-----Leap Year Skipped Detected - next will be 2100-----
                    if (np.mod(dateRange_padded[i,0]-1,100) == 0) and (np.mod(dateRange_padded[i,0]-1,400) != 0):
                        #Same alg as no leap year used here
                        dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+365]); #adjust to the year before
                    else:
                    #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                        dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1], dateRange_padded[i,1]+366); #adjust to the year before
                    #END IF
                        
                #-----No leap year detected-----
                else: #no leap year if this
                    dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+365]); #adjust to the year before
                #END IF
            #END WHILE
        #END FOR
            
        return(dateRange_padded); #return success
    
    else:
        #-----Anything else returns the default option of padding on both sides of the date range by the padNumber requested-----
        print("\n==============~Warning~==============");
        #print("Incorrect option number of {} input. Returning in default form of [dayNum, YR]\nAllowed options: 0 == [dayNum/YR], 1 == [YR/dayNum], 2 == [dayNum]\n".format(options)); #.format version
        print("Incorrect padDirection option number of "+str(padDirection)+" input. Returning in default form of padding both sides of the date range.\nAllowed options: 0 == pad both sides, 1 == pad forward, 2 == pad backwards\n"); #non-.format version for nopython calls
        
        dateRange_padded = np.zeros([len(dateRange[:,0])+2*padNumber,2],dtype="int16"); #preallocate    
        dateRange_padded[ (padNumber) : (len(dateRange[:,0])+padNumber) ,:] = dateRange; #insert original data
        
        cntrBefore = padNumber; #run seperate counters for before and after
        FLG_after = 1; #flag for after that facilitates support for padding multiple years forward sucessfully (1 == off shut up)
        for i in range(0,padNumber):
            #==============Tackle before padding==============
            dateRange_padded[i,:] = np.array([dateRange[0,0], dateRange[0,1] - cntrBefore]); #fill in before date
            cntrBefore = cntrBefore - 1; #subtract off of the counter
            
            while( dateRange_padded[i,1] < 1 ): #while supports padding multiple years if it ever came to that
                #NOTE: the year checks are run on the previous year since we're going back one!
                #-----Leap Year Detection-----
                if np.mod(dateRange_padded[i,0]-1,4) == 0: #leap year
                    #-----Leap Year Skipped Detected - next will be 2100-----
                    if (np.mod(dateRange_padded[i,0]-1,100) == 0) and (np.mod(dateRange_padded[i,0]-1,400) != 0):
                        #Same alg as no leap year used here
                        dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+365]); #adjust to the year before
                    else:
                    #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                        dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+366]); #adjust to the year before
                    #END IF
                        
                #-----No leap year detected-----
                else: #no leap year if this
                    dateRange_padded[i,:] = np.array([dateRange_padded[i,0]-1, dateRange_padded[i,1]+365]); #adjust to the year before
                #END IF
            #END WHILE
            
            #==============Tackle after padding==============
            cntrAfter = i+len(dateRange[:,0])+padNumber; #calc forward counter
            dateRange_padded[cntrAfter,:] = np.array([dateRange[-1,0], dateRange[-1,1] + i+1]); #fill in after date
            #NOW running a check on the current year to see if an increment is needed
            #-----Leap Year Detection-----
            if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                    #Same alg as no leap year used here
                    if( dateRange_padded[cntrAfter,1] > 365 ):
                        FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if( dateRange_padded[cntrAfter,1] > 366 ):
                        FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                    #END IF
                #END IF 
            #-----No leap year detected-----
            else: #no leap year if this
                if( dateRange_padded[cntrAfter,1] > 365 ):
                    FLG_after = 0; #flag for after that facilitates support for padding multiple years forward sucessfully
                #END IF
            #END IF
            
            while( FLG_after == 0 ): #while supports padding multiple years if it ever came to that
                #-----Leap Year Detection-----
                if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                    #-----Leap Year Skipped Detected - next will be 2100-----
                    if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                        #Same alg as no leap year used here
                        dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1, dateRange_padded[cntrAfter,1]-365]); #adjust to the year after
                        #-----Rerun the year check on the next to see if we need to go through this again-----
                        if np.mod(dateRange_padded[i,0],4) == 0: #leap year
                            #-----Leap Year Skipped Detected - next will be 2100-----
                            if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                                #Same alg as no leap year used here
                                if( dateRange_padded[cntrAfter,1] < 366 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            else:
                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                if( dateRange_padded[cntrAfter,1] < 367 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            #END IF
                            
                        #-----No leap year detected-----
                        else: #no leap year if this
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    else:
                    #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                        dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1, dateRange_padded[cntrAfter,1]-366]); #adjust to the year after
                        #-----Rerun the year check on the next to see if we need to go through this again-----
                        if np.mod(dateRange_padded[i,0],4) == 0: #leap year
                            #-----Leap Year Skipped Detected - next will be 2100-----
                            if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                                #Same alg as no leap year used here
                                if( dateRange_padded[cntrAfter,1] < 366 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            else:
                            #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                                if( dateRange_padded[cntrAfter,1] < 367 ):
                                    FLG_after = 1; #we're free
                                #END IF
                            #END IF
                        #-----No leap year detected-----
                        else: #no leap year if this
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    #END IF
                #-----No leap year detected-----
                else: #no leap year if this
                    dateRange_padded[cntrAfter,:] = np.array([dateRange_padded[cntrAfter,0]+1, dateRange_padded[cntrAfter,1]-365]); #adjust to the year after
                    #-----Rerun the year check on the next to see if we need to go through this again-----
                    if np.mod(dateRange_padded[cntrAfter,0],4) == 0: #leap year
                        #-----Leap Year Skipped Detected - next will be 2100-----
                        if (np.mod(dateRange_padded[cntrAfter,0],100) == 0) and (np.mod(dateRange_padded[cntrAfter,0],400) != 0):
                            #Same alg as no leap year used here
                            if( dateRange_padded[cntrAfter,1] < 366 ):
                                FLG_after = 1; #we're free
                            #END IF
                        else:
                        #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                            if( dateRange_padded[cntrAfter,1] < 367 ):
                                FLG_after = 1; #we're free
                            #END IF
                        #END IF
                    #-----No leap year detected-----
                    else: #no leap year if this
                        if( dateRange_padded[cntrAfter,1] < 366 ):
                            FLG_after = 1; #we're free
                        #END IF
                    #END IF
                #END IF
            #END WHILE
        #END FOR
        
        return(dateRange_padded); #return success?
        
    #END IF