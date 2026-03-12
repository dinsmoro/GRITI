#GOAL: Date to Day Number
#RD on 8/23/2018
#
#INPUT: date in [YRstart/M/D , YRend/M/D] format, as many as you like - must be [#dates , 3] sized
#OUTPUT: dateRange_dayNum as [YR/dayNum] format for each date input
#options!: 
#0 = Output as [dayNum/YR] 
#1 = Output as [YR/dayNum] [DEFAULT!]
#2 = Output as just [dayNum]
import numpy as np #import in here I dunno
#from numba import jit, prange

#@jit(parallel=True, fastmath=True) #try to activate some parallel speed up stuff (doesn't work really) #,nogil=True, parallel=True, fastmath=True
def subfun_date_to_dayNum(dateRange,options = 1):
    
#dateRange = np.array([[2014,4,7],[2014,4,7],[2014,4,8],[2014,4,8],[2014,4,9]],dtype="int16"); #for debugging
##dateRange = np.tile(np.array([[2014,7,30]],dtype="int16"), (100,1) ); #for debugging
##dateRange = np.tile(np.array([[2014],[7],[30]],dtype="int16"), (1,100) ); #for debugging
#options = 1; #for debugging
    
#==============Catch input format issues==============
    if( isinstance(dateRange,tuple) | isinstance(dateRange,list) ):
        dateRange = np.array(dateRange); #convert to numpy array
    #END IF
    if( dateRange.ndim == 1 ):
        dateRange = dateRange[..., np.newaxis]; #add on a new axis b/c code demands it
    #END IF

    if (len(dateRange[0,:]) != 3) and (len(dateRange[:,0]) == 3):
        #Catch where a [3,arbitrary] was sent but this needs to operate on a [arbitrary,3] sized array
        # print("\n==============~Warning~==============");
        # #print("Input size was [{},{}] and it is assumed that the dimensions are flipped -> adjusting to [{},{}] as [arbitrary,3] format is needed.\n".format( len(dateRange[0,:]),len(dateRange[:,0]),len(dateRange[:,0]),len(dateRange[0,:]) ) ); #.format version
        # print("Input size was ["+str(len(dateRange[0,:]))+","+str(len(dateRange[:,0]))+"] and it is assumed that the dimensions are flipped -> adjusting to ["+str(len(dateRange[:,0]))+","+str(len(dateRange[0,:]))+"] as [arbitrary,3] format is needed.\n"); #non-.format version for nopython calls
        dateRange = np.swapaxes(dateRange,0,1); #flip axes
        
    elif (len(dateRange[0,:]) != 3) and (len(dateRange[:,0]) != 3):
        #Catch where something formatted completely wrong was sent
        print("\n==============ERROR==============");
        # print("Input size was [{},{}] and it needs to be [arbitrary,3] - exiting date to dayNum fun!\n".format( len(dateRange[0,:]),len(dateRange[:,0]) ) ); #.format version
        print("Input size was ["+str(len(dateRange[0,:]))+","+str(len(dateRange[:,0]))+"}] and it needs to be [arbitrary,3] - exiting date to dayNum fun!\n"); #non-.format version for nopython calls
        # return("No"); #can I? I will
        import sys
        sys.crash(); #even better
    #END IF

#Split code base into the seperate if's for a possible speed-up
#==============Convert dates given in M/D/YR format to dayNum/YR format==============    
    if( options == 0 ):
    #-----Option 0 returns the default output of [dayNum/Yr] format-----
        
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) , 2 ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) , 2 ) ) ); #preallocate , numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1), : ] = np.array([tempMonth + dateRange_unique[i,2],dateRange_unique[i,0]]); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success
    
    elif( options == 1 ):
    #-----Option 1 returns the output in the [Yr/dayNum] format-----
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) , 2 ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) , 2 ) ) ); #preallocate, numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1), : ] = np.array([dateRange_unique[i,0],tempMonth + dateRange_unique[i,2]]); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success
    
    elif( options == 2 ):
    #-----Option 2 returns the output in the [dayNum] only format-----
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) ) ) ); #preallocate , numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1) ] = np.array( [tempMonth + dateRange_unique[i,2]] ); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success
    
    else:
        #-----Anything else returns the default output [dayNum/Yr] format-----
        print("\n==============~Warning~==============");
        #print("Incorrect option number of {} input. Returning in default form of [dayNum, YR]\nAllowed options: 0 == [dayNum/YR], 1 == [YR/dayNum], 2 == [dayNum]\n".format(options)); #.format version
        print("Incorrect option number of "+str(options)+" input. Returning in default form of [YR, dayNum]\nAllowed options: 0 == [dayNum/YR], 1 == [YR/dayNum], 2 == [dayNum]\n"); #non-.format version for nopython calls
        
        #dateRange_dayNum = np.zeros( (len(dateRange[:,0]) , 2 ), dtype=np.int16 ); #preallocate 
        dateRange_dayNum = np.int16(np.zeros( (len(dateRange[:,0]) , 2 ) ) ); #preallocate , numba happy
        tempMonth = np.int16(0); #prep to make numba happy
        dateRange_unique = np.unique(dateRange, axis = 0); #get the unique dates involved (helps if there are many repeats) 
    
        for i in range(0, len(dateRange_unique[:,0]) ):
            #-----Leap Year Detection-----
            if np.mod(dateRange_unique[i,0],4) == 0: #leap year
                #-----Leap Year Skipped Detected - next will be 2100-----
                if (np.mod(dateRange_unique[i,0],100) == 0) and (np.mod(dateRange_unique[i,0],400) != 0):
                    #Same alg as no leap year used here
                    if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 28
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(59); #prev month's contribution #31+28
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(90); #prev month's contribution #31+28+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                    #END IF
                else:
                #-----Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)-----
                    if dateRange_unique[i,1] == 1: #Jan 31
                        #Nothing happens, days count up here
                        tempMonth = np.int16(0); #no contribution
                    elif dateRange_unique[i,1] == 2: #Feb 29 - leap year
                        tempMonth = np.int16(31); #jan's contribution
                    elif dateRange_unique[i,1] == 3: #March 31
                        tempMonth = np.int16(60); #prev month's contribution #31+29
                    elif dateRange_unique[i,1] == 4: #April 30
                        tempMonth = np.int16(91); #prev month's contribution #31+29+31
                    elif dateRange_unique[i,1] == 5: #May 31
                        tempMonth = np.int16(121); #prev month's contribution #31+29+31+30
                    elif dateRange_unique[i,1] == 6: #June 30
                        tempMonth = np.int16(152); #prev month's contribution #31+29+31+30+31
                    elif dateRange_unique[i,1] == 7: #July 31
                        tempMonth = np.int16(182); #prev month's contribution  #31+29+31+30+31+30
                    elif dateRange_unique[i,1] == 8: #August 31
                        tempMonth = np.int16(213); #prev month's contribution #31+29+31+30+31+30+31)
                    elif dateRange_unique[i,1] == 9: #September 30
                        tempMonth = np.int16(244); #prev month's contribution #31+29+31+30+31+30+31+31
                    elif dateRange_unique[i,1] == 10: #Oct 31
                        tempMonth = np.int16(274); #prev month's contribution #31+29+31+30+31+30+31+31+30
                    elif dateRange_unique[i,1] == 11: #Nov 30
                        tempMonth = np.int16(305); #prev month's contribution #31+29+31+30+31+30+31+31+30+31
                    elif dateRange_unique[i,1] == 12: #Dec 31
                        tempMonth = np.int16(335); #prev month's contribution #31+29+31+30+31+30+31+31+30+31+30
                    #END IF
                #END IF
                    
            #-----No leap year detected-----
            else: #no leap year if this
                if dateRange_unique[i,1] == 1: #Jan 31
                    #Nothing happens, days count up here
                    tempMonth = np.int16(0); #no contribution
                elif dateRange_unique[i,1] == 2: #Feb 28
                    tempMonth = np.int16(31); #jan's contribution
                elif dateRange_unique[i,1] == 3: #March 31
                    tempMonth = np.int16(59); #prev month's contribution #31+28
                elif dateRange_unique[i,1] == 4: #April 30
                    tempMonth = np.int16(90); #prev month's contribution #31+28+31
                elif dateRange_unique[i,1] == 5: #May 31
                    tempMonth = np.int16(120); #prev month's contribution #31+28+31+30
                elif dateRange_unique[i,1] == 6: #June 30
                    tempMonth = np.int16(151); #prev month's contribution #31+28+31+30+31
                elif dateRange_unique[i,1] == 7: #July 31
                    tempMonth = np.int16(181); #prev month's contribution #31+28+31+30+31+30
                elif dateRange_unique[i,1] == 8: #August 31
                    tempMonth = np.int16(212); #prev month's contribution #31+28+31+30+31+30+31
                elif dateRange_unique[i,1] == 9: #September 30
                    tempMonth = np.int16(243); #prev month's contribution #31+28+31+30+31+30+31+31
                elif dateRange_unique[i,1] == 10: #Oct 31
                    tempMonth = np.int16(273); #prev month's contribution #31+28+31+30+31+30+31+31+30
                elif dateRange_unique[i,1] == 11: #Nov 30
                    tempMonth = np.int16(304); #prev month's contribution #31+28+31+30+31+30+31+31+30+31
                elif dateRange_unique[i,1] == 12: #Dec 31
                    tempMonth = np.int16(334); #prev month's contribution #31+28+31+30+31+30+31+31+30+31+30
                #END IF
            #END IF
            
            dateRange_dayNum[ np.all(dateRange == dateRange_unique[i,:],axis=1), : ] = np.array([dateRange_unique[i,0],tempMonth + dateRange_unique[i,2]]); #Adds the days up to the month plus the days in the month, puts them where they match
        #END FOR
        
        return(dateRange_dayNum); #return success?
    #END IF