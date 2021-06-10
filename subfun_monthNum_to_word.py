#GOAL: Date to Day Number
#RD on 8/23/2018
#
#INPUT: month number from 1 to 12
#OUTPUT: month name in string form
#options!: 
#0 = full names [DEFAULT!]
#1 = abreviations with . at end (Jan., May, etc.) - 4 or 3 letters long
#2 = abbreviations with no dot (Jan, May, etc.) - 3 or 4 letters, no dots tho
#3 = ALWAYS 3 letters long abbrev


def subfun_monthNum_to_word(monthNumber, options = 0):
    import numpy as np #import in here I dunno
#monthNumber = np.array([[4,4,4],[4,4,4]],dtype="int16"); #for debugging
#monthNumber = np.array([4,4,4,1],dtype="int16"); #for debugging
#monthNumber = np.reshape(monthNumber,(1,-1)); #for debugging
#monthNumber = 4; #for debugging
#options = 0; #for debugging

#==============Catch input format issues==============
    if( np.isscalar(monthNumber) == 1 ):
        monthNumber = np.atleast_1d(monthNumber); #convert to 1D array of length 1 if scalar, so rest can work
        #because python is weak
    #END IF
    
    if( (np.asarray(np.where(np.asarray(monthNumber.shape) == 1 )).size != 0) and (monthNumber.size  != 1) ): #this test checks to see if there is a dimension that is 1 while the whole size not being 1 (if whole size was 1 then it's a [1,] array yeah that is something)
        #tip toe around asking about dims that don't exist because python is -> weak
        print("\n==============~Warning~==============");
        print("Input size was [{},{}] and format [#,] is needed, adjusting to as [arbitrary,]\n".format( monthNumber.shape[0],monthNumber.shape[1],len(monthNumber[:,0]),len(monthNumber[0,:]) ) );
        monthNumber = np.squeeze(monthNumber); #remove the extra dimension so it goes from [arb,1] or [1,arb] to [arb,] because that's a thing
        print("went here size now {}\n".format(np.asarray(monthNumber.shape).size));
    #END IF
        
    if( np.asarray(monthNumber.shape).size > 1 ):
        #Catch where something formatted completely wrong was sent
        print("\n==============ERROR==============");
        print("Input size was [{},{}] and it needs to be [arbitrary,] - exiting month number to word fun!\n".format( monthNumber.shape[0],monthNumber.shape[1] ) );
        return "No"; #can I? I will
    #END IF
    
    if( np.issubdtype(monthNumber[0], np.int16) == 0 ):
        monthNumber = monthNumber.astype("int16"); #force integer16 formatting for dates because I demand it
    #END IF

#==============Convert month # to month word==============
    monthWord = []; #prep word list, it will be the size of monthNum requested, python doesn't really preallocate strings it seems so w/e
    for i in range(0, len(monthNumber) ):
        if( options == 0 ): #if options is 0, which is the default, full month names
            if( monthNumber[i] == 1 ):
                monthWord.append('January');
            elif( monthNumber[i] == 2 ):
                monthWord.append('February');
            elif( monthNumber[i] == 3 ):
                monthWord.append('March');
            elif( monthNumber[i] == 4 ):
                monthWord.append('April');
            elif( monthNumber[i] == 5 ):
                monthWord.append('May');
            elif( monthNumber[i] == 6 ):
                monthWord.append('June');
            elif( monthNumber[i] == 7 ):
                monthWord.append('July');
            elif( monthNumber[i] == 8 ):
                monthWord.append('August');
            elif( monthNumber[i] == 9 ):
                monthWord.append('September');
            elif( monthNumber[i] == 10 ):
                monthWord.append('October');
            elif( monthNumber[i] == 11 ):
                monthWord.append('November');
            elif( monthNumber[i] == 12 ):
                monthWord.append('December');
            else:
                print("\n==============ERROR==============");
                print("{} IS NOT A MONTH NUMBER WOW\n".format( monthNumber[i]) );
                monthWord.append('No.');
            #END IF
        elif( options == 1 ): #if options is 1, month names with . to denote abbreviation if needed
            if( monthNumber[i] == 1 ):
                monthWord.append('Jan.');
            elif( monthNumber[i] == 2 ):
                monthWord.append('Feb.');
            elif( monthNumber[i] == 3 ):
                monthWord.append('Mar.');
            elif( monthNumber[i] == 4 ):
                monthWord.append('Apr.');
            elif( monthNumber[i] == 5 ):
                monthWord.append('May');
            elif( monthNumber[i] == 6 ):
                monthWord.append('June');
            elif( monthNumber[i] == 7 ):
                monthWord.append('July');
            elif( monthNumber[i] == 8 ):
                monthWord.append('Aug.');
            elif( monthNumber[i] == 9 ):
                monthWord.append('Sept.');
            elif( monthNumber[i] == 10 ):
                monthWord.append('Oct.');
            elif( monthNumber[i] == 11 ):
                monthWord.append('Nov.');
            elif( monthNumber[i] == 12 ):
                monthWord.append('Dec.');
            else:
                print("\n==============ERROR==============");
                print("{} IS NOT A MONTH NUMBER WOW\n".format( monthNumber[i]) );
                monthWord.append('No.');
            #END IF
        elif( options == 2 ): #if options is 2, month names with no . so all are 3 letters long
            if( monthNumber[i] == 1 ):
                monthWord.append('Jan');
            elif( monthNumber[i] == 2 ):
                monthWord.append('Feb');
            elif( monthNumber[i] == 3 ):
                monthWord.append('Mar');
            elif( monthNumber[i] == 4 ):
                monthWord.append('Apr');
            elif( monthNumber[i] == 5 ):
                monthWord.append('May');
            elif( monthNumber[i] == 6 ):
                monthWord.append('June');
            elif( monthNumber[i] == 7 ):
                monthWord.append('July');
            elif( monthNumber[i] == 8 ):
                monthWord.append('Aug');
            elif( monthNumber[i] == 9 ):
                monthWord.append('Sept');
            elif( monthNumber[i] == 10 ):
                monthWord.append('Oct');
            elif( monthNumber[i] == 11 ):
                monthWord.append('Nov');
            elif( monthNumber[i] == 12 ):
                monthWord.append('Dec');
            else:
                print("\n==============ERROR==============");
                print("{} IS NOT A MONTH NUMBER WOW\n".format( monthNumber[i]) );
                monthWord.append('No.');
            #END IF
        #END IF
        elif( options == 3 ): #if options is 3, constant length of 3 is forced
            if( monthNumber[i] == 1 ):
                monthWord.append('Jan');
            elif( monthNumber[i] == 2 ):
                monthWord.append('Feb');
            elif( monthNumber[i] == 3 ):
                monthWord.append('Mar');
            elif( monthNumber[i] == 4 ):
                monthWord.append('Apr');
            elif( monthNumber[i] == 5 ):
                monthWord.append('May');
            elif( monthNumber[i] == 6 ):
                monthWord.append('Jun');
            elif( monthNumber[i] == 7 ):
                monthWord.append('Jul');
            elif( monthNumber[i] == 8 ):
                monthWord.append('Aug');
            elif( monthNumber[i] == 9 ):
                monthWord.append('Sep');
            elif( monthNumber[i] == 10 ):
                monthWord.append('Oct');
            elif( monthNumber[i] == 11 ):
                monthWord.append('Nov');
            elif( monthNumber[i] == 12 ):
                monthWord.append('Dec');
            else:
                print("\n==============ERROR==============");
                print("{} IS NOT A MONTH NUMBER WOW\n".format( monthNumber[i]) );
                monthWord.append('No.');
            #END IF
        #END IF
    #END FOR
    
    return(monthWord); #return success