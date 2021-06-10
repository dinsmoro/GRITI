#GOAL: Date to Day Number
#RD on 1/30/2019
#
#INPUT: month name in string form
#OUTPUT: month number from 1 to 12
#currently only does it one at a time

def subfun_monthWord_to_num(monthWord):
#monthWord = "Jan";

#==============Convert month word to month #==============
    if( (monthWord.lower() == "january") | (monthWord.lower() == "jan.")  | (monthWord.lower() == "jan") ):
        monthNumber = 1; #get the number
    elif( (monthWord.lower() == "february") | (monthWord.lower() == "feb.")  | (monthWord.lower() == "feb") ):
        monthNumber = 2; #get the number
    elif( (monthWord.lower() == "march") | (monthWord.lower() == "mar.")  | (monthWord.lower() == "mar") ):
        monthNumber = 3; #get the number
    elif( (monthWord.lower() == "april") | (monthWord.lower() == "apr.")  | (monthWord.lower() == "apr") ):
        monthNumber = 4; #get the number
    elif( (monthWord.lower() == "may") ):
        monthNumber = 5; #get the number
    elif( (monthWord.lower() == "june") | (monthWord.lower() == "jun.")  | (monthWord.lower() == "jun") ):
        monthNumber = 6; #get the number
    elif( (monthWord.lower() == "july") | (monthWord.lower() == "jul.")  | (monthWord.lower() == "jul") ):
        monthNumber = 7; #get the number
    elif( (monthWord.lower() == "august") | (monthWord.lower() == "aug.")  | (monthWord.lower() == "aug") ):
        monthNumber = 8; #get the number
    elif( (monthWord.lower() == "september") | (monthWord.lower() == "sep.")  | (monthWord.lower() == "sep") | (monthWord.lower() == "sept") | (monthWord.lower() == "sept.") ):
        monthNumber = 9; #get the number
    elif( (monthWord.lower() == "october") | (monthWord.lower() == "oct.")  | (monthWord.lower() == "oct") ):
        monthNumber = 10; #get the number
    elif( (monthWord.lower() == "november") | (monthWord.lower() == "nov.")  | (monthWord.lower() == "nov") ):
        monthNumber = 11; #get the number
    elif( (monthWord.lower() == "december") | (monthWord.lower() == "dec.")  | (monthWord.lower() == "dec") ):
        monthNumber = 12; #get the number
    else:
        print("\n==============ERROR==============");
        print("{} IS NOT A MONTH NAME WOW\n".format( monthWord) );
        monthNumber = 13; #big oof
    #END IF
    
    return(monthNumber); #return success