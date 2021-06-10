#GOAL: Find all occurances of a string B within a string A
#RD on 9/15/2020
#Greatly inspired by https://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
#Syntax looks like & acts like C
#
#INPUT: string to find (B) and string to search (A)
#OUTPUT: Indexes of found string (B) in string (A)

import numpy as np #import numpy
#tried parallel using simple numba + jit/prange and it just got worse - maybe explicit if it mattered would be better

#strstrNB means NO BACKSLASH support, which seems important in some cases
def strstrNB(strSearch, strFind):
    # strSearch = "%r"%strSearch; #this is a workaround to allow for finding of \ in \n [this is very important don't disable it]
    #this requires strSearch[1:-1] to be used
    
    cntr = 0; #counter
    strFound = []; #init list that holds found indexes per string
    while cntr != -1:
        #use [1:-1] for "%r" mode above
        # cntr = strSearch[1:-1].find(strFind, cntr); #search string for strFind (B) starting at cntr, which essentially moves the search up the strSearch (A) by 1 each loop
        cntr = strSearch.find(strFind, cntr); #search string for strFind (B) starting at cntr, which essentially moves the search up the strSearch (A) by 1 each loop
        if cntr != -1: #watch for when strFind (B) is not found
            strFound.append(cntr); #record cntr
            cntr = cntr + 1; #Increment to next char
            #To make it go faster, use +len(strFind) instead of +1 - but it will not find overlapping strings
        #END IF
    #END WHILE
    
    strFound = np.asarray(strFound); #I hate lists 
    return strFound; #time to exit when -1 hits as that means no string was found