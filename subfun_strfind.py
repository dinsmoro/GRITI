#GOAL: Search an array or list of strings for a single string
#RD on 8/27/2018
#Greatly inspired by https://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
#Syntax looks like & acts like C
#
#INPUT: array or list of strings to search (A), string to find in said array/list (B)
#OUTPUT: Boolean array of true/false, opt==1 causes the number of times it was found
#option opt, 0 == returns array size of strSearch with # occurances found, 1 == returns the total number of occurances found as 1 integer value

import numpy as np #import numpy
from subfun_strstr import strstr

def strfind(strSearch, strFind, opt=0):
    strSearch = np.array(strSearch); #convert to array
    
    strFound = np.zeros(strSearch.size,dtype=np.int64); #preallocate array to find
    for i in range(0,strSearch.size):
        strFound[i] = strstr(strSearch[i],strFind).size; #record how many times it was found
    #END FOR
    
    if(opt == 1):
        strFound = np.sum(strFound); #sum em up
    #END IF
    
    return strFound #returns how many were found where (or just how many were found)
#END DEF
    
from subfun_strstr import strstrNB

def strfindNB(strSearch, strFind, opt=0):
    strSearch = np.array(strSearch); #convert to array
    
    strFound = np.zeros(strSearch.size,dtype=np.int64); #preallocate array to find
    for i in range(0,strSearch.size):
        strFound[i] = strstrNB(strSearch[i],strFind).size; #record how many times it was found
    #END FOR
    
    if(opt == 1):
        strFound = np.sum(strFound); #sum em up
    #END IF
    
    return strFound #returns how many were found where (or just how many were found)
#END DEF