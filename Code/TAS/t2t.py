#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
time2time: t2t

Rely on datetime to get seconds since 19whatever

!! Not finished !!
"""

# offset time in seconds...
# POSIX -> no leap seconds
# UTC -> should have leap seconds
# TAI -> immune to leap seconds

import numpy as np
from os import path as ospath
from re import compile as recompile
import datetime
from subfun_date_to_dayNum import subfun_date_to_dayNum
from subfun_monthWord_to_num import subfun_monthWord_to_num

# J2000 epoch is in TT time; need to convert to TAI w/ formula
# J2000_epoch = datetime.datetime(2000, 1, 1, 11, 58, 55, 816000, tzinfo=datetime.timezone.utc); # J2000 epoch (afaik per https://aa.usno.navy.mil/faq/TT)


# Pulled from GRITI_import_ISR_Haystack.py
try: #try to open the leap second file
    leapSec = np.load('.leapsec.npy', allow_pickle=False)
    
    #now gotta make sure file we have is current enough
    leapCurrentDay = subfun_date_to_dayNum(np.int16(np.array(datetime.date.today().strftime("%Y %m %d").split(" "))))[0]; #get current day
    leapExpire = leapSec[0, 2:]; #create a number version
    
    if( leapExpire[0] < leapCurrentDay[0] ): #if the expiry year is less than the current one, easy catch
        raise ValueError("Raising an error so the except triggers and it goes and gets new data from the website, as the current data is old! Current date:\n"+str(leapCurrentDay)+"\nExpiring date:\n"+str(leapExpire));
    elif( leapExpire[0] == leapCurrentDay[0] ): #if the year is the same
        if( leapExpire[1] < leapCurrentDay[1] ): #check the days within them
            raise ValueError("Raising an error so the except triggers and it goes and gets new data from the website, as the current data is old! Current date:\n"+str(leapCurrentDay)+"\nExpiring date:\n"+str(leapExpire));
    #END IF
except: #fail and do this to get it from 
    from urllib.request import urlopen #only need it here
    regexr_exp = recompile(r'[0-9]+ \w+ [0-9]+$'); # Prep it
    regexr_line = recompile(r'[0-9]+\s+#\s[0-9]+\s\w+\s[0-9]+$'); # Prep it
    
    web_leapSecond = "https://data.iana.org/time-zones/data/leap-seconds.list"; #build site to go to
    web_page = urlopen(web_leapSecond); #get raw HTML
    web_htmlContent = web_page.read().decode("UTF-8").split("\n"); #read off the HTML from whatever the page holder is
    leapStart = -1; #prep flag/counter recorders
    leapEnd = -1;
    leapExpire = -1;
    for i in range(0,len(web_htmlContent)):
        if( len(web_htmlContent[i]) > 0 ): #prevent "none lines" from being read
            if( (web_htmlContent[i][0] != "#" ) & ( leapStart == -1 ) ): #catch the first instance of no #'s
                leapStart = i; #record the place where the comments end
            #END IF
            if( (web_htmlContent[i][0] == "#" ) & ( leapStart != -1 ) & ( leapEnd == -1 ) ): #catch where the #'s start again
                leapEnd = i-1; #record the place where the comments start again
            #END IF
            if( "File expires" in web_htmlContent[i] ): #catch the file expires line
                leapExpire = i; #record that line too
            # END IF
        #END IF
    #END FOR i
    
    leapSec = np.empty( (leapEnd-leapStart+2, 4), dtype=np.int16); # Prep the leap sec holder
    leapStrings = web_htmlContent[leapStart:(leapEnd+1)]; #pull out the leap strings
    leapExpire = regexr_exp.search(web_htmlContent[leapExpire]).group().split(' ');
    leapSec[0, 2:] = subfun_date_to_dayNum( [int(leapExpire[2]), subfun_monthWord_to_num( leapExpire[1] ), int(leapExpire[0])] )[0, :]; # Write in the expiry date
    for i in range(0, len(leapStrings)):
        bonk = regexr_line.search(leapStrings[i]).group().split(); # Bonk it
        
        leapSec[i+1, :] = [int(bonk[0]), int(bonk[0]), int(bonk[4]), subfun_date_to_dayNum( [int(bonk[4]), subfun_monthWord_to_num( bonk[3] ), int(bonk[2])] , options=2).item()]
    # END FOR i
    leapSec[1:, 1] = leapSec[1:,1] - leapSec[1,1]; # Reference to 1970 whatever
    
    np.save('.leapsec.npy', leapSec, allow_pickle=False)
#END TRY

#Go hard
#secondsSince1970 = 

secSince_yearDiff = np.arange(1970,dateRange_dayNum_zeroHr[0],1,dtype="int16") #get the year range from 1970 to dataYear-1 b/c slicing but it's OK here, I guess.
secSince_yearDiffDays = np.zeros(len(secSince_yearDiff),dtype="int16"); #use to record the days in the years
for i in range(0, len(secSince_yearDiff) ):
    #Leap Year Detection
    if np.mod(secSince_yearDiff[i],4) == 0: #leap year
        #LEAP YEAR! Possibly.
        if (np.mod(secSince_yearDiff[i],100) == 0) and (np.mod(secSince_yearDiff[i],400) != 0):
            #NO LEAP YEAR
            #Leap Year Skipped Detected - next will be 2100
            secSince_yearDiffDays[i] = 365; #no leap year
        else:
            #Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
            secSince_yearDiffDays[i] = 366; #it was a leap year
        #END IF
    else: #no leap year if this
        #no leap year
        secSince_yearDiffDays[i] = 365; #no leap year
    #END IF 
#END FOR i loop per year
    
secSince = (np.sum(secSince_yearDiffDays) + (dateRange_dayNum_zeroHr[1]-1))*86400; #sec, seconds since 1970 Jan 1st.
#-1 on the day so we don't count the seconds in the current day - only up to the current day!

#NOW ADD IN SOME LEAP SECONDS IT'S SUPER EASY
secSince_leapWhere = np.where(leapDates_dayNum[:,0] < dateRange_dayNum_zeroHr[0])[0]; #get leap seconds that are relevant for the date
secSince_leapWhereCheck = np.where(leapDates_dayNum[:,0] == dateRange_dayNum_zeroHr[0])[0]; #get leap seconds that are relevant for the date
for i in range(0, len(secSince_leapWhereCheck)):
    if( leapDates_dayNum[secSince_leapWhereCheck[i],1] <= dateRange_dayNum_zeroHr[1] ): #find days that are before or on the start day
        secSince_leapWhereCheck[i] = -1; #use this flag to remove them as those leap seconds happened
    #END IF
#END FOR i
    
secSince_leapWhereCheck = np.delete(secSince_leapWhereCheck,np.where(secSince_leapWhereCheck == -1)); #delete ones that are ok
if( len(secSince_leapWhereCheck) != 0 ): #if it's not empty, there's something to delete
    secSince_leapWhereDel = np.empty(len(secSince_leapWhereCheck),dtype=np.int64); #preallocate
    for i in range(0, len(secSince_leapWhereCheck)):
        secSince_leapWhereDel[i] = np.where(secSince_leapWhereCheck[i] == secSince_leapWhere)[0];
    #END FOR i
    secSince_leapWhere = np.delete(secSince_leapWhere,secSince_leapWhereDel); #delete the ones from the list
#END IF
    
secSince = secSince + len(secSince_leapWhere); #add in the number of leap seconds needed