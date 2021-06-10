#Function to import Kp data from the internet
#RD on 12/10/2018

#GOAL: Get 3 hour Kp indexes and stuff
#expecting: [year,day#min;year,day#max] numerical format. So, say [2013,126;2013,128] for 2013 Day 126 to Day 128. No chars pls
#expecting: monthOutputOption is either 0 or 1. 1 for month form, 0 for day form. 0 day form is default.
#NOTE: Kp is every 3 hours, 0-3, 3-6, 6-9, 9-12, 12-15, 15-18, 18-21, 21-24 UT

import numpy as np
from urllib.request import urlopen
from subfun_dayNum_to_date import subfun_dayNum_to_date
from subfun_strstrNB import strstrNB

#-----Testing variables-----
#import os
###Date range goes Month-Day-Year
##dateRange = np.array([[2013,5,8],[2013,5,10]],dtype="int16"); #dates are in int16 because they can be
#dateRange = np.array([[2013,5,6],[2013,5,8]],dtype="int16"); #dates are in int16 because they can be
##dateRange = np.array([[2014,12,31],[2015,3,1]],dtype="int16"); #for debug, check year success
##dates better go earlier -> later
##print("{}".format(dateRange))
#folder = [os.getcwd()]; #current working directory, leave it as this call usually
#folder.append('E:\Big Data'); #place to save data files to
##folder var structure: 0 = running folder, 1 = data folder
#from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units

def GRITI_import_Kp(dateRange_dayNum_full):

    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #get the full date range
    
    siteURL_base = 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/kp'; #base for Kp index, will tack on dates needed
    
    dateRange_yearRange = np.arange(dateRange_dayNum_full[0,0],dateRange_dayNum_full[-1,0]+1,1,dtype=np.int16); #get the full year range from min to max
    
    for yarr in range(0,dateRange_yearRange.size):
        temp = np.unique( dateRange_full[ (dateRange_yearRange[yarr] == dateRange_dayNum_full[:,0]), 1 ] ); #get number of unique months
        try: #try n fail
            dateRange_yearRange_numberOfMonths = np.vstack( (dateRange_yearRange_numberOfMonths, np.vstack( (np.tile(dateRange_yearRange[yarr],temp.size),temp) ).T ) ); #update the vector
        except(NameError):
            dateRange_yearRange_numberOfMonths = np.vstack( (np.tile(dateRange_yearRange[yarr],temp.size),temp) ).T; #create the vector
        #END IF
    #END FOR yarr
    
    Kp_output = np.zeros( (dateRange_dayNum_full.shape[0],8) ); #preallocate (8 kp's per day)
    cntr = 0; #preallocate
    for i in range(0,dateRange_yearRange_numberOfMonths.shape[0]):
        siteURL_prep = str(dateRange_yearRange_numberOfMonths[i,0]); #get the year
        siteURL = siteURL_base+siteURL_prep[2:len(siteURL_prep)+1]; #put in the last 2 year digits
        temp_YearMonthBase = siteURL_prep[2:len(siteURL_prep)+1]; #create a temp var that holds the year and month info
        siteURL_prep = str(dateRange_yearRange_numberOfMonths[i,1]); #get the month
        if( len(siteURL_prep) == 1 ): #if 1, tack on a 0
            siteURL_prep = '0'+siteURL_prep; #tack it on
        elif( len(siteURL_prep) > 2 ): #shouldn't happen
            import sys
            print('ERROR: Month length greater than 2 chars, reported month is '+siteURL_prep+'\nExiting Import Kp.');
            sys.exit(); #yolo out
        #END IF
        siteURL = siteURL+siteURL_prep+'.tab'; #put in the 2 month digits (padded w/ 0 if needed)
        temp_YearMonthBase = temp_YearMonthBase+siteURL_prep;
        Kp_raw = urlopen(siteURL); #download the data needed
        Kp_rawRead = Kp_raw.read(); #read it into a much more useful format
        charset = Kp_raw.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
        if( charset == None ):
            charset = 'utf-8'; #set it to this
        #END IF
        Kp_data = Kp_rawRead.decode(charset); #"decode" the HTML content so it's more legible
    
        numDaysInMonth = dateRange_full[ (dateRange_yearRange_numberOfMonths[i,0] == dateRange_dayNum_full[:,0]) & (dateRange_yearRange_numberOfMonths[i,1] == dateRange_full[:,1]), 2 ]; #long but figures out what days are in the month and year desired
        
        for j in range(0,numDaysInMonth.size): #run through each day
            temp = str(numDaysInMonth[j]); #get the day
            if( len(temp) == 1 ): #if 1, tack on a 0
                temp = '0'+temp; #tack it on
            elif( len(temp) > 2 ): #shouldn't happen
                import sys
                print('ERROR: Day length greater than 2 chars, reported day is '+temp+'\nExiting Import Kp.');
                sys.crash(); #yolo out
            #END IF
            temp_YearMonthDay = temp_YearMonthBase+temp; #create the Year/Month/Day string
            locale = strstrNB(Kp_data,temp_YearMonthDay)[0]; #get the locale data starts
            Kp_data_temp = Kp_data[locale+8:locale+33+1].split(); #get the Kp data for the day in question (+1 for python)
            for k in range(0,8):
                temp = Kp_data_temp[k]; #get the string data
                if( temp[1] == 'o' ):
                    Kp_output[cntr,k] = int(temp[0]); #write in the value read
                elif( temp[1] == '-' ):
                    Kp_output[cntr,k] = int(temp[0]) - 1/3; #write in the value read
                elif( temp[1] == '+' ):
                    Kp_output[cntr,k] = int(temp[0]) + 1/3; #write in the value read
                #END IF
            #END FOR k
            cntr = cntr + 1; #increment counter to next day line
            
        #END FOR j
    #END FOR i
    
    #finally, make Kp_time
    # Kp_time = np.arange(3,dateRange_dayNum_full.shape[0]*24+3,3)/24 + np.tile(dateRange_dayNum_full[0,1],dateRange_dayNum_full.shape[0]*8); #days, time for Kp measurements (0-3 is assumed to be relevant for hr 3, etc...); 
    Kp_time = np.arange(3,dateRange_dayNum_full.shape[0]*24+3,3)*3600 + np.tile(dateRange_dayNum_full[0,1],dateRange_dayNum_full.shape[0]*8)*86400; #days, time for Kp measurements (0-3 is assumed to be relevant for hr 3, etc...); 
    
    Kp_data = {
        'Kp':Kp_output,
        'time':Kp_time,
        }; #make a dict to hold it all
    
    # return Kp_output, Kp_time
    return Kp_data
#END DEF