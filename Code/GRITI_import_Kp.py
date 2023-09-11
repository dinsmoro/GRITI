#Function to import Kp data from the internet
#RD on 12/10/2018

#GOAL: Get 3 hour Kp indexes and stuff
#expecting: [year,day#min;year,day#max] numerical format. So, say [2013,126;2013,128] for 2013 Day 126 to Day 128. No chars pls
#expecting: monthOutputOption is either 0 or 1. 1 for month form, 0 for day form. 0 day form is default.
#NOTE: Kp is every 3 hours, 0-3, 3-6, 6-9, 9-12, 12-15, 15-18, 18-21, 21-24 UT

import numpy as np
import os
from datetime import datetime
import requests #requests doesn't seem to have cert errors
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date

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
#from Code.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units

def GRITI_import_Kp(dateRange_dayNum_full, settings_paths_data):

    dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #get the full date range
    
    # siteURL_base = 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/kp'; #base for Kp index, will tack on dates needed
    siteURL_base = 'https://datapub.gfz-potsdam.de/download/10.5880.Kp.0001/Kp_definitive/'; #FTP site funked up it seems
    
    path_base = os.path.join(settings_paths_data,'Kp'); #prep the base path
    if( not os.path.exists(path_base) ):
        print('NOTA BENE: GRITI_import_Kp - Created Kp data directory: '+path_base+'\n');
        os.makedirs(path_base); #create directory
    #END IF
    
    current_year = datetime.now().year; #get the current year (needed if Kp needs to be updated since the file isn't set yet)
    
    dateRange_yearRange = np.arange(dateRange_dayNum_full[0,0],dateRange_dayNum_full[-1,0]+1,1,dtype=np.int16); #get the full year range from min to max
        
    Kp_output = np.zeros( (dateRange_dayNum_full.shape[0],8) ); #preallocate (8 kp's per day)
    for i in range(0,dateRange_yearRange.size):
        Kp_fileName = 'Kp_def'+str(dateRange_yearRange[i])+'.wdc'; #get the expected file name
        path_year = os.path.join(path_base,Kp_fileName); #don't need to put Kp data in year folders, 1 file does it all
        
        #download if not downloaded, or if current year is requested year (data ain't done)
        if( (not os.path.exists(path_year)) | (dateRange_yearRange[i] == current_year) ):
            siteURL = siteURL_base+Kp_fileName; #get expected path
            Kp_raw = requests.get(siteURL, stream=True); #get the data, stream is key I think
            # Kp_data = Kp_raw.text; #get the text (requests is used up so can't get text also)
            with open(path_year, 'wb') as Kp_file:
                for data in Kp_raw.iter_content():
                    Kp_file.write(data)
                #END FOR data
            #END WITH
        #END IF
                
        with open(path_year, 'r') as Kp_file:
            Kp_dataText = []; #prep a list
            for line in Kp_file:
                if( line[0] != '#' ): #don't read commented lines as data
                    Kp_dataText.append(line.rstrip()); #.rstrip() ditches new line symbols
                #END IF
            #END FOR line
        #END WITH
                
        def WDC_reader(text):
            data = {}; #prep a dict
            for line in text:
                month = int(line[2:4]);
                day = int(line[4:6]);
                if( month not in data ):
                    data[month] = {}; #create a sub-dict
                #END IF
                if( day not in data[month] ):
                    data[month][day] = {}; #create a sub-sub-dict
                #END IF
                # data[month][day]['solar rotation #'] = int(line[6:10]); #Bartels solar rotation number (a sequence of 27-day intervals counted from February 8, 1832)
                # data[month][day]['solar rotation # day#'] = int(line[10:12]); #Number of day within the Bartels solar rotation ( 1-27)
                data[month][day]['Kp'] = line[12:28]; #Kp indices: the planetary three-hour index Kp for the intervals 0-3, 3-6, 6-9, 9-12, 12-15, 15-18, 18-21 and 21-24 UT
                data[month][day]['Kp'] = [data[month][day]['Kp'][j:j+2] for j in range(0, len(data[month][day]['Kp']), 2)]; #split line every 2
                Kp_line = np.empty(len(data[month][day]['Kp'])); #prep
                for j in range(0,len(data[month][day]['Kp'])):
                    if( data[month][day]['Kp'][j][0] != ' ' ):
                        Kp_line[j] = np.float64(data[month][day]['Kp'][j][0]);
                    else:
                        Kp_line[j] = 0.;
                    #END IF
                    if( data[month][day]['Kp'][j][1] == '3' ):
                        Kp_line[j] = Kp_line[j] + 1/3; #write in the value read
                    elif( data[month][day]['Kp'][j][1] == '7' ):
                        Kp_line[j] = Kp_line[j] + 2/3; #write in the value read
                    #END IF
                #END FOR j
                data[month][day]['Kp'] = np.copy(Kp_line); #load in the converted numpy array
                # Kp_line = 0.; #prep a number
                # data[month][day]['Kp daily sum'] = line[28:31]; #Daily Kp sum rounded to thirds, supplied only for historic reasons, use Ap for scientific studies
                # if( data[month][day]['Kp daily sum'][0:2] != '  ' ):
                #     Kp_line = np.float64(data[month][day]['Kp daily sum'][0:2]);
                # else:
                #     Kp_line = 0.;
                # #END IF
                # if( data[month][day]['Kp daily sum'][2] == '3' ):
                #     Kp_line = Kp_line + 1/3; #write in the value read
                # elif( data[month][day]['Kp daily sum'][2] == '7' ):
                #     Kp_line = Kp_line + 2/3; #write in the value read
                # #END IF
                # data[month][day]['Kp daily sum'] = Kp_line; #load in the converted
                # data[month][day]['Ap'] = line[31:55]; #ap indices: the three-hourly equivalent planetary amplitude for the intervals 0-3, 3-6, 6-9, 9-12, 12-15,  5-18, 18-21 and 21-24 UT, in units of 2 nT
                # data[month][day]['Ap'] = [data[month][day]['Ap'][j:j+3] for j in range(0, len(data[month][day]['Ap']), 3)]; #split line every 3
                # Kp_line = np.empty(len(data[month][day]['Ap']), dtype=np.int64); #prep
                # for j in range(0,len(data[month][day]['Ap'])):
                #     Kp_line[j] = np.int64(data[month][day]['Ap'][j])*2; #nT, Ap index
                # #END FOR j
                # data[month][day]['Ap'] = np.copy(Kp_line); #load in the converted numpy array
                # data[month][day]['Ap daily sum'] = int(line[55:58]); #Ap index: the daily equivalent planetary amplitude, the arithmetic mean of the 8 ap values, rounded to integer, in units of 2 nT
                # data[month][day]['Cp'] = float(line[58:61]); #Cp index: the daily planetary character figure, a qualitative estimate of the overall level of geomagnetic activity for this day determined from the sum of the eight ap amplitudes, ranging, in steps of 0.1, from 0.0 to 2.5
                # data[month][day]['C9'] = int(line[61]); #C9: the contracted scale for Cp with only 1 digit, from 0 to 9
            #END FOR line
            return data
        #END DEF
        
        Kp_data = WDC_reader(Kp_dataText); #convert to something useful
        
        daysInTheYear = np.where(dateRange_yearRange[i] == dateRange_dayNum_full[:,0])[0]; #indices of days in the year to load in
        
        for j in range(0,daysInTheYear.size):
            Kp_output[daysInTheYear[j],:] = Kp_data[dateRange_full[daysInTheYear[j],1]][dateRange_full[daysInTheYear[j],2]]['Kp']; #get out the Kps
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