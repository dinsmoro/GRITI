#Function to import OMNI data from the internet
#RD on 12/10/2018

'''
GOAL: Get OMNI data
expecting: 
 dateRange_dayNum_full - [year,day#min;year,day#inbetween;year,day#max] numerical format. So, say [2013,126;2013,127;2013,128] for 2013 Day 126 to Day 128. No chars pls
 folder - list with folder[0] being the running folder (where the code is) and folder[1] being where the data is kept (can be the same, or not)
 OMNI_columns_wanted [optional] - list or 1D-array of columns, refer to note below for column #s 
     can override default OMNI_columns_wanted = [0,1,2,3,16,21,27,37,41] with whatever columns you want
 FLG_overwrite [optional, default 0] - 1 downloads data again and overwrites old data, 0 just reads previously saved local data. Has no effect on new data.
NOTE: OMNI is every 1 minute.
NOTE: Reads OMNI data, which is 525600x46 (minutes in a year x 46 parameters reported, including 4 time parameters)
NOTE: https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/hroformat.txt describes OMNI data setup

'''

import numpy as np
from urllib.request import urlopen
import h5py
import os

##-----Testing variables-----
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
#from subfun_date_to_dayNum import subfun_date_to_dayNum
#from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units
#OMNI_columns_wanted = [0,1,2,3,16,21,27,37,41];

def GRITI_import_OMNI(dateRange_dayNum_full, folder, OMNI_columns_wanted = [0,1,2,3,16,18,19,21,27,37,41], FLG_overwrite = 0):
    
    siteURL_base = 'https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/'; #base for Kp index, will tack on dates needed
    #OMNI_columns_wanted = [0,1,2,3,16,21,27,37,41]; #columns wanted for testing
    #use index | OMNI index
    #0  | 0 - yr
    #1  | 1 - day#
    #2  | 2 - hr
    #3  | 3 - min
    #4  | 16 - Bz GSE (nT)
    #5  | 18 - Bz GSM (nT)
    #6  | 19 - RMS SD B scalar (nT)
    #7  | 21 - Flow speed (km/s)
    #8  | 27 - Flow pressure (nPa)
    #9  | 37 - AE Index (nT)
    #10 | 41 - SYM/H Index (nT)
    #see https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/hroformat.txt
    #for more documentation
    OMNI_dict = {'Year':0, #year,
                 'Day Num':1, #daynum, 
                 'Hour':2, #hour,
                 'Min':3, #min, 
                 'ID for IMF SC':4, #
                 'ID for SW Plasma SC':5, #
                 '# Pts in IMF Avg':6, #
                 '# Pts in Plasma Avg':7, #
                 '% Interp':8, #
                 'Time Shift':9, #sec,
                 'Timeshift RMS':10,
                 'Phase Front Normal RMS':11,
                 'Time Delta':12, #sec, time between observations
                 'Field Mag Avg':13, #nT
                 'Bx GSE & GSM':14, #nT, GSE and GSM share this axis
                 'By GSE':15, #nT,
                 'Bz GSE':16, #nT,
                 'By GSM':17, #nT,
                 'Bz GSM':18, #nT,
                 'SD B Scalar RMS':19, #nT,
                 'SD Field Vect RMS':20, #nT,
                 'Vsw':21, #km/s, 'flow speed' 'solar wind speed'
                 'VswX':22, #km/s in GSE coord
                 'VswY':23, #km/s in GSE coord
                 'VswZ':24, #km/s in GSE coord
                 'Proton Density':25, # #/cc,
                 'Temp':26, #K,
                 'Psw':27, #nPa, 'flow pressure' 'solar wind speed'
                 'Elec Field':28, #mV/m,
                 'Plasma Beta':29, #
                 'Alfven Mach Num':30, #
                 'SC X':31, #Re units I think, GSE coords
                 'SC Y':32, #Re units I think, GSE coords
                 'SC Z':33, #Re units I think, GSE coords
                 'BSN Loc X':34, #Re units I think, GSE coords, BSN:bow shock noise
                 'BSN Loc Y':35, #Re units I think, GSE coords, BSN:bow shock noise
                 'BSN Loc Z':36, #Re units I think, GSE coords, BSN:bow shock noise
                 'AE':37, #nT, index
                 'AL':38, #nT, index
                 'AU':39, #nT, index
                 'SYM/D':40, #nT, index
                 'SYM/H':41, #nT, index
                 'ASY/D':42, #nT, index
                 'ASY/H':43, #nT, index
                 'PC(N)':44, #nT, index
                 'Magnetosonic Mach Num':45}; #
    #this one's for plot labels
    OMNI_dictPlot = {0:'Year [Year]', #year,
                     1:'Day Num [Day Num]', #daynum, 
                     2:'Hour [Hour]', #hour,
                     3:'Min [Min]', #min, 
                     4:'ID for IMF SC', #
                     5:'ID for SW Plasma SC', #
                     6:'# Pts in IMF Avg', #
                     7:'# Pts in Plasma Avg', #
                     8:'% Interp', #
                     9:'Time Shift [Sec]', #sec,
                     10:'Timeshift RMS',
                     11:'Phase Front Normal RMS',
                     12:'Time Delta [Sec]', #sec, time between observations
                     13:'Field Mag Avg [nT]', #nT
                     14:'Bx GSE & GSM [nT]', #nT, GSE and GSM share this axis
                     15:'By GSE [nT]', #nT,
                     16:'Bz GSE [nT]', #nT,
                     17:'By GSM [nT]', #nT,
                     18:'Bz GSM [nT]', #nT,
                     19:'SD B Scalar RMS [nT]', #nT,
                     20:'SD Field Vect RMS [nT]', #nT,
                     21:'Vsw [km/s]', #km/s, 'flow speed'
                     22:'VswX [km/s]', #km/s in GSE coord
                     23:'VswY [km/s]', #km/s in GSE coord
                     24:'VswZ [km/s]', #km/s in GSE coord
                     25:'Proton Density [#/cc]', # #/cc,
                     26:'Temp [K]', #K,
                     27:'Psw [nPa]', #nPa, 'flow pressure'
                     28:'Elec Field [mV/m]', #mV/m,
                     29:'Plasma Beta', #
                     30:'Alfven Mach Num', #
                     31:'SC X [Re]', #Re units I think, GSE coords
                     32:'SC Y [Re]', #Re units I think, GSE coords
                     33:'SC Z [Re]', #Re units I think, GSE coords
                     34:'BSN X [Re]', #Re units I think, GSE coords, BSN:bow shock noise
                     35:'BSN Y [Re]', #Re units I think, GSE coords, BSN:bow shock noise
                     36:'BSN Z [Re]', #Re units I think, GSE coords, BSN:bow shock noise
                     37:'AE [nT]', #nT, index
                     38:'AL [nT]', #nT, index
                     39:'AU [nT]', #nT, index
                     40:'SYM/D [nT]', #nT, index
                     41:'SYM/H [nT]', #nT, index
                     42:'ASY/D [nT]', #nT, index
                     43:'ASY/H [nT]', #nT, index
                     44:'PC(N) [nT]', #nT, index
                     45:'Magnetosonic Mach Num'}; #
                 
                 
    dateRange_yearRange = np.arange(dateRange_dayNum_full[0,0],dateRange_dayNum_full[-1,0]+1,1,dtype=np.int16); #get the full year range from min to max
    
    
    #OMNI data cadence is 1 per minute
    #download all that OMNI data
    OMNI_output_raw = np.zeros( (dateRange_dayNum_full.shape[0]*1440,46) ); #preallocate, 46 total entries for the OMNI data
    cntr = 0; #preallocate
    FLG_dataAlreadyDownloaded = 0; #flag to show if data is already downloaded or not
    for i in range(0,dateRange_yearRange.size): #steps through each year
        
        #check to see if data is already downloaded
        OMNI_dataFilePath = folder[1] + '\\OMNI\\omni_min' + str(dateRange_yearRange[i]) + '.h5'; #create the expected data path
        if( os.path.isfile(OMNI_dataFilePath) == 1 ): #check if omni data file exists            
            try: #gonna try to read the file - if we fail, it's borked or something
                with h5py.File(OMNI_dataFilePath, 'r') as OMNIfile:
                    OMNI_data_raw = OMNIfile["/OMNIdata"][:,:].astype("float64"); #import satellite ID
    #                OMNIfile_folders = list(OMNIfile.keys()); #lists the data folders in the HDF5 file
                #END WITH
                FLG_dataAlreadyDownloaded = 1; #set the flag to not downloaded
            except OSError:
                #if this happens then the data is corrupted or something
                FLG_dataAlreadyDownloaded = 0; #set the flag to downloaded
            #END TRY
        #END IF
        
        if( FLG_dataAlreadyDownloaded == 0 ):
            siteURL_prep = 'omni_min'+str(dateRange_yearRange[i])+'.asc'; #put in the year
            
            siteURL = siteURL_base+siteURL_prep; #combine to make the site URL
        
            OMNI_raw = urlopen(siteURL); #download the data needed
            OMNI_rawRead = OMNI_raw.read(); #read it into a much more useful format
            charset = OMNI_raw.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
            if( charset == None ):
                charset = 'utf-8'; #set it to this
            #END IF
            OMNI_data = OMNI_rawRead.decode(charset); #"decode" the HTML content so it's more legible
            
            OMNI_data_split = OMNI_data.split('\n'); #split by the \n line seperator
            OMNI_data_raw = np.zeros( (len(OMNI_data_split)-1,46) ); #preallocate, 46 total entries for the OMNI data
            for j in range(0,len(OMNI_data_split)-1): #-1 cause last split is a blank line
                OMNI_data_raw[j,:] = np.float64(OMNI_data_split[j].split()); #make it numbers n such
            #END FOR j
            
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(OMNI_dataFilePath, 'w') as OMNIfile:
                OMNI_data_raw_record = OMNIfile.create_dataset("OMNIdata", OMNI_data_raw.shape, dtype='float64' ,compression="gzip"); #create dataset for the integers
                OMNI_data_raw_record[...] = OMNI_data_raw; #write that data
            #END WITH
        #END IF
        
        dateRange_full_currentYear = dateRange_dayNum_full[:,0] == dateRange_yearRange[i] ; #get the bits that are for the current year
        
        dateRange_full_allDays = np.isin(OMNI_data_raw[:,1],dateRange_dayNum_full[dateRange_full_currentYear,1]); #get the days of the current year needed
        
        OMNI_output_raw[cntr:cntr+(np.sum(dateRange_full_currentYear)*1440),:] = OMNI_data_raw[dateRange_full_allDays,:]; #copy in all relevant data
        
        cntr = cntr + np.sum(dateRange_full_currentYear)*1440; #increment the counter    
        FLG_dataAlreadyDownloaded = 0; #reset the flag
    #END FOR i
    
    #Weed out bad numbers and such (denoted by varying amounts of 9's)
    #Columns 7, 8, and 9 seem to be solid markers for data being bad - not sure for sure but doin it
    OMNI_output_badData = (OMNI_output_raw[:,6] == 999) | (OMNI_output_raw[:,7] == 999) | (OMNI_output_raw[:,8] == 999); #find the bad data, it's denoted by 999
    OMNI_output_raw = OMNI_output_raw[np.logical_not(OMNI_output_badData),:]; #remove the bad data
    
    OMNI_output = OMNI_output_raw; #get the data desired only, call it a day
    
    return OMNI_output, OMNI_dict, OMNI_dictPlot #finish it out