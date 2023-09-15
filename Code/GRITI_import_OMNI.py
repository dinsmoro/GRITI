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
#from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
#from Code.subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange
#dateRange_dayNum = subfun_date_to_dayNum(dateRange); #call function to get the date range into dayNumber form (easy to work with)
#(_, dateRange_dayNum_full) = subfun_dateORdayNum_to_fullRange(dateRange_dayNum); #call fun to get fully enumerated days between range
#dateRange_dayNum_zeroHr = dateRange_dayNum_full[np.int16( np.floor((len(dateRange_dayNum_full[:,0]) - 1)/2) ),:]; #choose day for when "0 hr" is - makes plots nice, no day units just hr units
#OMNI_columns_wanted = [0,1,2,3,16,21,27,37,41];

def GRITI_import_OMNI(dateRange_dayNum_full, settings_paths, OMNI_columns_wanted = [0,1,2,3,16,18,19,21,27,37,41], FLG_overwrite = 0):
    
    if( os.path.isdir(os.path.join(settings_paths['data'], 'OMNI')) == False ): #check if TEC folder exists
        #if not, make it
        os.makedirs(os.path.join(settings_paths['data'], 'OMNI'));
        print('NOTA BENE: In GRITI_import_OMNI - Created OMNI directory: '+os.path.join(settings_paths['data'], 'OMNI')+'\n');
    #END IF
    
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
                 'PC(N)':44, #mV/m, index
                 'Magnetosonic Mach Num':45,
                 'IMF clock angle':46, #this one's a derived product (angle between Bz and By (but accounts for real angle w/ Bx) w/ 0deg at Bz+ (geomag north))
                 'IMF clock mag':47,#this one's a derived product (angle between Bz and By (but accounts for real angle w/ Bx) w/ 0deg at Bz+ (geomag north))
                 'SYM/H diff':48, #derived delta of SYM/H
                 'Epsilon':49, #GW, per Akasofu, S.I. Energy coupling between the solar wind and the magnetosphere. Space Sci Rev 28, 121–190 (1981). https://doi.org/10.1007/BF00218810
                 'PC(S)':50, #mV/m, not actually included in OMNI but the code will attempt to load it in
                 }; 
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
                     25:'Proton Density [#/cc]', # #/cc, shorter alt: Proton ρ [#/cc]
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
                     44:'PCN [mV/m]', #mV/m, index
                     45:'Magnetosonic Mach Num',
                     46:'IMF ClkAngl [deg]', #this one's a derived product (angle between Bz and By (but accounts for real angle w/ Bx) w/ 0deg at Bz+ (geomag north))
                     47:'IMF ClkMag',
                     48:'SYM/H Diff [nT]',
                     49:'Epsilon',
                     50:'PCS [mV/m]',
                     }; #this one's a derived product (angle between Bz and By (but accounts for real angle w/ Bx) w/ 0deg at Bz+ (geomag north))
                 
                 
    dateRange_yearRange = np.arange(dateRange_dayNum_full[0,0],dateRange_dayNum_full[-1,0]+1,1,dtype=np.int16); #get the full year range from min to max
    
    
    #OMNI data cadence is 1 per minute
    #download all that OMNI data
    OMNI_output_raw = np.zeros( (dateRange_dayNum_full.shape[0]*1440,46) ); #preallocate, 46 total entries for the OMNI data
    cntr = 0; #preallocate
    FLG_dataAlreadyDownloaded = 0; #flag to show if data is already downloaded or not
    for i in range(0,dateRange_yearRange.size): #steps through each year
        
        #check to see if data is already downloaded
        OMNI_dataFilePath = settings_paths['data'] + '\\OMNI\\omni_min' + str(dateRange_yearRange[i]) + '.h5'; #create the expected data path
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
        
            try:
                OMNI_raw = urlopen(siteURL); #download the data needed
                OMNI_rawRead = OMNI_raw.read(); #read it into a much more useful format
                charset = OMNI_raw.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                if( charset == None ):
                    charset = 'utf-8'; #set it to this
                #END IF
                OMNI_data = OMNI_rawRead.decode(charset); #"decode" the HTML content so it's more legible
            except:
                import requests #requests doesn't seem to have cert errors, use it if urlopen is complainy
                OMNI_raw = requests.get(siteURL, stream=True); #get the data, stream is key I think
                # OMNI_raw.encoding = OMNI_raw.apparent_encoding; #not needed, UTF-8 detected correctly
                OMNI_data = OMNI_raw.text; #decode to text
            #END TRY
            
            OMNI_data_split = OMNI_data.split('\n'); #split by the \n line seperator
            OMNI_data_raw = np.zeros( (len(OMNI_data_split)-1,46) ); #preallocate, 46 total entries for the OMNI data
            for j in range(0,len(OMNI_data_split)-1): #-1 cause last split is a blank line
                OMNI_data_raw[j,:] = np.float64(OMNI_data_split[j].split()); #make it numbers n such
            #END FOR j
            
            #use the with thing because it prevents the file being "open in python" still in the event of a crash THERE IS NO WAY TO CLOSE A FILE?
            with h5py.File(OMNI_dataFilePath, 'w') as OMNIfile:
                OMNI_data_raw_record = OMNIfile.create_dataset("OMNIdata", OMNI_data_raw.shape, dtype='float64', compression='gzip', compression_opts=4, shuffle=True, fletcher32=True); #create dataset for the integers
                OMNI_data_raw_record[...] = OMNI_data_raw; #write that data
            #END WITH
        #END IF
        
        dateRange_full_currentYear = dateRange_dayNum_full[:,0] == dateRange_yearRange[i] ; #get the bits that are for the current year
        
        dateRange_full_allDays = np.isin(OMNI_data_raw[:,1],dateRange_dayNum_full[dateRange_full_currentYear,1]); #get the days of the current year needed
        
        OMNI_output_raw[cntr:cntr+(np.sum(dateRange_full_currentYear)*1440),:] = OMNI_data_raw[dateRange_full_allDays,:]; #copy in all relevant data
        
        cntr = cntr + np.sum(dateRange_full_currentYear)*1440; #increment the counter    
        FLG_dataAlreadyDownloaded = 0; #reset the flag
        
        #---very fast and lazy PCS import code, no way to make it automatic since---
        try:
            from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
            dateRange_full = subfun_dayNum_to_date(dateRange_dayNum_full); #can't use date dict here because date range may be artificially extended
            #try to import PCS if the PCS file is there from - wrapped in a try b/c even if it fails idk
            OMNI_PCFilePath = settings_paths['data'] + '\\OMNI\\pcnpcs_min' + str(dateRange_yearRange[i]) + '.txt'; #create the expected data path
            if( os.path.isfile(OMNI_PCFilePath) == 1 ): #check if omni data file exists
                import pandas as pd      
                PC_read = pd.read_csv(OMNI_PCFilePath, delim_whitespace=True, lineterminator='\n'); #read in the data
                
                dateRange_full_currentYear_where = np.where(dateRange_full_currentYear)[0]; #GET WHERE
                for j in range(0,dateRange_full_currentYear_where.size):
                    PC_read_currDay = PC_read.iloc[:,0] == str(dateRange_full[dateRange_full_currentYear_where[j],0]).zfill(4)+'-'+str(dateRange_full[dateRange_full_currentYear_where[j],1]).zfill(2)+'-'+str(dateRange_full[dateRange_full_currentYear_where[j],2]).zfill(2);
                    PCS_temp = PC_read['PCS'][PC_read_currDay].to_numpy(); #get the data we need
                    PCS_temp[PCS_temp >= 999] = np.nan; #nan the silly numbers
                    if( 'PCS' not in locals() ):
                        PCS = PCS_temp.copy(); #copy it over
                    else:
                        PCS = np.append(PCS,PCS_temp); #tack it on
                    #END IF
                #END FOR j
            #END IF
        except:
            pass; #whatever
        #END TRY
    #END FOR i
    

    
    #Weed out bad numbers and such (denoted by varying amounts of 9's)
    #Columns 7, 8, and 9 seem to be solid markers for data being bad - not sure for sure but doin it
    # turns out data can be good too with the bads (makes sense, multiple sources) - so remove 999 columns manually
    # OMNI_output_badData = (OMNI_output_raw[:,6] == 999) | (OMNI_output_raw[:,7] == 999) | (OMNI_output_raw[:,8] == 999) | (OMNI_output_raw[:,44] == 999.99); #find the bad data, it's denoted by 999
    # OMNI_output_raw = OMNI_output_raw[np.logical_not(OMNI_output_badData),:]; #remove the bad data
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,4] , 99),4] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,5] , 99),5] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,6] , 999),6] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,7] , 999),7] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,8] , 999),8] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,9] , 999999),9] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,10] , 999999),10] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,11] , 99.99),11] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,12] , 999999),12] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,13] , 9999.99),13] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,14] , 9999.99),14] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,15] , 9999.99),15] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,16] , 9999.99),16] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,17] , 9999.99),17] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,18] , 9999.99),18] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,19] , 9999.99),19] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,20] , 9999.99),20] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,21] , 99999.9),21] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,22] , 99999.9),22] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,23] , 99999.9),23] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,24] , 99999.9),24] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,25] , 999.99),25] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,26] , 9999999),26] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,27] , 99.99),27] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,28] , 999.99),28] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,29] , 999.99),29] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,30] , 999.9),30] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,31] , 9999.99),31] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,32] , 9999.99),32] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,33] , 9999.99),33] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,34] , 9999.99),34] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,35] , 9999.99),35] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,36] , 9999.99),36] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,44] , 999.99),44] = np.nan;
    OMNI_output_raw[np.isclose(OMNI_output_raw[:,45] , 99.9),45] = np.nan;
    
    OMNI_output = OMNI_output_raw; #get the data desired only, call it a day
    # #--- this calc for IMF clock angle includes the X axis by making an effective 2 axis system by combining the X and Y planes while aligning the orientation to the Y plane (so you can get it to go 360 - w/o orientation fixes it can only go to 180 b/c mag is always+) ---
    # mag_Bt = np.sqrt(OMNI_output[:,14]**2+OMNI_output[:,17]**2+OMNI_output[:,18]**2); #magnitude of IMF B vector
    # angle_BzToVect = np.arcsin(OMNI_output[:,18]/mag_Bt)*180/np.pi; #deg, calc the XYplane-to-Vect angle
    # k = (((OMNI_output[:,17]>=0)-(OMNI_output[:,17]<0))>0); #logical
    # angle_BzToVect[k] = 90-angle_BzToVect[k]; #deg, align to the Y plane & adjust to be Bz-to-Vect angle
    # angle_BzToVect[~k] += 270; #deg, align to the Y plane & adjust to be Bz-to-Vect angle
    # OMNI_output = np.append(OMNI_output,np.hstack(angle_BzToVect,mag_Bt)); #append on
    #--- this calc for IMF clock angle is just the Y and Z components with 0 at Bz+, it ignores X completely (https://www.sws.bom.gov.au/Category/Solar/Solar%20Conditions/Solar%20Wind%20Clock%20Angle/Solar%20Wind%20Clock%20Angle.php does this so doing it)---
    mag_Bclock = np.sqrt(OMNI_output[:,17]**2+OMNI_output[:,18]**2); #magnitude of IMF B vector
    # angle_BzToVect = np.arctan2(OMNI_output[:,18],OMNI_output[:,17])*180/np.pi; #deg, calc the Yplane-to-Vect angle (wrt Y plane) (do 90-... to flip around -180 and 180 AND add a 90 deg offset)
    # k = (90-angle_BzToVect)>=0; #logical
    # angle_BzToVect[k] = 90-angle_BzToVect[k]; #deg, align to the Y plane & adjust to be Bz-to-Vect angle
    # angle_BzToVect[~k] = 450-angle_BzToVect[~k]; #deg, align to the Y plane & adjust to be Bz-to-Vect angle
    # angle_BzToVect[angle_BzToVect > 180] += -360; #deg, convert from 0to360 to -180to180
    # angle_BzToVect = -angle_BzToVect; #flip it around -180 to 180
    # Per https://supermag.jhuapl.edu/info/data.php?page=swdata the clock angle may be = atan2(By,Bz)
    #--- per Akasofu, S.I. Energy coupling between the solar wind and the magnetosphere. Space Sci Rev 28, 121–190 (1981). https://doi.org/10.1007/BF00218810
    angle_BzToVect = np.arctan(np.abs(OMNI_output[:,17])/np.abs(OMNI_output[:,18]))*180/np.pi; #deg, calc per the paper
    angle_BzToVect[OMNI_output[:,18] < 0] = 180 - angle_BzToVect[OMNI_output[:,18] < 0];
    
    symh_diff = np.insert(np.diff(OMNI_output[:,41]),0,0); #change in SYM/H
    
    epsilonParam = ((4*np.pi/(1.25663706212*10**(-6)))*OMNI_output[:,21]*1000*((np.sqrt(OMNI_output[:,14]**2+OMNI_output[:,17]**2+OMNI_output[:,18]**2)*10**-9)**2)*(np.sin((angle_BzToVect*np.pi/180)/2)**4)*(7*6371*1000)**2)/10**9; #GW, epsilon parameter per Akasofu, S.I. Energy coupling between the solar wind and the magnetosphere. Space Sci Rev 28, 121–190 (1981). https://doi.org/10.1007/BF00218810
    
    try:
        OMNI_output = np.append(OMNI_output,np.vstack((angle_BzToVect,mag_Bclock,symh_diff,epsilonParam,PCS)).T,axis=1); #append on
    except:
        OMNI_output = np.append(OMNI_output,np.vstack((angle_BzToVect,mag_Bclock,symh_diff,epsilonParam)).T,axis=1); #append on
        OMNI_dict.pop('PC(S)', None); #remove PCS from the lists
        OMNI_dictPlot.pop(50, None);
    #END TRY
    return OMNI_output, OMNI_dict, OMNI_dictPlot #finish it out