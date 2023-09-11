"""
Creates movies
Movie types:
#0 = dTEC moving data points (Fastest, confusing to actually see what up)
#1 = dTEC stationary data points
#2 = dTEC stationary data points + Zenith ISR overlay on 2nd plot
#3 = dTEC stationary data points + MISA ISR overlay on 2nd plot
#4 = dTEC stationary data points + AMPERE data on same plot with time average to AMPERE data (every X min)
#5 = dTEC stationary data points + AMPERE data on same plot (no time average for TEC)
#6 = dTEC stationary data points + AMPERE data on same plot (no time average for TEC) + TEC plot with line for current time
#7 = dTEC stationary data points + AMPERE data on same plot (no time average for TEC) + 2 TEC plots with line for current time
#71 = dTEC stationary data points + 2 TEC plots with line for current time
#8 = dTEC stationary data points + OMNI Index of User Choice data on same plot (no time average for TEC)
#9 = dTEC stationary data points + AMPERE data on same plot (no time average for TEC) + OMNI Index of User Choice data on 2nd plot (no time average for TEC)
#10 = dTEC stationary data points + AMPERE data on same plot (no time average for TEC) + AMPERE side plot + TEC side plot
#11 = ONLY AMPERE data
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import ListedColormap
import cartopy
from cartopy.feature.nightshade import Nightshade
import os
import joblib
from time import time
from datetime import datetime
from copy import deepcopy
import pickle
from Code.GRITI_movieMaker_subfun_figMaker import GRITI_movieMaker_subfun_figMaker
from Code.GRITI_movieMaker_subfun_imageWriter import GRITI_movieMaker_subfun_imageWriter
from Code.GRITI_plotHelper_area_init import GRITI_plotHelper_area_init
from Code.GRITI_keo_keogrammer import GRITI_keo_keogrammer
from Code.subfun_timeMatch import subfun_timeMatch
from Code.GRITI_movieMaker_subfun_dataGridder import GRITI_movieMaker_subfun_dataGridder
from Code.GRITI_movieMaker_subfun_dayniteCalc import GRITI_movieMaker_subfun_dayniteCalc
from Code.subfun_textNice import textNice
from Code.subfun_figFitter import figFitter
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.subfun_monthNum_to_word import subfun_monthNum_to_word
        

def GRITI_movieMaker(data, dates, settings, FLG_parallelizer=True):

    movie_uses_TEC = [0,1,2,3,4,5,6,7,71,8,9,10];
    movie_uses_TEC_gridded = [1,2,3,4,5,6,7,71,8,9,10];
    movie_uses_TEC_keo = [6, 7, 71, 10];
    
    movie_uses_AMPERE = [4,5,6,7,9,10,11];
    movie_uses_AMPERE_gridded = [4, 5, 6, 7, 9, 10];
    
    movie_uses_OMNI = [8, 9];
    
    print('\nAnalysis: Starting Movie Making\n');
    tic = time(); #time some stuff
    # plt.ioff(); #disable showing the plot cause no point, can disable this to check it
    
    #----------------------PREP VALUES NEEDED----------------------
    time4mag = datetime(dates['date range zero hr'][0],dates['date range zero hr'][1],dates['date range zero hr'][2]); #date time object for aacgmv2   

    #----------------------CALCULATE SUNRISE SUNSET TIMES----------------------
    #SECOND STEP: CALC LOCAL SUNRISE/SUNSET TIMES
    #calcs done in UT/GMT
    # based on calc steps in https://www.mathworks.com/examples/matlab/community/21093-estimating-sunrise-and-sunset
    #Preallocate
    dayNite_Grid_size = 1000; #number to subdivide grid into
    dayNite_sunrise = np.zeros( (dayNite_Grid_size,dayNite_Grid_size,dates['date range full dayNum'].shape[0]) ,dtype=np.float64); #hr UT, prep sunrise time for each day in the lat/long grid
    dayNite_sunset = np.zeros( (dayNite_Grid_size,dayNite_Grid_size,dates['date range full dayNum'].shape[0]) ,dtype=np.float64); #hr UT, prep sunset time for each day in the lat/long grid
    [dayNite_Grid_Long,dayNite_Grid_Lat] = np.meshgrid(np.linspace(np.min(settings['map']['long range']),np.max(settings['map']['long range']),dayNite_Grid_size),np.linspace(np.min(settings['map']['lat range']),np.max(settings['map']['lat range']),dayNite_Grid_size)); #degc, make two matrixes that have all the corresponding points in a lat/long grid
    dayNite_long_corrected = 4*(dayNite_Grid_Long); #calc corrected longitude, for sunrise/sunset time
    
    for i in range(0,dates['date range full dayNum'].shape[0]):
    
        dayNite_B = 360*(dates['date range full dayNum'][i,1] - 81)/365*np.pi/180; #rad, some sort of angle based on days and stuff
        dayNite_EoT_corrected = 9.87*np.sin(2*dayNite_B) - 7.53*np.cos(dayNite_B) - 1.5*np.sin(dayNite_B); #eq for Time Correction
        dayNite_solar_corrected = dayNite_long_corrected + np.tile(dayNite_EoT_corrected,(dayNite_Grid_size,dayNite_Grid_size) ); #min, solar time correction - for noon
    
        dayNite_solar_declination = np.arcsin(np.sin(23.45*np.pi/180)*np.sin(360*(dates['date range full dayNum'][i,1] - 81)/365*np.pi/180)); #rad, solar declination
    
        dayNite_temp = -np.tan(dayNite_Grid_Lat*np.pi/180)*np.tan(dayNite_solar_declination); #calc some mid step
    #     dayNite_temp( dayNite_temp >= 1 ) = dayNite_temp( dayNite_temp >= 1 ) - -2*(1 - dayNite_temp( dayNite_temp >= 1 )); #attempt a flip
    #     dayNite_temp( dayNite_temp <= -1) = dayNite_temp( dayNite_temp <= -1) + -2*(1 + dayNite_temp( dayNite_temp <= -1) ); #attempt a flip
        k = (dayNite_temp[:,0] <= -1) | (dayNite_temp[:,0] >= 1); #prep to replace these
    #     kF = find(k == 1,1,'first');
        
        dayNite_sunrise[:,:,i] = 12 - np.real(np.arccos(dayNite_temp)*180/np.pi)/15 - dayNite_solar_corrected/60; #hr UT, sunrise time
    #     for(j = 1:dayNite_Grid_size) %interpolation didn't work - too precipitous
    #         dayNite_sunrise(k,j,i) = interp1(1:1:(kF-1),dayNite_sunrise(~k,j,i),kF:1:dayNite_Grid_size,'spline','extrap'); %interp for each long set (calc oofs out at specific lat)
    #     end
        dayNite_sunset[:,:,i] = 12 + np.real(np.arccos(dayNite_temp)*180/np.pi)/15 - dayNite_solar_corrected/60; #hr UT, sunrise time
        #Accurate to like 10 minutes or whatever (breaks at high latitudes...)
        
        dayNite_sunrise[k,:,i] = np.nan; #remove data that we can't calc with this alg
        dayNite_sunset[k,:,i] = np.nan; #remove data that we can't calc with this alg
        
        dayNite_sunrise[:,:,i] = dayNite_sunrise[:,:,i] + dates['date range zero hr hours'][i]; #adjust to make it align to the hourly schedule
        dayNite_sunset[:,:,i] = dayNite_sunset[:,:,i] + dates['date range zero hr hours'][i]; #adjust to make it align to the hourly schedule
    #EMD FOR i
    
    #yo I straight up just transpose and reshape in different ways till it works how I want, matlab's is so much easier
    dayNite_sunrise =  np.reshape(dayNite_sunrise.transpose(2,1,0), (dayNite_Grid_size*dates['date range full dayNum'].shape[0], dayNite_Grid_size) ).T; #reshape into one big thing since each day's 6 AM sunrise or whatever is now per-day (+/-24 etc)
    dayNite_sunset = np.reshape(dayNite_sunset.transpose(2,1,0), (dayNite_Grid_size*dates['date range full dayNum'].shape[0], dayNite_Grid_size) ).T; #reshape into one big thing since each day's 6 AM sunrise or whatever is now per-day (+/-24 etc)
    dayNite_Grid_Lat = np.tile(dayNite_Grid_Lat, (1, dates['date range full dayNum'].shape[0]) ); #copy this to match above 1:1
    dayNite_Grid_Long = np.tile(dayNite_Grid_Long, (1, dates['date range full dayNum'].shape[0]) ); #copy this to match above 1:1
    
    #------------------------CREATE THE FIGURE WE WILL USE---------------------
    #call it just to get the offsets so the dynamic fitter function isn't needed every time
    movie_dict = GRITI_movieMaker_subfun_figMaker(data, settings, time4mag, movie_figOffsets=None);
    # plt.close(movie_dict['fig']); #close the figure, just needed its offsets
        
    plt.rcParams['animation.ffmpeg_path'] = settings['paths']['ffmpeg']; #sets where the ffmpeg path is
    
    #---------------------FORCE GRID TO HAVE SQUARE GRID POINTS----------------
    #I don't know if it's a great idea. But it works
#    gif_Grid_TotalGrids = (settings['map']['lat range'][-1] - settings['map']['lat range'][0])*(settings['map']['long range'][-1] - settings['map']['long range'][0])*settings['movie']['grid divider']; #get an estimate for how many pixels we want   
    movie_grid_latSpaces = (settings['map']['lat range'][-1] - settings['map']['lat range'][0])*np.sqrt(settings['movie']['grid divider']); #get the lat degrees
    movie_grid_longSpaces = (settings['map']['long range'][-1] - settings['map']['long range'][0])*np.sqrt(settings['movie']['grid divider']); #get the long degrees
    if( 'stere' in settings['map']['projection name'] ):
        #polar plot so the square stuff won't work here, but we do know the shape has to be a circle (no oval!)
        #latitude is radius, longitude covers whole area so doesn't matter
        movie_grid_longSpaces = 2*np.pi*movie_grid_latSpaces; #circumnference
        movie_grid_longSpaces = np.int64(np.round(movie_grid_longSpaces)); #make gif grid total be about the right size
        movie_grid_latSpaces = np.int64(np.round(movie_grid_latSpaces)); #make gif grid total be about the right size
    else:
        gif_Grid_AdjustmentRatio = settings['movie']['fig size'][0]/settings['movie']['fig size'][1]; #get the adjustment ratio needed
        gif_Grid_TotalGrids = movie_grid_latSpaces*movie_grid_longSpaces; #get total grid points
        if( movie_grid_longSpaces > movie_grid_latSpaces):
            #adjust long spaces since it's the bigger number
            movie_grid_longSpaces = movie_grid_latSpaces*gif_Grid_AdjustmentRatio; #make it square (make it be lat spaces) then apply the ratio
            gif_Grid_TotalGrids_Adj = movie_grid_latSpaces*movie_grid_longSpaces; #get the adjusted total
            gif_Grid_TotalGrids_RatioSqrt = np.sqrt(gif_Grid_TotalGrids/gif_Grid_TotalGrids_Adj); #get the ratio of total grids, take sqrt b/c of how we'll adjust lat/long spaces
            movie_grid_longSpaces = np.int64(np.round(movie_grid_longSpaces*gif_Grid_TotalGrids_RatioSqrt)); #make gif grid total be about the right size
            movie_grid_latSpaces = np.int64(np.round(movie_grid_latSpaces*gif_Grid_TotalGrids_RatioSqrt)); #make gif grid total be about the right size
        else:
            #adjust lat spaces since it's the bigger number
            movie_grid_latSpaces = movie_grid_longSpaces/gif_Grid_AdjustmentRatio; #make it square (make it be long spaces) then apply the ratio
            gif_Grid_TotalGrids_Adj = movie_grid_latSpaces*movie_grid_longSpaces; #get the adjusted total
            gif_Grid_TotalGrids_RatioSqrt = np.sqrt(gif_Grid_TotalGrids/gif_Grid_TotalGrids_Adj); #get the ratio of total grids, take sqrt b/c of how we'll adjust lat/long spaces
            movie_grid_longSpaces = np.int64(np.round(movie_grid_longSpaces*gif_Grid_TotalGrids_RatioSqrt)); #make gif grid total be about the right size
            movie_grid_latSpaces = np.int64(np.round(movie_grid_latSpaces*gif_Grid_TotalGrids_RatioSqrt)); #make gif grid total be about the right size
        #END IF
    #END IF
    
    #----------------------PRIME THE LAT/LONG GRID-----------------------------
    movie_grid_lat = np.linspace(np.min(settings['map']['lat range']),np.max(settings['map']['lat range']),movie_grid_latSpaces+1); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
    movie_grid_long = np.linspace(np.min(settings['map']['long range']),np.max(settings['map']['long range']),movie_grid_longSpaces+1); #degc, create long points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
    movie_grid_latDelta = np.abs(movie_grid_lat[1] - movie_grid_lat[0]); #degc, lat delta
    movie_grid_longDelta = np.abs(movie_grid_long[1] - movie_grid_long[0]); #degc, long delta
    settings['movie']['grid lat'] = movie_grid_lat;
    settings['movie']['grid long'] = movie_grid_long;
    settings['movie']['grid lat delta'] = movie_grid_latDelta;
    settings['movie']['grid long delta'] = movie_grid_longDelta;
    settings['movie']['grid lat spaces'] = movie_grid_latSpaces
    settings['movie']['grid long spaces'] = movie_grid_longSpaces
    
    #---------------------PRIME THE TIME LIMIT RANGE---------------------------
    if( settings['movie']['movie type'] != 11 ):
        if( settings['movie']['time lim'] == 1 ):
            movie_timeRange = np.zeros(2,dtype=np.int64); #preallocate
            movie_timeRange[0] = np.where( np.abs(np.min(settings['movie']['time lim range']) - (data['TEC']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600) == np.min(np.abs(np.min(settings['movie']['time lim range']) - (data['TEC']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600)) )[0]; #get index where min is closest for the time step
            movie_timeRange[1] = np.where( np.abs(np.max(settings['movie']['time lim range']) - (data['TEC']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600) == np.min(np.abs(np.max(settings['movie']['time lim range']) - (data['TEC']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600)) )[0]; #get index where max is closest for the time step
        else:
            movie_timeRange = np.array([0,data['TEC']['time unique'].size]); #set to max time if settings['movie']['time lim'] isn't on
        #END IF
    else:
        if( settings['movie']['time lim'] == 1 ):
            movie_timeRange = np.zeros(2,dtype=np.int64); #preallocate
            movie_timeRange[0] = np.where( np.abs(np.min(settings['movie']['time lim range']) - (data['AMPERE']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600) == np.min(np.abs(np.min(settings['movie']['time lim range']) - (data['AMPERE']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600)) )[0]; #get index where min is closest for the time step
            movie_timeRange[1] = np.where( np.abs(np.max(settings['movie']['time lim range']) - (data['AMPERE']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600) == np.min(np.abs(np.max(settings['movie']['time lim range']) - (data['AMPERE']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600)) )[0]; #get index where max is closest for the time step
        else:
            movie_timeRange = np.array([0,data['AMPERE']['time unique'].size]); #set to max time if settings['movie']['time lim'] isn't on
        #END IF
    #END IF
    
    #---------------------PRIME THE POLAR CIRCULAR STUFF-----------------------
    # if( 'stere' in settings['map']['projection name']):
    #     #Constants needed
    #     ReAdj = Re*1000*np.pi/180; #convert Re to meters, toss in the pi/180 conversion needed for degrees as well
    
    #     #THIS IS FOR ROTATING THE EARTH (from ref J2000)
    #     #J2000 is
    #     #January 1, 2000, 11:58:55.816
    #     #Good luck
    #     dateRange_J2000 = np.zeros( (2,2) ); #preallocate
    #     dateRange_J2000[0,0] = 2000; #put in J2000's defining stuff
    #     dateRange_J2000[0,1] = 1; #put in J2000's defining stuff
    #     dateRange_J2000[1,0] = dates['date range zero hr dayNum'][0]; #get the year number for the reference day
    #     dateRange_J2000[1,1] = dates['date range zero hr dayNum'][1]; #get the day number for the reference day
    #     [_,dateRange_J2000] = subfun_dateORdayNum_to_fullRange(dateRange_J2000); #get the days involved
    #     secRange_J2000 = 3600*24*(dateRange_J2000.shape[0]-1) - (11*3600+58*60+55.816); #sec since J2000 epoch of Jan 1st, 2000 @ 11:58:55.816 from https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)#J2000.0
    #     # used full hour since that's what days do, made sense promise
    #     # Mean angular velocity of Earth (7.292115*10^-5) is from http://hpiers.obspm.fr/eop-pc/models/constants.html    
    #     radRange_J2000_EmpiricalConst = -156.42; #deg, found from seeing what it was and what it was supposed to be in stellarium
    #     if( settings['movie']['movie type'] != 11 ):
    #         radRange_J2000 = (secRange_J2000 + (data['TEC']['time unique'][movie_timeRange[0]:movie_timeRange[1]] - dates['date range zero hr dayNum'][1]*86400))*-7.292115*10**-5; #rad, converts time since J2000 to radians covered since J2000 through earth's angular velocity
    #     else:
    #         radRange_J2000 = (secRange_J2000 + (data['AMPERE']['time unique'][movie_timeRange[0]:movie_timeRange[1]] - dates['date range zero hr dayNum'][1]*86400))*-7.292115*10**-5; #rad, converts time since J2000 to radians covered since J2000 through earth's angular velocity
    #     #END IF
    #     radRange_J2000 = radRange_J2000 + radRange_J2000_EmpiricalConst*np.pi/180; #adjustment
    #     # i = find( min(abs((timeUnique - dates['date range zero hr'](2))*24 - 10)) == abs((timeUnique - dates['date range zero hr'](2))*24 - 10) );
    #     # radRange_J2000(i)*180/pi
    #     #Reports 84.93 deg, should be settings['map']['site coords'][0,1] -71.49 deg from using stellarium at 13:00 EDT (17 UT) where sun is directly overhead (more or less) - add -156.42 deg to fix
    #     #more accurate would be to relocate to equator I guess
    #     #Couldn't find where earth was pointing with J2000 with reference to the sun easily so yolo
    #     radRange_J2000 = np.mod(radRange_J2000,2*np.pi) - np.pi; #converts to -180 to 180
    #     # mod brings it down to 0<=angle<=2*pi b/c cyclic trig
    #     # - pi brings it down to -pi<=angle<=pi to help with understanding what it is
    #     #***************MAKE THIS NOT EMPIRICAL SOMEDAY******************** CLOSE!
    # #END IF
    if( 'stere' in settings['map']['projection name']):
        #today is that day #***************MAKE THIS NOT EMPIRICAL SOMEDAY******************** CLOSE!
        from Code.subfun_sunAlsoRises_location import sunAlsoRises_location
        import aacgmv2
        if( settings['movie']['movie type'] != 11 ):
            sunSubSolar_loc = np.unique(data['TEC']['time unique']//86400); #reuse var
            sunSubSolar_loc = np.vstack( (np.ones(sunSubSolar_loc.size, dtype=np.int32)*dates['date range zero hr dayNum'][0], sunSubSolar_loc) ).T; #stack em (assumes no year change, will need revamp)
            sunSubSolar_loc = sunAlsoRises_location(dates['date range full dayNum'],timeIndexes=data['TEC']['time unique'],timeZone='UTC');
            aacgmv2_alt = np.median(data['TEC']['pierceAlt']);
        else:
            sunSubSolar_loc = np.unique(data['AMPERE']['time unique']//86400); #reuse var
            sunSubSolar_loc = np.vstack( (np.ones(sunSubSolar_loc.size, dtype=np.int32)*dates['date range zero hr dayNum'][0], sunSubSolar_loc) ).T; #stack em (assumes no year change, will need revamp)
            sunSubSolar_loc = sunAlsoRises_location(sunSubSolar_loc, timeIndexes=data['AMPERE']['time unique'],timeZone='UTC');
            if( 'altitude' in settings['AMPERE'] ):
                aacgmv2_alt = settings['AMPERE']['altitude']; #use AMPERE altitude
            else:
                aacgmv2_alt = 120.; #default, great for auroral zone stuff (like field aligned currents)
            #END IF
        #END IF        
        if( settings['map']['coord type'] == 'mag' ):
            if( settings['movie']['movie type'] != 11 ):
                tempTimesSec = data['TEC']['time unique']; #get it as a var
            else:
                tempTimesSec = data['AMPERE']['time unique']; #get it as a var
            #END IF
            time4mag_hr = np.int32(np.mod(tempTimesSec, 86400)//3600); #get hours
            time4mag_min = np.int32(np.mod(tempTimesSec, 86400)//60-time4mag_hr*60); #get the minutes
            time4mag_sec = np.int32(np.mod(tempTimesSec, 86400)-time4mag_min*60-time4mag_hr*3600); #get the seconds
            sunSubSolar_loc['lat geo'] = np.copy(sunSubSolar_loc['lat']);
            sunSubSolar_loc['long geo'] = np.copy(sunSubSolar_loc['long']);
            for jj in range(0, tempTimesSec.size): #easier to loop for now
                if( settings['movie']['movie type'] != 11 ):
                    kk = np.where(np.int64(tempTimesSec[jj]/86400) == dates['date range full dayNum'][:,1])[0].item(); #get where the year is gonna be
                    time4mag = datetime(dates['date range full'][kk,0],dates['date range full'][kk,1],dates['date range full'][kk,2], \
                                                 hour = time4mag_hr[jj], minute = time4mag_min[jj], second = time4mag_sec[jj]); #date time object for aacgmv2
                else:
                    kk = np.where(np.int64(tempTimesSec[jj]/86400) == dates['date range full padded dayNum'][:,1])[0].item(); #get where the year is gonna be
                    time4mag = datetime(dates['date range full padded'][kk,0],dates['date range full padded'][kk,1],dates['date range full padded'][kk,2], \
                                                 hour = time4mag_hr[jj], minute = time4mag_min[jj], second = time4mag_sec[jj]); #date time object for aacgmv2
                #END IF
                #---nan result protection---
                incrementor = 0; #increments
                tempLat = np.nan; #set nan to start
                tempLong = np.nan;
                while( np.isnan(tempLat) | np.isnan(tempLong) ):
                    tempAlt = aacgmv2_alt+incrementor*100;
                    # if( tempAlt > 2000 ):
                    #     tempTrace = True;
                    # else:
                    #     tempTrace = False;
                    # #END IF
                    [tempLat, tempLong, _] = aacgmv2.convert_latlon(sunSubSolar_loc['lat'][jj], sunSubSolar_loc['long geo'][jj], tempAlt, time4mag, method_code='G2A|ALLOWTRACE'); #converts from geographic to geomagnetic (AACGMv2)
                    incrementor += 1; #increment
                    if( incrementor > 40 ):
                        break
                    #END IF
                #END WHILE
                sunSubSolar_loc['lat'][jj] = tempLat; #record, this way we avoid NaNs (but get spammed, oh well)
                sunSubSolar_loc['long'][jj] = tempLong;
            #END FOR jj
        #END IF
        radRange_J2000 = (sunSubSolar_loc['long']-90)*np.pi/180; #convert to radians and adjust by 90
    #END IF
    settings['movie']['sun loc'] = sunSubSolar_loc;
    settings['movie']['sun loc']['long rad'] = radRange_J2000;
    
    #----------------------PREP FOR MOVIE MAKING-------------------------------
    if( os.path.isdir(settings['movie']['save locale']) == 0 ): #check if TEC folder exists
        #if not, make it
        os.makedirs(settings['movie']['save locale']);
        print("NOTA BENE: Moving Making Function - Created plot directory: "+settings['movie']['save locale']+"\n");
    #END IF
    os.chdir(settings['movie']['save locale']); #move to the movie save folder
    
    #need to do same frame for a while to make it good
    movie_desiredFPS = np.round(1/settings['movie']['desired frame time']); #frames/sec, desired frames per sec from the frame time (how long a frame is on the screen)
    
    if( settings['movie']['mp4 mode'] == 1 ): #supports choosing video or GIF (video only right now)
        movie_actualFrameRate = 30; #set FPS to 30 to start off
        if( settings['movie']['movie type'] == 0 | settings['movie']['movie type'] == 1 | settings['movie']['movie type'] == 2 | settings['movie']['movie type'] == 3 | settings['movie']['movie type'] == 5 | settings['movie']['movie type'] == 6 | settings['movie']['movie type'] == 7 | settings['movie']['movie type'] == 8 | settings['movie']['movie type'] == 9): #supports doubling frame rate for the movies not using time-averaging (they do their own thing)
            if( ( data['TEC']['time unique'].size/30 > settings['movie']['desired max run time'] ) & ( settings['movie']['disable FPS shift'] == 0 ) ): #try to keep run time down
                movie_actualFrameRate = 60; #set the FPS to 60
            #END IF
        #END IF
    else:
        pass; #nothing special now
    #END IF
    
    if( settings['movie']['movie type'] in [6, 7, 71, 10] ):
        #build keo-specific stuff (add more options as needed)
        settings_keo1 = deepcopy(settings['TEC']['keo']);
        settings_keo1['keo angle'] = settings['movie']['side plots']['keo 1 angle'];
        settings_keo1['keo width orig'] = settings['movie']['side plots']['keo 1 width'];
        settings_keo1['keo N'] = settings['movie']['side plots']['keo 1 N'];
        settings_keo1['keo polar mode'] = settings['movie']['side plots']['keo 1 polar mode'];
        settings_keo1['keo 45 lat or long'] = settings['movie']['side plots']['keo 1 45 lat or long'];
        settings_keo1_map = deepcopy(settings['map']);
        settings_keo1_map['lat range'] = settings['movie']['side plots']['keo 1 lat range'];
        settings_keo1_map['long range'] = settings['movie']['side plots']['keo 1 long range'];
        
        settings_keo2 = deepcopy(settings['TEC']['keo']);
        settings_keo2['keo angle'] = settings['movie']['side plots']['keo 2 angle'];
        settings_keo2['keo width orig'] = settings['movie']['side plots']['keo 2 width'];
        settings_keo2['keo N'] = settings['movie']['side plots']['keo 2 N'];
        settings_keo2['keo polar mode'] = settings['movie']['side plots']['keo 2 polar mode'];
        settings_keo2['keo 45 lat or long'] = settings['movie']['side plots']['keo 2 45 lat or long'];
        settings_keo2_map = deepcopy(settings['map']);
        settings_keo2_map['lat range'] = settings['movie']['side plots']['keo 2 lat range'];
        settings_keo2_map['long range'] = settings['movie']['side plots']['keo 2 long range'];
    #END IF
        
    #prep movie naming
    movie_name = deepcopy(settings['movie']['name base']); #prep it up
    if( settings['movie']['movie type'] in movie_uses_TEC ):
        movie_name += '_ΔvTEC'; #TEC involved
    #END IF
    if( settings['movie']['movie type'] in [2] ):
        movie_name += '_ISR+Zenith';
    #END IF
    if( settings['movie']['movie type'] in [3] ):
        movie_name += '_ISR+MISA';
    #END IF
    if( settings['movie']['movie type'] in movie_uses_AMPERE ):
        movie_name += '_'+settings['AMPERE']['data type'];
    #END IF
    if( settings['movie']['movie type'] in movie_uses_OMNI ):
        movie_name += '_'+settings['OMNI']['data type'].replace('/','-');
    #END IF
    if( settings['movie']['movie type'] in [2,3,4] ):
        movie_name += '_TECtimeAvgd'; #TEC involved
    #END IF
    #same for all gif types
        
    #prime stuff for parallel
    parallel_numThreads = settings['config']['parallel num threads']; #it is CPU bound while drawing, but getting them out is what takes a bit it seems
    batches = np.append(np.arange(movie_timeRange[0],movie_timeRange[1],parallel_numThreads),movie_timeRange[1]);
    
    # precalc what can be
    time_cutout_range_delay_AMPERE_sec = settings['AMPERE']['delay wrt TEC']*3600;
    if( np.isclose(np.mod(time_cutout_range_delay_AMPERE_sec, 1), 0) ):
        time_cutout_range_delay_AMPERE_sec = np.int64(time_cutout_range_delay_AMPERE_sec); #integer it
    #END IF
    
    #frame repeat for proper playback speed
    #Stuff for repeating a frame - needed if want a constant 30 FPS at the desired frame time
    movie_actualFrameRate = 30; #default
    #check if exceed desired max run time at desired frame time
    if( settings['movie']['desired max run time'] < settings['movie']['desired frame time']*(movie_timeRange[1]-movie_timeRange[0]) ):
        #if exceeds, check if 1 frame works
        if( settings['movie']['desired max run time'] < (1/movie_actualFrameRate)*(movie_timeRange[1]-movie_timeRange[0]) ):
            #if exceeds, move to 60 fps
            movie_actualFrameRate = 60; #shift up to 60
            movie_repeatFrameNum = np.int64(np.round(settings['movie']['desired max run time']/(movie_timeRange[1]-movie_timeRange[0])*movie_actualFrameRate)); #how many frames to repeat
        else:
            movie_repeatFrameNum = np.int64(np.round(settings['movie']['desired max run time']/(movie_timeRange[1]-movie_timeRange[0])*movie_actualFrameRate)); #how many frames to repeat
        #END IF
    else:
        movie_repeatFrameNum = np.int64(np.round(settings['movie']['desired frame time']*movie_actualFrameRate)); #how many frames to repeat
    #END IF
    if( movie_repeatFrameNum < 1 ):
        movie_repeatFrameNum = 1; #just in case
    #END IF
    
    # begin movie making
    tic_batch = time();
    #------------------------Initialize Video------------------------------
    FFMpegWriter = manimation.writers['ffmpeg'];
    if( settings['movie']['mp4 mode'] == 1 ):
        movie_writer = FFMpegWriter(fps=movie_actualFrameRate, codec="h264", bitrate=25000); #create a movie writer
        #extra stuff , codec="h264", extra_args=["-y"]
        movie_name = movie_name+".mp4"; #add the file ending
    else: #otherwise it is a gif
        movie_writer = FFMpegWriter(fps=movie_actualFrameRate, codec="gif", bitrate=25000); #create a movie writer
        movie_name = movie_name+".gif"; #add the file ending
    #END IF
    with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
        with joblib.parallel_backend('loky'):
            with joblib.Parallel(n_jobs=parallel_numThreads,pre_dispatch=parallel_numThreads,batch_size=1) as parallel_arbiter:    
                for b in range(0,batches.size-1): #roll through each batch (batched by processor size)
                
                    # Create a list to hold all inputs for a batch
                    parallel_list = []; #prep a list
                    for i in range(batches[b], batches[b+1]):
                        
                        movie_data = {}; #prep holder
                        # TEC data prep
                        if( settings['movie']['movie type'] in movie_uses_TEC ):
                            movie_data['TEC'] = {}; #new dict
                            k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                            movie_data['TEC']['dTEC'] = data['TEC']['dTEC'][k]; #pull out the vTEC now
                            movie_data['TEC']['lat'] = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                            movie_data['TEC']['long'] = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            movie_data['TEC']['time unique'] = data['TEC']['time unique'][i]; #just it
                            movie_data['TEC']['pierceAlt'] = data['TEC']['pierceAlt'];
                        #END IF
                   
                        # AMPERE data prep matched to TEC
                        if( (settings['movie']['movie type'] in movie_uses_AMPERE) & (settings['movie']['movie type'] in movie_uses_TEC) ): 
                            movie_data['AMPERE'] = {}; #new dict
                            movie_data['AMPERE']['data info'] = {}; #new dict
                            if( settings['movie']['use time delays'] == 0 ):
                                k = np.where( (data['AMPERE']['time'] <= data['TEC']['time unique'][i]) & (data['AMPERE']['time'] >= (data['TEC']['time unique'][i]-data['AMPERE']['data rate'])) )[0]; #get where the time point is, make sure it is within the data rate window
                            else:
                                k = np.where( ( (data['AMPERE']['time']+time_cutout_range_delay_AMPERE_sec) <= data['TEC']['time unique'][i]) & ( (data['AMPERE']['time']+time_cutout_range_delay_AMPERE_sec) >= (data['TEC']['time unique'][i]-data['AMPERE']['data rate'])) )[0]; #get where the time point is, make sure it is within the data rate window
                            #END IF
                            movie_data['AMPERE'][settings['AMPERE']['data type']] = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                            movie_data['AMPERE']['lat'] = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                            movie_data['AMPERE']['long'] = data['AMPERE']['long'][k]; #degc, corresponding long values
                            movie_data['AMPERE']['data info']['lat delta'] = data['AMPERE']['data info']['lat delta']; #it's provided
                            movie_data['AMPERE']['data info']['long delta'] = data['AMPERE']['data info']['long delta']; #it's provided
                            
                            if( settings['AMPERE']['data type'] == 'JH' ):
                                #allows for low intensity JH to be ignored (JH can't be less than 0)
                                k = movie_data['AMPERE'][settings['AMPERE']['data type']] > np.min(settings['AMPERE']['plot lim']);
                                movie_data['AMPERE'][settings['AMPERE']['data type']] = movie_data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                                movie_data['AMPERE']['lat'] = movie_data['AMPERE']['lat'][k]; #degc, corresponding lat values
                                movie_data['AMPERE']['long'] = movie_data['AMPERE']['long'][k]; #degc, corresponding long values
                            #END IF
                        elif( (settings['movie']['movie type'] in movie_uses_AMPERE) & (settings['movie']['movie type'] not in movie_uses_TEC) ):
                            k =  np.where( data['AMPERE']['time'] == data['AMPERE']['time unique'][i])[0]; #gets during a time period
                            movie_data['AMPERE'][settings['AMPERE']['data type']] = data['AMPERE'][settings['AMPERE']['data type']][k];
                            movie_data['AMPERE']['lat'] = data['AMPERE']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                            movie_data['AMPERE']['long'] = data['AMPERE']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                            if( settings['AMPERE']['data type'] == 'JH' ):
                                #allows for low intensity JH to be ignored (JH can't be less than 0)
                                k = movie_data['AMPERE'][settings['AMPERE']['data type']] > np.min(settings['AMPERE']['plot lim']);
                                movie_data['AMPERE'][settings['AMPERE']['data type']] = movie_data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                                movie_data['AMPERE']['lat'] = movie_data['AMPERE']['lat'][k]; #degc, corresponding lat values
                                movie_data['AMPERE']['long'] = movie_data['AMPERE']['long'][k]; #degc, corresponding long values
                            #END IF
                        #END IF
                        
                        # General stuff
                        if( settings['movie']['day nite line'] == 2 ):
                            movie_data['movie'] = {}; #new dict
                            movie_data['movie']['sun loc'] = {}; #new dict
                            movie_data['movie']['sun loc']['lat'] = settings['movie']['sun loc']['lat'][i];
                            movie_data['movie']['sun loc']['long'] = settings['movie']['sun loc']['long'][i];
                            movie_data['movie']['sun loc']['long rad'] = settings['movie']['sun loc']['long rad'][i];
                            if( settings['map']['coord type'] == 'mag' ):
                                movie_data['movie']['sun loc']['lat geo'] = settings['movie']['sun loc']['lat geo'][i];
                                movie_data['movie']['sun loc']['long geo'] = settings['movie']['sun loc']['long geo'][i];
                            #END IF
                        #END IF
                        
                        # Load everything into the list
                        parallel_list.append([movie_data, dates, settings, movie_dict['fig offsets'], movie_dict['title offset']]);
                    #END FOR i
                                        
                    #------------------------Draw & Write Images------------------------------
                    # Parallel process the batch
                    if( FLG_parallelizer ):
                        parallel_figs = parallel_arbiter(joblib.delayed(GRITI_movieMaker_subfun_imageWriter)(j, k, l, m, n) for j, k, l, m, n in parallel_list); #will this not destroy the world?
                    else:
                        #misnomer, for testing
                        parallel_figs = []; #prep
                        for j, k, l, m, n in parallel_list:
                            parallel_figs.append(GRITI_movieMaker_subfun_imageWriter(j, k, l, m, n)); #just regular
                        #END FOR lots
                    #END IF                    
                    
                    # Save the figures in the batch
                    for j in range(0,len(parallel_figs)):                        
                        tempFig = pickle.loads(parallel_figs[j]); #change the figure on the fly - legal?
                        tempAxz = tempFig.axes; #get the new axis handle
                        currAxz = movie_dict['fig'].axes; #get the current one too
                        #yolooooo shout out to https://gist.github.com/salotz/8b4542d7fe9ea3e2eacc1a2eef2532c5 for the dark methods
                        currAxz_pos = [None for k in range(0,len(currAxz))]; #preallocate
                        for k in range(0,len(currAxz)):
                            #record positioning info
                            currAxz_pos[k] = deepcopy(currAxz[k].get_position()); #copy out the positioning info in the figure
                            #clean house
                            currAxz[k].remove(); #disconnect from current figure with .remove(), remove is NOT the python remove
                        #END FOR k
                        del currAxz #yeet
                        for k in range(0,len(tempAxz)):
                            tempAxz[k].remove(); #disconnect from temp figure with .remove(), remove is NOT the python remove
                        #END FOR k
                        for k in range(0,len(tempAxz)):
                            #attach to the current fig that I gotta use
                            tempAxz[k].figure = movie_dict['fig']; #attach to the fig I gotta use 
                            movie_dict['fig'].axes.append(tempAxz[k]); #"register" it
                            movie_dict['fig'].add_axes(tempAxz[k]); #"register" it in another way I guess
                            
                            #update the attached axis with the position (not needed since fig was pre-made good)
                            # tempAxz[k].set_position(currAxz_pos[k]); #set the position
                        #END FOR k
                        
                        plt.close(tempFig); #close it again
                          
                        movie_dict['fig'].canvas.draw(); #key for all instances
                        for j in range(0,movie_repeatFrameNum): #runs this as many times as needed to get 30 FPS or whatever
                            movie_writer.grab_frame(); #get the frame and save it
                        #END FOR j
                    #END FOR j
                    
                    if( np.mod(b,batches.size//10) == 0 ):
                        tic_now = time(); #sec, calc current toc
                        toc_now = tic_now - tic;
                        toc_batch = tic_now - tic_batch;
                        batch_2go = batches.size-(b+1);
                        print('At batch '+str(b+1)+' of '+str(batches.size)+' batches.');
                        print('Time spent so far: '+textNice(np.round(toc_now,2))+' sec/ '+textNice(np.round(toc_now/60,2))+' min/ '+textNice(np.round(toc_now/3600,2))+' hrs');
                        print('ETA: '+textNice(np.round((toc_batch/(b+1))*(batch_2go),2))+' sec/ '+textNice(np.round((toc_batch/(b+1))*(batch_2go)/60,2))+' min/ '+textNice(np.round((toc_batch/(b+1))*(batch_2go)/3600,2))+' hrs');
                    #END IF
                #END FOR b
            #END WITH
        #END WITH
        del parallel_list #save some mem
    #END WITH
    
    
    toc = time() - tic; #sec, calc current toc
    print("Time to run: "+textNice(np.round(toc,2))+" sec/ "+textNice(np.round(toc/60,2))+" min /"+textNice(np.round(toc/3600,2))+" hrs");
    plt.ion(); #re-enable plots cause yee
#END DEF