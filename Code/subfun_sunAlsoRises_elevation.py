#GOAL: Get sunrise and sunset times
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: sunrise suinset
#based on NOAA example https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#I implemented all bits cause I can, not all are used but they seemed like they might be useful at some point idk
#Note that timeIndexes can really be any form of time in seconds associated with the days (prolly)

import numpy as np #import in here I dunno
import timezonefinder
from datetime import datetime, timedelta
import pytz
from subfun_strstr import strstr
from subfun_dayNum_to_date import subfun_dayNum_to_date
from urllib.request import urlopen


def sunAlsoRises_elevation(dateRange_full,timeIndexes=None,dataRate=None,lat=None,long=None,timeZone=None):
    
    if( dateRange_full.shape[1] == 2 ):
        dateRange_full = subfun_dayNum_to_date(dateRange_full); #convert to dateRange_full from dateRange_dayNum_full that was actually sent in
    #END IF
    
    if( np.any((lat == None) | (long == None)) == False ):
        #--- Find time zone of lat and long ---
        tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
        dayNite_timeZoneID = tf.certain_timezone_at(lat=lat, lng=long); #use it to find the time zone
        if dayNite_timeZoneID is None:
            #use geonames site as a backup
            url = 'http://api.geonames.org/timezone?lat='+str(lat)+'&lng='+str(long)+'&username=razzluhdzuul'; #create link for lat/long
            webpage = urlopen(url).read(); #get the raw HTML and read it
            try:
                charset = webpage.headers.get_content_charset(); #get the charset from the page holder, w/o it doesn't work
                if( charset is None ):
                    charset = 'utf-8'; #assume utf-8
                #END IF
            except:
                charset = 'utf-8'; #assume utf-8
            #END TRY
            webpage = webpage.decode(charset); #"decode" the HTML content so it's legible
            # index_start = strstr(webpage,'<dstOffset>')[0]; #get where dstOffset is
            # index_end = strstr(webpage,'</dstOffset>')[0]; #get where dstOffset is
            # dayNite_DSToffset_str = webpage[index_start+11:index_end]; #get dst offset
            # dayNite_DSToffset = np.float64(dayNite_DSToffset_str); #and convert to number
            # index_start = strstr(webpage,'<gmtOffset>')[0]; #get where UT offset is
            # index_end = strstr(webpage,'</gmtOffset>')[0]; #get where UT offset is
            # dayNite_UToffset_str = webpage[index_start+11:index_end]; #get UT offset
            # dayNite_UToffset = np.float64(dayNite_UToffset_str); #and convert to number
            index_start = strstr(webpage,'<timezoneId>')[0]; #get where time zone ID is
            index_end = strstr(webpage,'</timezoneId>')[0]; #get where time zone ID is
            dayNite_timeZoneID = webpage[index_start+12:index_end]; #get time zone ID
        #END IF
        
        timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
        # timeZoneObj_UTC = pytz.timezone('UTC'); #make a timezone object
        
        #--- Calc time zones ---
        tzString = []; #prep a list
        tzInt = np.empty(dateRange_full.shape[0],dtype=np.int64); #preallocate
        for j in range(0,dateRange_full.shape[0]):
            tzString.append(datetime.strptime(str(dateRange_full[0,:]), '[%Y\t%m\t%d]').astimezone(timeZoneObj).strftime('%z').rstrip('0')); #get the timezone string
            if( tzString[j] != '+' ):
                tzInt[j] = np.int64(tzString[j]); #record as a number
            else:
                tzInt[j] = 0; #record as a number, UTC work-around
            #END IF
        #END FOR j
    elif( timeZone != None ):
        dayNite_timeZoneID = timeZone;
        timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
        # timeZoneObj_UTC = pytz.timezone('UTC'); #make a timezone object
        
        #--- Calc time zones ---
        tzString = []; #prep a list
        tzInt = np.empty(dateRange_full.shape[0],dtype=np.int64); #preallocate
        for j in range(0,dateRange_full.shape[0]):
            tzString.append(datetime.strptime(str(dateRange_full[0,:]), '[%Y\t%m\t%d]').astimezone(timeZoneObj).strftime('%z').rstrip('0')); #get the timezone string
            if( tzString[j] != '+' ):
                tzInt[j] = np.int64(tzString[j]); #record as a number
            else:
                tzInt[j] = 0; #record as a number, UTC work-around
            #END IF
        #END FOR j
    else:
        dayNite_timeZoneID = 'UTC';
        # timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
        timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object that's UTC b/c no location info provided
        
        #--- Calc time zones ---
        tzString = []; #prep a list
        tzInt = np.empty(dateRange_full.shape[0],dtype=np.int64); #preallocate
        for j in range(0,dateRange_full.shape[0]):
            tzString.append(datetime.strptime(str(dateRange_full[0,:]), '[%Y\t%m\t%d]').astimezone(timeZoneObj).strftime('%z').rstrip('0')); #get the timezone string
            if( tzString[j] != '+' ):
                tzInt[j] = np.int64(tzString[j]); #record as a number
            else:
                tzInt[j] = 0; #record as a number, UTC work-around
            #END IF
        #END FOR j
    #END IF
    
    if( np.any(timeIndexes != None) ):
        if( type(timeIndexes) is list ):
            localTime = timeIndexes; #assume timeIndexes came in just right to be a local time list
            FLG_outType = 'list';
        else:
            #make localTime list happen
            if( np.any(timeIndexes > 86400)):
                timeIndexes = np.mod(timeIndexes,86400); #convert to per day time, assume sec
            #END IF
            localTime = [None for j in range(0,dateRange_full.shape[0])]; #prep an array
            timeIndexes_splitz = np.append(np.insert(np.where(np.diff(timeIndexes) < 0)[0]+1,0,0),timeIndexes.size); #get where splitz occur
            for j in range(0,dateRange_full.shape[0]):
                localTime[j] = timeIndexes[timeIndexes_splitz[j]:timeIndexes_splitz[j+1]]/3600; #split it up, assume sec -> hr conversion
            #END FOR j
            FLG_outType = 'flat';
        #END IF
    elif( dataRate != None ):
        #calculate some time indexes based on the data rate
        localTime = [np.arange(0,86400,dataRate)/3600 for j in range(0,dateRange_full.shape[0])]; #hr, local time set to noon for sunrise/sunset calcs I guess
        FLG_outType = 'list';
    else:
        print('Error in sunAlsoRises_location: No time indexes (timeIndexes) or data rate (dataRate) provided. Cannot calculate sun locations w/o times. Returning an error.'); #report error
        sys.crash() #better
        return 'Error in sunAlsoRises_location' #very helpful
    #END IF
    
    sunSubSolar_elevation = {'elevation':[None for j in range(0,dateRange_full.shape[0])],
                   'zenith':[None for j in range(0,dateRange_full.shape[0])],
                   }; #prep a dict
    dateRange_full = np.int64(dateRange_full); #convert to int64 for math ez
    for j in range(0,dateRange_full.shape[0]):        
        #--- Julian Day from YR/M/D ---
        JD = 367*dateRange_full[j,0] - np.int64(7*(dateRange_full[j,0]+np.int64((dateRange_full[j,1]+9)/12))/4) - \
            np.int64(3*(np.int64((dateRange_full[j,0]+(dateRange_full[j,1]-9)/7)/100)+1)/4) + \
            np.int64(275*dateRange_full[j,1]/9) + dateRange_full[j,2] + 1721028.5; #JD, from https://scienceworld.wolfram.com/astronomy/JulianDate.html
        JD += localTime[j]/24 - np.float64(tzInt[j])/24; #JD, add on the desired local time (12 noon) and adjust that local time by the timezone offset from UTC to give it in a UTC time
        #--- Julian Century ---
        JC = (JD - 2451545)/36525; #JC, it's a julian century or something idk what this is
        #--- Geometric Mean Longitude & Anomoly of Sun ---
        sunLongGeomMean = np.mod(280.46646 + JC*(36000.76983 + JC*0.0003032),360); #deg, 
        sunAnomolyGeomMean = 357.52911 + JC*(35999.05029 - 0.0001537*JC); #deg,
        #--- Earth Orbit Eccentricity ---
        earthEccentricity = 0.016708634 - JC*(0.000042037+0.0000001267*JC); #eccentricity of Earth orbit around sun
        #--- Sun Eq of Ctr? ---
        sunEqOfCtr = np.sin(sunAnomolyGeomMean*np.pi/180)*(1.914602 - JC*(0.004817+0.000014*JC)) + \
            np.sin(2*sunAnomolyGeomMean*np.pi/180)*(0.019993-0.000101*JC) + np.sin(3*sunAnomolyGeomMean*np.pi/180)*0.000289;
        #--- Sun True Longitude & Sun True Anomoly ---
        sunLongTrue = sunLongGeomMean + sunEqOfCtr; #deg, 
        # sunAnomolyTrue = sunAnomolyGeomMean + sunEqOfCtr; #deg, 
        #--- Sun-to-Earth Radius ---
        # sunRadius = (1.000001018*(1 - earthEccentricity*earthEccentricity))/(1 + earthEccentricity*np.cos(sunAnomolyTrue*np.pi/180)); #AU,
        #--- Sun Appogee Longitude ---
        sunLongAppogee = sunLongTrue - 0.00569 - 0.00478*np.sin((125.04 - 1934.136*JC)*np.pi/180); #deg,
        #--- Mean Oblique Elliptic & Oblique Correction & "var y" ---
        obliqueMeanElliptic = 23 + (26 + ((21.448 - JC*(46.815 + JC*(0.00059 - JC*0.001813))))/60)/60; #deg,
        obliqueCorr = obliqueMeanElliptic + 0.00256*np.cos((125.04 - 1934.136*JC)*np.pi/180); #deg, 
        varY = np.tan((obliqueCorr/2)*np.pi/180)*np.tan((obliqueCorr/2)*np.pi/180);
        #--- Sun Right Ascension & Declination (latitude) ---
        # sunRightAscension = np.arctan2( np.cos(obliqueCorr*np.pi/180)*np.sin(sunLongAppogee*np.pi/180), np.cos(sunLongAppogee*np.pi/180) )*180/np.pi; #deg,
        sunDeclination = np.arcsin(np.sin(obliqueCorr*np.pi/180)*np.sin(sunLongAppogee*np.pi/180))*180/np.pi; #deg,
        #--- Eq of Time ---
        eqOfTime = 4*(varY*np.sin(2*sunLongGeomMean*np.pi/180) - 2*earthEccentricity*np.sin(sunAnomolyGeomMean*np.pi/180) + \
            4*earthEccentricity*varY*np.sin(sunAnomolyGeomMean*np.pi/180)*np.cos(2*sunLongGeomMean*np.pi/180) - \
            0.5*varY*varY*np.sin(4*sunLongGeomMean*np.pi/180) - \
            1.25*earthEccentricity*earthEccentricity*np.sin(2*sunAnomolyGeomMean*np.pi/180))*180/np.pi; #min,
        # #--- Sun Sub Solar Point (longitude)  ---
        # sunSubSolar = -15*(localTime[j] - tzInt[j] - 12 + eqOfTime/60); #deg, sub solar point (sun longitude pt where azimuth = 0 deg) [from wikipedia]
        # #--- Hour Angle of Sunrise ---
        # sunHA = (np.arccos(np.cos(90.833*np.pi/180)/(np.cos(lat*np.pi/180)*np.cos(sunDeclination*np.pi/180)) - \
        #     np.tan(lat*np.pi/180)*np.tan(sunDeclination*np.pi/180)))*180/np.pi; #deg,
        # #--- Local Solar Noon ---
        # sunSolarNoon = (720 - 4*long - eqOfTime+np.float64(tzInt[j])*60)/1440; #days, fraction of a day that is local solar noon
        # #--- Local Sunrise Time ---
        # sunRise_local = (sunSolarNoon*1440 - sunHA*4)/1440; #days, fraction of a day that is local sunrise
        # sunRise[j] = sunRise_local - tzInt[j]/24; #days, convert to UTC
        # #--- Local Sunset Time ---
        # sunSet_local = (sunSolarNoon*1440 + sunHA*4)/1440; #days, fraction of a day that is local sunrise
        # sunSet[j] = sunSet_local - tzInt[j]/24; #days, convert to UTC
        # #--- Local Sun Duration ---
        # sunDuration = 8*sunHA; #minutes
        #--- True Solar Time ---
        sunTrueSolarTime = np.mod(localTime/24*1440+eqOfTime+4*long-60*np.float64(tzInt[j]),1440); #min,
        #--- Hour Angle ---
        if( sunTrueSolarTime/4 < 0 ):
            HA = sunTrueSolarTime/4 + 180; #deg,
        else:
            HA = sunTrueSolarTime/4 - 180; #deg,
        #END IF
        #--- Solar Zenith Angle ---
        sunZenith = np.arccos(np.sin(lat*np.pi/180)*np.sin(sunDeclination*np.pi/180)+np.cos(lat*np.pi/180)*np.cos(sunDeclination*np.pi/180)*np.cos(HA*np.pi/180))*180/np.pi; #deg, 
        sunSubSolar_elevation['zenith'][j] = sunZenith; #record
        #--- Solar Elevation Angle ---
        sunElev = 90 - sunZenith; #deg, straight up is 90 deg here while with Zenith straight up is 0 deg
        sunSubSolar_elevation['elevation'][j] = sunElev; #record
        # #--- Approx. Atmospheric Refraction Effect in Deg ---
        # if( sunElev > 85 ):
        #     atmoRef = 0
        # else:
        #     if( sunElev > 5 ):
        #        atmoRef = 58.1/np.tan(sunElev*np.pi/180) - 0.07/np.tan(sunElev*np.pi/180)**3 + 0.000086/np.tan(sunElev*np.pi/180)**5;
        #     else:
        #         if( sunElev > -0.575 ):
        #             atmoRef = 1735 + sunElev*(-518.2 + sunElev*(103.4 + sunElev*(-12.79 + sunElev*0.711)));
        #         else:
        #             atmoRef = -20.772/np.tan(sunElev*np.pi/180);
        #         #END IF
        #     #END IF
        # #END IF
        # atmoRef = atmoRef/3600; #deg, divide by 3600 for all situations
        # #--- Corrected Solar Elevation Angle ---
        # sunElevCorr = sunElev + atmoRef; #deg, 
        # #--- Solar Azimuth Angle (clockwise from North) ---
        # if( HA > 0 ):
        #     sunAzimuth = np.mod( np.arccos(((np.sin(lat*np.pi/180)*np.cos(sunZenith*np.pi/180)) - \
        #         np.sin(sunDeclination*np.pi/180))/(np.cos(lat*np.pi/180)*np.sin(sunZenith*np.pi/180)))*180/np.pi + 180 , 360); #deg, 
        # else:
        #     sunAzimuth = np.mod( 540 - np.arccos(((np.sin(lat*np.pi/180)*np.cos(sunZenith*np.pi/180)) - \
        #         np.sin(sunDeclination*np.pi/180))/(np.cos(lat*np.pi/180)*np.sin(sunZenith*np.pi/180)))*180/np.pi , 360); #deg, 
        # #END IF
    #END FOR i
    if( FLG_outType == 'flat' ):
        sunSubSolar_elevation['elevation'] = np.concatenate(sunSubSolar_elevation['elevation']).ravel(); #flatten into an array
        sunSubSolar_elevation['zenith'] = np.concatenate(sunSubSolar_elevation['zenith']).ravel(); #flatten into an array
    #END IF
    
    return sunSubSolar_elevation
#END DEF
