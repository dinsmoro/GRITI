#GOAL: Get sunrise and sunset times
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: sunrise suinset
#based on NOAA example https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#I implemented all bits cause I can, not all are used but they seemed like they might be useful at some point idk

import numpy as np #import in here I dunno
import timezonefinder
from datetime import datetime, timedelta
import pytz
from Code.subfun_strstr import strstr
from urllib.request import urlopen

def sunAlsoRises(dateRange_full,lat,long,FLG_excess=False):
    
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
    
    #--- Pad dateRange_full to account for the timezone offset (a day/night change might happen the next local day, so pad for it) ---
    # if( tzString[0][0] == '+' ):
    tempTime = datetime.strptime(str(dateRange_full[-1,:]), '[%Y\t%m\t%d]') + timedelta(days=1);
    dateRange_fullPad = np.vstack(( dateRange_full, np.array( ( np.int16(tempTime.strftime('%Y')), np.int16(tempTime.strftime('%m')), np.int16(tempTime.strftime('%d'))) ) ));
    #END IF
    # if( tzString[-1][0] == '-' ):
    tempTime = datetime.strptime(str(dateRange_full[0,:]), '[%Y\t%m\t%d]') - timedelta(days=1);
    dateRange_fullPad = np.vstack(( np.array( ( np.int16(tempTime.strftime('%Y')), np.int16(tempTime.strftime('%m')), np.int16(tempTime.strftime('%d'))) ), dateRange_fullPad ));
    #END IF
    
    #--- Calc time zones again after updating the date range ---
    tzString = []; #prep a list
    tzInt = np.empty(dateRange_fullPad.shape[0],dtype=np.int64); #preallocate
    for j in range(0,dateRange_fullPad.shape[0]):
        tzString.append(datetime.strptime(str(dateRange_fullPad[0,:]), '[%Y\t%m\t%d]').astimezone(timeZoneObj).strftime('%z').rstrip('0')); #get the timezone string
        if( tzString[j] != '+' ):
            tzInt[j] = np.int64(tzString[j]); #record as a number
        else:
            tzInt[j] = 0; #record as a number, UTC work-around
        #END IF
    #END FOR j
    
    tzInt_unique, tzInt_unique_indexes = np.unique(tzInt,  return_inverse=True); #get unique time zones
    tzInt_unique_currentSiteArray = np.split(np.argsort(tzInt, kind='mergesort'), np.cumsum(np.bincount(tzInt_unique_indexes)[:-1])); #count how many indexes, cumulative sum, assign to sorted index as we go
    
    sunRise = np.empty( dateRange_fullPad.shape[0] , dtype=np.float64); #preallocate
    sunSet = np.empty( dateRange_fullPad.shape[0] , dtype=np.float64); #preallocate
    dateRange_fullPad = np.int64(dateRange_fullPad); #convert to int64 for math ez
    for i in range(0,tzInt_unique.size):
        j = tzInt_unique_currentSiteArray[i]; #get the indexes with the same timezone
        
        localTime = 12; #hr, local time set to noon for sunrise/sunset calcs I guess
        #--- Julian Day from YR/M/D ---
        JD = 367*dateRange_fullPad[j,0] - np.int64(7*(dateRange_fullPad[j,0]+np.int64((dateRange_fullPad[j,1]+9)/12))/4) - \
            np.int64(3*(np.int64((dateRange_fullPad[j,0]+(dateRange_fullPad[j,1]-9)/7)/100)+1)/4) + \
            np.int64(275*dateRange_fullPad[j,1]/9) + dateRange_fullPad[j,2] + 1721028.5; #JD, from https://scienceworld.wolfram.com/astronomy/JulianDate.html
        JD += localTime/24 - np.float64(tzInt[j])/24; #JD, add on the desired local time (12 noon) and adjust that local time by the timezone offset from UTC to give it in a UTC time
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
        #--- Sun Right Ascension & Declination (sun latitude) ---
        # sunRightAscension = np.arctan2( np.cos(obliqueCorr*np.pi/180)*np.sin(sunLongAppogee*np.pi/180), np.cos(sunLongAppogee*np.pi/180) )*180/np.pi; #deg,
        sunDeclination = np.arcsin(np.sin(obliqueCorr*np.pi/180)*np.sin(sunLongAppogee*np.pi/180))*180/np.pi; #deg,
        #--- Eq of Time ---
        eqOfTime = 4*(varY*np.sin(2*sunLongGeomMean*np.pi/180) - 2*earthEccentricity*np.sin(sunAnomolyGeomMean*np.pi/180) + \
            4*earthEccentricity*varY*np.sin(sunAnomolyGeomMean*np.pi/180)*np.cos(2*sunLongGeomMean*np.pi/180) - \
            0.5*varY*varY*np.sin(4*sunLongGeomMean*np.pi/180) - \
            1.25*earthEccentricity*earthEccentricity*np.sin(2*sunAnomolyGeomMean*np.pi/180))*180/np.pi; #min,
        #--- Sun Sub Solar Point (longitude)  ---
        # sunSubSolar = -15*(localTime - tzInt[j] - 12 + eqOfTime/60); #deg, sub solar point (sun longitude pt where azimuth = 0 deg) [from wikipedia]
        #--- Hour Angle of Sunrise ---
        sunHA = (np.arccos(np.cos(90.833*np.pi/180)/(np.cos(lat*np.pi/180)*np.cos(sunDeclination*np.pi/180)) - \
            np.tan(lat*np.pi/180)*np.tan(sunDeclination*np.pi/180)))*180/np.pi; #deg,
        #--- Local Solar Noon ---
        sunSolarNoon = (720 - 4*long - eqOfTime+np.float64(tzInt[j])*60)/1440; #days, fraction of a day that is local solar noon
        #--- Local Sunrise Time ---
        sunRise_local = (sunSolarNoon*1440 - sunHA*4)/1440; #days, fraction of a day that is local sunrise
        sunRise[j] = sunRise_local - tzInt[j]/24; #days, convert to UTC
        #--- Local Sunset Time ---
        sunSet_local = (sunSolarNoon*1440 + sunHA*4)/1440; #days, fraction of a day that is local sunrise
        sunSet[j] = sunSet_local - tzInt[j]/24; #days, convert to UTC
        # #--- Local Sun Duration ---
        # sunDuration = 8*sunHA; #minutes
        # #--- True Solar Time ---
        # sunTrueSolarTime = np.mod(localTime/24*1440+eqOfTime+4*long-60*np.float64(tzInt[j]),1440); #min,
        # #--- Hour Angle ---
        # if( sunTrueSolarTime/4 < 0 ):
        #     HA = sunTrueSolarTime/4 + 180; #deg,
        # else:
        #     HA = sunTrueSolarTime/4 - 180; #deg,
        # #END IF
        # #--- Solar Zenith Angle ---
        # sunZenith = np.arccos(np.sin(lat*np.pi/180)*np.sin(sunDeclination*np.pi/180)+np.cos(lat*np.pi/180)*np.cos(sunDeclination*np.pi/180)*np.cos(HA*np.pi/180))*180/np.pi; #deg, 
        # #--- Solar Elevation Angle ---
        # sunElev = 90 - sunZenith; #deg, straight up is 90 deg here while with Zenith straight up is 0 deg
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
    
    return sunRise, sunSet, dateRange_fullPad
#END DEF
