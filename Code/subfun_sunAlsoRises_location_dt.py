#GOAL: Get sunrise and sunset times
#RD on 3/22/19
#
#INPUT: datetime obj, can be in list or not. Aware of timezones but not picky if there aren't any (assumes you did it right with UTC)
#OUTPUT: sunrise suinset
#based on NOAA example https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#I implemented all bits cause I can, not all are used but they seemed like they might be useful at some point idk
#Note that timeIndexes can really be any form of time in seconds associated with the days (prolly)

import numpy as np #import in here I dunno
from datetime import timezone
import pytz


def sunAlsoRises_location_dt(datetime_list):
    
    FLG_inputNotListed = False;
    if( (not isinstance(datetime_list,list)) | (not isinstance(datetime_list,tuple)) ):
        datetime_list = [datetime_list]; #wrap
        FLG_inputNotListed = True;
    #END IF
    
    sunSubSolar_loc = {'lat':[None for j in range(0,len(datetime_list))],
                   'long':[None for j in range(0,len(datetime_list))],
                   }; #prep a dict
    for j in range(0,len(datetime_list)):        
        #--- Ensure datetime obj is in UTC ---
        if( datetime_list[j].tzinfo != None ):
            if( (datetime_list[j].tzinfo != pytz.utc) | (datetime_list[j].tzinfo != timezone.utc) ):
                datetime_list[j] = datetime_list[j].astimezone(timezone.utc); #convert as needed
            #END IF
        #END IF
        
        #--- Julian Day from YR/M/D ---
        dayOfHrs = (datetime_list[j].hour + datetime_list[j].minute/60 + datetime_list[j].second/3600); #hours in the day
        JD = 367*datetime_list[j].year - np.int64(7*(datetime_list[j].year+np.int64((datetime_list[j].month+9)/12))/4) - \
            np.int64(3*(np.int64((datetime_list[j].year+(datetime_list[j].month-9)/7)/100)+1)/4) + \
            np.int64(275*datetime_list[j].month/9) + datetime_list[j].day + 1721028.5; #JD, from https://scienceworld.wolfram.com/astronomy/JulianDate.html
        JD += dayOfHrs/24; #JD, in UTC so we all good - it needs to be in days, so hours/24
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
        sunSubSolar_loc['lat'][j] = sunDeclination; #record
        #--- Eq of Time ---
        eqOfTime = 4*(varY*np.sin(2*sunLongGeomMean*np.pi/180) - 2*earthEccentricity*np.sin(sunAnomolyGeomMean*np.pi/180) + \
            4*earthEccentricity*varY*np.sin(sunAnomolyGeomMean*np.pi/180)*np.cos(2*sunLongGeomMean*np.pi/180) - \
            0.5*varY*varY*np.sin(4*sunLongGeomMean*np.pi/180) - \
            1.25*earthEccentricity*earthEccentricity*np.sin(2*sunAnomolyGeomMean*np.pi/180))*180/np.pi; #min,
        #--- Sun Sub Solar Point (longitude)  ---
        sunSubSolar = -15*(dayOfHrs - 12 + eqOfTime/60); #deg, sub solar point (sun longitude pt where azimuth = 0 deg) [from wikipedia]
        sunSubSolar_loc['long'][j] = sunSubSolar; #record
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
    if( FLG_inputNotListed ):
        sunSubSolar_loc['lat'] = sunSubSolar_loc['lat'][0]; #unwrap b/c came in not as a list
        sunSubSolar_loc['long'] = sunSubSolar_loc['long'][0]; #unwrap b/c came in not as a list
    #END IF
    
    return sunSubSolar_loc
#END DEF
