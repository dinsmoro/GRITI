#GOAL: Get sunrise and sunset times
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: sunrise suinset
#based on NOAA example https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#I implemented all bits cause I can, not all are used but they seemed like they might be useful at some point idk
#Note that timeIndexes can really be any form of time in seconds associated with the days (prolly)

import numpy as np #import in here I dunno
from numba import jit

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def sunAlsoRises_zenith_numba(lat, long, year, month, day, hour, minute, sec, tzOffset = 0):
         
    #--- Compile local time ---
    localTime = hour + minute/60 + sec/3600; #hr, local time - tzOffset will convert it to UTC (tzOffset=0 if already UTC!)
    
    #--- Julian Day from YR/M/D ---
    JD = 367*year - np.int64(7*(year+np.int64((month+9)/12))/4) - \
        np.int64(3*(np.int64((year+(month-9)/7)/100)+1)/4) + \
        np.int64(275*month/9) + day + 1721028.5; #JD, from https://scienceworld.wolfram.com/astronomy/JulianDate.html
    JD += localTime/24 - tzOffset/24; #JD, add on the desired local time (12 noon) and adjust that local time by the timezone offset from UTC to give it in a UTC time
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
    # sunSubSolar = -15*(localTime - tzInt - 12 + eqOfTime/60); #deg, sub solar point (sun longitude pt where azimuth = 0 deg) [from wikipedia]
    # #--- Hour Angle of Sunrise ---
    # sunHA = (np.arccos(np.cos(90.833*np.pi/180)/(np.cos(lat*np.pi/180)*np.cos(sunDeclination*np.pi/180)) - \
    #     np.tan(lat*np.pi/180)*np.tan(sunDeclination*np.pi/180)))*180/np.pi; #deg,
    # #--- Local Solar Noon ---
    # sunSolarNoon = (720 - 4*long - eqOfTime+np.float64(tzInt)*60)/1440; #days, fraction of a day that is local solar noon
    # #--- Local Sunrise Time ---
    # sunRise_local = (sunSolarNoon*1440 - sunHA*4)/1440; #days, fraction of a day that is local sunrise
    # sunRise = sunRise_local - tzInt/24; #days, convert to UTC
    # #--- Local Sunset Time ---
    # sunSet_local = (sunSolarNoon*1440 + sunHA*4)/1440; #days, fraction of a day that is local sunrise
    # sunSet = sunSet_local - tzInt/24; #days, convert to UTC
    # #--- Local Sun Duration ---
    # sunDuration = 8*sunHA; #minutes
    #--- True Solar Time ---
    sunTrueSolarTime = np.mod(localTime/24*1440+eqOfTime+4*long-60*tzOffset,1440); #min,
    #--- Hour Angle ---
    if( sunTrueSolarTime/4 < 0 ):
        HA = sunTrueSolarTime/4 + 180; #deg,
    else:
        HA = sunTrueSolarTime/4 - 180; #deg,
    #END IF
    #--- Solar Zenith Angle ---
    sunZenith = np.arccos(np.sin(lat*np.pi/180)*np.sin(sunDeclination*np.pi/180)+np.cos(lat*np.pi/180)*np.cos(sunDeclination*np.pi/180)*np.cos(HA*np.pi/180))*180/np.pi; #deg, 
    #--- Solar Elevation Angle ---
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
    
    return sunZenith
#END DEF
