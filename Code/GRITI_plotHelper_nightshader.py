# Roughly based on Cartopy nightshade, shout out to them for shapley interaction

import numpy as np
from datetime import timezone
import shapely.geometry as sgeom
from cartopy.feature.__init__ import ShapelyFeature
import cartopy.crs as ccrs
# from Code.subfun_sunAlsoRises_location_dt import sunAlsoRises_location_dt

#!! LAT ADJUSTMENT IS NOT FINISHED, I JUST GOT OTHER STUFF TO DO SO OH WELL !!

class nightshader(ShapelyFeature):
    def __init__(self, date=None, stepSize=0.1, atmoRefrac=-0.833, latLongShift=None,
                 color='k', alpha=0.5, **kwargs):
        """
        Shade the darkside of the Earth, accounting for atmoRefrac.

        Parameters
        ----------
        date : datetime obj
            Datetime object of when you want the Sun's position, supports non-UTC via shifting to UTC as well
            Default: crash
            
        stepSize : float
            Step size in degrees to determine the resolution of the
            night polygon feature (``npts = 180 / stepSize``).
            
        atmoRefrac : float
            Atmospheric refraction angle in degrees, see https://gml.noaa.gov/grad/solcalc/calcdetails.html

        latLongShift : {'lat':latShift, 'long':longShift, 'rel':False}
            Shift lat/long by something, relative=False replaces long and adjusts latitude by the difference
            relative=True adds 'lat'/'long' dict values to Sun's true position [-180 to 180 enforcement is always on so anything is safe]

        Note
        ----
            Matplotlib keyword arguments can be used when drawing the feature.
            This allows standard Matplotlib control over aspects such as
            'color', 'alpha', etc.

        """
        if date == None:
            print('Error: Can\'t be no date, I don\'t fill it in for you. Yeeting ya!');
            np.crash(); #it doesn't exist
        #END IF

        # Get lat/long
        # latLong = sunAlsoRises_location_dt(date);
        # lat = latLong['lat'];
        # long = latLong['long'];
        lat, long = sunAlsoRises_location_dt(date);
        
        # if( latLongShift != None ):
        #     if( latLongShift['rel'] == True ):
        #         lat += latLongShift['lat']; #relatively adjust
        #         long += latLongShift['long']; #relatively adjust
        #     else:
        #         lat = lat + (latLongShift['lat'] - lat); #not relatively adjust
        #         long = long + (latLongShift['long'] - long); #not relatively adjust
        #     #END IF
        # #END IF
        
        
        
        if lat > 0:
            pole_lat = -90 + lat
            central_lon = 180
        else:
            pole_lat = 90 + lat
            central_lon = 0
        #END IF
        if( latLongShift != None ):
            if( latLongShift['rel'] == True ):
                pole_lat += latLongShift['lat']; #relatively adjust
                long += latLongShift['long']; #relatively adjust
            else:
                pole_lat = pole_lat + (latLongShift['lat'] - lat); #not relatively adjust
                pole_long = latLongShift['long']; #abrupt shift over
            #END IF
        else:
            pole_long = long
        #END IF
        pole_long = ((pole_long + 180) % 360) - 180; #wrap around enforcer, keeps it -180 to 180


        rotated_pole = ccrs.RotatedPole(pole_latitude=pole_lat,
                                        pole_longitude=pole_long,
                                        central_rotated_longitude=central_lon)

        npts = int(180 / stepSize)
        # longs = np.empty(npts * 2)
        lats = np.empty(npts * 2)
        
        # Solve the equation for sunrise/sunset:
        # https://en.wikipedia.org/wiki/Sunrise_equation#Generalized_equation
        # NOTE: In the generalized equation on Wikipedia,
        #       stepSize == 0. in the rotated pole coordinate system.
        #       Therefore, the max/min latitude is +/- (90+atmoRefrac)

        #EQ:
        #omega0 = acos( (sin(atmoRefract) - sin(obsLat)*sin(sunDec[sunLat])) / (cos(obsLat)*cos(sunDec[sunLat])) )
        # obsLat = 0 I guess for equator or something
        #omega0 = acos( sin(atmoRefract) / cos(sunDec[sunLat]) )

        # Fill latitudes up and then down
        lats[:npts] = np.linspace(-(90 + atmoRefrac), 90 + atmoRefrac, npts)
        lats[npts:] = lats[:npts][::-1]

        # Solve the generalized equation for omega0, which is the
        # angle of sunrise/sunset from solar noon
        # We need to clip the input to arccos to [-1, 1] due to floating
        # point precision and arccos creating nans for values outside
        # of the domain
        arccos_tmp = np.clip(np.sin(np.deg2rad(atmoRefrac)) /
                             np.cos(np.deg2rad(lats)), -1, 1)
        longs = np.rad2deg(np.arccos(arccos_tmp)) #this is omegaw on wiki page

        # Fill the longitude values from the offset for midnight.
        # This needs to be a closed loop to fill the polygon.
        # Negative longitudes
        longs[:npts] = -(180 - longs[:npts])
        # Positive longitudes
        longs[npts:] = 180 - longs[npts:]

        kwargs.setdefault('facecolor', color)
        kwargs.setdefault('alpha', alpha)

        geom = sgeom.Polygon(np.column_stack((longs, lats)))
        return super().__init__([geom], rotated_pole, **kwargs) #what even is this, yeah readable python ok
    #END DEF
#END CLASS


# #based on NOAA example https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
def sunAlsoRises_location_dt(dater):
        
    #--- Ensure datetime obj is in UTC ---
    if( dater.tzinfo != None ):
        # if( (date.tzinfo != pytz.utc) | (date.tzinfo != timezone.utc) ): #remove check for now
        dater = dater.astimezone(timezone.utc); #convert as needed
        #END IF
    #END IF
    
    #--- Julian Day from YR/M/D ---
    dayOfHrs = (dater.hour + dater.minute/60 + dater.second/3600); #hours in the day
    JD = 367*dater.year - np.int64(7*(dater.year+np.int64((dater.month+9)/12))/4) - \
        np.int64(3*(np.int64((dater.year+(dater.month-9)/7)/100)+1)/4) + \
        np.int64(275*dater.month/9) + dater.day + 1721028.5; #JD, from https://scienceworld.wolfram.com/astronomy/JulianDate.html
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
    #--- Sun True Longitude ---
    sunLongTrue = sunLongGeomMean + sunEqOfCtr; #deg, 
    #--- Sun Appogee Longitude ---
    sunLongAppogee = sunLongTrue - 0.00569 - 0.00478*np.sin((125.04 - 1934.136*JC)*np.pi/180); #deg,
    #--- Mean Oblique Elliptic & Oblique Correction & "var y" ---
    obliqueMeanElliptic = 23 + (26 + ((21.448 - JC*(46.815 + JC*(0.00059 - JC*0.001813))))/60)/60; #deg,
    obliqueCorr = obliqueMeanElliptic + 0.00256*np.cos((125.04 - 1934.136*JC)*np.pi/180); #deg, 
    varY = np.tan((obliqueCorr/2)*np.pi/180)*np.tan((obliqueCorr/2)*np.pi/180);
    #--- Sun Declination (latitude) ---
    lat = np.arcsin(np.sin(obliqueCorr*np.pi/180)*np.sin(sunLongAppogee*np.pi/180))*180/np.pi; #deg, sun declination
    #--- Eq of Time ---
    eqOfTime = 4*(varY*np.sin(2*sunLongGeomMean*np.pi/180) - 2*earthEccentricity*np.sin(sunAnomolyGeomMean*np.pi/180) + \
        4*earthEccentricity*varY*np.sin(sunAnomolyGeomMean*np.pi/180)*np.cos(2*sunLongGeomMean*np.pi/180) - \
        0.5*varY*varY*np.sin(4*sunLongGeomMean*np.pi/180) - \
        1.25*earthEccentricity*earthEccentricity*np.sin(2*sunAnomolyGeomMean*np.pi/180))*180/np.pi; #min,
    #--- Sun Sub Solar Point (longitude)  ---
    long = -15*(dayOfHrs - 12 + eqOfTime/60); #deg, sub solar point (sun longitude pt where azimuth = 0 deg) [from wikipedia]
    
    return lat, long
#END DEF
