#GOAL: Get sunrise and sunset times
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: sunrise suinset
#based on NOAA example https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#I implemented all bits cause I can, not all are used but they seemed like they might be useful at some point idk


import numpy as np #import in here I dunno
import timezonefinder
from datetime import datetime
import pytz
from Code.subfun_strstr import strstr
from urllib.request import urlopen

def sunAlsoRises_timeZone(latLong_ref, plotLatRange, plotLongRange, dateRange_zeroHr):
    #FIRST BATTLE: REGION WHERE LAT/LONG IS   
    if( (np.min(plotLatRange) <= latLong_ref[0][0]) & (np.max(plotLatRange) >= latLong_ref[0][0]) & \
       (np.min(plotLongRange) <= latLong_ref[0][1]) & (np.max(plotLongRange) >= latLong_ref[0][1]) ):
        # Only use latLong_ref if its within the plot area, otherwise ditch it b/c it's set wronk
        latLong_use = (latLong_ref[0][0],latLong_ref[0][1]);
    else:
        print('WARNING: latLong_ref[0] '+str(np.round(latLong_ref[0][0],2)).rstrip('0').rstrip('.')+' lat | '+str(np.round(latLong_ref[0][1],2)).rstrip('0').rstrip('.')+\
              ' long is not within lat range of '+str(np.min(plotLatRange))+' to '+str(np.max(plotLatRange))+\
              ' | long range of '+str(np.min(plotLongRange))+' to '+str(np.max(plotLongRange))+'. Ignoring latLong_ref[0] and using mean of plotLat/LongRange.'+\
              '('+str(np.round(np.mean(plotLatRange),2)).rstrip('0').rstrip('.')+' lat | '+str(np.round(np.mean(plotLongRange),2)).rstrip('0').rstrip('.')+' long)'); # Report a warning on a wrong setting
        latLong_use = (np.mean(plotLatRange),np.mean(plotLongRange)); # Set to mean of plotLat/LongRange
    #END IF    
    tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
    dayNite_timeZoneID = tf.certain_timezone_at(lat=latLong_use[0], lng=latLong_use[1]); #use it to find the time zone
    if dayNite_timeZoneID is None:
        #use geonames site as a backup
        url = 'http://api.geonames.org/timezone?lat='+str(latLong_use[0])+'&lng='+str(latLong_use[1])+'&username=razzluhdzuul'; #create link for lat/long
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

    timeZoneObj_zeroHr = timeZoneObj.localize(datetime.strptime(str(dateRange_zeroHr), '[%Y\t%m\t%d]')); #time zone info at zero hr

    dayNite_DSToffset_str = timeZoneObj_zeroHr.strftime('%z').replace('0',''); #remove the 0 that was extraneous

    # dayNite_DSToffset = np.int64(dayNite_DSToffset_str); #get the number version

    dayNite_timeZoneName = timeZoneObj_zeroHr.tzname(); #get the time zone name (like 'EST' or 'EDT' depending on standard or daylight savings time)
    # dateNite_DSTnUTCOffset = timeZoneObj_zeroHr.dst().total_seconds()/3600; #get the time offset
    # if( np.mod(dateNite_DSTnUTCOffset,1) == 0 ):
    #     dateNite_DSTnUTCOffset = np.int64(dateNite_DSTnUTCOffset); #convert to integer
    # #END IF
    
    return dayNite_timeZoneName, dayNite_DSToffset_str