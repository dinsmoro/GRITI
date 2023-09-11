#Function to calc day line real fast
#RD on 5/9/2019

#dayOrNite is the only option - 0 for day, 1 for nite

import numpy as np
from scipy.signal import savgol_filter

def GRITI_movieMaker_subfun_dayniteCalc(dayOrNite,dayNite_temp_sum,dayNite_Grid_Long,dayNite_Grid_Lat,dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig):
    #remember dayNite_Grid_Long is actually dayNite_Grid_Long[dayNite_temp] in the main space - for debugging
    
    if( dayOrNite == 0 ):
        dayNite_textRotation_default = -90; #deg, sets to -90 for day time
        dayNite_textRotation_angleAdjust = 0; #deg, adjusts angle calculated to face a certain direction (the original calc'd one is good for day)
    else:
        dayNite_textRotation_default = 90; #deg, sets to 90 for nite time
        dayNite_textRotation_angleAdjust = 180; #deg, adjusts angle calculated to face a certain direction (flips angle calc'd so text faces other way)
    #END IF
    
    #CODE FOR THE LINE LAT/LONG
    dayNite_savgolFiltLen = dayNite_savgolFiltLenOrig; #set it as a constant to start off
    if( dayNite_savgolFiltLen >= dayNite_temp_sum ):
        dayNite_savgolFiltLen = dayNite_temp_sum-1; #set it to this
        if( np.remainder(dayNite_savgolFiltLen,2) == 0 ): #if there's no remainder, it's even. Gotta add 1 cause Sovitsky-Golay filter needs an odd window length
            dayNite_savgolFiltLen = dayNite_savgolFiltLen - 1; #subtract 1 so it's odd
        #END IF
    #END IF
    if( dayNite_savgolFiltLen > 1): #if len is greater than 0 we savgol to SMOOTH
        dayNite_Long_line = savgol_filter(dayNite_Grid_Long,dayNite_savgolFiltLen,1); #filter it up
        dayNite_Lat_line = savgol_filter(dayNite_Grid_Lat,dayNite_savgolFiltLen,1); #filter it up
    else: #if it was 1, time to not savgol it - too little data
        dayNite_Long_line = dayNite_Grid_Long;
        dayNite_Lat_line = dayNite_Grid_Lat;
    #END IF
    
    #CODE FOR THE TEXT LAT/LONG AND ROTATION DEGREE
    dayNite_Long_text = dayNite_Grid_Long[ np.min(np.abs(dayNite_Grid_Lat - dayNite_textLatAbs)) == np.abs(dayNite_Grid_Lat - dayNite_textLatAbs) ][0]; #get the longitude closest to the desired latitude
    dayNite_Lat_text = dayNite_Grid_Lat[ np.min(np.abs(dayNite_Grid_Long - dayNite_Long_text)) == np.abs(dayNite_Grid_Long - dayNite_Long_text) ][0]; #get the latitude closest to the longitude that was found (since it might not be at desired latitude due to plot range)
    dayNite_textRotation = np.where(dayNite_Lat_text == dayNite_Grid_Lat)[0][0]; #use this to hold the index where the latitude is (longitude should be at the same since they're related)
    if( (dayNite_Grid_Long.size > (dayNite_textRotationLenOrig+dayNite_textRotation) ) & (0 <= (dayNite_textRotation-dayNite_textRotationLenOrig) )  ):
        dayNite_textRotationLen = dayNite_textRotationLenOrig; #use original length no problem
        dayNite_textRotation = np.array( (-dayNite_Grid_Long[dayNite_textRotation+dayNite_textRotationLen] + dayNite_Grid_Long[dayNite_textRotation-dayNite_textRotationLen], -dayNite_Grid_Lat[dayNite_textRotation+dayNite_textRotationLen] + dayNite_Grid_Lat[dayNite_textRotation-dayNite_textRotationLen]) ); #create a vector to use with a dot product with the x axis to get an angle off of the x axis
        dayNite_textRotation = dayNite_textRotation_angleAdjust+np.arctan2(dayNite_textRotation[1]*plot_ratio,dayNite_textRotation[0])*180/np.pi; #get the angle in degrees, the plot ratio adjusts so the angle shows up right on the screen
    elif( (dayNite_Grid_Long.size <= (dayNite_textRotationLenOrig+dayNite_textRotation) ) & (0 <= (dayNite_textRotation-dayNite_textRotationLenOrig) )  ):
        dayNite_textRotationLen = dayNite_Grid_Long.size-dayNite_textRotation-1; #calc new length -1 cause python
        if( dayNite_textRotationLen != 0 ):
            dayNite_textRotation = np.array( (-dayNite_Grid_Long[dayNite_textRotation+dayNite_textRotationLen] + dayNite_Grid_Long[dayNite_textRotation-dayNite_textRotationLen], -dayNite_Grid_Lat[dayNite_textRotation+dayNite_textRotationLen] + dayNite_Grid_Lat[dayNite_textRotation-dayNite_textRotationLen]) ); #create a vector to use with a dot product with the x axis to get an angle off of the x axis
            dayNite_textRotation = dayNite_textRotation_angleAdjust+np.arctan2(dayNite_textRotation[1]*plot_ratio,dayNite_textRotation[0])*180/np.pi; #get the angle in degrees, the plot ratio adjusts so the angle shows up right on the screen
        else:
            dayNite_textRotation = dayNite_textRotation_default; #just say it's vertical
        #END IF
    elif( (dayNite_Grid_Long.size > (dayNite_textRotationLenOrig+dayNite_textRotation) ) & (0 > (dayNite_textRotation-dayNite_textRotationLenOrig) )  ):
        dayNite_textRotationLen = dayNite_textRotation; #use new length
        if( dayNite_textRotationLen != 0 ):
            dayNite_textRotation = np.array( (-dayNite_Grid_Long[dayNite_textRotation+dayNite_textRotationLen] + dayNite_Grid_Long[dayNite_textRotation-dayNite_textRotationLen], -dayNite_Grid_Lat[dayNite_textRotation+dayNite_textRotationLen] + dayNite_Grid_Lat[dayNite_textRotation-dayNite_textRotationLen]) ); #create a vector to use with a dot product with the x axis to get an angle off of the x axis
            dayNite_textRotation = dayNite_textRotation_angleAdjust+np.arctan2(dayNite_textRotation[1]*plot_ratio,dayNite_textRotation[0])*180/np.pi; #get the angle in degrees, the plot ratio adjusts so the angle shows up right on the screen
        else:
            dayNite_textRotation = dayNite_textRotation_default; #just say it's vertical
        #END IF
    else:
        dayNite_textRotation = dayNite_textRotation_default; #just say it's vertical
    #END IF
    
    return dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation