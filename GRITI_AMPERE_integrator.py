#AMPERE integrator
# 0 = Integrate within plot lat/long range
# 1 = Integrate within plot long range with max latitude of 90
# 2 = Integrate within plot long range with user-defined max value
# 3 = Integrate entire hemisphere
import numpy as np

def GRITI_AMPERE_integrator(AMPERE_data, dates, settings_AMPERE, plotLatRange, plotLongRange, AMPERE_integrateMethod, AMPERE_integrateMethod_val, AMPERE_integrateMethod_log=0):
    #-----Unpack-----
    AMPERE_timeUnique = AMPERE_data['time unique'];
    # AMPERE_integrateMethod = settings_AMPERE['integrate method'];
    # AMPERE_integrateMethod_val = settings_AMPERE['integrate method lat val'];
    # AMPERE_integrateMethod_log = settings_AMPERE['integrate method log'];
    # plotLatRange = settings['map']['lat range'];
    # plotLongRange = settings['map']['long range'];
    AMPERE_dataType = settings_AMPERE['data type']; #get the current AMPERE data type
    
    #prep to plot
    if( AMPERE_integrateMethod == 0 ):
        #regular
        kInRange = (AMPERE_data['lat'] >= np.min(plotLatRange)) & (AMPERE_data['lat'] <= np.max(plotLatRange)) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 1 ):
        #max is 90 now
        kInRange = (AMPERE_data['lat'] >= np.min(plotLatRange)) & (AMPERE_data['lat'] <= 90) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 2 ):
        #max is user-defined now
        kInRange = (AMPERE_data['lat'] >= np.min(plotLatRange)) & (AMPERE_data['lat'] <= AMPERE_integrateMethod_val) & \
            (AMPERE_data['long'] >= np.min(plotLongRange)) & (AMPERE_data['long'] <= np.max(plotLongRange)); #get data in the range
    elif( AMPERE_integrateMethod == 3 ):
        if( (np.min(plotLatRange) <= 0) & (np.max(plotLatRange) >= 0) ):
            kInRange = np.ones(AMPERE_data['lat'].shape,dtype=np.bool_); #all are good to go, both hemisphere integration yolo
        else:
            if( (np.min(plotLatRange) >= 0) & (np.max(plotLatRange) >= 0) ):
                #northern hemisphere
                kInRange = AMPERE_data['lat'] >= 0; #get data in the range
            else:
                #southern hemisphere
                kInRange = AMPERE_data['lat'] <= 0; #get data in the range
            #END IF
        #END IF
    #END IF  
    
    AMPERE_integrate = np.zeros( AMPERE_timeUnique.size , dtype=np.float64); #prep integrated joule heating
    for i in range(AMPERE_timeUnique.size):
        k = AMPERE_timeUnique[i] == AMPERE_data['time']; #get the right time
        AMPERE_integrate[i] = np.sum(AMPERE_data[AMPERE_dataType][k&kInRange]); #ergs/(cm^2*sec), get the Joule Heating for the current time stamp
    #END FOR i
    if( AMPERE_integrateMethod_log == 1 ):
        #Can't log 0's, interpolate over them
        k = np.where( AMPERE_integrate == 0 )[0]; #get where 0's are
        AMPERE_integrate[k] = 1; #set to 1 to prevent error
        AMPERE_integrate = np.log10(AMPERE_integrate); #log it
    #END IF
    
    return AMPERE_integrate