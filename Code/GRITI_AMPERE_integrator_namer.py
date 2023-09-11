#AMPERE integrator
# 0 = Integrate within plot lat/long range
# 1 = Integrate within plot long range with max latitude of 90
# 2 = Integrate within plot long range with user-defined lat max value
# 3 = Integrate entire hemisphere
# 4 = Integrate within plot long range with user-defined lat min value and max of 90
# 5 = Integrate entire hemisphere within user-defined lat min value and max of 90 (so all long)
# 6 = Integrate entire hemisphere less than user-defined lat value to equator (and all long values)
# 7 = Integrate at a radius around a point (it's back!)
import numpy as np

def GRITI_AMPERE_integrator_namer(AMPERE_integrateMethod, AMPERE_integrateMethod_val, plotLatRange=None, FLG_partial=False):
    if( AMPERE_integrateMethod == 0 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' within keo area'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated within keo area'; #add to mecha title
        #END IF
    elif( AMPERE_integrateMethod == 1 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' within keo long & up to pole'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated within keo long & up to pole'; #add to mecha title
        #END IF
    elif( AMPERE_integrateMethod == 2 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' within keo long & up to '+str(AMPERE_integrateMethod_val)+' degc lat'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated within keo long & up to '+str(AMPERE_integrateMethod_val)+' degc lat'; #add to mecha title
        #END IF
    elif( AMPERE_integrateMethod == 3 ):
        if( np.any(plotLatRange == None) ):
            print('WARNING in GRITI_AMPERE_integrator_namer: plotLatRange was NOT supplied but needed. Assuming Northern Hemisphere, 33% chance!');
            if( FLG_partial == True ):
                AMPERE_integrateMethod_string = ' Northern Hemisphere'; #add to mecha title
            else:
                AMPERE_integrateMethod_string = 'AMPERE integrated within the Northern Hemisphere'; #add to mecha title
            #END IF
        else:
            if( (np.min(plotLatRange) <= 0) & (np.max(plotLatRange) >= 0) ):
                if( FLG_partial == True ):
                    AMPERE_integrateMethod_string = ' both Hemispheres'; #add to mecha title
                else:
                    AMPERE_integrateMethod_string = 'AMPERE integrated within both Hemispheres'; #add to mecha title
                #END IF
            else:
                if( (np.min(plotLatRange) >= 0) & (np.max(plotLatRange) >= 0) ):
                    #northern hemisphere
                    if( FLG_partial == True ):
                        AMPERE_integrateMethod_string = ' Northern Hemisphere'; #add to mecha title
                    else:
                        AMPERE_integrateMethod_string = 'AMPERE integrated within the Northern Hemisphere'; #add to mecha title
                    #END IF
                else:
                    #southern hemisphere
                    if( FLG_partial == True ):
                        AMPERE_integrateMethod_string = ' Southern Hemisphere'; #add to mecha title
                    else:
                        AMPERE_integrateMethod_string = 'AMPERE integrated within the Southern Hemisphere'; #add to mecha title
                    #END IF
                #END IF
            #END IF
        #END IF
    elif( AMPERE_integrateMethod == 4 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' within keo long & from '+str(AMPERE_integrateMethod_val)+' degc lat & up to pole'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated within keo long & from '+str(AMPERE_integrateMethod_val)+' degc lat & up to pole'; #add to mecha title
        #END IF
    elif( AMPERE_integrateMethod == 5 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' within entire hemisphere & from '+str(AMPERE_integrateMethod_val)+' degc lat & up to pole'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated within entire hemisphere & from '+str(AMPERE_integrateMethod_val)+' degc lat & up to pole'; #add to mecha title
        #END IF
    elif( AMPERE_integrateMethod == 6 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' within entire hemisphere & up to '+str(AMPERE_integrateMethod_val)+' degc lat'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated within entire hemisphere & up to '+str(AMPERE_integrateMethod_val)+' degc lat'; #add to mecha title
        #END IF
    elif( AMPERE_integrateMethod == 7 ):
        if( FLG_partial == True ):
            AMPERE_integrateMethod_string = ' around '+AMPERE_integrateMethod_val[0][0]+' lat, '+AMPERE_integrateMethod_val[0][1]+' long with a radius of '+AMPERE_integrateMethod_val[1]+' km'; #add to mecha title
        else:
            AMPERE_integrateMethod_string = 'AMPERE integrated around '+AMPERE_integrateMethod_val[0][0]+' lat, '+AMPERE_integrateMethod_val[0][1]+' long with a radius of '+AMPERE_integrateMethod_val[1]+' km'; #add to mecha title
        #END IF
    #END IF
    
    return AMPERE_integrateMethod_string
#END DEF