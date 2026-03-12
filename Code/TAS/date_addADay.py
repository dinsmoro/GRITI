#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def date_addADay(date, days2add):
    if( isinstance(date, (list, tuple)) ):
       origType = type(date); # Record
       date = np.asarray(date); # Convert to numpy
    else:
        origType = None; # Set
    # END IF
    if( date.size == 3 ):
        ref_monthDays = np.array( (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) ); # Preps number of days in a month
        ref_monthDays_leap = np.array( (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) ); # Preps number of days in a month
        
        if( date.ndim == 2 ):
            date = date.ravel(); # Flatten
        # END IF
        
        #-----Leap year detection-----
        if( ((np.mod(date[0],4) == 0) and (np.mod(date[0],100) != 0)) or (np.mod(date[0],400) == 0) ): # Leap year
            monthDays = ref_monthDays_leap; # Use the leap days in a month
        else:
            monthDays = ref_monthDays; # Use the regular days in a month
        # END IF
        
        date[-1] += days2add; # Bamzo
        if( days2add >= 0 ): # FORWARDS
            while( date[-1] > monthDays[date[1]-1] ):
                date[-1] -= monthDays[date[1]-1]; # Subtract off the days in the month
                date[1] += 1; # Move up to the next month
                if( date[1] > 12 ):
                    date[0] += 1; # Move up the year
                    date[1] = 1; # Set to January
                    # Re-evaluate if leap year
                    if( ((np.mod(date[0],4) == 0) and (np.mod(date[0],100) != 0)) or (np.mod(date[0],400) == 0) ): # Leap year
                        monthDays = ref_monthDays_leap; # Use the leap days in a month
                    else:
                        monthDays = ref_monthDays; # Use the regular days in a month
                    # END IF
                # END IF
            # END WHILE
        else: # BACKWARDS
            while( date[-1] < 1 ):
                if( date[1]-2 >= 0 ):
                    date[-1] += monthDays[date[1]-2]; # Add on the days in the previous month
                    date[1] -= 1; # Move up to the previous month
                else:
                    date[-1] += monthDays[11]; # Add on the days in December (only Feb is afflicted with leaping)
                    date[1] = 12; # Set to December
                    date[0] -= 1; # Move back a year
                    # Re-evaluate if leap year
                    if( ((np.mod(date[0],4) == 0) and (np.mod(date[0],100) != 0)) or (np.mod(date[0],400) == 0) ): # Leap year
                        monthDays = ref_monthDays_leap; # Use the leap days in a month
                    else:
                        monthDays = ref_monthDays; # Use the regular days in a month
                    # END IF
                # END IF
            # END WHILE
        # END IF
    elif( date.size == 2 ):
        if( date.ndim == 2 ):
            date = date.ravel(); # Flatten
        # END IF
        
        date[-1] += days2add; # Bamzo
        if( days2add >= 0 ): # FORWARDS
            #-----Leap year detection-----
            if( ((np.mod(date[0],4) == 0) and (np.mod(date[0],100) != 0)) or (np.mod(date[0],400) == 0) ): # Leap year
                yearDays = 366; # Use the leap days in a year
            else:
                yearDays = 365; # Use the regular days in a year
            # END IF
            while( date[-1] > yearDays ):
                date[-1] -= yearDays; # Subtract off the days in the year
                date[0] += 1; # Move up to the next year
    
                # Re-evaluate if leap year
                if( ((np.mod(date[0],4) == 0) and (np.mod(date[0],100) != 0)) or (np.mod(date[0],400) == 0) ): # Leap year
                    yearDays = 366; # Use the leap days in a year
                else:
                    yearDays = 365; # Use the regular days in a year
                # END IF
            # END WHILE
        else: # BACKWARDS
            while( date[-1] < 0 ):
                date[0] -= 1; # Move up to the previous year
    
                # Re-evaluate if leap year
                if( ((np.mod(date[0],4) == 0) and (np.mod(date[0],100) != 0)) or (np.mod(date[0],400) == 0) ): # Leap year
                    yearDays = 366; # Use the leap days in a year
                else:
                    yearDays = 365; # Use the regular days in a year
                # END IF
                
                date[-1] += yearDays; # Add on the days in the previous year
            # END WHILE
        # END IF
    else:
        raise Exception('ERROR in date_addADay: date provided is NOT size 3 (yr/mon/day) or size 2 (yr/dayNum). Ditching.');
    # END IF
    
    if( origType is None ):
        return date
    else:
        return origType(date) # Re-convert
    # END IF
# END DEF