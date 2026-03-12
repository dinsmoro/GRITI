#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RIGHT now THIS cannot WORK because of LOCAL functions BEING NOT HERE
It prob still won't work on ones you import, but there's a glimmer of hope it might
Currently copy into your function and bring any features bakc here
"""

from multiprocessing import Pool, cpu_count
from os import nice as osnice
from platform import system as platformsystem
import numpy as np

def parallel_starmap_helper(fn, args, kwargs=None):
    #for multiprocess starmap with or without kwargs
    # MUST be outside of the function that calls it or it just hangs
    if( kwargs is None ):
        return fn(*args)
    else:
        return fn(*args, **kwargs)
    # END IF
# END DEF
    
def beNice(niceness=0): # Make nice level niceness if it's not already
    if( niceness != 0 ):
        nice = osnice(0) # Add zero to nice and return it 
        if nice < niceness: 
            osnice(niceness-nice) # Increment nice to niceness
        # END IF
    # END IF
# END DEF

def parallel_orchastrator(parallel_list, parallel_CPUnum=None, parallel_CPUmargin=0, parallel_niceNess=0): # parallel_CPUmargin is in whole percentages like 10 == 10%, NOT 0.1 == 10%
    # Lowers the amount of parallel helper code needed in the main function
    
    # --- Prep Parallel ---
    # Get the CPU num to use
    coresInstalled = cpu_count(); # Get the number of CPU cores here (probably includes hyperthreads on Intel systems and SMT threads on AMD systems)
    if( parallel_CPUnum is None ):
        if( parallel_CPUmargin == 0 ):
            parallel_CPUnum = coresInstalled; # Use multiprocessing to get the cpu count, it includes tiny cores and does not figure out who is who yolo
        else:
            parallel_CPUnum = int(np.ceil(coresInstalled*(100-parallel_CPUmargin)/100));
            if( (parallel_CPUnum == coresInstalled) and (parallel_CPUnum > 1) ):
                parallel_CPUnum -= 1; # Make sure one core is free if the % didn't manage to get a core free - margin implies you don't want ALL
            # END IF
        # END IF
    else:
        if( coresInstalled < parallel_CPUnum ):
            # from warnings import warn
            # warn('WARNING: Number of parallel processes requested `'+str(parallel_CPUnum)+'` exceeds the number of CPU cores found `'+str(coresInstalled)+'`. Limiting to '+str(coresInstalled)+' parallel processes to not oversubscribe the CPU.');
            parallel_CPUnum = coresInstalled; # Directly limit so cores not oversubscribed
        # END IF
        if( parallel_CPUmargin != 0 and (parallel_CPUnum > int(np.ceil(coresInstalled*(100-parallel_CPUmargin)/100))) ):
            parallel_CPUnum = int(np.ceil(coresInstalled*(100-parallel_CPUmargin)/100)); # Limit it to enforce the margin requested
            if( (parallel_CPUnum == coresInstalled) and (parallel_CPUnum > 1) ):
                parallel_CPUnum -= 1; # Make sure one core is free if the % didn't manage to get a core free - margin implies you don't want ALL
            # END IF
        # END IF
    # END IF
    if( parallel_CPUnum > len(parallel_list) ):
        parallel_CPUnum = len(parallel_list); # Limit parallel CPU number to the number of things to work on
    # END IF
        
    # --- Get if Linux to Activate Nice ---
    if( platformsystem() == 'Linux' ): # Only Linux has a working implementation of the scheduling system "nice"
        #--- Execute function on list of inputs ---
        with Pool(processes=parallel_CPUnum, initializer=beNice, initargs=(parallel_niceNess,)) as executor:
             results = executor.starmap(parallel_starmap_helper, parallel_list); #function you want is replaced with; parallel_starmap_helper helps starmap distribute everything right
        #END WITH
    else:
        #--- Execute function on list of inputs ---
        with Pool(processes=parallel_CPUnum) as executor:
             results = executor.starmap(parallel_starmap_helper, parallel_list); #function you want is replaced with; parallel_starmap_helper helps starmap distribute everything right
        #END WITH   
    # END WITH
    
    return results
# END DEF