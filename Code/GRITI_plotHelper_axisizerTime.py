"""
GOAL: Initiate a Cartopy plot
RD on 6/11/21

INPUT: 
    timez as you'd like to plot
    tickNumGoal [default 19] - number of tick marks on bottom axis goal, it won't exceed this but will try to get near it
    tickReducer [default 3] - for each extra character on the ticks, reduce the number of ticks by this amount
    tickMin [default 4] - minimum number of ticks to have
    FLG_force24max [default False] - ignores automagic sizing and forces 24 hour ticks even if they may not fit [only for unit='hr']
    FLG_force60max [default False] - ignores automagic sizing and forces 60 minute/sec ticks even if they may not fit [only for unit='min'/'sec']
    FLG_manualLims [default None] - allows for manual time limits to be set (for time cutout stuff prolly), must be a list/tuple/numpy array of 2 values with [0] being the lesser
    FLG_removeLabels [default False] - True removes labels on the x axis (for subplots)
    FLG_tickDirIn [default False] - True turns tick direction inward instead of default outward
OUTPUT:
    timez axis ticks
    timez axis lims
"""

import numpy as np #import in here I dunno
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.ticker as tick
# import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated

def GRITI_plotHelper_axisizerTime(timez,ax=None,unit='hr',tickNumGoal=28,tickReducer=3,tickMin=4,FLG_force24max=False,FLG_force60max=False,FLG_manualLims=None,FLG_removeLabels=False,FLG_tickDirIn=False):
        
    #----- Axis timez Ticks -----
    #character number estimtor
    charEst = np.ones(timez.shape,dtype=np.int64); #preallocate
    kj = (np.round(timez,0) != 0) & (~np.isinf(timez)) & (~np.isnan(timez)); #anywhere that is 0 has 1 digit, so ones already takes care of it
    charEst[kj] = np.int64(np.log10(np.abs(np.round(timez[kj],0))))+1; #get the number of characters
    charEst[timez < 0] = charEst[timez < 0]+1; #xtra character for negative sign
    charEst = charEst.max(); #get the most characters we're dealing with
    
    tickNumGoal += -charEst*tickReducer; #reduce the ticks by the number of characters we're dealing with
    
    #support manual time limits
    if( np.any(FLG_manualLims == None) ):
        FLG_manualLims = np.array( (np.min(timez),np.max(timez)) ); #populate the limits to use as the min/max of the time if none provided
    #END IF
    
    timez_autoTick = (np.ceil(FLG_manualLims[1]) - np.floor(FLG_manualLims[0]))/tickNumGoal; #tries to split the time range into tickNumGoal # of times
    if( (unit.lower() == 'hr') | (unit.lower() == 'hour') | (unit.lower() == 'h') ):
        if( timez_autoTick > 12 ):
            if( FLG_force24max == False ):
                if( timez_autoTick < 24 ):
                    timez_autoTick = 24;
                else:
                    timez_autoTick = 24*(timez_autoTick//24); #sets the tick setting to 15 arcdegrees per tick
                #END IF
            else:
                timez_autoTick = 24; #forced to 24 otherwise
            #END IF
        elif( timez_autoTick > 6 ):
            timez_autoTick = 12; #sets the tick setting to 15 arcdegrees per tick
        elif( timez_autoTick > 4 ):
            timez_autoTick = 6; #sets the tick setting to 10 arcdegrees per tick
        elif( timez_autoTick > 2 ):
            timez_autoTick = 4; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick > 1 ):
            timez_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick >= 0.5 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            timez_autoTick = 1; # sets the tick setting to 1 arcdegree per tick   
        elif( timez_autoTick >= 0.1 ):
            timez_autoTick = 0.5;
        else:
            timez_autoTick = (FLG_manualLims[1] - FLG_manualLims[0])/tickNumGoal; #just goes for it if it's a super tiny range
        #END IF
    else:
        if( timez_autoTick > 60 ):
            if( FLG_force60max == False ):
                if( timez_autoTick < 60 ):
                    timez_autoTick = 60;
                else:
                    timez_autoTick = 60*(timez_autoTick//60); #sets the tick setting to 15 arcdegrees per tick
                #END IF
            else:
                timez_autoTick = 60; #forced to 60 otherwise
            #END IF
        elif( timez_autoTick > 30 ):
            timez_autoTick = 60; #sets the tick setting to 15 arcdegrees per tick
        elif( timez_autoTick > 20 ):
            timez_autoTick = 30; #sets the tick setting to 10 arcdegrees per tick
        elif( timez_autoTick > 15 ):
            timez_autoTick = 20; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick > 10 ):
            timez_autoTick = 15; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick > 5 ):
            timez_autoTick = 10; #sets the tick setting to 5 arcdegrees per tick   
        elif( timez_autoTick > 3 ):
            timez_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick > 2 ):
            timez_autoTick = 3; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick > 1 ):
            timez_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
        elif( timez_autoTick >= 0.5 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            timez_autoTick = 1; # sets the tick setting to 1 arcdegree per tick
        elif( timez_autoTick >= 0.1 ):
            timez_autoTick = 0.5;
        else:
            timez_autoTick = (FLG_manualLims[1] - FLG_manualLims[0])/tickNumGoal; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    #uses mod to force sensible numbers
    timezAxisTicks = np.arange( (np.round(FLG_manualLims[0]) - np.mod(np.round(FLG_manualLims[0]),timez_autoTick)) , \
            (np.round(FLG_manualLims[1]) - np.mod(np.round(FLG_manualLims[1]),timez_autoTick)) + timez_autoTick , \
            timez_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
        
    #----- Axis timez Limits -----
    #makes plotting look better
    timezAxisLims = np.array((FLG_manualLims[0] , FLG_manualLims[1])); #get those xAxis limits, ensures everything is aligned
    if( (np.round(FLG_manualLims[0]) - FLG_manualLims[0]) < 5*10**-2 ): #fix x axis stuff (timez)
        timezAxisLims[0] = np.round(FLG_manualLims[0]); #force the rounded value
    #END IF
    if( (np.round(FLG_manualLims[1]) - FLG_manualLims[1]) < 5*10**-2 ): #fix x axis stuff (timez)
        timezAxisLims[1] = np.round(FLG_manualLims[1]); #force the rounded value
    #END IF
    
    #----- Apply to axes if provided -----
    if( ax is not None ):
        ax.set_xticks(timezAxisTicks); #set x axis ticks
        ax.set_xlim(timezAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
    #END IF
    
    #----- Remove axis labels if reqd -----
    if( FLG_removeLabels == True ):
        # ax.set_xticklabels([]); #remove the labels
        ax.tick_params(labelbottom=False); #remove the labels but better
    #END IF
    
    #----- Set tick dir inward if reqd -----
    if( FLG_tickDirIn == True ):
        ax.tick_params(axis='x',direction='in'); #set tick dir in
    #END IF

    return timezAxisTicks, timezAxisLims
#END DEF