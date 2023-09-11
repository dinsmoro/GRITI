"""
GOAL: Initiate a Cartopy plot
RD on 6/11/21

INPUT: 
    latLong as you'd like to plot
    tickNumGoal [default 19] - number of tick marks on bottom axis goal, it won't exceed this but will try to get near it
    tickReducer [default 3] - for each extra character on the ticks, reduce the number of ticks by this amount
    tickMin [default 4] - minimum number of ticks to have
    FLG_force24 [default False] - ignores automagic sizing and forces 24 hour ticks even if they may not fit
OUTPUT:
    latLong axis ticks
    latLong axis lims
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.ticker as tick
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
from Code.subfun_textNice import textNice

def GRITI_plotHelper_axisizerLatLong(latLong,ax=None,axDir=None,tickNumGoal=28,tickReducer=3.5,tickMin=4, \
                                     FLG_forceTick=False,FLG_extendLims=False,FLG_removeLabels=False,FLG_labelUnits=False,FLG_tickDirIn=False):
        
    #----- Prep -----
    if( isinstance(latLong, list) ):
        latLong = np.asarray(latLong); #convert
    #END IF
    
    if( (ax is not None) & (axDir is None) ):
        axDir = 'x'; #assume x dir
    #END IF
    
    if( axDir == 'y' ):
        tickReducer = 0; #no need since not consecutive
        if( tickNumGoal > 17 ):
            tickNumGoal = 17; #reduced
        #END IF
    #END IF
    
    #----- Axis latLong Ticks -----
    if( FLG_forceTick == False ):
        #character number estimtor
        charEst = np.ones(latLong.shape,dtype=np.int64); #preallocate
        kj = np.round(latLong,0) != 0; #anywhere that is 0 has 1 digit, so ones already takes care of it
        charEst[kj] = np.int64(np.log10(np.abs(np.round(latLong[kj],0))))+1; #get the number of characters
        charEst[latLong < 0] = charEst[latLong < 0]+1; #xtra character for negative sign
        charEst = charEst.max(); #get the most characters we're dealing with
        
        tickNumGoal += -np.round(charEst*tickReducer); #reduce the ticks by the number of characters we're dealing with
        
        latLong_delta = np.ceil(np.max(latLong)) - np.floor(np.min(latLong)); # difference to cover
        if( np.isclose(np.mod(latLong_delta,1), 0) ):
            FLG_latLong_delta_isInteger = True;
            latLong_delta = np.int64(latLong_delta); #integer it
        else:
            FLG_latLong_delta_isInteger = False;
        # END IF
        latLong_autoTick = (latLong_delta)/tickNumGoal; #tries to split the time range into tickNumGoal # of times
        
        latLong_autoTick_values =    np.array((360, 180, 90, 60, 45, 30, 15, 10, 5, 3, 2, 1)); #create acceptable array
        latLong_autoTick_threshold = np.array((270, 140, 76, 54, 38, 25, 10, 5, 3, 2, 1, 0.6)); #greater than the threshold # yields the _values tick value
    
        if( FLG_latLong_delta_isInteger == True ):
            latLong_autoTick_values_less = np.isclose(np.mod(latLong_delta/latLong_autoTick_values, 1), 0) & (latLong_delta/latLong_autoTick_values > tickNumGoal*.25); #get evenly dividing values with enough ticks
            if( latLong_autoTick_values_less.sum() > 1 ):
                latLong_autoTick_values_less[-1] = False; #remove 1 as an option if others exist
            # END IF
            if( (latLong_autoTick_values_less.sum() == 1) | ((np.sum(latLong_autoTick >= latLong_autoTick_threshold[latLong_autoTick_values_less]) == 0) & (latLong_autoTick_values_less.sum() > 0)) ):
                latLong_autoTick_where = np.where(latLong_autoTick_threshold[latLong_autoTick_values_less][0] == latLong_autoTick_threshold)[0]; #where ourselves out of the subindexing
            else:
                latLong_autoTick_where = np.where(latLong_autoTick >= latLong_autoTick_threshold[latLong_autoTick_values_less])[0]; #get where they happen
                if( latLong_autoTick_where.size != 0 ):
                    latLong_autoTick_where = np.where(latLong_autoTick_threshold[latLong_autoTick_values_less][latLong_autoTick_where[0]] == latLong_autoTick_threshold)[0]; #where ourselves out of the subindexing
                else:
                    # No limits if not integer, use all options
                    latLong_autoTick_where = np.where(latLong_autoTick >= latLong_autoTick_threshold)[0];
                #END IF
            #END IF
        else:
            # No limits if not integer, use all options
            latLong_autoTick_where = np.where(latLong_autoTick >= latLong_autoTick_threshold)[0];
        # END IF
        if( latLong_autoTick_where.size != 0 ):
            latLong_autoTick = latLong_autoTick_values[latLong_autoTick_where[0]]; #get the 1st match
        else:
            tickNumGoal += np.round(charEst*tickReducer); #undo the subtraction
            tickNumGoal += -np.round((charEst+3)*tickReducer); #reduce the ticks by the number of characters we're dealing with w/ included decimal+2decimalplaces
            latLong_autoTick = (np.max(latLong) - np.min(latLong))/tickNumGoal; #just goes for it if it's a super tiny range
        #END IF
    else:
        latLong_autoTick = FLG_forceTick; #forced to desired tick otherwise
    #END IF
    
    if( np.mod(np.round(np.max(latLong))-np.round(np.min(latLong)),latLong_autoTick) == 0 ):
        #if it fits perfectly, use it
        latLongAxisTicks = np.round(np.arange( np.round(np.min(latLong)) , \
                np.round(np.max(latLong)) + latLong_autoTick , \
                latLong_autoTick),2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    else:
        #uses mod to force sensible numbers
        latLongAxisTicks = np.round(np.arange( (np.round(np.min(latLong)) - np.mod(np.round(np.min(latLong)),latLong_autoTick)) , \
                (np.round(np.max(latLong)) - np.mod(np.round(np.max(latLong)),latLong_autoTick)) + latLong_autoTick , \
                latLong_autoTick),2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    #END IF
        
    #----- Axis latLong Limits -----
    #makes plotting look better
    latLongAxisLims = np.array((np.min(latLong) , np.max(latLong))); #get those xAxis limits, ensures everything is aligned
    if( np.abs(np.round(np.min(latLong)) - np.min(latLong)) < 5*10**-2 ): #fix x axis stuff (latLong)
        latLongAxisLims[0] = np.round(np.min(latLong)); #force the rounded value
    #END IF
    if( np.abs(np.round(np.max(latLong)) - np.max(latLong)) < 5*10**-2 ): #fix x axis stuff (latLong)
        latLongAxisLims[1] = np.round(np.max(latLong)); #force the rounded value
    #END IF
    if( (FLG_extendLims == True) & (np.abs(latLongAxisTicks[0]-latLongAxisLims[0]) <= latLong_autoTick/2) & (latLongAxisTicks[0] < latLongAxisLims[0]) ):
        latLongAxisLims[0] = latLongAxisTicks[0]; #set the yAxisLims to be the tick if its worth making the tick happen
    #END IF
    if( (FLG_extendLims == True) & (np.abs(latLongAxisTicks[-1]-latLongAxisLims[1]) <= latLong_autoTick/2) & (latLongAxisTicks[-1] > latLongAxisLims[1]) ):
        latLongAxisLims[1] = latLongAxisTicks[-1]; #set the yAxisLims to be the tick if its worth making the tick happen
    #END IF
    
    #----- Apply to axes if provided -----
    if( (ax is not None) & (axDir is not None) ):
        if( axDir == 'x' ):
            ax.set_xticks(latLongAxisTicks); #set x axis ticks
            ax.set_xlim(latLongAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
        else:
            ax.set_yticks(latLongAxisTicks); #set x axis ticks
            ax.set_ylim(latLongAxisLims); #set the xlims now [gotta do it after ticks in case they're a little too far]
        #END IF
        
        #----- Add units to each tick text -----
        if( FLG_labelUnits == True ):
            if( axDir == 'x' ):
                labelz = ax.get_xticklabels(); #get them labels
            else:
                labelz = ax.get_yticklabels(); #get them labels
            #END IF
            #--- edit ticks ---
            for jk in range (0,len(labelz)):
                middleLad = labelz[jk].get_text()+'°'; #add the unit on!
                labelz[jk].set_text(middleLad);
            #END FOR jk
            #--- applky ticks---
            if( axDir == 'x' ):
                ax.set_xticklabels(labelz); #apply new labels
            else:
                ax.set_yticklabels(labelz); #apply new labels
            #END IF
        #END IF
        
        #----- Remove axis labels if reqd -----
        if( FLG_removeLabels == True ):
            if( axDir == 'x' ):
                # ax.set_xticklabels([]); #remove the labels
                ax.tick_params(labelbottom=False); #remove the labels but better
            else:
                # ax.set_yticklabels([]); #remove the labels
                ax.tick_params(labelleft=False); #remove the labels but better
            #END IF
        #END IF
        
        #----- Set tick dir inward if reqd -----
        if( FLG_tickDirIn == True ):
            if( axDir == 'x' ):
                ax.tick_params(axis='x',direction='in'); #set tick dir in
            else:
                ax.tick_params(axis='y',direction='in'); #set tick dir in
            #END IF
        #END IF
    #END IF
            
    return latLongAxisTicks, latLongAxisLims
#END DEF