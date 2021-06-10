#GOAL: Plot only TEC keogram
#RD on 3/22/19
#
#INPUT: buncha stuff
#OUTPUT: no vars, just a plot is made

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

def GRITI_Mag_keo_plot(data, dates, settings):

    # #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    # if( np.isscalar(TEC_plotLimValu) == 1 ):
    #     TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    # #END IF

    #-----Plot TEC results as a Keogram-----
    #Plot just the TEC
    fig, ax = plt.subplots(); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    divider = make_axes_locatable(ax); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    
    pltHelprX, pltHelprY = np.meshgrid( data['Mag']['keo time hr'], \
                settings['Mag']['keo plot range chunks']); 
        #,vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu)
    im = ax.pcolormesh(pltHelprX, pltHelprY,  data['Mag']['keo'].T ,cmap=settings['Mag']['keo colorbar']); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    if( settings['Mag']['keo plot lim value'] != None ): 
        cax.yaxis.set_ticks(np.linspace(np.min(settings['Mag']['keo plot lim value']),np.max(settings['Mag']['keo plot lim value']),5)); #create useful tick marks
    #END IF
    # cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
    cbar.set_label(settings['Mag']['keo plot name']); #tabel the colorbar
    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
    if( settings['Mag']['keo plot lim value'] != None ): 
        #cbar.set_clim(vmin=np.min(settings['Mag']['keo plot lim value']), vmax=np.max(settings['Mag']['keo plot lim value'])); #they changed how the code works, this doesn't work anymore
        cbar.mappable.set_clim(vmin=np.min(settings['Mag']['keo plot lim value']), vmax=np.max(settings['Mag']['keo plot lim value'])); #now it's this
    #END IF
    cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
    
    #    string_Title = 'TEC Averaged on Angle of '+str(np.round(avg_anyAngle*180/np.pi,2))+' deg and Width of '+ \
    #        str(np.round(avg_anyAngle_Width,2))+' arcdeg, Avg Step # = '+str(avg_anyAngle_N)+ \
    #        ' arcdeg, Line Shows '+settings['Mag']['keo plot name']+' of Millstone Hill Zenith Beam'; #create mecha title
    if( settings['Mag']['delta method'] == 'mean' ):
        string_titleTemp = '0 Mean Mag'; #set the specific word to describe the data
    elif( settings['Mag']['delta method'] == 'savgol' ):
        string_titleTemp = 'Delta-Mag'; #set the specific word to describe the data
    else: #otherwise do nothing
        string_titleTemp = 'Mag'; #set the specific word to describe the data
    #END IF
    if( settings['Mag']['keo normalize'] == 1 ):
        string_titleTemp = string_titleTemp+' Norm\'d'; #set the specific word to describe the data
    #END IF
    string_Title = string_titleTemp+' Keo Averaged on Angle of '+str(settings['Mag']['keo angle orig'])+' deg and Width of '+ \
        str(np.round(settings['Mag']['keo width'],2))+' arcdeg'; #create mecha title
    ax.set_title(string_Title,fontproperties=settings['plot']['font title FM']); #set the title
    ax.set_xlabel('Time in UT - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
    ax.set_ylabel(settings['Mag']['keo plot range name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
    
    xAxisTicks = np.arange( (np.round(np.min(data['Mag']['keo time hr'])) - np.mod(np.round(np.min(data['Mag']['keo time hr'])),2)) , \
            (np.round(np.max(data['Mag']['keo time hr'])) - np.mod(np.round(np.max(data['Mag']['keo time hr'])),2)) + 4 , \
            4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
    ax.set_xticks(xAxisTicks); #set x axis ticks
    
    kep_range_chunks_long_plot_autoTick = (np.ceil(np.max(settings['Mag']['keo plot range chunks'])) - np.floor(np.min(settings['Mag']['keo plot range chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    if( kep_range_chunks_long_plot_autoTick > 25 ):
        kep_range_chunks_long_plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    elif( kep_range_chunks_long_plot_autoTick > 10 ):
        kep_range_chunks_long_plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    elif( kep_range_chunks_long_plot_autoTick > 5 ):
        kep_range_chunks_long_plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    elif( kep_range_chunks_long_plot_autoTick > 2 ):
        kep_range_chunks_long_plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    elif( kep_range_chunks_long_plot_autoTick > 1 ):
        kep_range_chunks_long_plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    elif( kep_range_chunks_long_plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
        kep_range_chunks_long_plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    else:
        if(settings['Mag']['keo plot range name'] == 'Latitude'): #if Y axis is latitude, use latitude
            kep_range_chunks_long_plot_autoTick = (np.max(settings['Mag']['lat range']) - np.min(settings['Mag']['lat range']))/13; #just goes for it if it's a super tiny range
        elif(settings['Mag']['keo plot range name'] == 'Longitude'): #if Y axis is longitude, use longitude
            kep_range_chunks_long_plot_autoTick = (np.max(settings['Mag']['long range']) - np.min(settings['Mag']['long range']))/13; #just goes for it if it's a super tiny range
        #END IF
    #END IF
    yAxisTicks = np.round(np.arange( np.floor(np.min(settings['Mag']['keo plot range chunks'])),np.ceil(np.max(settings['Mag']['keo plot range chunks'])),kep_range_chunks_long_plot_autoTick ),2); #creates y ticks automagically
    ax.set_yticks(yAxisTicks); #set x axis ticks
    
    # #Now drawing line of interest
    # if( settings['Mag']['keo plot range name'] == 'Longitude' ): #if true, longitude
    #     if( (np.min(plotLongRange) <= longMillstone) & (np.max(plotLongRange) >= longMillstone) ): #only plot if it's in the long range specified
    #         ax.plot( np.linspace(np.min((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),np.max((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),10,endpoint=True) , #X time hr
    #                 np.tile(longMillstone,10) , #Y latitude OR longitude arcdeg
    #                 c='xkcd:black',linewidth=1); #plots a point with a black line
    #     #END IF
    # else: #else latitude
    #     if( (np.min(plotLatRange) <= latMillstone) & (np.max(plotLatRange) >= latMillstone) ): #only plot if it's in the lat range specified
    #         ax.plot( np.linspace(np.min((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),np.max((TEC_timeUnique - np.floor(np.mean(dateRange_dayNum_zeroHr[1])))*24),10,endpoint=True) , #X time hr
    #                 np.tile(latMillstone,10) , #Y latitude OR longitude arcdeg
    #                 c='xkcd:black',linewidth=1); #plots a point with a black line
    #     #END IF
    # #END IF
    
    fig.subplots_adjust(left = 0.065, right = 0.935, top = 0.96, bottom = 0.070); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up
       
#no return, plot is the return!




