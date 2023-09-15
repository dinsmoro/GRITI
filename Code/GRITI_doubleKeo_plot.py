import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import StrMethodFormatter
from Code.GRITI_plotHelper_axisizerLatLong import GRITI_plotHelper_axisizerLatLong
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.GRITI_AMPERE_integrator_namer import GRITI_AMPERE_integrator_namer
# from Code.subfun_timeMatch import subfun_timeMatch
from Code.subfun_figFitter import figFitter
from Code.subfun_filter import subfun_filter
from Code.subfun_sunAlsoRises import sunAlsoRises
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_dictDelver import dictDelver_get
from Code.subfun_strfind import strfind
from Code.subfun_textNice import textNice
from copy import deepcopy
# from time import sleep as nappytime
        
def GRITI_doubleKeo_plot(doubleKeo, overlayName, dates, settings, dataAll=None, timeCutout = None, FLG_fancyPlot=0):
    
    if( FLG_fancyPlot >= 1 ):
        if( np.any(timeCutout == None) ):
            print('MAKING FANCY PLOT: doubleKeo_'+settings['double keo']['keo data type']+'&'+overlayName+'_keo IN fancyPlot FOLDER'); #report since you won't see anything
        else:
            print('MAKING FANCY PLOT: doubleKeo_'+settings['double keo']['keo data type']+'&'+overlayName+'_keo_cutout'+' IN fancyPlot FOLDER'); #report since you won't see anything
        #END IF
    #END IF
    
    #-----Unpack-----
    doubleKeo_keo = doubleKeo['keo'];
    doubleKeo_timeUnique = doubleKeo['time unique'];
    doubleKeo_dataRate = doubleKeo['data rate'];
    doubleKeo_overlay_integrated = doubleKeo[overlayName+' integrated'];
    # doubleKeo_overlay_time = doubleKeo[overlayName+' time'];
    doubleKeo_overlay_timeHr = doubleKeo[overlayName+' time hr'];
    #from settings now
    doubleKeo_latLong = settings['double keo']['lat long'];
    doubleKeo_latLongComb = settings['double keo']['lat long combo'];
    doubleKeo_angle = settings['double keo']['angle'];
    doubleKeo_width = settings['double keo']['width'];
    doubleKeo_alignments = settings['double keo']['plot alignments'];
    doubleKeo_plotSpacing = settings['double keo']['plot spacing'];
    doubleKeo_plotSpacingName = settings['double keo']['plot spacing name'];
    doubleKeo_niteTimes = settings['double keo']['nite times'];
    doubleKeo_arrowTimes = settings['double keo']['arrow times'];
    doubleKeo_overlay_integrateMethod = settings['double keo'][overlayName+' integrate method'];
    doubleKeo_overlay_integrateMethod_val = settings['double keo'][overlayName+' integrate val'];
    # doubleKeo_overlay_coordType = settings['double keo'][overlayName+' integrate coord type'];
    # doubleKeo_overlay_radiusNloc = settings['double keo'][overlayName+' radius n loc'];
    # doubleKeo_overlay_filtMethod = settings['double keo'][overlayName+' filter method'];
    doubleKeo_overlay_plotLim = settings['double keo'][overlayName+' plot lim'];
    doubleKeo_overlay_timeDelay = settings['double keo'][overlayName+' time delay'];
    doubleKeo_overlay_latAlign = settings['double keo'][overlayName+' alignments'];
    doubleKeo_overlay_plotLabel = settings['double keo'][overlayName+' plot label'];
    doubleKeo_overlay_plotLabel_noUnits = settings['double keo']['AMPERE plot label no units'];
    doubleKeo_extraLine = settings['double keo']['extra line'];
    doubleKeo_extraLine_delayHr = settings['double keo']['extra line delay hr'];

    if( 'degree label' in settings['map'] ):
        latlong_unitName = settings['map']['degree label']; #use supplied label
    else:
        latlong_unitName = '°';
    #END IF
    if( latlong_unitName != '' ):
        latlong_unitName_bracketed = settings['plot']['unit L']+latlong_unitName+settings['plot']['unit R']; #bracket it
    else:
        latlong_unitName_bracketed = ''; #nada
    #END IF
    doubleKeo_plotLatLong_dirAdder = [None for i in range(0,len(doubleKeo_latLong))]; #preallocate
    for i in range(0,len(doubleKeo_latLong)):
        if( 'indicate direction' in settings['map'] ):
            if( settings['map']['indicate direction'][0] == True ):
                if(doubleKeo_plotSpacingName[i] == 'Latitude'):
                    doubleKeo_plotLatLong_dirAdder[i] = settings['map']['indicate direction'][1]['lat']; #get the lat
                else:
                    doubleKeo_plotLatLong_dirAdder[i] = settings['map']['indicate direction'][1]['long']; #get the long
                #END IF
            else:
                doubleKeo_plotLatLong_dirAdder[i] = ''; #nada
            #END IF
        else:
            doubleKeo_plotLatLong_dirAdder[i] = ''; #nada
        #END IF
    #END IF
    
    

    #-----Plot TEC results as a Keogram w/ overlay integrated line as well-----
    #Prep the plot
    
    if( np.all(doubleKeo_extraLine[0] == True) ):
        FLG_extraLine = True;
        doubleKeo_extraLine = deepcopy(doubleKeo_extraLine[1:]); #cut off the front
        if( np.all(dataAll == None) ):
            print('ERROR in GRITI_doubleKeo_plot: dataAll was not provided, which should be a dict that holds the data at the path provided by doubleKeo_extraLine ('+textNice(doubleKeo_extraLine)+'). Crashing cause it\'ll crash later anyway.');
            import sys
            sys.crash();
        #END IF
        dataAlias = [None for i in range(0,len(doubleKeo_latLong))]; #preallocate
        dataAlias_timeUnique = [None for i in range(0,len(doubleKeo_latLong))]; #preallocate
        for i in range(0,len(doubleKeo_latLong)):
            extraOptions = strfind(doubleKeo_extraLine[i]['dict path'],'$|',opt=0);
            if( np.any(extraOptions > 0) ):
                FLG_splittr = True;
                kk = np.where(extraOptions)[0][-1];
                splittr_splitz = doubleKeo_extraLine[i]['dict path'][kk].split('$|'); #split the sub-dict from the sub-vect
                doubleKeo_extraLine[i]['dict path'][kk] = splittr_splitz[0]; #remove the extra bit for dict delving
                if( (splittr_splitz[1] == '0') | (splittr_splitz[1].lower() == 'x') ):
                    splittr_splitz = 0;
                    # splittr_addz = '$\mathregular{_X}$';
                elif( (splittr_splitz[1] == '1') | (splittr_splitz[1].lower() == 'y') ):
                    splittr_splitz = 1;
                    # splittr_addz = '$\mathregular{_Y}$';
                elif( (splittr_splitz[1] == '2') | (splittr_splitz[1].lower() == 'z') ):
                    splittr_splitz = 2;
                    # splittr_addz = '$\mathregular{_Z}$';
                elif( (splittr_splitz[1] == '3') | (splittr_splitz[1].lower() == 't') ):
                    splittr_splitz = 3;
                    # splittr_addz = '$\mathregular{_T}$'; #(By^2 + Bz^2)^(1/2) - not implemented yet
                #END IF
            else:
                FLG_splittr = False;
                # splittr_addz = '';
            #END IF
            dataAlias[i] = dictDelver_get(dataAll, doubleKeo_extraLine[i]['dict path'], FLG_crashOnFail=True); #get that data
            dataAlias_timeUnique_path = [doubleKeo_extraLine[i]['dict path'][0],'time unique'];
            # dataAlias_timeUnique_path.append('time unique');
            dataAlias_timeUnique[i] = dictDelver_get(dataAll, dataAlias_timeUnique_path, FLG_crashOnFail=True); #get that data
            if( FLG_splittr == True ):
                splittr_dimz = dataAlias[i].shape;
                if( splittr_dimz[0] == dataAlias_timeUnique[i].size ):
                    dataAlias[i] = dataAlias[i][:,splittr_splitz]; #automated way to get the right bit to slice on
                else:
                    dataAlias[i] = dataAlias[i][splittr_splitz,:];
                #END IF
            #END IF 
            dataAlias[i] = subfun_filter( dataAlias[i], doubleKeo_extraLine[i]['filter'], dataTime = dataAlias_timeUnique[i], dataRate = None, settings_spectra = settings['spectra'], reduceWindow = 0, FLG_reportNaNs = False); #filter same as AMPERE if AMPERE is filtered here
            if( isinstance(dataAlias[i], tuple) ): #catch time return if it happens
                dataAlias[i], dataAlias_timeUnique[i] = dataAlias[i]; #split it up
            #END IF
        #END FOR i
        FLG_extraLine_cntr = 1; #prep cntr
    else:
        FLG_extraLine = False;
    #END IF
    
    if( np.any(timeCutout == None) ):
        lineWidth = settings['plot']['line width']['regular'];
    else:
        lineWidth = settings['plot']['line width']['double plus'];
    #END IF
    
    if( FLG_fancyPlot == 0 ):
        fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(14,8.5),dpi=settings['plot']['journal dpi']); #use instead of fig because it inits an axis too (I think I dunno)
    #END IF
    axTwinx = []; #prep this list, holds the twinx axes, will grow as we plot
    caxTwinx = []; #prep
    divider = make_axes_locatable(ax); #prep to add an axis
    if( FLG_extraLine ):
        padVal = 0.15; #no right axis stuff
    else:
        padVal = 1.75; #make room for right axis stuff
    #END IF
    cax = divider.append_axes('right', size='2.0%', pad=padVal); #make a color bar axis, append [0.35]
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    for i in range(0,len(doubleKeo_latLong)):
        #-----Plot TEC results as a keogram-----
        pltHelprX, pltHelprY = np.meshgrid( (np.append(doubleKeo_timeUnique,doubleKeo_timeUnique[-1]+doubleKeo_dataRate) - dates['date range zero hr dayNum'][1]*86400)/3600, \
                    doubleKeo_plotSpacing[i]);
        im = ax.pcolormesh(pltHelprX, pltHelprY,  doubleKeo_keo[i].T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap=settings['TEC']['colormap'],zorder=i); # pseudocolor plot "stretched" to the grid
    #END FOR i
        
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar.set_label(settings['TEC']['name']+settings['TEC']['units']); #tabel the colorbar
    for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
        tick.label2.set_fontproperties(settings['plot']['font axis tick FM']); #yee
    #END FOR tick
    #cbar.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim'])); #they changed how the code works, this doesn't work anymore
    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim'])); #now it's this
    # cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
    cax.yaxis.set_ticks(np.arange(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim'])+0.1,0.1)); #create useful tick marks
    # cax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
    cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
    
    #----- Draw lines of separation between the keograms -----
    # doubleKeo_alignments_limits = np.empty(doubleKeo_alignments.size+2); #preallocate
    for i in range(0,doubleKeo_alignments.size):
        ax.plot( np.array( (np.min((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600)) ), \
                np.ones(2)*doubleKeo_alignments[i],
                c='xkcd:white',linewidth=settings['plot']['line width']['thicc'],linestyle='-'); #plots a point with a black line
        # plt.pause(0.0001); #another pause needed to make matplotlib work [really broke ioff plots]
        # nappytime(0.0001); #another pause needed to make matplotlib work
        plt.draw();
        # doubleKeo_alignments_limits[i+1] = ax.transLimits.transform((1,doubleKeo_alignments[i]))[1]; #get the position of the split in axis units [this was NOT working on ioff plots, barely on ion]
    #END FOR i
    # doubleKeo_alignments_limits[0] = 1; #these are known
    # doubleKeo_alignments_limits[-1] = 0; #these are known
    #manual workaround
    doubleKeo_alignments_limits = np.concatenate( ([np.max(doubleKeo_latLongComb),],doubleKeo_alignments,[np.min(doubleKeo_latLongComb),]) )-np.min(doubleKeo_latLongComb);
    doubleKeo_alignments_limits = doubleKeo_alignments_limits/np.max(doubleKeo_alignments_limits);
    
    #!! hindsight: solve this with 2 separate plots that have no spacing between them !!    
    for i in range(0,len(doubleKeo_latLong)):
        # -----Plot overlay results as a 1D line -----
        axTwinx.append(ax.twinx()); #add on the new twinx axis
        doubleKeo_overlay_integrated_plot = np.copy(doubleKeo_overlay_integrated[i]); #copy
        if( doubleKeo_overlay_plotLim != False ):
            doubleKeo_overlay_integrated_plot[ doubleKeo_overlay_integrated_plot > doubleKeo_overlay_plotLim ] = doubleKeo_overlay_plotLim; #cap the max val
        #END IF
        axTwinx[-1].plot( doubleKeo_overlay_timeHr[i], doubleKeo_overlay_integrated_plot , linewidth=lineWidth , color='#B129FF',zorder=i); #plot xkcd:violet
        if( (i == len(doubleKeo_latLong)-1) & (FLG_extraLine == False) ):
            axTwinx[-1].set_ylabel(doubleKeo_overlay_plotLabel,fontproperties=settings['plot']['font axis label FM']); #set the y axis label
        #END IF
        
        #now adjust the ylim so the lines are about where we want them to be
        # #--- Adjust yLim max to be at the right place --- (turns out can't do both at once since they're linked - so gotta... solve)
        # twinx_yFactor = 1 - (doubleKeo_alignments_limits[i] - axTwinx[-1].transLimits.transform((1,np.max(doubleKeo_overlay_integrated[i])))[1]); #get a factor
        # twinx_yLims = np.array( axTwinx[-1].get_ylim() ); #get the ylims
        # twinx_yLims[1] = np.max(doubleKeo_overlay_integrated[i])/doubleKeo_alignments_limits[i]-(twinx_yLims[0]*(1-doubleKeo_alignments_limits[i]))*2; #solving this is nasty I'm not big brain here
        # axTwinx[-1].set_ylim(twinx_yLims); #apply the new ylim because it adjusts the next thing
        
        # #--- Adjust yLim min to align the data on the line of interest ---
        # lineOfInterestLoc = (doubleKeo_alignments_limits[i]-doubleKeo_alignments_limits[i+1])*(doubleKeo_overlay_latAlign[i]-np.min(doubleKeo_latLongComb[i]))/(np.max(doubleKeo_latLongComb[i])-np.min(doubleKeo_latLongComb[i]))+doubleKeo_alignments_limits[i+1]; #get line of interest location on TEC keogram in a form of 1 to 0 (scaled by the keogram plot-within-plots)
        # twinx_valAtLineOfInterest = np.mean(doubleKeo_overlay_integrated[i]) + 1.5*np.std(doubleKeo_overlay_integrated[i]); #this is the value we want to align with the line of interest
        # twinx_minToAlign = (lineOfInterestLoc*twinx_yLims[1] - twinx_valAtLineOfInterest)/(lineOfInterestLoc - 1); #calc the min to get that value at the LOI
        # axTwinx[-1].set_ylim( twinx_minToAlign , twinx_yLims[1] ); #set y axis limits
        
        #--- I had to do actual algebra for this (combined two equations above) ---
        twinx_yLims = np.array( axTwinx[-1].get_ylim() ); #get the ylims
        lim1 = doubleKeo_alignments_limits[i]; #get the upper limit in axis units
        yMax = np.nanmax(doubleKeo_overlay_integrated_plot); #get the desired yMax
        # yMin = twinx_yLims[0]; #get the current yMin
        lim2 = (doubleKeo_alignments_limits[i]-doubleKeo_alignments_limits[i+1])*(doubleKeo_overlay_latAlign[i]-np.min(doubleKeo_latLongComb[i]))/(np.max(doubleKeo_latLongComb[i])-np.min(doubleKeo_latLongComb[i]))+doubleKeo_alignments_limits[i+1]; #get line of interest location on TEC keogram in a form of 1 to 0 (scaled by the keogram plot-within-plots)
        yGoal = np.nanmean(doubleKeo_overlay_integrated_plot) + 1.5*np.nanstd(doubleKeo_overlay_integrated_plot); #this is the value we want to align with the line of interest
        
        #--- Calc the yLim Max and Min ---
        twinx_yLims = np.array( axTwinx[-1].get_ylim() ); #get the ylims
        twinx_yLims[1] = (yMax/lim1 + 2*yGoal/(lim2-1)*(1-lim1))/(1+2*lim2/(lim2-1)*(1-lim1)); #calc max limit
        twinx_yLims[0] = (lim2*twinx_yLims[1] - yGoal)/(lim2-1); #calc min limit
        if( twinx_yLims[0] < np.max(doubleKeo_overlay_integrated_plot) ):
            axTwinx[-1].set_ylim( twinx_yLims ); #set y axis limits
        #END IF
        # #alignment check
        # ax.transLimits.transform((1,doubleKeo_overlay_latAlign[i]))[1]
        # axTwinx[-1].transLimits.transform((1,twinx_valAtLineOfInterest))[1]
        
        #this is to keep plotting similar
        dividerTwinx = make_axes_locatable(axTwinx[-1]); #prep to add an axis
        caxTwinx.append(dividerTwinx.append_axes('right', size='2.0%', pad=padVal)); #make a color bar axis, append
        caxTwinx[-1].set_visible(False); #mkae it invisible so it matches the other plots in width
        
        # xAxisTicksStep = 4; #hr, hour steps between each x axis tick
        # xAxisTicks = np.arange( (np.round((np.min(doubleKeo_timeUnique)-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((np.min(doubleKeo_timeUnique)-dates['date range zero hr dayNum'][1]*86400)/3600), xAxisTicksStep)) , \
        #     (np.round((np.max(doubleKeo_timeUnique)-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((np.max(doubleKeo_timeUnique)-dates['date range zero hr dayNum'][1]*86400)/3600),xAxisTicksStep)) + xAxisTicksStep , \
        #     xAxisTicksStep); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
        # ax.set_xticks(xAxisTicks); #set x axis ticks
        # ax.set_xlim( np.min(xAxisTicks) , np.max(xAxisTicks) ); #set x axis limits
        # axTwinx[-1].set_xlim( np.min(xAxisTicks) , np.max(xAxisTicks) ); #set x axis limits
        
        if( FLG_fancyPlot == 0 ):
            if( np.any(timeCutout == None) ):
                GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1]); #automagic time ticks here
            else:
                GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1],FLG_manualLims=timeCutout); #automagic time ticks here
            #END IF
        else:
            if( np.any(timeCutout == None) ):
                GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1],tickNumGoal=22); #automagic time ticks here
            else:
                GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1],FLG_manualLims=timeCutout,tickNumGoal=22); #automagic time ticks here
            #END IF
        #END IF
        
        if( FLG_extraLine == False ):
            #adjust the tick marks
            # tickzLim = 40*np.median(doubleKeo_overlay_integrated[i]); #get the median*40 of the integrated JH
            tickz = axTwinx[-1].get_yticks(); #get the ticks
            # tickz = tickz[ (tickz < tickzLim) & (tickz >= 0) ]; #get only the tickz we want
            tickzPos = axTwinx[-1].transLimits.transform(np.vstack( (np.ones(tickz.size), tickz) ).T)[:,1]; #get the tickz positions on the plot
            tickz = tickz[ (tickz >= 0) & (doubleKeo_alignments_limits[i] >= tickzPos) & (doubleKeo_alignments_limits[i+1] <= tickzPos) ]; #make sure the ticks fall within the keogram area and they're positive
            if( (tickz.size == 1) & (doubleKeo_overlay_plotLim != False) ):
                tickz = np.append(tickz,np.floor(np.max(doubleKeo_overlay_integrated_plot)/10**np.int64(np.log10(np.max(doubleKeo_overlay_integrated_plot))))*10**np.int64(np.log10(np.max(doubleKeo_overlay_integrated_plot)))); #make sure it's not just 0 on the axis labels
            #END IF
            axTwinx[-1].set_yticks(tickz); #set the ticks
        else:
            tickz = axTwinx[-1].get_yticklabels(); #get the tick labels
            for jj in range(0,len(tickz)):
                tickz[jj].set_text(''); #remove the label
            #END FOR jj
            axTwinx[-1].set_yticklabels(tickz); #set the tick labels
            axTwinx[-1].tick_params(top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False); #remove the tick marks
        #END IF
        
        if( FLG_extraLine == True ):
            #more plotting to do
            # -----Plot overlay results as a 1D line -----
            axTwinx.append(ax.twinx()); #add on the new twinx axis
            # doubleKeo_overlay_integrated_plot = np.copy(dataAlias[i]); #copy
            # if( doubleKeo_overlay_plotLim != False ):
            #     doubleKeo_overlay_integrated_plot[ doubleKeo_overlay_integrated_plot > doubleKeo_overlay_plotLim ] = doubleKeo_overlay_plotLim; #cap the max val
            # #END IF doubleKeo_overlay_timeDelay[i]-38/60
            axTwinx[-1].plot( (dataAlias_timeUnique[i]-dates['date range zero hr dayNum'][1]*86400)/3600+doubleKeo_extraLine_delayHr[i], dataAlias[i] , linewidth=settings['plot']['line width']['double plus'], linestyle=settings['plot']['line style'][FLG_extraLine_cntr] , color='xkcd:gunmetal',zorder=i); #plot xkcd:violet
            
            #--- I had to do actual algebra for this (combined two equations above) ---
            twinx_yLims = np.array( axTwinx[-1].get_ylim() ); #get the ylims
            lim1 = doubleKeo_alignments_limits[i]; #get the upper limit in axis units
            yMax = np.nanmax(dataAlias[i]); #get the desired yMax
            # yMin = twinx_yLims[0]; #get the current yMin
            lim2 = (doubleKeo_alignments_limits[i]-doubleKeo_alignments_limits[i+1])*(doubleKeo_overlay_latAlign[i]-np.min(doubleKeo_latLongComb[i]))/(np.max(doubleKeo_latLongComb[i])-np.min(doubleKeo_latLongComb[i]))+doubleKeo_alignments_limits[i+1]; #get line of interest location on TEC keogram in a form of 1 to 0 (scaled by the keogram plot-within-plots)
            yGoal = np.nanmean(dataAlias[i]) + 1.5*np.nanstd(dataAlias[i]); #this is the value we want to align with the line of interest
            
            #--- Calc the yLim Max and Min ---
            twinx_yLims = np.array( axTwinx[-1].get_ylim() ); #get the ylims
            twinx_yLims[1] = (yMax/lim1 + 2*yGoal/(lim2-1)*(1-lim1))/(1+2*lim2/(lim2-1)*(1-lim1)); #calc max limit
            twinx_yLims[0] = (lim2*twinx_yLims[1] - yGoal)/(lim2-1); #calc min limit
            if( twinx_yLims[0] < np.nanmax(dataAlias[i]) ):
                axTwinx[-1].set_ylim( twinx_yLims ); #set y axis limits
            #END IF
            
            #this is to keep plotting similar
            dividerTwinx = make_axes_locatable(axTwinx[-1]); #prep to add an axis
            caxTwinx.append(dividerTwinx.append_axes('right', size='2.0%', pad=padVal)); #make a color bar axis, append
            caxTwinx[-1].set_visible(False); #mkae it invisible so it matches the other plots in width
            
            if( FLG_fancyPlot == 0 ):
                if( np.any(timeCutout == None) ):
                    GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1]); #automagic time ticks here
                else:
                    GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1],FLG_manualLims=timeCutout); #automagic time ticks here
                #END IF
            else:
                if( np.any(timeCutout == None) ):
                    GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1],tickNumGoal=22); #automagic time ticks here
                else:
                    GRITI_plotHelper_axisizerTime((doubleKeo_timeUnique-dates['date range zero hr dayNum'][1]*86400)/3600,ax=axTwinx[-1],FLG_manualLims=timeCutout,tickNumGoal=22); #automagic time ticks here
                #END IF
            #END IF
            
            # #adjust the tick marks
            # tickz = axTwinx[-1].get_yticklabels(); #get the tick labels
            # for jj in range(0,len(tickz)):
            #     tickz[jj].set_text(''); #remove the label
            # #END FOR jj
            # axTwinx[-1].set_yticklabels(tickz); #set the tick labels
            axTwinx[-1].tick_params(top=False, bottom=False, left=False, right=False,
                labelright=False, labelbottom=False); #remove the tick marks
            
            # FLG_extraLine_cntr += 1; #increment cntr
        #END IF
        axTwinx[-1].format_coord = plot_cursor_relabeler(axTwinx[-1], ax); #hopefully make it work
        

        #----- Plot Shading to represent 'night' -----
        if( doubleKeo_niteTimes[0] == False ):
            if( doubleKeo_plotSpacingName[i] == 'Latitude' ): #if true, latitude
                if( doubleKeo_overlay_latAlign[i] > 45 ):
                    latToUse = 45; #cap at 50 deg lat to keep it from getting too zesty at the poles
                else:
                    latToUse = doubleKeo_overlay_latAlign[i];
                #END IF
                (doubleKeo_niteTimes_sunRise, doubleKeo_niteTimes_sunSet, doubleKeo_dateRange_fullPad) = sunAlsoRises(dates['date range full'],latToUse,np.mean(doubleKeo_latLong[i][1])); #call sunrise/set function
            else:
                (doubleKeo_niteTimes_sunRise, doubleKeo_niteTimes_sunSet, doubleKeo_dateRange_fullPad) = sunAlsoRises(dates['date range full'],np.mean(doubleKeo_latLong[i][0]),doubleKeo_overlay_latAlign[i]); #call sunrise/set function
            #END IF
            doubleKeo_dateRange_dayNum_fullPad = subfun_date_to_dayNum(doubleKeo_dateRange_fullPad); #convert to dayNum
            doubleKeo_niteTimes_sunRise = (doubleKeo_niteTimes_sunRise + doubleKeo_dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*24; #hrs, center around zero hr and convert ot hrs
            doubleKeo_niteTimes_sunSet = (doubleKeo_niteTimes_sunSet + doubleKeo_dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*24; #hrs, center around zero hr and convert ot hrs
            doubleKeo_niteTimes_sunRise = doubleKeo_niteTimes_sunRise[1:]; #remove 1st
            doubleKeo_niteTimes_sunSet = doubleKeo_niteTimes_sunSet[:-1]; #remove last
            #FIFTH STEP: PLOT THIS STUFF
            for j in range(0,doubleKeo_niteTimes_sunSet.size):
                # ax.axvspan(doubleKeo_overlay_niteTimes_array[j,0], doubleKeo_overlay_niteTimes_array[j,1], alpha=0.25, color='xkcd:black');
                if(doubleKeo_plotSpacingName[i] == 'Latitude'): #if Y axis is latitude, use latitude
                    recta = Rectangle((doubleKeo_niteTimes_sunSet[j], np.min(doubleKeo_latLong[i][0])), doubleKeo_niteTimes_sunRise[j]-doubleKeo_niteTimes_sunSet[j], np.diff(doubleKeo_latLong[i][0]).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.25); #make that patch
                    ax.add_patch(recta); #add on that patch
                else: #otherwise longitude
                    recta = Rectangle((doubleKeo_niteTimes_sunSet[j], np.min(doubleKeo_latLong[i][1])), doubleKeo_niteTimes_sunRise[j]-doubleKeo_niteTimes_sunSet[j], np.diff(doubleKeo_latLong[i][1]).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.25); #make that patch
                    ax.add_patch(recta); #add on that patch
                #END IF
            #END FOR j
        else:
            for j in range(0,dates['date range zero hr hours'].size):
                doubleKeo_niteTimes_per = dates['date range zero hr hours'][j] + np.asarray(doubleKeo_niteTimes[i]); #get the nite time ranges
                if(doubleKeo_plotSpacingName[i] == 'Latitude'): #if Y axis is latitude, use latitude
                    recta = Rectangle((doubleKeo_niteTimes_per[0], np.min(doubleKeo_latLong[i][0])), np.diff(doubleKeo_niteTimes_per).item(), np.diff(doubleKeo_latLong[i][0]).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.25); #make that patch
                    ax.add_patch(recta); #add on that patch
                else: #otherwise longitude
                    recta = Rectangle((doubleKeo_niteTimes_per[0], np.min(doubleKeo_latLong[i][1])), np.diff(doubleKeo_niteTimes_per).item(), np.diff(doubleKeo_latLong[i][1]).item(), edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.25); #make that patch
                    ax.add_patch(recta); #add on that patch
                #END IF
            #END FOR j
        #END IF
        
        #----- Now drawing line of interest -----
        if( doubleKeo_plotSpacingName[i] == 'Latitude' ): #if true, latitude
            if( (np.min(doubleKeo_latLong[i][0]) <= doubleKeo_overlay_latAlign[i]) & (np.max(doubleKeo_latLong[i][0]) >= doubleKeo_overlay_latAlign[i]) ): #only plot if it's in the lat range specified
                ax.plot( np.linspace(np.min((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(doubleKeo_overlay_latAlign[i],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=settings['plot']['line width']['smol']); #plots a point with a black line
            #END IF
        else:
            if( (np.min(doubleKeo_latLong[i][1]) <= doubleKeo_overlay_latAlign[i]) & (np.max(doubleKeo_latLong[i][1]) >= doubleKeo_overlay_latAlign[i]) ): #only plot if it's in the lat range specified
                ax.plot( np.linspace(np.min((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                    np.tile(doubleKeo_overlay_latAlign[i],10) , #Y latitude OR longitude arcdeg
                    c='xkcd:black',linewidth=settings['plot']['line width']['smol']); #plots a point with a black line
            #END IF
        #END IF
    #END FOR i
    
    #----- Label and make it look nice -----
    if( FLG_fancyPlot == 0 ): #different aspect ratios require different spacing assumptions
        string_title = textNice(len(doubleKeo_latLong))+' '+settings['TEC']['name']+' Keograms on '+textNice(np.round(doubleKeo_angle[i],2))+' ° Angle & Widths of '+ \
            textNice(np.round(doubleKeo_width,2))+' arc° \n Integrated '+doubleKeo_overlay_plotLabel_noUnits; #create mecha title
        string_title = string_title + GRITI_AMPERE_integrator_namer(doubleKeo_overlay_integrateMethod[i], doubleKeo_overlay_integrateMethod_val[i], plotLatRange=doubleKeo_latLong[i][0], FLG_partial=True); #if new integration methods added only need to adjust the function now so smart (lol no I ain't using a method)
        string_title = string_title + ' w/ time delays '+textNice(doubleKeo_overlay_timeDelay)+' hrs';
        ax.set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
    #END IF
    if( i == len(doubleKeo_latLong)-1 ):
        ax.set_xlabel('Time in UT - 0 Hr on Day '+textNice(dates['date range zero hr dayNum'][1])+', '+textNice(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
    #END IF
    ax.set_ylabel(doubleKeo_plotSpacingName[i]+' '+latlong_unitName_bracketed+doubleKeo_plotLatLong_dirAdder[i],fontproperties=settings['plot']['font axis label FM']); #set the y axis label
    
    # autoTick = (np.ceil(np.max(doubleKeo_plotSpacing)) - np.floor(np.min(doubleKeo_plotSpacing)))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
    # if( autoTick > 25 ):
    #     autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
    # elif( autoTick > 10 ):
    #     autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
    # elif( autoTick > 5 ):
    #     autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
    # elif( autoTick > 2 ):
    #     autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
    # elif( autoTick > 1 ):
    #     autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
    # elif( autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    #     autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
    # else:
    #     autoTick = (np.max(doubleKeo_latLongComb) - np.min(doubleKeo_latLongComb))/13; #just goes for it if it's a super tiny range
    # #END IF
    # yAxisTicks = np.round(np.arange( np.floor(np.min(doubleKeo_plotSpacing)),np.ceil(np.max(doubleKeo_plotSpacing))+autoTick,autoTick ),2); #creates y ticks automagically
    # ax.set_yticks(yAxisTicks); #set x axis ticks
    
    if( FLG_fancyPlot == 0 ): #different aspect ratios require different spacing assumptions
        GRITI_plotHelper_axisizerLatLong(doubleKeo_plotSpacing,ax=ax,axDir='y',tickNumGoal=17);
    else:
        GRITI_plotHelper_axisizerLatLong(doubleKeo_plotSpacing,ax=ax,axDir='y',tickNumGoal=13);
    #END IF
    
    #-----Draw line from 1st TEC event to 2nd-----
    if( np.any(timeCutout == None) ):
        if( np.any(doubleKeo_arrowTimes[0] != False) ):
            for jj in range(0,len(doubleKeo_arrowTimes)):
                if( (doubleKeo_arrowTimes[jj][0] >= np.min((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600)) & \
                    (doubleKeo_arrowTimes[jj][1] <= np.max((doubleKeo_timeUnique - dates['date range zero hr dayNum'][1]*86400)/3600)) ):
                    
                    con = ConnectionPatch(xyA=(doubleKeo_arrowTimes[jj][0],doubleKeo_overlay_latAlign[0]), coordsA=ax.transData, #-11.81 #-10.54
                                          xyB=(doubleKeo_arrowTimes[jj][1],doubleKeo_overlay_latAlign[1]), coordsB=ax.transData, #-11.13 #-9.96
                                          arrowstyle="-|>", shrinkA=5, shrinkB=5,mutation_scale=20, fc="xkcd:pink",
                                          color='xkcd:pink', linewidth=settings['plot']['line width']['plus']); #prep a line between plots
                    fig.add_artist(con); #draw the line            
                #END IF
            #END FOR jj
        #END IF
    else:
        if( np.any(doubleKeo_arrowTimes[0] != False) ):
            for jj in range(0,len(doubleKeo_arrowTimes)):
                if( (doubleKeo_arrowTimes[jj][0] >= np.min(timeCutout)) & \
                    (doubleKeo_arrowTimes[jj][1] <= np.max(timeCutout)) ):
                    
                    con = ConnectionPatch(xyA=(doubleKeo_arrowTimes[jj][0],doubleKeo_overlay_latAlign[0]), coordsA=ax.transData, #-11.81 #-10.54
                                          xyB=(doubleKeo_arrowTimes[jj][1],doubleKeo_overlay_latAlign[1]), coordsB=ax.transData, #-11.13 #-9.96
                                          arrowstyle="-|>", shrinkA=5, shrinkB=5,mutation_scale=20, fc="xkcd:pink",
                                          color='xkcd:pink', linewidth=settings['plot']['line width']['plus']); #prep a line between plots
                    fig.add_artist(con); #draw the line   
                #END IF
            #END FOR jj
        #END IF
    #END IF
    
    figFitter(fig); #fit the fig fast
    if( FLG_fancyPlot != 0 ):
        if( np.any(timeCutout == None) ):
            fig.savefig(os.path.join(settings['paths']['fancyPlots'],'doubleKeo_'+settings['double keo']['keo data type']+'&'+overlayName+'_keo'+settings['plot']['save file type'])); #save the figure
        else:
            fig.savefig(os.path.join(settings['paths']['fancyPlots'],'doubleKeo_'+settings['double keo']['keo data type']+'&'+overlayName+'_keo_cutout'+settings['plot']['save file type'])); #save the figure
        #END IF
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF

def plot_cursor_relabeler(ax_top, ax_under): #inspired by https://stackoverflow.com/a/21585524/2403531
    #ax_under needs to be a list, even if just [ax0]
    #ax_names should be in the order of [top, most under, less under, least under] (e.g. 1 longer than ax_under)
    def format_coord(x, y):
        # x, y are data coordinates, they appear mystically
        # returnString = '';
        display_coord = ax_top.transData.transform((x,y)); # convert to display coords
        ax_coord = ax_under.transData.inverted().transform(display_coord); #fire up the inverter & convert back to data coords with respect to ax
        # coords = []; #prep list
        # for j in range(0,len(ax_under)):
        #     ax_coord = ax_under[j].transData.inverted().transform(display_coord); #fire up the inverter & convert back to data coords with respect to ax
        #     print(str(ax_coord))
        #     print(str(len(ax_under)))
        #     print(str(ax_names))
        #     # coords.append(ax_coord); #create list of coords
        #     returnString = returnString + ax_names[j]+' ({:.3f}, {:.3f}) | '.format(ax_coord[0], ax_coord[1]); #tack it on
        # #END FOR j
        # returnString = ' | '+ax_names[-1]+': ({:.3f}, {:.3f})'.format(x, y); #end of the string        
        # return (returnString)
        return ('time={:.3f}, lat/long={:.3f}'.format(ax_coord[0], ax_coord[1]))
    #END DEF
    return format_coord
#END DEF