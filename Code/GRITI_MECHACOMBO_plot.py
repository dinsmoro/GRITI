"""
GOAL: Plot Kp and OMNI and SuperMAG
RD on 5/11/19

INPUT: buncha Kp stuff
OUTPUT: no vars, just a plot is made

opt: 0 (default) - plots regular
    1 - plots with vertical lines dividing the days
    2 - plots with vertical lines dividing the days, and the dates written just above the bottom axis
    3 - does the same as above, but uses day number format instead of dates
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import copy
from Code.subfun_figFitter import figFitter
from Code.GRITI_import_SuperMAG_stations import GRITI_import_SuperMAG_stations
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.subfun_strstr import strstrNB as strstr
from Code.subfun_strfind import strfind
from Code.subfun_filter import subfun_filter


def GRITI_MECHACOMBO_plot(FLG_MECHACOMBO_names, data, dates, settings, opt=0, FLG_ignoreUnitMismatch=False, FLG_MECHACOMBO_plot_timeLim=[False], FLG_MECHACOMBO_plot_localTime=False, FLG_MECHACOMBO_plot_singleColumn=False, FLG_fancyPlot=0):
    #----Start plotting-----
    if( FLG_fancyPlot >= 1 ):
        print('MAKING FANCY PLOT: MECHACOMBO_fancyPlot IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    
    #Unpack major sub-stuff
    settings_plot = settings['plot'];
    settings_paths = settings['paths'];
    
    #Unpack everything but line widths
    journal_dpi = settings_plot['journal dpi'];
    FONT_grandioseFM = settings_plot['font grandiose FM'];
    FONT_axisLabelFM = settings_plot['font axis label FM'];
    
    #Unpack line widths
    # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    
    if( FLG_MECHACOMBO_plot_localTime == True ):
        #do the timezone calcs
        #make it local time, some legwork required
        import timezonefinder
        from datetime import datetime
        import pytz
        from Code.subfun_textNice import textNice
        
        if( settings['map']['coord type'] == 'geo' ):
            localLat = settings['map']['site coords'][0][0];
            localLong = settings['map']['site coords'][0][1];
        else:
            localLat = settings['map']['site coords geo'][0][0];
            localLong = settings['map']['site coords geo'][0][1];
        #END IF
        
        #--- Find time zone of lat and long ---
        tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
        dayNite_timeZoneID = tf.certain_timezone_at(lat=localLat, lng=localLong); #use it to find the time zone
        if dayNite_timeZoneID is None:
            from urllib.request import urlopen
            #use geonames site as a backup
            url = 'http://api.geonames.org/timezone?lat='+str(localLat)+'&lng='+str(localLong)+'&username=razzluhdzuul'; #create link for lat/long
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
        timeZoneObj_zeroHr = timeZoneObj.localize(datetime.strptime(str(dates['date range zero hr']), '[%Y\t%m\t%d]')); #time zone info at zero hr
        dayNite_DSToffset_str = timeZoneObj_zeroHr.strftime('%z').replace('0',''); #remove the 0 that was extraneous
        dayNite_DSToffset = np.int64(dayNite_DSToffset_str); #get the number version
        dayNite_timeZoneName = timeZoneObj_zeroHr.tzname(); #get the time zone name (like 'EST' or 'EDT' depending on standard or daylight savings time)
        dateNite_DSTnUTCOffset = timeZoneObj_zeroHr.dst().total_seconds()/3600; #get the time offset
        if( np.mod(dateNite_DSTnUTCOffset,1) == 0 ):
            dateNite_DSTnUTCOffset = np.int64(dateNite_DSTnUTCOffset); #convert to integer
        #END IF
    else:
        dayNite_DSToffset = 0; #set so things work
    #END IF
    if( FLG_MECHACOMBO_plot_timeLim[0] == None ):
        FLG_manualLims = dates['date range zero hr hour bounds'];
    else:
        FLG_manualLims = np.asarray(FLG_MECHACOMBO_plot_timeLim); #set directly
        if( FLG_MECHACOMBO_plot_localTime == True ):
            FLG_manualLims -= dayNite_DSToffset; #adjust to UTC b/c these are assumed to be in local time if local time is on
        #END IF
    #END IF
    
    #Prime the plot
    if( FLG_fancyPlot == 0 ):
        fig = plt.figure();
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig = plt.figure(figsize=(14,14.5),dpi=journal_dpi);
    #END IF
    #Make sub-axes
    #various supported ez splits
    FLG_sameLineUnits = True; #show them units on the same line
    if( (FLG_MECHACOMBO_plot_singleColumn == False) & (len(FLG_MECHACOMBO_names) == 6) ):
        # letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
        # letteringPositionY = 0.90; #set the X position of the lettering (e.g., a. b. c. ...)
        plt_nrows = 3;
        plt_ncols = 2;
        plt_topIndexr = [0,1];
        plt_bottomIndexr = [4,5];
    elif( (FLG_MECHACOMBO_plot_singleColumn == False) & (len(FLG_MECHACOMBO_names) == 8) ):
        # letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
        # letteringPositionY = 0.90; #set the X position of the lettering (e.g., a. b. c. ...)
        plt_nrows = 4;
        plt_ncols = 2;
        plt_topIndexr = [0,1];
        plt_bottomIndexr = [6,7];
    elif( (FLG_MECHACOMBO_plot_singleColumn == False) & (len(FLG_MECHACOMBO_names) == 10) ):
        # letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
        # letteringPositionY = 0.90; #set the X position of the lettering (e.g., a. b. c. ...)
        plt_nrows = 5;
        plt_ncols = 2;
        plt_topIndexr = [0,1];
        plt_bottomIndexr = [8,9];
    elif( (FLG_MECHACOMBO_plot_singleColumn == False) & (len(FLG_MECHACOMBO_names) == 12) ):
        # letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
        # letteringPositionY = 0.90; #set the X position of the lettering (e.g., a. b. c. ...)
        plt_nrows = 6;
        plt_ncols = 2;
        plt_topIndexr = [0,1];
        plt_bottomIndexr = [5,11];
        FLG_sameLineUnits = False; #too dense for same line units
    else:
        if( len(FLG_MECHACOMBO_names) > 5 ):
            #places it inside the plot
            # letteringPositionX = 0.01; #set the X position of the lettering (e.g., a. b. c. ...)
            # letteringPositionY = 0.85; #set the X position of the lettering (e.g., a. b. c. ...)
            FLG_sameLineUnits = False; #too dense for same line units
        # else:
        #     letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
        #     letteringPositionY = 0.90; #set the X position of the lettering (e.g., a. b. c. ...)
        #END IF
        plt_nrows = len(FLG_MECHACOMBO_names);
        plt_ncols = 1;
        plt_topIndexr = [0];
        plt_bottomIndexr = [len(FLG_MECHACOMBO_names)-1];
    #END IF
    gridr = gridspec.GridSpec(nrows=plt_nrows, ncols=plt_ncols, figure=fig); #make a grid (used in case need larger than 1 "grid" plots)
    gridr.update(hspace=0.05,wspace=0.30); # set the spacing between axes. 
    # fig.add_subplot(gridr[0:1]);
    for i in range(0,len(FLG_MECHACOMBO_names),1):
        fig.add_subplot(gridr[i:i+1]); #dynamically add axes
    #END FOR i
    ax = fig.axes; #get a list of the axes
    
    #Plot each part needed
    secretLabel = []; #to hold secret labels
    secretLabel_string = [];
    secretLabel_color = []; 
    OMNI_keyz = list(data['OMNI'].keys()); #get OMNI keys
    SuperMAG_keyz = list(data['SuperMAG'].keys()); #get SuperMAG keys
    try:
        MagCAN_keyz = list(data['MagCAN'].keys()); #get MagCAN keys
    except:
        MagCAN_keyz = []; #empty list if MagCAN isn't available
    #END TRY
    for i in range(0,len(FLG_MECHACOMBO_names)):
        ax[i].set_aspect('auto'); #set to auto for all axes
        
        #this and the extra loop allow for stuff to be plotted on the same plot
        if( isinstance(FLG_MECHACOMBO_names[i],list) | isinstance(FLG_MECHACOMBO_names[i],tuple) ):
            plotNames_now = copy.deepcopy(FLG_MECHACOMBO_names[i]); #use it
        else:
            plotNames_now = [FLG_MECHACOMBO_names[i]]; #wrap it
        #END IF
        plotNames_units = [[] for j in range(0,len(plotNames_now))]; #preallocate
        plotNames_imHandles = [[] for j in range(0,len(plotNames_now))]; #preallocate
        for j in range(0,len(plotNames_now)):
            #catch decorators
            FLG_use_AMPERE = None; #activates AMPERE data type if not none
            if( (strfind(plotNames_now[j].lower(),'amp:',opt=1) > 0) | (strfind(plotNames_now[j].lower(),'ampere:',opt=1) > 0) ):
                if( (strfind(plotNames_now[j].lower(),':int:',opt=1) > 0) | (strfind(plotNames_now[j].lower(),':integrate:',opt=1) > 0) ):
                    plotNames_now[j] = plotNames_now[j][plotNames_now[j].rfind(':')+1:]; #remove the decorators
                    from Code.GRITI_AMPERE_integrator import GRITI_AMPERE_integrator
                    settings_AMPERE_copy = copy.deepcopy(settings['AMPERE']);
                    settings_AMPERE_copy['data type'] = plotNames_now[j];
                    FLG_use_AMPERE = GRITI_AMPERE_integrator(data['AMPERE'], dates, settings_AMPERE_copy, settings['map']['lat range'], settings['map']['long range'], settings_AMPERE_copy['integrate method'], settings_AMPERE_copy['integrate method lat val'], 
                                                AMPERE_integrateMethod_coordType=settings_AMPERE_copy['integrate method coord type'], AMPERE_integrateMethod_coordType_global=settings['map']['coord type'], GRITI_import_AMPERE=None, AMPERE_desired_latLongSteps=settings['AMPERE']['lat long steps'],
                                                AMPERE_import_AMPERE_hemi=settings['map']['hemi'], AMPERE_import_AMPERE_login=settings['config']['login AMPERE'], settings_paths=settings['paths'],
                                                AMPERE_integrateMethod_log=settings_AMPERE_copy['integrate method log'], AMPERE_integrateMethod_radiusLoc=settings_AMPERE_copy['integrate method radius n loc'][1], AMPERE_integrateMethod_radius=settings_AMPERE_copy['integrate method radius n loc'][0]);
                #END IF
            #END IF            
            
            #prep the plot name
            requestLoc = plotNames_now[j].find(' &R '); #get request loc (for data import stuff)
            filtLoc = plotNames_now[j].find(' &F '); #get filt loc (goes to filter)
            
            if( requestLoc > -1 ):
                requestMethod = plotNames_now[j].split(' &R ')[1]; #get mostly just the request string
                if( filtLoc > requestLoc ):
                    requestMethod = requestMethod[:requestMethod.find( ' &F ')]; #remove any filter stuff after
                #END IF
            else:
                requestMethod = ''; #more general so make it clear abs nothing here
            #END IF
            if( filtLoc > -1 ):
                filtMethod = plotNames_now[j].split(' &F ')[1]; #get mostly just the request string
                if( requestLoc > filtLoc ):
                    filtMethod = filtMethod[:requestMethod.find( ' &R ')]; #remove any filter stuff after
                #END IF
            else:
                filtMethod = ''; #filt alg works with empty strings just fine
            #END IF
            if( (filtLoc == -1) & (requestLoc == -1) ):
                pass; #saves work
            elif( (filtLoc > -1) & (requestLoc == -1) ):
                plotNames_now[j] = plotNames_now[j].split(' &F ')[0]; #only filt involved, split at it and keep the bit before to get just the data type
            elif( (filtLoc == -1) & (requestLoc > -1) ):
                plotNames_now[j] = plotNames_now[j].split(' &R ')[0]; #only request involved, split at it and keep the bit before to get just the data type
            else:
                if( filtLoc > requestLoc ):
                    plotNames_now[j] = plotNames_now[j].split(' &R ')[0]; #request 1st so split at it to get just the data type
                else:
                    plotNames_now[j] = plotNames_now[j].split(' &F ')[0]; #filt 1st so split at it to get just the data type
                #END IF
            #END IF

            if( strfind(plotNames_now[j],'-',opt=1) > 0 ): #check for filt info
                #only for use with Mag data types that have x/y/z components
                componentNum = plotNames_now[j][plotNames_now[j].find('-')+1:]; #get the filt requirements
                plotNames_now[j] = plotNames_now[j][:plotNames_now[j].find('-')]; #remove the filt info
                if( componentNum.lower() == 'x' ):
                    componentNum = 0; #component index
                    componentName = '$\mathregular{'+plotNames_now[j]+'_X}$';
                    FLG_components = True;
                elif( componentNum.lower() == 'y' ):
                    componentNum = 1; #component index
                    componentName = '$\mathregular{'+plotNames_now[j]+'_Y}$';
                    FLG_components = True;
                elif( componentNum.lower() == 'z' ):
                    componentNum = 2; #component index
                    componentName = '$\mathregular{'+plotNames_now[j]+'_Z}$';
                    FLG_components = True;
                elif( componentNum.lower() == 'f' ):
                    componentName = '$\mathregular{'+plotNames_now[j]+'_F}$';
                    FLG_components = False;           
                #END IF
            else:
                componentNum = 2; #component index
                componentName = '$\mathregular{'+plotNames_now[j]+'_Z}$';
                FLG_components = True;
            #END IF
            
            #Kp plotting
            if( plotNames_now[j] == 'Kp' ):
                #-----PREP Kp STUFF-----
                if( filtMethod != '' ):
                    pltr_data = subfun_filter( np.copy(data['Kp']['Kp']), filtMethod, dataTime = data['Kp']['time'], dataRate = None, settings_spectra = settings['spectra'], reduceWindow = 0, FLG_reportNaNs = False);
                else:
                    pltr_data = data['Kp']['Kp'];
                #END IF
                Kp_data_plot = np.repeat(pltr_data.flatten(),2); #replicate teh Kp for plotting
                Kp_time_plot = (np.repeat(data['Kp']['time'],2) - dates['date range zero hr dayNum'][1]*86400)/3600; #hrs at 0 hr, replicate the hours for plotting purposes
                Kp_time_plot = np.hstack( (np.nanmin(Kp_time_plot)-3,Kp_time_plot[0:len(Kp_time_plot)-1]) ); #readjust so the hour range matches what is real (e.g. Kp lasts 3 hr, so 0-3 hr is same Kp value)
                plotNames_imHandles[j], = ax[i].plot( Kp_time_plot, Kp_data_plot,linewidth=PLOT_lineWidthDoublePlus, antialiased=True, color=settings_plot['color'][j], linestyle=settings_plot['line style'][j]); #plot
                ylabelString = 'Kp Index';
                ax[i].set_ylim( 0 , np.nanmax(pltr_data)+0.5 ); #set y axis limits
                
                #declare general stuff for afterwards
                pltr_yLims = (0 , np.nanmax(pltr_data)+0.5); #get the ylims
                pltr_time = (data['Kp']['time'] - dates['date range zero hr dayNum'][1]*86400)/3600; #set the time to use for x axis stuff
                plotNames_units[j] = '';
            #OMNI&SuperMAG plotting
            elif( (plotNames_now[j] in OMNI_keyz) | (plotNames_now[j] in SuperMAG_keyz) ):
                #get the data to use
                if( plotNames_now[j] in OMNI_keyz ):
                    pltr_time = (data['OMNI']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                    pltr_data = np.copy(data['OMNI'][plotNames_now[j]]);
                    pltr_name = 'OMNI'; #for reading labels
                else:
                    pltr_time = (data['SuperMAG']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                    pltr_data = np.copy(data['SuperMAG'][plotNames_now[j]]);
                    pltr_name = 'SuperMAG'; #for reading labels
                #END IF
                #apply any filter requested
                if( filtMethod != '' ):
                    pltr_data = subfun_filter( pltr_data, filtMethod, dataTime = pltr_time*3600, dataRate = None, settings_spectra = settings['spectra'], reduceWindow = 0, FLG_reportNaNs = False);
                #END IF
                #make sure no nan gaps
                k = np.isnan(pltr_data) | (pltr_time > FLG_manualLims[1]) | (pltr_time < FLG_manualLims[0]); #get data to yeet
                pltr_data = np.delete(pltr_data,k,axis=0); #delete em
                pltr_time = np.delete(pltr_time,k,axis=0); #delete em
                
                plotNames_imHandles[j], = ax[i].plot( pltr_time, pltr_data , linewidth=PLOT_lineWidthRegular, antialiased=True, color=settings_plot['color'][j], linestyle=settings_plot['line style'][j]); #plot
                
                #--- special coloring for certain data types ---
                if( (plotNames_now[j] == 'IMF clock angle') ):
                    kk = data['OMNI']['Bz GSM'][~k] > 0; #get where Bz GSM is positive (to remove it)
                    kk2 = np.copy(pltr_data); #copy this over
                    kk2[kk] = np.nan;
                    ax[i].plot( pltr_time, kk2 , linewidth=PLOT_lineWidthRegular, color='xkcd:brick red', antialiased=True); #plot
                #END IF
                if( j == 0 ):
                    if( FLG_sameLineUnits ):
                        ylabelString = settings[pltr_name]['labels'][plotNames_now[j]]+settings[pltr_name]['units'][plotNames_now[j]]; #set the y axis label
                        plotNames_units[j] = settings[pltr_name]['units'][plotNames_now[j]];
                    else:
                        ylabelString = settings[pltr_name]['labels'][plotNames_now[j]]+settings[pltr_name]['units'][plotNames_now[j]].replace(' ','\n'); #set the y axis label
                        plotNames_units[j] = settings[pltr_name]['units'][plotNames_now[j]].replace(' ','\n');
                    #END IF
                else:
                    if( j == 1 ):
                        ylabelString = [ylabelString]; #to list for many things
                        ylabelColor = [settings_plot['color'][j-1]]; #record the color
                    #END IF
                    ylabelString.append(' & '); #split
                    ylabelColor.append('xkcd:black');
                    ylabelColor.append(settings_plot['color'][j]);
                    if( FLG_sameLineUnits ):
                        ylabelString.append(settings[pltr_name]['labels'][plotNames_now[j]]+settings[pltr_name]['units'][plotNames_now[j]]); #set the y axis label
                        plotNames_units[j] = settings[pltr_name]['units'][plotNames_now[j]];
                    else:
                        ylabelString.append(settings[pltr_name]['labels'][plotNames_now[j]]+settings[pltr_name]['units'][plotNames_now[j]].replace(' ','\n')); #set the y axis label
                        plotNames_units[j] = settings[pltr_name]['units'][plotNames_now[j]].replace(' ','\n');
                    #END IF
                #END IF
                
                if( j == 0 ):
                    pltr_yLims = (np.nanmin(pltr_data) , np.nanmax(pltr_data)); #get the ylims
                else:
                    ylim_curr = ax[i].get_ylim(); #use prev set values
                    pltr_yLims = (np.min((np.min(ylim_curr),np.nanmin(pltr_data))) , np.max((np.max(ylim_curr),np.nanmax(pltr_data)))); #get the ylims
                #END IF
                ax[i].set_ylim( pltr_yLims ); #set y axis limits
                            
                # if( 98+i < 123): #avoids weird stuff in a weird case
                #     if( (i == 0) | (i == 1) ):
                #         #adjusted y positioning for b. and c.
                #         ax[i].text( letteringPositionX, letteringPositionY+.05, chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i].transAxes); #print the text labelling the lettering a. b. c. ect.
                #         #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
                #     else:
                        # ax[i].text( letteringPositionX, letteringPositionY, chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i].transAxes); #print the text labelling the lettering a. b. c. ect.
                #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
                #     #END IF
                # #END IF
                
                #declare general stuff for afterwards
            elif( plotNames_now[j] in MagCAN_keyz ):
                #get the data to use
                # if( plotNames_now[j] in MagCAN_keyz ):
                pltr_time = (data['MagCAN']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                if( FLG_components == True ):
                    pltr_data = np.copy(data['MagCAN'][plotNames_now[j]]['mag'][:,componentNum]);
                else:
                    pltr_data = np.copy(data['MagCAN'][plotNames_now[j]]['magF']);
                #END IF
                pltr_name = 'MagCAN'; #for reading labels
                # #END IF
                #apply any filter
                if( filtMethod != '' ):
                    pltr_data = subfun_filter( pltr_data, filtMethod, dataTime = pltr_time*3600, dataRate = None, settings_spectra = settings['spectra'], reduceWindow = 0, FLG_reportNaNs = True);
                #END IF
                #make sure no nan gaps
                k = np.isnan(pltr_data) | (pltr_time > FLG_manualLims[1]) | (pltr_time < FLG_manualLims[0]); #get data to yeet
                pltr_data = np.delete(pltr_data,k,axis=0); #delete em
                pltr_time = np.delete(pltr_time,k,axis=0); #delete em
                
                plotNames_imHandles[j], = ax[i].plot( pltr_time, pltr_data , linewidth=PLOT_lineWidthRegular, antialiased=True, color=settings_plot['color'][j], linestyle=settings_plot['line style'][j]); #plot
                
                #--- special coloring for certain data types ---
                if( (plotNames_now[j] == 'IMF clock angle') ):
                    kk = data['OMNI']['Bz GSM'][~k] > 0; #get where Bz GSM is positive (to remove it)
                    kk2 = np.copy(pltr_data); #copy this over
                    kk2[kk] = np.nan;
                    ax[i].plot( pltr_time, kk2 , linewidth=PLOT_lineWidthRegular, color='xkcd:brick red', antialiased=True); #plot
                #END IF
                if( j == 0 ):
                    if( FLG_sameLineUnits ):
                        ylabelString = componentName+' [nT]'; #set the y axis label
                        plotNames_units[j] = ' [nT]';
                    else:
                        ylabelString = componentName+'\n[nT]'; #set the y axis label
                        plotNames_units[j] = '\n[nT]';
                    #END IF
                else:
                    if( j == 1 ):
                        ylabelString = [ylabelString]; #to list for many things
                        ylabelColor = [settings_plot['color'][j-1]]; #record the color
                    #END IF
                    ylabelString.append(' & '); #split
                    ylabelColor.append('xkcd:black');
                    ylabelColor.append(settings_plot['color'][j]);
                    if( FLG_sameLineUnits ):
                        ylabelString.append(componentName+' [nT]'); #set the y axis label
                        plotNames_units[j] = ' [nT]';
                    else:
                        ylabelString.append(componentName+'\n[nT]'); #set the y axis label
                        plotNames_units[j] = '\n[nT]';
                    #END IF
                #END IF
                
                if( j == 0 ):
                    pltr_yLims = (np.nanmin(pltr_data) , np.nanmax(pltr_data)); #get the ylims
                else:
                    ylim_curr = ax[i].get_ylim(); #use prev set values
                    pltr_yLims = (np.min((np.min(ylim_curr),np.nanmin(pltr_data))) , np.max((np.max(ylim_curr),np.nanmax(pltr_data)))); #get the ylims
                #END IF
                ax[i].set_ylim( pltr_yLims ); #set y axis limits
                            
                # if( 98+i < 123): #avoids weird stuff in a weird case
                #     if( (i == 0) | (i == 1) ):
                #         #adjusted y positioning for b. and c.
                #         ax[i].text( letteringPositionX, letteringPositionY+.05, chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i].transAxes); #print the text labelling the lettering a. b. c. ect.
                #         #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
                #     else:
                        # ax[i].text( letteringPositionX, letteringPositionY, chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i].transAxes); #print the text labelling the lettering a. b. c. ect.
                #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
                #     #END IF
                # #END IF
            elif( np.all(FLG_use_AMPERE != None) ):
                #get the data to use
                # if( plotNames_now[j] in MagCAN_keyz ):
                pltr_time = (data['AMPERE']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                pltr_data = FLG_use_AMPERE; #alias
                pltr_name = 'AMPERE'; #for reading labels
                # #END IF
                #apply any filter
                if( filtMethod != '' ):
                    pltr_data = subfun_filter( pltr_data, filtMethod, dataTime = pltr_time*3600, dataRate = None, settings_spectra = settings['spectra'], reduceWindow = 0, FLG_reportNaNs = True);
                #END IF
                #make sure no nan gaps
                k = np.isnan(pltr_data) | (pltr_time > FLG_manualLims[1]) | (pltr_time < FLG_manualLims[0]); #get data to yeet
                pltr_data = np.delete(pltr_data,k,axis=0); #delete em
                pltr_time = np.delete(pltr_time,k,axis=0); #delete em
                
                plotNames_imHandles[j], = ax[i].plot( pltr_time, pltr_data , linewidth=PLOT_lineWidthRegular, antialiased=True, color=settings_plot['color'][j], linestyle=settings_plot['line style'][j]); #plot
                
                #--- special coloring for certain data types ---
                if( j == 0 ):
                    if( FLG_sameLineUnits ):
                        ylabelString = plotNames_now[j]+settings['AMPERE']['units'][plotNames_now[j]]; #set the y axis label
                        plotNames_units[j] = settings['AMPERE']['units'][plotNames_now[j]];
                    else:
                        ylabelString = plotNames_now[j]+'\n'+settings['AMPERE']['units'][plotNames_now[j]].replace(' ',''); #set the y axis label
                        plotNames_units[j] = '\n'+settings['AMPERE']['units'][plotNames_now[j]].replace(' ','');
                    #END IF
                else:
                    if( j == 1 ):
                        ylabelString = [ylabelString]; #to list for many things
                        ylabelColor = [settings_plot['color'][j-1]]; #record the color
                    #END IF
                    ylabelString.append(' & '); #split
                    ylabelColor.append('xkcd:black');
                    ylabelColor.append(settings_plot['color'][j]);
                    if( FLG_sameLineUnits ):
                        ylabelString.append(componentName+' [nT]'); #set the y axis label
                        plotNames_units[j] = ' [nT]';
                    else:
                        ylabelString.append(componentName+'\n[nT]'); #set the y axis label
                        plotNames_units[j] = '\n[nT]';
                    #END IF
                #END IF
                
                if( j == 0 ):
                    pltr_yLims = (np.nanmin(pltr_data) , np.nanmax(pltr_data)); #get the ylims
                else:
                    ylim_curr = ax[i].get_ylim(); #use prev set values
                    pltr_yLims = (np.min((np.min(ylim_curr),np.nanmin(pltr_data))) , np.max((np.max(ylim_curr),np.nanmax(pltr_data)))); #get the ylims
                #END IF
                ax[i].set_ylim( pltr_yLims ); #set y axis limits
            
                #declare general stuff for afterwards
            else:
                try:
                    if( 'detrend' in requestMethod.lower() ):
                        FLG_detrended = True;
                    else:
                        FLG_detrended = False;
                    #END IF
                    if( 'overwrite' in requestMethod.lower().replace(' ','').replace('-','') ):
                        FLG_overwrite = True;
                    else:
                        FLG_overwrite = False;
                    #END IF
                    data['SuperMAG stations'] = GRITI_import_SuperMAG_stations(dates['date range full dayNum'], plotNames_now[j], settings['paths'], settings['config'], coordType='mag', FLG_overwrite=FLG_overwrite, FLG_detrended=FLG_detrended); #import SuperMAG station data from internet [per station, see https://supermag.jhuapl.edu/mag/ for station names]
                except AttributeError:
                    print('ERROR in GRITI_MECHACOMBO_plot: Data type "'+plotNames_now[j]+'" is not in Kp/OMNI/SuperMAG/MagCAN/SuperMAGstations. Fix that. Crashing.');
                    import sys
                    sys.crash()
                #END TRY
                #get the data to use
                pltr_time = (data['SuperMAG stations']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                if( FLG_components == True ):
                   pltr_data = np.copy(data['SuperMAG stations'][plotNames_now[j]]['mag'][:,componentNum]);
                else:
                    pltr_data = np.copy(data['SuperMAG stations'][plotNames_now[j]]['magF']);
                #END IF
                pltr_name = 'SuperMAG stations'; #for reading labels
                # #END IF
                #apply any filter
                if( filtMethod != '' ):
                    pltr_data = subfun_filter( pltr_data, filtMethod, dataTime = pltr_time*3600, dataRate = None, settings_spectra = settings['spectra'], reduceWindow = 0, FLG_reportNaNs = True);
                #END IF
                #make sure no nan gaps
                k = np.isnan(pltr_data) | (pltr_time > FLG_manualLims[1]) | (pltr_time < FLG_manualLims[0]); #get data to yeet
                pltr_data = np.delete(pltr_data,k,axis=0); #delete em
                pltr_time = np.delete(pltr_time,k,axis=0); #delete em
                
                plotNames_imHandles[j], = ax[i].plot( pltr_time, pltr_data , linewidth=PLOT_lineWidthRegular, antialiased=True, color=settings_plot['color'][j], linestyle=settings_plot['line style'][j]); #plot
                
                #--- special coloring for certain data types ---
                if( (plotNames_now[j] == 'IMF clock angle') ):
                    kk = data['OMNI']['Bz GSM'][~k] > 0; #get where Bz GSM is positive (to remove it)
                    kk2 = np.copy(pltr_data); #copy this over
                    kk2[kk] = np.nan;
                    ax[i].plot( pltr_time, kk2 , linewidth=PLOT_lineWidthRegular, color='xkcd:brick red', antialiased=True); #plot
                #END IF
                if( j == 0 ):
                    if( FLG_sameLineUnits ):
                        ylabelString = componentName+' [nT]'; #set the y axis label
                        plotNames_units[j] = ' [nT]';
                    else:
                        ylabelString = componentName+'\n[nT]'; #set the y axis label
                        plotNames_units[j] = '\n[nT]';
                    #END IF
                else:
                    if( j == 1 ):
                        ylabelString = [ylabelString]; #to list for many things
                        ylabelColor = [settings_plot['color'][j-1]]; #record the color
                    #END IF
                    ylabelString.append(' & '); #split
                    ylabelColor.append('xkcd:black');
                    ylabelColor.append(settings_plot['color'][j]);
                    if( FLG_sameLineUnits ):
                        ylabelString.append(componentName+' [nT]'); #set the y axis label
                        plotNames_units[j] = ' [nT]';
                    else:
                        ylabelString.append(componentName+'\n[nT]'); #set the y axis label
                        plotNames_units[j] = '\n[nT]';
                    #END IF
                #END IF
                
                if( j == 0 ):
                    pltr_yLims = (np.nanmin(pltr_data) , np.nanmax(pltr_data)); #get the ylims
                else:
                    ylim_curr = ax[i].get_ylim(); #use prev set values
                    pltr_yLims = (np.min((np.min(ylim_curr),np.nanmin(pltr_data))) , np.max((np.max(ylim_curr),np.nanmax(pltr_data)))); #get the ylims
                #END IF
                ax[i].set_ylim( pltr_yLims ); #set y axis limits
                            
                # if( 98+i < 123): #avoids weird stuff in a weird case
                #     if( (i == 0) | (i == 1) ):
                #         #adjusted y positioning for b. and c.
                #         ax[i].text( letteringPositionX, letteringPositionY+.05, chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i].transAxes); #print the text labelling the lettering a. b. c. ect.
                #         #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
                #     else:
                        # ax[i].text( letteringPositionX, letteringPositionY, chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, transform=ax[i].transAxes); #print the text labelling the lettering a. b. c. ect.
                #uses ASCII coded letters to be able to letter b. -> c. -> d. in a loop
                #     #END IF
                # #END IF
                
                #declare general stuff for afterwards
                #END IF
        #END FOR j
        
        #--- make sure units are the same ---
        if( (FLG_ignoreUnitMismatch == False) & (not all(plotNames_units[0] == j for j in plotNames_units)) ):
            print('ERROR in GRITI_COMBO_KpOMNISuperMAGMag_plot: Values plotted on the same plot "'+plotNames_now+'" do not have the same units ('+plotNames_units+'). Fix that or set FLG_ignoreUnitMismatch to True. Crashing.');
            import sys
            sys.crash()
        #END IF
        
        #--- prep ylabel ---
        if( len(plotNames_now) > 1 ):
            if( all(plotNames_units[0] == j for j in plotNames_units) ):
                #check for redundant units and remove them
                for k in range(0,len(ylabelString)):
                    kk = strstr(ylabelString[k],plotNames_units[0]); #find redundant unit spots
                    if( kk.size != 0 ):
                        ylabelString[k] = ylabelString[k][:kk.item()] + ylabelString[k][kk.item()+len(plotNames_units[0]):]; #remove redundant units
                    #END IF
                #END FOR k
                if( strstr(plotNames_units[0],'\n').size > 0  ):
                    ylabelString.append('\n'); #new line for the new units
                    ylabelColor.append('xkcd:black');
                #END IF
                if( len(plotNames_units) <= 3 ):
                    ylabelString.append(plotNames_units[0].replace('\n','')); #append
                    ylabelColor.append('xkcd:black');
                else:
                    try:
                        ylabelColor.insert(np.where(strfind(ylabelString,plotNames_now[len(plotNames_units)//2]))[0][0]+1,'xkcd:black'); #insert in the middle
                        ylabelString.insert(np.where(strfind(ylabelString,plotNames_now[len(plotNames_units)//2]))[0][0]+1,plotNames_units[0].replace('\n','')); #insert in the middle
                    except:
                        ylabelString.append(plotNames_units[0].replace('\n','')); #append anyway
                        ylabelColor.append('xkcd:black');
                    #END TRY
                #END IF
                
            #END IF
            len(plotNames_units)//2
            ylabelString_len = len(''.join(ylabelString));
            for k in range(0,len(ylabelString)-1):
                if( (ylabelString_len > 10) & (ylabelString[k].find(' & ') == 0) ):
                    ylabelString[k] = '\n'; #replace concatenator with a new line
                #END IF
            #END FOR k
        #     if( FLG_sameLineUnits == False ):
        #         #check for similarly named things and compact them
                
        #         for k in range(0,len(ylabelString)):
        #             for kk in range(0,len(ylabelString[0])):
        #                 ylabelString_common = '';
        #                 if( ylabelString[0][kk] == ylabelString[k][kk] ):
        #                     ylabelString_common += ylabelString[0][kk];
        #                 #END IF
        #             #END FOR kk
        #         #END FOR k
                
        #         ylabelString_split = ylabelString.split(' & '); #split
        #             too hard I gave up
        #         ylabelString_common = ''.join(ylabelString_split[0][i] if ylabelString_split[0][i] == ylabelString_split[k][i] else '*^*' for i in range(0,len(ylabelString_split[0]))); #get the commonalities
        #         ylabelString_common = ylabelString_common[:ylabelString_common.find('*^*')];
        #         for j in range(0,len(plotNames_now)):
        #             plotNames_now[j]
        #         #END FOR j
        #     #END IF
        # #END IF
        
        #--- make legend if needed ---
        if( len(plotNames_now) > 1 ):
            # ax[i].legend(plotNames_imHandles, [x for x in ylabelString if((x != ' & ') & (x != '\n')) ],ncol=len(plotNames_now),loc='lower right',framealpha=0.5); #add the legend
            secretLabel.append(ax[i].set_ylabel(''.join(ylabelString),color='xkcd:red',alpha=0.5,fontproperties=FONT_axisLabelFM)); #set the y axis label
            secretLabel_string.append(ylabelString);
            secretLabel_color.append(ylabelColor);
        else:
            #otherwise ylabel it
            ax[i].set_ylabel(ylabelString,fontproperties=FONT_axisLabelFM); #set the y axis label
        #END IF
        
        #--- apply common axis stuff ---
        ax[i].grid(visible=True, which='major', axis='both', color='xkcd:grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
        
        #--- draw on big diffs to the time ref ---
        if( ( ((np.nanmin(dates['date range zero hr hour bounds'])*3600 + dates['date range zero hr dayNum'][1]*86400) - np.nanmin(data['time ref'])) <= -0.25*86400) & ((settings_plot['time ref'] != 'Kp') & (settings_plot['time ref'] != 'OMNI') & (settings_plot['time ref'] != 'SuperMAG')) ): #as long as min Kp time is 15 min diff or more from the other time reference, plot where the time ref begins (not Kp tho)
            ax[i].plot( np.repeat( (np.nanmin(data['time ref']) - dates['date range zero hr dayNum'][1]*86400)/3600 , 2) , np.linspace(np.nanmin(pltr_yLims),np.nanmax(pltr_yLims),num=2), linewidth=PLOT_lineWidthPlus, color='r', antialiased=True); #plot red lines showing ISR data time
        #END IF
        if( ( ((np.nanmin(dates['date range zero hr hour bounds'])*3600 + dates['date range zero hr dayNum'][1]*86400) - np.nanmin(data['time ref'])) >= 0.25*86400) & ((settings_plot['time ref'] != 'Kp') & (settings_plot['time ref'] != 'OMNI') & (settings_plot['time ref'] != 'SuperMAG')) ): #as long as max Kp time is 15 min diff or more from the other time reference, plot where the time ref ends (not Kp tho)
            ax[i].plot( np.repeat( (np.nanmax(data['time ref']) - dates['date range zero hr dayNum'][1]*86400)/3600 , 2) , np.linspace(np.nanmin(pltr_yLims),np.nanmax(pltr_yLims),num=2), linewidth=PLOT_lineWidthPlus, color='r', antialiased=True); #plot red lines showing ISR data time
        #END IF
        
        #--- draw optional lines and dates on ---
        if( (opt >= 1) & (np.abs(np.diff(FLG_manualLims)) >= 18) ):
            for j in range(0,dates['date range full'].shape[0]):
                #run through each day, print the line
                if( j != 0 ):
                    ax[i].plot( np.repeat( dates['date range zero hr hours'][j]-dayNite_DSToffset, 2) , np.linspace(np.nanmin(pltr_yLims),np.nanmax(pltr_yLims),num=2), linewidth=PLOT_lineWidthRegularPlus, color='k', linestyle='--', antialiased=True); #plot black dotted lines showing the date separation
                #END IF
                if( (i in plt_topIndexr) & (opt >= 2) ): #only print on the top plot
                    # topOfPlot = ax[i].transLimits.transform((0,1));
                    if( opt == 2 ):
                        ax[i].text(dates['date range zero hr hours'][j]+12-dayNite_DSToffset, np.nanmax(pltr_yLims), str(dates['date range full'][j,1])+'/'+str(dates['date range full'][j,2])+'/'+str(dates['date range full'][j,0]), fontproperties=FONT_axisLabelFM, horizontalalignment='center', verticalalignment='bottom'); #print the text saying the days
                    elif( opt == 3 ):
                        ax[i].text(dates['date range zero hr hours'][j]+12-dayNite_DSToffset, np.nanmax(pltr_yLims), str(dates['date range full dayNum'][j,1])+', '+str(dates['date range full dayNum'][j,0]), fontproperties=FONT_axisLabelFM, horizontalalignment='center', verticalalignment='bottom' ); #print the text saying the days
                    #END IF
                    # if( opt == 2 ):
                    #     ax[i].text(dates['date range zero hr hours'][j]+12, np.nanmax(pltr_yLims)-(np.nanmax(pltr_yLims)-np.nanmin(pltr_yLims))*.15, str(dates['date range full'][j,1])+'/'+str(dates['date range full'][j,2])+'/'+str(dates['date range full'][j,0]), fontproperties=FONT_axisLabelFM, horizontalalignment='center', verticalalignment='bottom'); #print the text saying the days
                    # elif( opt == 3 ):
                    #     ax[i].text(dates['date range zero hr hours'][j]+12, np.nanmax(pltr_yLims)-(np.nanmax(pltr_yLims)-np.nanmin(pltr_yLims))*.15, str(dates['date range full dayNum'][j,1])+', '+str(dates['date range full dayNum'][j,0]), fontproperties=FONT_axisLabelFM, horizontalalignment='center', verticalalignment='bottom'); #print the text saying the days
                    # #END IF
                #END IF
            #END FOR j
            if( dayNite_DSToffset == 0 ):
                pass; #speeds it up a tiny bit
            elif( dayNite_DSToffset < 0 ):
                ax[i].plot( np.repeat( dates['date range zero hr hours'][0]-dayNite_DSToffset, 2) , np.linspace(np.nanmin(pltr_yLims),np.nanmax(pltr_yLims),num=2), linewidth=PLOT_lineWidthRegularPlus, color='k', linestyle='--', antialiased=True); #plot black dotted lines showing the date separation
            else:
                ax[i].plot( np.repeat( dates['date range zero hr hours'][-1]-dayNite_DSToffset, 2) , np.linspace(np.nanmin(pltr_yLims),np.nanmax(pltr_yLims),num=2), linewidth=PLOT_lineWidthRegularPlus, color='k', linestyle='--', antialiased=True); #plot black dotted lines showing the date separation
            #END IF
        #END IF
                
        #--- deal with x axis stuff ---
        if( FLG_fancyPlot == 0 ):
            if( i in plt_bottomIndexr ):
                GRITI_plotHelper_axisizerTime(pltr_time,ax=ax[i],tickNumGoal=28//plt_ncols,tickReducer=2,FLG_manualLims=FLG_manualLims,FLG_removeLabels=False,FLG_tickDirIn=True); #automagic time ticks here
            else:
                GRITI_plotHelper_axisizerTime(pltr_time,ax=ax[i],tickNumGoal=28//plt_ncols,tickReducer=2,FLG_manualLims=FLG_manualLims,FLG_removeLabels=True,FLG_tickDirIn=True); #automagic time ticks here
            #END IF
        else:
            if( i in plt_bottomIndexr ):
                GRITI_plotHelper_axisizerTime(pltr_time,ax=ax[i],tickNumGoal=21//plt_ncols,tickReducer=2,FLG_manualLims=FLG_manualLims,FLG_removeLabels=False,FLG_tickDirIn=True); #automagic time ticks here
            else:
                GRITI_plotHelper_axisizerTime(pltr_time,ax=ax[i],tickNumGoal=21//plt_ncols,tickReducer=2,FLG_manualLims=FLG_manualLims,FLG_removeLabels=True,FLG_tickDirIn=True); #automagic time ticks here
            #END IF
        #END IF
    #END FOR i
    
    #label the bottom
    if( FLG_MECHACOMBO_plot_localTime == False ):
        if( plt_ncols == 2 ):
            fig.supxlabel('Time in UT [hr] - 0 Hr on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr day post fix']+' | Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        else:
            ax[i].set_xlabel('Time in UT [hr] - 0 Hr on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr day post fix']+' | Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        #END IF
    else:
        #make it local time, some legwork required        
        #mod the labels now
        # fig.draw();
        fig.canvas.draw_idle(); #key or everything fails
        # plt.draw()
        fig.canvas.flush_events(); 
        xTickLabels = copy.deepcopy(ax[i].get_xticklabels()); #get the labels
        for j in range(0,len(xTickLabels)):
            xTickLabels[j].set_text(textNice(xTickLabels[j].get_position()[0]+dayNite_DSToffset)); #adjust
        #END FOR j
        ax[i].set_xticklabels(xTickLabels); #set the labels
        plt.close(); #it opens a 2nd window because idk, so just close it???
        
        if( plt_ncols == 2 ):
            fig.supxlabel('Local Time [hr] | '+dayNite_timeZoneName+' = UT '+dayNite_DSToffset_str+' hr | 0 hr on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr day post fix']+' (Day '+str(dates['date range zero hr dayNum'][1])+') '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        else:
            ax[i].set_xlabel('Local Time [hr] | '+dayNite_timeZoneName+' = UT '+dayNite_DSToffset_str+' hr | 0 hr on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr day post fix']+' (Day '+str(dates['date range zero hr dayNum'][1])+') '+str(dates['date range zero hr dayNum'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        #END IF
    #END IF
    
    #highlight based on Bz being orientated south (negative)
    if( any('Bz GSM' in strang for strang in FLG_MECHACOMBO_names) or any('Bz GSE' in strang for strang in FLG_MECHACOMBO_names) ):
        if( any('Bz GSM' in strang for strang in FLG_MECHACOMBO_names) ):
            pltr_time = (data['OMNI']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
            pltr_data = np.copy(data['OMNI']['Bz GSM']);
            axIndx = ['Bz GSM' in strang for strang in FLG_MECHACOMBO_names].index(True); #get the axis index
        else:
            pltr_time = (data['OMNI']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
            pltr_data = np.copy(data['OMNI']['Bz GSE']);
            axIndx = ['Bz GSE' in strang for strang in FLG_MECHACOMBO_names].index(True); #get the axis index
        #END IF
        k = np.isnan(pltr_data) | (pltr_time > FLG_manualLims[1]) | (pltr_time < FLG_manualLims[0]); #get data to yeet
        pltr_data = np.delete(pltr_data,k,axis=0); #delete em
        pltr_time = np.delete(pltr_time,k,axis=0); #delete em
        
        k = pltr_data >= 0; #remove positive stuff
        pltr_data[k] = np.nan; #nan it for the plot
        pltr_time[k] = np.nan; #nan it for the plot
        ax[axIndx].plot( pltr_time, pltr_data , linewidth=PLOT_lineWidthRegular, antialiased=True, color='xkcd:brick red', linestyle=settings_plot['line style'][0]); #plot
        pltr_data = np.delete(pltr_data,k,axis=0); #delete em
        pltr_time = np.delete(pltr_time,k,axis=0); #delete em

        
        k = np.isclose(np.diff(pltr_time), data['OMNI']['data rate']/3600); #find contiguous time bits
        kj = np.where(~k)[0]+1; #add 1 for index b/c diff
        if(kj[0] != 0):
            kj = np.insert(kj, 0, 0); #insert 0
        #END IF
        if(kj[-1] != (pltr_time.size-1) ):
            kj = np.append(kj, pltr_time.size-1); #insert end size
        #END IF
        for i in range(0,len(FLG_MECHACOMBO_names)):
            for j in range(0,kj.size-1):
                ax[i].axvspan(pltr_time[kj[j]], pltr_time[kj[j+1]-1], ymin=0, ymax=1, alpha=0.45, edgecolor='none', facecolor='xkcd:brick red');
            #END FOR j
        #END FOR i      
    #END IF

    figFitter(fig); #fit the fig fast
    #remakes labels with color if needed - must be called after figFitter othwrise the plot isn't "settled" enough to work
    
    #make sure at least 3 values are shown on the graphs after the fig fitting
    for i in range(0,len(FLG_MECHACOMBO_names)):
        limz = ax[i].get_ylim();
        tickz = ax[i].get_yticks();
        tickz_within_limz = (tickz <= np.max(limz)) & (tickz >= np.min(limz)); #get which tickz are within the limz
        if( np.sum(tickz_within_limz) < 3 ):
            tickz_within_limz_where = np.where(tickz_within_limz)[0]; #get where the tickz are within the limz
            tickz_within_limz_whereNot = np.where(~tickz_within_limz)[0]; #get where the tickz are within the limz
            if( np.sum(tickz_within_limz) == 1 ):
                limz = (tickz[np.min(tickz_within_limz_where)-1], tickz[np.max(tickz_within_limz_where)+1]); #expand both limz ends
            else:
                maxCheck = np.abs(tickz[~tickz_within_limz] - np.max(limz)); #use where +/- changes to check for if bounds are reasonable
                minCheck = np.abs(tickz[~tickz_within_limz] - np.min(limz)); #use where +/- changes to check for if bounds are reasonable
                if( np.min(minCheck) < np.min(maxCheck) ): #choose whichever is closer to 0 (thus closer to adding a new pt)
                    limz = (tickz[np.min(tickz_within_limz_where)-1], np.max(limz)); #set the limz to be the next tick mark value
                else:
                    limz = (np.min(limz), tickz[np.max(tickz_within_limz_where)+1]); #set the limz to be the next tick mark value
                #END IF
            #END IF
            ax[i].set_ylim(limz); #set the new ylimz
        #END IF
    #END FOR i
    
    #prep secretLabel_string for alignment
    # for i in range(0,len(secretLabel_string)):
    #     jk = np.where(strfind(secretLabel_string[i],' '))[0]; #get strings with spaces in them
    #     for j in range(0,jk.size):
    #         if( secretLabel_string[i][jk[j]][0] == ' ' ): #only do this on stuff supposed to be spaced together, not words with spaces in them
    #             tempLen = len(secretLabel_string[i][jk[j]]); #hold this
    #             if( '$\\mathregular{' in secretLabel_string[i][jk[j]-1] ):
    #                 secretLabel_string[i][jk[j]] = ''.ljust(len(secretLabel_string[i][jk[j]-1].replace('$\\mathregular{','').replace('_','').replace('}$','')))+secretLabel_string[i][jk[j]]; #also offset spaced string with spaces so everything lines up
    #             else:
    #                 secretLabel_string[i][jk[j]] = ''.ljust(len(secretLabel_string[i][jk[j]-1]))+secretLabel_string[i][jk[j]]; #also offset spaced string with spaces so everything lines up
    #             #END IF
    #             secretLabel_string[i][jk[j]-1] = secretLabel_string[i][jk[j]-1]+''.ljust(tempLen+1); #offset 1st string with spaces
    #         #END IF
    #     #END FOR j
    # #END FOR i
    cntr = 0;
    FLG_spacer = False; #keep false usually
    for i in range(0,len(FLG_MECHACOMBO_names)):
        #this and the extra loop allow for stuff to be plotted on the same plot
        if( isinstance(FLG_MECHACOMBO_names[i],list) | isinstance(FLG_MECHACOMBO_names[i],tuple) ):
            # multicolor_ylabel(ax[i],ylabelString,ylabelColor,axis='y',anchorpad=0,fontproperties=FONT_axisLabelFM,position=secretLabel.get_position())
            # axisLabel_magical(fig, ax[i], ylabelString, ylabelColor, (0,0), axis='y', fontproperties=FONT_axisLabelFM);
            plt.draw();
            secretLabel[cntr].draw(fig.canvas.get_renderer());
            txtPos = secretLabel[cntr].get_window_extent().transformed(ax[i].transAxes.inverted()); #get bbox for the secret label
            txtPos = [txtPos.x0, (txtPos.y0+txtPos.y1)/2]; #place it where the ylabel actually is
            txtPosOrig = copy.deepcopy(txtPos); #keep an original copy
            for jk in range(0,len(secretLabel_color[cntr])):
                if( secretLabel_string[cntr][jk] != '\n' ):
                    if( secretLabel_string[cntr][jk][0] != ' ' ):
                        #if text write the text and adjust for new text to the "right" of it
                        text = ax[i].text(txtPos[0],txtPos[1],secretLabel_string[cntr][jk],color=secretLabel_color[cntr][jk], transform=ax[i].transAxes, \
                            rotation=90, va='center', ha='left', fontproperties=FONT_axisLabelFM);
                        text.draw(fig.canvas.get_renderer());
                        txtPos[1] += text.get_window_extent().transformed(ax[i].transAxes.inverted()).height;
                    else:
                        #two spaced words need lots of black magicks
                        text.set_alpha(0); #hide the 1st word that was already made
                        txtPos[1] = txtPosOrig[1]; #reset
                        textSpacer = ax[i].text(txtPos[0],txtPos[1],secretLabel_string[cntr][jk-1]+secretLabel_string[cntr][jk],color='xkcd:bright red', transform=ax[i].transAxes, \
                            rotation=90, va='center', ha='left', fontproperties=FONT_axisLabelFM);
                        textSpacer.draw(fig.canvas.get_renderer());
                        textSpacerPos = textSpacer.get_window_extent().transformed(ax[i].transAxes.inverted()); #get bbox for the secret label
                        txtPos[1] = textSpacerPos.y0; #align to the bottom
                        #1st word
                        text = ax[i].text(txtPos[0],txtPos[1],secretLabel_string[cntr][jk-1],color=secretLabel_color[cntr][jk-1], transform=ax[i].transAxes, \
                            rotation=90, va='bottom', ha='left', fontproperties=FONT_axisLabelFM); #note va='bottom'
                        text.draw(fig.canvas.get_renderer());
                        txtPos[1] += text.get_window_extent().transformed(ax[i].transAxes.inverted()).height;
                        
                        #2nd word
                        text = ax[i].text(txtPos[0],txtPos[1],secretLabel_string[cntr][jk],color=secretLabel_color[cntr][jk], transform=ax[i].transAxes, \
                            rotation=90, va='bottom', ha='left', fontproperties=FONT_axisLabelFM); #note va='bottom'
                        text.draw(fig.canvas.get_renderer());
                        txtPos[1] += text.get_window_extent().transformed(ax[i].transAxes.inverted()).height;
                        
                        textSpacerPosWidth = copy.deepcopy(textSpacer.get_window_extent().transformed(ax[i].transAxes.inverted()).width);
                        textSpacer.set_alpha(0); #hide the placeholder words
                        FLG_spacer = True;
                    #END IF
                else:
                    #if a new line undo the "right" adjustment and instead move to a new line for the next text
                    # txtPos[1] -= text.get_window_extent().transformed(ax[i].transAxes.inverted()).height; #undo
                    txtPos[1] = txtPosOrig[1]; #reset
                    if( FLG_spacer == False ):
                        txtPos[0] += text.get_window_extent().transformed(ax[i].transAxes.inverted()).width; #move right
                    else:
                        txtPos[0] += textSpacerPosWidth; #move right wholistically [purposely spellt]
                    #END IF
                    FLG_spacer = False; #keep false usually
                #END IF
            #END FOR jk
            secretLabel[cntr].set_alpha(0); #required after, if starts invisible it's ignored
            cntr += 1; #increment
        #END IF
        # letteringPositionTrans = ax[i].transLimits.inverted().transform((letteringPositionX,letteringPositionY));
        ax[i].text( np.min(ax[i].get_xlim())+np.diff(ax[i].get_xlim())*.005, np.max(ax[i].get_ylim()), chr(97+i)+'.', color='r', fontproperties=FONT_grandioseFM, horizontalalignment='left', verticalalignment='top'); #print the text labelling the lettering a. b. c. ect.
    #END FOR i
    figFitter(fig); #fit the fig fast AGAIN
    if( FLG_fancyPlot != 0 ):
        fig.savefig(os.path.join(settings_paths['fancyPlots'],'MECHACOMBO_fancyPlot'+settings_plot['save file type'])); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF