"""
GOAL: Plot only ISR POPL HP RTI
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
# from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from os.path import join as jointer
# from Code.subfun_monthNum_to_word import subfun_monthNum_to_word
from Code.subfun_textNice import textNice
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.subfun_figFitter import figFitter

def GRITI_ISR_Haystack_plot_POPL_HP(Zenith_time,Zenith_height,Zenith_POPL_hp,MISA_time,MISA_height,MISA_POPL_hp, \
                                    filter_cutoffPeriod,ISR_RTI_heightLimValues,ISR_POPL_plotLimValu, \
                                    dates, settings_plot, settings_paths, \
                                    time_cutout_range=None, FLG_dayNite=0, settings_map=None, settings_config=None, FLG_fancyPlot=0):
    #----- Declare -----
    if( (FLG_dayNite >= 1) & np.all(settings_map == None) ):
        print('WARNING IN GRITI_ISR_Haystack_plot_POPL_HP: FLG_dayNite ON ('+str(FLG_dayNite)+') but no settings_map provided. It needs it for the geographic area.\nSend it in to fix! Disabling FLG_dayNite for now.');
        FLG_dayNite = 0; #disable
    #END IF
    if( FLG_fancyPlot >= 1 ):
        fileName = 'ISR_RTI_Z&M_POPL+HP'; #prep
        if( np.all(time_cutout_range != None) ):
            fileName = fileName+'_cut'+textNice(np.min(time_cutout_range))+'to'+textNice(np.max(time_cutout_range)); #tack on hrs
        if( FLG_dayNite >= 1 ):
            fileName = fileName+'_wDayNite'; #tack on wDayNite
            from Code.subfun_sunAlsoRises import sunAlsoRises
            from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
            import timezonefinder
            from urllib.request import urlopen
            from Code.subfun_strstr import strstr
            import pytz
            from datetime import datetime
        #END IF
        print('MAKING FANCY PLOT: '+fileName+' IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    ISR_m3tocc = 100**3; #1 m^3 is 100^3 cm^3 and 1 cm^3 is 1 cc
    ISR_POPL_plotLimValu = ISR_POPL_plotLimValu/ISR_m3tocc; #adjust this plot limit value
    Zenith_POPL_hp = Zenith_POPL_hp/ISR_m3tocc; #adjust POPL values
    MISA_POPL_hp = MISA_POPL_hp/ISR_m3tocc; #adjust POPL values
    yAxisTicks = np.arange(100,ISR_RTI_heightLimValues[1]+50,50); #get the y axis tick locations
    
    letteringPositionX = -0.105; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.88; #set the X position of the lettering (e.g., a. b. c. ...)
    
    #----- Plot ISR POPL HP results as a RTI -----
    #unpack
    # dateRange = dates['date range'];
    dateRange_zeroHr = dates['date range zero hr'];
    dateRange_dayNum_zeroHr = dates['date range zero hr dayNum'];
    
    
    #----- Plot just the ISR POPL HP results -----
    if( FLG_dayNite == 0 ):
        figNumRows = 2;
    else:
        figNumRows = 3; #xtra for the dayNite
    #END IF
    if( FLG_fancyPlot == 0 ):
        if( figNumRows == 2 ):
            fig, ax = plt.subplots(nrows=figNumRows, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        else:
            import matplotlib.gridspec as gridspec
            fig = plt.figure();
            # gridr = gridspec.GridSpec(nrows=7, ncols=1, figure=fig);
            # fig.subplots_adjust(hspace=0.75)
            gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
            gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=20, ncols=1, subplot_spec = gridr[0], hspace = 6.0)
            gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
            # gridr.update(hspace=0.05); # set the spacing between axes.
            fig.add_subplot(gridr1[0:10]); #RTI plots are 2 tall
            # gridr.update(hspace=0.80); # set the spacing between axes.
            fig.add_subplot(gridr1[10:20]); #RTI plots are 2 tall
            # gridr.update(hspace=0.90); # set the spacing between axes.
            fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
            ax = fig.axes; #get a list of the axes
        #END IF
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
    else:
        plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        if( figNumRows == 2 ):
            fig, ax = plt.subplots(nrows=figNumRows, ncols=1, figsize=(14,8.5), dpi=settings_plot['journal dpi']); #use instead of fig because it inits an axis too (I think I dunno)
        else:
            import matplotlib.gridspec as gridspec
            fig = plt.figure(figsize=(14,10.5),dpi=settings_plot['journal dpi']);
            # gridr = gridspec.GridSpec(nrows=7, ncols=1, figure=fig);
            # fig.subplots_adjust(hspace=0.75)
            gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
            gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=20, ncols=1, subplot_spec = gridr[0], hspace = 6.0)
            gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
            # gridr.update(hspace=0.05); # set the spacing between axes.
            fig.add_subplot(gridr1[0:10]); #RTI plots are 2 tall
            # gridr.update(hspace=0.80); # set the spacing between axes.
            fig.add_subplot(gridr1[10:20]); #RTI plots are 2 tall
            # gridr.update(hspace=0.90); # set the spacing between axes.
            fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
            ax = fig.axes; #get a list of the axes
        #END IF
    #END IF
    #Remove the aspect ratio from the basemap so it fills the screen better
    for i in range(0,len(ax)):
        ax[i].set_aspect('auto');
    #END FOR i
    
    #----- prep if time_cutout_range is active -----
    if( np.all(time_cutout_range != None) ):
        time_cutout_indexes = np.array( ( np.where(np.min(np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) )) == np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) ) )[0][0] , \
            np.where(np.min(np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) )) == np.abs( (Zenith_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) ) )[0][0] ) ); #get the indexes for that time cutout range
    
        Zenith_time = np.copy(Zenith_time[time_cutout_indexes[0]:time_cutout_indexes[1]+1]); #protect memory w/ copy
        Zenith_POPL_hp = np.copy(Zenith_POPL_hp[:,time_cutout_indexes[0]:time_cutout_indexes[1]+1]);
        
        time_cutout_indexes = np.array( ( np.where(np.min(np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) )) == np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.min(time_cutout_range) ) )[0][0] , \
            np.where(np.min(np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) )) == np.abs( (MISA_time-dateRange_dayNum_zeroHr[1]*86400)/3600 - np.max(time_cutout_range) ) )[0][0] ) ); #get the indexes for that time cutout range
                                        
        MISA_time = np.copy(MISA_time[time_cutout_indexes[0]:time_cutout_indexes[1]+1]);
        MISA_POPL_hp = np.copy(MISA_POPL_hp[:,time_cutout_indexes[0]:time_cutout_indexes[1]+1]);
    #END IF
    
    #----- ZENITH STUFF ----- 
    divider = make_axes_locatable(ax[0]); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    pltHelprX, pltHelprY = np.meshgrid( (Zenith_time - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                Zenith_height);
    im = ax[0].pcolormesh(pltHelprX , pltHelprY , Zenith_POPL_hp , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu) , cmap='gray' , shading='gouraud', antialiased=True); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar.mappable.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
    #cbar.ax.set_yticklabels(np.round(np.linspace(-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu,5), len(str(ISR_POPL_plotLimValu).split('.')[1])+1 )); #create useful tick marks
    # cbar.set_label("POPL [$e^-/m^3$]"); #tabel the colorbar
    cbar.set_label("N$_{e^-}$ [$e^-/cc$]"); #tabel the colorbar
    cbar.ax.tick_params(labelsize=settings_plot['font axis tick']);
    cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    if( FLG_fancyPlot == 0 ):
        string_Title = 'Zenith POPL with High-pass Filter '+textNice(np.round(filter_cutoffPeriod/3600,2))+' hr Cutoff'; #create mecha title
        if( np.all(time_cutout_range != None) ):
            string_Title = string_Title + ' - Hours '+textNice(np.min(time_cutout_range))+' to '+textNice(np.max(time_cutout_range));
        #END IF
        ax[0].set_title(string_Title,fontproperties=settings_plot['font title FM']); #set the title
    #END IF
    # ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=settings_plot['font axis label FM']); #set the x axis label
    ax[0].set_ylabel('Height [km]',fontproperties=settings_plot['font axis label FM']); #set the y axis label
    
    
    xAxis_ticks, xAxis_lims = GRITI_plotHelper_axisizerTime(Zenith_time/3600,ax=ax[0],unit='hr',FLG_manualLims=((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(Zenith_time)-dateRange_dayNum_zeroHr[1]*86400)/3600),FLG_removeLabels=False);
    ax[0].set_yticks(yAxisTicks); #set y axis ticks
    ax[0].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[0].transAxes); #print the text saying the day or nite
    
    #----- MISA STUFF ----- 
    divider = make_axes_locatable(ax[1]); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    pltHelprX, pltHelprY = np.meshgrid( (MISA_time - dateRange_dayNum_zeroHr[1]*86400)/3600, \
                MISA_height);
    im = ax[1].pcolormesh(pltHelprX , pltHelprY , MISA_POPL_hp , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu) , cmap='gray', shading='gouraud', antialiased=True); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar.mappable.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
    #cbar.ax.set_yticklabels(np.round(np.linspace(-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu,5), len(str(ISR_POPL_plotLimValu).split('.')[1])+1 )); #create useful tick marks
    # cbar.set_label("POPL [$e^-/m^3$]"); #tabel the colorbar
    cbar.set_label("N$_{e^-}$ [$e^-/cc$]"); #tabel the colorbar
    #    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
    cbar.ax.tick_params(labelsize=settings_plot['font axis tick']);
    cax.yaxis.label.set_font_properties(settings_plot['font axis label FM']);
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    if( FLG_fancyPlot == 0 ):
        string_Title = 'MISA POPL with High-pass Filter '+textNice(np.round(filter_cutoffPeriod/3600,2))+' hr Cutoff'; #create mecha title
        if( np.all(time_cutout_range != None) ):
            string_Title = string_Title + ' - Hours '+textNice(np.min(time_cutout_range))+' to '+textNice(np.max(time_cutout_range));
        #END IF
        ax[1].set_title(string_Title,fontproperties=settings_plot['font title FM']); #set the title
    #END IF
    ax[1].set_xlabel('Time in UT [hr] - 0 Hr on '+dates['date range zero hr month name']+' '+str(dateRange_zeroHr[2])+dates['date range zero hr day post fix']+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=settings_plot['font axis label FM']); #set the x axis label
    ax[1].set_ylabel('Height [km]',fontproperties=settings_plot['font axis label FM']); #set the y axis label
    
    # GRITI_plotHelper_axisizerTime(MISA_time/3600,ax=ax[1],unit='hr',FLG_manualLims=((np.min(MISA_time)-dateRange_dayNum_zeroHr[1]*86400)/3600 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1]*86400)/3600),FLG_removeLabels=False);
    ax[1].set_xticks(xAxis_ticks); #set x axis ticks
    ax[1].set_xlim( xAxis_lims ); #set x axis limits (based off of Zenith)
    ax[1].set_yticks(yAxisTicks); #set y axis ticks
    ax[1].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[1].transAxes); #print the text saying the day or nite
    
    
    #----- DAYNITE STUFF IF NEEDED ----- 
    if( FLG_fancyPlot >= 1 ):       
        #FIRST: GET SUNRISE/SUNSET TIMES
        (dayNite_sunrise, dayNite_sunset, daynites_dateRange_fullPad) = sunAlsoRises(dates['date range full'],settings_map['site coords'][0][0],settings_map['site coords'][0][1]); #call sunrise/set function
    
        daynite_dateRange_dayNum_fullPad = subfun_date_to_dayNum(daynites_dateRange_fullPad); #convert to dayNum
        dayNite_sunrise = (dayNite_sunrise + daynite_dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hrs, center around zero hr and convert ot hrs
        dayNite_sunset = (dayNite_sunset + daynite_dateRange_dayNum_fullPad[:,1] - dateRange_dayNum_zeroHr[1])*24; #hrs, center around zero hr and convert ot hrs
        dayNite_sunrise = dayNite_sunrise[1:]; #remove 1st
        dayNite_sunset = dayNite_sunset[:-1]; #remove last
        
        #SECOND: GET TIME ZONE FOR LOCAL TIME
        tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
        dayNite_timeZoneID = tf.certain_timezone_at(lat=settings_map['site coords'][0][0], lng=settings_map['site coords'][0][1]); #use it to find the time zone
        if dayNite_timeZoneID is None:
            #use geonames site as a backup
            url = 'http://api.geonames.org/timezone?lat='+str(settings_map['site coords'][0][0])+'&lng='+str(settings_map['site coords'][0][1])+'&username='+settings_config['login GeoNames TimeZone']['user']; #create link for lat/long
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
            # timeZone_offset_str = webpage[index_start+11:index_end]; #get dst offset
            # timeZone_offset = np.float64(timeZone_offset_str); #and convert to number
            # index_start = strstr(webpage,'<gmtOffset>')[0]; #get where UT offset is
            # index_end = strstr(webpage,'</gmtOffset>')[0]; #get where UT offset is
            # dayNite_UToffset_str = webpage[index_start+11:index_end]; #get UT offset
            # dayNite_UToffset = np.float64(dayNite_UToffset_str); #and convert to number
            index_start = strstr(webpage,'<timezoneId>')[0]; #get where time zone ID is
            index_end = strstr(webpage,'</timezoneId>')[0]; #get where time zone ID is
            dayNite_timeZoneID = webpage[index_start+12:index_end]; #get time zone ID
        #END IF
        
        #THIRD: TIME ZONE STUFF
        timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
        timeZoneObj_UTC = pytz.timezone('UTC'); #make a timezone object
        
        timeZone_datetime = timeZoneObj_UTC.localize(datetime.strptime(str(dates['date range zero hr'][0])+'-'+str(dates['date range zero hr'][1])+'-'+str(dates['date range zero hr'][2])+'T00:00:00.000', \
            '%Y-%m-%dT%H:%M:%S.%f')).astimezone(timeZoneObj); #create datetime object, set it to the UTC time zone (which it is, datetime just doesn't know), then convert it to local time zone
        timeZone_offset_str = timeZone_datetime.isoformat()[strstr(timeZone_datetime.isoformat(),'-')[-1]:]; #get the DST time offset
        if( strstr(timeZone_offset_str,':').size > 0 ):
            if( timeZone_offset_str[strstr(timeZone_offset_str,':')[0]+1:] == '00' ):
                timeZone_offset_str = timeZone_offset_str[:strstr(timeZone_offset_str,':')[0]]; #remove the :00 if it's just that
                if( timeZone_offset_str[1] == '0' ):
                    timeZone_offset_str = timeZone_offset_str.replace('0',''); #remove the 0 that was extraneous
                #END IF
                timeZone_offset = np.int64(timeZone_offset_str); #get the number version
            #END IF
            else:
                timeZone_offset = np.float64(timeZone_offset_str[:strstr(timeZone_offset_str,':')[0]]) + \
                    np.int64(timeZone_offset_str[strstr(timeZone_offset_str,':')[0]+1:])/60; #convert to an hour decimal
            #END IF
        #END IF
        timeZone_name = timeZone_datetime.tzname(); #get the time zone name (like 'EST' or 'EDT' depending on standard or daylight savings time)
        timeZone_DSTnUTCOffset = timeZone_datetime.dst().total_seconds()/3600; #get the time offset
        if( np.mod(timeZone_DSTnUTCOffset,1) == 0 ):
            timeZone_DSTnUTCOffset = np.int64(timeZone_DSTnUTCOffset); #convert to integer
        #END IF
            
        #FOURTH STEP: PREP FOR PLOTTING BY ALIGNING TIMES, MAKING PLOT VARIABLES
        dayNite_sunrise += timeZone_offset; #make local time
        dayNite_sunset += timeZone_offset; #make local time
    
        xTime = np.sort( np.concatenate( (dayNite_sunrise,dayNite_sunset) ) ); #hr, xvar to plot against
        yDayNite = np.ones(xTime.shape); #prep if day or night, set all to day
        for i in range(0,dayNite_sunset.size):
            yDayNite[xTime == dayNite_sunset[i]] = 0; #set sunset times to sunset
        #END FOR i
        xTime = xTime.repeat(2); #interleave repeated values
        yDayNite = np.roll(yDayNite.repeat(2),1) #interleave repeated values and circular shift by 1
        yDayNite[0] = yDayNite[1]; #set that to match (for plotting niceness)
        yDayNite[-1] = yDayNite[-2]; #set that to match (for plotting niceness)
        
        # xTime[ xTime < np.min(xAxis_lims)+timeZone_offset ] = np.min(xAxis_lims)+timeZone_offset; #limit the xTime to the plotted times
        # xTime[ xTime > np.max(xAxis_lims)+timeZone_offset ] = np.max(xAxis_lims)+timeZone_offset;
        # xTime[0] = xTime[1]; #set that to match (for plotting niceness)
        # xTime[-1] = xTime[-2]; #set that to match (for plotting niceness)
        yDayNite = np.delete(yDayNite, (xTime > np.max(xAxis_lims)+timeZone_offset) | (xTime < np.min(xAxis_lims)+timeZone_offset) ); #remove out of bounds stuff
        xTime = np.delete(xTime, (xTime > np.max(xAxis_lims)+timeZone_offset) | (xTime < np.min(xAxis_lims)+timeZone_offset) ); #remove out of bounds stuff
        #fill in edges so they just end
        yDayNite = np.append(np.insert(yDayNite,0,yDayNite[0]),yDayNite[-1]);
        xTime = np.append(np.insert(xTime,0,np.min(xAxis_lims)+timeZone_offset),np.max(xAxis_lims)+timeZone_offset);
        
        
        #FIFTH STEP: ACTUALLY PLOTTING
        divider2 = make_axes_locatable(ax[2]); #prep to add an axis
        cax2 = divider2.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
        cax2.set_visible(False); #mkae it invisible so it matches the other plots in width
        ax[2].plot(xTime,yDayNite,color='xkcd:black',linewidth=settings_plot['line width']['thicc'], antialiased=True);
        if( timeZone_DSTnUTCOffset != 0 ):
            strang = timeZone_offset_str+' (Daylight Savings) '+timeZone_name+' Time Zone';
        else:
             strang = timeZone_offset_str+' '+timeZone_name+' Time Zone';
        #END IF
        # if( yDayNite[-1] == yDayNite[-2] ):
        #     #hacky offset so day or night isn't printed off the plot (or skipped)
        #     offset = 3;
        # else:
        #     offset = 1;
        # #END IF
        for i in range(0,xTime.size,2):
            if( yDayNite[i] == 1 ):
                ax[2].text( (xTime[i+1]-xTime[i])/2+xTime[i]-1.3, \
                    0.22,'Day', color='k', fontproperties=settings_plot['font title FM']); #print the text saying the day or nite
            else:
                ax[2].text( (xTime[i+1]-xTime[i])/2+xTime[i]-1.9, \
                   0.28, 'Night', color='k', fontproperties=settings_plot['font title FM']); #print the text saying the day or nite
            #END IF
        #END FOR i
        ax[2].set_xlabel('Local Time [hr] | '+strang,fontproperties=settings_plot['font axis label FM'])
        ax[2].set_xticks(xAxis_ticks+timeZone_offset); #set x axis ticks
        # ax[2].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
        ax[2].set_yticklabels([]); #remove y labels for day nite plot
        ax[2].set_yticks([]); #remove y tick marks for day nite plot
        ax[2].set_xlim( (xAxis_lims[0]+timeZone_offset,xAxis_lims[1]+timeZone_offset) ); #set x axis limits
        ax[2].set_ylim( (0,1) ); #set y lims to 0 and 1
        ax[2].spines['left'].set_visible(False); #turn off box lines
        ax[2].spines['right'].set_visible(False); #turn off box lines
        ax[2].spines['top'].set_visible(False); #turn off box lines
        # ax[2].grid(b=True, which='major', axis='x', color='xkcd:light grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
        ax[2].text( letteringPositionX, letteringPositionY, 'c.', color='r', fontproperties=settings_plot['font grandiose FM'], transform=ax[2].transAxes); #print the text saying the day or nite
    #END IF
    
    #----- FINISH UP ----- 
    figFitter(fig); #fit that fig fast
    # fig.subplots_adjust(left = 0.045, right = 0.95, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
    if( FLG_fancyPlot == 0 ):
        plt.show(); #req to make plot show up
    else:
        figFileName = jointer(settings_paths['fancyPlots'],fileName+settings_plot['save file type']);
        if( settings_plot['save file type'].lower() != '.pdf' ):
            fig.savefig(figFileName); #save the figure
        else:
            from matplotlib.backends.backend_pdf import PdfPages
            pp = PdfPages(figFileName);
            pp.savefig(fig);
            pp.close();
        #END IF
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF