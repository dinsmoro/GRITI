"""
GOAL: Plot only ISR POPL HP RTI
RD on 6/4/20

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made

"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tick
import timezonefinder
from datetime import datetime
import pytz
from astroplan import Observer
import astropy.units as astroUnits
from astropy.time import Time
from urllib.request import urlopen
from subfun_strstr import strstr
from subfun_dateORdayNum_to_fullRange import subfun_dateORdayNum_to_fullRange


def GRITI_ISR_Haystack_fancyPlot_POPL_HP_wDayNite(Zenith_time,Zenith_height,Zenith_POPL_hp, \
        MISA_time,MISA_height,MISA_POPL_hp,filter_cutoffPeriod,ISR_RTI_heightLimValues,ISR_POPL_plotLimValu, \
        dateRange_full,dateRange_dayNum_full,dateRange_dayNum_zeroHr,dateRange_zeroHr, dateRange_zeroHr_monthName, \
        dateRange_zeroHr_dayPostfix, latLong_ref,FONT_grandioseFM, FONT_titleFM, FONT_axisLabelFM, FONT_axisTick, \
        PLOT_lineWidth, folder, journal_width_2C,journal_height_max,journal_dpi):
    print('MAKING FANCY PLOT: ISR_Haystack_fancyPlot_POPL_HP_wDayNite IN fancyPlot FOLDER'); #report since you won't see anything
    
    ISR_m3tocc = 100**3; #1 m^3 is 100^3 cm^3 and 1 cm^3 is 1 cc
    ISR_POPL_plotLimValu = ISR_POPL_plotLimValu/ISR_m3tocc; #adjust this plot limit value
    Zenith_POPL_hp = Zenith_POPL_hp/ISR_m3tocc; #adjust POPL values
    MISA_POPL_hp = MISA_POPL_hp/ISR_m3tocc; #adjust POPL values
    
    letteringPositionX = -0.09; #set the X position of the lettering (e.g., a. b. c. ...)
    letteringPositionY = 0.88; #set the X position of the lettering (e.g., a. b. c. ...)
    
    # def cbarFormatter(x, y):
    #     cbarVal = '{:1.0e}'.format(x).replace('+',''); #get it in basic scientific format
    #     if( cbarVal[cbarVal.find('e')+1] == '0' ):
    #         cbarVal = cbarVal[:cbarVal.find('e')+1]+cbarVal[cbarVal.find('e')+2:]; #remove the 0 in '2e05' or something like that
    #     #END IF
    #     if( cbarVal[0] == '0' ):
    #         cbarVal = '0'; #just set to 0
    #     #END IF
    #     return cbarVal
    # #END DEF
    
    xAxisLims = ((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24); #get those xAxis limits, ensures everything is aligned
    ISR_RTI_heightLimValues = np.array(ISR_RTI_heightLimValues); #convert to array so it's not weak python tuples
    ISR_RTI_heightLimValues[0] = np.min(Zenith_height); #move the height limit to the actual lower-limit
    
    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth['thicc']; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth['double plus']; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth['plus']; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth['regular plus']; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth['regular']; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth['smol']; #get the line widths
    
    #prep the plot
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig = plt.figure(figsize=(14,10.5),dpi=journal_dpi);
    # gridr = gridspec.GridSpec(nrows=7, ncols=1, figure=fig);
    # fig.subplots_adjust(hspace=0.75)
    gridr = gridspec.GridSpec(nrows=2, ncols=1, height_ratios = [20, 1], figure=fig); #prep for a nested gridspec (it's 2) and note the ratios of the plots (8 to 1)
    gridr1 = gridspec.GridSpecFromSubplotSpec(nrows=20, ncols=1, subplot_spec = gridr[0], hspace = 5.0)
    gridr2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec = gridr[1], hspace = 0.0)
    # gridr.update(hspace=0.05); # set the spacing between axes.
    fig.add_subplot(gridr1[0:10]); #RTI plots are 2 tall
    # gridr.update(hspace=0.80); # set the spacing between axes.
    fig.add_subplot(gridr1[10:20]); #RTI plots are 2 tall
    # gridr.update(hspace=0.90); # set the spacing between axes.
    fig.add_subplot(gridr2[0]); #dayNite plot is 1 tall
    ax = fig.axes; #get a list of the axes
    
    #Remove the aspect ratio so it fills the screen better
    ax[0].set_aspect('auto'); #set to auto for all axes
    ax[1].set_aspect('auto'); #set to auto for all axes
    ax[2].set_aspect('auto'); #set to auto for all axes
    
    #!!!-----Plot ISR POPL HP results as a RTI-----!!!
    #Plot just the ISR POPL HP results
    #ZENITH STUFF
    divider0 = make_axes_locatable(ax[0]); #prep to add an axis
    cax0 = divider0.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    pltHelprX, pltHelprY = np.meshgrid( (Zenith_time - dateRange_dayNum_zeroHr[1])*24, Zenith_height);
    im0 = ax[0].pcolormesh(pltHelprX , pltHelprY , Zenith_POPL_hp , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu), \
        cmap='gray' , shading='gouraud', antialiased=True); # pseudocolor plot "stretched" to the grid
    cbar0 = fig.colorbar(im0, cax=cax0, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar0.mappable.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
    cbar0.set_label(r"N$_{e^-}$ [$e^-/cc$]"); #tabel the colorbar
    cbar0.ax.tick_params(labelsize=FONT_axisTick);
    # cax0.yaxis.set_major_formatter(tick.FuncFormatter(cbarFormatter)); #force a rounded format
    cax0.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    # string_Title = 'Zenith POPL with High-pass Filter '+str(filter_cutoffPeriod)+' hr Cutoff'; #create mecha title
    # ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    # ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('Height [km]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 4 , \
            4); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    # ax[0].set_xticklabels([]); #if statement to remove x axis labels except for the last line
    ax[0].set_xlim( xAxisLims ); #set x axis limits
    yAxisTicks = np.arange(100,ISR_RTI_heightLimValues[1]+50,50); #get the y axis tick locations
    ax[0].set_yticks(yAxisTicks); #set x axis ticks
    ax[0].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text saying the day or nite
    
    
    #MISA STUFF
    divider1 = make_axes_locatable(ax[1]); #prep to add an axis
    cax1 = divider1.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    
    pltHelprX, pltHelprY = np.meshgrid( (MISA_time - dateRange_dayNum_zeroHr[1])*24, MISA_height);
    im1 = ax[1].pcolormesh(pltHelprX , pltHelprY , MISA_POPL_hp , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu), \
        cmap='gray', shading='gouraud', antialiased=True); # pseudocolor plot "stretched" to the grid
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar1.mappable.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
    #cbar.ax.set_yticklabels(np.round(np.linspace(-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu,5), len(str(ISR_POPL_plotLimValu).split('.')[1])+1 )); #create useful tick marks
    cbar1.set_label("N$_{e^-}$ [$e^-/cc$]"); #tabel the colorbar
    cbar1.ax.tick_params(labelsize=FONT_axisTick);
    cax1.yaxis.label.set_font_properties(FONT_axisLabelFM);
    # cax1.yaxis.set_major_formatter(tick.FuncFormatter(cbarFormatter)); #force a rounded format
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks    
    
    # string_Title = 'MISA POPL with High-pass Filter '+str(filter_cutoffPeriod)+' hr Cutoff'; #create mecha title
    # ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[1].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel('Height [km]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    # xAxisTicks = np.arange( (np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
    #         (np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
    #         2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks
    # ax[1].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    ax[1].set_xlim( xAxisLims ); #set to Zenith to match
    # yAxisTicks = np.arange(100,ISR_RTI_heightLimValues[1]+50,50); #get the y axis tick locations
    ax[1].set_yticks(yAxisTicks); #set x axis ticks
    ax[1].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=FONT_grandioseFM, transform=ax[1].transAxes); #print the text saying the day or nite
    
    
    #!!!DAY NITE PLOT STUFF!!!
    Xaxisvar_min = (Zenith_time[0] - dateRange_dayNum_zeroHr[1])*24; #UT hr, min time to compare to
    Xaxisvar_max = (Zenith_time[-1] - dateRange_dayNum_zeroHr[1])*24; #UT hr, max time to compare to
#         #below is alt. full time
#         Xaxisvar_min = (min(timeUnique) - dateRange_zeroHr(2))*24; #UT hr, min time to compare to
#         Xaxisvar_max = (max(timeUnique) - dateRange_zeroHr(2))*24; #UT hr, min time to compare to

    #FIRST BATTLE: REGION WHERE LAT/LONG IS   
    tf = timezonefinder.TimezoneFinder(); #prep the time zone finder function thing
    dayNite_timeZoneID = tf.certain_timezone_at(lat=latLong_ref[0], lng=latLong_ref[1]); #use it to find the time zone
    if dayNite_timeZoneID is None:
        #use geonames site as a backup
        url = 'http://api.geonames.org/timezone?lat='+str(latLong_ref[0])+'&lng='+str(latLong_ref[1])+'&username=razzluhdzuul'; #create link for lat/long
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
    
    # dateRange_dates = sFUN_dayNumber_to_Date_MULTIPLE(dateRange); #convert day# range to yr/mon/day range with all inbetween filled in
    # dateRange_extended = sFUN_dayNumber_to_Date_MULTIPLE(dateRange,1); #convert day# range to yr/day# flushed out
    # dayNite_DSTactive = isdst(datetime(dateRange_dates[1,1],dateRange_dates[1,2],dateRange_dates[:,3],'TimeZone',dayNite_timeZoneID)); #get DST active or not for each day
    
    #SECOND STEP: CALC LOCAL SUNRISE/SUNSET TIMES
    #will adjust for DST later - calcs done in UT/GMT
    # based on calc steps in https://www.mathworks.com/examples/matlab/community/21093-estimating-sunrise-and-sunset
    # long_corrected = 4*(latLong_ref[1] - 15*dayNite_UToffset); #calc corrected longitude, for sunrise/sunset time

    # dayNite_B = 360*(dateRange_dayNum_full[:,1] - 81)/365; #some sort of angle based on days and stuff
    # dayNite_EoT_corrected = 9.87*np.sin(2*dayNite_B*np.pi/180) - 7.53*np.cos(dayNite_B*np.pi/180) - 1.5*np.sin(dayNite_B*np.pi/180); #eq for Time Correction
    # dayNite_solar_corrected = long_corrected + dayNite_EoT_corrected; #min, solar time correction - for noon

    # dayNite_solar_declination = np.arcsin(np.sin(23.45*np.pi/180)*np.sin(360*(dateRange_dayNum_full[:,1] - 81)/365)*np.pi/180); #deg, solar declination

    # dayNite_sunrise = 12 - np.arccos(-np.tan(latLong_ref[0]*np.pi/180)*np.tan(dayNite_solar_declination*np.pi/180))/15 - dayNite_solar_corrected/60; #hr, sunrise time
    # dayNite_sunset = 12 + np.arccos(-np.tan(latLong_ref[0]*np.pi/180)*np.tan(dayNite_solar_declination*np.pi/180))/15 - dayNite_solar_corrected/60; #hr, sunrise time

    ISR = Observer(longitude=latLong_ref[1]*astroUnits.deg, latitude=latLong_ref[0]*astroUnits.deg, \
        elevation=0*astroUnits.m, name='ISR', timezone=dayNite_timeZoneID); #make an atroplan observer at the ISR location
    timeZoneObj = pytz.timezone(dayNite_timeZoneID); #make a timezone object
    timeZoneObj_UTC = pytz.timezone('UTC'); #make a timezone object
    #Daylight savings local time fix
    dayNite_sunrise = np.zeros( (dateRange_full.shape[0],) ); #preallocate
    dayNite_sunset = np.zeros( (dateRange_full.shape[0],) ); #preallocate
    for i in range(0,dateRange_full.shape[0]):
        dateRange_timeObj = Time(str(dateRange_full[i,0])+'-'+str(dateRange_full[i,1])+'-'+str(dateRange_full[i,2])+'T12:00:00', format='isot', scale='utc'); #make an astropy time object
        
        #Run through the days
        dayNite_sunriseTemp = ISR.sun_rise_time(dateRange_timeObj, which='nearest', horizon=0*astroUnits.deg, n_grid_points=150).to_value('isot'); #get the sunrise time
        dayNite_sunriseObj = timeZoneObj_UTC.localize(datetime.strptime(dayNite_sunriseTemp, '%Y-%m-%dT%H:%M:%S.%f')).astimezone(timeZoneObj); #create datetime object, set it to the UTC time zone (which it is, datetime just doesn't know), then convert it to local time zone
        #convert to decimal hour
        dayNite_sunriseTemp = dayNite_sunriseObj.isoformat(); #convert to a string that's easily readable
        index_start = strstr(dayNite_sunriseTemp,'T')[0]; #get where the hour starts
        index_colon = strstr(dayNite_sunriseTemp,':'); #get where the colons are
        index_dash = strstr(dayNite_sunriseTemp,'-')[-1]; #get where the last -, end of time
        dayNite_sunrise[i] = dateRange_dayNum_full[i,1] + np.int64(dayNite_sunriseTemp[index_start+1:index_colon[0]])/24 + \
            np.int64(dayNite_sunriseTemp[index_colon[0]+1:index_colon[1]])/1440 + np.float64(dayNite_sunriseTemp[index_colon[1]+1:index_dash])/86400; #make into day units
        
        #Run through the days
        dayNite_sunsetTemp = ISR.sun_set_time(dateRange_timeObj, which='nearest', horizon=0*astroUnits.deg, n_grid_points=150).to_value('isot'); #get the sunrise time
        dayNite_sunsetObj = timeZoneObj_UTC.localize(datetime.strptime(dayNite_sunsetTemp, '%Y-%m-%dT%H:%M:%S.%f')).astimezone(timeZoneObj); #create datetime object, set it to the UTC time zone (which it is, datetime just doesn't know), then convert it to local time zone
        #convert to decimal hour
        dayNite_sunsetTemp = dayNite_sunsetObj.isoformat(); #convert to a string that's easily readable
        index_start = strstr(dayNite_sunsetTemp,'T')[0]; #get where the hour starts
        index_colon = strstr(dayNite_sunsetTemp,':'); #get where the colons are
        index_dash = strstr(dayNite_sunsetTemp,'-')[-1]; #get where the last -, end of time
        dayNite_sunset[i] = dateRange_dayNum_full[i,1] + np.int64(dayNite_sunsetTemp[index_start+1:index_colon[0]])/24 + \
            np.int64(dayNite_sunsetTemp[index_colon[0]+1:index_colon[1]])/1440 + np.float64(dayNite_sunsetTemp[index_colon[1]+1:index_dash])/86400; #make into day units
    #END FOR i
    dayNite_DSToffset_str = dayNite_sunriseTemp[index_dash:]; #get the DST time offset
    if( strstr(dayNite_DSToffset_str,':').size > 0 ):
        if( dayNite_DSToffset_str[strstr(dayNite_DSToffset_str,':')[0]+1:] == '00' ):
            dayNite_DSToffset_str = dayNite_DSToffset_str[:strstr(dayNite_DSToffset_str,':')[0]]; #remove the :00 if it's just that
            if( dayNite_DSToffset_str[1] == '0' ):
                dayNite_DSToffset_str = dayNite_DSToffset_str.replace('0',''); #remove the 0 that was extraneous
            #END IF
            dayNite_DSToffset = np.int64(dayNite_DSToffset_str); #get the number version
        #END IF
        else:
            dayNite_DSToffset = np.float64(dayNite_DSToffset_str[:strstr(dayNite_DSToffset_str,':')[0]]) + \
                np.int64(dayNite_DSToffset_str[strstr(dayNite_DSToffset_str,':')[0]+1:])/60; #convert to an hour decimal
        #END IF
    #END IF
    dayNite_timeZoneName = dayNite_sunsetObj.tzname(); #get the time zone name (like 'EST' or 'EDT' depending on standard or daylight savings time)
    dateNite_DSTnUTCOffset = dayNite_sunsetObj.dst().total_seconds()/3600; #get the time offset
    if( np.mod(dateNite_DSTnUTCOffset,1) == 0 ):
        dateNite_DSTnUTCOffset = np.int64(dateNite_DSTnUTCOffset); #convert to integer
    #END IF

    #THIRD STEP: PREP FOR PLOTTING BY ALIGNING TIMES, MAKING PLOT VARIABLES
    Xaxisvar_min = Xaxisvar_min + dayNite_DSToffset; #hr local, UT time of -12 conv. to local
    Xaxisvar_max = Xaxisvar_max + dayNite_DSToffset; #hr local, UT time of -12 conv. to local

    dayNite_sunrise = (dayNite_sunrise - dateRange_dayNum_zeroHr[1])*24; #hr, convert so 0 hr is in the middle (and convert from days to hours)
    dayNite_sunset = (dayNite_sunset - dateRange_dayNum_zeroHr[1])*24; #hr, convert so 0 hr is in the middle (and convert from days to hours)

    xTime = np.sort( np.concatenate( (dayNite_sunrise,dayNite_sunset) ) ); #hr, xvar to plot against
    yDayNite = np.ones(xTime.shape); #prep if day or night, set all to day
    for i in range(0,dayNite_sunset.size):
        yDayNite[xTime == dayNite_sunset[i]] = 0; #set sunset times to sunset
    #END FOR i
    xTime = xTime.repeat(2); #interleave repeated values
    yDayNite = np.roll(yDayNite.repeat(2),1) #interleave repeated values and circular shift by 1
    yDayNite[0] = yDayNite[1]; #set that to match (for plotting niceness)
    yDayNite[-1] = yDayNite[-2]; #set that to match (for plotting niceness)
    
    xTime[ xTime < np.min(xAxisLims)+dayNite_DSToffset ] = np.min(xAxisLims)+dayNite_DSToffset; #limit the xTime to the plotted times
    xTime[ xTime > np.max(xAxisLims)+dayNite_DSToffset ] = np.max(xAxisLims)+dayNite_DSToffset;
    # xTime[0] = xTime[1]; #set that to match (for plotting niceness)
    # xTime[-1] = xTime[-2]; #set that to match (for plotting niceness)
    
    #FOURTH STEP: ACTUALLY PLOTTING
    divider2 = make_axes_locatable(ax[2]); #prep to add an axis
    cax2 = divider2.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    cax2.set_visible(False); #mkae it invisible so it matches the other plots in width
    ax[2].plot(xTime,yDayNite,color='xkcd:black',linewidth=PLOT_lineWidthThicc, antialiased=True);
    if( dateNite_DSTnUTCOffset != 0 ):
        strang = dayNite_DSToffset_str+' (Daylight Savings) '+dayNite_timeZoneName+' Time Zone';
    else:
         strang = dayNite_DSToffset_str+' '+dayNite_timeZoneName+' Time Zone';
    #END IF
    # if( yDayNite[-1] == yDayNite[-2] ):
    #     #hacky offset so day or night isn't printed off the plot (or skipped)
    #     offset = 3;
    # else:
    #     offset = 1;
    # #END IF
    for i in range(1,xTime.size-2,2):
        if( yDayNite[i] == 1 ):
            ax[2].text( (xTime[i+1]-xTime[i])/2+xTime[i]-1.3, \
                0.25,'Day', color='k', fontproperties=FONT_titleFM); #print the text saying the day or nite
        else:
            ax[2].text( (xTime[i+1]-xTime[i])/2+xTime[i]-1.9, \
               0.25, 'Night', color='k', fontproperties=FONT_titleFM); #print the text saying the day or nite
        #END IF
    #END FOR i

    ax[2].set_xlabel('Local Time [hr] | '+strang,fontproperties=FONT_axisLabelFM)
    ax[2].set_xticks(xAxisTicks+dayNite_DSToffset); #set x axis ticks
    # ax[2].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    ax[2].set_yticklabels([]); #remove y labels for day nite plot
    ax[2].set_yticks([]); #remove y tick marks for day nite plot
    ax[2].set_xlim( (xAxisLims[0]+dayNite_DSToffset,xAxisLims[1]+dayNite_DSToffset) ); #set x axis limits
    ax[2].set_ylim( (0,1) ); #set y lims to 0 and 1
    ax[2].spines['left'].set_visible(False); #turn off box lines
    ax[2].spines['right'].set_visible(False); #turn off box lines
    ax[2].spines['top'].set_visible(False); #turn off box lines
    # ax[2].grid(b=True, which='major', axis='x', color='xkcd:light grey',linewidth=PLOT_lineWidthSmol); #sets major axis grid lines to be on
    ax[2].text( letteringPositionX, letteringPositionY, 'c.', color='r', fontproperties=FONT_grandioseFM, transform=ax[2].transAxes); #print the text saying the day or nite
    
    fig.subplots_adjust(left = 0.08, right = 0.910, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig(folder[3]+'\\ISR_RTI_Z&M_wDayNite.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff