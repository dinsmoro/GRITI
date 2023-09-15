"""
Writes the image to a picture
"""
import numpy as np
from datetime import datetime
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import cartopy
from cartopy.feature.nightshade import Nightshade
from Code.GRITI_plotHelper_nightshader import nightshader
from Code.GRITI_plotHelper_area_init import GRITI_plotHelper_area_init
from Code.GRITI_keo_keogrammer import GRITI_keo_keogrammer
from Code.subfun_timeMatch import subfun_timeMatch
from Code.GRITI_movieMaker_subfun_figMaker import GRITI_movieMaker_subfun_figMaker
from Code.GRITI_movieMaker_subfun_dataGridder import GRITI_movieMaker_subfun_dataGridder
from Code.GRITI_movieMaker_subfun_dayniteCalc import GRITI_movieMaker_subfun_dayniteCalc
from Code.subfun_textNice import textNice
from Code.subfun_figFitter import figFitter
from Code.subfun_time_to_dateRange import subfun_time_to_dateRange
from Code.subfun_dayNum_to_date import subfun_dayNum_to_date
from Code.subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_movieMaker_subfun_imageWriter(data, dates, settings, movie_writer, movie_figOffsets, movie_title_yOffset):
    #-------------Multiple Movie Algs for Difference Scenarios-----------------
    #================================================================= MOVIE TYPE 0 =================================================================
    if( settings['movie']['movie type'] == 0 ): #moving gif implementation (scatter points)
        
        #-------------------------Start Making Pictures------------------------
        for i in range(0,data['TEC']['time unique'].size):
    
            #----------------Corral the data to the right place----------------
            k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
            
            #----------------------------Tack on Title-------------------------
            string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
            
            ax.set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                
            #-------------------Starting the Plotting--------------------------
            if( i == 0 ): #first run preps axes, color bars, etc.
                
                TEC_latLongMapped = geoMap(data['TEC']['long'][k],data['TEC']['lat'][k]); #convert the lat/long arcdeg to the current map coordinates
                
                #Do the TEC plotting
                imTEC = ax.scatter(TEC_latLongMapped[0],TEC_latLongMapped[1],s=settings['map']['TEC scatter size'],c=data['TEC']['dTEC'][k],cmap='jet', vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']), zorder=5);
                cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));

                if( settings['movie']['day nite line'] == 1 ):
                    movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                    #Plot the sunrise/sunset terminators
                    dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                    hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                    #constants will use all the time - only for plotting of day/nite line so minor importance
                    bbox = ax.get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
                    plot_ratio = bbox.width/bbox.height; #get the plot ratio, will use it to fix up the angle
                    dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                    dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                    dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                    if( hLgd_FLG_day > 0 ): #only do work if it is there
                        #calc all that day/nite stuff in one function to keep it from getting cluttered
                        dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                        dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                        imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                        if( settings['movie']['day nite text'] == 1 ):
                            dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                            textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                        #END IF
                    #END IF
                    
                    dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                    hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                    if( hLgd_FLG_nite > 0 ): #only do work if it is there
                        #calc all that day/nite stuff in one function to keep it from getting cluttered
                        dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                        dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                        imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                        if( settings['movie']['day nite text'] == 1 ):
                            dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                            textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                        #EMD IF
                    #END IF
                #END IF
                
                #Now drawing line of interest
#                imMillstone = plot(fig1Axes,settings['map']['site coords'][0,1],settings['map']['site coords'][0,0],settings['map']['site marker type'],'Color',settings['map']['site marker color'],'MarkerSize',settings['map']['site marker size']); #plots a point with a red big *
                
                Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                imMillstone = ax.plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                
            else:
                #Just plotting now - prep done, hopefully speeds it!
                TEC_latLongMapped = geoMap(data['TEC']['long'][k],data['TEC']['lat'][k]); #convert the lat/long arcdeg to the current map coordinates

                #Do the TEC plotting
                imTEC = ax.scatter(TEC_latLongMapped[0],TEC_latLongMapped[1],s=settings['map']['TEC scatter size'],c=data['TEC']['dTEC'][k],cmap='jet', vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']), zorder=5);

                if( settings['movie']['day nite line'] == 1 ):
                    #Plot the sunrise/sunset terminators
                    dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                    hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                    if( hLgd_FLG_day > 0 ): #only do work if it is there
                        #calc all that day/nite stuff in one function to keep it from getting cluttered
                        dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                        dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                        imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                        if( settings['movie']['day nite text'] == 1 ):
                            dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                            textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                        #END IF
                    #END IF
                    
                    dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                    hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                    if( hLgd_FLG_nite > 0 ): #only do work if it is there
                        #calc all that day/nite stuff in one function to keep it from getting cluttered
                        dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                        dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                        imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                        if( settings['movie']['day nite text'] == 1 ):
                            dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                            textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                        #END IF
                    #END IF
                #END IF
                
                #Now drawing line of interest
#                imMillstone = ax.plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], markerfacecolor=settings['map']['site marker color'] ,markersize=settings['map']['site marker size']); #plot this
            #END IF
        
        
            #-----------------------Create Movie/GIF---------------------------
            #Makes the gif now
            plt.draw();
            
            movie_writer.grab_frame(); #get the frame and save it
                        
            #-------------------Clean up for re-use----------------------------
            #if forget one (like hOverlay) slows it way down after many plots
            imTEC.remove();
#            imOverlay.pop(0).remove();
#            imMillstone.pop(0).remove();
            if( settings['movie']['day nite line'] == 1 ):
                if(hLgd_FLG_day > 0): #only delete if it is there
                    imDayNite_day.pop(0).remove();
                    if( settings['movie']['day nite text'] == 1 ):
                        textDayNite_day.remove();
                    #END IF
                #END IF
                if(hLgd_FLG_nite > 0): #only delete if it is there
                    imDayNite_nite.pop(0).remove();
                    if( settings['movie']['day nite text'] == 1 ):
                        textDayNite_nite.remove();
                    #END IF
                #END IF
            #END IF
            
        #END FOR i
    #================================================================= MOVIE TYPE 1 =================================================================
    elif(settings['movie']['movie type'] == 1): #stationary data points (through averaging)
        if( np.isscalar(settings['TEC']['plot lim']) == 1 ):
            gif_TEC_plotLimValu = np.array( (-settings['TEC']['plot lim'],settings['TEC']['plot lim']) ); #make it a vector
        else:
            gif_TEC_plotLimValu = settings['TEC']['plot lim']; #keep it the same
        #END IF
    
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
        
            #-------------------------Start Making Pictures------------------------
            for i in range(0,data['TEC']['time unique'].size):
        
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
                
                ax.set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                    
                #-------------------Starting the Plotting--------------------------
                if( i == 0 ): #first run preps axes, color bars, etc.
                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax.pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap=settings['TEC']['colormap'],zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label(settings['TEC']['name']+settings['TEC']['units']); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),6)); #create useful tick marks
    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax.get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i]/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax, zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                            hSun = ax.add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax.set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest
    #                imMillstone = plot(fig1Axes,settings['map']['site coords'][0,1],settings['map']['site coords'][0,0],settings['map']['site marker type'],'Color',settings['map']['site marker color'],'MarkerSize',settings['map']['site marker size']); #plots a point with a red big *
                    
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax.plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                    
                    figFitter(fig); #make sure everything is fit well
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax.pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap=settings['TEC']['colormap'],zorder=5); # pseudocolor plot "stretched" to the grid
    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i]/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax, zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                            hSun = ax.add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax.set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest
    #                imMillstone = ax.plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], markerfacecolor=settings['map']['site marker color'] ,markersize=settings['map']['site marker size']); #plot this
                #END IF
            
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                movie_writer.grab_frame(); #get the frame and save it
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
    #            imOverlay.pop(0).remove();
    #            imMillstone.pop(0).remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                elif( settings['movie']['day nite line'] == 2): #only delete if it is there
                    for middleman in dayNite_shade.collections: #this needs an extra helping hand
                        middleman.remove();
                #END IF
                if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
                    hSun.remove();
                #END IF
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 2 =================================================================
    elif(settings['movie']['movie type'] == 2): #stationary data points (through averaging) + Zenith ISR overlay
        
        #-----------------------Clip Time to ISR Zenith------------------------
        gif_Zenith_time_start = np.where( np.min(np.abs(np.min(data['ISR']['mill']['zenith']['time']) - data['TEC']['time unique'])) == np.abs(np.min(data['ISR']['mill']['zenith']['time']) - data['TEC']['time unique']) )[0][0]; #match Zenith start time
        gif_Zenith_time_end = np.where( np.min(np.abs(np.max(data['ISR']['mill']['zenith']['time']) - data['TEC']['time unique'])) == np.abs(np.max(data['ISR']['mill']['zenith']['time']) - data['TEC']['time unique']) )[0][0]; #match Zenith end time
        
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
        
            #-------------------------Start Making Pictures------------------------
            for i in range(gif_Zenith_time_start,gif_Zenith_time_end+1):
        
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                                
                #-------------------Starting the Plotting--------------------------
                if( i == gif_Zenith_time_start ): #first run preps axes, color bars, etc.
                                   
                    #tack on 2nd title that doesn't change
                    string_title = 'Zenith SNR '+str(settings['spectra']['filter cutoff period'])+' Hr Highpass Filtered'; #create mecha title
                
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    #Do the ISR plotting
                    pltHelprX, pltHelprY = np.meshgrid( (data['ISR']['mill']['zenith']['time'] - dates['date range zero hr dayNum'][1]*86400)/3600, data['ISR']['mill']['zenith']['height']);
                    imISR = ax[1].pcolormesh(pltHelprX , pltHelprY , data['ISR']['mill']['zenith']['SNR hp'] , vmin=-settings['ISR']['mill']['plot lim SNR'] , vmax=settings['ISR']['mill']['plot lim SNR'] , cmap='gray' , shading='gouraud', zorder = 1); # pseudocolor plot "stretched" to the grid
                    cbar2 = movie_dict['fig'].colorbar(imISR, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cbar2.mappable.set_clim(vmin=-settings['ISR']['mill']['plot lim SNR'], vmax=settings['ISR']['mill']['plot lim SNR']);
    #                cbar2.ax.set_yticklabels(np.round(np.linspace(-settings['ISR']['mill']['plot lim SNR'],settings['ISR']['mill']['plot lim SNR'],5), len(str(settings['ISR']['mill']['plot lim SNR']).split('.')[1])+1 )); #create useful tick marks
                    cax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar2.set_label("SNR (unitless)"); #tabel the colorbar
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    ax[1].set_xlabel(subfun_monthNum_to_word(dates['date range'][0,1])[0]+" "+str(dates['date range'][0,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][0,0])+" to "+subfun_monthNum_to_word(dates['date range'][1,1])[0]+" "+str(dates['date range'][1,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][1,0]),fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[1].set_ylabel('Height (km)',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((np.min(data['ISR']['mill']['zenith']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((np.min(data['ISR']['mill']['zenith']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((np.max(data['ISR']['mill']['zenith']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((np.max(data['ISR']['mill']['zenith']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    ax[1].set_xlim( ((np.min(data['ISR']['mill']['zenith']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600 , (np.max(data['ISR']['mill']['zenith']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600) ); #set x axis limits
                    
                    ax[1].set_ylim( settings['ISR']['mill']['height lim'] ); #set y axis limits
                    
                    #Draw time line on bottom ISR plot
                    imISRLine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( np.min(settings['ISR']['mill']['height lim']),np.max(settings['ISR']['mill']['height lim']))) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        bbox = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
                        plot_ratio = bbox.width/bbox.height; #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest               
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
    
                    #Draw time line on bottom ISR plot
                    imISRLine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( np.min(settings['ISR']['mill']['height lim']),np.max(settings['ISR']['mill']['height lim']))) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    #END IF
                    
                #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                movie_writer.grab_frame(); #get the frame and save it
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imISRLine.pop(0).remove();
    #            imOverlay.pop(0).remove();
    #            imMillstone.pop(0).remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 3 =================================================================
    elif(settings['movie']['movie type'] == 3): #stationary data points (through averaging) + MISA ISR overlay
        
        #-----------------------Clip Time to ISR MISA------------------------
        gif_MISA_time_start = np.where( np.min(np.abs(np.min(data['ISR']['mill']['MISA']['time']) - data['TEC']['time unique'])) == np.abs(np.min(data['ISR']['mill']['MISA']['time']) - data['TEC']['time unique']) )[0][0]; #match MISA start time
        gif_MISA_time_end = np.where( np.min(np.abs(np.max(data['ISR']['mill']['MISA']['time']) - data['TEC']['time unique'])) == np.abs(np.max(data['ISR']['mill']['MISA']['time']) - data['TEC']['time unique']) )[0][0]; #match MISA end time
        
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
        
            #-------------------------Start Making Pictures------------------------
            for i in range(gif_MISA_time_start,gif_MISA_time_end+1):
        
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                                
                #-------------------Starting the Plotting--------------------------
                if( i == gif_MISA_time_start ): #first run preps axes, color bars, etc.
                                   
                    #tack on 2nd title that doesn't change
                    string_title = 'MISA SNR '+str(settings['spectra']['filter cutoff period'])+' Hr Highpass Filtered'; #create mecha title
                
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    #Do the ISR plotting
                    pltHelprX, pltHelprY = np.meshgrid( (data['ISR']['mill']['MISA']['time'] - dates['date range zero hr dayNum'][1]*86400)/3600, data['ISR']['mill']['MISA']['height']);
                    imISR = ax[1].pcolormesh(pltHelprX , pltHelprY , data['ISR']['mill']['MISA']['SNR hp'] , vmin=-settings['ISR']['mill']['plot lim SNR'] , vmax=settings['ISR']['mill']['plot lim SNR'] , cmap='gray' , shading='gouraud', zorder = 1); # pseudocolor plot "stretched" to the grid
                    cbar2 = movie_dict['fig'].colorbar(imISR, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cbar2.mappable.set_clim(vmin=-settings['ISR']['mill']['plot lim SNR'], vmax=settings['ISR']['mill']['plot lim SNR']);
    #                cbar2.ax.set_yticklabels(np.round(np.linspace(-settings['ISR']['mill']['plot lim SNR'],settings['ISR']['mill']['plot lim SNR'],5), len(str(settings['ISR']['mill']['plot lim SNR']).split('.')[1])+1 )); #create useful tick marks
                    cax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar2.set_label("SNR (unitless)"); #tabel the colorbar
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    ax[1].set_xlabel(subfun_monthNum_to_word(dates['date range'][0,1])[0]+" "+str(dates['date range'][0,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][0,0])+" to "+subfun_monthNum_to_word(dates['date range'][1,1])[0]+" "+str(dates['date range'][1,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][1,0]),fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[1].set_ylabel('Height (km)',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((np.min(data['ISR']['mill']['MISA']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((np.min(data['ISR']['mill']['MISA']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((np.max(data['ISR']['mill']['MISA']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((np.max(data['ISR']['mill']['MISA']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    ax[1].set_xlim( ((np.min(data['ISR']['mill']['MISA']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600 , (np.max(data['ISR']['mill']['MISA']['time'])-dates['date range zero hr dayNum'][1]*86400)/3600) ); #set x axis limits
                    
                    ax[1].set_ylim( settings['ISR']['mill']['height lim'] ); #set y axis limits
                    
                    #Draw time line on bottom ISR plot
                    imISRLine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( np.min(settings['ISR']['mill']['height lim']),np.max(settings['ISR']['mill']['height lim']))) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        bbox = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
                        plot_ratio = bbox.width/bbox.height; #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest               
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
    
                    #Draw time line on bottom ISR plot
                    imISRLine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( np.min(settings['ISR']['mill']['height lim']),np.max(settings['ISR']['mill']['height lim']))) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    #END IF
                    
                #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                movie_writer.grab_frame(); #get the frame and save it
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imISRLine.pop(0).remove();
    #            imOverlay.pop(0).remove();
    #            imMillstone.pop(0).remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 4 =================================================================
    elif(settings['movie']['movie type'] == 4): #stationary data points + AMPERE data on same plot with time average to AMPERE data (every 10 min)
        
        #---------------Prep for Averaging & Frame Multiplying-----------------
        gif_TimeAvg = np.int64(np.round(np.mean(np.diff(data['AMPERE']['time unique']))*24*60)); #min, force AVG to match the AMPERE time step (only for reporting at this point)
        
        timeUnique_TimeAvg = data['AMPERE']['time unique']; #hr, set up the time steps to be exactly what the AMPERE data time steps are
        if( timeUnique_TimeAvg[-1] < data['TEC']['time unique'][-1] ):
            timeUnique_TimeAvg = np.hstack( (timeUnique_TimeAvg,data['TEC']['time unique'][-1]) ); #hr, force last time range to be added if it didn't fit in well
        #END IF
        #makes coding easier
        
        #Stuff for repeating a frame - needed if want a constant 30 FPS
        #This needs to go in settings['movie']['movie type'] == 2 area
        gif_DesiredTotalTime = (timeUnique_TimeAvg.size-1)/gif_DesiredFPS; #s, total time from desired FPS
        gif_DesiredReqFrameNum = gif_DesiredTotalTime*30; #frames, number of frames needed to go for the time wanted
        gif_RepeatFrame = gif_DesiredReqFrameNum/(timeUnique_TimeAvg.size-1); #number of frames to repeat every time (exact w/ decimal)
        if( gif_RepeatFrame > 0.5 ):
            gif_RepeatFrame = np.int64(np.round(gif_RepeatFrame)); #round it
        else: #if less than 0.5 will round to 0 - so checking to see if there's a ton of frames and need to just boost it to 60 to keep run time in check
            if( ( (timeUnique_TimeAvg.size-1)/30 > settings['movie']['desired max run time'] ) & ( settings['movie']['disable FPS shift'] == 0 ) ): #try to keep run time down
                gif_vidFrameRate = 60; #set the FPS to 60
            #END IF
            gif_RepeatFrame = 1; #I guess to make sure it works
        #EMD IF
        
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
        
            #-------------------------Start Making Pictures------------------------
            for i in range(0,timeUnique_TimeAvg.size-1):
                
                #------------Corral the AMPERE data to the right place-------------
                k = np.where(data['AMPERE']['time'] == timeUnique_TimeAvg[i])[0]; #get where the time point is
                AMPERE_data_portion = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                AMPERE_lat_portion = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                AMPERE_long_portion = data['AMPERE']['long'][k]; #degc, corresponding long values
                
                #----------------Corral the data to the right place----------------
                k = np.where( ( data['TEC']['time'] >= timeUnique_TimeAvg[i]) & ( data['TEC']['time'] < timeUnique_TimeAvg[i+1] ) )[0]; #get where time points fall within the averaging boundaries
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot | AMPERE Time Avg\'d to Every '+str(gif_TimeAvg)+' Min | Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr'][0]); #create mecha title
                
                ax.set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                    
                #-------------------Starting the Plotting--------------------------
                if( i == 0 ): #first run preps axes, color bars, etc.
                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax.pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    #Do the AMPERE plotting
                    AMPERE_colorMap = ListedColormap( np.hstack(( np.array( ( (np.linspace(1,0.492063492063492,128)),(np.linspace(1,0.507936507936508,128)),(np.linspace(1,1,128)) ) ) , np.array( ( (np.linspace(0.492063492063492,1,128)) , (np.linspace(0.507936507936508,0,128)) , (np.linspace(1,1,128)) ) ) )).T ); #white to purpleblue to pink (based off of 'cool')
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax.scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    cbar2 = movie_dict['fig'].colorbar(imAMP, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    #                cbar2.ax.set_yticklabels(np.linspace(np.min(settings['AMPERE']['plot lim']),np.max(settings['AMPERE']['plot lim']),5)); #create useful tick marks
                    cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    cax2.yaxis.set_ticks_position('left'); #move it to the left
                    cax2.yaxis.set_label_position('left'); #move it to the left
                    cbar2.set_label('AMPERE Joule Heating '+r'$(erg/cm^{2}·sec)$'); #tabel the colorbar
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar2.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
                                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        bbox = ax.get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
                        plot_ratio = bbox.width/bbox.height; #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    #END IF
                    
                    #Now drawing point of interest
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax.plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax.pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
    
                    #Do the AMPERE plotting
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax.scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    #END IF
                #END IF
            
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();

                for j in range(0,gif_RepeatFrame): #runs this as many times as needed to get 30 FPS or whatever
                    movie_writer.grab_frame(); #get the frame and save it
                #END FOR j
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imAMP.remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                #END IF
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 5 =================================================================
    elif(settings['movie']['movie type'] == 5): #stationary data points + AMPERE data on same plot (no time average for TEC)
        
        # time alignment
        dateNow_dayNum = subfun_time_to_dateRange(data['TEC']['time unique'], dates_zeroHr = dates['date range zero hr'], FLG_timesWRTzeroHr = False, options = 0)[0]; #get the date the starting time index is on
        dateNow = subfun_dayNum_to_date(dateNow_dayNum)[0]; #convert to date from dayNum
        
        dateNow_secInDay = data['TEC']['time unique'] - np.int64(dateNow_dayNum[1])*86400; #get seconds in the day
        dateNow_hour = dateNow_secInDay//3600; #hours that fit in
        dateNow_min = (dateNow_secInDay - dateNow_hour*3600)//60; #minutes that fit in
        dateNow_sec = dateNow_secInDay - dateNow_min*60 - dateNow_hour*3600; #left over seconds
                
        time4mag = datetime(dateNow[0], dateNow[1], dateNow[2],\
            hour = dateNow_hour, minute = dateNow_min, second = dateNow_sec); #date time object for aacgmv2    
            
        # create the frame
        movie_dict = GRITI_movieMaker_subfun_figMaker(data, settings, time4mag, movie_figOffsets=movie_figOffsets);
        fig = movie_dict['fig']; #unpack
        ax = movie_dict['ax'];
        cax = movie_dict['cax'];
        cax2 = movie_dict['cax2'];
    
        #pre-calc what can be
        AMPERE_lat_delta = data['AMPERE']['data info']['lat delta']; #it's provided
        AMPERE_long_delta = data['AMPERE']['data info']['long delta']; #it's provided
        gif_Grid_Lat_AMPERE = np.arange(np.min(settings['map']['lat range']),np.max(settings['map']['lat range'])+AMPERE_lat_delta,AMPERE_lat_delta); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
        gif_Grid_Long_AMPERE = np.arange(np.min(settings['map']['long range']),np.max(settings['map']['long range'])+AMPERE_long_delta,AMPERE_long_delta); #degc, create long points specilized for AMPERE
        time_cutout_range_delay_AMPERE_sec = settings['AMPERE']['delay wrt TEC']*3600;
        if( np.isclose(time_cutout_range_delay_AMPERE_sec, np.int64(time_cutout_range_delay_AMPERE_sec)) ):
            time_cutout_range_delay_AMPERE_sec = np.int64(time_cutout_range_delay_AMPERE_sec); #integer it
        #END IF
        pltHelprX_TEC, pltHelprY_TEC = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
        pltHelprX_AMPERE, pltHelprY_AMPERE = np.meshgrid( gif_Grid_Long_AMPERE, gif_Grid_Lat_AMPERE); #helps the pcolor work
    
        # AMPERE_dataRate = data['AMPERE']['data rate']; #sec, get the median data rate for AMPERE data (avoids outliers)
        if( 'stere' in settings['map']['projection name'] ):
            if( settings['movie']['spin'] != 0):
                #this is a constant sun since the plot is continually spun to match
                x = (1.05*0.5*np.cos(np.pi/2))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.sin(np.pi/2))+0.5; #geoMap coordinate 
                circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                hSunStatic = ax.add_artist(circleSun); #plot the sun
            #END IF
        #END IF
        
        gif_Grid = GRITI_movieMaker_subfun_dataGridder( data['TEC']['lat'], data['TEC']['long'], data['TEC']['dTEC'],settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
        #call a numba'd function that makes the movie quicker by crunching the numbers gooder
        
        gif_Grid_AMPERE = GRITI_movieMaker_subfun_dataGridder(data['AMPERE']['lat'],data['AMPERE']['long'],data['AMPERE'][settings['AMPERE']['data type']],gif_Grid_Lat_AMPERE,gif_Grid_Long_AMPERE,gif_Grid_Lat_AMPERE.size-1,gif_Grid_Long_AMPERE.size-1,AMPERE_lat_delta,AMPERE_long_delta,settings['movie']['data reject ratio'],101,settings['movie']['data reject ratio max']).T; #101 disables the data rejection stuff b/c AMPERE doesn't need it
        
        #----------------------------Tack on Title-------------------------
        curr_time_mod = np.mod(data['TEC']['time unique'],86400);
        curr_time_date = subfun_dayNum_to_date((dates['date range zero hr dayNum'][0],data['TEC']['time unique']//86400))[0]; #if year changes needs overhaul everywhere
        string_title = '{0:.2f}'.format(np.round((data['TEC']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
            ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0])+\
            ' | ('+str(curr_time_date[0])+'/'+str(curr_time_date[1])+'/'+str(curr_time_date[2])+', '+str(curr_time_mod//3600).zfill(2)+\
            ':'+str((curr_time_mod-curr_time_mod//3600*3600)//60).zfill(2)+':'+str(curr_time_mod-(curr_time_mod-curr_time_mod//3600*3600)//60*60-curr_time_mod//3600*3600).zfill(2)+')'; #create mecha title
        if( settings['movie']['use time delays'] == 1 ):
            string_title += ' [Time Shifted by '+settings['AMPERE']['delay wrt TEC']+' hrs]'; 
        #END IF
        
        ax.set_title(string_title,fontproperties=settings['plot']['font title FM'],y=movie_title_yOffset); #set the title, properties always needed
            
        #-------------------Starting the Plotting--------------------------
                           
        #Do the TEC plotting
        imTEC = ax.pcolormesh(pltHelprX_TEC, pltHelprY_TEC,  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=8, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid
        cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
        cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
        cbar.set_label(settings['TEC']['name']+settings['TEC']['units']); #tabel the colorbar
        for tick in cbar.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
            tick.label2.set_fontproperties(settings['plot']['font axis tick FM']); #yee
        #END FOR tick
        cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));

        #Do the AMPERE plotting
        if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
            imAMP = ax.pcolormesh(pltHelprX_AMPERE, pltHelprY_AMPERE, gif_Grid_AMPERE, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), cmap=settings['AMPERE']['colormap'], zorder=6, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid
        else:
            imAMP = ax.pcolormesh(pltHelprX_AMPERE, pltHelprY_AMPERE, gif_Grid_AMPERE, cmap=settings['AMPERE']['colormap'],zorder=6, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid
        #END IF
        cbar2 = movie_dict['fig'].colorbar(imAMP, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
#                cbar2.ax.set_yticklabels(np.linspace(np.min(settings['AMPERE']['plot lim']),np.max(settings['AMPERE']['plot lim']),5)); #create useful tick marks
        cax2.yaxis.set_ticks_position('left'); #move it to the left
        cax2.yaxis.set_label_position('left'); #move it to the left        
        cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
        cbar2.set_label(settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]); #tabel the colorbar
        # cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
        # for tick in cbar2.ax.yaxis.get_major_ticks(): #if only I could apply the font manager directly
        #     tick.label.set_fontproperties(settings['plot']['font axis tick FM']); #yee
        #     tick.label2.set_fontproperties(settings['plot']['font axis tick FM']); #yee
        # #END FOR tick
        # cbar2.ax.set_yticklabels(cbar2.ax.get_yticks(), settings['plot']['font axis tick FM']); #the others won't stick -maybe this?
        for label in cbar2.ax.get_yticklabels(): #if only I could apply the font manager directly
            label.set_fontproperties(settings['plot']['font axis tick FM']); #yee <- this one worked for cbar2 only, idk why
        #END FOR label
        cbar2.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
        
        if( settings['movie']['day nite line'] == 1 ):
            print('not supported sorry bring it back if you want it ')
#             movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
#             #Plot the sunrise/sunset terminators
#             dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
#             hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
#             #constants will use all the time - only for plotting of day/nite line so minor importance
#             #this stuff makes the text angle plotted mostly correct most of the time
#             bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
#             bboxAx0 = ax.get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
# #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
#             plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
#             dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
#             dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
#             dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
#             if( hLgd_FLG_day > 0 ): #only do work if it is there
#                 #calc all that day/nite stuff in one function to keep it from getting cluttered
#                 dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
#                 dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
#                 imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
#                 if( settings['movie']['day nite text'] == 1 ):
#                     dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
#                     textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
#                 #END IF
#             #END IF
            
#             dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
#             hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
#             if( hLgd_FLG_nite > 0 ): #only do work if it is there
#                 #calc all that day/nite stuff in one function to keep it from getting cluttered
#                 dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
#                 dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
#                 imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
#                 if( settings['movie']['day nite text'] == 1 ):
#                     dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
#                     textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
#                 #EMD IF
#             #END IF
        elif( settings['movie']['day nite line'] == 2 ):
            if( settings['map']['coord type'] == 'geo' ):
                shifted = data['TEC']['time unique']; #no shift
                latLongShift = None;
            else:
                #this ain't perfect, but it's plenty good for general where not sun
                # shifter = (np.mod(data['movie']['sun loc']['long'],360) - np.mod(data['movie']['sun loc']['long geo'],360))*240 #sec, effective time to shift by so Sun is in right spot; 240 = 24*3600/360
                # shifted = data['TEC']['time unique'] - shifter; #shift the time so Sun is in the right spot
                shifted = data['TEC']['time unique']; #no shift
                latLongShift = {'lat':data['movie']['sun loc']['lat'],'long':data['movie']['sun loc']['long'],'rel':False};
            #END IF
            dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(shifted/86400)),ndmin=2))[0]; #get the current yr/month/day
            dayNite_currentHr = np.int64((shifted/86400-np.int64(shifted/86400))*24); #hr, get the current hour
            dayNite_currentMin = np.int64( ((shifted/86400-np.int64(shifted/86400))*24 - dayNite_currentHr)*60); #min, get the current min
            dayNite_currentSec = np.int64( (((shifted/86400-np.int64(shifted/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
            dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
            # dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax, zorder=25); #basemap version
            dayNite_shade = ax.add_feature(nightshader(dayNite_currentTime, latLongShift=latLongShift, alpha=0.25),zorder=25); #nighttime shading, only relevant for geographic
        #END IF
        
        if( 'stere' in settings['map']['projection name'] ):
            if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc']['long']))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc']['long']+np.pi))+0.5; #geoMap coordinate 
                x = (1.05*0.5*np.cos(data['movie']['sun loc']['long rad']))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.sin(data['movie']['sun loc']['long rad']))+0.5; #geoMap coordinate 
                #hSun = ax[0].text(x,y,'SUN\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                hSun = ax.add_artist(circleSun); #plot the sun
            else:
                hSun = ax.set_theta_offset(settings['movie']['sun loc']['long rad']); #turn the whole plot so top is where the sun is
                #I'm not sure if I can get this working easily - not a lot of optons with cartopy.
            #END IF
        #END IF
        
        #Now drawing point of interest
        # Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
        if( (settings['map']['site coords'][0,0] <= np.max(settings['map']['lat range'])) & (settings['map']['site coords'][0,0] >= np.min(settings['map']['lat range'])) & (settings['map']['site coords'][0,1] <= np.max(settings['map']['long range'])) & (settings['map']['site coords'][0,1] >= np.min(settings['map']['long range'])) ):
            imMillstone = ax.plot(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], transform=cartopy.crs.PlateCarree(), zorder=50); #plot this, 50 always on top
        else:
            if( 'stere' in settings['map']['projection name'] ):
                diameter_offset = 0.09;
                x = ((1+diameter_offset)*0.5*np.cos((settings['map']['site coords'][0,1][0]-90)*np.pi/180))+0.5; #geoMap coordinate 
                y = ((1+diameter_offset)*0.5*np.sin((settings['map']['site coords'][0,1][0]-90)*np.pi/180))+0.5; #geoMap coordinate
                imMillstone = ax.arrow(x, y, (diameter_offset/2)*np.cos((settings['map']['site coords'][0,1][0]+90)*np.pi/180), (diameter_offset/2)*np.sin((settings['map']['site coords'][0,1][0]+90)*np.pi/180), width=0.007, head_width=0.007*3.5, head_length=0.007*3.5*1.00, length_includes_head=True, color=settings['map']['site marker color'], clip_on=False, transform=ax.transAxes); #plot this, 50 always on top
            #END IF
            #if needed, code something for rectangular plots but idk what it should be exactly rn
        #END IF
        
        #-----------------------Create Movie/GIF---------------------------
        #Makes the gif now
        plt.draw();
        
        # #-------------------Clean up for re-use----------------------------
        # #if forget one (like hOverlay) slows it way down after many plots
        # imTEC.remove();
        # imAMP.remove();
        # if( settings['movie']['day nite line'] == 1 ):
        #     if(hLgd_FLG_day > 0): #only delete if it is there
        #         imDayNite_day.pop(0).remove();
        #         if( settings['movie']['day nite text'] == 1 ):
        #             textDayNite_day.remove();
        #         #END IF
        #     #END IF
        #     if(hLgd_FLG_nite > 0): #only delete if it is there
        #         imDayNite_nite.pop(0).remove();
        #         if( settings['movie']['day nite text'] == 1 ):
        #             textDayNite_nite.remove();
        #         #END IF
        #     #END IF
        # elif( settings['movie']['day nite line'] == 2): #only delete if it is there
        #     # for middleman in dayNite_shade.collections: #this needs an extra helping hand
        #     #     middleman.remove();
        #     dayNite_shade.remove(); #ditch it
        # #END IF
        # if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
        #     hSun.remove();
        # #END IF

    #================================================================= MOVIE TYPE 6 =================================================================
    elif(settings['movie']['movie type'] == 6): #stationary data points + AMPERE data on same plot (no time average for TEC) + TEC keogram with line for current time
        #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
        if( np.isscalar(settings['TEC']['plot lim']) == 1 ):
            gif_TEC_plotLimValu = np.array( (-settings['TEC']['plot lim'],settings['TEC']['plot lim']) ); #make it a vector
        else:
            gif_TEC_plotLimValu = settings['TEC']['plot lim']; #keep it the same
        #END IF
        
        #calc the TEC keogram first, only need to calc it once
        movie_TEC_keo1, settings_keo1 = GRITI_keo_keogrammer(data['TEC']['dTEC'] ,data['TEC']['time'], data['TEC']['lat'], data['TEC']['long'], data['TEC']['time unique'], data['TEC']['time unique'], dates,
                settings_keo1, settings['paths'], settings_keo1_map, settings['plot'],
                FLG_disablePlot=2,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True);
        
        AMPERE_dataRate = np.median(np.diff(data['AMPERE']['time unique'])); #days, get the median data rate for AMPERE data (avoids outliers)
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
                        
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                
                #------------Corral the AMPERE data to the right place-------------
                k = np.where( (data['AMPERE']['time'] <= data['TEC']['time unique'][i]) & (data['AMPERE']['time'] >= (data['TEC']['time unique'][i]-AMPERE_dataRate)) )[0]; #get where the time point is, make sure it is within the data rate window
                AMPERE_data_portion = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                AMPERE_lat_portion = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                AMPERE_long_portion = data['AMPERE']['long'][k]; #degc, corresponding long values
                
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot | Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                    
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label("delta-vTEC [TECU]",y=0.87); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    #Do the AMPERE plotting
                    AMPERE_colorMap = ListedColormap( np.hstack(( np.array( ( (np.linspace(1,0.492063492063492,128)),(np.linspace(1,0.507936507936508,128)),(np.linspace(1,1,128)) ) ) , np.array( ( (np.linspace(0.492063492063492,1,128)) , (np.linspace(0.507936507936508,0,128)) , (np.linspace(1,1,128)) ) ) )).T ); #white to purpleblue to pink (based off of 'cool')
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    cbar2 = movie_dict['fig'].colorbar(imAMP, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    #                cbar2.ax.set_yticklabels(np.linspace(np.min(settings['AMPERE']['plot lim']),np.max(settings['AMPERE']['plot lim']),5)); #create useful tick marks
                    cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    cax2.yaxis.set_ticks_position('left'); #move it to the left
                    cax2.yaxis.set_label_position('left'); #move it to the left
                    cbar2.set_label('AMPERE Joule Heating [$erg/cm^{2}·sec$]'); #tabel the colorbar
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar2.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i])),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    #Now drawing point of interest
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax[0].text(x,y,'SUN\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #------Now plot the keogram--------
                    pltHelprX, pltHelprY = np.meshgrid( (data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600, \
                                settings_keo1['keo plot latlong chunks']);
                    im = ax[1].pcolormesh(pltHelprX, pltHelprY,  movie_TEC_keo1.T ,vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu),cmap='jet'); # pseudocolor plot "stretched" to the grid
                    cbarKeo = movie_dict['fig'].colorbar(im, cax=caxKeo, orientation='vertical'); #create a colorbar using the prev. defined cax
                    caxKeo.yaxis.set_ticks(np.linspace(np.min(gif_TEC_plotLimValu),np.max(gif_TEC_plotLimValu),5)); #create useful tick marks
                    caxKeo.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbarKeo.set_label('delta-vTEC [TECU]'); #tabel the colorbar
                    cbarKeo.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbarKeo.mappable.set_clim(vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu));
                    caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    #    string_title = 'TEC Averaged on Angle of '+str(np.round(keo,2))+' deg and Width of '+ \
                    #        str(np.round(settings['TEC']['keo']['keo width'],2))+' arcdeg, Avg Step # = '+str(keo_N)+ \
                    #        ' arcdeg, Line Shows '+settings['TEC']['keo']['keo plot latlong name']+' of Millstone Hill Zenith Beam'; #create mecha title
                    string_title = 'delta-vTEC Keogram of '+settings['movie']['side plots']['keo 1 name']; #create mecha title
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[1].set_xlabel('Time in UT - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[1].set_ylabel(settings_keo1['keo plot latlong name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            8); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    
                    keo_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(settings_keo1['keo plot latlong chunks'])) - np.floor(np.min(settings_keo1['keo plot latlong chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
                    if( keo_Range_Chunks_Long_Plot_autoTick > 25 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 10 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 5 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 2 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 1 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                        keo_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
                    else:
                        if(settings_keo1['keo plot latlong name'] == 'Latitude'): #if Y axis is latitude, use latitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 1 lat range']) - np.min(settings['movie']['side plots']['keo 1 lat range']))/13; #just goes for it if it's a super tiny range
                        elif(settings_keo1['keo plot latlong name'] == 'Longitude'): #if Y axis is longitude, use longitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 1 long range']) - np.min(settings['movie']['side plots']['keo 1 long range']))/13; #just goes for it if it's a super tiny range
                        #END IF
                    #END IF
                    yAxisTicks = np.round(np.arange( np.floor(np.min(settings_keo1['keo plot latlong chunks'])),np.ceil(np.max(settings_keo1['keo plot latlong chunks']))+keo_Range_Chunks_Long_Plot_autoTick,keo_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
                    ax[1].set_yticks(yAxisTicks); #set x axis ticks
                    ax[1].set_ylim( np.round(ax[1].get_ylim()) ); #avoid weird rounding errors
                    
                    #Now drawing line of interest
                    if( settings_keo1['keo plot latlong name'] == 'Longitude' ): #if true, longitude
                        if( (np.min(settings['movie']['side plots']['keo 1 long range']) <= settings['map']['site coords'][0,1]) & (np.max(settings['movie']['side plots']['keo 1 long range']) >= settings['map']['site coords'][0,1]) ): #only plot if it's in the long range specified
                            ax[1].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,1],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    else: #else latitude
                        if( (np.min(settings['movie']['side plots']['keo 1 lat range']) <= settings['map']['site coords'][0,0]) & (np.max(settings['movie']['side plots']['keo 1 lat range']) >= settings['map']['site coords'][0,0]) ): #only plot if it's in the lat range specified
                            ax[1].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,0],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    #END IF
                    
                    #-----plot vertical line showing the time on the keogram-----
                    hVline = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
    
                    #Do the AMPERE plotting
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i])),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #----plot vertical line showing the time on the keogram-----
                    hVline = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                #END IF
                
                #Dynamically show and hide the day/nite labels
    #            if( (hLgd_FLG_day > 0) & (hLgd_FLG_nite > 0) ):
    #                hLgd = ax.legend( (imDayNite_day , imDayNite_nite) , ('Sunrise','Sunset') , loc=4 , fontsize = settings['plot']['font axis tick'] ); #make that legend
    #            elif( (hLgd_FLG_day > 0) & (hLgd_FLG_nite == 0) ):
    #                hLgd = ax.legend( (imDayNite_day) , ('Sunrise') , loc=4 , fontsize = settings['plot']['font axis tick'] , handlelength=.75 ); #make that legend
    #            elif( (hLgd_FLG_day == 0) & (hLgd_FLG_nite > 0) ):
    #                hLgd = ax.legend( (imDayNite_nite,) , ('Sunset',) , loc=4 , fontsize = settings['plot']['font axis tick'] ); #make that legend
    #            #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                movie_writer.grab_frame(); #get the frame and save it
                    
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imAMP.remove();
                hVline.remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                elif( settings['movie']['day nite line'] == 2): #only delete if it is there
                    for middleman in dayNite_shade.collections: #this needs an extra helping hand
                        middleman.remove();
                #END IF
                if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
                    hSun.remove();
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 7 =================================================================
    elif(settings['movie']['movie type'] == 7): #stationary data points + AMPERE data on same plot (no time average for TEC) + TEC keogram with line for current time
        #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
        if( np.isscalar(settings['TEC']['plot lim']) == 1 ):
            gif_TEC_plotLimValu = np.array( (-settings['TEC']['plot lim'],settings['TEC']['plot lim']) ); #make it a vector
        else:
            gif_TEC_plotLimValu = settings['TEC']['plot lim']; #keep it the same
        #END IF
                
        #calc the TEC keogram first, only need to calc it once
        movie_TEC_keo1, settings_keo1 = GRITI_keo_keogrammer(data['TEC']['dTEC'] ,data['TEC']['time'], data['TEC']['lat'], data['TEC']['long'], data['TEC']['time unique'], data['TEC']['time unique'], dates,
                settings_keo1, settings['paths'], settings_keo1_map, settings['plot'],
                FLG_disablePlot=2,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True);
        
        movie_TEC_keo2, settings_keo2 = GRITI_keo_keogrammer(data['TEC']['dTEC'] ,data['TEC']['time'], data['TEC']['lat'], data['TEC']['long'], data['TEC']['time unique'], data['TEC']['time unique'], dates,
                settings_keo2, settings['paths'], settings_keo2_map, settings['plot'],
                FLG_disablePlot=2,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True);
        
        AMPERE_dataRate = np.median(np.diff(data['AMPERE']['time unique'])); #days, get the median data rate for AMPERE data (avoids outliers)
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
                        
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                
                #------------Corral the AMPERE data to the right place-------------
                if( settings['movie']['use time delays'] == 1 ):
                    k = np.where( ((data['AMPERE']['time']+settings['AMPERE']['delay wrt TEC']/24) <= data['TEC']['time unique'][i]) & ((data['AMPERE']['time']+settings['AMPERE']['delay wrt TEC']/24) >= (data['TEC']['time unique'][i]-AMPERE_dataRate)) )[0]; #get where the time point is, make sure it is within the data rate window
                else:
                    k = np.where( (data['AMPERE']['time'] <= data['TEC']['time unique'][i]) & (data['AMPERE']['time'] >= (data['TEC']['time unique'][i]-AMPERE_dataRate)) )[0]; #get where the time point is, make sure it is within the data rate window
                #END IF
                AMPERE_data_portion = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                AMPERE_lat_portion = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                AMPERE_long_portion = data['AMPERE']['long'][k]; #degc, corresponding long values
                
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot | Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title, properties always needed
                    
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    # cbar.set_label("delta-vTEC [TECU]",y=0.5); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    #Do the AMPERE plotting
                    AMPERE_colorMap = ListedColormap( np.hstack(( np.array( ( (np.linspace(1,0.492063492063492,128)),(np.linspace(1,0.507936507936508,128)),(np.linspace(1,1,128)) ) ) , np.array( ( (np.linspace(0.492063492063492,1,128)) , (np.linspace(0.507936507936508,0,128)) , (np.linspace(1,1,128)) ) ) )).T ); #white to purpleblue to pink (based off of 'cool')
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    cbar2 = movie_dict['fig'].colorbar(imAMP, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    #                cbar2.ax.set_yticklabels(np.linspace(np.min(settings['AMPERE']['plot lim']),np.max(settings['AMPERE']['plot lim']),5)); #create useful tick marks
                    if( settings['movie']['use time delays'] == 1 ):
                        cbar2.set_label('AMPERE '+settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]+' (Delayed by '+str(settings['AMPERE']['delay wrt TEC'])+' hrs)'); #tabel the colorbar
                    else:
                        cbar2.set_label('AMPERE '+settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]); #tabel the colorbar
                    #END IF
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar2.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
                    if( np.all(np.mod(cbar2.get_ticks(),1) == 0) ):
                        cax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f')); #force a rounded format
                    else:
                        cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    #END IF
                    cax2.yaxis.set_ticks_position('left'); #move it to the left
                    cax2.yaxis.set_label_position('left'); #move it to the left
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i])),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    #Now drawing point of interest
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax[0].text(x,y,'SUN\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #------Now plot the keogram--------
                    pltHelprX, pltHelprY = np.meshgrid( (data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600, \
                                settings_keo1['keo plot latlong chunks']);
                    im = ax[1].pcolormesh(pltHelprX, pltHelprY,  movie_TEC_keo1.T ,vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu),cmap='jet'); # pseudocolor plot "stretched" to the grid
                    cbarKeo = movie_dict['fig'].colorbar(im, cax=caxKeo, orientation='vertical'); #create a colorbar using the prev. defined cax
                    caxKeo.yaxis.set_ticks(np.linspace(np.min(gif_TEC_plotLimValu),np.max(gif_TEC_plotLimValu),5)); #create useful tick marks
                    caxKeo.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbarKeo.set_label('delta-vTEC [TECU]'); #tabel the colorbar
                    cbarKeo.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbarKeo.mappable.set_clim(vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu));
                    caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    #    string_title = 'TEC Averaged on Angle of '+str(np.round(keo,2))+' deg and Width of '+ \
                    #        str(np.round(settings['TEC']['keo']['keo width'],2))+' arcdeg, Avg Step # = '+str(keo_N)+ \
                    #        ' arcdeg, Line Shows '+settings['TEC']['keo']['keo plot latlong name']+' of Millstone Hill Zenith Beam'; #create mecha title
                    string_title = 'delta-vTEC Keogram of '+settings['movie']['side plots']['keo 1 name']; #create mecha title
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[1].set_xlabel('Time in UT - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[1].set_ylabel(settings_keo1['keo plot latlong name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            8); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    
                    keo_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(settings_keo1['keo plot latlong chunks'])) - np.floor(np.min(settings_keo1['keo plot latlong chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
                    if( keo_Range_Chunks_Long_Plot_autoTick > 25 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 10 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 5 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 2 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 1 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                        keo_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
                    else:
                        if(settings_keo1['keo plot latlong name'] == 'Latitude'): #if Y axis is latitude, use latitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 1 lat range']) - np.min(settings['movie']['side plots']['keo 1 lat range']))/13; #just goes for it if it's a super tiny range
                        elif(settings_keo1['keo plot latlong name'] == 'Longitude'): #if Y axis is longitude, use longitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 1 long range']) - np.min(settings['movie']['side plots']['keo 1 long range']))/13; #just goes for it if it's a super tiny range
                        #END IF
                    #END IF
                    yAxisTicks = np.round(np.arange( np.floor(np.min(settings_keo1['keo plot latlong chunks'])),np.ceil(np.max(settings_keo1['keo plot latlong chunks']))+keo_Range_Chunks_Long_Plot_autoTick,keo_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
                    ax[1].set_yticks(yAxisTicks); #set x axis ticks
                    ax[1].set_ylim( np.round(ax[1].get_ylim()) ); #avoid weird rounding errors
                    
                    #Now drawing line of interest
                    if( settings_keo1['keo plot latlong name'] == 'Longitude' ): #if true, longitude
                        if( (np.min(settings['movie']['side plots']['keo 1 long range']) <= settings['map']['site coords'][0,1]) & (np.max(settings['movie']['side plots']['keo 1 long range']) >= settings['map']['site coords'][0,1]) ): #only plot if it's in the long range specified
                            ax[1].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,1],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    else: #else latitude
                        if( (np.min(settings['movie']['side plots']['keo 1 lat range']) <= settings['map']['site coords'][0,0]) & (np.max(settings['movie']['side plots']['keo 1 lat range']) >= settings['map']['site coords'][0,0]) ): #only plot if it's in the lat range specified
                            ax[1].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,0],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    #END IF
                    
                    #-----plot vertical line showing the time on the keogram-----
                    hVline = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    #------Now plot the keogram2--------
                    pltHelprX, pltHelprY = np.meshgrid( (data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600, \
                                settings_keo2['keo plot latlong chunks']);
                    imKeo2 = ax[2].pcolormesh(pltHelprX, pltHelprY,  movie_TEC_keo2.T ,vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu),cmap='jet'); # pseudocolor plot "stretched" to the grid
                    cbarKeo2 = movie_dict['fig'].colorbar(imKeo2, cax=caxKeo2, orientation='vertical'); #create a colorbar using the prev. defined cax
                    caxKeo2.yaxis.set_ticks(np.linspace(np.min(gif_TEC_plotLimValu),np.max(gif_TEC_plotLimValu),5)); #create useful tick marks
                    caxKeo2.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbarKeo2.set_label('delta-vTEC [TECU]'); #tabel the colorbar
                    cbarKeo2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbarKeo2.mappable.set_clim(vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu));
                    caxKeo2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    #    string_title = 'TEC Averaged on Angle of '+str(np.round(keo,2))+' deg and Width of '+ \
                    #        str(np.round(settings['TEC']['keo']['keo width'],2))+' arcdeg, Avg Step # = '+str(keo_N)+ \
                    #        ' arcdeg, Line Shows '+settings['TEC']['keo']['keo plot latlong name']+' of Millstone Hill Zenith Beam'; #create mecha title
                    string_title = 'delta-vTEC Keogram of '+settings['movie']['side plots']['keo 2 name']; #create mecha title
                    ax[2].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[2].set_xlabel('Time in UT - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[2].set_ylabel(settings_keo2['keo plot latlong name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            8); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[2].set_xticks(xAxisTicks); #set x axis ticks
                    
                    keo_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(settings_keo2['keo plot latlong chunks'])) - np.floor(np.min(settings_keo2['keo plot latlong chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
                    if( keo_Range_Chunks_Long_Plot_autoTick > 25 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 10 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 5 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 2 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 1 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                        keo_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
                    else:
                        if(settings_keo2['keo plot latlong name'] == 'Latitude'): #if Y axis is latitude, use latitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 2 lat range']) - np.min(settings['movie']['side plots']['keo 2 lat range']))/13; #just goes for it if it's a super tiny range
                        elif(settings_keo2['keo plot latlong name'] == 'Longitude'): #if Y axis is longitude, use longitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 2 long range']) - np.min(settings['movie']['side plots']['keo 2 long range']))/13; #just goes for it if it's a super tiny range
                        #END IF
                    #END IF
                    yAxisTicks = np.round(np.arange( np.floor(np.min(settings_keo2['keo plot latlong chunks'])),np.ceil(np.max(settings_keo2['keo plot latlong chunks']))+keo_Range_Chunks_Long_Plot_autoTick,keo_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
                    ax[2].set_yticks(yAxisTicks); #set x axis ticks
                    ax[2].set_ylim( np.round(ax[2].get_ylim()) ); #avoid weird rounding errors
                    
                    #Now drawing line of interest
                    if( settings_keo2['keo plot latlong name'] == 'Longitude' ): #if true, longitude
                        if( (np.min(settings['movie']['side plots']['keo 2 long range']) <= settings['map']['site coords'][0,1]) & (np.max(settings['movie']['side plots']['keo 2 long range']) >= settings['map']['site coords'][0,1]) ): #only plot if it's in the long range specified
                            ax[2].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,1],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    else: #else latitude
                        if( (np.min(settings['movie']['side plots']['keo 2 lat range']) <= settings['map']['site coords'][0,0]) & (np.max(settings['movie']['side plots']['keo 2 lat range']) >= settings['map']['site coords'][0,0]) ): #only plot if it's in the lat range specified
                            ax[2].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,0],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    #END IF
                    
                    #-----plot vertical line showing the time on the keogram2-----
                    hVline2 = ax[2].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
    
                    #Do the AMPERE plotting
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i])),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]-np.int64(data['TEC']['time unique'][i]))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #----plot vertical line showing the time on the keogram-----
                    hVline = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    #----plot vertical line showing the time on the keogram2-----
                    hVline2 = ax[2].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                #END IF
                
                #Dynamically show and hide the day/nite labels
    #            if( (hLgd_FLG_day > 0) & (hLgd_FLG_nite > 0) ):
    #                hLgd = ax.legend( (imDayNite_day , imDayNite_nite) , ('Sunrise','Sunset') , loc=4 , fontsize = settings['plot']['font axis tick'] ); #make that legend
    #            elif( (hLgd_FLG_day > 0) & (hLgd_FLG_nite == 0) ):
    #                hLgd = ax.legend( (imDayNite_day) , ('Sunrise') , loc=4 , fontsize = settings['plot']['font axis tick'] , handlelength=.75 ); #make that legend
    #            elif( (hLgd_FLG_day == 0) & (hLgd_FLG_nite > 0) ):
    #                hLgd = ax.legend( (imDayNite_nite,) , ('Sunset',) , loc=4 , fontsize = settings['plot']['font axis tick'] ); #make that legend
    #            #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                movie_writer.grab_frame(); #get the frame and save it
                    
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imAMP.remove();
                hVline.remove();
                hVline2.remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                elif( settings['movie']['day nite line'] == 2): #only delete if it is there
                    for middleman in dayNite_shade.collections: #this needs an extra helping hand
                        middleman.remove();
                #END IF
                if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
                    hSun.remove();
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 71 =================================================================
    elif(settings['movie']['movie type'] == 71): #stationary data points (no time average for TEC) + TEC keogram with line for current time
        #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
        if( np.isscalar(settings['TEC']['plot lim']) == 1 ):
            gif_TEC_plotLimValu = np.array( (-settings['TEC']['plot lim'],settings['TEC']['plot lim']) ); #make it a vector
        else:
            gif_TEC_plotLimValu = settings['TEC']['plot lim']; #keep it the same
        #END IF
        
        #calc the TEC keogram first, only need to calc it once
        movie_TEC_keo1, settings_keo1 = GRITI_keo_keogrammer(data['TEC']['dTEC'] ,data['TEC']['time'], data['TEC']['lat'], data['TEC']['long'], data['TEC']['time unique'], data['TEC']['time unique'], dates,
                settings_keo1, settings['paths'], settings_keo1_map, settings['plot'],
                FLG_disablePlot=2,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True);
        
        movie_TEC_keo2, settings_keo2 = GRITI_keo_keogrammer(data['TEC']['dTEC'] ,data['TEC']['time'], data['TEC']['lat'], data['TEC']['long'], data['TEC']['time unique'], data['TEC']['time unique'], dates,
                settings_keo2, settings['paths'], settings_keo2_map, settings['plot'],
                FLG_disablePlot=2,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True);
        
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
                        
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                                
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot | Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title, properties always needed
                    
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    # cbar.set_label("delta-vTEC [TECU]",y=0.5); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i]/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax[0].text(x,y,'SUN\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #Now drawing point of interest
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                                        
                    #------Now plot the keogram--------
                    pltHelprX, pltHelprY = np.meshgrid( (data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600, \
                                settings_keo1['keo plot latlong chunks']);
                    im = ax[1].pcolormesh(pltHelprX, pltHelprY,  movie_TEC_keo1.T ,vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu),cmap='jet'); # pseudocolor plot "stretched" to the grid
                    cbarKeo = movie_dict['fig'].colorbar(im, cax=caxKeo, orientation='vertical'); #create a colorbar using the prev. defined cax
                    caxKeo.yaxis.set_ticks(np.linspace(np.min(gif_TEC_plotLimValu),np.max(gif_TEC_plotLimValu),5)); #create useful tick marks
                    caxKeo.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbarKeo.set_label('delta-vTEC [TECU]'); #tabel the colorbar
                    cbarKeo.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbarKeo.mappable.set_clim(vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu));
                    caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    #    string_title = 'TEC Averaged on Angle of '+str(np.round(keo,2))+' deg and Width of '+ \
                    #        str(np.round(settings['TEC']['keo']['keo width'],2))+' arcdeg, Avg Step # = '+str(keo_N)+ \
                    #        ' arcdeg, Line Shows '+settings['TEC']['keo']['keo plot latlong name']+' of Millstone Hill Zenith Beam'; #create mecha title
                    string_title = 'delta-vTEC Keogram of '+settings['movie']['side plots']['keo 1 name']; #create mecha title
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[1].set_xlabel('Time in UT - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[1].set_ylabel(settings_keo1['keo plot latlong name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            8); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    
                    keo_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(settings_keo1['keo plot latlong chunks'])) - np.floor(np.min(settings_keo1['keo plot latlong chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
                    if( keo_Range_Chunks_Long_Plot_autoTick > 25 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 10 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 5 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 2 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 1 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                        keo_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
                    else:
                        if(settings_keo1['keo plot latlong name'] == 'Latitude'): #if Y axis is latitude, use latitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 1 lat range']) - np.min(settings['movie']['side plots']['keo 1 lat range']))/13; #just goes for it if it's a super tiny range
                        elif(settings_keo1['keo plot latlong name'] == 'Longitude'): #if Y axis is longitude, use longitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 1 long range']) - np.min(settings['movie']['side plots']['keo 1 long range']))/13; #just goes for it if it's a super tiny range
                        #END IF
                    #END IF
                    yAxisTicks = np.round(np.arange( np.floor(np.min(settings_keo1['keo plot latlong chunks'])),np.ceil(np.max(settings_keo1['keo plot latlong chunks']))+keo_Range_Chunks_Long_Plot_autoTick,keo_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
                    ax[1].set_yticks(yAxisTicks); #set x axis ticks
                    ax[1].set_ylim( np.round(ax[1].get_ylim()) ); #avoid weird rounding errors
                    
                    #Now drawing line of interest
                    if( settings_keo1['keo plot latlong name'] == 'Longitude' ): #if true, longitude
                        if( (np.min(settings['movie']['side plots']['keo 1 long range']) <= settings['map']['site coords'][0,1]) & (np.max(settings['movie']['side plots']['keo 1 long range']) >= settings['map']['site coords'][0,1]) ): #only plot if it's in the long range specified
                            ax[1].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,1],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    else: #else latitude
                        if( (np.min(settings['movie']['side plots']['keo 1 lat range']) <= settings['map']['site coords'][0,0]) & (np.max(settings['movie']['side plots']['keo 1 lat range']) >= settings['map']['site coords'][0,0]) ): #only plot if it's in the lat range specified
                            ax[1].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,0],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    #END IF
                    
                    #-----plot vertical line showing the time on the keogram-----
                    hVline = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    #------Now plot the keogram2--------
                    pltHelprX, pltHelprY = np.meshgrid( (data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600, \
                                settings_keo2['keo plot latlong chunks']);
                    imKeo2 = ax[2].pcolormesh(pltHelprX, pltHelprY,  movie_TEC_keo2.T ,vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu),cmap='jet'); # pseudocolor plot "stretched" to the grid
                    cbarKeo2 = movie_dict['fig'].colorbar(imKeo2, cax=caxKeo2, orientation='vertical'); #create a colorbar using the prev. defined cax
                    caxKeo2.yaxis.set_ticks(np.linspace(np.min(gif_TEC_plotLimValu),np.max(gif_TEC_plotLimValu),5)); #create useful tick marks
                    caxKeo2.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbarKeo2.set_label('delta-vTEC [TECU]'); #tabel the colorbar
                    cbarKeo2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbarKeo2.mappable.set_clim(vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu));
                    caxKeo2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    #    string_title = 'TEC Averaged on Angle of '+str(np.round(keo,2))+' deg and Width of '+ \
                    #        str(np.round(settings['TEC']['keo']['keo width'],2))+' arcdeg, Avg Step # = '+str(keo_N)+ \
                    #        ' arcdeg, Line Shows '+settings['TEC']['keo']['keo plot latlong name']+' of Millstone Hill Zenith Beam'; #create mecha title
                    string_title = 'delta-vTEC Keogram of '+settings['movie']['side plots']['keo 2 name']; #create mecha title
                    ax[2].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[2].set_xlabel('Time in UT - 0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0])+' [hr]',fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[2].set_ylabel(settings_keo2['keo plot latlong name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            8); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[2].set_xticks(xAxisTicks); #set x axis ticks
                    
                    keo_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(settings_keo2['keo plot latlong chunks'])) - np.floor(np.min(settings_keo2['keo plot latlong chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
                    if( keo_Range_Chunks_Long_Plot_autoTick > 25 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 10 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 5 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 2 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 1 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                        keo_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
                    else:
                        if(settings_keo2['keo plot latlong name'] == 'Latitude'): #if Y axis is latitude, use latitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 2 lat range']) - np.min(settings['movie']['side plots']['keo 2 lat range']))/13; #just goes for it if it's a super tiny range
                        elif(settings_keo2['keo plot latlong name'] == 'Longitude'): #if Y axis is longitude, use longitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 2 long range']) - np.min(settings['movie']['side plots']['keo 2 long range']))/13; #just goes for it if it's a super tiny range
                        #END IF
                    #END IF
                    yAxisTicks = np.round(np.arange( np.floor(np.min(settings_keo2['keo plot latlong chunks'])),np.ceil(np.max(settings_keo2['keo plot latlong chunks']))+keo_Range_Chunks_Long_Plot_autoTick,keo_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
                    ax[2].set_yticks(yAxisTicks); #set x axis ticks
                    ax[2].set_ylim( np.round(ax[2].get_ylim()) ); #avoid weird rounding errors
                    
                    #Now drawing line of interest
                    if( settings_keo2['keo plot latlong name'] == 'Longitude' ): #if true, longitude
                        if( (np.min(settings['movie']['side plots']['keo 2 long range']) <= settings['map']['site coords'][0,1]) & (np.max(settings['movie']['side plots']['keo 2 long range']) >= settings['map']['site coords'][0,1]) ): #only plot if it's in the long range specified
                            ax[2].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,1],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    else: #else latitude
                        if( (np.min(settings['movie']['side plots']['keo 2 lat range']) <= settings['map']['site coords'][0,0]) & (np.max(settings['movie']['side plots']['keo 2 lat range']) >= settings['map']['site coords'][0,0]) ): #only plot if it's in the lat range specified
                            ax[2].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,0],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    #END IF
                    
                    #-----plot vertical line showing the time on the keogram2-----
                    hVline2 = ax[2].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    figFitter(fig); #fit the fig fast one last time
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid

                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i]/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #----plot vertical line showing the time on the keogram-----
                    hVline = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    #----plot vertical line showing the time on the keogram2-----
                    hVline2 = ax[2].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                #END IF
                
                #Dynamically show and hide the day/nite labels
    #            if( (hLgd_FLG_day > 0) & (hLgd_FLG_nite > 0) ):
    #                hLgd = ax.legend( (imDayNite_day , imDayNite_nite) , ('Sunrise','Sunset') , loc=4 , fontsize = settings['plot']['font axis tick'] ); #make that legend
    #            elif( (hLgd_FLG_day > 0) & (hLgd_FLG_nite == 0) ):
    #                hLgd = ax.legend( (imDayNite_day) , ('Sunrise') , loc=4 , fontsize = settings['plot']['font axis tick'] , handlelength=.75 ); #make that legend
    #            elif( (hLgd_FLG_day == 0) & (hLgd_FLG_nite > 0) ):
    #                hLgd = ax.legend( (imDayNite_nite,) , ('Sunset',) , loc=4 , fontsize = settings['plot']['font axis tick'] ); #make that legend
    #            #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                movie_writer.grab_frame(); #get the frame and save it
                    
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                hVline.remove();
                hVline2.remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                elif( settings['movie']['day nite line'] == 2): #only delete if it is there
                    for middleman in dayNite_shade.collections: #this needs an extra helping hand
                        middleman.remove();
                #END IF
                if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
                    hSun.remove();
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 8 =================================================================
    elif(settings['movie']['movie type'] == 8): #stationary data points (through averaging) + OMNI overlay
            
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
            
            #-------------------------OMNI prep stuff-------------------------
            OMNI_timeUnique_hr = (data['OMNI']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
            
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                                
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
    
                    #Do the OMNI plotting
                    imOMNI = ax[1].plot( OMNI_timeUnique_hr, data['OMNI'][settings['OMNI']['data type']], linewidth=0.8 , zorder=1); #plot
            
                    if( np.mod(np.round(np.min(OMNI_timeUnique_hr)),2) == 0 ):
                        OMNI_time_hr_axis_min = np.round(np.min(OMNI_timeUnique_hr)); #is even, good to go
                    else:
                        OMNI_time_hr_axis_min = np.round(np.min(OMNI_timeUnique_hr))+1; #is odd, make even
                    #END IF
                    if( np.mod(np.round(np.max(OMNI_timeUnique_hr)),2) == 0 ):
                        OMNI_time_hr_axis_max = np.round(np.max(OMNI_timeUnique_hr)); #is even, good to go
                    else:
                        OMNI_time_hr_axis_max = np.round(np.max(OMNI_timeUnique_hr))-1; #is odd, make even
                    #END IF
                    
                    xAxisTicks = np.arange(OMNI_time_hr_axis_min,OMNI_time_hr_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    ax[1].set_xlim( OMNI_time_hr_axis_min , OMNI_time_hr_axis_max ); #set y axis limits
                    ax[1].set_ylabel(settings['OMNI']['labels'][settings['OMNI']['data type']]+settings['OMNI']['units'][settings['OMNI']['data type']],fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    ax[1].set_ylim( np.min(data['OMNI'][settings['OMNI']['data type']]) , np.max(data['OMNI'][settings['OMNI']['data type']]) ); #set y axis limits
                    ax[1].grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on
                    
                    #tack on 2nd title that doesn't change
                    OMNI_plot_scargle_label = settings['OMNI']['labels'][settings['OMNI']['data type']]+settings['OMNI']['units'][settings['OMNI']['data type']]; #get the label
                    OMNI_plot_scargle_labelNoUnits = OMNI_plot_scargle_label[0:OMNI_plot_scargle_label.find('(')-1]; #remove the (units)
                    string_title = 'OMNI '+OMNI_plot_scargle_labelNoUnits+' Index'; #create mecha title
                
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[1].set_xlabel(subfun_monthNum_to_word(dates['date range'][0,1])[0]+" "+str(dates['date range'][0,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][0,0])+" to "+subfun_monthNum_to_word(dates['date range'][1,1])[0]+" "+str(dates['date range'][1,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][1,0]),fontproperties=settings['plot']['font axis label FM']); #set the x axis label                
                    
                    #Draw time line on bottom OMNI plot
                    OMNI_data_min = np.min(data['OMNI'][settings['OMNI']['data type']]); #pre-calc the min
                    OMNI_data_max = np.max(data['OMNI'][settings['OMNI']['data type']]); #pre-calc the max
                    imOMNILine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400) , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400) )) , np.array(( OMNI_data_min,OMNI_data_max)) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest               
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
    
                    #Draw time line on bottom OMNI plot
                    imOMNILine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( OMNI_data_min,OMNI_data_max)) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    #END IF
                    
                #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                movie_writer.grab_frame(); #get the frame and save it
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imOMNILine.pop(0).remove();
    #            imOverlay.pop(0).remove();
    #            imMillstone.pop(0).remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 9 =================================================================
    elif(settings['movie']['movie type'] == 9): #stationary data points (through averaging) + AMPERE data on same plot (no time average for TEC) + OMNI overlay
        
        AMPERE_dataRate = np.median(np.diff(data['AMPERE']['time unique'])); #days, get the median data rate for AMPERE data (avoids outliers)
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
            
            #-------------------------OMNI prep stuff-------------------------
            OMNI_timeUnique_hr = (data['OMNI']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
            
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                
                #------------Corral the AMPERE data to the right place-------------
                k = np.where( (data['AMPERE']['time'] <= data['TEC']['time unique'][i]) & (data['AMPERE']['time'] >= (data['TEC']['time unique'][i]-AMPERE_dataRate)) )[0]; #get where the time point is, make sure it is within the data rate window
                AMPERE_data_portion = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                AMPERE_lat_portion = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                AMPERE_long_portion = data['AMPERE']['long'][k]; #degc, corresponding long values
                
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title, properties always needed
                                
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),5)); #create useful tick marks
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
                    
                    #Do the AMPERE plotting
                    AMPERE_colorMap = ListedColormap( np.hstack(( np.array( ( (np.linspace(1,0.492063492063492,128)),(np.linspace(1,0.507936507936508,128)),(np.linspace(1,1,128)) ) ) , np.array( ( (np.linspace(0.492063492063492,1,128)) , (np.linspace(0.507936507936508,0,128)) , (np.linspace(1,1,128)) ) ) )).T ); #white to purpleblue to pink (based off of 'cool')
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    cbar2 = movie_dict['fig'].colorbar(imAMP, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    #                cbar2.ax.set_yticklabels(np.linspace(np.min(settings['AMPERE']['plot lim']),np.max(settings['AMPERE']['plot lim']),5)); #create useful tick marks
                    cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    cax2.yaxis.set_ticks_position('left'); #move it to the left
                    cax2.yaxis.set_label_position('left'); #move it to the left
                    cbar2.set_label('AMPERE Joule Heating '+r'$(erg/cm^{2}·sec)$'); #tabel the colorbar
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar2.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
    
                    #Do the OMNI plotting
                    imOMNI = ax[1].plot( OMNI_timeUnique_hr, data['OMNI'][settings['OMNI']['data type']] , linewidth=0.8 , zorder=1); #plot
            
                    if( np.mod(np.round(np.min(OMNI_timeUnique_hr)),2) == 0 ):
                        OMNI_time_hr_axis_min = np.round(np.min(OMNI_timeUnique_hr)); #is even, good to go
                    else:
                        OMNI_time_hr_axis_min = np.round(np.min(OMNI_timeUnique_hr))+1; #is odd, make even
                    #END IF
                    if( np.mod(np.round(np.max(OMNI_timeUnique_hr)),2) == 0 ):
                        OMNI_time_hr_axis_max = np.round(np.max(OMNI_timeUnique_hr)); #is even, good to go
                    else:
                        OMNI_time_hr_axis_max = np.round(np.max(OMNI_timeUnique_hr))-1; #is odd, make even
                    #END IF
                    
                    xAxisTicks = np.arange(OMNI_time_hr_axis_min,OMNI_time_hr_axis_max+4,4); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    ax[1].set_xlim( OMNI_time_hr_axis_min , OMNI_time_hr_axis_max ); #set y axis limits
                    ax[1].set_ylabel(settings['OMNI']['labels'][settings['OMNI']['data type']]+settings['OMNI']['units'][settings['OMNI']['data type']],fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    ax[1].set_ylim( np.min(data['OMNI'][settings['OMNI']['data type']]) , np.max(data['OMNI'][settings['OMNI']['data type']]) ); #set y axis limits
                    ax[1].grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on
                    
                    #tack on 2nd title that doesn't change
                    OMNI_plot_scargle_label = settings['OMNI']['labels'][settings['OMNI']['data type']]+settings['OMNI']['units'][settings['OMNI']['data type']]; #get the label
                    OMNI_plot_scargle_labelNoUnits = OMNI_plot_scargle_label[0:OMNI_plot_scargle_label.find('(')-1]; #remove the (units)
                    string_title = 'OMNI '+OMNI_plot_scargle_labelNoUnits+' Index'; #create mecha title
                
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[1].set_xlabel(subfun_monthNum_to_word(dates['date range'][0,1])[0]+" "+str(dates['date range'][0,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][0,0])+" to "+subfun_monthNum_to_word(dates['date range'][1,1])[0]+" "+str(dates['date range'][1,2])+" (Day "+str(dates['date range dayNum'][0,1])+"), "+str(dates['date range'][1,0]),fontproperties=settings['plot']['font axis label FM']); #set the x axis label                
                    
                    #Draw time line on bottom OMNI plot
                    OMNI_data_min = np.min(data['OMNI'][settings['OMNI']['data type']]); #pre-calc the min
                    OMNI_data_max = np.max(data['OMNI'][settings['OMNI']['data type']]); #pre-calc the max
                    imOMNILine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( OMNI_data_min,OMNI_data_max)) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest               
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    
                    #Do the AMPERE plotting
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=AMPERE_colorMap, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    
                    #Draw time line on bottom OMNI plot
                    imOMNILine = ax[1].plot( np.array(( (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 , (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600 )) , np.array(( OMNI_data_min,OMNI_data_max)) , color=settings['map']['site marker color'] , zorder = 10);
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    #END IF
                    
                #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                movie_writer.grab_frame(); #get the frame and save it
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imAMP.remove();
                imOMNILine.pop(0).remove();
    #            imOverlay.pop(0).remove();
    #            imMillstone.pop(0).remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                #END IF
                
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 10 =================================================================
    elif(settings['movie']['movie type'] == 10): #stationary data points (through averaging) + AMPERE data on same plot (no time average for TEC) + AMPERE side plot + TEC side plot
        if( np.isscalar(settings['TEC']['plot lim']) == 1 ):
            gif_TEC_plotLimValu = np.array( (-settings['TEC']['plot lim'],settings['TEC']['plot lim']) ); #make it a vector
        else:
            gif_TEC_plotLimValu = settings['TEC']['plot lim']; #keep it the same
        #END IF
        # AMPERE_dataRate = np.median(np.diff(data['AMPERE']['time unique'])); #days, get the median data rate for AMPERE data (avoids outliers)
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
            
            #calc the TEC keogram first, only need to calc it once      
            movie_TEC_keo2, settings_keo2 = GRITI_keo_keogrammer(data['TEC']['dTEC'] ,data['TEC']['time'], data['TEC']['lat'], data['TEC']['long'], data['TEC']['time unique'], data['TEC']['time unique'], dates,
                    settings_keo2, settings['paths'], settings_keo2_map, settings['plot'],
                    FLG_disablePlot=2,FLG_disableText=0,FLG_disableCache=0,FLG_useRightExact=1,FLG_raytraceMethod='cos',FLG_raytraceExactSqrt=True);
            
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                
                #------------Corral the AMPERE data to the right place-------------
                if( settings['movie']['use time delays'] == 1 ):
                    k = np.where( ((data['AMPERE']['time']+settings['AMPERE']['delay wrt TEC']*3600) <= data['TEC']['time unique'][i]) & ((data['AMPERE']['time']+settings['AMPERE']['delay wrt TEC']*3600) >= (data['TEC']['time unique'][i]-data['AMPERE']['data rate'])) )[0]; #get where the time point is, make sure it is within the data rate window
                else:
                    k = np.where( (data['AMPERE']['time'] <= data['TEC']['time unique'][i]) & (data['AMPERE']['time'] >= (data['TEC']['time unique'][i]-data['AMPERE']['data rate'])) )[0]; #get where the time point is, make sure it is within the data rate window
                #END IF
                AMPERE_data_portion = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                AMPERE_lat_portion = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                AMPERE_long_portion = data['AMPERE']['long'][k]; #degc, corresponding long values
                
                #----------------Corral the data to the right place----------------
                k =  np.where( data['TEC']['time'] == data['TEC']['time unique'][i])[0]; #gets during a time period
                vTEC_portion = data['TEC']['dTEC'][k]; #pull out the vTEC now
                pplat_portion = data['TEC']['lat'][k]; #get the pplat (pierce-point lat) at the time required
                pplong_portion = data['TEC']['long'][k]; #get the pplong (pierce-point long) at the time required
                            
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(pplat_portion,pplong_portion,vTEC_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat spaces'],settings['movie']['grid long spaces'],settings['movie']['grid lat delta'],settings['movie']['grid long delta'],settings['movie']['data reject ratio'],settings['movie']['data reject perc lim'],settings['movie']['data reject ratio max']);
                #call a numba'd function that makes the movie quicker by crunching the numbers gooder
                
                #----------------------------Tack on Title-------------------------
                string_title = 'TEC Global Plot, Time =  '+'{0:.2f}'.format(np.round((data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0]); #create mecha title
                
                ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.025); #set the title, properties always needed
                                
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                                   
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    cbar = movie_dict['fig'].colorbar(imTEC, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
                    # cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbar.mappable.set_clim(vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']));
                    cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),6)); #create useful tick marks
                    
                    #Do the AMPERE plotting
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
                        k2 = AMPERE_data_portion >= np.min(settings['AMPERE']['plot lim']); #make sure only to plot the values above the min
                        imAMP = ax[0].scatter(AMPERE_latLongMapped[0][k2],AMPERE_latLongMapped[1][k2],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion[k2],cmap=settings['AMPERE']['colormap'], vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    else:
                        imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=settings['AMPERE']['colormap'], zorder=7);
                    #END IF
                    cbar2 = movie_dict['fig'].colorbar(imAMP, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    #                cbar2.ax.set_yticklabels(np.linspace(np.min(settings['AMPERE']['plot lim']),np.max(settings['AMPERE']['plot lim']),5)); #create useful tick marks
                    # cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    cbar2.set_label(settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]); #tabel the colorbar
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    
                    if( settings['movie']['use time delays'] == 1 ):
                        cbar2.set_label(settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]+' (Delayed '+textNice(settings['AMPERE']['delay wrt TEC'])+' hr)'); #tabel the colorbar
                    else:
                        cbar2.set_label(settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]); #tabel the colorbar
                    #END IF
                    cbar2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
                        cbar2.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
                    #END IF
                    if( np.all(np.mod(cbar2.get_ticks(),1) == 0) ):
                        cax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f')); #force a rounded format
                    else:
                        cax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    #END IF    
                    cax2.yaxis.set_ticks_position('left'); #move it to the left
                    cax2.yaxis.set_label_position('left'); #move it to the left
    
                    #Do the AMPERE plotting
                    #----- Integrate AMPERE Data -----    
                    # AMPERE_integrated = GRITI_AMPERE_integrator(data['AMPERE'], dates, settings['AMPERE'], settings['map']['lat range'], settings['map']['long range'], settings['AMPERE']['integrate method'], settings['AMPERE']['integrate method lat val'], AMPERE_integrateMethod_log=settings['AMPERE']['integrate method log']); #integrate with the integrator function
                    AMPERE_integrated = data['AMPERE']['integrated']; #alias
                    
                    #--- Time match to 6 minutes if needed ---
                    if( np.isclose(data['AMPERE']['data rate'],360.) == False ):
                        sixMin_timeUnique = np.arange(dates['date range zero hr hour bounds'][0]*3600,dates['date range zero hr hour bounds'][1]*3600,360); #sec, arange time stamps in 6 minute steps
                        AMPERE_integrated, AMPERE_timeUnique_matched = subfun_timeMatch(AMPERE_integrated, data['AMPERE']['time unique'], sixMin_timeUnique, timeMatch_delta=data['AMPERE']['data rate'], FLG_useSum=1); #time match alg to align to 6 minute cadence, add because it's a count (?)
                        AMPERE_timeUnique_hr = (AMPERE_timeUnique_matched-dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                    else:
                        AMPERE_timeUnique_hr = (data['AMPERE']['time unique']-dates['date range zero hr dayNum'][1]*86400)/3600; #hr, convert to hr with 0 hr at specified day
                    #END IF
                    
                    #--- Prepare axis limits ---
                    if( np.mod(np.round(np.min(AMPERE_timeUnique_hr)),2) == 0 ):
                        AMPERE_time_hr_axis_min = np.round(np.min(AMPERE_timeUnique_hr)); #is even, good to go
                    else:
                        AMPERE_time_hr_axis_min = np.round(np.min(AMPERE_timeUnique_hr))+1; #is odd, make even
                    #END IF
                    if( np.mod(np.round(np.max(AMPERE_timeUnique_hr)),2) == 0 ):
                        AMPERE_time_hr_axis_max = np.round(np.max(AMPERE_timeUnique_hr)); #is even, good to go
                    else:
                        AMPERE_time_hr_axis_max = np.round(np.max(AMPERE_timeUnique_hr))-1; #is odd, make even
                    #END IF         
                    
                    #-----BEGIN THE PLOTTING!------
                    AMPERE_plot_label = settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]; #get the label
                    AMPERE_plot_label_noUnits = settings['AMPERE']['labels'][settings['AMPERE']['data type']]; #remove the (units)
                    
                    imAMPplt = ax[1].plot( AMPERE_timeUnique_hr+settings['AMPERE']['delay wrt TEC'], AMPERE_integrated , linewidth=settings['plot']['line width']['regular'] ); #plot
                    
                    xAxisTicks = np.arange(AMPERE_time_hr_axis_min,AMPERE_time_hr_axis_max+8,8); #sets the start hr, stop hr, and the step size between (in this case, 4 hr)
                    ax[1].set_xticks(xAxisTicks); #set x axis ticks
                    
                    ax[1].set_xlim( AMPERE_time_hr_axis_min , AMPERE_time_hr_axis_max ); #set y axis limits
                    
                    ax[1].set_ylabel(settings['AMPERE']['units'][settings['AMPERE']['data type']],fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    ax[1].yaxis.set_label_position('right');
                    
                    ax[1].set_ylim( np.min(AMPERE_integrated) , np.max(AMPERE_integrated) ); #set y axis limits
                    # ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.0e')); #force a rounded format
                    
                    ax[1].grid(b=True, which='major', axis='both', color='xkcd:light grey'); #sets major axis grid lines to be on        
                    
                    string_title = 'Int. '+AMPERE_plot_label_noUnits; #prep the title
                    if( settings['AMPERE']['integrate method'] == 0 ):
                        string_title = string_title + ' within keo area'; #add to mecha title
                    elif( settings['AMPERE']['integrate method'] == 1 ):
                        string_title = string_title + ' within keo long & up to pole'; #add to mecha title
                    elif( settings['AMPERE']['integrate method'] == 2 ):
                        string_title = string_title + ' within keo long & up to '+str(settings['AMPERE']['integrate method lat val'])+' degc lat'; #add to mecha title
                    elif( settings['AMPERE']['integrate method'] == 3 ):
                        if( (np.min(settings['map']['lat range']) < 0) & (np.max(settings['map']['lat range']) >= 0) ):
                            string_title = string_title + ' both Hemis.'; #add to mecha title
                        else:
                            if( (np.min(settings['map']['lat range']) >= 0) & (np.max(settings['map']['lat range']) >= 0) ):
                                #northern hemisphere
                                string_title = string_title + ' Northern Hemi.'; #add to mecha title
                            else:
                                #southern hemisphere
                                string_title = string_title + ' Southern Hemi.'; #add to mecha title
                            #END IF
                        #END IF
                    #END IF
                    ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    
                    # ax[1].set_xlabel('Time in UT [hr] - 0 Hr on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr']_dayPostfix+' | Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=settings['plot']['font axis label FM']); #set the x axis label          
                    ax[1].set_xlabel('Time in UT [hr]'); #just a bit
                    
                    #Draw time line on bottom AMPERE plot
                    imAMPpltLine = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    #------Now plot the keogram2--------
                    pltHelprX, pltHelprY = np.meshgrid( (data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600, \
                                settings_keo2['keo plot latlong chunks']);
                    imKeo2 = ax[2].pcolormesh(pltHelprX, pltHelprY,  movie_TEC_keo2.T ,vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu),cmap='jet'); # pseudocolor plot "stretched" to the grid
                    cbarKeo2 = movie_dict['fig'].colorbar(imKeo2, cax=caxKeo2, orientation='vertical'); #create a colorbar using the prev. defined cax
                    cbarKeo2.set_label('delta-vTEC [TECU]'); #tabel the colorbar
                    cbarKeo2.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    cbarKeo2.mappable.set_clim(vmin=np.min(gif_TEC_plotLimValu), vmax=np.max(gif_TEC_plotLimValu));
                    caxKeo2.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    caxKeo2.yaxis.set_ticks(np.linspace(np.min(gif_TEC_plotLimValu),np.max(gif_TEC_plotLimValu),6)); #create useful tick marks
                    caxKeo2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                    
                    #    string_title = 'TEC Averaged on Angle of '+str(np.round(keo,2))+' deg and Width of '+ \
                    #        str(np.round(settings['TEC']['keo']['keo width'],2))+' arcdeg, Avg Step # = '+str(keo_N)+ \
                    #        ' arcdeg, Line Shows '+settings['TEC']['keo']['keo plot latlong name']+' of Millstone Hill Zenith Beam'; #create mecha title
                    string_title = 'delta-vTEC Keogram of '+settings['movie']['side plots']['keo 2 name']; #create mecha title
                    ax[2].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
                    ax[2].set_xlabel('0 Hr on Day '+str(dates['date range zero hr dayNum'][1])+', '+str(dates['date range zero hr dayNum'][0]),fontproperties=settings['plot']['font axis label FM']); #set the x axis label
                    ax[2].set_ylabel(settings_keo2['keo plot latlong name']+' [arcdeg]',fontproperties=settings['plot']['font axis label FM']); #set the y axis label
                    
                    xAxisTicks = np.arange( (np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][0]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) , \
                            (np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600) - np.mod(np.round((data['TEC']['time unique'][-1]-dates['date range zero hr dayNum'][1]*86400)/3600),2)) + 2 , \
                            8); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
                    ax[2].set_xticks(xAxisTicks); #set x axis ticks
                    
                    keo_Range_Chunks_Long_Plot_autoTick = (np.ceil(np.max(settings_keo2['keo plot latlong chunks'])) - np.floor(np.min(settings_keo2['keo plot latlong chunks'])))/13; #tries to split the latitude range into 13 parts (based off of 180/15+1)
                    if( keo_Range_Chunks_Long_Plot_autoTick > 25 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 30; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 10 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 15; #sets the tick setting to 15 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 5 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 10; #sets the tick setting to 10 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 2 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 5; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick > 1 ):
                        keo_Range_Chunks_Long_Plot_autoTick = 2; #sets the tick setting to 5 arcdegrees per tick
                    elif( keo_Range_Chunks_Long_Plot_autoTick >= 0.6 ): #0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                        keo_Range_Chunks_Long_Plot_autoTick = 1; #                                                        sets the tick setting to 1 arcdegree per tick
                    else:
                        if(settings_keo2['keo plot latlong name'] == 'Latitude'): #if Y axis is latitude, use latitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 2 lat range']) - np.min(settings['movie']['side plots']['keo 2 lat range']))/13; #just goes for it if it's a super tiny range
                        elif(settings_keo2['keo plot latlong name'] == 'Longitude'): #if Y axis is longitude, use longitude
                            keo_Range_Chunks_Long_Plot_autoTick = (np.max(settings['movie']['side plots']['keo 2 long range']) - np.min(settings['movie']['side plots']['keo 2 long range']))/13; #just goes for it if it's a super tiny range
                        #END IF
                    #END IF
                    yAxisTicks = np.round(np.arange( np.floor(np.min(settings_keo2['keo plot latlong chunks'])),np.ceil(np.max(settings_keo2['keo plot latlong chunks']))+keo_Range_Chunks_Long_Plot_autoTick,keo_Range_Chunks_Long_Plot_autoTick ),2); #creates y ticks automagically
                    ax[2].set_yticks(yAxisTicks); #set x axis ticks
                    ax[2].set_ylim( np.round(ax[2].get_ylim()) ); #avoid weird rounding errors
                    
                    #Now drawing line of interest
                    if( settings_keo2['keo plot latlong name'] == 'Longitude' ): #if true, longitude
                        if( (np.min(settings['movie']['side plots']['keo 2 long range']) <= settings['map']['site coords'][0,1]) & (np.max(settings['movie']['side plots']['keo 2 long range']) >= settings['map']['site coords'][0,1]) ): #only plot if it's in the long range specified
                            ax[2].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,1],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    else: #else latitude
                        if( (np.min(settings['movie']['side plots']['keo 2 lat range']) <= settings['map']['site coords'][0,0]) & (np.max(settings['movie']['side plots']['keo 2 lat range']) >= settings['map']['site coords'][0,0]) ): #only plot if it's in the lat range specified
                            ax[2].plot( np.linspace(np.min((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),np.max((data['TEC']['time unique'] - dates['date range zero hr dayNum'][1]*86400)/3600),10,endpoint=True) , #X time hr
                                    np.tile(settings['map']['site coords'][0,0],10) , #Y latitude OR longitude arcdeg
                                    c='xkcd:black',linewidth=1); #plots a point with a black line
                        #END IF
                    #END IF
                    
                    #-----plot vertical line showing the time on the keogram2-----
                    hVline2 = ax[2].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax[0].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i]/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax[0].text(x,y,'SUN\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #Now drawing line of interest               
                    Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    imMillstone = ax[0].plot(Millstone_latLongMapped[0],Millstone_latLongMapped[1],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], zorder=50); #plot this, 50 always on top
                                    
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the TEC plotting
                    pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    TEC_latLongMapped = geoMap(pltHelprX,pltHelprY); #convert the lat/long arcdeg to the current map coordinates
                    imTEC = ax[0].pcolormesh(TEC_latLongMapped[0], TEC_latLongMapped[1],  gif_Grid.T ,vmin=np.min(settings['TEC']['plot lim']), vmax=np.max(settings['TEC']['plot lim']),cmap='jet',zorder=5); # pseudocolor plot "stretched" to the grid
                    
                    #Do the AMPERE plotting
                    AMPERE_latLongMapped = geoMap(AMPERE_long_portion,AMPERE_lat_portion); #convert the lat/long arcdeg to the current map coordinates
                    if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
                        k2 = AMPERE_data_portion >= np.min(settings['AMPERE']['plot lim']); #make sure only to plot the values above the min
                        imAMP = ax[0].scatter(AMPERE_latLongMapped[0][k2],AMPERE_latLongMapped[1][k2],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion[k2],cmap=settings['AMPERE']['colormap'], vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), zorder=7);
                    else:
                        imAMP = ax[0].scatter(AMPERE_latLongMapped[0],AMPERE_latLongMapped[1],s=settings['map']['AMPERE scatter size'],c=AMPERE_data_portion,cmap=settings['AMPERE']['colormap'], zorder=7);
                    #END IF
                    
                    #Draw time line on bottom OMNI plot
                    imAMPpltLine = ax[1].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    hVline2 = ax[2].axvline( x=(data['TEC']['time unique'][i] - dates['date range zero hr dayNum'][1]*86400)/3600 ,c='xkcd:black',linewidth=1, linestyle='--'); #plot a vertical line to show the current time
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['TEC']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax[0].text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(data['TEC']['time unique'][i]/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((data['TEC']['time unique'][i]/86400-np.int64(data['TEC']['time unique'][i]/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax[0], zorder=2);
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax[0].transAxes); #make a sun figure
                            hSun = ax[0].add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax[0].set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                movie_writer.grab_frame(); #get the frame and save it
                            
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imTEC.remove();
                imAMP.remove();
                imAMPpltLine.remove();
                hVline2.remove();
                # imOMNILine.pop(0).remove();
    #            imOverlay.pop(0).remove();
    #            imMillstone.pop(0).remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                elif( settings['movie']['day nite line'] == 2): #only delete if it is there
                    for middleman in dayNite_shade.collections: #this needs an extra helping hand
                        middleman.remove();
                #END IF
                if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
                    hSun.remove();
                #END IF
            #END FOR i
        #END WITH
    #================================================================= MOVIE TYPE 11 =================================================================
    elif(settings['movie']['movie type'] == 11): #ONLY AMPERE data on same plot
        AMPERE_lat_delta = data['AMPERE']['data info']['lat delta']; #it's provided
        AMPERE_long_delta = data['AMPERE']['data info']['long delta']; #it's provided
        
        #Stuff for repeating a frame - needed if want a constant 30 FPS
        #This needs to go in settings['movie']['movie type'] == 2 area
        gif_DesiredTotalTime = (data['AMPERE']['time unique'].size-1)/gif_DesiredFPS; #s, total time from desired FPS
        gif_DesiredReqFrameNum = gif_DesiredTotalTime*30; #frames, number of frames needed to go for the time wanted
        gif_RepeatFrame = gif_DesiredReqFrameNum/(data['AMPERE']['time unique'].size-1); #number of frames to repeat every time (exact w/ decimal)
        if( gif_RepeatFrame > 0.5 ):
            gif_RepeatFrame = np.int64(np.round(gif_RepeatFrame)); #round it
        else: #if less than 0.5 will round to 0 - so checking to see if there's a ton of frames and need to just boost it to 60 to keep run time in check
            if( ( (data['AMPERE']['time unique'].size-1)/30 > settings['movie']['desired max run time'] ) & ( settings['movie']['disable FPS shift'] == 0 ) ): #try to keep run time down
                gif_vidFrameRate = 60; #set the FPS to 60
            #END IF
            gif_RepeatFrame = 1; #I guess to make sure it works
        #EMD IF
        
        with movie_writer.saving(movie_dict['fig'], movie_name, settings['movie']['fig PPI']): #figure to save, file name, and PPI (DPI in documentation)
                        
            #------------Prep reused vars-------------
            settings['movie']['grid lat'] = np.arange(np.min(settings['map']['lat range']),np.max(settings['map']['lat range'])+AMPERE_lat_delta,AMPERE_lat_delta); #degc, create lat points (+1 lets us use delta to edge wanted range - yields correct # of spaces)
            settings['movie']['grid long'] = np.arange(np.min(settings['map']['long range']),np.max(settings['map']['long range'])+AMPERE_long_delta,AMPERE_long_delta); #degc, create long points specilized for AMPERE
            pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
            if( 'stere' in settings['map']['projection name'] ):
                if( settings['movie']['spin'] != 0):
                    #this is a constant sun since the plot is continually spun to match
                    x = (1.05*0.5*np.cos(90))+0.5; #geoMap coordinate 
                    y = (1.05*0.5*np.sin(90))+0.5; #geoMap coordinate 
                    circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                    hSunStatic = ax.add_artist(circleSun); #plot the sun
                #END IF
            #END IF
                    
            #-------------------------Start Making Pictures------------------------
            for i in range(movie_timeRange[0],movie_timeRange[1]):
                
                #------------Corral the AMPERE data to the right place-------------
                if( settings['movie']['use time delays'] == 1 ):
                    k = np.where( ((data['AMPERE']['time']+settings['AMPERE']['delay wrt TEC']*3600) <= data['AMPERE']['time unique'][i]) & ((data['AMPERE']['time']+settings['AMPERE']['delay wrt TEC']*3600) > (data['AMPERE']['time unique'][i]-data['AMPERE']['data rate'])) )[0]; #get where the time point is, make sure it is within the data rate window
                else:
                    k = np.where( (data['AMPERE']['time'] <= data['AMPERE']['time unique'][i]) & (data['AMPERE']['time'] > (data['AMPERE']['time unique'][i]-data['AMPERE']['data rate'])) )[0]; #get where the time point is, make sure it is within the data rate window
                #END IF
                AMPERE_data_portion = data['AMPERE'][settings['AMPERE']['data type']][k]; #ergs/(cm^2*sec), get the Joule Heating for the current time step
                AMPERE_lat_portion = data['AMPERE']['lat'][k]; #degc, corresponding lat values
                AMPERE_long_portion = data['AMPERE']['long'][k]; #degc, corresponding long values
                
                gif_Grid = GRITI_movieMaker_subfun_dataGridder(AMPERE_lat_portion,AMPERE_long_portion,AMPERE_data_portion,settings['movie']['grid lat'],settings['movie']['grid long'],settings['movie']['grid lat'].size-1,settings['movie']['grid long'].size-1,AMPERE_lat_delta,AMPERE_long_delta,settings['movie']['data reject ratio'],101,settings['movie']['data reject ratio max']).T; #101 disables the data rejection stuff b/c AMPERE doesn't need it
                # gif_Grid = np.vstack((np.hstack((np.roll(AMPERE_data_portion.reshape(settings['movie']['grid lat'].size-1,settings['movie']['grid long'].size-1),60,axis=1),np.nan*np.ones((settings['movie']['grid lat'].size-1,1)))),np.nan*np.ones((1,settings['movie']['grid long'].size)))); # replicate what gif_Grid would do with reshaping b/c there's no need to avg pts
                
                #----------------------------Tack on Title-------------------------
                curr_time_mod = np.mod(data['AMPERE']['time unique'][i],86400);
                curr_time_date = subfun_dayNum_to_date((dates['date range zero hr dayNum'][0],data['AMPERE']['time unique'][i]//86400))[0]; #if year changes needs overhaul everywhere
                string_title = '{0:.2f}'.format(np.round((data['AMPERE']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600,2))+\
                    ', 0 UT on '+dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+', '+str(dates['date range zero hr dayNum'][0])+\
                    ' | ('+str(curr_time_date[0])+'/'+str(curr_time_date[1])+'/'+str(curr_time_date[2])+', '+str(curr_time_mod//3600).zfill(2)+\
                    ':'+str((curr_time_mod-curr_time_mod//3600*3600)//60).zfill(2)+':'+str(curr_time_mod-(curr_time_mod-curr_time_mod//3600*3600)//60*60-curr_time_mod//3600*3600).zfill(2)+')'; #create mecha title
                if( settings['movie']['use time delays'] == 1 ):
                    string_title += '[Time Shifted by '+settings['AMPERE']['delay wrt TEC']+' hrs]'; 
                #END IF
                
                ax.set_title(string_title,fontproperties=settings['plot']['font title FM'],y=movie_title_yOffset); #set the title, properties always needed
                                
                #-------------------Starting the Plotting--------------------------
                if( i == movie_timeRange[0] ): #first run preps axes, color bars, etc.
                                                   
                    #Do the AMPERE plotting
                    if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
                        imAMP = ax.pcolormesh(pltHelprX, pltHelprY,  gif_Grid, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), cmap=settings['AMPERE']['colormap'],zorder=5, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid, data is pre .T'd
                    else:
                        imAMP = ax.pcolormesh(pltHelprX, pltHelprY,  gif_Grid, cmap=settings['AMPERE']['colormap'],zorder=5, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid, data is pre .T'd
                    #END IF
                    cbar = movie_dict['fig'].colorbar(imAMP, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
                    # cax.yaxis.set_major_formatter(FsormatStrFormatter('%.2f')); #force a rounded format
                    cbar.set_label(settings['AMPERE']['labels'][settings['AMPERE']['data type']]+settings['AMPERE']['units'][settings['AMPERE']['data type']]); #tabel the colorbar
                    cbar.ax.tick_params(labelsize=settings['plot']['font axis tick']);
                    if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
                        cbar.mappable.set_clim(vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']));
                    #END IF
                    if( np.all(np.mod(cbar.get_ticks(),1) == 0) ):
                        cax.yaxis.set_major_formatter(FormatStrFormatter('%.0f')); #force a rounded format
                    else:
                        cax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    #END IF 
                    # cax.yaxis.set_ticks(np.linspace(np.min(settings['TEC']['plot lim']),np.max(settings['TEC']['plot lim']),6)); #create useful tick marks
                    
                    if( settings['movie']['day nite line'] == 1 ):
                        movie_dict['fig'].canvas.flush_events(); #this is req. to get get_window_extent() to get the current window size
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['AMPERE']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        #constants will use all the time - only for plotting of day/nite line so minor importance
                        #this stuff makes the text angle plotted mostly correct most of the time
                        bboxFig = movie_dict['fig'].get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the entire figure dimensions
                        bboxAx0 = ax.get_window_extent().transformed(movie_dict['fig'].dpi_scale_trans.inverted()); #get the plot dimensions
    #                    plot_ratio = (bboxAx0.width/bboxAx0.height)/( (bboxAx0.width/bboxAx0.height)/(bboxmovie_dict['fig'].width/bboxmovie_dict['fig'].height) )**2; #get the plot ratio, will use it to fix up the angle
                        plot_ratio = (bboxAx0.width/bboxAx0.height); #get the plot ratio, will use it to fix up the angle
                        dayNite_textRotationLenOrig = 35; #length to go up and down the sunrise/sunset line to estimate an angle
                        dayNite_savgolFiltLenOrig = 101; #set it as a constant to start off
                        dayNite_textLatAbs = (np.max(settings['map']['lat range'])-np.min(settings['map']['lat range']))*.9 + np.min(settings['map']['lat range']); #calc like 80% of the max latitude
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax[0].plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['AMPERE']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #EMD IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        if( settings['map']['coord type'] == 'geo' ):
                            shifted = data['AMPERE']['time unique'][i]; #no shift
                        else:
                            #this ain't perfect, but it's plenty good for general where not sun
                            shifter = (np.mod(settings['movie']['sun loc']['long'][i],360) - np.mod(settings['movie']['sun loc']['long geo'][i],360))*240 #sec, effective time to shift by so Sun is in right spot; 240 = 24*3600/360
                            shifted = data['AMPERE']['time unique'][i] - shifter; #shift the time so Sun is in the right spot
                        #END IF
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(shifted/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((shifted/86400-np.int64(shifted/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((shifted/86400-np.int64(shifted/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((shifted/86400-np.int64(shifted/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        # dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax, zorder=25); #basemap version
                        dayNite_shade = ax.add_feature(Nightshade(dayNite_currentTime, alpha=0.25),zorder=25); #nighttime shading, only relevant for geographic
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax[0].text(x,y,'SUN\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                            hSun = ax.add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax.set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                    #Now drawing point of interest
                    # Millstone_latLongMapped = geoMap(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0]); #convert the lat/long arcdeg to the current map coordinates
                    if( (settings['map']['site coords'][0,0] <= np.max(settings['map']['lat range'])) & (settings['map']['site coords'][0,0] >= np.min(settings['map']['lat range'])) & (settings['map']['site coords'][0,1] <= np.max(settings['map']['long range'])) & (settings['map']['site coords'][0,1] >= np.min(settings['map']['long range'])) ):
                        imMillstone = ax.plot(settings['map']['site coords'][0,1],settings['map']['site coords'][0,0],marker=settings['map']['site marker type'], color=settings['map']['site marker color'], markersize=settings['map']['site marker size'], transform=cartopy.crs.PlateCarree(), zorder=50); #plot this, 50 always on top
                    else:
                        if( 'stere' in settings['map']['projection name'] ):
                            diameter_offset = 0.09;
                            x = ((1+diameter_offset)*0.5*np.cos((settings['map']['site coords'][0,1][0]-90)*np.pi/180))+0.5; #geoMap coordinate 
                            y = ((1+diameter_offset)*0.5*np.sin((settings['map']['site coords'][0,1][0]-90)*np.pi/180))+0.5; #geoMap coordinate
                            imMillstone = ax.arrow(x, y, (diameter_offset/2)*np.cos((settings['map']['site coords'][0,1][0]+90)*np.pi/180), (diameter_offset/2)*np.sin((settings['map']['site coords'][0,1][0]+90)*np.pi/180), width=0.007, head_width=0.007*3.5, head_length=0.007*3.5*1.00, length_includes_head=True, color=settings['map']['site marker color'], clip_on=False, transform=ax.transAxes); #plot this, 50 always on top
                        #END IF
                        #if needed, code something for rectangular plots but idk what it should be exactly rn
                    #END IF
                    
                    # figFitter(fig); #fit that fig fast [already been done, slowly]
                else:
                    #Just plotting now - prep done, hopefully speeds it!
                    #Do the AMPERE plotting
                    # pltHelprX, pltHelprY = np.meshgrid( settings['movie']['grid long'], settings['movie']['grid lat']); #helps the pcolor work
                    if( ~np.any(np.isinf(settings['AMPERE']['plot lim'])) ):
                        imAMP = ax.pcolormesh(pltHelprX, pltHelprY,  gif_Grid, vmin=np.min(settings['AMPERE']['plot lim']), vmax=np.max(settings['AMPERE']['plot lim']), cmap=settings['AMPERE']['colormap'],zorder=5, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid
                    else:
                        imAMP = ax.pcolormesh(pltHelprX, pltHelprY,  gif_Grid, cmap=settings['AMPERE']['colormap'],zorder=5, transform=cartopy.crs.PlateCarree()); # pseudocolor plot "stretched" to the grid
                    #END IF

                    if( np.all(np.mod(cbar.get_ticks(),1) == 0) ):
                        cax.yaxis.set_major_formatter(FormatStrFormatter('%.0f')); #force a rounded format
                    else:
                        cax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')); #force a rounded format
                    #END IF 
                    if( settings['movie']['day nite line'] == 1 ):
                        #Plot the sunrise/sunset terminators
                        dayNite_temp = np.abs(dayNite_sunrise - (data['AMPERE']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunrise locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_day = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_day > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(0,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_day = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='g',linewidth=1.75,zorder=6); #plots a line to show the sunrise time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_day = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunrise", rotation=dayNite_textRotation, color='g', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                        
                        dayNite_temp = np.abs(dayNite_sunset - (data['AMPERE']['time unique'][i]-dates['date range zero hr dayNum'][1]*86400)/3600) <= 10/3600; #gets the sunset locations that are within 10 sec of the current TEC time (2 min had way too many hits, zigzags)
                        hLgd_FLG_nite = np.sum(dayNite_temp); #plotting flag, shows when legend is active
                        if( hLgd_FLG_nite > 0 ): #only do work if it is there
                            #calc all that day/nite stuff in one function to keep it from getting cluttered
                            dayNite_Lat_line, dayNite_Long_line, dayNite_Lat_text, dayNite_Long_text, dayNite_textRotation = GRITI_movieMaker_subfun_dayniteCalc(1,np.sum(dayNite_temp),dayNite_Grid_Long[dayNite_temp],dayNite_Grid_Lat[dayNite_temp],dayNite_textLatAbs,plot_ratio,dayNite_savgolFiltLenOrig,dayNite_textRotationLenOrig);
                            dayNite_latLongMapped = geoMap(dayNite_Long_line,dayNite_Lat_line); #convert the lat/long arcdeg to the current map coordinates
                            imDayNite_nite = ax.plot(dayNite_latLongMapped[0],dayNite_latLongMapped[1],color='b',linewidth=1.75,zorder=6); #plots a line to show the sunset time
                            if( settings['movie']['day nite text'] == 1 ):
                                dayNite_latLongMapped = geoMap(dayNite_Long_text,dayNite_Lat_text); #convert the lat/long arcdeg to the current map coordinates
                                textDayNite_nite = ax.text(dayNite_latLongMapped[0], dayNite_latLongMapped[1], "Sunset", rotation=dayNite_textRotation, color='b', fontsize=settings['plot']['font axis tick'], zorder=6);
                            #END IF
                        #END IF
                    elif( settings['movie']['day nite line'] == 2 ):
                        if( settings['map']['coord type'] == 'geo' ):
                            shifted = data['AMPERE']['time unique'][i]; #no shift
                        else:
                            #this ain't perfect, but it's plenty good for general where not sun
                            shifter = (np.mod(settings['movie']['sun loc']['long'][i],360) - np.mod(settings['movie']['sun loc']['long geo'][i],360))*240 #sec, effective time to shift by so Sun is in right spot; 240 = 24*3600/360
                            shifted = data['AMPERE']['time unique'][i] - shifter; #shift the time so Sun is in the right spot
                        #END IF
                        dayNite_currentDay = subfun_dayNum_to_date( np.array((dates['date range zero hr dayNum'][0],np.int64(shifted/86400)),ndmin=2))[0]; #get the current yr/month/day
                        dayNite_currentHr = np.int64((shifted/86400-np.int64(shifted/86400))*24); #hr, get the current hour
                        dayNite_currentMin = np.int64( ((shifted/86400-np.int64(shifted/86400))*24 - dayNite_currentHr)*60); #min, get the current min
                        dayNite_currentSec = np.int64( (((shifted/86400-np.int64(shifted/86400))*24 - dayNite_currentHr)*60 - dayNite_currentMin)*60); #min, get the current min
                        dayNite_currentTime = datetime(dayNite_currentDay[0], dayNite_currentDay[1], dayNite_currentDay[2], dayNite_currentHr, dayNite_currentMin,dayNite_currentSec);
                        # dayNite_shade = geoMap.nightshade(dayNite_currentTime, color='k', delta=0.25, alpha=0.25, ax=ax, zorder=25);
                        dayNite_shade = ax.add_feature(Nightshade(dayNite_currentTime, alpha=0.25),zorder=25); #nighttime shading, only relevant for geographic
                    #END IF
                    
                    if( 'stere' in settings['map']['projection name'] ):
                        if( settings['movie']['spin'] == 0):
#                            x = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
#                            y = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]+np.pi))+0.5; #geoMap coordinate 
                            x = (1.05*0.5*np.cos(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            y = (1.05*0.5*np.sin(settings['movie']['sun loc'][i]))+0.5; #geoMap coordinate 
                            #hSun = ax.text(x,y,'SUN\N{DEGREE SIGN}',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM']);
                            circleSun = plt.Circle((x, y), radius=0.015, color='xkcd:sun yellow', clip_on=False,transform=ax.transAxes); #make a sun figure
                            hSun = ax.add_artist(circleSun); #plot the sun
                        else:
                            hSun = ax.set_theta_offset(settings['movie']['sun loc'][i]); #turn the whole plot so top is where the sun is
                            #I'm not sure if I can get this working easily - not a lot of optons.
                        #END IF
                    #END IF
                    
                #END IF
            
                #-----------------------Create Movie/GIF---------------------------
                #Makes the gif now
                plt.draw();
                
                for j in range(0,gif_RepeatFrame): #runs this as many times as needed to get 30 FPS or whatever
                    movie_writer.grab_frame(); #get the frame and save it
                #END FOR j
                
                #-------------------Clean up for re-use----------------------------
                #if forget one (like hOverlay) slows it way down after many plots
                imAMP.remove();
                if( settings['movie']['day nite line'] == 1 ):
                    if(hLgd_FLG_day > 0): #only delete if it is there
                        imDayNite_day.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_day.remove();
                        #END IF
                    #END IF
                    if(hLgd_FLG_nite > 0): #only delete if it is there
                        imDayNite_nite.pop(0).remove();
                        if( settings['movie']['day nite text'] == 1 ):
                            textDayNite_nite.remove();
                        #END IF
                    #END IF
                elif( settings['movie']['day nite line'] == 2): #only delete if it is there
                    # for middleman in dayNite_shade.collections: #this needs an extra helping hand
                    #     middleman.remove();
                    dayNite_shade.remove(); #ditch it
                #END IF
                if( 'stere' in settings['map']['projection name'] ): #only delete if it is there
                    hSun.remove();
                #END IF
            #END FOR i
        #END WITH
    #END IF

    figDump = pickle.dumps(fig, protocol=pickle.HIGHEST_PROTOCOL);
    plt.close(fig); #close it down
    return figDump
#END DEF