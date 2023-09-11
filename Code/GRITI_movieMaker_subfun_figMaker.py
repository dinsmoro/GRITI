"""
Creates the figure to paint the movie on
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Code.GRITI_plotHelper_area_init import GRITI_plotHelper_area_init
from Code.subfun_figFitter import figFitter

def GRITI_movieMaker_subfun_figMaker(data, settings, time4mag, movie_figOffsets=None):
    movie_dict = {}; #prep it
    if( settings['movie']['movie type'] in [0, 1, 4, 5, 11] ): #creates figure with 1 plot
                    
        if( settings['movie']['movie type'] in [4, 5] ): #creates figure with 1 plot (geomap) and 2 colorbars (left is AMPERE)
            FLG_double_cax = True;
        else:
            FLG_double_cax = False;
        #END IF
                
        if( settings['movie']['movie type'] != 11 ):
            alt4mag = np.median(data['TEC']['pierceAlt']);
        else:
            if( 'altitude' in settings['AMPERE'] ):
                alt4mag = settings['AMPERE']['altitude']; #use AMPERE altitude
            else:
                alt4mag = 120.; #default, great for auroral zone stuff (like field aligned currents)
            #END IF
        #END IF
        # fig = plt.figure(figsize = settings['movie']['fig size']/settings['movie']['fig PPI']); #use fig because we gotta init an axis later
        figsizer_override = settings['movie']['fig size']/settings['movie']['fig PPI']; #calc the size desired
        figPPI_override = settings['movie']['fig PPI'];
        fig, ax, cax = GRITI_plotHelper_area_init(settings['map']['lat range'], settings['map']['long range'], settings['map'], settings['plot'], 1, \
                                                  time4mag=time4mag, alt4mag=alt4mag, FLG_double_cax=FLG_double_cax, \
                                                  FLG_latLabel=True,FLG_longLabel=True, FLG_forcePolarLatLabels=True, \
                                                  figsizer_override=figsizer_override, figPPI_override=figPPI_override); #init a cartopy thing
        if( settings['movie']['movie type'] in [4, 5] ):
            cax2 = cax[1]; #unpack
            cax = cax[0]; #unpack (hope cax2 isn't donezo)
        #END IF
        
        string_title = 'Prep Title'; #create mecha title
        if( ('stere' in settings['map']['projection name']) ):
            movie_title_yOffset = 1.128; #more b/c words
        else:
            movie_title_yOffset = 1.035; #less, it chill
        #END IF
        ax.set_title(string_title,fontproperties=settings['plot']['font title FM'],y=movie_title_yOffset); #set the title
        
        if( np.all(movie_figOffsets == None) ):
            movie_figOffsets = figFitter(fig, fast=False, returnOffsets=True); #fit that fig slow, cause it's gonna be used lots
        else:
            fig.subplots_adjust(left = movie_figOffsets[0], right = movie_figOffsets[1], bottom = movie_figOffsets[2], top = movie_figOffsets[3]); #apply known ones instead
        #END IF
        if( ('stere' in settings['map']['projection name']) ):
            cax.yaxis.label.set_text(''); #clear out the text in case it isn't used
            cax2.yaxis.label.set_text(''); #clear out the text in case it isn't used
        #END IF
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
        
        # Load in references
        movie_dict['fig'] = fig;
        movie_dict['fig offsets'] = movie_figOffsets;
        movie_dict['ax'] = ax;
        movie_dict['cax'] = cax;
        if( settings['movie']['movie type'] in [4, 5] ):
            movie_dict['cax2'] = cax2;
        #END IF
        movie_dict['fig size'] = settings['movie']['fig size']; #it might be diff?
        movie_dict['title offset'] = movie_title_yOffset;
        
    elif( settings['movie']['movie type'] in [2, 3] ): #creates figure with 2 plots (TEC on top, RTI on bottom)
    
        fig, ax = plt.subplots(figsize = settings['movie']['fig size']/settings['movie']['fig PPI'],nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        divider = make_axes_locatable(ax[0]); #prep to add an axis (top one is TEC so it gets it)
        cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        divider2 = make_axes_locatable(ax[1]); #prep to add an axis (bottom is ISR so it gets one too)
        cax2 = divider2.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        #mill for square Mercator style
        #robin for oval shape
        geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
            resolution = 'i', area_thresh = 10000, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
            llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
            urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[0].set_aspect('auto');
        ax[1].set_aspect('auto');
        
        #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
        geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,settings['map']['long autotick']),2), 
            labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
        geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
            labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
        #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
        #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
        string_title = 'Prep Title'; #create mecha title
        ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the upper title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the lower title
        
        figFitter(fig, fast=False); #fit that fig slow, cause it's gonna be used lots
        # fig.subplots_adjust(left = 0.040, right = 0.95, top = 0.96, bottom = 0.055, hspace = 0.15); #sets padding to small numbers for minimal white space  
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
                
    if( (settings['movie']['movie type'] == 6) ): #creates figure with 1 plot (geomap) and 2 colorbars (left is AMPERE) and smaller plot next to geomap as well for a keogram or something
        fig = plt.figure(figsize = settings['movie']['fig size']/settings['movie']['fig PPI']); #use fig because we gotta init an axis later
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        ax = [plt.subplot2grid((1, 3), (0, 0), colspan=2,fig=fig) , plt.subplot2grid((7, 3), (2, 2),rowspan=3, colspan=1,fig=fig) ]; #init an axis
        #the subplot2grid is weird, 2nd one makes a 7x3 grid and puts a subplot that spans 3 of the 7 rows 2 rows down and in the 3rd column
        #1st one does a 1x3 grid and puts a subplot that spans 2 columns in the 1st column - they don't overlap, but we get weird sizes

        divider = make_axes_locatable(ax[0]); #prep to add an axis
        dividerKeo = make_axes_locatable(ax[1]); #prep to add an axis

        #mill for square Mercator style
        #robin for oval shape
        #npstere/spstere for polar
        if( ('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'north') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.80); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.80); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.min(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
                #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(90,np.min(settings['map']['lat range']),-15),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(90,np.min(settings['map']['lat range']),-15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.min(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j

            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoCircle.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.92, bottom = 0.045); #sets padding to small numbers for minimal white space
        elif( ('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'south') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.80); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.80); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.max(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                 resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),2), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(-90,np.max(settings['map']['lat range']),15),0), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(-90,np.max(settings['map']['lat range']),15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.max(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoMap.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
            
            fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.96, bottom = 0.045); #sets padding to small numbers for minimal white space
        else:
            cax = divider.append_axes('right', size='1.5%', pad=0.35); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.65); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
                resolution = 'i', area_thresh = 10000, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
                llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
                urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,np.int64(settings['map']['long autotick']*2)),2), 
                labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
                labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('auto');
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
            
            fig.subplots_adjust(left = 0.04, right = 0.95, top = 0.96, bottom = 0.045); #sets padding to small numbers for minimal white space
        #END IF
        
        #now draw some coastlines
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #now for the other plot (dividerKeo, ax[1])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[1].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo = dividerKeo.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        
        caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        figFitter(fig, fast=False); #fit that fig slow, cause it's gonna be used lots
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
        
    if( (settings['movie']['movie type'] == 7) ): #creates figure with 1 plot (geomap) and 2 colorbars (left is AMPERE) and smaller plot next to geomap as well for a keogram or something
        fig = plt.figure(figsize = settings['movie']['fig size']/settings['movie']['fig PPI']); #use fig because we gotta init an axis later
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        ax = [plt.subplot2grid((1, 3), (0, 0), colspan=2,fig=fig) , plt.subplot2grid((9, 3), (0, 2),rowspan=4, colspan=1,fig=fig), plt.subplot2grid((9, 3), (5, 2),rowspan=4, colspan=1,fig=fig) ]; #init an axis
        #the subplot2grid is weird, 2nd one makes a 7x3 grid and puts a subplot that spans 3 of the 7 rows 0 rows down and in the 3rd column, 3rd does similar but 5 rows down
        #1st one does a 1x3 grid and puts a subplot that spans 2 columns in the 1st column - they don't overlap, but we get weird sizes

        divider = make_axes_locatable(ax[0]); #prep to add an axis
        dividerKeo = make_axes_locatable(ax[1]); #prep to add an axis
        dividerKeo2 = make_axes_locatable(ax[2]); #prep to add an axis

        #mill for square Mercator style
        #robin for oval shape
        #npstere/spstere for polar
        if( ('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'north') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.60); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.60); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.min(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
                #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(90,np.min(settings['map']['lat range']),-15),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(90,np.min(settings['map']['lat range']),-15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.min(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j

            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoCircle.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            # fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.935, bottom = 0.072); #sets padding to small numbers for minimal white space
        elif(('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'south') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.60); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.60); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.max(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                 resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),2), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(-90,np.max(settings['map']['lat range']),15),0), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(-90,np.max(settings['map']['lat range']),15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.max(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoMap.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            # fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.96, bottom = 0.063); #sets padding to small numbers for minimal white space
        else:
            cax = divider.append_axes('right', size='1.5%', pad=0.15); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.90); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
                resolution = 'i', area_thresh = 10000, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
                llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
                urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,np.int64(settings['map']['long autotick']*2)),2), 
                labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
                labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('auto');
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
            
            fig.subplots_adjust(wspace=0.30); #pads the wspace to keep the labels from overlapping
            # fig.subplots_adjust(left = 0.045, right = 0.947, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #END IF
        
        #now draw some coastlines
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #now for the other plot (dividerKeo, ax[1])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[1].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo = dividerKeo.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        
        caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                
        #now for the other other plot (dividerKeo2, ax[2])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[2].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[2].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo2 = dividerKeo2.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        
        caxKeo2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        figFitter(fig, fast=False); #fit that fig slow, cause it's gonna be used lots
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
        
    if( (settings['movie']['movie type'] == 71) ): #creates figure with 1 plot (geomap) and 2 colorbars (left is AMPERE) and smaller plot next to geomap as well for a keogram or something
        fig = plt.figure(figsize = settings['movie']['fig size']/settings['movie']['fig PPI']); #use fig because we gotta init an axis later
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        ax = [plt.subplot2grid((1, 3), (0, 0), colspan=2,fig=fig) , plt.subplot2grid((9, 3), (0, 2),rowspan=4, colspan=1,fig=fig), plt.subplot2grid((9, 3), (5, 2),rowspan=4, colspan=1,fig=fig) ]; #init an axis
        #the subplot2grid is weird, 2nd one makes a 7x3 grid and puts a subplot that spans 3 of the 7 rows 0 rows down and in the 3rd column, 3rd does similar but 5 rows down
        #1st one does a 1x3 grid and puts a subplot that spans 2 columns in the 1st column - they don't overlap, but we get weird sizes

        divider = make_axes_locatable(ax[0]); #prep to add an axis
        dividerKeo = make_axes_locatable(ax[1]); #prep to add an axis
        dividerKeo2 = make_axes_locatable(ax[2]); #prep to add an axis

        #mill for square Mercator style
        #robin for oval shape
        #npstere/spstere for polar
        if( ('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'north') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.60); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.60); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.min(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
                #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(90,np.min(settings['map']['lat range']),-15),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(90,np.min(settings['map']['lat range']),-15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.min(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j

            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoCircle.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            # fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.935, bottom = 0.072); #sets padding to small numbers for minimal white space
        elif(('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'south') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.60); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.max(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                 resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),2), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(-90,np.max(settings['map']['lat range']),15),0), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(-90,np.max(settings['map']['lat range']),15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.max(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoMap.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            # fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.96, bottom = 0.063); #sets padding to small numbers for minimal white space
        else:
            cax = divider.append_axes('right', size='1.5%', pad=0.15); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            if( (np.diff(settings['map']['long range']).item() > 40 ) | (np.diff(settings['map']['lat range']).item() > 40 ) ):
                # geoMap_res = 
                geoMap_thresh = 10000;
            else:
                geoMap_thresh = 1000;
            #END IF
            
            geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
                resolution = 'i', area_thresh = geoMap_thresh, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
                llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
                urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,np.int64(settings['map']['long autotick']*2)),2), 
                labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
                labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('auto');
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
            
            # fig.subplots_adjust(wspace=0.30); #pads the wspace to keep the labels from overlapping
            fig.subplots_adjust(wspace=0.60); #pads the wspace to keep the labels from overlapping
            # fig.subplots_adjust(left = 0.045, right = 0.947, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #END IF
        
        #now draw some coastlines
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #now for the other plot (dividerKeo, ax[1])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[1].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo = dividerKeo.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        
        caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                
        #now for the other other plot (dividerKeo2, ax[2])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[2].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[2].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo2 = dividerKeo2.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        
        caxKeo2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        figFitter(fig, fast=False); #fit that fig slow, cause it's gonna be used lots
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
        
    elif( (settings['movie']['movie type'] == 8) ): #creates figure with 2 plots (TEC on top, time series on bottom)
        
        fig, ax = plt.subplots(figsize = settings['movie']['fig size']/settings['movie']['fig PPI'],nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        divider = make_axes_locatable(ax[0]); #prep to add an axis (top one is TEC so it gets it)
        cax = divider.append_axes('right', size='1.5%', pad=0.35); #make a color bar axis
        cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        #mill for square Mercator style
        #robin for oval shape
        geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
            resolution = 'i', area_thresh = 10000, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
            llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
            urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[0].set_aspect('auto');
        ax[1].set_aspect('auto');
        
        #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
        geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,settings['map']['long autotick']),2), 
            labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
        geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
            labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
        #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
        #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
        string_title = 'Prep Title'; #create mecha title
        ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the upper title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the lower title
        
        figFitter(fig, fast=False); #fit that fig slow, cause it's gonna be used lots
        # fig.subplots_adjust(left = 0.045, right = 0.945, top = 0.96, bottom = 0.065, hspace = 0.15); #sets padding to small numbers for minimal white space
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
        
    if( (settings['movie']['movie type'] == 9) ): #creates figure with 1 plot (geomap) and 2 colorbars (left is AMPERE) AND 1 plot time series (OMNI data)
    
        fig, ax = plt.subplots(figsize = settings['movie']['fig size']/settings['movie']['fig PPI'],nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        divider = make_axes_locatable(ax[0]); #prep to add an axis
        cax = divider.append_axes('right', size='1.5%', pad=0.35); #make a color bar axis
        cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        cax2 = divider.append_axes('left', size='1.5%', pad=0.85); #make a color bar axis
        cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        #mill for square Mercator style
        #robin for oval shape
        geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
            resolution = 'i', area_thresh = 10000, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
            llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
            urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[0].set_aspect('auto');
        ax[1].set_aspect('auto');
        
        #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
        geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,settings['map']['long autotick']),2), 
            labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
        geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
            labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
        #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
        #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
        string_title = 'Prep Title'; #create mecha title
        ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the upper title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the lower title
        
        figFitter(fig, fast=False); #fit that fig slow, cause it's gonna be used lots
        # fig.subplots_adjust(left = 0.045, right = 0.945, top = 0.96, bottom = 0.065, hspace = 0.15); #sets padding to small numbers for minimal white space
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
        
    if( (settings['movie']['movie type'] == 10) ): #creates figure with 1 plot (geomap) and 2 colorbars (left is AMPERE) and smaller plot next to geomap as well for a keogram or something
        fig = plt.figure(figsize = settings['movie']['fig size']/settings['movie']['fig PPI']); #use fig because we gotta init an axis later
        figManager = fig.canvas.manager; #req to maximize
        figManager.window.showMaximized(); #force maximized
        ax = [plt.subplot2grid((1, 3), (0, 0), colspan=2,fig=fig) , plt.subplot2grid((9, 3), (0, 2),rowspan=4, colspan=1,fig=fig), plt.subplot2grid((9, 3), (5, 2),rowspan=4, colspan=1,fig=fig) ]; #init an axis
        #the subplot2grid is weird, 2nd one makes a 7x3 grid and puts a subplot that spans 3 of the 7 rows 0 rows down and in the 3rd column, 3rd does similar but 5 rows down
        #1st one does a 1x3 grid and puts a subplot that spans 2 columns in the 1st column - they don't overlap, but we get weird sizes

        divider = make_axes_locatable(ax[0]); #prep to add an axis
        dividerKeo = make_axes_locatable(ax[1]); #prep to add an axis
        dividerKeo2 = make_axes_locatable(ax[2]); #prep to add an axis

        #mill for square Mercator style
        #robin for oval shape
        #npstere/spstere for polar
        if( ('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'north') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.60); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.60); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.min(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
                #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(90,np.min(settings['map']['lat range']),-15),0), 
                labels=[0,0,0,0], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(90,np.min(settings['map']['lat range']),-15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.min(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j

            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoCircle.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            # fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.935, bottom = 0.072); #sets padding to small numbers for minimal white space
        elif(('stere' in settings['map']['projection name']) & (settings['map']['hemi'] == 'south') ):
            cax = divider.append_axes('right', size='1.5%', pad=0.60); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.60); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'],boundinglat=np.max(settings['map']['lat range']),lon_0=np.mean(settings['map']['long range']),
                 resolution = 'i', area_thresh = 10000, ax=ax[0], round=True);
            if( settings['movie']['day nite line'] == 1 ): # == 2 is OK as it's better than my rough calcs
                settings['movie']['day nite line'] = 0; #override, they do weird things on polar
                settings['movie']['day nite text'] = 0; #override, they do weird things on polar
            #END IF
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,30),2), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(-90,np.max(settings['map']['lat range']),15),0), 
                labels=[True,True,True,True], labelstyle='+/-', color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
        
            for j in np.arange(0,360,30): #longitude labels
                x = (1.05*0.5*np.sin(np.deg2rad(j)))+0.5; #geoMap coordinate 
                y = (1.05*0.5*np.cos(np.deg2rad(j+180)))+0.5;
                if( j > 180 ):
                    angle = j-360; #deg, flip to negative
                else:
                    angle = j; #deg, angle is OK
                #END IF
                if( angle == 180):
                    y = y - 0.01; #small nudge, cause this one is too close
                #END IF
                if( angle == -60):
                    x = x - 0.003; #small nudge, cause this one is too close
                #END IF
                ax[0].text(x,y,str(angle)+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            latPts = np.roll(np.arange(-90,np.max(settings['map']['lat range']),15),-1); #degc, latitude points to note (this extra work is to replace the 90 w/ the last latitude value)
            latPts[-1] = np.max(settings['map']['lat range']); #degc, last is the final latitude value
            latPts_mapped = geoMap(np.tile(180,(latPts.size,)),latPts); #convert to geoMap values
            for j in range(0,latPts.size): #latitude labels
                x = latPts_mapped[0][j]/latPts_mapped[1][-1] + 0.025; #geoMap coordinate 
                y = latPts_mapped[1][j]/latPts_mapped[1][-1] - 0.0085; #geoMap coordinate 
                ax[0].text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',transform=ax[0].transAxes,horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
                #ax.text(x,y,str(latPts[j])+'\N{DEGREE SIGN}',horizontalalignment='center',verticalalignment='center',fontproperties=settings['plot']['font axis tick FM'])
            #END FOR j
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('equal');
            geoCircle = geoMap.drawmapboundary(linewidth=2, color='k'); #polar circle is clipped, so this draws it then makes sure it isn't
            geoMap.set_clip_on(False); #prevent weird things where the circle is clipped off
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM'],y=1.035); #set the title
            
            # fig.subplots_adjust(left = 0.02, right = 0.94, top = 0.96, bottom = 0.063); #sets padding to small numbers for minimal white space
        else:
            cax = divider.append_axes('right', size='1.5%', pad=0.15); #make a color bar axis
            cax.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            cax2 = divider.append_axes('left', size='1.5%', pad=0.90); #make a color bar axis
            cax2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
            
            geoMap = Basemap(projection=settings['map']['projection name'], lat_0=np.mean(settings['map']['lat range']), lon_0=np.mean(settings['map']['long range']), #projection type, and I think lat_0/lon_0 are the centers?
                resolution = 'i', area_thresh = 10000, ax=ax[0], #resolutions I know are l, i, h - i seems good. area_thresh being big prevents it drawing lil lakes, 0.1 makes everything
                llcrnrlon=np.float32(settings['map']['long range'][0]), llcrnrlat=np.float32(settings['map']['lat range'][0]), #lower left corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
                urcrnrlon=np.float32(settings['map']['long range'][1]), urcrnrlat=np.float32(settings['map']['lat range'][1])); #upper right corner lat/long - MUST BE FLOAT cause CODED BY THE BEST
        
            #np.linspace(startlat,endlat,5) # 5 = number of "ticks"
            geoMap.drawmeridians( np.round(np.arange(np.floor(np.min(settings['map']['long range'])),np.ceil(np.max(settings['map']['long range']))+1,np.int64(settings['map']['long autotick']*2)),2), 
                labels=[True,False,False,True], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            geoMap.drawparallels(np.round(np.arange(np.floor(np.min(settings['map']['lat range'])),np.ceil(np.max(settings['map']['lat range']))+1,settings['map']['lat autotick']),2), 
                labels=[True,False,True,False], labelstyle='+/-', dashes=[6,15000], color='black' ); #adds the labels but keeps the lines invisible
            #dashes[6,15000] seems to be enough to keep the repeating dash from happening on a world plot (I think it's how wide a dash is and the dash spacing). 900 was good for 22 degc longitude, for ref.
            #labels=[left,right,top,bottom] is the order. true/false to turn them on and off. not sure why there's a top/bottom for the parallels but w/e people do it
            
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[0].set_aspect('auto');
            
            string_title = 'Prep Title'; #create mecha title
            ax[0].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
            
            fig.subplots_adjust(wspace=0.30); #pads the wspace to keep the labels from overlapping
            # fig.subplots_adjust(left = 0.045, right = 0.947, top = 0.96, bottom = 0.065); #sets padding to small numbers for minimal white space
        #END IF
        
        #now draw some coastlines
        geoMap.drawcoastlines(zorder=25); #always on top
        if( settings['map']['world color'] == 1 ): #color in stuff if this is on
            #geoMap.drawcountries();
            #geoMap.drawmapboundary();
            geoMap.fillcontinents(color=settings['map']['land color'],lake_color=settings['map']['water color'],zorder=1);
            geoMap.drawmapboundary(fill_color=settings['map']['water color'],zorder=0);
        #END IF
        
        #now for the other plot (dividerKeo, ax[1])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[1].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[1].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo = dividerKeo.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        caxKeo.set_visible(False); #mkae it invisible so it matches the other plots in width
        
        # caxKeo.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
                
        #now for the other other plot (dividerKeo2, ax[2])
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax[2].set_aspect('auto');
        
        string_title = 'Prep Title'; #create mecha title
        ax[2].set_title(string_title,fontproperties=settings['plot']['font title FM']); #set the title
        caxKeo2 = dividerKeo2.append_axes('right', size='3%', pad=0.30); #make a color bar axis
        
        caxKeo2.yaxis.label.set_font_properties(settings['plot']['font axis label FM']);
        
        
        
        #update settings['movie']['fig size'] to now be the size of the movie plot (it'll be smaller if there are other plots on the screen)
        bbox = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted()); #get the size of the movie plot currently
        width, height = bbox.width, bbox.height; #break out the width/height in inches or something
        settings['movie']['fig size'] = np.int64(np.round(np.array( (width*fig.dpi, height*fig.dpi) ))); #save that size, convert from inches to pixels using DPI
    #END IF
    
    
    return movie_dict #holds everything that needs to leave
#END DEF