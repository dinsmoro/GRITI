
import numpy as np
from scipy import stats
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from matplotlib.cm import get_cmap
from Code.GRITI_plotHelper_axisizerTime import GRITI_plotHelper_axisizerTime
from Code.subfun_strfind import strfind
from Code.subfun_sunAlsoRises import sunAlsoRises
from Code.subfun_textNice import textNice
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum
from Code.subfun_figFitter import figFitter

#corrRet is a list of corrRets that correspond to each data1types string
#data1types are a list of names that are the same length as corrRet

#data2types is a list of names that are the same length as every sub-list of corrRet for each data1types

def subfun_correlator_walking_plotter(corrRet, data1types, data2types, \
                                      time_cutout_range_walking, time2span, \
                                      time2step, time2lim, time2shiftDir, \
                                      settings_plot, settings_paths, dates, \
                                      FLG_showNiteTimes = False, showNiteTimesDict = None, \
                                      reportDivisor = [3600,'hr'], redZone=None, FLG_fancyPlot = False):
    
    #--- Declare & Unpack ---
    if( FLG_fancyPlot >= 1 ):
        print('MAKING FANCY PLOT: corr_walking_timeRange_'+str(data1types).replace('[\'','').replace('\']','').replace('\', \'',' & ').replace(' ','_')+'_data2typeHERE'+' IN fancyPlot FOLDER'); #report since you won't see anything
    #END IF
    # PLOT_lineWidthThicc = settings_plot['line width']['thicc']; #get the line widths
    # PLOT_lineWidthDoublePlus = settings_plot['line width']['double plus']; #get the line widths
    # PLOT_lineWidthPlus = settings_plot['line width']['plus']; #get the line widths
    # PLOT_lineWidthRegularPlus = settings_plot['line width']['regular plus']; #get the line widths
    # PLOT_lineWidthRegular = settings_plot['line width']['regular']; #get the line widths
    # PLOT_lineWidthSmol = settings_plot['line width']['smol']; #get the line widths
    # PLOT_lineWidthSmoller = settings_plot['line width']['smoller']; #get the line widths
    FONT_titleFM = settings_plot['font title FM'];
    # FONT_axisTick = settings_plot['font axis tick'];
    FONT_axisLabelFM = settings_plot['font axis label FM'];
    journal_dpi = settings_plot['journal dpi'];
    
    data2types_num = len(data2types); #number of data 2 (2nd set of inputs into function) types involved
    
    txt_every_few = (time_cutout_range_walking[:,0][-1] - time_cutout_range_walking[:,0][0])/time2step;
    txt_every_few = np.ceil(txt_every_few/30).astype(np.int64); #get spacing for text later on, assume 40 texts are the max
    
    #pre-extract the data
    corrRet_combo = [None for i in range(0,len(corrRet))]; #prep list holder
    timeShift_combo = [None for i in range(0,len(corrRet))]; #prep list holder
    for i in range(0,len(corrRet)):
        corrRet_combo[i] = np.empty( (data2types_num, len(corrRet[i])) ); #prep
        timeShift_combo[i] = np.empty( (data2types_num, len(corrRet[i])) );
        for j in range(0,data2types_num):
            for k in range(0,len(corrRet[i])):
                corrRet_combo[i][j, k] = corrRet[i][k][j]['corr'].item(); #get the value out
                timeShift_combo[i][j, k] = corrRet[i][k][j]['time shift'].item();
            #END FOR j
        #END FOR k
    #END FOR i
    
    for j in range(0,data2types_num):
        #--- Prep Plot ---
        if( FLG_fancyPlot == 0 ):
            # fig, ax = plt.subplots(nrows=len(corrRet), ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
            fig = plt.figure(); #fire up a figure only
            griddr = GridSpec(len(corrRet), 51, figure=fig);
            figManager = fig.canvas.manager; #req to maximize
            figManager.window.showMaximized(); #force maximized
            
            ax = [None for i in range(0,len(corrRet))]; #prep list
            for i in range(0,len(corrRet)):
                ax[i] = fig.add_subplot(griddr[i, 0:49]); #PLOTS FOR THE PLOT GOD
            #END FOR i
            cax = fig.add_subplot(griddr[:, 50]); #cax for the colorbar
        else:
            plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
            # fig, ax = plt.subplots(nrows=len(corrRet), ncols=1,figsize=(14,8.5),dpi=journal_dpi); #use instead of fig because it inits an axis too (I think I dunno)
            fig = plt.figure(figsize=(14,8.5), dpi=journal_dpi); #fire up a figure only
            griddr = GridSpec(len(corrRet), 51, figure=fig);
            ax = [None for i in range(0,len(corrRet))]; #prep list
            for i in range(0,len(corrRet)):
                ax[i] = fig.add_subplot(griddr[i, 0:49]); #PLOTS FOR THE PLOT GOD
            #END FOR i
            cax = fig.add_subplot(griddr[:, 50]); #cax for the colorbar
        #END IF
        fig.subplots_adjust(hspace=0.12);
        
        for i in range(0,len(corrRet)):
            #Remove the aspect ratio from the basemap so it fills the screen better
            ax[i].set_aspect('auto');
                    
            if( time2shiftDir == 'both' ):
                walking_vmin = -time2lim/60;
                walking_vmax = time2lim/60;
                walking_cmap = 'bwr';
            elif( time2shiftDir == 'pos' ):
                walking_vmin = 0;
                walking_vmax = time2lim/60;
                walking_cmap = ListedColormap( get_cmap('bwr')(np.linspace(0.5, 1, 300)) ); #just the red side
            else:
                walking_vmin = -time2lim/60;
                walking_vmax = 0;
                walking_cmap = ListedColormap( get_cmap('bwr')(np.linspace(0, 0.5, 300)) ); #just the blue side
            #END IF
            
            #--- Actual Plotting ---
            #first is to plot the corr coeff
            h_lgnd, = ax[i].plot(time_cutout_range_walking[:,0]/reportDivisor[0], corrRet_combo[i][j], linewidth=settings_plot['line width']['thicc'], color=settings_plot['color'][i], zorder=3);
            #second is to represent the time shift
            scat = ax[i].scatter(time_cutout_range_walking[:,0]/reportDivisor[0], corrRet_combo[i][j], s=settings_plot['line width']['thicc']*40, c=timeShift_combo[i][j]/60, cmap=walking_cmap, edgecolors='black', vmin=walking_vmin, vmax=walking_vmax, zorder=4);
            txt_every_few_cntr = 1; #prep cntr
            for k in range(1,timeShift_combo[i][j].size-1):
                if( txt_every_few == txt_every_few_cntr ): #only go a few times
                    if( corrRet_combo[i][j][k] >= 0 ):
                        va_pos = 'bottom'; #U+1F97A
                    else:
                        va_pos = 'top';
                    #END IF
                    ax[i].text(time_cutout_range_walking[:,0][k]/reportDivisor[0], corrRet_combo[i][j][k], textNice(np.round(timeShift_combo[i][j][k]/60,2)), ha='center', va=va_pos, zorder=5); #add the text on
                    txt_every_few_cntr = 1; #reset cntr
                else:
                    txt_every_few_cntr += 1; #increment cntr otherwise
                #END IF
            #END FOR k
    
            #--- extra bonus shading to show weak stuff ---
            if( redZone == None ):
                samples = np.diff(corrRet[i][0][j]['time range']).item()//corrRet[i][0][j]['data rate']; #calc # of samples
                confidence = 0.99; #use this confidence interval
                tval_crit = stats.t.ppf(1-(1-confidence)/2, samples-2); #equiv to Excel T.INV.2T(alpha, dof) or TINV(alpha, dof) (alpha = 1-confidence, dof = samples-2)
                redZone_now = tval_crit/np.sqrt(samples - 2 + tval_crit**2); #corr coeff values below this are insignificant            
            else:
                redZone_now = redZone; #use it manually to represent insiginificant region
            #END IF
            ax[i].axhspan(-redZone_now, redZone_now, facecolor='xkcd:brick red', alpha=0.5, zorder=1); #makes it easier to see what is a bad val
            #--- force -1 to 1 limts ---
            ax[i].set_ylim( (-1, 1) );
            
            #----- Plot Shading to represent 'night' -----
            if( FLG_showNiteTimes == True ):
                niteTimes = showNiteTimesDict['nite times']; #unpack
                if( niteTimes[0] == False ):
                    niteTimes_latAlign = showNiteTimesDict['lat align'][i];
                    niteTimes_latLongRange = showNiteTimesDict['lat long range'][i];
                    if( showNiteTimesDict['lat or long'][i].lower() == 'latitude' ): #if true, latitude
                        if( niteTimes_latAlign > 45 ):
                            latToUse = 45; #cap at 50 deg lat to keep it from getting too zesty at the poles
                        else:
                            latToUse = niteTimes_latAlign;
                        #END IF
                        (niteTimes_sunRise, niteTimes_sunSet, dateRange_fullPad) = sunAlsoRises(dates['date range full'],latToUse,np.mean(niteTimes_latLongRange[1])); #call sunrise/set function
                    else:
                        (niteTimes_sunRise, niteTimes_sunSet, dateRange_fullPad) = sunAlsoRises(dates['date range full'],np.mean(niteTimes_latLongRange[0]),niteTimes_latAlign); #call sunrise/set function
                    #END IF
                    dateRange_dayNum_fullPad = subfun_date_to_dayNum(dateRange_fullPad); #convert to dayNum
                    niteTimes_sunRise = (niteTimes_sunRise + dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*24; #hrs, center around zero hr and convert ot hrs
                    niteTimes_sunSet = (niteTimes_sunSet + dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*24; #hrs, center around zero hr and convert ot hrs
                    niteTimes_sunRise = niteTimes_sunRise[1:]; #remove 1st
                    niteTimes_sunSet = niteTimes_sunSet[:-1]; #remove last
                    #FIFTH STEP: PLOT THIS STUFF
                    for jk in range(0,niteTimes_sunSet.size):
                        ax[i].axvspan(niteTimes_sunSet[jk], niteTimes_sunRise[jk], edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.25, zorder=2); #make that patch
                    #END FOR jk
                else:
                    for jk in range(0,dates['date range zero hr hours'].size):
                        niteTimes_per = dates['date range zero hr hours'][jk] + np.asarray(niteTimes[i]); #get the nite time ranges
                        ax[i].axvspan(niteTimes_per[0], niteTimes_per[-1], edgecolor='xkcd:black', facecolor='xkcd:black', alpha=0.25, zorder=2); #make that patch
                    #END FOR jk
                #END IF
            #END IF
            
            #--- Axis and Titles and Stuff ---
            if( FLG_fancyPlot == 0 ): #only title non-fancy plot
                string_Title = 'Corr. Coeff Time Span '+textNice(time2span/reportDivisor[0])+' '+reportDivisor[1]+ \
                    ' | Time Increment '+textNice(time2step/reportDivisor[0])+' '+reportDivisor[1]+ \
                    ' | Lag Limit '+textNice(time2lim/reportDivisor[0])+' '+reportDivisor[1];
                ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            #END IF
            
            ax[i].set_ylabel(data1types[i],fontproperties=FONT_axisLabelFM); #set the y axis label
            
            #nice axis ticks
            if( FLG_fancyPlot == 0 ):
                if( i != (len(corrRet)-1) ):
                    GRITI_plotHelper_axisizerTime(time_cutout_range_walking[:,0]/reportDivisor[0],ax=ax[i],unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True);
                else:
                    GRITI_plotHelper_axisizerTime(time_cutout_range_walking[:,0]/reportDivisor[0],ax=ax[i],unit=reportDivisor[1],FLG_removeLabels=False,FLG_tickDirIn=True);
                #END IF
            else:
                if( i != (len(corrRet)-1) ):
                    GRITI_plotHelper_axisizerTime(time_cutout_range_walking[:,0]/reportDivisor[0],ax=ax[i],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=False,FLG_tickDirIn=True);
                else:
                    GRITI_plotHelper_axisizerTime(time_cutout_range_walking[:,0]/reportDivisor[0],ax=ax[i],unit=reportDivisor[1],tickNumGoal=23,FLG_removeLabels=False,FLG_tickDirIn=True);
                #END IF
            #END IF
            
            #legend creation
            #&| means sub-dict, $| means sub-array
            if( (strfind(data2types[j],'&|',opt=1) > 0) | (strfind(data2types[j],'$|',opt=1) > 0) ):
                data2_adder = '$\mathregular{_F}$'; #default to F
                data2_setNames_subSegs = data2types[j].split('&|'); #get the sub-dicts
                for d2s in range(0,len(data2_setNames_subSegs)):
                    if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) != 0 ): #check for sub-array for the sub-dict
                        data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                        if( data2_setNames_subSegs_subVects[1] == '0' ):
                            data2_adder = '$\mathregular{_X}$';
                            data2_namer = '-X';
                        elif( data2_setNames_subSegs_subVects[1] == '1' ):
                            data2_adder = '$\mathregular{_Y}$';
                            data2_namer = '-Y';
                        elif( data2_setNames_subSegs_subVects[1] == '2' ):
                            data2_adder = '$\mathregular{_Z}$';
                            data2_namer = '-Z';
                        elif( data2_setNames_subSegs_subVects[1] == '3' ):
                            data2_adder = '$\mathregular{_T}$'; #(By^2 + Bz^2)^(1/2) - not implemented yet
                            data2_namer = '-T';
                        #END IF
                    #END IF
                #END FOR d2s
                ax[i].legend([h_lgnd], [data2_setNames_subSegs[0]+data2_adder], loc='center right');
                data2types_name = data2_setNames_subSegs[0]+data2_namer; #build it
                data2types_name = data2types_name.replace('[\'','').replace('\']','').replace('\', \'',' & ').replace(' ','_').replace('/', '-'); #fix it up
            else:
                ax[i].legend([h_lgnd], [data2types[j]], loc='center right');
                data2types_name = data2types[j].replace('[\'','').replace('\']','').replace('\', \'',' & ').replace(' ','_').replace('/', '-'); #use it
            #END IF
            
            ax[i].text( np.min(ax[i].get_xlim())+np.diff(ax[i].get_xlim())*.005, np.max(ax[i].get_ylim()), chr(97+i)+'.', color='r', fontproperties=settings_plot['font grandiose FM'], horizontalalignment='left', verticalalignment='top'); #print the text labelling the lettering a. b. c. ect.
        #END FOR i
        
        cbar_autoTick = (np.ceil(walking_vmax) - np.floor(walking_vmin))/13; #tries to split the time range into tickNumGoal # of times
        if( cbar_autoTick > 100*.9 ):
            cbar_autoTick = 100;
        elif( cbar_autoTick > 60*.9 ):
            cbar_autoTick = 60;
        elif( cbar_autoTick > 50*.9 ):
            cbar_autoTick = 50;
        elif( cbar_autoTick > 30*.9 ):
            cbar_autoTick = 30;
        elif( cbar_autoTick > 20*.9 ):
            cbar_autoTick = 20;
        elif( cbar_autoTick > 15*.9 ):
            cbar_autoTick = 15;
        elif( cbar_autoTick > 10*.9 ):
            cbar_autoTick = 10;
        elif( cbar_autoTick > 5*.9 ):
            cbar_autoTick = 5;
        elif( cbar_autoTick > 3*.9 ):
            cbar_autoTick = 3;
        elif( cbar_autoTick > 2*.9 ):
            cbar_autoTick = 2;
        elif( cbar_autoTick > 1*.9 ):
            cbar_autoTick = 1;
        else:
            cbar_autoTick = (walking_vmax - walking_vmin)/13; #just goes for it if it's a super tiny range
        #END IF
        cbar_tickz = np.arange( (np.round(walking_vmin) - np.mod(np.round(walking_vmin),cbar_autoTick)) , \
                (np.round(walking_vmax) - np.mod(np.round(walking_vmax),cbar_autoTick)) + cbar_autoTick , \
                cbar_autoTick); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
        fig.colorbar(scat, cax=cax, ticks=cbar_tickz, label='Lag [min]'); #add in colorbar
        
        ax[-1].set_xlabel('Time ['+reportDivisor[1]+'] | Zero Hr on '+ \
            dates['date range zero hr month name']+' '+str(dates['date range zero hr'][2])+dates['date range zero hr day post fix']+', '+str(dates['date range zero hr'][0]),fontproperties=FONT_axisLabelFM); #set the x axis label
        
        figFitter(fig); #fit the fig fast
        if( FLG_fancyPlot != 0 ):
            fig.savefig(os.path.join(settings_paths['fancyPlots'], 'Correlations', 'corr_walking_timeRange_'+str(data1types).replace('[\'','').replace('\']','').replace('\', \'',' & ').replace(' ','_').replace('/', '-')+'_'+data2types_name+settings_plot['save file type'])); #save the figure
            plt.close(); #close figure b/c it lurks apparently
            plt.ion(); #re-enable it for later stuff
        #END IF
    #END FOR j
#END DEF