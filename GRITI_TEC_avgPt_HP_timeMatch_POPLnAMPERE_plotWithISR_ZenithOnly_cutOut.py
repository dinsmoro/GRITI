"""
GOAL: Plot only ISR POPL HP RTI
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_TEC_avgPt_HP_timeMatch_POPLnAMPERE_plotWithISR_ZenithOnly_cutOut(avgPt_vTEC_timeMatch_HP_cutOut,avgPt_vTEC_timeMatch_time_cutOut, \
        Zenith_time_cutOut,Zenith_height,Zenith_POPL_hp_cutOut,Zenith_POPL_hp_altAvgd_cutOut,MISA_time_cutOut,MISA_height,MISA_POPL_hp_altAvgd_cutOut, \
        AMPERE_jouleHeating_integrate_cutOut, AMPERE_timeUnique_hr_cutOut, AMPERE_plot_label, AMPERE_plot_label_noUnits, FLG_AMPERE_log, FLG_AMPERE_hp, \
        filter_cutoffPeriod,ISR_RTI_heightLimValues,ISR_POPL_plotLimValu,avgPt_coords,avgPt_ISRavgAlt,pointAltitude,time_cutout_range, time_cutout_range_delay, \
        dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,PLOT_lineWidth):
    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth[0]; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth[1]; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth[2]; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth[3]; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth[4]; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth[5]; #get the line widths
    
    #-----Plot ISR POPL HP results as a RTI-----
    #Plot just the ISR POPL HP results
    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    
    avgPt_pwr = np.sqrt(1/avgPt_vTEC_timeMatch_HP_cutOut.size*np.sum(avgPt_vTEC_timeMatch_HP_cutOut**2)); #estimate power of signal
    Zenith_pwr = np.sqrt(1/Zenith_POPL_hp_altAvgd_cutOut.size*np.sum(Zenith_POPL_hp_altAvgd_cutOut**2)); #estimate power of signal
    MISA_pwr = np.sqrt(1/MISA_POPL_hp_altAvgd_cutOut.size*np.sum(MISA_POPL_hp_altAvgd_cutOut**2)); #estimate power of signal
    AMPERE_pwr = np.sqrt(1/AMPERE_jouleHeating_integrate_cutOut.size*np.sum(AMPERE_jouleHeating_integrate_cutOut**2)); #estimate power of signal
    line1, = ax[0].plot( (avgPt_vTEC_timeMatch_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , 1/avgPt_pwr*avgPt_vTEC_timeMatch_HP_cutOut , color='xkcd:fire engine red' ,linewidth=PLOT_lineWidthRegularPlus, linestyle='-', marker='o'); # plot some time series
    line2, = ax[0].plot( (Zenith_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , 1/Zenith_pwr*Zenith_POPL_hp_altAvgd_cutOut, color='xkcd:electric blue' ,linewidth=PLOT_lineWidthRegularPlus, linestyle='--'); # plot some time series
    line3, = ax[0].plot( (MISA_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , 1/MISA_pwr*MISA_POPL_hp_altAvgd_cutOut, color='xkcd:vivid green' ,linewidth=PLOT_lineWidthRegularPlus, linestyle ='-.'); # plot some time series
    line4, = ax[0].plot( AMPERE_timeUnique_hr_cutOut , 1/AMPERE_pwr*AMPERE_jouleHeating_integrate_cutOut , color='xkcd:purple' ,linewidth=PLOT_lineWidthDoublePlus, linestyle =':'); # plot some time series
    
    AMPERE_legend = AMPERE_plot_label_noUnits; #prep legend
    if( FLG_AMPERE_log == 1 ):
        AMPERE_legend = 'Log10 '+AMPERE_legend; #add log10 info
    #END IF
    if( FLG_AMPERE_hp == 1 ):
        AMPERE_legend = AMPERE_legend + ' High-passed'; #add high-pass info
    #END IF
    if( time_cutout_range_delay != 0 ):
        AMPERE_legend = AMPERE_legend + ' (time delayed '+str(time_cutout_range_delay)+' hr)'; #add time delay info
    #END IF
    ax[0].legend((line1, line4, line2, line3), ('delta-vTEC ISR Interval Matched AVG''d High-passed',AMPERE_legend,'Zenith POPL High-passed AVG''d', 'MISA POPL High-passed AVG''d') ,loc='lower center', ncol=2);
    string_Title = 'delta-vTEC ISR Interval Matched AVG''d High-passed AND Zenith/MISA POPL HP AVG''d +/-'+str(avgPt_ISRavgAlt)+' km around '+str(pointAltitude)+' km at '+str(np.round(avgPt_coords[0],2))+' deg lat/'+str(np.round(avgPt_coords[1],2))+'deg long'; #create mecha title
    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('Amplitude [Norm\'d Power]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            1); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    ax[0].set_xlim( ((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    #ZENITH RTI PLOT   
#    divider = make_axes_locatable(ax[1]); #prep to add an axis
#    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    pltHelprX, pltHelprY = np.meshgrid( (Zenith_time_cutOut - dateRange_dayNum_zeroHr[1])*24, \
                Zenith_height);
    im = ax[1].pcolormesh(pltHelprX , pltHelprY , Zenith_POPL_hp_cutOut , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu) , cmap='gray', shading='gouraud'); # pseudocolor plot "stretched" to the grid
#    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
#    cbar.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
#    #cbar.ax.set_yticklabels(np.round(np.linspace(np.min(ISR_POPL_plotLimValu),np.max(ISR_POPL_plotLimValu),5), len(str(np.min(ISR_POPL_plotLimValu)).split('.')[1])+1 )); #create useful tick marks
#    #cbar.set_label("POPL (unitless)"); #tabel the colorbar
#    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
#    cbar.ax.tick_params(labelsize=FONT_axisTick);
#    cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
#    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
#    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
#    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    string_Title = 'Zenith POPL with High-pass Filter '+str(filter_cutoffPeriod)+' hr Cutoff - Hours '+str(np.min(time_cutout_range))+' to '+str(np.max(time_cutout_range)); #create mecha title
    ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[1].set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' (hr)',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel('Height (km)',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            1); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks
    ax[1].set_xlim( ((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    ax[1].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    
    fig.subplots_adjust(left = 0.055, right = 0.98, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up