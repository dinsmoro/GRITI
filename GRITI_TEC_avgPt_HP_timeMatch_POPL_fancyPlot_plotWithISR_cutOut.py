"""
GOAL: Plot only ISR POPL HP RTI
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_TEC_avgPt_HP_timeMatch_POPL_fancyPlot_plotWithISR_cutOut(avgPt_vTEC_timeMatch_HP_cutOut,avgPt_vTEC_timeMatch_time_cutOut, \
        Zenith_time_cutOut,Zenith_height,Zenith_POPL_hp_cutOut,Zenith_POPL_hp_altAvgd_cutOut,MISA_time_cutOut,MISA_height, \
        MISA_POPL_hp_cutOut,MISA_POPL_hp_altAvgd_cutOut,filter_cutoffPeriod,ISR_RTI_heightLimValues,ISR_POPL_plotLimValu, \
        avgPt_coords,avgPt_ISRavgAlt,pointAltitude,time_cutout_range,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,
        dateRange_zeroHr, dateRange_zeroHr_monthName, dateRange_zeroHr_dayPostfix, \
        FONT_grandioseFM, FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,PLOT_lineWidth, journal_width_2C,journal_height_max,journal_dpi):
    print('MAKING FANCY PLOT: TEC_avgPt_HP_timeMatch_POPL_fancyPlot_plotWithISR_cutOut IN fancyPlot FOLDER'); #report since you won't see anything

    ISR_m3tocc = 100**3; #1 m^3 is 100^3 cm^3 and 1 cm^3 is 1 cc
    ISR_POPL_plotLimValu = ISR_POPL_plotLimValu/ISR_m3tocc; #adjust this plot limit value
    Zenith_POPL_hp_cutOut = Zenith_POPL_hp_cutOut/ISR_m3tocc; #adjust POPL values
    Zenith_POPL_hp_altAvgd_cutOut = Zenith_POPL_hp_altAvgd_cutOut/ISR_m3tocc; #adjust POPL values
    MISA_POPL_hp_cutOut = MISA_POPL_hp_cutOut/ISR_m3tocc; #adjust POPL values
    MISA_POPL_hp_altAvgd_cutOut = MISA_POPL_hp_altAvgd_cutOut/ISR_m3tocc; #adjust POPL values

    letteringPositionX = -0.105; #set the X position of the lettering (e.g., a. b. c. ...)
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
    
    ISR_RTI_heightLimValues = np.array(ISR_RTI_heightLimValues); #convert to array so it's not weak python tuples
    ISR_RTI_heightLimValues[0] = np.min(Zenith_height); #move the height limit to the actual lower-limit
    
    #Unpack line widths
    PLOT_lineWidthThicc = PLOT_lineWidth[0]; #get the line widths
    PLOT_lineWidthDoublePlus = PLOT_lineWidth[1]; #get the line widths
    PLOT_lineWidthPlus = PLOT_lineWidth[2]; #get the line widths
    PLOT_lineWidthRegularPlus = PLOT_lineWidth[3]; #get the line widths
    PLOT_lineWidthRegular = PLOT_lineWidth[4]; #get the line widths
    PLOT_lineWidthSmol = PLOT_lineWidth[5]; #get the line widths
    
    #-----Plot ISR POPL HP results as a RTI-----
    #Plot just the ISR POPL HP results
    plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
    fig = plt.figure(figsize=(14,10.5),dpi=journal_dpi);
    gridr = gridspec.GridSpec(nrows=3, ncols=1, figure=fig);
    gridr.update(hspace=0.15); # set the spacing between axes. 
    fig.add_subplot(gridr[0]); #RTI plots are 2 tall
    fig.add_subplot(gridr[1]); #RTI plots are 2 tall
    fig.add_subplot(gridr[2]); #dayNite plot is 1 tall
    ax = fig.axes; #get a list of the axes
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    ax[2].set_aspect('auto');
    
    #----PLOT TIME SERIES DATA----
    divider0 = make_axes_locatable(ax[0]); #prep to add an axis
    cax0 = divider0.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    cax0.set_visible(False); #mkae it invisible so it matches the other plots in width
    
    avgPt_pwr = np.sqrt(1/avgPt_vTEC_timeMatch_HP_cutOut.size*np.sum(avgPt_vTEC_timeMatch_HP_cutOut**2)); #estimate power of signal
    Zenith_pwr = np.sqrt(1/Zenith_POPL_hp_altAvgd_cutOut.size*np.sum(Zenith_POPL_hp_altAvgd_cutOut**2)); #estimate power of signal
    MISA_pwr = np.sqrt(1/MISA_POPL_hp_altAvgd_cutOut.size*np.sum(MISA_POPL_hp_altAvgd_cutOut**2)); #estimate power of signal
    line1, = ax[0].plot( (avgPt_vTEC_timeMatch_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , 1/avgPt_pwr*avgPt_vTEC_timeMatch_HP_cutOut , color='xkcd:fire engine red' ,linewidth=PLOT_lineWidthRegularPlus, linestyle='-', marker='o'); # plot some time series
    line2, = ax[0].plot( (Zenith_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , 1/Zenith_pwr*Zenith_POPL_hp_altAvgd_cutOut , color='xkcd:electric blue' ,linewidth=PLOT_lineWidthRegularPlus, linestyle='--'); # plot some time series
    line3, = ax[0].plot( (MISA_time_cutOut - dateRange_dayNum_zeroHr[1])*24 , 1/MISA_pwr*MISA_POPL_hp_altAvgd_cutOut , color='xkcd:vivid green' ,linewidth=PLOT_lineWidthRegularPlus, linestyle ='-'); # plot some time series
    
    ax[0].legend((line2, line3, line1), ('Zenith N$_{e^-}$', 'MISA N$_{e^-}$', 'delta-vTEC') ,loc='lower right', ncol=3, prop=FONT_axisLabelFM);
    # string_Title = 'delta-vTEC ISR Interval Matched AVG''d High-passed AND Zenith/MISA POPL HP AVG''d +/-'+str(avgPt_ISRavgAlt)+' km around '+str(pointAltitude)+' km at '+str(np.round(avgPt_coords[0],2))+' deg lat/'+str(np.round(avgPt_coords[1],2))+'deg long'; #create mecha title
    # ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
#    ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('Amplitude\n[Power Normalized]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            1); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    ax[0].set_xlim( ((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    ax[0].text( letteringPositionX, letteringPositionY, 'a.', color='r', fontproperties=FONT_grandioseFM, transform=ax[0].transAxes); #print the text saying the day or nite
    ax[0].spines['right'].set_visible(False); #turn off box lines
    ax[0].spines['top'].set_visible(False); #turn off box lines
    
    #-----ZENITH RTI PLOT-----
    divider1 = make_axes_locatable(ax[1]); #prep to add an axis
    cax1 = divider1.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    pltHelprX, pltHelprY = np.meshgrid( (Zenith_time_cutOut - dateRange_dayNum_zeroHr[1])*24, \
                Zenith_height);
    im1 = ax[1].pcolormesh(pltHelprX , pltHelprY , Zenith_POPL_hp_cutOut , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu) , cmap='gray', shading='gouraud'); # pseudocolor plot "stretched" to the grid
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar1.mappable.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
    cbar1.set_label("N$_{e^-}$ [$e^-/cc$]"); #tabel the colorbar
    cbar1.ax.tick_params(labelsize=FONT_axisTick);
    cax1.yaxis.label.set_font_properties(FONT_axisLabelFM);
    # cax1.yaxis.set_major_formatter(tick.FuncFormatter(cbarFormatter)); #force a rounded format
#    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
#    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
#    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    # string_Title = 'Zenith POPL with High-pass Filter '+str(filter_cutoffPeriod)+' hr Cutoff - Hours '+str(np.min(time_cutout_range))+' to '+str(np.max(time_cutout_range)); #create mecha title
    # ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
#    ax[1].set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' (hr)',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel('Height [km]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    # xAxisTicks = np.arange( (np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) , \
    #         (np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
    #         1); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks
    ax[1].set_xlim( ((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    ax[1].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    ax[1].text( letteringPositionX, letteringPositionY, 'b.', color='r', fontproperties=FONT_grandioseFM, transform=ax[1].transAxes); #print the text saying the day or nite
    
    #----MISA RTI PLOT-----
    divider2 = make_axes_locatable(ax[2]); #prep to add an axis
    cax2 = divider2.append_axes('right', size='2.0%', pad=0.15); #make a color bar axis
    pltHelprX, pltHelprY = np.meshgrid( (MISA_time_cutOut - dateRange_dayNum_zeroHr[1])*24, 
                MISA_height);
    im2 = ax[2].pcolormesh(pltHelprX , pltHelprY , MISA_POPL_hp_cutOut , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu) , cmap='gray', shading='gouraud'); # pseudocolor plot "stretched" to the grid
    cbar2 = fig.colorbar(im2, cax=cax2, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar2.mappable.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
    cbar2.set_label("N$_{e^-}$ [$e^-/cc$]"); #tabel the colorbar
    cbar2.ax.tick_params(labelsize=FONT_axisTick);
    cax2.yaxis.label.set_font_properties(FONT_axisLabelFM);
    # cax2.yaxis.set_major_formatter(tick.FuncFormatter(cbarFormatter)); #force a rounded format
#    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
#    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
#    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    # string_Title = 'MISA POPL with High-pass Filter '+str(filter_cutoffPeriod)+' hr Cutoff - Hours '+str(np.min(time_cutout_range))+' to '+str(np.max(time_cutout_range)); #create mecha title
    # ax[2].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[2].set_xlabel('Time in UT [hr] - 0 Hr on '+dateRange_zeroHr_monthName+' '+str(dateRange_zeroHr[2])+dateRange_zeroHr_dayPostfix+' | Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[2].set_ylabel('Height [km]',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    # xAxisTicks = np.arange( (np.round((np.min(MISA_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(MISA_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) , 
    #         (np.round((np.max(MISA_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(MISA_time_cutOut)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , 
    #         1); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[2].set_xticks(xAxisTicks); #set x axis ticks
    # ax[2].set_xlim( ((np.min(MISA_time_cutOut)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    ax[2].set_xlim( ((np.min(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time_cutOut)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits

    ax[2].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    ax[2].text( letteringPositionX, letteringPositionY, 'c.', color='r', fontproperties=FONT_grandioseFM, transform=ax[2].transAxes); #print the text saying the day or nite
    
    fig.subplots_adjust(left = 0.095, right = 0.910, top = 0.985, bottom = 0.070); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    fig.savefig('fancyPlots//avgPt_HP_timeMatch&POPL_cutOut.png'); #save the figure
    plt.close(); #close figure b/c it lurks apparently
    plt.ion(); #re-enable it for later stuff