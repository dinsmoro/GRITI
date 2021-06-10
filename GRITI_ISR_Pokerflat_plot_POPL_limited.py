"""
GOAL: Plot only ISR POPL limited RTI
RD on 4/11/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_ISR_Pokerflat_plot_POPL_limited(PFISR_time,PFISR_height,PFISR_POPL,PFISR_el,PFISR_az,ISR_cutoffPeriod,ISR_RTI_heightLimValues,ISR_POPL_plotLimValu,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):
    
    for j in range(0,len(PFISR_time)):
    
        #-----Plot ISR POPL results as a RTI-----
        #Plot just the ISR POPL results
        fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
        figManager = plt.get_current_fig_manager(); #req to maximize
        figManager.window.showMaximized(); #force maximized
        divider = make_axes_locatable(ax); #prep to add an axis
        cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax.set_aspect('auto');
        
        pltHelprX, pltHelprY = np.meshgrid( (PFISR_time[j] - dateRange_dayNum_zeroHr[1])*24, \
                    PFISR_height[j]);
        im = ax.pcolormesh(pltHelprX , pltHelprY , PFISR_POPL[j].T , vmin=np.min(ISR_POPL_plotLimValu) , vmax=np.max(ISR_POPL_plotLimValu) , cmap='gray' , shading='gouraud'); # pseudocolor plot "stretched" to the grid
        cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
        cbar.set_clim(vmin=np.min(ISR_POPL_plotLimValu), vmax=np.max(ISR_POPL_plotLimValu));
        #cbar.ax.set_yticklabels(np.round(np.linspace(-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu,5), len(str(ISR_POPL_plotLimValu).split('.')[1])+1 )); #create useful tick marks
        cbar.set_label("POPL [$e^-/m^3$]"); #tabel the colorbar
        cbar.ax.tick_params(labelsize=FONT_axisTick);
        cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
        #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
        #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
        #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
        
        string_Title = 'PFISR POPL - No Filters, Amplitude Limited | El = '+str(np.median(PFISR_el[j]))+' & Az = '+str(np.median(PFISR_az[j])); #create mecha title
        ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    #    ax.set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    #    ax.set_ylabel('Height (km)',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        xAxisTicks = np.arange( (np.round((np.min(PFISR_time[j])-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(PFISR_time[j])-dateRange_dayNum_zeroHr[1])*24),2)) , \
                (np.round((np.max(PFISR_time[j])-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(PFISR_time[j])-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
                2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
        ax.set_xticks(xAxisTicks); #set x axis ticks
        ax.set_xlim( ((np.min(PFISR_time[j])-dateRange_dayNum_zeroHr[1])*24 , (np.max(PFISR_time[j])-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
        
        ax.set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
                
        ax.set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' (hr)',fontproperties=FONT_axisLabelFM); #set the x axis label
        ax.set_ylabel('Height (km)',fontproperties=FONT_axisLabelFM); #set the y axis label
        
        fig.subplots_adjust(left = 0.045, right = 0.95, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
        #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
        plt.show(); #req to make plot show up
    
    #END FOR j