"""
GOAL: Plot only ISR ION VELOCITY RTI
RD on 8/13/19

INPUT: buncha stuff
OUTPUT: no vars, just a plot is made
"""

import numpy as np #import in here I dunno
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Code.subfun_monthNum_to_word import subfun_monthNum_to_word

def GRITI_ISR_Haystack_plot_ionVel(Zenith_time,Zenith_height,Zenith_vel,MISA_time,MISA_height,MISA_vel,ISR_RTI_heightLimValues,ISR_ionVel_plotLimValu_noFilt,dateRange,dateRange_dayNum,dateRange_dayNum_zeroHr,FONT_titleFM,FONT_axisTick,FONT_axisLabelFM):
    #-----Plot ISR POPL HP results as a RTI-----
    #Plot just the ISR POPL HP results
    fig, ax = plt.subplots(nrows=2, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    divider = make_axes_locatable(ax[0]); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax[0].set_aspect('auto');
    ax[1].set_aspect('auto');
    
    pltHelprX, pltHelprY = np.meshgrid( (Zenith_time - dateRange_dayNum_zeroHr[1])*24, \
                Zenith_height);
    im = ax[0].pcolormesh(pltHelprX , pltHelprY , Zenith_vel , vmin=np.min(ISR_ionVel_plotLimValu_noFilt) , vmax=np.max(ISR_ionVel_plotLimValu_noFilt) , cmap='gray' , shading='gouraud'); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar.mappable.set_clim(vmin=np.min(ISR_ionVel_plotLimValu_noFilt), vmax=np.max(ISR_ionVel_plotLimValu_noFilt));
    #cbar.ax.set_yticklabels(np.round(np.linspace(-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu,5), len(str(ISR_POPL_plotLimValu).split('.')[1])+1 )); #create useful tick marks
    cbar.set_label("Ion Velocity [$m/s$]"); #tabel the colorbar
    cbar.ax.tick_params(labelsize=FONT_axisTick);
    cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    string_Title = 'Zenith LOS Ion Velocity'; #create mecha title
    ax[0].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[0].set_xlabel(subfun_monthNum_to_word(dateRange[0,1])[0]+" "+str(dateRange[0,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[0,0])+" to "+subfun_monthNum_to_word(dateRange[1,1])[0]+" "+str(dateRange[1,2])+" (Day "+str(dateRange_dayNum[0,1])+"), "+str(dateRange[1,0]),fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[0].set_ylabel('Height (km)',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[0].set_xticks(xAxisTicks); #set x axis ticks
    ax[0].set_xlim( ((np.min(Zenith_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(Zenith_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    ax[0].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    
    #MISA STUFF
    divider = make_axes_locatable(ax[1]); #prep to add an axis
    cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
    
    pltHelprX, pltHelprY = np.meshgrid( (MISA_time - dateRange_dayNum_zeroHr[1])*24, \
                MISA_height);
    im = ax[1].pcolormesh(pltHelprX , pltHelprY , MISA_vel , vmin=np.min(ISR_ionVel_plotLimValu_noFilt) , vmax=np.max(ISR_ionVel_plotLimValu_noFilt) , cmap='gray', shading='gouraud'); # pseudocolor plot "stretched" to the grid
    cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
    cbar.mappable.set_clim(vmin=np.min(ISR_ionVel_plotLimValu_noFilt), vmax=np.max(ISR_ionVel_plotLimValu_noFilt));
    #cbar.ax.set_yticklabels(np.round(np.linspace(-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu,5), len(str(ISR_POPL_plotLimValu).split('.')[1])+1 )); #create useful tick marks
    cbar.set_label("Ion Velocity [$m/s$]"); #tabel the colorbar
#    cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); #force a rounded format
    cbar.ax.tick_params(labelsize=FONT_axisTick);
    cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
    #cbar.ax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cax.set_ylim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    #cbar.set_clim( (-ISR_POPL_plotLimValu,ISR_POPL_plotLimValu) ); #create useful tick marks
    
    string_Title = 'MISA LOS Ion Velocity'; #create mecha title
    ax[1].set_title(string_Title,fontproperties=FONT_titleFM); #set the title
    ax[1].set_xlabel('Time in UT - 0 Hr on Day '+str(dateRange_dayNum_zeroHr[1])+', '+str(dateRange_dayNum_zeroHr[0])+' (hr)',fontproperties=FONT_axisLabelFM); #set the x axis label
    ax[1].set_ylabel('Height (km)',fontproperties=FONT_axisLabelFM); #set the y axis label
    
    xAxisTicks = np.arange( (np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) , \
            (np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) - np.mod(np.round((np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24),2)) + 2 , \
            2); #sets the start hr, stop hr, and the step size between (in this case, 2 hr)
    ax[1].set_xticks(xAxisTicks); #set x axis ticks
    ax[1].set_xlim( ((np.min(MISA_time)-dateRange_dayNum_zeroHr[1])*24 , (np.max(MISA_time)-dateRange_dayNum_zeroHr[1])*24) ); #set x axis limits
    
    ax[1].set_ylim( ISR_RTI_heightLimValues ); #set y axis limits
    
    fig.subplots_adjust(left = 0.045, right = 0.95, top = 0.96, bottom = 0.065 , hspace = 0.225); #sets padding to small numbers for minimal white space
    #fig.tight_layout(); #function for a tight layout, doesn't seem to do much here
    plt.show(); #req to make plot show up