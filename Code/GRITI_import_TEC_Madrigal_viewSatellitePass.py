import numpy as np
import matplotlib.pyplot as plt
from subfun_figFitter import subfun_figFitter

def GRITI_import_TEC_Madrigal_viewSatellitePass(currentsTEC_singleSatLine, currentvTEC_singleSatLine, currentPolyYvals, current_deltavTEC, \
        currentElv_singleSatLine, currentSat_singleSatLine, currentSite_singleSatLine, currentYear_singleSatLine, currentDayNum_singleSatLine, \
        currentHour_singleSatLine, currentMin_singleSatLine, currentSec_singleSatLine, currentLat_singleSatLine, currentLong_singleSatLine, \
        TEC_deltaLim, minElevation, folder, savedVersion = False):
    
    #----- Interactive Version -----
    fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
    figManager = plt.get_current_fig_manager(); #req to maximize
    figManager.window.showMaximized(); #force maximized
    
    #twinx time
    axr = ax.twinx(); #get that other x
    
    #Remove the aspect ratio from the basemap so it fills the screen better
    ax.set_aspect('auto');
    axr.set_aspect('auto');
    
    custElv = np.copy(currentElv_singleSatLine); #copy it over
    custElv[np.where(currentElv_singleSatLine == np.max(currentElv_singleSatLine))[0].item():] = custElv[np.where(currentElv_singleSatLine == np.max(currentElv_singleSatLine))[0].item():] + 2*(np.max(currentElv_singleSatLine) - currentElv_singleSatLine[np.where(currentElv_singleSatLine == np.max(currentElv_singleSatLine))[0].item():]); #make it always increase
    custElv -= np.max(currentElv_singleSatLine); #0 it
    p1, = ax.plot( custElv, currentsTEC_singleSatLine ,color='xkcd:turquoise' ,linewidth=3.00, linestyle='-',label='sTEC'); #plot
    p2, = ax.plot( custElv , currentvTEC_singleSatLine ,color='xkcd:fire engine red' ,linewidth=3.00, linestyle='--',label='vTEC'); #plot
    p3, = ax.plot( custElv , currentPolyYvals ,color='xkcd:electric blue' ,linewidth=3.00, linestyle=':',label='Savitzky-Golay Fit'); #plot
    
    p4, = axr.plot( custElv , current_deltavTEC ,color='xkcd:purple' ,linewidth=3.00, linestyle='-.',label='delta-vTEC'); #plot
    axr.set_ylim( (-TEC_deltaLim,TEC_deltaLim) ); 
    axr.set_yticks( np.arange(-TEC_deltaLim,TEC_deltaLim+.25,.25) );
    axr.set_ylabel('delta-vTEC [TECU]'); #set the y axis label
    
    ax.set_xlabel("Elevation Angle [deg]"); #set the x axis label
    ax.set_ylabel('TEC [TECU]'); #set the y axis label
    ax.legend((p1, p2, p3, p4), ('sTEC', 'vTEC', 'Savitzky-Golay Fit', 'delta-vTEC')); #throw in the legend #loc="center right"
    
    ax.set_xlim( (np.ceil(minElevation-np.max(currentElv_singleSatLine)),np.floor(np.abs(minElevation-np.max(currentElv_singleSatLine)))) ); 
    ax.set_xticks( np.arange(np.ceil(minElevation-np.max(currentElv_singleSatLine)),np.floor(np.abs(minElevation-np.max(currentElv_singleSatLine)))+10,10) );
    
    #switch out the xlabels for other labels (elevation labels!)
    xlabels = ax.get_xticklabels(); #get the xlabels
    xlabelsTru = np.hstack( (np.arange(minElevation,np.round(np.max(currentElv_singleSatLine),-1),10),np.round(np.float64(np.max(currentElv_singleSatLine)),2),np.flip(np.arange(minElevation,np.round(np.max(currentElv_singleSatLine),-1),10))) ); #make the true labels
    for jk in range(0,len(xlabels)):
        xlabels[jk].set_text(str(xlabelsTru[jk]).rstrip('0').rstrip('.'))
    #END FOR jk
    ax.set_xticklabels(xlabels); #get the xlabels
    
    string_Title = 'Satellite G'+str(currentSat_singleSatLine[0])+' over site '+currentSite_singleSatLine[0].decode('UTF-8')+ \
        ' for '+str(currentYear_singleSatLine[0])+', '+str(currentDayNum_singleSatLine[0])+' from '+ \
        str(currentHour_singleSatLine[0])+':'+str(currentMin_singleSatLine[0])+':'+str(currentSec_singleSatLine[0])+' to '+ \
        str(currentHour_singleSatLine[-1])+':'+str(currentMin_singleSatLine[-1])+':'+str(currentSec_singleSatLine[-1])+ \
        ' at '+str(np.round(currentLat_singleSatLine[currentElv_singleSatLine == np.max(currentElv_singleSatLine)].item(),2))+' lat, '+ \
        str(np.round(currentLong_singleSatLine[currentElv_singleSatLine == np.max(currentElv_singleSatLine)].item(),2))+' long'; #create mecha title
    ax.set_title(string_Title); #set the title
    print(string_Title); #print it also
    
    subfun_figFitter(fig); #fit that fig fast
    plt.show();
    
    
    if( savedVersion == True ):
        #----- Non-Interactive, Saved Version -----
        plt.ioff(); #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        fig, ax = plt.subplots(figsize=(14,8.5),dpi=300); #use instead of fig because it inits an axis too (I think I dunno)
        
        #twinx time
        axr = ax.twinx(); #get that other x
        
        #Remove the aspect ratio from the basemap so it fills the screen better
        ax.set_aspect('auto');
        axr.set_aspect('auto');
        
        custElv = np.copy(currentElv_singleSatLine); #copy it over
        custElv[np.where(currentElv_singleSatLine == np.max(currentElv_singleSatLine))[0].item():] = custElv[np.where(currentElv_singleSatLine == np.max(currentElv_singleSatLine))[0].item():] + 2*(np.max(currentElv_singleSatLine) - currentElv_singleSatLine[np.where(currentElv_singleSatLine == np.max(currentElv_singleSatLine))[0].item():]); #make it always increase
        custElv -= np.max(currentElv_singleSatLine); #0 it
        p1, = ax.plot( custElv, currentsTEC_singleSatLine ,color='xkcd:turquoise' ,linewidth=3.00, linestyle='-',label='sTEC'); #plot
        p2, = ax.plot( custElv , currentvTEC_singleSatLine ,color='xkcd:fire engine red' ,linewidth=3.00, linestyle='--',label='vTEC'); #plot
        p3, = ax.plot( custElv , currentPolyYvals ,color='xkcd:electric blue' ,linewidth=3.00, linestyle=':',label='Savitzky-Golay Fit'); #plot
        
        p4, = axr.plot( custElv , current_deltavTEC ,color='xkcd:purple' ,linewidth=3.00, linestyle='-.',label='delta-vTEC'); #plot
        axr.set_ylim( (-TEC_deltaLim,TEC_deltaLim) ); 
        axr.set_yticks( np.arange(-TEC_deltaLim,TEC_deltaLim+.25,.25) );
        axr.set_ylabel('delta-vTEC [TECU]'); #set the y axis label
        
        ax.set_xlabel("Elevation Angle [deg]"); #set the x axis label
        ax.set_ylabel('TEC [TECU]'); #set the y axis label
        ax.legend((p1, p2, p3, p4), ('sTEC', 'vTEC', 'Savitzky-Golay Fit', 'delta-vTEC')); #throw in the legend #loc="center right"
        
        ax.set_xlim( (np.ceil(minElevation-np.max(currentElv_singleSatLine)),np.floor(np.abs(minElevation-np.max(currentElv_singleSatLine)))) ); 
        ax.set_xticks( np.arange(np.ceil(minElevation-np.max(currentElv_singleSatLine)),np.floor(np.abs(minElevation-np.max(currentElv_singleSatLine)))+10,10) );
        
        #switch out the xlabels for other labels (elevation labels!)
        xlabels = ax.get_xticklabels(); #get the xlabels
        xlabelsTru = np.hstack( (np.arange(minElevation,np.round(np.max(currentElv_singleSatLine),-1),10),np.round(np.float64(np.max(currentElv_singleSatLine)),2),np.flip(np.arange(minElevation,np.round(np.max(currentElv_singleSatLine),-1),10))) ); #make the true labels
        for jk in range(0,len(xlabels)):
            xlabels[jk].set_text(str(xlabelsTru[jk]).rstrip('0').rstrip('.'))
        #END FOR jk
        ax.set_xticklabels(xlabels); #get the xlabels
        
        # string_Title = 'Satellite G'+str(currentSat_singleSatLine[0])+' over site '+currentSite_singleSatLine[0].decode('UTF-8')+ \
        #     ' for '+str(currentYear_singleSatLine[0])+', '+str(currentDayNum_singleSatLine[0])+' from '+ \
        #     str(currentHour_singleSatLine[0])+':'+str(currentMin_singleSatLine[0])+':'+str(currentSec_singleSatLine[0])+' to '+ \
        #     str(currentHour_singleSatLine[-1])+':'+str(currentMin_singleSatLine[-1])+':'+str(currentSec_singleSatLine[-1])+ \
        #     ' at '+str(np.round(currentLat_singleSatLine[currentElv_singleSatLine == np.max(currentElv_singleSatLine)].item(),2))+' lat, '+ \
        #     str(np.round(currentLong_singleSatLine[currentElv_singleSatLine == np.max(currentElv_singleSatLine)].item(),2))+' long'; #create mecha title
        # ax.set_title(string_Title); #set the title
        # print(string_Title); #print it also
        
        subfun_figFitter(fig); #fit that fig fast
        
        fig.savefig(folder[3]+'\\debug_import_TEC_Madrigal_viewSatellitePass.png'); #save the figure
        plt.close(); #close figure b/c it lurks apparently
        plt.ion(); #re-enable it for later stuff
    #END IF
#END DEF