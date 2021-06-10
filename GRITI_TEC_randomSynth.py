"""
Replaces delta-vTEC data with random data OR with random data and an embedded synthetic wave (or two)
"""
import numpy as np

def GRITI_TEC_randomSynth(deltavTECsize,TEC_lat,TEC_long,TEC_time, \
    noise_background_mean,noise_background_stdev,Re,dateRange_zeroHr, \
    plotLatRange,plotLongRange,plotLatRange_autoTick,plotLongRange_autoTick, \
    wave_latRange,wave_longRange,wave_N,wave_angle,wave_phase,wave_waveLength,wave_period,wave_amp, \
    FONT_titleFM,FONT_axisTick,FONT_axisLabelFM,TEC_plotLimValu,FLG_TEC_noise,FLG_plotStuff=0):
    
    #prep to plot, if the colorbar limit is 1 value, make it 2 because it's meant to be a +/-# situation.
    if( np.isscalar(TEC_plotLimValu) == 1 ):
        TEC_plotLimValu = np.array( (-TEC_plotLimValu,TEC_plotLimValu) ); #make it a vector
    #END IF
    
    #replace TEC values with noise, keep the time/location the same
    if( FLG_TEC_noise == 1):
        #just flat noise, no synthetic wave injected into it
        deltavTECrando = np.random.normal(loc=noise_background_mean,scale=noise_background_stdev,size=deltavTECsize); #fill in random noise into the vTEC
    elif(FLG_TEC_noise == 2):
        #add a synth wave in
        # CALC N CONVERT AS NEEDED
        wave_latRange_lin = np.linspace(np.min(wave_latRange),np.max(wave_latRange),num=wave_N); #arcdeg, get a vector of pts along this
        wave_longRange_lin = np.linspace(np.min(wave_longRange),np.max(wave_longRange),num=wave_N); #arcdeg, get a vector of pts along this

        [wave_longRange_mesh, wave_latRange_mesh] = np.meshgrid(wave_longRange_lin,wave_latRange_lin); #arcdeg, meshgrid of X (long) and Y (lat) for plotting and stuff

        wave_angle = wave_angle*np.pi/180; #rad, convert
        wave_phase = wave_phase*np.pi/180; #rad, convert

        wave_waveLength_km = wave_waveLength; #km, record
        wave_waveLength = wave_waveLength/Re*(180/np.pi); #arcdeg, convert to arcdeg convention from km arc
        wave_period = wave_period/3600; #sec -> hr
        wave_freq = 1/wave_period; #1/hr, freq of wave
        wave_speed = wave_waveLength/wave_period; #arcdeg/hr
        wave_waveNumber = 2*np.pi/wave_waveLength; #rad/arcdeg, wave number
        wave_freqAngular = 2*np.pi*wave_freq; #rad/hr, angular freq
        
        if( FLG_plotStuff == 1 ):
            import matplotlib.pyplot as plt
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            
            # TEST WAVEFORM MADE
            grating = np.flipud( wave_amp[0]*np.sin( wave_waveNumber[0]*np.cos(wave_angle[0])*wave_longRange_mesh + wave_waveNumber[0]*np.sin(wave_angle[0])*wave_latRange_mesh + wave_freqAngular[0]*0 + wave_phase[0]) 
            + wave_amp[1]*np.sin( wave_waveNumber[1]*np.cos(wave_angle[1])*wave_longRange_mesh + wave_waveNumber[1]*np.sin(wave_angle[1])*wave_latRange_mesh + wave_freqAngular[1]*0 + wave_phase[1])     
            );
            
            #Plot the entire world to see
            fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
            figManager = plt.get_current_fig_manager(); #req to maximize
            figManager.window.showMaximized(); #force maximized
            divider = make_axes_locatable(ax); #prep to add an axis
            cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
            
            im = ax.pcolormesh( wave_longRange_lin , wave_latRange_lin , np.fliplr(grating), cmap='gray');# pseudocolor plot "stretched" to the grid
            cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
            cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
            cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
            cbar.ax.tick_params(labelsize=FONT_axisTick);
            cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
            string_Title = 'Synthetic Wave - ';
            if(wave_amp[1] == 0):
                string_Title = string_Title+' WaveParams:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr';
            else:
                string_Title = string_Title+' WaveParams: W1:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr & W2:'+str(wave_angle[1]*180/np.pi)+'deg+'+str(wave_waveLength_km[1])+'km+'+str(wave_period[1])+'hr';
            #END IF
            ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            ax.set_xlabel('Longitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the x axis label
            ax.set_ylabel('Latitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
            xAxisTicks = np.arange( -180,180+15,15 ); #preps the x ticks
            ax.set_xticks(xAxisTicks); #set x axis ticks
            yAxisTicks = np.arange( -90,90+15,15 ); #preps the x ticks
            ax.set_yticks(yAxisTicks); #set x axis ticks
    
            #plot a local shot
            fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
            figManager = plt.get_current_fig_manager(); #req to maximize
            figManager.window.showMaximized(); #force maximized
            divider = make_axes_locatable(ax); #prep to add an axis
            cax = divider.append_axes('right', size='2.0%', pad=0.35); #make a color bar axis
            
            im = ax.pcolormesh( wave_longRange_lin , wave_latRange_lin , np.fliplr(grating), cmap='gray');# pseudocolor plot "stretched" to the grid
            cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
            cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
            cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
            cbar.ax.tick_params(labelsize=FONT_axisTick);
            cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
            string_Title = 'Synthetic Wave - ';
            if(wave_amp[1] == 0):
                string_Title = string_Title+' WaveParams:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr';
            else:
                string_Title = string_Title+' WaveParams: W1:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr & W2:'+str(wave_angle[1]*180/np.pi)+'deg+'+str(wave_waveLength_km[1])+'km+'+str(wave_period[1])+'hr';
            #END IF
            ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            ax.set_xlabel('Longitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the x axis label
            ax.set_ylabel('Latitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
            xAxisTicks = np.arange(np.min(plotLongRange),np.max(plotLongRange)+plotLongRange_autoTick,plotLongRange_autoTick); #creates x ticks automagically
            yAxisTicks = np.arange(np.min(plotLatRange),np.max(plotLatRange)+plotLatRange_autoTick,plotLatRange_autoTick); #creates y ticks automagically
            ax.set_xticks(xAxisTicks); #set x axis ticks
            ax.set_yticks(yAxisTicks); #set y axis ticks
            ax.set_xlim( [np.min(plotLongRange), np.max(plotLongRange)] ); #limit the x to the plot range
            ax.set_ylim( [np.min(plotLatRange), np.max(plotLatRange)] ); #limit the y to the plot range
    
    
            #Plot just world with TEC noise
            fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
            figManager = plt.get_current_fig_manager(); #req to maximize
            figManager.window.showMaximized(); #force maximized
            # imagesc( grating + normrnd(noise_background_mean,noise_Background_STDEV,size(grating)) , [-1 1] );
            im = ax.pcolormesh( wave_longRange_lin , wave_latRange_lin , np.fliplr(grating + np.random.normal(loc=noise_background_mean,scale=noise_background_stdev,size=np.shape(grating))), cmap='gray');# pseudocolor plot "stretched" to the grid
            cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
            cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
            cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
            cbar.ax.tick_params(labelsize=FONT_axisTick);
            cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
            string_Title = 'Synthetic Wave with Noise - ';
            if(wave_amp[1] == 0):
                string_Title = string_Title+' WaveParams:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr';
            else:
                string_Title = string_Title+' WaveParams: W1:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr & W2:'+str(wave_angle[1]*180/np.pi)+'deg+'+str(wave_waveLength_km[1])+'km+'+str(wave_period[1])+'hr';
            #END IF
            string_Title = string_Title+ ' - NoiseParams: Normal+Mean:'+str(np.round(noise_background_mean+2))+'+Stdev:'+str(np.round(noise_background_stdev+2))+'';
            ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            ax.set_xlabel('Longitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the x axis label
            ax.set_ylabel('Latitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
            #local shot
            fig, ax = plt.subplots(nrows=1, ncols=1); #use instead of fig because it inits an axis too (I think I dunno)
            figManager = plt.get_current_fig_manager(); #req to maximize
            figManager.window.showMaximized(); #force maximized
            # imagesc( grating + normrnd(noise_background_mean,noise_Background_STDEV,size(grating)) , [-1 1] );
            im = ax.pcolormesh( wave_longRange_lin , wave_latRange_lin , np.fliplr(grating + np.random.normal(loc=noise_background_mean,scale=noise_background_stdev,size=np.shape(grating))), cmap='gray');# pseudocolor plot "stretched" to the grid
            cbar = fig.colorbar(im, cax=cax, orientation='vertical'); #create a colorbar using the prev. defined cax
            cbar.mappable.set_clim(vmin=np.min(TEC_plotLimValu), vmax=np.max(TEC_plotLimValu));
            cbar.set_label("delta-vTEC [TECU]"); #tabel the colorbar
            cbar.ax.tick_params(labelsize=FONT_axisTick);
            cax.yaxis.label.set_font_properties(FONT_axisLabelFM);
            string_Title = 'Synthetic Wave with Noise - ';
            if(wave_amp[1] == 0):
                string_Title = string_Title+' WaveParams:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr';
            else:
                string_Title = string_Title+' WaveParams: W1:'+str(wave_angle[0]*180/np.pi)+'deg+'+str(wave_waveLength_km[0])+'km+'+str(wave_period[0])+'hr & W2:'+str(wave_angle[1]*180/np.pi)+'deg+'+str(wave_waveLength_km[1])+'km+'+str(wave_period[1])+'hr';
            #END IF
            string_Title = string_Title+ ' - NoiseParams: Normal+Mean:'+str(np.round(noise_background_mean+2))+'+Stdev:'+str(np.round(noise_background_stdev+2))+'';
            ax.set_title(string_Title,fontproperties=FONT_titleFM); #set the title
            ax.set_xlabel('Longitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the x axis label
            ax.set_ylabel('Latitude [arcdeg]',fontproperties=FONT_axisLabelFM); #set the y axis label
            xAxisTicks = np.arange(np.min(plotLongRange),np.max(plotLongRange)+plotLongRange_autoTick,plotLongRange_autoTick); #creates x ticks automagically
            yAxisTicks = np.arange(np.min(plotLatRange),np.max(plotLatRange)+plotLatRange_autoTick,plotLatRange_autoTick); #creates y ticks automagically
            ax.set_xticks(xAxisTicks); #set x axis ticks
            ax.set_yticks(yAxisTicks); #set y axis ticks
            ax.set_xlim( [np.min(plotLongRange), np.max(plotLongRange)] ); #limit the x to the plot range
            ax.set_ylim( [np.min(plotLatRange), np.max(plotLatRange)] ); #limit the y to the plot range
        #END IF
        
        #DANGER ADD SEPERATELY IMPLEMENTED BECAUSE RANDOM NUMBERS JUST WEREN'T BEING ADDED IF ALL IN ONE GO (WHY?)
        deltavTECrando = wave_amp[0]*np.sin( wave_waveNumber[0]*np.cos(wave_angle[0])*TEC_long + wave_waveNumber[0]*np.sin(wave_angle[0])*TEC_lat + wave_freqAngular[0]*(TEC_time - dateRange_zeroHr[1]*86400)/3600 + wave_phase[0]) 
        + wave_amp[1]*np.sin( wave_waveNumber[1]*np.cos(wave_angle[1])*TEC_long + wave_waveNumber[1]*np.sin(wave_angle[1])*TEC_lat + wave_freqAngular[1]*(TEC_time - dateRange_zeroHr[1]*86400)/3600 + wave_phase[1])     
        deltavTECrando = deltavTECrando + np.random.normal(loc=noise_background_mean,scale=noise_background_stdev,size=(deltavTECsize,)); #delta_vTEC, calculate synthetic vTEC
    #END IF
    

    return deltavTECrando